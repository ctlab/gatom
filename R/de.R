convertPvalDT <- function(de.pvals, map) {
    from <- key(map)
    to <- setdiff(colnames(map), from)
    res.pvals <- map[de.pvals, nomatch=0]
    # setnames(res.pvals, c("ID", to), c(from, "ID"))
    res.pvals <- res.pvals[order(pval)]
    res.pvals <- res.pvals[!duplicated(res.pvals[, to, with=F])]
    setkeyv(res.pvals, to)
    setcolorder(res.pvals, c(to,
                             setdiff(colnames(res.pvals),
                                     c(to, from)),
                             from))
    res.pvals
}

prepareDEColumn <- function(gene.de, columnName, from) {
    if (is.character(from)) {
        setnames(gene.de, old=from, new=columnName)
    } else {
        gene.de[, (columnName) := eval(from, envir = gene.de)]
    }
}

#' Makes data.table with differential expression results containing
#' all columns required for gatom in the expected format
#' @export
prepareDE <- function(de.raw, de.meta) {
    if (is.null(de.raw)) {
        return(NULL)
    }
    if (!is(de.raw, "data.table")) {
        de.raw <- as.data.table(
            as.data.frame(de.raw),
            keep.rownames = !is.numeric(attr(de.raw,
                                             "row.names")))
    }

    de <- copy(de.raw)

    for (columnName in names(de.meta$columns)) {
        prepareDEColumn(de,
                        columnName,
                        de.meta$columns[[columnName]])
    }
    de[]
}

findColumn <- function(de, names) {
     candidates <-  na.omit(match(tolower(names),
                                  tolower(colnames(de))))
     if (length(candidates) == 0) {
         return(NA)
     }
     return(colnames(de)[candidates[1]])
}

findIdColumn <- function(de, idsList,
                         sample.size=1000,
                         match.threshold=0.6,
                         remove.ensembl.revisions=TRUE) {
    # first looking for column with base IDs
    de.sample <- if (nrow(de) < sample.size) {
        copy(de)
    } else {
        de[sample(seq_len(nrow(de)), sample.size), ]
    }
    columnSamples <- lapply(de.sample, as.character)


    if (remove.ensembl.revisions) {
        columnSamples <- lapply(columnSamples, gsub,
                                pattern="(ENS\\w*\\d*)\\.\\d*",
                                replacement="\\1")
    }

    ss <- sapply(columnSamples,
                 .intersectionSize, idsList[[1]])

    if (max(ss) / nrow(de.sample) >= match.threshold) {
        # we found a good column with base IDs
        return(list(column=colnames(de)[which.max(ss)],
                    type=names(idsList)[1]))
    }

    z <- .pairwiseCompare(.intersectionSize,
                         columnSamples,
                         idsList)

    bestMatch <- which(z == max(z), arr.ind = TRUE)[1,]
    return(list(column=colnames(de)[bestMatch["row"]],
                type=names(idsList)[bestMatch["col"]]))
}

idsListFromAnnotation <- function(org.gatom.anno) {
    res <- c(list(org.gatom.anno$genes$gene),
                  lapply(names(org.gatom.anno$mapFrom),
                         function(n) org.gatom.anno$mapFrom[[n]][[n]])
    )
    names(res) <- c(org.gatom.anno$baseId, names(org.gatom.anno$mapFrom))
    res
}

#' Finds columns in gene differential expression table
#' required for gatom analysis
#' @param baseMeanColumn could be NULL (automatic), NA (no such column),
#'                       character (coumn name)
#' @param signalColumn could be NULL (automatic), character (coumn name)
#'                    function (evaluated in a scope of original data frame)
#' @export
getGeneDEMeta <- function(gene.de.raw, org.gatom.anno,
                          idColumn=NULL,
                          idType=NULL,
                          pvalColumn=NULL,
                          logPvalColumn=NULL,
                          log2FCColumn=NULL,
                          baseMeanColumn=NULL,
                          signalColumn=NULL,
                          signalRankColumn=NULL
                          ) {

    if (is.null(idColumn) != is.null(idType)) {
        stop("Either both or none of idColumn and idType can be specified")
    }

    if (is.null(idColumn)) {
        idsList <- idsListFromAnnotation(org.gatom.anno)
        idColumnInfo <- findIdColumn(gene.de.raw, idsList)
        idColumn <- idColumnInfo$column
        idType <- idColumnInfo$type
    }

    if (is.null(pvalColumn)) {
        pvalColumn <- findColumn(gene.de.raw,
                                 c("pval", "p.value", "pvalue"))
    }

    if (is.null(logPvalColumn)) {
        logPvalColumn <- findColumn(gene.de.raw,
                                    c("logpval"))
        if (is.na(logPvalColumn)) {
            logPvalColumn <- quote(log(pval))
        }
    }

    if (is.null(log2FCColumn)) {
        log2FCColumn <- findColumn(gene.de.raw,
                                   c("log2FC", "log2foldchange",
                                     "logfc"))
    }

    if (is.null(baseMeanColumn)) {
        baseMeanColumn <- findColumn(gene.de.raw,
                                     c("baseMean", "aveexpr"))
    }

    if (is.null(signalColumn)) {
        signalColumn <- findColumn(gene.de.raw,
                                    c("signal", "probe"))
        if (is.na(signalColumn)) {
            signalColumn <- quote(paste0(pval, "_", log2FC))
        }
    }

    if (is.null(signalRankColumn)) {
        signalRankColumn <- findColumn(gene.de.raw,
                                  c("signalRank", "rank"))
        if (is.na(signalRankColumn)) {
            signalRankColumn <- quote({
                signalLevels <- setNames(baseMean, signal)[!duplicated(signal)]
                signalRanks <- setNames(rank(-signalLevels), names(signalLevels))
                signalRanks[signal]
                })
        }
    }


    list(idType=idType,
         columns=list(
             ID=idColumn,
             pval=pvalColumn,
             logPval=logPvalColumn,
             log2FC=log2FCColumn,
             baseMean=baseMeanColumn,
             signal=signalColumn,
             signalRank=signalRankColumn
             ))
}

#' Finds columns in differential expression table for metabolites
#' required for gatom analysis
#' @export
getMetDEMeta <- function(met.de.raw, met.db,
                         idColumn=NULL,
                         idType=NULL,
                         pvalColumn=NULL,
                         logPvalColumn=NULL,
                         log2FCColumn=NULL,
                         baseMeanColumn=NULL,
                         signalColumn=NULL,
                         signalRankColumn=NULL
) {
    if (is.null(met.de.raw)) {
        return(NULL)
    }

    if (is.null(idColumn) != is.null(idType)) {
        stop("Either both or none of idColumn and idType can be specified")
    }

    if (is.null(idColumn)) {
        idsList <- c(list(met.db$metabolites$metabolite),
                     lapply(names(met.db$mapFrom),
                            function(n) met.db$mapFrom[[n]][[n]])
        )
        names(idsList) <- c(met.db$baseId, names(met.db$mapFrom))
        idColumnInfo <- findIdColumn(met.de.raw, idsList)
        idColumn <- idColumnInfo$column
        idType <- idColumnInfo$type
    }

    if (is.null(pvalColumn)) {
        pvalColumn <- findColumn(met.de.raw,
                                 c("pval", "p.value", "pvalue"))
    }

    if (is.null(logPvalColumn)) {
        logPvalColumn <- findColumn(met.de.raw,
                                    c("logpval"))
        if (is.na(logPvalColumn)) {
            logPvalColumn <- quote(log(pval))
        }
    }

    if (is.null(log2FCColumn)) {
        log2FCColumn <- findColumn(met.de.raw,
                                   c("log2FC", "log2foldchange",
                                     "logfc"))
    }

    if (is.null(baseMeanColumn)) {
        baseMeanColumn <- findColumn(met.de.raw,
                                     c("baseMean", "aveexpr"))
    }

    if (is.null(signalColumn)) {
        signalColumn <- findColumn(met.de.raw,
                                  c("signal", "ion"))
        if (is.na(signalColumn)) {
            signalColumn <- quote(paste0(pval, "_", log2FC))
        }
    }

    if (is.null(signalRankColumn)) {
        signalRankColumn <- findColumn(met.de.raw,
                                      c("signalRank", "rank"))
        if (is.na(signalRankColumn)) {
            signalRankColumn <- quote({
                signalLevels <- setNames(baseMean, signal)[!duplicated(signal)]
                signalRanks <- setNames(rank(-signalLevels), names(signalLevels))
                signalRanks[signal]
            })
        }
    }


    list(idType=idType,
         columns=list(
             ID=idColumn,
             pval=pvalColumn,
             logPval=logPvalColumn,
             log2FC=log2FCColumn,
             baseMean=baseMeanColumn,
             signal=signalColumn,
             signalRank=signalRankColumn))
}
