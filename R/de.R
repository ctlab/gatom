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

prepareDeColumn <- function(gene.de, columnName, from) {
    if (is.character(from)) {
        setnames(gene.de, old=from, new=columnName)
    } else {
        gene.de[, (columnName) := eval(from, envir = gene.de)]
    }
}

#' Makes data.table with gene expression table containing
#' all columns required for gatom in the expected format
#' @export
prepareGeneDE <- function(gene.de.raw, gene.de.meta) {
    if (!is(gene.de.raw, "data.table")) {
        gene.de.raw <- as.data.table(
            as.data.frame(gene.de.raw),
            keep.rownames = !is.numeric(attr(de,
                                             "row.names")))
    }

    gene.de <- copy(gene.de.raw)

    for (columnName in names(gene.de.meta$columns)) {
        prepareDeColumn(gene.de,
                        columnName,
                        gene.de.meta$columns[[columnName]])
    }
    gene.de[]
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
        de[sample(seq_len(nrow(de)), sample.size)]
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

#' Finds columns in gene differential expression tables
#' required for gatom analysis
#' @param baseMeanColumn could be NULL (automatic), NA (no such column),
#'                       character (coumn name)
#' @param probeColumn could be NULL (automatic), character (coumn name)
#'                    function (evaluated in a scope of original data.frame)
#' @export
getGeneDEMeta <- function(gene.de.raw, org.gatom.anno,
                          idColumn=NULL,
                          idType=NULL,
                          pvalColumn=NULL,
                          logPvalColumn=NULL,
                          log2FCColumn=NULL,
                          baseMeanColumn=NULL,
                          probeColumn=NULL,
                          probeRankColumn=NULL
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

    if (is.null(probeColumn)) {
        probeColumn <- findColumn(gene.de.raw,
                                    c("probe"))
        if (is.na(probeColumn)) {
            probeColumn <- quote(paste0(pval, "_", log2FC))
        }
    }

    if (is.null(probeRankColumn)) {
        probeRankColumn <- findColumn(gene.de.raw,
                                  c("probeRank", "rank"))
        if (is.na(probeRankColumn)) {
            probeRankColumn <- quote({
                probeLevels <- setNames(baseMean, probe)[!duplicated(probe)]
                probeRanks <- setNames(rank(-probeLevels), names(probeLevels))
                probeRanks[probe]
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
             probe=probeColumn,
             probeRank=probeRankColumn))
}
