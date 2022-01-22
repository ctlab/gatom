splitIds <- function(ids, split=" */// *") {
    res <- strsplit(ids, split)
    res[lengths(res) == 0] <- ids[lengths(res) == 0]
    res
}

convertPvalDT <- function(de.pvals, map,
                          idColumn=ifelse(haskey(de.pvals), key(de.pvals), "ID"),
                          removeDuplicates=TRUE, removeGeneVersions=FALSE,
                          split=" */// *") {
    from <- key(map)
    to <- setdiff(colnames(map), from)

    if (length(split) != 0) {
        newIds <- splitIds(de.pvals[[idColumn]])

        if (any(lengths(newIds) != 1)) {
            de.pvals.new <- de.pvals[rep(seq_len(nrow(de.pvals)), lengths(newIds))]
            de.pvals.new[, c(idColumn) := unlist(newIds)]
            de.pvals <- de.pvals.new
        }
    }

    if (removeGeneVersions) {
        de.pvals <- copy(de.pvals)
        de.pvals[, c(idColumn) := sub(pattern="\\.\\d*$", replacement="", x=de.pvals[[idColumn]]) ]
    }

    res.pvals <- merge.data.table(map, de.pvals, by.x=key(map), by.y=idColumn)
    # setnames(res.pvals, c("ID", to), c(from, "ID"))
    res.pvals <- res.pvals[order(pval)]
    if (removeDuplicates) {
        res.pvals <- res.pvals[!duplicated(res.pvals[, to, with=F])]
    }
    setkeyv(res.pvals, to)
    setcolorder(res.pvals, c(to,
                             setdiff(colnames(res.pvals),
                                     c(to, from)),
                             from))
    res.pvals
}

prepareDEColumn <- function(gene.de, columnName, from) {
    if (is.character(from)) {
        if (columnName %in% colnames(gene.de) && from != columnName) {
            columnName1 <- columnName
            while (columnName1 %in% colnames(gene.de)) {
                columnName1 <- paste0(columnName1, "_")
            }
            setnames(gene.de, old=columnName, new=columnName1)
        }
        setnames(gene.de, old=from, new=columnName)
    } else {
        gene.de[, (columnName) := eval(from, envir = gene.de)]
    }

    # :ToDo: add types in meta?
    if (columnName == "ID") {
        gene.de[, ID := as.character(ID)]
        setkey(gene.de, ID)
    }

    if (columnName %in% c("baseMean", "pval", "log2FC", "signalRank")) {
        x <- gene.de[[columnName]]
        gene.de[, (columnName) := as.numeric(x)]
    }

    if (columnName == "pval") {
        if (any(gene.de$pval == 0, na.rm = T)) {
            mpval <- min(gene.de$pval[gene.de$pval != 0])
            gene.de[pval == 0, pval := mpval]
            warning("Some of your p-values are equal to zero. Replaced it by the minimum non-zero p-value.")
        }
    }
}

#' Makes data.table with differential expression results containing
#' all columns required for gatom and in the expected format
#' based on metadata object
#' @param de.raw Table with defferential expression results, an object
#'        convertable to data.frame
#' @param de.meta Object with differential expression table metadata
#'        acquired with getGeneDEMeta or getMetDEMeta functions
#' @return data.table object with converted differential expressiond table
#' @export
#' @importFrom methods is
#' @examples
#' data("org.Mm.eg.gatom.annoEx")
#' data("gene.de.rawEx")
#' de.meta <- getGeneDEMeta(gene.de.rawEx, org.gatom.anno = org.Mm.eg.gatom.annoEx)
#' de <- prepareDE(gene.de.rawEx, de.meta)
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
                         remove.gene.versions=TRUE,
                         split=" */// *") {
    # first looking for column with base IDs
    de.sample <- if (nrow(de) < sample.size) {
        copy(de)
    } else {
        de[sample(seq_len(nrow(de)), sample.size), ]
    }
    columnSamples <- lapply(de.sample, as.character)

    if (length(split) != 0) {
        columnSamples <- lapply(columnSamples, function(ids) unlist(splitIds(ids, split=split)))
    }

    if (remove.gene.versions) {
        columnSamples <- lapply(columnSamples, gsub,
                                pattern="\\.\\d*$",
                                replacement="")
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



    # Check whether there is column called ID and it matches something
    if ("id" %in% tolower(rownames(z))) {
        z1 <- z[tolower(rownames(z)) == "id", , drop=FALSE]
        if (max(z1) / nrow(de.sample) >= match.threshold) {
            bestMatch <- which(z1 == max(z1), arr.ind = TRUE)[1,]
            return(list(column=rownames(z1)[bestMatch["row"]],
                        type=colnames(z1)[bestMatch["col"]]))
        }
    }

    # No column named ID, but may be something with unique values

    idCandidates <- colnames(de)[sapply(de,
                                        pryr::compose(all, `!`, duplicated))]

    z1 <- z[idCandidates, , drop=FALSE]
    if (max(z1) / nrow(de.sample) >= match.threshold) {
        bestMatch <- which(z1 == max(z1), arr.ind = TRUE)[1,]
        return(list(column=rownames(z1)[bestMatch["row"]],
                    type=colnames(z1)[bestMatch["col"]]))
    }

    # Whatever else matches best

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
#'
#' Default values for all columns are NULL which mean they are
#' determined automatically.
#'
#' @param gene.de.raw A table with differential expression results, an object
#'                    convertable to data.frame.
#' @param org.gatom.anno Organsim-specific annotation obtained from
#'                       makeOrgGatomAnnotation function.
#' @param idColumn Specifies column name with gene identifiers.
#' @param idType Specifies type of gene IDs (one of the supported by annotation).
#' @param pvalColumn Specifies column with p-values.
#' @param logPvalColumn Specifies column with log p-values, if there is no such
#'                      column one will be generated automatically.
#' @param log2FCColumn Specifies column with log2-fold changes.
#' @param baseMeanColumn Specifies column with average expression across samples.
#' @param signalColumn Specifies column with identifier of the measured entity
#'      (such as gene ID for RNA-seq and probe ID for microarrays).
#'      Could be NULL (automatic, set from based on pval and log2FC columns),
#'      character (column name), or function (evaluated in a scope of original data frame)
#' @param signalRankColumn Specifies how the genes are ranked from highly to lowly expressed,
#'      used in `addHighlyExpressedEdgues` function.
#'      Could be NULL (automatic), character (column name)
#'      function (evaluated in a scope of original data frame).
#' @export
#' @examples
#' data("org.Mm.eg.gatom.annoEx")
#' data("gene.de.rawEx")
#' de.meta <- getGeneDEMeta(gene.de.rawEx, org.gatom.anno = org.Mm.eg.gatom.annoEx)
#' de <- prepareDE(gene.de.rawEx, de.meta)
getGeneDEMeta <- function(gene.de.raw,
                          org.gatom.anno,
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
                                 c("pval", "p.value", "pvalue", "p-value"))
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
            signalColumn <- quote(paste0("gs", as.integer(as.factor(paste0(pval, "_", log2FC)))))
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
                         signalColumn=NULL
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
                                 c("pval", "p.value", "pvalue", "p-value"))
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
                                     "logfc", "log2(fc)", "log(fc)"))
    }

    if (is.null(baseMeanColumn)) {
        baseMeanColumn <- findColumn(met.de.raw,
                                     c("baseMean", "aveexpr"))
    }

    if (is.null(signalColumn)) {
        signalColumn <- findColumn(met.de.raw,
                                  c("signal", "ionIdx", "ion"))
        # :TODO: add check for (pval, log2FC) uniqueness for found signal column
        if (is.na(signalColumn)) {
            signalColumn <- quote(paste0("ms", as.integer(as.factor(paste0(pval, "_", log2FC)))))
        }
    }

    list(idType=idType,
         columns=list(
             ID=idColumn,
             pval=pvalColumn,
             logPval=logPvalColumn,
             log2FC=log2FCColumn,
             baseMean=baseMeanColumn,
             signal=signalColumn))
}
