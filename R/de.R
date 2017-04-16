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


prepareGeneDE <- function(gene.de, gene.de.meta)
{
    if (!is(de, "data.table")) {
        de <- as.data.table(as.data.frame(de), keep.rownames = !is.numeric(attr(de,
                                                                                "row.names")))
    }
    rename.smart(de, ID = c("gene", "entrez", "rn", "symbol"),
                 pval = c("p.value", "pvalue"), log2FC = c("log2foldchange",
                                                           "logfc"))
    de <- de[!duplicated(de$ID), ]
    de <- de[!is.na(de$pval), ]
    de <- as.data.table(de[order(de$pval), ])
    de
}

findColumn <- function(de, names) {
     candidates <-  na.omit(match(tolower(names),
                                  tolower(colnames(de))))
     if (length(candidates) == 0) {
         return(NA)
     }
     return(colnames(de)[candidates[1]])
}

findIdColumn <- function(de, ids.list,
                         sample.size=1000,
                         match.threshold=0.6,
                         remove.ensembl.revisions=TRUE) {
    # first looking for column with base IDs
    gene.de.raw.sample <- if (nrow(gene.de.raw) < sample.size) {
        copy(gene.de.raw)
    } else {
        gene.de.raw[sample(seq_len(nrow(gene.de.raw)), sample.size)]
    }
    columnSamples <- lapply(gene.de.raw.sample, as.character)


    if (remove.ensembl.revisions) {
        columnSamples <- lapply(columnSamples, gsub,
                                pattern="(ENS\\w*\\d*)\\.\\d*",
                                replacement="\\1")
    }

    ss <- sapply(columnSamples,
                 pryr::compose(length, intersect), ids.list[[1]])

    if (max(ss) / nrow(gene.de.raw.sample) >= match.threshold) {
        # we found a good column with base IDs
        return(list(idColumn=colnames(de)[which.max(ss)],
                    idType=names(ids.list)[1]))
    }

    z <- pairwiseCompare(pryr::compose(length, intersect),
                         columnSamples,
                         ids.list)

    bestMatch <- which(z == max(z), arr.ind = TRUE)[1,]
    return(list(idColumn=colnames(de)[bestMatch["row"]],
                idType=names(ids.list)[bestMatch["col"]]))
}

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
                          probeColumn=NULL
                          ) {

    if (is.null(idColumn) != is.null(idType)) {
        stop("Either both or none of idColumn and idType can be specified")
    }

    if (is.null(idColumn)) {
        ids.list <- c(list(org.gatom.anno$genes$gene),
                      lapply(names(org.gatom.anno$mapFrom),
                             function(n) org.gatom.anno$mapFrom[[n]][[n]])
                      )
        names(ids.list) <- c(org.gatom.anno$baseId, names(org.gatom.anno$mapFrom))
        idColumnInfo <- findIdColumn(gene.de.raw, ids.list)
        idColumn <- idColumnInfo$idColumn
        idType <- idColumnInfo$idType
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


    list(idColumn=idColumn,
         idType=idType,
         pvalColumn=pvalColumn,
         logPvalColumn=logPvalColumn,
         baseMeanColumn=baseMeanColumn,
         probeColumn=probeColumn)
}
