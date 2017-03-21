#' @import KEGGREST
#' @export
makeOrgGatomAnnotation <- function(org.db,
                                   idColumns=
                                       c("Entrez"="ENTREZID",
                                         "RefSeq"="REFSEQ",
                                         "Ensembl"="ENSEMBL",
                                         "Symbol"="SYMBOL"),
                                   nameColumn="SYMBOL",
                                   enzymeColumn="ENZYME",
                                   appendEnzymesFromKegg=TRUE,
                                   keggOrgCode=NULL
) {
    if (appendEnzymesFromKegg) {
        stopifnot(require(KEGGREST))
        if (is.null(keggOrgCode)) {
            meta <- AnnotationDbi::metadata(org.db)
            organismName <- meta$value[match("ORGANISM", meta$name)]
            keggOrgCodes <- data.table(keggList("organism"))
            keggOrgCode <- keggOrgCodes[grep(organismName, species), organism]
        }
    }

    baseColumn <- idColumns[1]

    org.gatom.anno <- list()
    org.gatom.anno$genes <- data.table(gene=keys(org.db, keytype = baseColumn))
    setkey(org.gatom.anno$genes, gene)
    org.gatom.anno$genes[, symbol := mapIds(org.db, keys=gene,
                                            keytype = baseColumn,
                                            column = nameColumn)]
    org.gatom.anno$baseId <- names(baseColumn)
    org.gatom.anno$gene2enzyme <- data.table(select(org.db, keys=org.gatom.anno$genes$gene, columns = c("ENZYME")))
    setnames(org.gatom.anno$gene2enzyme, c(baseColumn, enzymeColumn), c("gene", "enzyme"))
    org.gatom.anno$gene2enzyme <- org.gatom.anno$gene2enzyme[!is.na(enzyme)]

    if (appendEnzymesFromKegg) {
        kegg.gene2enzyme <- keggLink("enzyme", keggOrgCode)
        kegg.gene2enzyme <- data.table(gene=gsub(paste0(keggOrgCode, ":"), "",
                                                 names(kegg.gene2enzyme)),
                                       enzyme=gsub("ec:", "", kegg.gene2enzyme))

        org.gatom.anno$gene2enzyme <- rbind(org.gatom.anno$gene2enzyme, kegg.gene2enzyme)
        org.gatom.anno$gene2enzyme <- unique(org.gatom.anno$gene2enzyme)
        setkey(org.gatom.anno$gene2enzyme, gene)
    }


    org.gatom.anno$mapFrom <- list()
    for (i in tail(seq_along(idColumns), -1)) {
        otherColumn <- idColumns[i]
        n <- names(otherColumn)
        org.gatom.anno$mapFrom[[n]] <- data.table(select(org.db,
                                                         keys=org.gatom.anno$genes$gene,
                                                         columns = otherColumn))
        setnames(org.gatom.anno$mapFrom[[n]],
                 c(baseColumn, otherColumn),
                 c("gene", names(otherColumn)))
        setkeyv(org.gatom.anno$mapFrom[[n]], n)
    }

    org.gatom.anno
}
