#' Create an organism annotation object for network analysis
#'
#' @param org.db Bioconductor org.db object, e.g. org.Mm.eg.db
#' @param idColumns vector of column names from `org.db` object to cread ID mappings.
#'                  First ID will be used as a base identifier, should be compatible
#'                  with KEGG and Reactome databses.
#' @param nameColumn column with a human readable gene symbol. Default to "SYMBOL".
#' @param enzymeColumn column with an Enzyme Commission ID. Defatult to "ENZYME".
#' @param appendEnzymesFromKegg if TRUE, KEGG databases will be sued to extend gene to
#'     enzyme mappings obtained from org.db package.
#' @param keggOrgCode KEGG organizm code, e.g. "mmu". If set to NULL, the code is determined
#'     automatically.
#'
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
                                   keggOrgCode=NULL) {

    stopifnot(requireNamespace("AnnotationDbi"))
    if (appendEnzymesFromKegg) {
        stopifnot(requireNamespace("KEGGREST"))
        if (is.null(keggOrgCode)) {
            meta <- AnnotationDbi::metadata(org.db)
            organismName <- meta$value[match("ORGANISM", meta$name)]
            keggOrgCodes <- data.table(KEGGREST::keggList("organism"))
            keggOrgCode <- keggOrgCodes[grep(organismName, species), organism]
        }
    }

    baseColumn <- idColumns[1]

    org.gatom.anno <- list()
    org.gatom.anno$genes <- data.table(gene=keys(org.db, keytype = baseColumn))
    setkey(org.gatom.anno$genes, gene)
    org.gatom.anno$genes[, symbol := AnnotationDbi::mapIds(org.db, keys=gene,
                                            keytype = baseColumn,
                                            column = nameColumn)]
    org.gatom.anno$baseId <- names(baseColumn)
    org.gatom.anno$gene2enzyme <- data.table(
        AnnotationDbi::select(org.db,
                              keys=org.gatom.anno$genes$gene,
                              columns = c("ENZYME"),
                              keytype=baseColumn))
    setnames(org.gatom.anno$gene2enzyme,
             c(baseColumn, enzymeColumn),
             c("gene", "enzyme"))

    # sometime mapping is not direct (e.g. for Sc.sgd), keeping only gene & enzyme
    org.gatom.anno$gene2enzyme <- org.gatom.anno$gene2enzyme[, list(gene, enzyme)]

    org.gatom.anno$gene2enzyme <- org.gatom.anno$gene2enzyme[!is.na(enzyme)]

    if (appendEnzymesFromKegg) {
        kegg.gene2enzyme <- KEGGREST::keggLink("enzyme", keggOrgCode)
        kegg.gene2enzyme <- data.table(gene=gsub(paste0(keggOrgCode, ":"), "",
                                                 names(kegg.gene2enzyme)),
                                       enzyme=gsub("ec:", "", kegg.gene2enzyme))

        org.gatom.anno$gene2enzyme <- rbind(org.gatom.anno$gene2enzyme,
                                            kegg.gene2enzyme)
        org.gatom.anno$gene2enzyme <- unique(org.gatom.anno$gene2enzyme)
        setkey(org.gatom.anno$gene2enzyme, gene)
    }


    org.gatom.anno$mapFrom <- list()
    for (i in tail(seq_along(idColumns), -1)) {
        otherColumn <- idColumns[i]
        n <- names(otherColumn)
        org.gatom.anno$mapFrom[[n]] <- data.table(
            AnnotationDbi::select(org.db,
                                  keys=org.gatom.anno$genes$gene,
                                  columns = otherColumn))
        org.gatom.anno$mapFrom[[n]] <- na.omit(org.gatom.anno$mapFrom[[n]])
        setnames(org.gatom.anno$mapFrom[[n]],
                 c(baseColumn, otherColumn),
                 c("gene", names(otherColumn)))
        setkeyv(org.gatom.anno$mapFrom[[n]], n)
    }

    # todo: support reactome pathways for non-entrez base IDs
    pathways <- getMetabolicPathways(org.gatom.anno$genes$gene, keggOrgCode,
                                     includeReactome=(org.gatom.anno$baseId == "ENTREZID"))

    # keeping only pathways with at least half of enzymatic genes
    pathwaysEnzymeRatio <-
      lengths(lapply(pathways, intersect, y=unique(org.gatom.anno$gene2enzyme$gene)))/
      lengths(pathways)
    pathways <- pathways[pathwaysEnzymeRatio >= 0.5]

    org.gatom.anno$pathways <- pathways


    org.gatom.anno
}

#' Generate list of metabolic pathways from Reactome and KEGG databases
#'
#' @param universe list of gene
#' @param keggOrgCode KEGG organism code, like mmu or hsa
#' @param includeReactome whether to include Reactome pathways (only works for Entrez ID universe)
#' @param includeKEGG wheter to include KEGG pahtways and modules
getMetabolicPathways <- function(universe,
                                 keggOrgCode,
                                 includeReactome=TRUE, includeKEGG=TRUE){

  if (includeReactome) {
    reactomepath <- na.omit(AnnotationDbi::select(reactome.db::reactome.db, universe, "PATHID", "ENTREZID"))
    reactomepath <- split(reactomepath$ENTREZID, reactomepath$PATHID)

    reactomepathway2name <- data.table::as.data.table(na.omit(
      AnnotationDbi::select(reactome.db::reactome.db,
                            names(reactomepath),
                            c("PATHNAME"), 'PATHID')))
    reactomepathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)]

  } else {
    reactomepath <- NULL
    reactomepathway2name <- data.table(PATHNAME=character(0), PATHID=character(0))
  }

  if (includeKEGG) {
    keggmodule <- KEGGREST::keggLink(keggOrgCode, "module")
    keggmodule <- gsub(paste0(keggOrgCode, ":"), "", keggmodule)
    names(keggmodule) <- gsub("md:", "", names(keggmodule))
    keggmodule <- split(keggmodule, names(keggmodule))

    # keggmdnames <- KEGGREST::keggList("module", keggOrgCode) # 404 after September, 2019
    keggmdnames <- KEGGREST::keggList("module")
    keggmd2name <- data.table::as.data.table(keggmdnames, keep.rownames=T)
    keggmd2name$rn <- gsub("md:", "", keggmd2name$rn)
    data.table::setnames(keggmd2name, c("rn","keggmdnames"), c("PATHID","PATHNAME"))
    keggmd2name$PATHID <- paste0(keggOrgCode, "_", keggmd2name$PATHID)

    keggpathway <- KEGGREST::keggLink(keggOrgCode, "pathway")
    keggpathway <- gsub(paste0(keggOrgCode, ":"), "", keggpathway)
    names(keggpathway) <- gsub("path:", "", names(keggpathway))
    keggpathway <- split(keggpathway, names(keggpathway))
    keggpathway <- lapply(keggpathway, unname)
    keggpathnames <- KEGGREST::keggList("pathway", keggOrgCode)

    keggpath2name <- data.table::as.data.table(keggpathnames, keep.rownames=T)
    keggpath2name$rn <- gsub("path:", "", keggpath2name$rn)
    keggpath2name$keggpathnames <- gsub(" - [^-]*$", "", keggpath2name$keggpathnames)
    data.table::setnames(keggpath2name, c("rn","keggpathnames"), c("PATHID","PATHNAME"))
  } else {
    keggmodule <- NULL
    keggmdy2name <- data.table(PATHNAME=character(0), PATHID=character(0))
    keggpathway <- NULL
    keggpath2name <- data.table(PATHNAME=character(0), PATHID=character(0))
  }



  pathways <- c(reactomepath, keggmodule, keggpathway)
  pathways <- lapply(pathways, intersect, y=universe)
  pathways <- lapply(pathways, unique)
  pathways <- pathways[lengths(pathways) >= 1]
  pathway2name <- do.call("rbind", list(reactomepathway2name,
                                        keggmd2name,
                                        keggpath2name))
  pthws2exclude <- c("Metabolism",
                     "Metabolic pathways",
                     "Carbon metabolism",
                     "Fatty acid metabolism",
                     "Metabolism of lipids",
                     "Metabolism of carbohydrates",
                     "2-Oxocarboxylic acid metabolism",
                     "Biosynthesis of amino acids")
  ids2exclude <- which(names(pathways) %in%
                           pathway2name$PATHID[match(pthws2exclude, pathway2name$PATHNAME)])
  # adding non-metabolic KEGG pathways to ids2exclude:
  ids2exclude <- c(ids2exclude, grep("(mmu|hsa)0[2-9]", names(pathways)))
  pathways <- pathways[-ids2exclude]

  names(pathways) <- sapply(
      seq_along(names(pathways)),
      function(i) paste(names(pathways)[[i]],
                        pathway2name$PATHNAME[match(names(pathways)[[i]], pathway2name$PATHID)],
                        sep = ": "))
  pathways
}
