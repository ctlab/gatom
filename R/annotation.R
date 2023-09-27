#' Create an organism annotation object for network analysis
#'
#' @param org.db Bioconductor org.db object, e.g. org.Mm.eg.db
#' @param idColumns vector of column names from `org.db` object to creat ID mappings.
#'                  First ID will be used as a base identifier, should be compatible
#'                  with KEGG and Reactome databases.
#' @param nameColumn column with a human readable gene symbol. Default to "SYMBOL".
#' @param enzymeColumn column with an Enzyme Commission ID. Default to "ENZYME".
#' @param appendEnzymesFromKegg if TRUE, KEGG databases will be sued to extend gene to
#'     enzyme mappings obtained from org.db package.
#' @param appendOrthologiesFromKegg if TRUE, KEGG database will be sued to extend gene to
#'     orthology mappings obtained from org.db package
#' @param filterNonSpecificEnzymes if TRUE, will filter out non-specific enzymes from
#'     gene to enzyme mappings obtained from org.db package
#' @param keggOrgCode KEGG organism code, e.g. "mmu". If set to NULL, the code is determined
#'     automatically.
#'
#' @return organism annotation object that will be used for network analysis
#'
#' @examples
#' library(org.Mm.eg.db)
#' org.Mm.eg.gatom.anno <- makeOrgGatomAnnotation(org.db = org.Mm.eg.db)
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
                                   appendOrthologiesFromKegg=TRUE,
                                   filterNonSpecificEnzymes=TRUE,
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
    org.gatom.anno$genes[is.na(symbol), symbol := gene]
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

    if (appendOrthologiesFromKegg) {
        kegg.gene2orthology <- KEGGREST::keggLink("orthology", keggOrgCode)
        kegg.gene2orthology <- data.table(gene=gsub(paste0(keggOrgCode, ":"), "",
                                                    names(kegg.gene2orthology)),
                                          enzyme=gsub("ko:", "", kegg.gene2orthology))

        org.gatom.anno$gene2enzyme <- rbind(org.gatom.anno$gene2enzyme,
                                            kegg.gene2orthology)
        org.gatom.anno$gene2enzyme <- unique(org.gatom.anno$gene2enzyme)
        setkey(org.gatom.anno$gene2enzyme, gene)
    }

    if (filterNonSpecificEnzymes) {
        org.gatom.anno$gene2enzyme <- org.gatom.anno$gene2enzyme[!(like(enzyme, "-"))]
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
    org.gatom.anno$pathways <- getMetabolicPathways(universe=org.gatom.anno$genes$gene,
                                                    metGenes=unique(org.gatom.anno$gene2enzyme$gene),
                                                    keggOrgCode=keggOrgCode,
                                                    includeReactome=(org.gatom.anno$baseId == "Entrez"))

    org.gatom.anno
}

#' Generate list of metabolic pathways from Reactome and KEGG databases
#'
#' @param universe list of genes
#' @param metGenes list of metabolic genes
#' @param keggOrgCode KEGG organism code, like mmu or hsa
#' @param threshold threshold for Fisher test to filter out non-metabolic pathways
#' @param includeReactome whether to include Reactome pathways (only works for Entrez ID universe)
#' @param includeKEGG whether to include KEGG pathways and modules
#' @return list of metabolic pathways for given organism and list of genes
getMetabolicPathways <- function(universe,
                                 metGenes,
                                 keggOrgCode,
                                 threshold=0.01,
                                 includeReactome=TRUE, includeKEGG=TRUE){


  if (includeReactome) {
    reactomepath <- na.omit(AnnotationDbi::select(reactome.db::reactome.db, universe, "PATHID", "ENTREZID"))
    reactomepath <- split(reactomepath$ENTREZID, reactomepath$PATHID)

    reactomepathway2name <- data.table::as.data.table(na.omit(
      AnnotationDbi::select(reactome.db::reactome.db,
                            names(reactomepath),
                            c("PATHNAME"), 'PATHID')))
    reactomepathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)]

    # keeping only metabolic pathways (by Reactome tree and metabolic genes enrichment):
    ids.reactome.relation <- read.table("https://reactome.org/download/current/ReactomePathwaysRelation.txt")
    colnames(ids.reactome.relation) <- c("from", "to")

    g <- igraph::graph_from_data_frame(ids.reactome.relation, directed=TRUE)
    ids <- reactomepathway2name$PATHID[which(reactomepathway2name$PATHNAME %in% c(
        "Metabolism"))]

    res <- igraph::bfs(g, root=ids, mode="out", unreachable=FALSE, dist=TRUE)
    ids.reactome.Metabolism <- names(which(res$dist > 0))

    reactomepath <- reactomepath[ids.reactome.Metabolism]
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
    keggmd2name <- data.table::as.data.table(keggmdnames, keep.rownames=TRUE)
    keggmd2name$rn <- gsub("md:", "", keggmd2name$rn)
    data.table::setnames(keggmd2name, c("rn","keggmdnames"), c("PATHID","PATHNAME"))
    keggmd2name$PATHID <- paste0(keggOrgCode, "_", keggmd2name$PATHID)

    keggpathway <- KEGGREST::keggLink(keggOrgCode, "pathway")
    keggpathway <- gsub(paste0(keggOrgCode, ":"), "", keggpathway)
    names(keggpathway) <- gsub("path:", "", names(keggpathway))
    keggpathway <- split(keggpathway, names(keggpathway))
    keggpathway <- lapply(keggpathway, unname)
    keggpathnames <- KEGGREST::keggList("pathway", keggOrgCode)

    keggpath2name <- data.table::as.data.table(keggpathnames, keep.rownames=TRUE)
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

#' Abbreviate lipid labels for lipid module
#'
#' @param module Module to prepare
#' @param orig.names whether to use original names from the dataset
#' @param abbrev.names whether to use abbreviated names for all lipids
#'
#' @return module object with abbreviated labels
#'
#' @export
abbreviateLabels <- function(module,
                             orig.names,
                             abbrev.names){
    if (is.null(module)) {
        return(NULL)
    }

    if (abbrev.names) {
        if (orig.names) {
            V(module)$tmp <- V(module)$label
            V(module)$SpecialSpeciesLabelColumn <- V(module)$Species
            V(module)$label <- V(module)$SpecialSpeciesLabelColumn
            V(module)$SpecialSpeciesLabelColumn <- V(module)$tmp
            module <- delete_vertex_attr(module, "tmp")

            V(module)$label[is.na(V(module)$label)] <- V(module)$Abbreviation_LipidMaps[is.na(V(module)$label)]
        } else {
            V(module)$SpecialSpeciesLabelColumn <- V(module)$label
            V(module)$label <- V(module)$Abbreviation_LipidMaps
        }
        V(module)$label[is.na(V(module)$label)] <- V(module)$Abbreviation_SwissLipids[is.na(V(module)$label)]
        V(module)$label[V(module)$label == "-"] <- V(module)$Abbreviation_SwissLipids[V(module)$label == "-"]
        V(module)$label[is.na(V(module)$label)] <- V(module)$SpecialSpeciesLabelColumn[is.na(V(module)$label)]
        V(module)$label[V(module)$label == "-"] <- V(module)$SpecialSpeciesLabelColumn[V(module)$label == "-"]
        module <- delete_vertex_attr(module, "SpecialSpeciesLabelColumn")
    } else if (orig.names) {
        V(module)$tmp <- V(module)$label
        V(module)$SpecialSpeciesLabelColumn <- V(module)$Species
        V(module)$label <- V(module)$SpecialSpeciesLabelColumn
        V(module)$SpecialSpeciesLabelColumn <- V(module)$tmp
        module <- delete_vertex_attr(module, "tmp")
        V(module)$label[is.na(V(module)$label)] <- V(module)$SpecialSpeciesLabelColumn[is.na(V(module)$label)]
        V(module)$label[V(module)$label == "-"] <- V(module)$SpecialSpeciesLabelColumn[V(module)$label == "-"]
        module <- delete_vertex_attr(module, "SpecialSpeciesLabelColumn")
    }

    module
}
