#' @import data.table
.makeEdgeTable <- function(network, org.gatom.anno, gene.de, gene.de.meta,
                           gene2reaction.extra=NULL,
                           keepReactionsWithoutEnzymes=FALSE) {
    if (is.null(gene.de)) {
        gene.pvals <- data.table(gene=org.gatom.anno$genes$gene,
                                 pval=NA)
        setkey(gene.pvals, gene)
    } else {
        gene.pvals <- gene.de[, list(ID=ID, pval=pval, origin=seq_len(nrow(gene.de)))]
        gene.pvals <- gene.pvals[order(pval)][!duplicated(ID)]
        setkey(gene.pvals, ID)

        if (gene.de.meta$idType != org.gatom.anno$baseId) {
            gene.pvals <- convertPvalDT(gene.pvals,
                                        org.gatom.anno$mapFrom[[gene.de.meta$idType]],
                                        removeGeneVersions = TRUE)
        } else {
            setnames(gene.pvals, "ID", "gene")
        }
    }

    gene.pvals <- org.gatom.anno$genes[gene.pvals]

    enzyme.pvals <- convertPvalDT(gene.pvals,
                                  org.gatom.anno$gene2enzyme)

    if (keepReactionsWithoutEnzymes) {
        enzyme.pvals <- rbind(enzyme.pvals, list(enzyme="-", symbol="-", gene="-"), fill=TRUE)
        setkeyv(enzyme.pvals, "enzyme")
    }

    reaction.pvals <- convertPvalDT(enzyme.pvals,
                                    network$enzyme2reaction)

    if (!is.null(gene2reaction.extra)) {
        stopifnot(c("gene", "reaction") %in% colnames(gene2reaction.extra))
        if (!identical(key(gene2reaction.extra), "gene")) {
            gene2reaction.extra <- as.data.table(gene2reaction.extra)
            setkey(gene2reaction.extra, "gene")
        }
        reaction.pvals.extra <- convertPvalDT(gene.pvals,
                                              gene2reaction.extra)
        reaction.pvals.extra[, enzyme := "-"]

        reaction.pvals <- rbind(reaction.pvals, reaction.pvals.extra)
        reaction.pvals <- reaction.pvals[order(reaction, pval, gene == "-", enzyme == "-"), ]
        reaction.pvals <- reaction.pvals[!duplicated(reaction), ]
    }

    reaction.pvals <- merge(reaction.pvals, network$reactions)

    align.pvals <- convertPvalDT(reaction.pvals, network$reaction2align)

    edge.table <- copy(align.pvals)
    if (!is.null(gene.de)) {
        edge.table[, colnames(gene.de) := gene.de[origin]]
        edge.table[, ID := NULL]
    }
    setnames(edge.table, "symbol", "label")
    setnames(edge.table, "reaction_url", "url")


    # if there are parallel edges select one with the best gene by p-value
    edge.table[atom.x > atom.y, c("atom.x", "atom.y") := list(atom.y, atom.x)]
    edge.table <- edge.table[order(pval),]
    edge.table <- edge.table[!duplicated(edge.table[, list(atom.x, atom.y)]), ]

    edge.table[]
}

#' @import data.table
.makeVertexTable <- function(network, atoms, met.db, met.de, met.de.meta) {
    if (is.null(met.de)) {
        metabolite.pvals <- data.table(metabolite=character(0), pval=numeric(0), origin=integer(0))
        setkeyv(metabolite.pvals, "metabolite")
    } else {
        metabolite.pvals <- met.de[, list(ID=ID, pval=pval, origin=seq_len(nrow(met.de)))]
        if (met.de.meta$idType != met.db$baseId) {
            metabolite.pvals <- convertPvalDT(metabolite.pvals,
                                              met.db$mapFrom[[met.de.meta$idType]])
        } else {
            setnames(metabolite.pvals, "ID", "metabolite")
        }
    }



    # Extending p-values to anomers
    base_metabolite.pvals <- convertPvalDT(metabolite.pvals, met.db$anomers$metabolite2base_metabolite)
    base_metabolite.pvals[, metabolite := NULL]
    inferred_metabolite.pvals <- convertPvalDT(base_metabolite.pvals, met.db$anomers$base_metabolite2metabolite)
    inferred_metabolite.pvals[, base_metabolite := NULL]

    metabolite.pvals <- rbind(metabolite.pvals, inferred_metabolite.pvals)
    metabolite.pvals <- metabolite.pvals[!duplicated(metabolite)]
    setkey(metabolite.pvals, metabolite)

    vertex.table <- data.table(atom=atoms)
    setkey(vertex.table, atom)
    vertex.table <- network$atoms[vertex.table]
    vertex.table <- merge(vertex.table,
                          met.db$metabolites,
                          by="metabolite",
                          all.x=TRUE)
    vertex.table <- merge(vertex.table, metabolite.pvals, all.x=T)
    if (!is.null(met.de)) {
        vertex.table[, colnames(met.de) := met.de[origin]]
        vertex.table[, ID := NULL]
    }

    setcolorder(vertex.table, c("atom", setdiff(colnames(vertex.table), "atom")))
    setnames(vertex.table, "metabolite_name", "label")
    setnames(vertex.table, "metabolite_url", "url")
    vertex.table[]
}

#' Creates metabolic graph based on specified data
#' @param network Network object
#' @param topology Way to determine network vertices
#' @param org.gatom.anno Organism annotation object
#' @param gene.de Table with the differential gene expression, set to NULL if absent
#' @param gene.de.meta Annotation of `gene.de` table
#' @param gene.keep.top Only the `gene.keep.top` of the most expressed genes will be kept for the network
#' @param met.db Metabolite database
#' @param met.de Table with the differential expression for metabolites, set to NULL if absent
#' @param met.de.meta Annotation of `met.de` table
#' @param met.to.filter List of metabolites to filter from the network
#' @param largest.component If TRUE, only the largest connected component is returned
#' @param gene2reaction.extra Additional gene to reaction mappings. Should be a data.table with `gene` and `reaction` columns
#' @param keepReactionsWithoutEnzymes If TRUE, keep reactions that have no annotated enzymes, thus expanding the network but
#'                                    including some reactions which are not possible in the considered species.
#' @export
#' @import igraph
makeMetabolicGraph <- function(network,
                          topology=c("atoms", "metabolites"),
                          org.gatom.anno,
                          gene.de,
                          gene.de.meta=getGeneDEMeta(gene.de, org.gatom.anno),
                          gene.keep.top=12000,
                          met.db,
                          met.de,
                          met.de.meta=getMetDEMeta(met.de, met.db),
                          met.to.filter=fread(system.file("mets2mask.lst", package="gatom"))$ID,
                          gene2reaction.extra=NULL,
                          keepReactionsWithoutEnzymes=FALSE,
                          largest.component=TRUE) {

    topology <- match.arg(topology)

    if (!is.null(gene.de)) {
        .messagef("Found DE table for genes with %s IDs", gene.de.meta$idType)
        gene.de <- prepareDE(gene.de, gene.de.meta)
        gene.de <- gene.de[signalRank <= gene.keep.top]
    }

    if (!is.null(met.de)) {
        .messagef("Found DE table for metabolites with %s IDs", met.de.meta$idType)
    }

    met.de <- prepareDE(met.de, met.de.meta)


    edge.table <- .makeEdgeTable(network=network,
                                 org.gatom.anno=org.gatom.anno,
                                 gene.de=gene.de,
                                 gene.de.meta=gene.de.meta,
                                 gene2reaction.extra=gene2reaction.extra,
                                 keepReactionsWithoutEnzymes=keepReactionsWithoutEnzymes)
    all.atoms <- union(edge.table$atom.x, edge.table$atom.y)

    if (topology == "metabolites"){

        # This construction is used to keep metaboliteas as two FIRST columns
        edge.table <- cbind(network$atoms$metabolite[match(edge.table$atom.x, network$atoms$atom)],
                            network$atoms$metabolite[match(edge.table$atom.y, network$atoms$atom)],
                            edge.table)
        colnames(edge.table)[c(1,2)] <- c("metabolite.x", "metabolite.y")
        edge.table[, atom.x := NULL]
        edge.table[, atom.y := NULL]

        # removing parallel edges by selecting only the best gene by p-value
        stopifnot(!is.unsorted(edge.table$pval, na.rm = TRUE))
        edge.table <- edge.table[!duplicated(edge.table[, list(metabolite.x, metabolite.y)]), ]
    }

    vertex.table <- .makeVertexTable(network=network,
                                     atoms=all.atoms,
                                     met.db=met.db,
                                     met.de=met.de,
                                     met.de.meta=met.de.meta)

    if (topology == "metabolites"){
        # During conversion to graph, first column will be renamed to "name";
        # we will still need the column called "metabolite" for further functions
        vertex.table[, 1] <- vertex.table$metabolite
        vertex.table <- vertex.table[!duplicated(vertex.table), ]
    }

    if ("name" %in% colnames(vertex.table)) {
        setnames(vertex.table, "name", "name.orig")
    }

    g <- graph.data.frame(edge.table, directed=FALSE, vertices = vertex.table)

    if (!is.null(met.to.filter)) {
        nodes.to.del <- V(g)[metabolite %in% met.to.filter]
        if (length(nodes.to.del) == 0) {
            warning("Found no metabolites to mask")
        } else {
            g <- delete_vertices(g, v = V(g)[metabolite %in% met.to.filter])
        }
    }

    if (largest.component) {
        gc <- components(g)
        g <- induced.subgraph(g, gc$membership == which.max(gc$csize))
    }

    g
}

.reversefdrThreshold <- function(pt, fb){
    pihat <- fb$lambda + (1 - fb$lambda) * fb$a
    (pihat * pt) / (-fb$lambda * pt^fb$a + pt^fb$a + fb$lambda * pt)
}

#' Score metabolic graph
#'
#' @param g Metabolic graph obtained with makeMetabolic graph function
#' @param k.gene Number of gene signals to be scored positively, the higher is the number,
#'     the larger will be the resulting module. If set to NULL, genes will not be used for scoring.
#' @param k.met Number of metabolite signals to be scored positively, the higher is the number,
#'     the larger will be the resulting module. If set to NULL, metabolites will not be used for scoring.
#' @import BioNet
#' @importFrom mwcsr normalize_sgmwcs_instance
#' @export
scoreGraph <- function(g, k.gene, k.met,
                       vertex.threshold.min=0.1,
                       edge.threshold.min=0.1,
                       met.score.coef=1,
                       show.warnings=TRUE,
                       raw=FALSE) {
    if (show.warnings) {
        warnWrapper <- identity
    } else {
        warnWrapper <- suppressWarnings
    }

    vertex.table <- data.table(as_data_frame(g, what="vertices"))
    edge.table <- data.table(as_data_frame(g, what="edges"))
    if (!is.null(k.met)) {
        pvalsToFit <- vertex.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]

        warnWrapper(vertex.bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
        if (vertex.bum$a > 0.5) {
            V(g)$score <- 0
            warning("Vertex scores have been assigned to 0 due to an inappropriate p-value distribution")

        } else {
            vertex.threshold <- if (k.met > length(pvalsToFit)) 1 else {
            sort(pvalsToFit)[k.met]
                }

            vertex.threshold <- min(vertex.threshold,
                                BioNet::fdrThreshold(vertex.threshold.min, vertex.bum))

            met.fdr <- .reversefdrThreshold(vertex.threshold, vertex.bum)

            .messagef("Metabolite p-value threshold: %f", vertex.threshold)
            .messagef("Metabolite BU alpha: %f", vertex.bum$a)
            .messagef("FDR for metabolites: %f", met.fdr)

            V(g)$score <- with(vertex.table,
                               (vertex.bum$a - 1) *
                               (log(.replaceNA(pval, 1)) - log(vertex.threshold)))
            V(g)$score <- V(g)$score * met.score.coef
            }
        }
    else {
        V(g)$score <- 0
        V(g)$signal <- ""

    }

    if (!is.null(k.gene)) {
        pvalsToFit <- edge.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]

        warnWrapper(edge.bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F))
        if(edge.bum$a > 0.5) {
            E(g)$score <- 0
            warning("Edge scores have been assigned to 0 due to an inappropriate p-value distribution")

        } else {
            edge.threshold <- if (k.gene > length(pvalsToFit)) 1 else {
            sort(pvalsToFit)[k.gene]
                }

            edge.threshold <- min(edge.threshold,
                              BioNet::fdrThreshold(edge.threshold.min, edge.bum))

            gene.fdr <- .reversefdrThreshold(edge.threshold, edge.bum)

            .messagef("Gene p-value threshold: %f", edge.threshold)
            .messagef("Gene BU alpha: %f", edge.bum$a)
            .messagef("FDR for genes: %f", gene.fdr)

            E(g)$score <- with(edge.table,
                               (edge.bum$a - 1) *
                               (log(.replaceNA(pval, 1)) - log(edge.threshold)))
        }
    }
    else {
        E(g)$score <- 0
        E(g)$signal <- ""
    }
    g
    if (raw) {
        return(g)
    }

    res <- normalize_sgmwcs_instance(g,
                                     nodes.weight.column = "score",
                                     edges.weight.column = "score",
                                     nodes.group.by = "signal",
                                     edges.group.by = "signal",
                                     group.only.positive = TRUE)
    res
}

#' Connect atoms belonging to the same metabolite with edges
#' @param m Metabolic module
#' @export
connectAtomsInsideMetabolite <- function(m) {
    t <- data.frame(v=V(m)$name, met=V(m)$metabolite, stringsAsFactors=F)
    toCollapse <- merge(t, t, by="met")
    toCollapse <- toCollapse[(toCollapse$v.x < toCollapse$v.y), ]

    z <- matrix(c(toCollapse$v.x, toCollapse$v.y),
                nrow=2,
                byrow=T)
    res <- igraph::add.edges(m, t(as.matrix(toCollapse[, c("v.x", "v.y")])))
    res
}

#' Collapse atoms belonging to the same metabolite into one vertex
#' @param m Metabolic module
#' @export
collapseAtomsIntoMetabolites <- function(m) {
    vertex.table <- data.table(as_data_frame(m, what="vertices"))
    edge.table <- data.table(as_data_frame(m, what="edges"))

    atom2metabolite <- vertex.table[, setNames(metabolite, name)]

    vertex.table.c <- copy(vertex.table)
    vertex.table.c[, name := atom2metabolite[name]]
    vertex.table.c <- vertex.table.c[!duplicated(name)]

    edge.table.c <- copy(edge.table)
    edge.table.c[, from :=  atom2metabolite[from]]
    edge.table.c[, to :=  atom2metabolite[to]]
    edge.table.c[from > to, c("from", "to") := list(to, from)]
    edge.table.c <- edge.table.c[!duplicated(edge.table.c[, list(from, to)])]


    res <- graph.data.frame(edge.table.c, directed=FALSE, vertices=vertex.table.c)
    res
}

#' Add reactions without highly changing genes but with high average expression
#' @param m Metabolic module
#' @param g Scored graph
#' @param top Maximum rank value for the gene to be considered highly expressed
#' @import igraph
#' @import plyr
#' @export
addHighlyExpressedEdges <- function(m, g, top=3000) {
    if (!"signalRank" %in% list.edge.attributes(g)) {
        warning("No signalRank edge attribute, returing graph as is")
        return(m)
    }

    vertex.table <- as_data_frame(m, what=c("vertices"))
    edge.table <- as_data_frame(m, what=c("edges"))



    toAdd.edge.table <-
        as_data_frame(g, what=c("edges")) %>%
        subset(from %in% vertex.table$name) %>%
        subset(to %in% vertex.table$name) %>%
        subset(signalRank <= top)

    res.edge.table <-
        rbind.fill(edge.table, toAdd.edge.table) %>%
        subset(!duplicated(paste0(from, "_", to)))

    res <- graph.data.frame(res.edge.table, directed=FALSE, vertices=vertex.table)
}
