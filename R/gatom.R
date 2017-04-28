#' @import data.table
.makeEdgeTable <- function(network, org.gatom.anno, gene.de, gene.de.meta) {
    gene.pvals <- gene.de[, list(ID=ID, pval=pval, origin=seq_len(nrow(gene.de)))]
    gene.pvals <- gene.pvals[order(pval)][!duplicated(ID)]
    setkey(gene.pvals, ID)

    if (gene.de.meta$idType != org.gatom.anno$baseId) {
        gene.pvals <- convertPvalDT(gene.pvals,
                                    org.gatom.anno$mapFrom[[gene.de.meta$idType]])
    } else {
        setnames(gene.pvals, "ID", "gene")
    }

    gene.pvals <- org.gatom.anno$genes[gene.pvals]

    enzyme.pvals <- convertPvalDT(gene.pvals,
                                  org.gatom.anno$gene2enzyme)
    reaction.pvals <- convertPvalDT(enzyme.pvals,
                                    network$enzyme2reaction)
    reaction.pvals <- merge(reaction.pvals, network$reactions)
    rpair.pvals <- convertPvalDT(reaction.pvals,
                                 network$reaction2rpair)

    align.pvals <- convertPvalDT(rpair.pvals, network$rpair2align)

    edge.table <- copy(align.pvals)
    edge.table[, colnames(gene.de) := gene.de[origin]]
    edge.table[, ID := NULL]
    setnames(edge.table, "symbol", "label")
    setnames(edge.table, "reaction_url", "url")


    edge.table[]
}

#' @import data.table
.makeVertexTable <- function(network, atoms, met.db, met.de, met.de.meta) {
    metabolite.pvals <- met.de[, list(ID=ID, pval=pval, origin=seq_len(nrow(met.de)))]
    if (met.de.meta$idType != met.db$baseId) {
        metabolite.pvals <- convertPvalDT(metabolite.pvals,
                                          met.db$mapFrom[[met.de.meta$idType]])
    } else {
        setnames(metabolite.pvals, "ID", "metabolite")
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
    vertex.table[, colnames(met.de) := met.de[origin]]
    vertex.table[, ID := NULL]
    setcolorder(vertex.table, c("atom", setdiff(colnames(vertex.table), "atom")))
    setnames(vertex.table, "metabolite_name", "label")
    setnames(vertex.table, "metabolite_url", "url")
    vertex.table[]
}

#' @export
#' @import igraph
makeAtomGraph <- function(network,
                          org.gatom.anno, gene.de, gene.de.meta,
                          met.db, met.de, met.de.meta) {
    edge.table <- .makeEdgeTable(network=network,
                                 org.gatom.anno=org.gatom.anno,
                                 gene.de=gene.de,
                                 gene.de.meta=gene.de.meta)
    all.atoms <- union(edge.table$atom.x, edge.table$atom.y)
    vertex.table <- .makeVertexTable(network=network,
                                     atoms=all.atoms,
                                     met.db=met.db,
                                     met.de=met.de,
                                     met.de.meta=met.de.meta)
    g <- graph.data.frame(edge.table, directed=FALSE, vertices = vertex.table)
    gc <- components(g)
    g <- induced.subgraph(g, gc$membership == which.max(gc$csize))
    g
}

#' @import BioNet
#' @export
scoreGraph <- function(g, k.gene, k.met,
                       vertex.threshold.min=0.1,
                       edge.threshold.min=0.1,
                       met.score.coef=1) {
    vertex.table <- data.table(as_data_frame(g, what="vertices"))
    edge.table <- data.table(as_data_frame(g, what="edges"))
    if (!is.null(met.de)) {
        pvalsToFit <- vertex.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]

        vertex.bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F)

        vertex.threshold <- if (k.met > length(pvalsToFit)) 1 else {
            sort(pvalsToFit)[k.met]
        }

        vertex.threshold <- min(vertex.threshold,
                                BioNet::fdrThreshold(vertex.threshold.min, vertex.bum))
        .messagef("Metabolite p-value threshold: %f", vertex.threshold)
        .messagef("Metabolite BU alpha: %f", vertex.bum$a)
        V(g)$score <- with(vertex.table,
                           (vertex.bum$a - 1) *
                               (log(.replaceNA(pval, 1)) - log(vertex.threshold)))
        V(g)$score <- V(g)$score * met.score.coef
    } else {
        V(g)$score <- 0
    }

    if (!is.null(gene.de)) {
        pvalsToFit <- edge.table[!is.na(pval)][!duplicated(signal), setNames(pval, signal)]

        edge.bum <- BioNet::fitBumModel(pvalsToFit[pvalsToFit > 0], plot = F)

        edge.threshold <- if (k.gene > length(pvalsToFit)) 1 else {
            sort(pvalsToFit)[k.gene]
        }

        edge.threshold <- min(edge.threshold,
                              BioNet::fdrThreshold(edge.threshold.min, edge.bum))
        .messagef("Gene p-value threshold: %f", edge.threshold)
        .messagef("Gene BU alpha: %f", edge.bum$a)
        E(g)$score <- with(edge.table,
                           (edge.bum$a - 1) *
                               (log(.replaceNA(pval, 1)) - log(edge.threshold)))

    } else {
        E(g)$score <- 0
    }
    g
}

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

#' @export
collapseAtomsIntoMetabolites <- function(m) {
    vertex.table <- data.table(as_data_frame(m, what="vertices"))
    edge.table <- data.table(as_data_frame(m, what="edges"))

    atom2metabolite <- vertex.table[, setNames(metabolite, name)]

    vertex.table.c <- copy(vertex.table)
    vertex.table.c[, name := atom2metabolite[name]]
    vertex.table.c <- unique(vertex.table.c)

    edge.table.c <- copy(edge.table)
    edge.table.c[, from :=  atom2metabolite[from]]
    edge.table.c[, to :=  atom2metabolite[to]]
    edge.table.c[from > to, c("from", "to") := list(to, from)]
    edge.table.c <- unique(edge.table.c)


    res <- graph.data.frame(edge.table.c, directed=FALSE, vertices=vertex.table.c)
    res
}

#' @import igraph
#' @import plyr
#' @export
addHighlyExpressedEdges <- function(m, g, top=3000) {
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
