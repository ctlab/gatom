context("Overall pipeline")

data("networkEx")
data("org.Mm.eg.gatom.annoEx")
data("met.kegg.dbEx")

data("met.de.rawEx")
data("gene.de.rawEx")

data("gEx")
data("gsEx")
data("mEx")

test_that("overall pipeline works", {
    g <- makeAtomGraph(network=networkEx,
                       org.gatom.anno=org.Mm.eg.gatom.annoEx,
                       gene.de=gene.de.rawEx,
                       met.db=met.kegg.dbEx,
                       met.de=met.de.rawEx)
    expect_is(g, "igraph")
    expect_true("Idh1" %in% E(g)$label)

    gs <- scoreGraph(g, k.gene = 25, k.met=25, show.warnings = FALSE)

    expect_true(E(gs)[match("Idh1", label)]$score > 0)
    expect_true(E(gs)[match("Gapdh", label)]$score < 0)


    set.seed(42)
    m <- solveSgmwcsRandHeur(gs, max.iterations = 2000)

    expect_true("Idh1" %in% E(m)$label)

    m.ext <- addHighlyExpressedEdges(m, g, top=200)

    addedGenes <- setdiff(E(m.ext)$label, E(m)$label)
    expect_true(all(E(m.ext)[label %in% addedGenes]$signalRank <= 200))

    m.collapsed <- collapseAtomsIntoMetabolites(m.ext)
    expect_equivalent(unique(V(m.ext)$metabolite), V(m.collapsed)$metabolite)
    expect_true(all(!duplicated(V(m.collapsed)$metabolite)))

    m.connected <- connectAtomsInsideMetabolite(m.ext)
    expect_equivalent(V(m.ext)$metabolite, V(m.connected)$metabolite)
    components <- split(V(m)$name, V(m.connected)$metabolite)
    components <- components[sapply(components, length) > 1]
    for (component in components) {
        for (v in component) {
            expect_equivalent(setdiff(component, V(m.connected)[nei(v)]$name), v)
        }
    }


})

test_that("overall pipeline works with data only for genes", {
    g <- makeAtomGraph(network=networkEx,
                       org.gatom.anno=org.Mm.eg.gatom.annoEx,
                       gene.de=gene.de.rawEx,
                       met.db=met.kegg.dbEx,
                       met.de=NULL)
    expect_is(g, "igraph")
    expect_true("Idh1" %in% E(g)$label)

    gs <- scoreGraph(g, k.gene=25, k.met=NULL, show.warnings = FALSE)

    expect_true(E(gs)[match("Idh1", label)]$score > 0)
    expect_true(E(gs)[match("Gapdh", label)]$score < 0)


    set.seed(42)
    m <- solveSgmwcsRandHeur(gs, max.iterations = 2000)

    expect_true("Idh1" %in% E(m)$label)

    m.ext <- addHighlyExpressedEdges(m, g, top=200)

    addedGenes <- setdiff(E(m.ext)$label, E(m)$label)
    expect_true(all(E(m.ext)[label %in% addedGenes]$signalRank <= 200))

    m.collapsed <- collapseAtomsIntoMetabolites(m.ext)
    expect_equivalent(unique(V(m.ext)$metabolite), V(m.collapsed)$metabolite)

    m.connected <- connectAtomsInsideMetabolite(m.ext)
    expect_equivalent(V(m.ext)$metabolite, V(m.connected)$metabolite)
})

test_that("overall pipeline works with data only for metabolites", {
    g <- makeAtomGraph(network=networkEx,
                       org.gatom.anno=org.Mm.eg.gatom.annoEx,
                       gene.de=NULL,
                       met.db=met.kegg.dbEx,
                       met.de=met.de.rawEx)
    expect_is(g, "igraph")
    expect_true("Idh1" %in% E(g)$label)

    gs <- scoreGraph(g, k.gene=NULL, k.met=25, show.warnings = FALSE)

    expect_true(V(gs)[match("Isocitrate", label)]$score > 0)
    expect_true(V(gs)[match("Acetyl-CoA", label)]$score < 0)


    set.seed(42)
    m <- solveSgmwcsRandHeur(gs, max.iterations = 2000)

    # expect_true("Idh1" %in% E(m)$label)

    # no gene expression, no adding edges
    expect_warning(m.ext <- addHighlyExpressedEdges(m, g, top=200))

    expect_equivalent(E(m)$label, E(m.ext)$label)

    m.collapsed <- collapseAtomsIntoMetabolites(m.ext)
    expect_equivalent(unique(V(m.ext)$metabolite), V(m.collapsed)$metabolite)

    m.connected <- connectAtomsInsideMetabolite(m.ext)
    expect_equivalent(V(m.ext)$metabolite, V(m.connected)$metabolite)
})



test_that(".makeVertexTable works with null DE", {
    all.atoms <- networkEx$atoms$atom
    vt <- .makeVertexTable(network=networkEx,
                           atoms=all.atoms,
                           met.db=met.kegg.dbEx,
                           met.de=NULL,
                           met.de.meta=NULL)
})

test_that(".makeVertexTable works with null DE", {
    et <- .makeEdgeTable(network=networkEx,
                         org.gatom.anno=org.Mm.eg.gatom.annoEx,
                         gene.de=NULL,
                         gene.de.meta=NULL)
})

test_that("scoreGraph shows warning on bad distribution", {
    g <- gEx
    tryCatch(gs <- scoreGraph(g, k.gene=NULL, k.met=25),
             warning=fail)

    V(g)$pval <- 1
    expect_warning(gs <- scoreGraph(g, k.gene=NULL, k.met=25))
})

test_that("scoreGraph assigns vertex scores to 0 if its p-value distribution has an inappropriate type", {
    g <- makeAtomGraph(network=networkEx,
                       org.gatom.anno=org.Mm.eg.gatom.annoEx,
                       gene.de=gene.de.rawEx,
                       met.db=met.kegg.dbEx,
                       met.de=met.de.rawEx)

    vertex_attr(g)$pval <- runif(n = length(vertex_attr(g)$pval))

    expect_warning(gs <- scoreGraph(g, k.gene = 25, k.met = 25, show.warnings = F))
})

test_that("scoreGraph assigns edge scores to 0 if its p-value distribution has an inappropriate type", {
    g <- makeAtomGraph(network=networkEx,
                       org.gatom.anno=org.Mm.eg.gatom.annoEx,
                       gene.de=gene.de.rawEx,
                       met.db=met.kegg.dbEx,
                       met.de=met.de.rawEx)

    edge_attr(g)$pval <- runif(n = length(edge_attr(g)$pval))

    expect_warning(gs <- scoreGraph(g, k.gene = 25, k.met = 25, show.warnings = F))
})

test_that("makeAtomGraph works on already good DE tables", {
    gene.de <- data.table(gene.de.rawEx)
    gene.de <- convertPvalDT(gene.de, org.Mm.eg.gatom.annoEx$mapFrom$RefSeq)

    met.de <- data.table(met.de.rawEx)
    met.de <- convertPvalDT(met.de, met.kegg.dbEx$mapFrom$HMDB)

    g <- makeAtomGraph(network=networkEx,
                       org.gatom.anno=org.Mm.eg.gatom.annoEx,
                       gene.de=gene.de,
                       met.db=met.kegg.dbEx,
                       met.de=met.de)
    expect_is(g, "igraph")
    expect_true("Idh1" %in% E(g)$label)
})

test_that("makeAtomGraph notifies if none of the metabolites was masked", {
    expect_warning(g <- makeAtomGraph(network=networkEx,
                                      org.gatom.anno=org.Mm.eg.gatom.annoEx,
                                      gene.de=gene.de.rawEx,
                                      met.db=met.kegg.dbEx,
                                      met.de=met.de.rawEx,
                                      met.to.filter="abrakadabra"))
})

