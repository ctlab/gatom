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
    g <- makeMetabolicGraph(network=networkEx,
                            topology="atoms",
                            org.gatom.anno=org.Mm.eg.gatom.annoEx,
                            gene.de=gene.de.rawEx,
                            met.db=met.kegg.dbEx,
                            met.de=met.de.rawEx)
    expect_is(g, "igraph")
    expect_true("Idh1" %in% E(g)$label)

    gs <- scoreGraph(g, k.gene = 25, k.met=25, show.warnings = FALSE)

    expect_true(E(gs)[match("Idh1", label)]$score > 0)
    expect_true(E(gs)[match("Gapdh", label)]$score < 0)


    vhsolver <- mwcsr::rnc_solver()
    m <- mwcsr::solve_mwcsp(vhsolver, gs)$graph

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
    g <- makeMetabolicGraph(network=networkEx,
                            topology="atoms",
                            org.gatom.anno=org.Mm.eg.gatom.annoEx,
                            gene.de=gene.de.rawEx,
                            met.db=met.kegg.dbEx,
                            met.de=NULL)
    expect_is(g, "igraph")
    expect_true("Idh1" %in% E(g)$label)

    gs <- scoreGraph(g, k.gene=25, k.met=NULL, show.warnings = FALSE)

    expect_true(E(gs)[match("Idh1", label)]$score > 0)
    expect_true(E(gs)[match("Gapdh", label)]$score < 0)


    vhsolver <- mwcsr::rnc_solver()
    m <- mwcsr::solve_mwcsp(vhsolver, gs)$graph


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
    g <- makeMetabolicGraph(network=networkEx,
                            topology="atoms",
                            org.gatom.anno=org.Mm.eg.gatom.annoEx,
                            gene.de=NULL,
                            met.db=met.kegg.dbEx,
                            met.de=met.de.rawEx)
    expect_is(g, "igraph")
    expect_true("Idh1" %in% E(g)$label)

    gs <- scoreGraph(g, k.gene=NULL, k.met=25, show.warnings = FALSE)

    expect_true(V(gs)[match("Isocitrate", label)]$score > 0)
    expect_true(V(gs)[match("Acetyl-CoA", label)]$score < 0)


    vhsolver <- mwcsr::rnc_solver()
    m <- mwcsr::solve_mwcsp(vhsolver, gs)$graph

    # expect_true("Idh1" %in% E(m)$label)

    # no gene expression, no adding edges
    expect_warning(m.ext <- addHighlyExpressedEdges(m, g, top=200))

    expect_equivalent(E(m)$label, E(m.ext)$label)

    m.collapsed <- collapseAtomsIntoMetabolites(m.ext)
    expect_equivalent(unique(V(m.ext)$metabolite), V(m.collapsed)$metabolite)

    m.connected <- connectAtomsInsideMetabolite(m.ext)
    expect_equivalent(V(m.ext)$metabolite, V(m.connected)$metabolite)
})

test_that("overall pipeline works with nodes as metabolites", {
    g <- makeMetabolicGraph(network=networkEx,
                            topology="atoms",
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
    vhsolver <- mwcsr::rnc_solver()
    res <- mwcsr::solve_mwcsp(vhsolver, gs)
    m <- res$graph

    # expect_true("Idh1" %in% E(m)$label)

    # no gene expression, no adding edges
    expect_warning(m.ext <- addHighlyExpressedEdges(m, g, top=200))

    expect_equivalent(E(m)$label, E(m.ext)$label)
})

test_that(".makeVertexTable works with null DE", {
    all.atoms <- networkEx$atoms$atom
    vt <- .makeVertexTable(network=networkEx,
                           atoms=all.atoms,
                           met.db=met.kegg.dbEx,
                           met.de=NULL,
                           met.de.meta=NULL)
    expect_true(!is.null(vt))
})

test_that(".makeVertexTable works with null DE", {
    et <- .makeEdgeTable(network=networkEx,
                         org.gatom.anno=org.Mm.eg.gatom.annoEx,
                         gene.de=NULL,
                         gene.de.meta=NULL)
    expect_true(!is.null(et))
})

test_that("scoreGraph shows warning on bad distribution", {
    g <- gEx
    tryCatch(gs <- scoreGraph(g, k.gene=NULL, k.met=25),
             warning=fail)

    V(g)$pval <- 1
    expect_warning(gs <- scoreGraph(g, k.gene=NULL, k.met=25))
})

test_that("scoreGraph assigns vertex scores to 0 if its p-value distribution has an inappropriate type", {
    g <- makeMetabolicGraph(network=networkEx,
                            topology="atoms",
                            org.gatom.anno=org.Mm.eg.gatom.annoEx,
                            gene.de=gene.de.rawEx,
                            met.db=met.kegg.dbEx,
                            met.de=met.de.rawEx)

    set.seed(42)
    vertex_attr(g)$pval <- rbeta(n = length(vertex_attr(g)$pval),
                                 shape1 = 2, shape2 = 1)

    expect_warning(gs <- scoreGraph(g, k.gene = 25, k.met = 25, show.warnings = T))
})

test_that("scoreGraph assigns edge scores to 0 if its p-value distribution has an inappropriate type", {
    g <- makeMetabolicGraph(network=networkEx,
                            topology="atoms",
                            org.gatom.anno=org.Mm.eg.gatom.annoEx,
                            gene.de=gene.de.rawEx,
                            met.db=met.kegg.dbEx,
                            met.de=met.de.rawEx)

    set.seed(42)
    edge_attr(g)$pval <- rbeta(n = length(edge_attr(g)$pval),
                               shape1 = 2, shape2 = 1)

    expect_warning(gs <- scoreGraph(g, k.gene = 25, k.met = 25, show.warnings = T))
})

test_that("makeMetabolicGraph works on already good DE tables", {
    gene.de <- data.table(gene.de.rawEx)
    gene.de <- convertPvalDT(gene.de, org.Mm.eg.gatom.annoEx$mapFrom$RefSeq)

    met.de <- data.table(met.de.rawEx)
    met.de <- convertPvalDT(met.de, met.kegg.dbEx$mapFrom$HMDB)

    g <- makeMetabolicGraph(network=networkEx,
                            topology="atoms",
                            org.gatom.anno=org.Mm.eg.gatom.annoEx,
                            gene.de=gene.de,
                            met.db=met.kegg.dbEx,
                            met.de=met.de)
    expect_is(g, "igraph")
    expect_true("Idh1" %in% E(g)$label)
})

test_that("makeMetabolicGraph notifies if none of the metabolites was masked", {
    expect_warning(g <- makeMetabolicGraph(network=networkEx,
                                           topology="atoms",
                                           org.gatom.anno=org.Mm.eg.gatom.annoEx,
                                           gene.de=gene.de.rawEx,
                                           met.db=met.kegg.dbEx,
                                           met.de=met.de.rawEx,
                                           met.to.filter="abrakadabra"))
})

test_that("makeMetabolicGraph fails if topology is misspelled", {
    expect_error(g <- makeMetabolicGraph(network=networkEx,
                                         topology="abrakadabra",
                                         org.gatom.anno=org.Mm.eg.gatom.annoEx,
                                         gene.de=gene.de.rawEx,
                                         met.db=met.kegg.dbEx,
                                         met.de=met.de.rawEx))
})
