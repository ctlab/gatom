context("Overall pipeline")

data("networkEx")
data("org.Mm.eg.gatom.annoEx")
data("met.kegg.dbEx")

data("met.de.rawEx")
data("gene.de.rawEx")


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

    m.collapsed <- collapseAtomsIntoMetabolites(m)
    expect_equivalent(unique(V(m)$metabolite), V(m.collapsed)$metabolite)
    expect_true(all(!duplicated(V(m.collapsed)$metabolite)))

    m.connected <- connectAtomsInsideMetabolite(m)
    expect_equivalent(V(m)$metabolite, V(m.connected)$metabolite)
    components <- split(V(m)$name, V(m.connected)$metabolite)
    components <- components[sapply(components, length) > 1]
    for (component in components) {
        for (v in component) {
            expect_equivalent(setdiff(component, V(m.connected)[nei(v)]$name), v)
        }
    }
})

test_that("overall pipeline works without met data", {
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
})

test_that(".makeVertexTable works with null DE", {
    all.atoms <- networkEx$atoms$atom
    vt <- .makeVertexTable(network=networkEx,
                           atoms=all.atoms,
                           met.db=met.kegg.dbEx,
                           met.de=NULL,
                           met.de.meta=NULL)
})
