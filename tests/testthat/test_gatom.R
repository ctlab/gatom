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

    gs <- scoreGraph(g, k.gene = 25, k.met=25)

    expect_true(E(gs)[match("Idh1", label)]$score > 0)
    expect_true(E(gs)[match("Gapdh", label)]$score < 0)


    set.seed(42)
    m <- solveSgmwcsRandHeur(gs, max.iterations = 2000)

    expect_true("Idh1" %in% E(m)$label)
})

