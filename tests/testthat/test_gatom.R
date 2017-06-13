context("Overall pipeline")

data("networkEx")
data("org.Mm.eg.gatom.annoEx")
data("met.kegg.dbEx")

data("met.de.rawEx")
data("gene.de.rawEx")


test_that("gatom works", {
    g <- makeAtomGraph(network=networkEx,
                       org.gatom.anno=org.Mm.eg.gatom.annoEx,
                       gene.de=gene.de.rawEx,
                       met.db=met.kegg.dbEx,
                       met.de=met.de.rawEx)
    expect_is(g, "igraph")
    expect_true("Idh1" %in% E(g)$label)
})

