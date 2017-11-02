context("IO")
data("mEx")

test_that("saveModuleToXgmml works", {
    f <- tempfile()
    saveModuleToXgmml(mEx, f)
    expect_true(file.exists(f))
})

test_that("saveModuleToDot works", {
    f <- tempfile()
    saveModuleToDot(mEx, f)
    expect_true(file.exists(f))
})

test_that("getDotNodeStyleAttributes and getDotEdgeStyleAttributes work", {
    g <- makeAtomGraph(network=networkEx,
                       org.gatom.anno=org.Mm.eg.gatom.annoEx,
                       gene.de=gene.de.rawEx,
                       met.db=met.kegg.dbEx,
                       met.de=NULL)
    gs <- scoreGraph(g, k.gene=25, k.met=NULL, show.warnings = FALSE)
    set.seed(42)
    m <- solveSgmwcsRandHeur(gs, max.iterations = 2000)
    options(stringsAsFactors = FALSE)
    produceNodeAttrs <- getDotNodeStyleAttributes(as_data_frame(m, what = "vertices"))


    g <- makeAtomGraph(network=networkEx,
                       org.gatom.anno=org.Mm.eg.gatom.annoEx,
                       gene.de=NULL,
                       met.db=met.kegg.dbEx,
                       met.de=met.de.rawEx)
    gs <- scoreGraph(g, k.gene=NULL, k.met=25, show.warnings = FALSE)
    set.seed(42)
    m <- solveSgmwcsRandHeur(gs, max.iterations = 2000)
    options(stringsAsFactors = FALSE)
    produceEdgeAttrs <- getDotEdgeStyleAttributes(as_data_frame(m))
})
