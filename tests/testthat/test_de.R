context("Differential expression utils")

load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rda"))

test_that("findIdColumn works", {
    de <- fread(system.file("extdata", "de.samples/Ctrl.vs.MandLPSandIFNg.gene.de.tsv", package="gatom"))

    idsList <- idsListFromAnnotation(org.gatom.anno = org.Mm.eg.gatom.anno)

    idColumn <- findIdColumn(de, idsList)
    expect_true(idColumn$column == "ID")
    expect_true(idColumn$type == "RefSeq")
})


test_that("getGeneDeMeta works", {
    de <- fread(system.file("extdata", "de.samples/Ctrl.vs.MandLPSandIFNg.gene.de.tsv", package="gatom"))

    de.meta <- getGeneDEMeta(de, org.gatom.anno = org.Mm.eg.gatom.anno)
    expect_true(de.meta$idType == "RefSeq")
    expect_true(de.meta$columns$ID == "ID")
    expect_true(de.meta$columns$pval == "pval")
    expect_true(de.meta$columns$baseMean == "baseMean")
    expect_true(de.meta$columns$log2FC == "log2FC")
})

test_that("prepareDE works", {
    de.raw <- fread(system.file("extdata", "de.samples/Ctrl.vs.MandLPSandIFNg.gene.de.tsv", package="gatom"))

    de.meta <- getGeneDEMeta(de.raw, org.gatom.anno = org.Mm.eg.gatom.anno)

    de <- prepareDE(de.raw, de.meta)

    expect_equal(de$logPval, log(de.raw$pval))
})

