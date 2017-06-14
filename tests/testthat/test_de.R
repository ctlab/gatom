context("Differential expression utils")

data("org.Mm.eg.gatom.annoEx")
data("gene.de.rawEx")
data("met.kegg.dbEx")

test_that("findIdColumn works", {
    de <- gene.de.rawEx

    idsList <- idsListFromAnnotation(org.gatom.anno = org.Mm.eg.gatom.annoEx)

    idColumn <- findIdColumn(de, idsList)
    expect_true(idColumn$column == "ID")
    expect_true(idColumn$type == "RefSeq")
})


test_that("getGeneDeMeta works", {
    de <- gene.de.rawEx

    de.meta <- getGeneDEMeta(de, org.gatom.anno = org.Mm.eg.gatom.annoEx)
    expect_true(de.meta$idType == "RefSeq")
    expect_true(de.meta$columns$ID == "ID")
    expect_true(de.meta$columns$pval == "pval")
    expect_true(de.meta$columns$baseMean == "baseMean")
    expect_true(de.meta$columns$log2FC == "log2FC")
})

test_that("prepareDE works", {
    de.raw <- gene.de.rawEx

    de.meta <- getGeneDEMeta(de.raw, org.gatom.anno = org.Mm.eg.gatom.annoEx)

    de <- prepareDE(de.raw, de.meta)

    expect_equal(de$logPval, log(de.raw$pval))
})


test_that("getMetDEMeta works for NULL", {
    nullMeta <- getMetDEMeta(NULL, met.db=met.kegg.dbEx)
    expect_is(nullMeta, "NULL")
})

test_that("prepareDE works for NULL", {
    de <- prepareDE(de.raw=NULL, de.meta=NULL)
    expect_is(de, "NULL")
})


