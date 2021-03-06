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

test_that("findIdColumn works with sampling", {
    de <- gene.de.rawEx

    idsList <- idsListFromAnnotation(org.gatom.anno = org.Mm.eg.gatom.annoEx)

    idColumn <- findIdColumn(de, idsList, sample.size=100)
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

test_that("prepareDEColumn renames conflicting names", {
    de <- data.table(ID=c("a", "b"), symbol=c("x", "y"))
    prepareDEColumn(de, "ID", "symbol")
    expect_true(all(!duplicated(names(de))))
})


test_that("findIdColumns have order of preferences", {
    de <- data.table(v1=c("a", "b", "c"),
                     v2=c("x", "x", "y"),
                     v3=c("qwdas", "as", "qq"),
                     v4=c("2", "1", "3"))

    idsList <- list(
        "d"=c("1", "2", "3"),
        "A"=c("a", "b", "c"),
        "X"=c("x", "y", "z"))


    de1 <- rename(de, c(v2="id"))

    # if there is column with base IDs, we selecting it
    expect_equal(findIdColumn(de1, idsList)$column, "v4")

    de2 <- de1; de2$v4 <- NULL
    # if there is ID column it's preferred when good enough
    expect_equal(findIdColumn(de2, idsList)$column, "id")

    # no ID, no base IDs, next priority is columns with unique values
    de3 <- rename(de2, c(id="bla"))
    expect_equal(findIdColumn(de3, idsList)$column, "v1")

    de4 <- de3; de4$v1 <- NULL
    # last call: whatever matches
    expect_equal(findIdColumn(de4, idsList)$column, "bla")
})
