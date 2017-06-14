context("Annotation")

test_that("makeOrgGatomAnnotation works", {
    require("org.Mm.eg.db")

    anno <- makeOrgGatomAnnotation(org.Mm.eg.db,
                                   idColumns=
                                       c("Entrez"="ENTREZID",
                                         "Symbol"="SYMBOL"))
    expect_equal(anno$baseId, "Entrez")
})

