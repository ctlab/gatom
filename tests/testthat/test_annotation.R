context("Annotation")

test_that("makeOrgGatomAnnotation works", {
    require("org.Mm.eg.db")

    anno <- makeOrgGatomAnnotation(org.Mm.eg.db,
                                   idColumns=
                                       c("Entrez"="ENTREZID",
                                         "Symbol"="SYMBOL"))
    expect_equal(anno$baseId, "Entrez")
})

test_that("getPathways2annotation works", {
    require("org.Mm.eg.db")

    anno <- makeOrgGatomAnnotation(org.Mm.eg.db,
                                   idColumns=
                                       c("Entrez"="ENTREZID",
                                         "Symbol"="SYMBOL"))
    anno$pathways <- getPathways2annotation(anno, organism="mmu")
    expect_true("mmu00010: Glycolysis / Gluconeogenesis" %in% names(anno$pathways))
})
