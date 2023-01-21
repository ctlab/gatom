context("Annotation")

test_that("makeOrgGatomAnnotation works", {
    require("org.Mm.eg.db")

    anno <- makeOrgGatomAnnotation(org.Mm.eg.db,
                                   idColumns=
                                       c("Entrez"="ENTREZID",
                                         "Symbol"="SYMBOL"))
    expect_equal(anno$baseId, "Entrez")
})

test_that("getMetabolicPathways works", {
    require("org.Mm.eg.db")

    universe <- keys(org.Mm.eg.db, "ENTREZID")
    metGenes <- mappedkeys(org.Mm.egENZYME)
    pathways <- getMetabolicPathways(universe,
                                     metGenes,
                                     keggOrgCode="mmu")

    expect_true("mmu00010: Glycolysis / Gluconeogenesis" %in% names(pathways))
})
