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

test_that("saveModuleToPdf works", {
    f <- tempfile()
    saveModuleToPdf(mEx, f, n_iter=100, force=1e-5)
    expect_true(file.exists(f))
})

test_that("saveModuleToHtml works", {
    f <- tempfile()
    saveModuleToHtml(mEx, f)
    expect_true(file.exists(f))
})

test_that("Node attributes are created", {
    vertex.table <- as_data_frame(mEx, what="vertices")
    res <- getJsNodeStyleAttributes(vertex.table)
    expect_true(!is.null(res$tooltip))
})

test_that("Edge attributes are created", {
    edge.table <- as_data_frame(mEx, what="edges")
    res <- getJsEdgeStyleAttributes(edge.table)
    expect_true(!is.null(res$tooltip))
})

