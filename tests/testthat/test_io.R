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
    saveModuleToPdf(mEx, n_iter=100, force=1e-5, seed=1, f)
    expect_true(file.exists(f))
})
