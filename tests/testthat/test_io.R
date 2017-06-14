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
