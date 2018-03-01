testthat::context("OM_init_doc")

name <- 'TEST'

testthat::test_that("OM_init works", {
  testthat::expect_error(OMinit(name, overwrite = TRUE), NA)
  testthat::expect_error(OMinit(name, DLMtool::Albacore, overwrite = TRUE), NA)
  testthat::expect_error(OMinit(name, DLMtool::Generic_DecE, overwrite = TRUE), NA)
  testthat::expect_error(OMinit(name, DLMtool::Generic_Obs, overwrite = TRUE), NA)
  testthat::expect_error(OMinit(name, DLMtool::Overages, overwrite = TRUE), NA)
  testthat::expect_error(OMinit(name, DLMtool::testOM, overwrite = TRUE), NA)
})


testthat::test_that("OM_doc works from XL", {
  testthat::expect_error(OMdoc(name, openFile=FALSE, quiet=TRUE), NA)
})

testthat::test_that("XL2OM works", {
  testthat::expect_error(OM <<- XL2OM(name), NA)
})


testthat::test_that("OM_doc works with OM", {
  testthat::expect_error(OMdoc(OM, openFile=FALSE, quiet=TRUE), NA)
})



# file clean-up
file.remove(paste0(name, c(".xlsx", ".rmd", ".html")))
unlink("build", recursive = TRUE)
unlink("images", recursive = TRUE)
unlink("robustness", recursive = TRUE)



if (basename(getwd()) != "DLMtool") {
  unlink("data", recursive = TRUE)
  unlink("docs", recursive = TRUE)
}
