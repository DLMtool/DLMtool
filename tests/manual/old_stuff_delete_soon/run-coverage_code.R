
library(testthat)
if (!suppressWarnings(require(covr))) {
  install.packages("covr")
  library(covr)
}

tests <- covr::package_coverage(type="none", code="testthat::test_dir('tests/manual/test-code')",
                                quiet = FALSE)

covr::report(tests, file='tests/manual/DLMtool-report.html')
