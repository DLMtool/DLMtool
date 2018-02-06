testthat::context("Import Data objects")

fls <- list.files("../DLMDev/DLMtool_Data/Data")
for (fl in fls) {
  testthat::test_that(paste("new('Data') works with ", fl), {
    testthat::expect_error(Data <- new("Data", file.path("../DLMDev/DLMtool_Data/Data", fl)), NA)
    testthat::expect_error(summary(Data), NA)
  })
}
