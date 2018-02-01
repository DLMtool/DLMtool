library(DLMtool)

# Test that runMSE works with all builtin OMs and default MPs
testthat::test_file("tests/manual/test_02a_runMSE_all_builtin_OMs.r")

# Test plotting functions (need to add more functions - currently only OM and component objects
testthat::test_file("tests/manual/test_03a_testPlotting.r")

# Test runMSE with all OMs and all MPs
testthat::test_file("tests/manual/test_04a_runMSE_allOM_allMP.r")

# Test Can, Cant, etc with all Data objects
testthat::test_file("tests/manual/test_05a_Can_Cant_allData.r")


