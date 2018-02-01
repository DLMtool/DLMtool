library(DLMtool)

# Test that runMSE works with testOM and default MPs
testthat::test_file("tests/manual/test_01_runMSE_testOM.r")

# Test that runMSE works with 9 random builtin OMs and default MPs
testthat::test_file("tests/manual/test_02_runMSE_9_builtin_OMs.r")












