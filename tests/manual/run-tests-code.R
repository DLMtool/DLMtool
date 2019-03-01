library(DLMtool)

options(testthat.output_file = "test-out.xml")
testthat::test_dir('tests/manual/test-code', reporter = "junit")

testthat::test_file("tests/manual/test-code/test-Data_Functions.R") #

testthat::test_file("tests/manual/test-code/test-Data_Plotting.R") #

testthat::test_file("tests/manual/test-code/test-Fease_Functions.R") #

testthat::test_file("tests/manual/test-code/test-Import_Data.R") #

testthat::test_file("tests/manual/test-code/test-MPs.R") #

testthat::test_file("tests/manual/test-code/test-MSE_functions.R") #

testthat::test_file("tests/manual/test-code/test-MSE_Plotting.R") #

testthat::test_file("tests/manual/test-code/test-OM_functions.R") #

testthat::test_file("tests/manual/test-code/test-OM_init_doc.R") #

testthat::test_file("tests/manual/test-code/test-OM_Plotting.R") # takes a while - some fail

testthat::test_file("tests/manual/test-code/test-runMSE.R") # to check still

testthat::test_file("tests/manual/test-code/test-slotDescription.R") #

testthat::test_file("tests/manual/test-code/test-cpars.R") #

testthat::test_file("tests/manual/test-code/test-PMobjects.R") #

testthat::test_file("tests/manual/test-code/test-RealIndices.R") #

testthat::test_file("tests/manual/test-code/test-Data2csv.R") #

testthat::test_file("tests/manual/test-code/test-checkPopdyn.R") #




