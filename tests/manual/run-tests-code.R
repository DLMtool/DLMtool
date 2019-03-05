library(DLMtool)

options(testthat.output_file = "test-out.xml")
testthat::test_dir('tests/manual/test-code', reporter = "junit")

# testthat::test_file("tests/manual/test-code/test-Data_Functions.R") # Ok 

# testthat::test_file("tests/manual/test-code/test-Data_Plotting.R") # Ok

# testthat::test_file("tests/manual/test-code/test-Fease_Functions.R") # Ok 

# testthat::test_file("tests/manual/test-code/test-Import_Data.R") # Ok

# testthat::test_file("tests/manual/test-code/test-MPs.R") # Ok

# testthat::test_file("tests/manual/test-code/test-MSE_functions.R") # Ok

# testthat::test_file("tests/manual/test-code/test-MSE_Plotting.R") # Ok

# testthat::test_file("tests/manual/test-code/test-OM_functions.R") # Ok 

# testthat::test_file("tests/manual/test-code/test-OM_init_doc.R") # Ok

# testthat::test_file("tests/manual/test-code/test-OM_Plotting.R") # Ok

# testthat::test_file("tests/manual/test-code/test-runMSE.R") # to check still

# testthat::test_file("tests/manual/test-code/test-slotDescription.R") # Ok

# testthat::test_file("tests/manual/test-code/test-cpars.R") # Ok 

# testthat::test_file("tests/manual/test-code/test-PMobjects.R") # Ok 

# testthat::test_file("tests/manual/test-code/test-RealIndices.R") # Ok

# testthat::test_file("tests/manual/test-code/test-Data2csv.R") # Ok

# testthat::test_file("tests/manual/test-code/test-checkPopdyn.R") # Ok




