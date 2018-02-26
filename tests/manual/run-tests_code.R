library(DLMtool)
library(testthat)

testthat::test_file("tests/manual/test-code/test-Data_Functions.R") # OK 

testthat::test_file("tests/manual/test-code/test-Data_Plotting.R") # OK 

testthat::test_file("tests/manual/test-code/test-Import_Data.R") # OK 

testthat::test_file("tests/manual/test-code/test-MPs.R") # OK 

testthat::test_file("tests/manual/test-code/test-MSE_functions.R") # OK 

testthat::test_file("tests/manual/test-code/test-MSE_Plotting.R") # OK 

testthat::test_file("tests/manual/test-code/test-OM_functions.R") # OK

testthat::test_file("tests/manual/test-code/test-OM_init_doc.R") # OK 

testthat::test_file("tests/manual/test-code/test-OM_Plotting.R") # OK 

testthat::test_file("tests/manual/test-code/test-runMSE.R") # 

testthat::test_file("tests/manual/test-code/test-slotDescription.R") # 



testthat::test_dir('tests/manual/test-code')
                   

# source("tests/manual/run-coverage_code.R")


