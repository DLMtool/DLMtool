library(DLMtool)

options(testthat.output_file = "test-out.xml")
testthat::test_dir('tests/manual/test-code', reporter = "junit")

options(testthat.output_file = "test-1.rda")
tt <- testthat::test_dir('tests/manual/test-1') #, reporter=ListReporter)


t1 <- readRDS(file.path(getwd(), '/tests/manual/test-1/test-1.rda'))


options(testthat.output_file = "tests-result.xml")
tt <- testthat::test_file("tests/manual/test-code/test-checkPopdyn.R"  ) # ok

t1 <- testthat::test_file("tests/manual/test-code/test-cpars.R", reporter = RstudioReporter) # ok

t1



# testthat::test_file("tests/manual/test-code/test-checkPopdyn.R") # ok

# testthat::test_file("tests/manual/test-code/test-cpars.R") # ok

# testthat::test_file("tests/manual/test-code/test-Data_Functions.R") # ok

# testthat::test_file("tests/manual/test-code/test-Data_Plotting.R") # ok

# testthat::test_file("tests/manual/test-code/test-Fease_Functions.R") # ok  

# testthat::test_file("tests/manual/test-code/test-Import_Data.R") # ok 

# testthat::test_file("tests/manual/test-code/test-MPs.R") # ok 

# testthat::test_file("tests/manual/test-code/test-MSE_functions.R") # ok

# testthat::test_file("tests/manual/test-code/test-MSE_Plotting.R") # ok

# testthat::test_file("tests/manual/test-code/test-OM_functions.R") # ok

# testthat::test_file("tests/manual/test-code/test-OM_init_doc.R") # ok

# testthat::test_file("tests/manual/test-code/test-OM_Plotting.R") # ok

# testthat::test_file("tests/manual/test-code/test-runMSE.R") # ok

# testthat::test_file("tests/manual/test-code/test-slotDescription.R") # ok

# testthat::test_file("tests/manual/test-code/test-PMobjects.R") # ok

# testthat::test_file("tests/manual/test-code/test-Data2csv.R") # ok

# testthat::test_file("tests/manual/test-code/test-RealIndices.R") # not used








