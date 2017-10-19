library(DLMtool)


# Test that MPs run on the Data files 
testthat::test_file("tests/manual/test-Data.r")

# Test plotting functions (need to add more functions - currently only OM and component objects
testthat::test_file("tests/manual/test-plotting.r")


# Test all OMs with all MPs 
testthat::test_file("tests/manual/test-MPs.r") # takes a while


# Test built-in objects 
testthat::test_file("tests/manual/test-builtin.r") # takes about 9 hours 


# Test all OMs including DLMextra
testthat::test_file("tests/manual/test-overnight.r") # 
