library(DLMtool)
library(testthat)

# Test that MPs run on the Data files 
test_file("manual/test-Data.r")

# Test plotting functions (need to add more functions - currently only OM and component objects
test_file("manual/test-plotting.r")

# Test built-in objects 
test_file("manual/test-builtin.r") # takes about 9 hours 

# Test all OMs with all MPs 
test_file("manual/test-overnight.r") # takes a long time!