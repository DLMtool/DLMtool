testthat::context("Test PM Functions")

PMs <- avail("PM")

MSE <- runMSE()

for (pm in PMs) {
  testthat::test_that(paste("Test ", pm, " works with 6 MPs"), {
    testthat::expect_error( get(pm)(MSE), NA)
    pmout <- get(pm)(MSE)
    testthat::expect_equal(length(pmout@Mean), MSE@nMPs)
  })
}

MSE2 <- Sub(MSE, MP=1)

for (pm in PMs) {
  testthat::test_that(paste("Test ", pm, " works with 1 MP"), {
    testthat::expect_error( get(pm)(MSE2), NA)
    pmout <- get(pm)(MSE2)
    testthat::expect_equal(length(pmout@Mean), MSE2@nMPs)
  })
}


# AAVY(MSE)
# AAVY(MSE2)
# 
# LTY(MSE)
# LTY(MSE2)
# 
# P10(MSE)
# P10(MSE2)
# 
# P100(MSE)
# P100(MSE2)
# 
# P50(MSE)
# P50(MSE2)
# 
# PNOF(MSE)
# PNOF(MSE2)
# 
# STY(MSE)
# STY(MSE2)
# 
# Yield(MSE)
# Yield(MSE2)

