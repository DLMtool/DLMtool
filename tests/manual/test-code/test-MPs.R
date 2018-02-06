testthat::context("test MPs with testOM")

OM <- DLMtool::testOM 
OM@nsim <- 6

MPs <- avail("MP")
output <- avail('Output')
input <- avail('Input')

mixed <- avail('Mixed')
reference <- avail('Reference')

testthat::test_that("avail MP returns correct length", {
  testthat::expect_equal(length(c(output, input, mixed, reference)), length(MPs))
})

                    
for (mm in MPs) {
  testthat::test_that(paste("testOM works with ", mm), {
    testthat::expect_is(runMSE(OM, MPs=mm, silent=TRUE), "MSE")
  })
}
 
                   
                    

