context("Plotting functions")

GetMoreData(TRUE)
library(DLMdata)
# setup()


test_that("plotStock works with all available Stock objects", {
  objs <- avail('Stock')
  
  for (i in seq_along(objs)) {
    graphics.off()
    obj <- get(objs[i])
    expect_error(plot(obj), NA, info=objs[i])
    expect_warning(plot(obj), NA, info=objs[i])
    graphics.off()
  }
})

test_that("plotFleet works with all available Fleet objects", {
  objs <- avail('Fleet')
  
  for (i in seq_along(objs)) {
    graphics.off()
    obj <- get(objs[i])
    stock <- sample(avail("Stock"),1)
    info <- paste(objs[i], stock)
    expect_error(plotFleet(obj, get(stock)), NA, info=info)
    # expect_warning(plotFleet(obj, get(stock)), NA, info=info)
    graphics.off()
  }
})

test_that("plotImp works with all available Imp objects", {
  objs <- avail('Imp')
  
  for (i in seq_along(objs)) {
    graphics.off()
    obj <- get(objs[i])
    expect_error(plot(obj), NA, info=objs[i])
    # expect_warning(plot(obj), NA, info=objs[i])
    graphics.off()
  }
})

test_that("plotObs works with all available Obs objects", {
  objs <- avail('Obs')
  
  for (i in seq_along(objs)) {
    graphics.off()
    obj <- get(objs[i])
    expect_error(plot(obj), NA, info=objs[i])
    # expect_warning(plot(obj), NA, info=objs[i])
    graphics.off()
  }
})

test_that("plotOM works with all OMs", {
  objs <- avail('OM')
  
  for (i in seq_along(objs)) {
    graphics.off()
    obj <- get(objs[i])
    expect_error(plot(obj), NA, info=objs[i])
    # expect_warning(plot(obj), NA, info=objs[i])
    graphics.off()
  }
})




