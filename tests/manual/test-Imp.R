
nsim <- 9
setup()


# --- Output imp -----------------------
context("Test TAC Implementation Error model")

OM<-new('OM',Albacore)
OM@seed<-ceiling(runif(1, 1, 1000))
OM@nsim <- nsim
info <- paste("seed=", OM@seed)
MPs<-c("AvC","FMSYref","DCAC")
nMPs<-length(MPs)

test1<-runMSEnomsg(OM,MPs) # Base

OM@TACFrac<-c(1.1,1.2)
test2<-runMSEnomsg(OM,MPs) # 10 - 20% overage

OM@TACFrac<-c(1.2,1.4)
test3<-runMSEnomsg(OM,MPs) # 20 - 40% overage

OM@TACFrac<-c(0.8,0.9)
test4<-runMSEnomsg(OM,MPs) # 10 - 20% underage

OM@TACSD<-c(0.4,0.5)
test5<-runMSEnomsg(OM,MPs) # Implementation uncertainty


test_that("Overage TAC results in lower PNOF", {
  expect_true(sum(NOAA_plot(test1)$PNOF>NOAA_plot(test2)$PNOF)==nMPs, info=info)
})

test_that("More overage should be even lower PNOF", {
  expect_true(sum(NOAA_plot(test2)$PNOF>NOAA_plot(test3)$PNOF)==nMPs, info=info)
})

test_that("Underage should be higher PNOF", {
  expect_true(sum(NOAA_plot(test4)$PNOF>NOAA_plot(test1)$PNOF)==nMPs, info=info)
})

test_that("Variability should be lower VY metric", {
  expect_true(sum(NOAA_plot(test3)$VY>NOAA_plot(test5)$VY)==nMPs, info=info)
})

test_that("Prob Biomass should be higher", {
  expect_true(sum(NOAA_plot(test1)$B50>NOAA_plot(test2)$B50)==nMPs, info=info)
})

test_that("Prob Biomass should be higher", {
  expect_true(sum(NOAA_plot(test2)$B50>NOAA_plot(test3)$B50)==nMPs, info=info)
})

test_that("Prob Biomass should be higher", {
  expect_true(sum(NOAA_plot(test4)$B50>=NOAA_plot(test1)$B50)==nMPs, info=info)
})



# --- Size limitimp -------------------
context("Test Size Limit Implementation Error model")
OM<-new('OM',Albacore)
OM@seed<-ceiling(runif(1, 1, 1000))
OM@nsim <- nsim
info <- paste("seed=", OM@seed)
MPs<-c("matlenlim","matlenlim2")
nMPs<-length(MPs)

test1<-runMSEnomsg(OM,MPs)

OM@SizeLimFrac<-c(1.2,1.2)

test2<-runMSEnomsg(OM,MPs)

test_that("Prob Not Overfishing should be higher when size limit is positively biased", {
  expect_true(sum(NOAA_plot(test2)$PNOF>NOAA_plot(test1)$PNOF)==nMPs, info=info)
})

OM@SizeLimFrac<-c(0.8,0.8)

test2<-runMSEnomsg(OM,MPs)

test_that("Prob Not Overfishing should be lower when size limit is negatively biased", {
  expect_true(sum(NOAA_plot(test2)$PNOF<NOAA_plot(test1)$PNOF)==nMPs, info=info)
})


# --- Effort imp ----------------------
context("Test Effort Implementation Error model")
OM<-new('OM',Albacore)
OM@seed<-ceiling(runif(1, 1, 1000))
OM@nsim <- nsim
info <- paste("seed=", OM@seed)
MPs<-c("matlenlim","matlenlim2")
nMPs<-length(MPs)

test1<-runMSEnomsg(OM,MPs)

OM@EFrac <- c(1.2, 1.4)
test2<-runMSEnomsg(OM,MPs)

test_that("Prob Not Overfishing should be lower when Effort implementation is positively biased", {
  expect_true(sum(NOAA_plot(test2)$PNOF<NOAA_plot(test1)$PNOF)==nMPs, info=info)
})

OM@EFrac <- c(0.6, 0.8)
test2<-runMSEnomsg(OM,MPs)
test_that("Prob Not Overfishing should be higher when size limit is negatively biased", {
  expect_true(sum(NOAA_plot(test2)$PNOF>NOAA_plot(test1)$PNOF)==nMPs, info=info)
})



