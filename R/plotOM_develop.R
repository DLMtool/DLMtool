# 
# stocks <- avail("Stock")
# 
# i <- 1
# 
# 
# Stock <- get(stocks[i])
# plot_NaturalMortality <- function(Stock, tabs=FALSE, nsamp=3) {
#   if (!inherits(Stock, "Stock")) stop("Object must be class 'Stock' or 'OM'", call. = FALSE)
#   
#   if (class(Stock) == "OM") {
#     cpars <- Stock@cpars
#     nsim <- Stock@nsim 
#     proyears <- Stock@proyears
#     nyears <- Stock@nyears
#   } 
#   set.seed(ifelse(class(Stock) == "OM", Stock@seed, seed))
#   
#   SampCpars <- list()
#   if(length(cpars)>0) SampCpars <- SampleCpars(cpars, nsim, msg=!silent)
#   StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars, msg=!silent)
#   
#   
#   plotPars <- list(col="gray", axes=FALSE, breaks=10, cex.main=1, lwd=2)
#   
#   its <- sample(1:nsim, nsamp)
#   rmarkdown::render('inst/Rmd/Stock_NaturalMortality.Rmd', params=list(title=Stock@Name,
#                                                                        StockPars=StockPars, 
#                                                                        plotPars=plotPars,
#                                                                        tabs=tabs,
#                                                                        its=its,
#                                                                        nyears=nyears,
#                                                                        proyears=proyears))
# }
# 
# 
# plot_NaturalMortality(Stock)
# plot_NaturalMortality(Stock, tabs=TRUE)
# 
# 
# plot_Growth <- function(Stock, tabs=FALSE, nsamp=3) {
#   if (!inherits(Stock, "Stock")) stop("Object must be class 'Stock' or 'OM'", call. = FALSE)
#   
#   if (class(Stock) == "OM") {
#     cpars <- Stock@cpars
#     nsim <- Stock@nsim 
#     proyears <- Stock@proyears
#     nyears <- Stock@nyears
#   } 
#   set.seed(ifelse(class(Stock) == "OM", Stock@seed, seed))
#   
#   SampCpars <- list()
#   if(length(cpars)>0) SampCpars <- SampleCpars(cpars, nsim, msg=!silent)
#   StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars, msg=!silent)
#   
#   
#   plotPars <- list(col="gray", axes=FALSE, breaks=10, cex.main=1, lwd=2)
#   
#   its <- sample(1:nsim, nsamp)
#   rmarkdown::render('inst/Rmd/Stock_Growth.Rmd', params=list(title=Stock@Name,
#                                                              StockPars=StockPars, 
#                                                              plotPars=plotPars,
#                                                              tabs=tabs,
#                                                              its=its,
#                                                              nyears=nyears,
#                                                              proyears=proyears))
#   
# }
# 
# plot_Growth(Stock)
# plot_Growth(Stock, tabs=TRUE)
# 
# cpars=NULL; 
# nsim=480; nyears=50; proyears=50; seed=101; silent=TRUE
# 
# plot_Stock <- function(Stock, cpars=NULL, nsamp=3, nsim=48, nyears=50, proyears=50, seed=101, silent=TRUE) {
#   if (!inherits(Stock, "Stock")) stop("Object must be class 'Stock' or 'OM'", call. = FALSE)
#   
#   if (class(Stock) == "OM") {
#     cpars <- Stock@cpars
#     nsim <- Stock@nsim 
#     proyears <- Stock@proyears
#     nyears <- Stock@nyears
#   } 
#   set.seed(ifelse(class(Stock) == "OM", Stock@seed, seed))
#   
#   SampCpars <- list()
#   if(length(cpars)>0) SampCpars <- SampleCpars(cpars, nsim, msg=!silent)
#   StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars, msg=!silent)
#   
#   
#   plotPars <- list(col="gray", axes=FALSE, breaks=10, cex.main=1)
#   
#   its <- sample(1:nsim, nsamp)
#   rmarkdown::render('inst/Rmd/plotStock.Rmd', params=list(title=Stock@Name,
#                                                           StockPars=StockPars, 
#                                                           plotPars=plotPars,
#                                                           tabs=TRUE,
#                                                           its=its,
#                                                           nyears=nyears,
#                                                           proyears=proyears))
#   
#   
# }
# 
# plot_Stock(Stock)
# 
# plotStock_Existing <- function(x, nsamp=3, nsim=500, nyears=50, proyears=28, 
#                       col="darkgray", breaks=10, lwd=2, ask=FALSE, incVB=TRUE, ...) {
#   Stock <- x 
#   SampCpars <- list() # empty list 
#   if (class(Stock) == "OM") {
#     Stock <- updateMSE(x) 
#     if (is.finite(Stock@nyears)) nyears <- Stock@nyears
#     if (is.finite(Stock@proyears)) proyears <- Stock@proyears
#     if (is.finite(Stock@nsim)) nsim <- Stock@nsim	
#     
#     if(length(Stock@cpars)>0){ # custom parameters exist - sample and write to list
#       #ncparsim<-cparscheck(Stock@cpars)   # check each list object has the same length and if not stop and error report
#       SampCpars <- SampleCpars(Stock@cpars, nsim, msg=FALSE) 
#     }
#     Stock <- SubOM(Stock)
#   }
#   its <- sample(1:nsim, nsamp)
#   
# }
# 
# 
# plot_Growth(Stock)
# plot_Stock(Stock)
# 
# 
# # Make table of plotting functions
# 
