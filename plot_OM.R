library(DLMtool)
DLMextra()
library(DLMextra)

OMs <- avail("OM")

OM <- get(OMs[1])


Stock <- SubOM(OM, "Stock")

Object <- Stock

render_plot <- function(Object, RMD=NULL, nsamp=3, nsim=200, nyears=50, proyears=28, 
                        output_file=NULL, output_dir=getwd(), quiet=TRUE,
                        tabs=TRUE, title=NULL, date=NULL,
                        plotPars =NULL) {
  
  SampCpars <- list() # empty list
  
  if (is.null(plotPars)) plotPars <- list(breaks=50, col="grey", axes=TRUE, cex.main=1, lwd=2)
  
  its <- sample(1:nsim, nsamp)
  
  Class <- class(Object)
  if (Class == "Stock") {
    Pars <- SampleStockPars(Object, nsim, nyears, proyears, SampCpars, msg=FALSE)
    Pars$Name <- Object@Name
  } else if (Class == "Fleet") {
   
  } else if (Class == "Obs") {
    
  } else if (Class == "Imp") {
    
  } else {
    stop("Object must be class 'Stock', 'Fleet', 'Obs', or 'Imp'", call.=FALSE)  
  }
  
  Params <- list(
    title = title,
    Pars = Pars,
    plotPars=plotPars,
    tabs = tabs,
    its = its,
    nyears=nyears,
    proyears=proyears,
    date=NULL
  )
  
  outname <- paste0("_", RMD, ".html")
  if (is.null(output_file)) output_file <- paste0(Pars$Name, outname)
  message("Rendering HTML file")
  
  RMD <- paste0(RMD, ".Rmd")
  input <- file.path(system.file(package = 'DLMtool'),'Rmd', Class, RMD)
  
  rend <- try(rmarkdown::render(input, params=Params,
                    output_file=output_file,
                    output_dir=output_dir,
                    quiet=quiet), silent=TRUE)
  
  if (class(rend) == "try-error") {
    print(rend)
  } else {
    message("Rendered ", output_file, " in ", output_dir)
  }
  
}

plot_M <- function(Stock=NULL, nsamp=3, nsim=200, nyears=50, proyears=28, 
                   output_file=NULL, output_dir=getwd(), quiet=TRUE, tabs=TRUE, title=NULL,
                   date=NULL, plotPars=NULL) {
  
  SampCpars <- list() # empty list
  if (class(Stock) == "OM") {
    Stock <- updateMSE(Stock) 
    if (is.finite(Stock@nyears)) nyears <- Stock@nyears
    if (is.finite(Stock@proyears)) proyears <- Stock@proyears
    if (is.finite(Stock@nsim)) nsim <- Stock@nsim	
    
    if(length(Stock@cpars)>0){ # custom parameters exist - sample and write to list
      #ncparsim<-cparscheck(Stock@cpars)   # check each list object has the same length and if not stop and error report
      SampCpars <- SampleCpars(Stock@cpars, nsim, msg=FALSE) 
    }
    Stock <- SubOM(Stock)
  } 
  if (class(Stock) != "Stock") stop("Object must be class 'Stock' or 'OM'", call.=FALSE)
  
  render_plot(Object=Stock, RMD='NaturalMortality', nsamp=nsamp, nsim=nsim, 
              nyears=nyears, proyears=proyears,
              output_file=output_file, output_dir=output_dir, quiet=quiet,
              tabs=tabs, title=title, date=date,
              plotPars=plotPars)

  
}

plot_M(Stock)





plot_Growth <- function(Stock=NULL, nsamp=3, nsim=200, nyears=50, proyears=28, 
                        output_file=NULL, output_dir=getwd(), quiet=TRUE) {
  SampCpars <- list() # empty list
  plotPars <- list(breaks=50, col="grey", axes=TRUE, cex.main=1, lwd=2)
  
  if (class(Stock) == "OM") {
    Stock <- updateMSE(Stock) 
    if (is.finite(Stock@nyears)) nyears <- Stock@nyears
    if (is.finite(Stock@proyears)) proyears <- Stock@proyears
    if (is.finite(Stock@nsim)) nsim <- Stock@nsim	
    
    if(length(Stock@cpars)>0){ # custom parameters exist - sample and write to list
      #ncparsim<-cparscheck(Stock@cpars)   # check each list object has the same length and if not stop and error report
      SampCpars <- SampleCpars(Stock@cpars, nsim, msg=FALSE) 
    }
    Stock <- SubOM(Stock)
  } 
  if (class(Stock) == "Stock") {
    StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars, msg=FALSE)
    StockPars$Name <- Stock@Name
    its <- sample(1:nsim, nsamp)
    Params <- list(
      title = NULL,
      StockPars = StockPars,
      plotPars=plotPars,
      tabs = TRUE,
      its = its,
      nyears=nyears,
      proyears=proyears,
      date=NULL
    )
  } 
  
  if (class(Stock) == "list") {
    # check that list is ok 
    if(!all(c("title", "StockPars", "plotPars", "tabs", "its", "nyears","proyears") %in% names(Stock)))
      stop("Object must be class 'Stock' or 'OM', or a named list", call.=FALSE)
    Params <- Stock 
  }
  
  if (is.null(output_file)) output_file <- paste0(StockPars$Name, "_Growth.html")
  message("Rendering ", output_file, " in ", output_dir)
  
  rmarkdown::render("inst/Rmd/Stock/Growth.Rmd", params=Params,
                    output_file=output_file,
                    output_dir=output_dir,
                    quiet=quiet)
  
}


plot_Growth(Stock, nsamp=5)
plot_M(Params)






plot_Stock <- function(Stock, nsamp=3, nsim=500, nyears=50, proyears=28) {
  
  SampCpars <- list() # empty list
  
  if (class(Stock) == "OM") {
    Stock <- updateMSE(x) 
    if (is.finite(Stock@nyears)) nyears <- Stock@nyears
    if (is.finite(Stock@proyears)) proyears <- Stock@proyears
    if (is.finite(Stock@nsim)) nsim <- Stock@nsim	
    
    if(length(Stock@cpars)>0){ # custom parameters exist - sample and write to list
      #ncparsim<-cparscheck(Stock@cpars)   # check each list object has the same length and if not stop and error report
      SampCpars <- SampleCpars(Stock@cpars, nsim, msg=FALSE) 
    }
    Stock <- SubOM(Stock)
  }
  if (class(Stock) != "Stock") stop("Object must be class 'Stock' or 'OM'", call.=FALSE)
  
  
  # --- Sample Stock Parameters ----
  StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars, msg=FALSE)
  
  plotPars <- list(breaks=50, col="grey", axes=TRUE, cex.main=1, lwd=2)
  

  its <- 1:3
  Params <- list(
    title = "Title" ,
    StockPars = StockPars,
    plotPars=plotPars,
    tabs = TRUE,
    its = its,
    nyears=nyears,
    proyears=proyears
    
  )

  rmarkdown::render("inst/Rmd/Stock_NaturalMortality.Rmd", params=Params)
  
  
  rmarkdown::render("inst/Rmd/test.Rmd", params=list(StockPars))
  
  
  
}



plotOM <- function(OM) {
  
  
  
  
  
  
  
}
