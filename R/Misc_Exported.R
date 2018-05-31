
#' What objects of this class are available
#' 
#' Generic class finder
#' 
#' Finds objects of the specified class in the global environment or the
#' DLMtool package.
#' 
#' @param classy A class of object (character string, e.g. 'Fleet')
#' @examples
#' avail("OM")
#' @author T. Carruthers
#' @seealso \link{Can} \link{Cant} \link{avail}
#' @examples 
#' Stocks <- avail("Stock")
#' Fleets <- avail("Fleet")
#' MPs <- avail("MP")
#' @export 
avail <- function(classy) {
  temp <- try(class(classy), silent=TRUE)
  if (class(temp) == "try-error") classy <- deparse(substitute(classy))
  if (temp == "function") classy <- deparse(substitute(classy))
  
  if (classy %in% c('Output', 'Input', "Mixed", "Reference")) {
    MPs <- avail('MP')
    gettype <- MPtype(MPs)
    temp <- gettype[gettype[,2] %in% classy,1]
    if (length(temp) < 1) stop("No MPs of type '", classy, "' found", call. = FALSE)
    return(temp)
    
  } else {
    temp <- c(ls("package:DLMtool")[vapply(ls("package:DLMtool"), getclass, logical(1), classy = classy)], 
              ls(envir = .GlobalEnv)[vapply(ls(envir = .GlobalEnv), getclass, logical(1), classy = classy)])
    pkgs <- search()
    if ("package:DLMextra" %in% pkgs) {
      temp_extra <- ls("package:DLMextra")[vapply(ls("package:DLMextra"), getclass, logical(1), classy = classy)]
      temp <- c(temp, temp_extra)
    }
    
    if (classy == "Observation") message("Class 'Observation' has been re-named 'Obs'")	
    if (length(temp) < 1) stop("No objects of class '", classy, "' found", call. = FALSE)
    return(unique(temp))
  }
}


#' Directory of the installed package on your computer
#' 
#' A way of locating where the package was installed so you can find example
#' data files and code etc.
#' 
#' @usage DLMDataDir(stock=NA)
#' @param stock Character string representing the name of a .csv file e.g.
#' 'Snapper', 'Rockfish'
#' @author T. Carruthers
#' @examples
#' \dontrun{
#' tilefish_location <- DLMDataDir("Gulf_blue_tilefish")
#' tilefish_Data <- new("Data", tilefish_location)
#' }
#' @export DLMDataDir
DLMDataDir <- function(stock = NA) {
  if (is.na(stock)) {
    system.file(package = "DLMtool")
  } else {
    system.file(paste0(stock, ".csv"), package = "DLMtool", mustWork = TRUE)
  }
}


#' Load more data from DLMextra package
#'
#' Downloads the DLMextra package from GitHub 
#' @param silent Logical. Should messages to printed?
#' @export
#'
DLMextra <- function(silent=FALSE) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!silent) message("\nDownloading 'DLMextra' from GitHub")
  tt <- devtools::install_github("DLMtool/DLMextra", quiet=TRUE)
  if (tt) {
    if (!silent) message("Use 'library(DLMextra)' to load additional data into workspace")
  } else {
    if (!silent) message("Package 'DLMextra' already up to date\n Use 'library(DLMextra)' to load additional data into workspace")
  }
  
}

#' Convert a OM object to one without observation or process error
#' 
#' Takes an existing OM object and converts it to one without any observation
#' error, and very little process error.  Used for debugging and testing that
#' MPs perform as expected under perfect conditions.
#' 
#' 
#' @param OMin An object of class \code{OM}
#' @param except An optional vector of slot names in the OM that will not be
#' changed (not tested perfectly so watch out!)
#' @return A new \code{OM} object
#' @author A. Hordyk
#' @examples 
#' OM_noerror <- makePerf(DLMtool::testOM)
#' @export 
makePerf <- function(OMin, except = NULL) {
  nms <- slotNames(OMin)
  # exceptions
  if (is.null(except)) except <- "EVERYTHING"
  exclude <- unique(grep(paste(except, collapse = "|"), nms, value = FALSE))
  
  vars <- c("grad", "cv", "sd", "inc")
  ind <- unique(grep(paste(vars, collapse = "|"), nms, value = FALSE))
  ind <- ind[(!(nms[ind] %in% exclude))]
  for (X in seq_along(ind)) {
    n <- length(slot(OMin, nms[ind[X]]))
    if (n == 0) n <- 2
    slot(OMin, nms[ind[X]]) <- rep(0, n)
  }
  
  if (!("Cobs" %in% exclude)) 
    OMin@Cobs <- c(0, 0)
  if (!("Perr" %in% exclude)) 
    OMin@Perr <- c(0, 0)
  if (!("Iobs" %in% exclude)) 
    OMin@Iobs <- c(0, 0)
  if (!("AC" %in% exclude)) 
    OMin@AC <- c(0, 0)
  if (!("Btbiascv" %in% exclude)) 
    OMin@Btbiascv <- c(1, 1)
  if (!("CAA_ESS" %in% exclude)) 
    OMin@CAA_ESS <- c(1000, 1000)
  if (!("CAA_nsamp" %in% exclude)) 
    OMin@CAA_nsamp <- c(2000, 2000)
  if (!("CAL_ESS" %in% exclude)) 
    OMin@CAL_ESS <- c(1000, 1000)
  if (!("CAL_nsamp" %in% exclude)) 
    OMin@CAL_nsamp <- c(2000, 2000)
  if (!("beta" %in% exclude)) 
    OMin@beta <- c(1, 1)
  return(OMin)
}


#' Management Procedure Type
#'
#' @param MPs A vector of MP names. If none are provided function is run on all available MPs
#'
#' @return A data.frame with MP names, management type (e.g "Input", "Output") and management recommendations returned by the MP
#' (e.g, TAC (total allowable catch), TAE (total allowable effort), SL (size-selectivity), and/or or Spatial)
#' @export
#' @seealso \link{Required}
#' @examples 
#' MPtype(c("AvC", "curE", "matlenlim", "MRreal", "FMSYref"))
#' 
MPtype <- function(MPs=NA) {
  if (any(is.na(MPs))) MPs <- avail("MP")
  
  Data <- DLMtool::SimulatedData
  
  runMPs <- applyMP(Data, MPs, reps = 2, nsims=1, silent=TRUE)
  recs <- runMPs[[1]]
  
  type <- rep("NA", length(MPs))
  rec <- rep("", length(MPs))
  rectypes <- c("TAE", "Spatial", "SL")
  for (mm in seq_along(recs)) {
    Effort <- Spatial <- Selectivity <- FALSE
    output <- length(recs[[mm]]$TAC) > 0 
    names <- names(recs[[mm]])
    names <- names[!names %in% c("TAC", "Spatial")]
    input <- sum(unlist(lapply(Map(function(x) recs[[mm]][[x]], names), length))) > 0
    if (all(!is.na(recs[[mm]]$Spatial))) input <- TRUE
    if (output) {
      type[mm] <- "Output"
      thisrec <- "TAC"
    }
    if (input) {
      # what recommentations have been made?
      if (any(is.finite(recs[[mm]]$Effort))) Effort <- TRUE
      if (any(is.finite(recs[[mm]]$Spatial))) Spatial <- TRUE
      if (any(is.finite(recs[[mm]]$LR5)) | any(is.finite(recs[[mm]]$LFR)) | any(is.finite(recs[[mm]]$HS)) |
          any(is.finite(recs[[mm]]$Rmaxlen)) | any(is.finite(recs[[mm]]$L5)) | any(is.finite(recs[[mm]]$LFS)) |
          any(is.finite(recs[[mm]]$Vmaxlen))) Selectivity <- TRUE
      
      dorecs <- rectypes[c(Effort, Spatial, Selectivity)]
      thisrec <- dorecs
      type[mm] <- "Input"
      
    }
    if (input & output) {
      type[mm] <- "Mixed"
      thisrec <- c("TAC", thisrec)
    }
    if (length(thisrec)>1)  {
      rec[mm] <- paste(thisrec, collapse=", ")
    } else {
      rec[mm] <- thisrec
    }
  }
  type[grep("ref", MPs)] <- "Reference"
  
  df <- data.frame(MP=MPs, Type=type, Recs=rec, stringsAsFactors = FALSE)
  df[order(df$Type),]
  
}

#' Is a value NA or zero.
#' 
#' As title
#' 
#' 
#' @param x A numeric value.
#' @return TRUE or FALSE 
#' @author T. Carruthers
#' @keywords internal
#' @export
NAor0 <- function(x) {
  if (length(x) == 0) 
    return(TRUE)
  if (length(x) > 0) 
    return(is.na(x[1]))
}


#' Print out plotting functions
#' 
#' This function prints out the available plotting functions for objects of
#' class MSE or Data
#' 
#' @param class Character string. Prints out the plotting functions for objects
#' of this class.
#' @param msg Logical. Should the functions be printed to screen?
#' @note Basically the function looks for any functions in the DLMtool that
#' have the word `plot` in them.  There is a chance that some plotting
#' functions are missed. Let us know if you find any and we will add them.
#' @author A. Hordyk
#' @export 
plotFun <- function(class = c("MSE", "Data"), msg = TRUE) {
  class <- match.arg(class)
  tt <- lsf.str("package:DLMtool")
  p <- p2 <- rep(FALSE, length(tt))
  for (X in seq_along(tt)) {
    temp <- grep("plot", tolower(tt[[X]]))
    if (length(temp) > 0) p[X] <- TRUE
    temp2 <- grep(class, paste(format(match.fun(tt[[X]])), collapse = " "))
    if (length(temp2) > 0)  p2[X] <- TRUE
  }
  if (msg) 
    message("DLMtool functions for plotting objects of class ", class, 
            " are:")
  out <- sort(tt[which(p & p2)])
  
  out <- out[!out %in% c('plotFleet', 'plotStock', 'plotFun',
                         "COSEWIC_plot", "DFO_hist",
                         'plotOFL', "boxplot",
                         'boxplot.Data','plotDep','plotM2',
                         'plotGrowth', 'plotMat','plotRec', 'plot.OM')]
  
  if (class == "MSE") {
    out <- c(out, "barplot", "VOI", "VOI2", "DFO_proj",
             "PWhisker")
    out <- sort(out)
  }
  if (class == "Data") {
    out <- c(out, "boxplot", "Sense")
    out <- sort(out)
  }
  if (msg) {
    if (length(out) > 5) {
      sq <- seq(from = 1, to = length(out), by = 5)
      for (x in seq_along(sq)) {
        cat(out[sq[x]:min(sq[x + 1] - 1, length(out), na.rm = TRUE)])
        cat("\n")
      }
    } else cat(out)
    cat("\n")
  }
  invisible(out)
}


#' What management procedures need what data
#' 
#' A function that finds all the MPs and searches the
#' function text for slots in the Data object
#' 
#' @param funcs A character vector of management procedures
#' @param noCV Logical. Should the CV slots be left out? 
#' @author T. Carruthers
#' @return A matrix of MPs and their required data in terms of Data slotnames
#' @examples 
#' Required(c("DCAC", "AvC"))
#' Required() # For all MPs
#' @seealso \link{Can} \link{Cant} \link{Needed} \link{MPtype} \linkS4class{Data}
#' @export 
Required <- function(funcs = NA, noCV=FALSE) {
  
  if (class(funcs) != 'logical' & class(funcs) != "character") stop("first argument must be character with MP name")
  
  if (all(is.na(funcs))) funcs <- avail("MP")
  for (x in 1:length(funcs)) {
    tt <- try(get(funcs[x]))
    if (class(tt) != "MP") stop(funcs[x], " is not class 'MP'")
  } 
  
  
  builtin <- funcs[funcs %in% ReqData$MP]
  custom <- funcs[!funcs %in% ReqData$MP]
  
  df <- ReqData[match(builtin, ReqData$MP),]
  
  funcs1 <- custom
  
  temp <- lapply(funcs1, function(x) paste(format(match.fun(x)), collapse = " "))
  repp <- vapply(temp, match_slots, character(1))
  #repp[!nzchar(repp)] <- "No data needed for this MP."
  
  df2 <- data.frame(MP=funcs1, Data=repp, stringsAsFactors = FALSE)
  
  dfout <- rbind(df, df2)

  if (noCV) {
    for (rr in 1:nrow(dfout)) {
      tt <- unlist(strsplit(dfout$Data[rr], ","))
      tt <- tt[!grepl("CV_", tt)]
      tt <- sort(trimws(sort(tt)))
      dfout$Data[rr] <- paste(tt, collapse = ", ")
    }
  }
  
  as.matrix(dfout)
}



#' Get help topic URL
#'
#' @param topic Name of the functions
#' @param url URL for the help documentation
#' @param nameonly Logical. Help file name only?
#'
#' @return file path to help file
#' @export
#'
#' @keywords internal
MPurl <- function(topic, url='https://dlmtool.github.io/DLMtool/reference/',
                    nameonly=FALSE) {
  
  paths <- file.path(.libPaths()[1], "DLMtool")
  
  res <- character()
  for (p in paths) {
    if (file.exists(f <- file.path(p, "help", "aliases.rds"))) 
      al <- readRDS(f)
    else if (file.exists(f <- file.path(p, "help", "AnIndex"))) {
      foo <- scan(f, what = list(a = "", b = ""), sep = "\t", 
                  quote = "", na.strings = "", quiet = TRUE)
      al <- structure(foo$b, names = foo$a)
    }
    else next
    f <- al[topic]
    if (is.na(f)) 
      next
    res <- c(res, file.path(p, "help", f))
   
  }
  
  if(nameonly) {
    return(basename(res))
  } else{
    return(paste0(url, basename(res), ".html"))
  }
  
}


#' Setup parallel processing
#'
#' Sets up parallel processing using the snowfall package
#'
#' @param cpus number of CPUs 
#' @param ... other arguments passed to 'snowfall::sfInit'
#' @examples
#' \dontrun{
#' setup() # set-up half the available processors
#' setup(6) # set-up 6 processors
#' }
#' @export 
setup <- function(cpus=parallel::detectCores()*0.5, ...) {
  if(snowfall::sfIsRunning()) snowfall::sfStop()
  snowfall::sfInit(parallel=TRUE,cpus=cpus, ...)
  sfLibrary("DLMtool", character.only = TRUE, verbose=FALSE)
}



#' Open the DLMtool User Guide
#'
#' Opens the DLMtool User Guide website (requires internet connection)
#' 
#' @export
#' @examples
#' \dontrun{
#' userguide()
#' }
userguide <- function() {
  utils::browseURL("https://dlmtool.github.io/DLMtool/userguide/introduction.html")
}


RepmissingVal <- function(object, name, vals=NA) {
  miss <- FALSE
  if (!name %in% slotNames(object)) return(object)
  if (!.hasSlot(object,name)) miss <- TRUE
  if (!miss) {
    if (length(slot(object, name))==0) miss <- TRUE
    if (all(is.na(slot(object, name)))) miss <- TRUE 
  }
  if (miss) slot(object, name) <- vals
  return(object)
}

#' @describeIn checkMSE Updates an existing MSE object (class MSE) from a previous version of the
#' DLMtool to include slots new to the lastest version. Also works with Stock, 
#' Fleet, Obs, Imp, and Data objects. The new slots will be empty, 
#' but avoids the 'slot doesn't exist' error that sometimes occurs. 
#' Returns an object of class matching class(MSEobj)
#' @export 
updateMSE <- function(MSEobj) {
  slots <- slotNames(MSEobj)
  for (X in seq_along(slots)) {
    classDef <- getClassDef(class(MSEobj))
    slotTypes <- classDef@slots
    tt <- try(slot(MSEobj, slots[X]), silent = TRUE)
    if (class(tt) == "try-error") {
      fun <- get(as.character(slotTypes[X]))
      if(as.character(slotTypes[X]) == "vector") {
        slot(MSEobj, slots[X]) <- fun("numeric", length=0)
      } else slot(MSEobj, slots[X]) <- fun(0)
    }
  }
  MSEobj <- RepmissingVal(MSEobj, 'Mexp', c(0,0))
  MSEobj <- RepmissingVal(MSEobj, 'LenCV', c(0.08,0.15))
  MSEobj <- RepmissingVal(MSEobj, 'LR5', c(0,0))
  MSEobj <- RepmissingVal(MSEobj, 'LFR', c(0,0))
  MSEobj <- RepmissingVal(MSEobj, 'Rmaxlen', c(1,1))
  MSEobj <- RepmissingVal(MSEobj, 'DR', c(0,0))
  MSEobj <- RepmissingVal(MSEobj, 'Fdisc', c(0,0))
  MSEobj <- RepmissingVal(MSEobj, 'nareas', 2)
  
  MSEobj
}




#' Calculate CV from vector of values 
#' 
#' 
#' @param x vector of numeric values 
#' @author T. Carruthers
#' @return numeric
#' @keywords internal
#' @export
cv <- function(x) sd(x)/mean(x)


#' Get parameters of lognormal distribution from mean and standard deviation in normal
#' space
#' 
#' @param m mean in normal space 
#' @param sd standard deviation in normal space
#' @author T. Carruthers
#' @return numeric
#' @describeIn sdconv Returns sigma of lognormal distribution
#' @keywords internal
#' @export
sdconv <- function(m, sd) (log(1 + ((sd^2)/(m^2))))^0.5


#' @describeIn sdconv Returns mu of lognormal distribution
#' @export
mconv <- function(m, sd) log(m) - 0.5 * log(1 + ((sd^2)/(m^2)))

#' Calculate parameters for beta distribution from mean and standard deviation in
#' normal space
#' 
#' @param m mean 
#' @param sd standard deviation
#' @author T. Carruthers
#' @return numeric
#' @describeIn alphaconv Returns alpha of beta distribution
#' @keywords internal
#' @export
alphaconv <- function(m, sd) m * (((m * (1 - m))/(sd^2)) - 1)


#' @describeIn alphaconv Returns beta of beta distribution
#' @export 
betaconv <- function(m, sd) (1 - m) * (((m * (1 - m))/(sd^2)) - 1)

#' Lognormal distribution for DLMtool 
#' 
#' Variant of rlnorm which returns the mean when reps = 1.
#' 
#' @param reps number of random numbers 
#' @param mu mean 
#' @param cv coefficient of variation
#' @param x vector 
#' @author T. Carruthers
#' @return numeric
#' @describeIn trlnorm Generate log-normally distributed random numbers
#' @keywords internal 
#' @export 
trlnorm <- function(reps, mu, cv) {
  if (all(is.na(mu))) return(rep(NA, reps))
  if (all(is.na(cv))) return(rep(NA, reps))
  if (reps == 1)  return(mu)
  return(rlnorm(reps, mconv(mu, mu * cv), sdconv(mu, mu * cv)))
}



#' @describeIn trlnorm Calculate density of log-normally distributed random numbers 
#' @export 
tdlnorm <- function(x, mu, cv) dlnorm(x, mconv(mu, mu * cv), sdconv(mu, mu * cv))




#' Depletion and F estimation from mean length of catches
#' 
#' A highly dubious means of getting very uncertain estimates of current stock
#' biomass and (equilibrium) fishing mortality rate from growth, natural
#' mortality rate, recruitment and fishing selectivity.
#' 
#' @param OM An object of class 'OM'
#' @param ML A estimate of current mean length of catches
#' @param nsim Number of simulations
#' @param ploty Produce a plot of depletion and F
#' @param Dlim Limits on the depletion that is returned as a fraction of
#' unfished biomass.
#' @return An object of class 'OM' with 'D' slot populated 
#' @author T. Carruthers
#' @export 
ML2D <- function(OM, ML, nsim = 100, ploty = T, Dlim = c(0.05, 0.6)) {
  
  maxage <- OM@maxage
  M <- runif(nsim, OM@M[1], OM@M[2])  # Natural mortality rate
  h <- runif(nsim, OM@h[1], OM@h[2])  # Steepness
  Linf <- runif(nsim, OM@Linf[1], OM@Linf[2])  # Maximum length
  K <- runif(nsim, OM@K[1], OM@K[2])  # Maximum growth rate
  t0 <- runif(nsim, OM@t0[1], OM@t0[2])  # Theorectical length at age zero
  
  if (OM@isRel == "0" | OM@isRel == "FALSE" | OM@isRel == FALSE) {
    if (max(OM@LFS) > 0) {
      LFS <- runif(nsim, OM@LFS[1], OM@LFS[2])
    } else {
      LFS <- runif(nsim, mean(OM@LFSLower), mean(OM@LFSUpper))
    }
  } else {
    if (max(OM@LFS) > 0) {
      LFS <- runif(nsim, OM@LFS[1], OM@LFS[2]) * mean(OM@L50)
    } else {
      LFS <- runif(nsim, mean(OM@LFSLower), mean(OM@LFSUpper)) * 
        mean(OM@L50)
    }
  }
  AFS <- L2A(t0, Linf, K, LFS, maxage)
  
  L5 <- runif(nsim, OM@L5[1], OM@L5[2]) * mean(OM@L50)
  age05 <- L2A(t0, Linf, K, L5, maxage)
  age05[age05<0] <- 0
  
  Vmaxage <- runif(nsim, OM@Vmaxlen[1], OM@Vmaxlen[2])  #runif(BT_fleet@Vmaxage[1],BT_fleet@Vmaxage[2]) # selectivity of oldest age class
  
  LM <- runif(nsim, OM@L50[1], OM@L50[2])
  AM <- L2A(t0, Linf, K, LM, maxage)
  
  # age at maturity
  a <- OM@a  # length-weight parameter a
  b <- OM@b  # length-weight parameter b
  
  mod <- AFS  # the age at modal (or youngest max) selectivity
  
  # deriv <- getDNvulnS(mod, age05, Vmaxage, maxage, nsim)  # The vulnerability schedule
  # vuln <- deriv[[1]]
  
  srs <- (maxage - AFS) / ((-log(Vmaxage,2))^0.5) # selectivity parameters are constant for all years
  sls <- (AFS - age05) /((-log(0.05,2))^0.5)
  
  vuln <- t(sapply(1:nsim, getsel, lens=matrix(1:maxage, nrow=nsim, ncol=maxage, byrow=TRUE), 
                   lfs=AFS, sls=sls, srs=srs))
  
  Agearray <- array(rep(1:maxage, each = nsim), c(nsim, maxage))
  mat <- 1/(1 + exp((AM - (Agearray))/(AM * 0.1)))  # Maturity at age array
  
  nyears <- 100
  # bootfun<-function(dat,ind)mean(dat[ind]) MLo<-boot(MLt,bootfun,nsim)
  # ML<-MLo$t
  out <- CSRA(M, h, Linf, K, t0, AM, a, b, vuln, mat, ML = rep(ML, nsim), 
              NA, NA, maxage, nyears)
  cond <- out[, 1] > Dlim[1] & out[, 1] < Dlim[2] & out[, 2] < 2.5  # Stock levels are unlikely to be above 80% unfished, F is unlikely to be above 2.5
  
  if (ploty & sum(cond) > 0) {
    par(mfrow = c(1, 2))
    plot(density(out[cond, 1], from = 0, adj = 0.4), main = "Depletion")
    plot(density(out[cond, 2], from = 0, adj = 0.4), main = "Fishing mortality rate")
    OM@D <- quantile(out[cond, 1], c(0.05, 0.95))
  }
  if (sum(cond) == 0) {
    message("All estimates of Depletion outside bounds of Dlim")
    message("Operating Model object not updated")
  }
  if (!ploty) 
    message("Operating Model object not updated")
  OM
}

# Composition stock reduction analysis


#' Catch at size reduction analysis
#' 
#' What depletion level and corresponding equlibrium F arise from data
#' regarding mean length of current catches, natural mortality rate, steepness
#' of the stock recruitment curve, maximum length, maximum growth rate, age at
#' maturity, age based vulnerability, maturity at age, maximum age and number
#' of historical years of fishing.
#' 
#' 
#' @usage CSRA(M,h,Linf,K,t0,AM,a,b,vuln,mat,ML,CAL,CAA,maxage,nyears)
#' @param M A vector of natural mortality rate estimates
#' @param h A vector of sampled steepness (Beverton-Holt stock recruitment)
#' @param Linf A vector of maximum length (von Bertalanffy growth)
#' @param K A vector of maximum growth rate (von Bertalanffy growth)
#' @param t0 A vector of theoretical age at length zero (von Bertalanffy
#' growth)
#' @param AM A vector of age at maturity
#' @param a Length-weight conversion parameter a (W=aL^b)
#' @param b Length-weight conversion parameter b (W=aL^b)
#' @param vuln A matrix nsim x nage of the vulnerabilty at age (max 1) to
#' fishing.
#' @param mat A matrix nsim x nage of the maturity at age (max 1)
#' @param ML A vector of current mean length estimates
#' @param CAL A catch-at-length matrix nyears x (1 Linf unit) length bins
#' @param CAA A catch-at-age matrix nyears x maximum age
#' @param maxage Maximum age
#' @param nyears Number of historical years of fishing
#' @author T. Carruthers
#' @export CSRA
#' @keywords internal
CSRA <- function(M, h, Linf, K, t0, AM, a, b, vuln, mat, ML, CAL, CAA, 
                 maxage, nyears) {
  nsim <- length(M)
  Dep <- rep(NA, nsim)
  Fm <- rep(NA, nsim)
  for (i in 1:nsim) {
    fit <- optimize(CSRAfunc, log(c(1e-04, 5)), Mc = M[i], hc = h[i], 
                    maxage, nyears, Linfc = Linf[i], Kc = K[i], t0c = t0[i], AMc = AM[i], 
                    ac = a, bc = b, vulnc = vuln[i, ], matc = mat[i, ], MLc = ML[i], 
                    CAL = NA, CAA = NA, opt = T)
    
    
    out <- CSRAfunc(fit$minimum, Mc = M[i], hc = h[i], maxage, nyears, 
                    Linfc = Linf[i], Kc = K[i], t0c = t0[i], AMc = AM[i], ac = a, 
                    bc = b, vulnc = vuln[i, ], matc = mat[i, ], MLc = ML[i], CAL = NA, 
                    CAA = NA, opt = 3)
    
    Dep[i] <- out[1]
    Fm[i] <- out[2]
    
    
  }
  cbind(Dep, Fm)
}

# The function that CSRA operates on

#' Optimization function for CSRA
#' 
#' What depletion level and corresponding equlibrium F arise from data
#' regarding mean length of current catches, natural mortality rate, steepness
#' of the stock recruitment curve, maximum length, maximum growth rate, age at
#' maturity, age based vulnerability, maturity at age, maximum age and number
#' of historical years of fishing.
#' 
#' 
#' @param lnF A proposed value of current instantaneous fishing mortality rate
#' @param Mc Natural mortality rate estimates
#' @param hc Steepness (Beverton-Holt stock recruitment)
#' @param maxage Maximum age
#' @param nyears Number of historical years of fishing
#' @param AFSc Age at full selection
#' @param AFCc Age at first capture
#' @param Linfc Maximum length (von Bertalanffy growth)
#' @param Kc Maximum growth rate (von Bertalanffy growth)
#' @param t0c Theoretical age at length zero (von Bertalanffy growth)
#' @param AMc Age at maturity
#' @param ac Length-weight conversion parameter a (W=aL^b)
#' @param bc Length-weight conversion parameter b (W=aL^b)
#' @param vulnc A vector (nage long) of the vulnerabilty at age (max 1) to
#' fishing.
#' @param matc A vector (nage long) of the maturity at age (max 1)
#' @param MLc A current mean length estimates
#' @param CAL A catch-at-length matrix nyears x (1 Linf unit) length bins
#' @param CAA A catch-at-age matrix nyears x maximum age
#' @param opt Should the measure of fit be returned?
#' @param meth Are we fitting to mean length or catch composition?
#' @author T. Carruthers
#' @keywords internal
CSRAfunc <- function(lnF, Mc, hc, maxage, nyears, AFSc, AFCc, Linfc, Kc, 
                     t0c, AMc, ac, bc, vulnc, matc, MLc, CAL, CAA, opt = T, meth = "ML") {
  
  Fm <- exp(lnF)
  Fc <- vulnc * Fm
  Lac <- Linfc * (1 - exp(-Kc * ((1:maxage) - t0c)))
  Wac <- ac * Lac^bc
  N <- exp(-Mc * ((1:maxage) - 1))
  SSN <- matc * N  # Calculate initial spawning stock numbers
  Biomass <- N * Wac
  SSB <- SSN * Wac  # Calculate spawning stock biomass
  
  B0 <- sum(Biomass)
  SSB0 <- sum(SSB)
  SSN0 <- SSN
  SSBpR <- sum(SSB)  # Calculate spawning stock biomass per recruit
  SSNpR <- SSN
  Zc <- Fc + Mc
  CN <- array(NA, dim = c(nyears, maxage))
  HR <- rep(0, maxage)
  pen <- 0
  for (y in 1:nyears) {
    VB <- Biomass * vulnc * exp(-Mc)
    CN[y, ] <- N * (1 - exp(-Zc)) * (Fc/Zc)
    N[2:maxage] <- N[1:(maxage - 1)] * exp(-Zc[1:(maxage - 1)])  # Total mortality
    N[1] <- (0.8 * hc * sum(SSB))/(0.2 * SSBpR * (1 - hc) + (hc - 0.2) * 
                                     sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    Biomass <- N * Wac
    SSN <- N * matc
    SSB <- SSN * Wac
  }  # end of year
  
  pred <- sum((CN[nyears, ] * Lac))/sum(CN[nyears, ])
  fobj <- (pred - MLc)^2  # Currently a least squares estimator. Probably not worth splitting hairs WRT likelihood functions!
  if (opt == 1) {
    return(fobj)
  } else {
    c(sum(SSB)/sum(SSB0), Fm)
  }
}

# Stochastic inverse growth curve used to back-calculate age at first
# capture from length at first capture
#' Calculate age at first capture from length at first capture and growth
#' 
#' As title.
#' 
#' @param t0c A vector of theoretical age at length zero (von Bertalanffy
#' growth)
#' @param Linfc A vector of maximum length (von Bertalanffy growth)
#' @param Kc A vector of maximum growth rate (von Bertalanffy growth)
#' @param LFC A vector of length at first capture
#' @param maxage Maximum age
#' @author T. Carruthers
#' @keywords internal
getAFC <- function(t0c, Linfc, Kc, LFC, maxage) {
  nsim <- length(t0c)
  agev <- c(1e-04, 1:maxage)
  agearray <- matrix(rep(agev, each = nsim), nrow = nsim)
  Larray <- Linfc * (1 - exp(-Kc * (agearray - t0c)))
  matplot(agev, t(Larray), type = "l")
  abline(h = LFC, col = "#ff000030", lwd = 2)
  AFC <- (log(1 - (LFC/Linfc))/-Kc) + t0c
  abline(v = AFC, col = "#0000ff30", lwd = 2)
  AFC
}



#' Length to age conversion
#' 
#' Simple deterministic length to age conversion given inverse von Bertalanffy
#' growth.
#' 
#' @param t0c Theoretical age at length zero
#' @param Linfc Maximum length
#' @param Kc Maximum growth rate
#' @param Len Length
#' @param maxage Maximum age
#' @return An age (vector of ages, matrix of ages) corresponding with Len
#' @author T. Carruthers
#' @keywords internal
L2A <- function(t0c, Linfc, Kc, Len, maxage) {
  nsim <- length(t0c)
  agev <- c(1e-04, 1:maxage)
  agearray <- matrix(rep(agev, each = nsim), nrow = nsim)
  Larray <- Linfc * (1 - exp(-Kc * (agearray - t0c)))
  matplot(agev, t(Larray), type = "l")
  abline(h = Len, col = "#ff000030", lwd = 2)
  age <- (log(1 - (Len/Linfc))/-Kc) + t0c
  abline(v = age, col = "#0000ff30", lwd = 2)
  age
}



#' Determine optimal number of cpus
#'
#' @param nsim Numeric. Number of simulations.
#' @param thresh Recommended n cpus is what percent of the fastest time?
#' @param plot Logical. Show the plot?
#' @param msg Logical. Should messages be printed to console?
#' @param maxn Optional. Maximum number of cpus. Used for demo purposes
#'
#' @export
#' @seealso \link{setup}
#' @examples
#' \dontrun{
#' optCPU()
#' }
#' @author A. Hordyk
optCPU <- function(nsim=96, thresh=5, plot=TRUE, msg=TRUE, maxn=NULL) {
  cpus <- 1:parallel::detectCores()
  if (!is.null(maxn)) cpus <- 1:maxn
  
  time <- NA
  OM <- DLMtool::testOM
  OM@nsim <- nsim
  for (n in cpus) {
    if (msg) message('Running MSE with ', nsim, ' simulations and ', n, ' of ', max(cpus), ' cpus')
    if (n == 1) {
      snowfall::sfStop()
      st <- Sys.time()
      tt <- runMSE(OM, silent = TRUE)
      time[n] <- difftime(Sys.time(), st, units='secs')
    } else{
      if (msg) {
        setup(cpus=n)
      } else {
        sink('temp')
        suppressMessages(setup(cpus=n))
        sink()
      }
      st <- Sys.time()
      tt <- runMSE(OM, silent=TRUE, parallel=TRUE)  
      
      time[n] <- difftime(Sys.time(), st, units='secs')
      
    }
  } 
  df <- data.frame(ncpu=cpus, time=time)
  df$time <- round(df$time,2)
  rec <- min(which(time < min(time) * (1 + thresh/100)))
  if (plot) {
    plot(df, type='b', ylab="time (seconds)", xlab= "# cpus", bty="l", lwd=2)
    points(rec, df[rec,2], cex=2, pch=16, col="blue")
  }
  return(df)
}

#' DLMenv blank environment
#' 
#' An environment allocated for MPs to print model output during the
#' management strategy evaluation. Is blank at the beginning of each call to \code{runMSE}.
#' 
#' @seealso \link{runMSE}
#' @export
DLMenv <- new.env()


