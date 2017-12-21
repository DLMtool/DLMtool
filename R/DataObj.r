# A collection of functions which operate on DLM data object (class
# Data) which is used as input object for the management
# procedures.

# These functions are used either internally in the MSE (runMSE.r) or
# used manually to apply a MP (or MPs) to a particular data object.

DLMdiag <- function(Data, command = "available", reps = 5, timelimit = 1, funcs1=NA) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  if (all(is.na(funcs1))) funcs1 <- avail("MP")
  # mpclasses <- MPtype(funcs1)[,2]
  # funcs1 <- funcs1[mpclasses %in% c("Output", "Input", "Mixed", "Reference")]
  good <- rep(TRUE, length(funcs1))
  report <- rep("Worked fine", length(funcs1))
  test <- new("list")
  timey <- new("list")
  on.exit(options(show.error.messages = TRUE))
  for (y in 1:length(funcs1)) {
    # First check for required data slots that are not NA
    time1 <- Sys.time() # initialise time1
    temp <- format(match.fun(funcs1[y])) # Extract the dependencies from function 
    temp <- paste(temp[1:(length(temp))], collapse = " ")
    
    slots <- slotNames(Data)
    rr <- try(slot(Data, "Misc"), silent = TRUE)
    if (class(rr) == "try-error") Data@Misc <- list()
    slotnams <- paste("Data@", slotNames(Data), sep = "")
    # slotnams <- slotnams[!grepl("OM", slotnams)]
    
    chk <- rep(FALSE, length(slotnams))
    for (j in 1:length(slotnams)) {
      if (grepl(slotnams[j], temp) & all(is.na(slot(Data, slots[j])))) chk[j] <- TRUE
    }
    
    if (sum(chk) > 0) {
      test[[y]] <- "missing data" # if some required slots are NA or NULL - return error
      class(test[[y]]) <- "try-error"
    } else {  # continue if all slots are valid 
      time1 <- Sys.time()
      suppressWarnings({
        setTimeLimit(timelimit * 1.5)
        options(show.error.messages = FALSE)
        test[[y]] <- try(do.call(funcs1[y], list(x = 1, Data = Data, reps = 5)), silent = T)
        options(show.error.messages = TRUE)
        if (class(test[[y]]) == "list") test[[y]] <- test[[y]][[1]]
        setTimeLimit(Inf)
      })
      
      # if (class(match.fun(funcs1[y])) == "Output") {
      #   time1 <- Sys.time()
      #   suppressWarnings({
      #     setTimeLimit(timelimit * 1.5)
      #     options(show.error.messages = FALSE)
      #     test[[y]] <- try(do.call(funcs1[y], list(x = 1, Data = Data, reps = 5)), silent = T)
      #     options(show.error.messages = TRUE)
      #     if (class(test[[y]]) == "list") test[[y]] <- test[[y]][[1]]
      #     setTimeLimit(Inf)
      #   })
      # } else {
      #   time1 <- Sys.time()
      #   suppressWarnings({
      #     setTimeLimit(timelimit * 1.5)
      #     options(show.error.messages = FALSE)
      #     test[[y]] <- try(do.call(funcs1[y], list(x = 1, Data = Data)), silent = T)
      #     options(show.error.messages = TRUE)
      #     setTimeLimit(Inf)
      #   })
      # }      
    }
    time2 <- Sys.time()
    timey[[y]] <- time2 - time1
    if (class(test[[y]]) == "try-error") {
      report[[y]] <- "Insufficient data"
      good[[y]] <- FALSE
    } else if (class(test[[y]]) == "Rec") {
      slts <- slotNames(test[[y]])
      tt <- rep(FALSE, length(slts))
      for (x in seq_along(slts)) tt[x] <- NAor0(slot(test[[y]], slts[x]))
      if (sum(tt) < 1) {
        report[[y]] <- "Produced all NA scores"
        good[[y]] <- FALSE
      }  
    } else if (sum(is.na(test[[y]])) == length(test[[y]])) {
      report[[y]] <- "Produced all NA scores"
      good[[y]] <- FALSE
    }
    if (timey[[y]] > timelimit) {
      report[[y]] <- "Exceeded the user-specified time limit"
      good[[y]] <- FALSE
    }
  }  # end of funcs
  options(show.error.messages = TRUE)
  if (command == "available") return(funcs1[good])
  if (command == "not available") return(cbind(funcs1[!good], report[!good]))
  if (command == "needed")  return(needed(Data, funcs = funcs1[!good]))
}

needed <- function(Data, funcs = NA) {
  if (is.na(funcs[1]))  funcs <- avail("Output")
  slots <- slotNames("Data")
  rr <- try(slot(Data, "Misc"), silent = TRUE)
  if (class(rr) == "try-error") Data@Misc <- list()
  
  slotnams <- paste("Data@", slotNames(Data), sep = "")
  repp <- rep("", length(funcs))
  Data@Misc <- list()
  for (i in 1:length(funcs)) {
    temp <- format(match.fun(funcs[i]))
    temp <- paste(temp[1:(length(temp))], collapse = " ")
    rec <- ""
    for (j in 1:length(slotnams)) {
      if (grepl(slotnams[j], temp) & NAor0(slot(Data, slots[j]))) rec <- c(rec, slots[j])
    }
    repp[i] <- paste(funcs[i], ": ", paste(rec[2:length(rec)], collapse = ", "), sep = "")
    if (length (rec) == 1) 
      repp[i] <- paste0(funcs[i], ": All data available. Model returns NA")
  }
  repp
}

parallelMPs <- function(x, Data, reps, MPs, ss) sapply(ss, MPs[x], 
  Data, reps = reps)

#' Make stochastic variables certain for only one rep
#' 
#' As title.
#'
#' @param Data An object of class Data that has been run though TAC()
#' @author T. Carruthers
#' @keywords internal 
#' @export OneRep
OneRep <- function(Data) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  Data@CV_Cat = Data@CV_Dt = Data@CV_AvC = Data@CV_Ind = Data@CV_Mort = Data@CV_FMSY_M = Data@CV_BMSY_B0 = Data@CV_Cref = Data@CV_Bref = Data@CV_Iref = Data@CV_Rec = Data@CV_Dep = Data@CV_Abun = Data@CV_L50 = Data@CV_vbK = Data@CV_vbLinf = Data@CV_vbt0 = Data@CV_LFC = Data@CV_LFS = Data@CV_wla = Data@CV_wlb = Data@CV_steep = Data@sigmaL = tiny
  Data
}


#' A generic OFL plot for NOAA use
#' 
#' As title.
#' 
#' 
#' @usage plotOFL(Data,xlims=NA,perc=0.5)
#' @param Data An object of class Data that has been run though TAC()
#' @param xlims x axis limits
#' @param perc The percentile of the OFL distribution to be plotted
#' @return A table of performance metrics.
#' @author T. Carruthers
#' @export plotOFL
plotOFL <- function(Data, xlims = NA, perc = 0.5) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  cols <- rep(c("black", "red", "green", "blue", "orange", "brown", "purple", 
    "dark grey", "violet", "dark red", "pink", "dark blue", "grey"), 
    4)
  ltys <- rep(1:4, each = 13)
  
  funcs <- Data@MPs
  nMPs <- length(funcs)
  
  if (is.na(xlims[1]) | length(xlims) != 2) {
    xlims <- quantile(Data@TAC, c(0.005, 0.9), na.rm = T)
    if (xlims[1] < 0) 
      xlims[1] <- 0
  }
  if (!NAor0(Data@Ref)) {
    if (xlims[1] > Data@Ref) 
      xlims[1] <- max(0, 0.98 * Data@Ref)
    if (xlims[2] < Data@Ref) 
      xlims[2] <- 1.02 * Data@Ref
    if (xlims[2] > Data@Ref * 2) 
      xlims[2] <- 2 * Data@Ref
  }
  ylims <- c(0, 1)
  
  plot(NA, NA, xlim = xlims, ylim = ylims, main = "", xlab = "", ylab = "", 
    col = "white", lwd = 3, type = "l")
  abline(h = 0)
  if (!NAor0(Data@Ref)) {
    abline(v = Data@Ref, col = "light grey", lwd = 3)
    if (!NAor0(Data@Ref_type[1])) 
      legend("bottomright", Data@Ref_type, text.col = "grey", 
        bty = "n")
  }
  
  for (m in 1:nMPs) {
    
    if (sum(!is.na(Data@TAC[m, , 1])) > 10) {
      # only plot if there are sufficient non-NA TAC samples
      x <- density(Data@TAC[m, , 1], from = 0, na.rm = T)$x
      y <- density(Data@TAC[m, , 1], from = 0, na.rm = T)$y
      y <- y/max(y)
      lines(x, y, col = cols[m])
    } else {
      print(paste("Method ", funcs[m], " produced too many NA TAC values for plotting densities", 
        sep = ""))
    }
    if (!is.na(perc[1])) 
      abline(v = quantile(Data@TAC[m, , 1], p = perc, na.rm = T), 
        col = cols[m], lty = 2)
  }
  
  cind <- 1:nMPs
  legend("topright", funcs, text.col = cols[cind], col = cols[cind], 
    lty = 1, bty = "n", cex = 0.75)
  
  mtext(paste("OFL (", Data@Units, ")", sep = ""), 1, outer = F, 
    line = 2.6)
  mtext(paste("Standardized relative frequency", sep = ""), 2, outer = F, 
    line = 2.6)
  # mtext(paste('OFL calculation for
  # ',Data@Name,sep=''),3,outer=F,line=1)
  
}

# Primary functions


#' What data-limited methods can be applied to this Data object?
#' 
#' An diagnostic tool that looks up the slot requirements of each method and
#' compares this to the data available to limit the analysis to methods that
#' have the correct data, do not produce errors and run within a time limit.
#' Time limit is the maximum time taken to carry out five reps (stochastic
#' samples) of a given method and is in units of seconds.
#' 
#' 
#' @usage Can(Data, timelimit = 1)
#' @param Data A data-limited methods data object (class Data)
#' @param timelimit The maximum time (seconds) taken for a method to undertake
#' 10 reps (this filters out methods that are too slow)
#' @param MPs Optional list of MP names
#' @export Can
Can <- function(Data, timelimit = 1, MPs=NA) {
  DLMdiag(Data, "available",  timelimit = timelimit, funcs1=MPs)
}


#' What methods can't be applied to this DLM data object
#' 
#' The methods that don't have sufficient data, lead to errors or don't run in
#' time along with a list of their data requirments.
#' 
#' 
#' @usage Cant(Data, timelimit = 1)
#' @param Data A data-limited methods data object (class Data)
#' @param timelimit The maximum time (seconds) taken for a method to undertake
#' 10 reps (this filters out methods that are too slow)
#' @export Cant
Cant <- function(Data, timelimit = 1) {
  DLMdiag(Data, "not available", timelimit = timelimit)
}

#' Data needed to get MPs running
#' 
#' Wrapper function for DLMdiag that lists what data are needed to run
#' data-limited methods that are current not able to run given a DLM_cdata
#' object
#' 
#' 
#' @usage Needed(Data, timelimit=1)
#' @param Data A data-limited methods data object
#' @param timelimit The maximum time (seconds) taken to complete 10 reps
#' @author T. Carruthers
#' @export Needed
Needed <- function(Data, timelimit = 1) {
  DLMdiag(Data, "needed", timelimit = timelimit)
}

# A function that determines the inputs for a given data-limited method
# of class Output and then analyses the sensitivity of TAC
# estimates to marginal differences in each input.



#' Sensitivity analysis
#' 
#' A function that determines the inputs for a given data-limited method of
#' class Output and then analyses the sensitivity of TAC estimates to 
#' marginal differences in each input. The range used for sensitivity is based 
#' on the user-specified CV for that input (e.g. CV_Mort, Mort)
#' 
#' 
#' @usage Sense(Data, MP, nsense = 6, reps = 100, perc = c(0.05, 0.5, 0.95), ploty = T)
#' @param Data A data-limited methods data object
#' @param MP A character string representing an MP applied in calculating the TAC recommendations in the DLM object
#' @param nsense The number of points over which to calculate the TAC (resolution)
#' @param reps The number of samples of the quota taken for the calculation of the TAC
#' @param perc The percentile of the sample TAC
#' @param ploty A logical switch, (T/F, should a plot be drawn?)
#'
#' @author T. Carruthers
#' @export Sense
Sense <- function(Data, MP, nsense = 6, reps = 100, perc = c(0.05, 0.5, 0.95), 
                  ploty = T) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  
  DLM_data2 <- Data
  nm <- deparse(substitute(DLM_data2))
  refTAC <- quantile(getTAC(DLM_data2, MP, reps)[[1]], perc, na.rm = T)
  
  Data <- DLM_data2
  reqs <- Required(MP)  #read.csv(paste(getwd(),'/Data/Data requirements.csv',sep=''),header=T)
  ind <- (1:nrow(reqs))[reqs[, match(MP, names(reqs))] == "Y"]
  # for(i in 1:length(reqs))
  
  
  slotsCV <- slotNames("Data")[grep("CV_", slotNames("Data"))]
  slots <- rep("", length(slotsCV))
  for (i in 1:length(slotsCV)) slots[i] <- substr(slotsCV[i], 4, nchar(slotsCV[i]))
  
  ind <- slots %in% unlist(strsplit(reqs[2], ", "))
  slots <- slots[ind]
  slotsCV <- slotsCV[ind]
  sname <- slots
  nslots <- length(slots)
  
  nrep <- nslots * nsense
  Data <- replic8(Data, nrep)
  pss <- seq(0, 1, length.out = nsense + 2)[2:(nsense + 1)]
  vals <- array(NA, dim = c(nslots, nsense))
  
  for (i in 1:nslots) {
    ind <- (((i - 1) * nsense + 1):(i * nsense))
    mn <- attr(Data, slots[i])[1]
    cv <- attr(Data, slotsCV[i])[1] * 2  # twice the CV of the variable specified in the DLM object
    if (class(attr(Data, slots[i])) == "numeric") {
      if (mn > 0) {
        attr(Data, slots[i])[ind] <- qlnorm(pss, mconv(mn, 
          cv * mn), sdconv(mn, cv * mn))
        vals[i, ] <- qlnorm(pss, mconv(mn, cv * mn), sdconv(mn, 
          cv * mn))
      } else {
        attr(Data, slots[i])[ind] <- -qlnorm(pss, mconv(-mn, 
          cv * -mn), sdconv(-mn, cv * -mn))
        vals[i, ] <- -qlnorm(pss, mconv(-mn, cv * -mn), sdconv(-mn, 
          cv * -mn))
      }
    } else {
      cv <- attr(Data, slotsCV[i])[1]
      attr(Data, slots[i])[ind, ] <- attr(Data, slots[i])[ind, 
        ] * qlnorm(pss, mconv(1, cv), sdconv(1, cv))
      vals[i, ] <- qlnorm(pss, mconv(1, cv), sdconv(1, cv))
    }
  }
  
  TACa <- getTAC(Data, MPs = MP, reps = reps)[[1]]
  TACa <- apply(TACa, 3, quantile, p = perc, na.rm = T)
  LB <- ((1:nslots) - 1) * 4 + 1
  UB <- (1:nslots) * 4
  sense <- matrix(data = NA, nrow = 4 * nslots, ncol = nsense + 1)
  
  for (i in 1:nslots) {
    ind <- ((i - 1) * nsense + 1):(i * nsense)
    dat <- TACa[, ind]
    
    sense[LB[i], 2:(nsense + 1)] <- vals[i, ]
    sense[(LB[i] + 1):UB[i], 2:(nsense + 1)] <- dat
    sense[LB[i], 1] <- slots[i]
    sense[(LB[i] + 1):UB[i], 1] <- perc
  }
  
  DLM_data2@Sense <- sense
  
  if (ploty) {
    ylimy <- range(TACa)
    # dev.new2(width=10,height=0.5+3*ceiling(nslots/2))
    par(mfrow = c(ceiling(nslots/2), 2), mai = c(0.4, 0.4, 0.01, 0.01), 
      omi = c(0.4, 0.4, 0.4, 0.01))
    for (i in 1:nslots) {
      ind <- (((i - 1) * nsense + 1):(i * nsense))
      dat <- TACa[, ind]
      xlimy <- range(vals[i, ])
      plot(xlimy, rep(refTAC[2], 2), ylim = ylimy, xlim = xlimy, 
        type = "l", col = "#99999960", main = "", xlab = "", ylab = "")
      abline(h = refTAC[c(1, 3)], col = "#99999960", lty = 2)
      abline(v = slot(DLM_data2, slots[i]), col = "#99999960", lty = 2)
      lines(vals[i, ], dat[2, ], col = "red", lwd = 1.5)
      lines(vals[i, ], dat[1, ], col = "red", lty = 2, lwd = 1.5)
      lines(vals[i, ], dat[3, ], col = "red", lty = 2, lwd = 1.5)
      legend("top", legend = sname[i], text.col = "blue", bty = "n")
    }
    
    mtext(paste("Output control (", Data@Units, ")", sep = ""), 
      2, outer = T, line = 0.5)
    mtext("Parameter / variable input level", 1, outer = T, line = 0.5)
    mtext(paste("Sensitivity analysis for ", Data@Name, ": ", MP, 
      sep = ""), 3, outer = T, line = 0.5)
  }
  # assign(nm,DLM2,envir=.GlobalEnv)
  DLM_data2
}

# Replicates position 1 data to multiple positions for sensitivity
# testing etc


#' Enlarge (replicate) a DLM data object to create an additional dimension for
#' simulation / sensitivity testing
#' 
#' Replicates position 1 data to multiple positions for sensitivity testing etc
#' 
#' 
#' @usage replic8(Data, nrep)
#' @param Data A data-limited methods data object
#' @param nrep The number of positions to expand the DLM object to
#' @author T. Carruthers
#' @export replic8
replic8 <- function(Data, nrep) {
  
  slotnam <- slotNames(Data)
  slotnam <- slotnam[slotnam != "Ref" & slotnam != "OM" & slotnam != 
    "MaxAge" & slotnam != "CAL_bins" & slotnam != "Year"]
  
  for (sl in 1:length(slotnam)) {
    slt <- attr(Data, slotnam[sl])
    if (class(slt) == "matrix") {
      attr(Data, slotnam[sl]) <- matrix(rep(slt, each = nrep), 
        nrow = nrep, ncol = ncol(slt))
    } else if (class(slt) == "numeric") {
      attr(Data, slotnam[sl]) <- rep(slt, nrep)
    } else if (class(slt) == "array") {
      attr(Data, slotnam[sl]) <- array(rep(slt, each = nrep), 
        dim = c(nrep, dim(slt)[2:3]))
    }
  }
  Data
}

# A function that returns the stochastic TAC recommendations from a
# vector of data-limited MPs (Output) given a data-limited data
# object Data


#' Calculate TAC recommendations for more than one MP
#' 
#' A function that returns the stochastic TAC recommendations from a vector of
#' data-limited MPs (Output) given a data-limited data object Data
#' 
#' 
#' @usage TAC(Data, MPs = NA, reps = 100, timelimit = 1)
#' @param Data A data-limited methods data object
#' @param MPs optional vector of MP names
#' @param reps Number of repititions
#' @param timelimit The maximum time (seconds) taken to complete 10 reps
#' @author T. Carruthers
#' @export TAC
TAC <- function(Data, MPs = NA, reps = 100, timelimit = 1) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  nm <- deparse(substitute(Data))
  PosMPs <- Can(Data, timelimit = timelimit)
  PosMPs <- PosMPs[PosMPs %in% avail("Output")]
  Data@PosMPs <- PosMPs
  if (!is.na(MPs[1])) 
    Data@MPs <- MPs[MPs %in% PosMPs]
  if (is.na(MPs[1])) 
    Data@MPs <- PosMPs
  funcs <- Data@MPs
  
  if (length(funcs) == 0) {
    stop("None of the methods 'MPs' are possible given the data available")
  } else {
    temp <- getTAC(Data, MPs = funcs, reps)
    TACa <- temp[[1]]
    Data <- temp[[2]]
    Data@TAC <- TACa
    return(Data)
    # assign(nm,DLM,envir=.GlobalEnv)
  }
  
}

getTAC <- function(Data, MPs = NA, reps = 100) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  nsims <- length(Data@Mort)
  nMPs <- length(MPs)
  TACa <- array(NA, dim = c(nMPs, reps, nsims))
  
  if (!sfIsRunning()) {
    for (ff in 1:nMPs) {
      temp <- sapply(1:nsims, MPs[ff], Data = Data, reps = reps)
      if (mode(temp) == "numeric") 
        TACa[ff, , ] <- temp
      if (mode(temp) == "list") {
        TACa[ff, , ] <- unlist(temp[1, ])
        for (x in 1:nsims) Data@Misc[[x]] <- temp[2, x][[1]]
      }
    }
  } else {
    sfExport("Data") 
    if (nsims < 8) {
      sfExport(list = c("MPs", "reps"))
      for (ss in 1:nsims) {
        temp <- (sfSapply(1:length(MPs), parallelMPs, Data = Data, 
          reps = reps, MPs = MPs, ss = ss))
        if (mode(temp) == "list") {
          Lens <- unlist(lapply(temp, length))
          for (X in 1:length(Lens)) {
          Classes <- unlist(lapply(temp[, X][[1]], class))
          if (length(unique(Classes)) == 1) {
            # no Misc object
            TACa[X, , ss] <- unlist(temp[, X])
          } else {
            # a Misc object is include
            ind <- which(Classes == "list")
            TACa[X, , ss] <- unlist(temp[, X][[1]][1:(ind - 1), 
            ])
            Data@Misc[[ss]] <- temp[, X][[1]][ind, ]
          }
          }
        } else {
          temp <- matrix(temp, nrow = nMPs, ncol = reps, byrow = TRUE)
          TACa[, , ss] <- temp
        }
      }
    } else {
      for (ff in 1:nMPs) {
        temp <- sfSapply(1:nsims, MPs[ff], Data = Data, 
          reps = reps)
        if (mode(temp) == "numeric") 
          TACa[ff, , ] <- temp
        if (mode(temp) == "list") {
          TACa[ff, , ] <- unlist(temp[1, ])
          for (x in 1:nsims) Data@Misc[[x]] <- temp[2, x][[1]]
        }
      }
    }
  }
  for (ff in 1:nMPs) {
    if (sum(is.na(TACa[ff, , ])) > sum(!is.na(TACa[ff, , ]))) {
      # only plot if there are sufficient non-NA TAC samples
      print(paste("Method ", MPs[ff], " produced greater than 50% NA values", 
        sep = ""))
    }
  }
  out <- list(TACa, Data)
  return(out)
}



#' Conduct stock assessment
#' 
#' A wrapper function that gets the OFL recommendation in cases where a method
#' of DLM quota has been specified
#' 
#' 
#' @usage Sam(Data, MPs = NA, reps = 100, perc = 0.5)
#' @param Data A data-limited methods data object
#' @param MPs A character vector of methods of DLM quota, DLM space or DLM size
#' @param reps The number of samples of quota recommendations by method
#' @param perc quantile of the recommendation to use
#' @author T. Carruthers
#' @export Sam
Sam <- function(Data, MPs = NA, reps = 100, perc = 0.5) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  nm <- deparse(substitute(DLM))
  Data@PosMPs <- MPs
  funcs <- Data@PosMPs
  nMPs <- length(funcs)
  Data@MPs <- funcs
  temp <- getTAC(Data, MPs = funcs, reps)
  TACa <- temp[[1]]
  Data <- temp[[2]]
  nsim <- length(Data@Mort)
  ref <- array(rep(Data@Ref, nMPs), c(nsim, nMPs))
  TACm <- apply(TACa, c(3, 1), quantile, p = perc, na.rm = T)
  TACbias <- (TACm - ref)/ref * 100
  POF <- round(apply(TACbias > 0, 2, sum)/length(Data@Mort) * 100, 
    1)
  Data@TAC <- TACa
  Data@TACbias <- TACbias
  Data
}


#' Boxplot of TAC recommendations
#' 
#' @param x An object of class MSE
#' @param upq Upper quantile of TACs for max ylim
#' @param lwq Lower quantile of TACs for min ylim
#' @param outline Logical. Include outliers in plot?
#' @param ...  Optional additional arguments passed to \code{boxplot}
#' @return Returns a data frame containing the information shown in the plot
#' @author A. Hordyk
#' @export
boxplot.Data <- function(x, upq=0.9, lwq=0.1, outline = FALSE, ...) {
  Data <- updateMSE(x)
  if (class(Data) != "Data")  stop("Object must be of class 'Data'")
  tacs <- t(Data@TAC[, , 1])
  if (all(is.na(tacs))) {
    message("Nothing found in TAC slot")
    return(invisible(NULL))
  }
  units <- TRUE
  if (length(nchar(x@Units)) < 1) units <- FALSE
  MPs <- Data@MPs
  ind <- grep("ref", MPs)
  if (length(ind) > 0) {
    tacs <- tacs[, -ind]
    MPs <- MPs[-ind]
  }
  if (nrow(tacs) > 1) {
    ord <- order(apply(tacs, 2, median, na.rm = TRUE))
    MPs <- MPs[ord]
    tacs <- tacs[, ord]
    cols <- rainbow(30)
    # ymax <- quantile(apply(tacs, 2, quantile, upq, na.rm = TRUE), upq)
    ymax <- max(apply(tacs, 2, quantile, upq, na.rm = TRUE))
    # ymin <- quantile(apply(tacs, 2, quantile, lwq, na.rm = TRUE), lwq)
    ymin <- min(apply(tacs, 2, quantile, lwq, na.rm = TRUE))
    ylim <- c(ymin, ymax)
    Median <- round(apply(tacs, 2, median, na.rm = TRUE), 2)
    SD <- round(apply(tacs, 2, sd, na.rm = TRUE), 2)
  } else {
    ylim <- c(quantile(tacs, lwq),quantile(tacs, upq))
    Median <- median(tacs)
    SD <- sd(tacs)
    tacs <- as.numeric(tacs)
    cols <- "black"
  }
  
  par(mfrow = c(1, 1), oma = c(2, 4, 1, 0), mar = c(3, 3, 0, 0))
  boxplot(tacs, names = MPs, las = 1, col = cols, outline = outline, 
    frame = FALSE, ylim = ylim, horizontal = TRUE, ...)
  
  if (units) mtext(paste("TAC (", Data@Units, ")", sep = ""), side = 1, outer = T, 
                   line = 0.5, cex = 1.25)
  if (!units) mtext("TAC (no units supplied)", side = 1, outer = T, 
                    line = 0.5, cex = 1.25)
  mtext(side = 2, "Management Procedures", outer = TRUE, line = 3, cex = 1.25)
  mtext(paste("TAC calculation for ", Data@Name, sep = ""), 3, outer = T, 
    line = -0.5, cex = 1.25)
  
  if (units) data.frame(MP = MPs, Median = Median, SD = SD, Units = Data@Units)
  if (!units) data.frame(MP = MPs, Median = Median, SD = SD)
  
}

#' Function to run a set of input control methods
#' 
#' Runs a set of input control methods are returns the output in a single table
#' 
#' @param Data A Data object
#' @param MPs A list of input MPs, if NA all available input MPs are run
#' @param reps Number of repetitions (for those methods that use them)
#' @param timelimit Maximum timelimit to run MP (in seconds)
#' @param CheckMPs Logical, the Can function is run if this is TRUE
#' @param msg Logical. Should messages be printed?
#' @author A. Hordyk
#' @export Input
Input <- function(Data, MPs = NA, reps = 100, timelimit = 10, CheckMPs = TRUE, 
                  msg=TRUE) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  if (msg) message("Checking which MPs can be run")

  if (CheckMPs) PosMPs <- Can(Data, timelimit = timelimit)
  if (!CheckMPs) PosMPs <- MPs
  PosMPs <- PosMPs[PosMPs %in% avail("Input")]
  if (!is.na(MPs[1])) Data@MPs <- MPs[MPs %in% PosMPs]
  if (is.na(MPs[1])) Data@MPs <- PosMPs
  funcs <- Data@MPs
  
  if (length(funcs) == 0) {
    stop("None of the methods 'MPs' are possible given the data available")
  } else {
    nareas <- Data@nareas
    areacol <- paste("Area", 1:nareas)
    Out <- matrix(NA, nrow = length(funcs), ncol = 4+nareas)
    colnames(Out) <- c("Effort", "LR5", "LFR", "Harvest Slot Limit", areacol)
    rownames(Out) <- funcs
    
    for (mm in 1:length(funcs)) {
      if (msg) message("Running ", mm, " of ", length(funcs), " - ", funcs[mm])
      
      runIn <- runInMP(Data, MPs = funcs[mm], reps = reps)[[1]][[1]]
      if (length(runIn$Effort) > 0) Out[mm, 1] <- runIn$Effort
      if (length(runIn$LR5) > 0) Out[mm, 2] <- runIn$LR5
      if (length(runIn$LFR) > 0) Out[mm, 3] <- runIn$LFR
      if (length(runIn$HS) > 0) Out[mm, 4] <- runIn$HS
      Out[mm, 5:ncol(Out)] <- runIn$Spatial
    }
  }
  round(Out,2)
  
}


