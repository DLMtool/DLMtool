#' Label class union for performance metric objects
#' 
#' @description Used internally. Nothing to see here!
#' @keywords internal
#' @export
setClassUnion(name="label.class", members=c("call", "character", "function"))

#' Prob class union for performance metric objects
#' 
#' @description Used internally. Nothing to see here!
#' @keywords internal  
#' @export
setClassUnion(name="prob.class", members=c("matrix", "numeric", "data.frame"))


# ---- Data Class ----

#' Class \code{'Data'}
#' 
#' An object for storing fishery data for analysis
#' 
#' 
#' @name Data-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Data', stock)} 
#' @slot Name The name of the Data object. Single value. Character string  
#' @slot Common_Name Common name of the species. Character string
#' @slot Species Scientific name of the species. Genus and species name. Character string
#' @slot Region Name of the general geographic region of the fishery. Character string
#' @slot LHYear The last historical year of the simulation (before projection). Single value. Positive integer 
#' @slot MPrec The previous recommendation of a management procedure. Vector of length nsim. Positive real numbers   
#' @slot Units Units of the catch/absolute abundance estimates. Single value. Character string
#' @slot MPeff The current level of effort. Vector of length nsim. Positive real numbers 
#' @slot nareas Number of fishing areas. Vector of length nsim. Non-negative integer 
#' 
#' @slot MaxAge Maximum age. Vector nsim long. Positive integer
#' @slot Mort Natural mortality rate. Vector nsim long. Positive real numbers 
#' @slot CV_Mort Coefficient of variation in natural mortality rate. Vector nsim long. Positive real numbers 
#' @slot vbLinf Maximum length. Vector nsim long. Positive real numbers
#' @slot CV_vbLinf Coefficient of variation in maximum length. Vector nsim long. Positive real numbers 
#' @slot vbK The von Bertalanffy growth coefficient K. Vector nsim long. Positive real numbers
#' @slot CV_vbK Coefficient of variation in the von Bertalanffy K parameter. Vector nsim long. Positive real numbers 
#' @slot vbt0 Theoretical age at length zero. Vector nsim long. Non-positive real numbers
#' @slot CV_vbt0 Coefficient of variation in age at length zero. Vector nsim long. Positive real numbers 
#' @slot wla Weight-Length parameter alpha. Vector nsim long. Positive real numbers 
#' @slot CV_wla Coefficient of variation in weight-length parameter a. Vector nsim long. Positive real numbers
#' @slot wlb Weight-Length parameter beta. Vector nsim long. Positive real numbers 
#' @slot CV_wlb Coefficient of variation in weight-length parameter b. Vector nsim long. Positive real numbers
#' @slot steep Steepness of stock-recruitment relationship. Vector nsim long. Value in the range of one-fifth to 1 
#' @slot CV_steep Coefficient of variation in steepness. Vector nsim long. Positive real numbers   
#' @slot sigmaR Recruitment variability. Vector nsim long. Positive real numbers
#' @slot CV_sigmaR Coefficient of variation in recruitment variability. Vector nsim long. Positive real numbers 
#' @slot L50 Length at 50 percent maturity. Vector nsim long. Positive real numbers 
#' @slot CV_L50 Coefficient of variation in length at 50 per cent maturity. Vector nsim long. Positive real numbers 
#' @slot L95 Length at 95 percent maturity. Vector nsim long. Positive real numbers 
#' @slot LenCV Coefficient of variation of length-at-age (assumed constant for all age classes). Vector nsim long. Positive real numbers 
#' @slot LFC Length at first capture. Vector nsim long. Positive real numbers 
#' @slot CV_LFC Coefficient of variation in length at first capture. Vector nsim long. Positive real numbers 
#' @slot LFS Shortest length at full selection.  Vector nsim long. Positive real numbers 
#' @slot CV_LFS Coefficient of variation in length at full selection. Vector nsim long. Positive real numbers 
#' @slot Vmaxlen Vulnerability of individuals at asymptotic length. Vector nsim long. Real number between 0 and 1.
#' 
#' @slot Year Years that corresponding to catch and relative abundance data. Vector nyears long. Positive integer
#' @slot Cat Total annual catches. Matrix of nsim rows and nyears columns. Non-negative real numbers 
#' @slot CV_Cat Coefficient of variation in annual catches. Matrix nsim rows and either 1 or nyear columns.
#'  Positive real numbers. Note: built-in MPs use only the first value of `CV_Cat` for all years.
#' @slot Effort Annual fishing effort. Matrix of nsim rows and nyears columns. Non-negative real numbers 
#' @slot CV_Effort Coefficient of variation in annual effort. Matrix nsim rows and either 1 or nyear columns.
#'  Positive real numbers. Note: built-in MPs use only the first value of `CV_Effort` for all years.
#' @slot Ind Relative total abundance index. Matrix of nsim rows and nyears columns. Non-negative real numbers
#' @slot CV_Ind Coefficient of variation in the relative total abundance index. Matrix nsim rows and either 1 or nyear columns.
#'  Positive real numbers. Note: built-in MPs use only the first value of `CV_Ind` for all years
#'  
#' @slot SpInd Relative spawning abundance index. Matrix of nsim rows and nyears columns. Non-negative real numbers
#' @slot CV_SpInd Coefficient of variation in the relative spawning abundance index. Matrix nsim rows and either 1 or nyear columns. Positive real numbers. 
#' 
#' @slot VInd Relative vulnerable abundance index. Matrix of nsim rows and nyears columns. Non-negative real numbers
#' @slot CV_VInd Coefficient of variation in the relative vulnerable abundance index. Matrix nsim rows and either 1 or nyear columns.
#'  Positive real numbers. 
#'  
#' @slot AddInd Optional additional indices. Array of dimensions `nsim`, n additional indices, and `nyears` (length `Year`).  
#' @slot CV_AddInd Coefficient of variation for additional indices. Array of same dimensions as `AddInd`
#' @slot AddIndV Vulnerability-at-age schedules for the additional indices. Array with dimensions: `nsim`, n additional indices, `MaxAge`.  
#' @slot AddIunits Units for the additional indices - biomass (1; default) or numbers (0). Numeric vector length n.ind.
#' @slot AddIndType Index calculated from total stock (1, default), spawning stock (2), or vulnerable stock (3). Numeric vector of length n.ind
#' 
#' @slot Rec Recent recruitment strength. Matrix of nsim rows and nyears columns. Non-negative real numbers
#' @slot CV_Rec Log-normal CV for recent recruitment strength.  Matrix nsim rows and either 1 or nyear columns.
#'  Positive real numbers. Note: built-in MPs use only the first value of `CV_Rec` for all years.
#' @slot ML Mean length time series. Matrix of nsim rows and nyears columns. Non-negative real numbers
#' @slot Lc Modal length of catches. Matrix of nsim rows and nyears columns. Positive real numbers 
#' @slot Lbar Mean length of catches over Lc. Matrix of nsim rows and nyears columns. Positive real numbers  
#' @slot Vuln_CAA Optional vulnerability-at-age schedule for catch-at-age samples. Used to condition OM for closed-loop
#' simulation testing. Replaces the fleet selectivity schedule in the OM used to generate CAA samples. Matrix
#' with dimensions `nsim` x `MaxAge`.
#' @slot CAA Catch at Age data (numbers). Array of dimensions nsim x nyears x MaxAge. Non-negative integers
#' @slot Vuln_CAL Optional vulnerability-at-length schedule for catch-at-length samples. Used to condition OM for closed-loop
#' simulation testing. Replaces the fleet selectivity schedule in the OM used to generate CAL samples. Matrix
#' with dimensions `nsim` x `length(CAL_mids)`.
#' @slot CAL_bins The values delimiting the length bins for the catch-at-length data. Vector. Non-negative real numbers
#' @slot CAL_mids The values of the mid-points of the length bins. Optional, calculated from `CAL_bins` if not entered. Vector. Non-negative real numbers.
#' @slot CAL Catch-at-length data. An array with dimensions nsim x nyears x length(CAL_mids). Non-negative integers  
#'    
#' @slot Dep Stock depletion SSB(current)/SSB(unfished). Vector nsim long. Fraction.  
#' @slot CV_Dep Coefficient of variation in current stock depletion. Vector nsim long. Positive real numbers
#' @slot Abun An estimate of absolute current vulnerable abundance. Vector nsim long. Positive real numbers    
#' @slot CV_Abun Coefficient of variation in estimate of absolute current stock size. Vector nsim long. Positive real numbers 
#' @slot SpAbun An estimate of absolute current spawning stock abundance. Vector nsim long. Positive real numbers
#' @slot CV_SpAbun Coefficient of variation in estimate of absolute spawning current stock size. Vector nsim long. Positive real numbers  
#' 
#' @slot FMSY_M An assumed ratio of FMSY to M. Vector nsim long. Positive real numbers  
#' @slot CV_FMSY_M Coefficient of variation in the ratio in FMSY/M. Vector nsim long. Positive real numbers 
#' @slot BMSY_B0 The most productive stock size relative to unfished. Vector nsim long. Fraction 
#' @slot CV_BMSY_B0 Coefficient of variation in the position of the most productive stock size relative to unfished. Vector nsim long. Positive real numbers 
#' @slot Cref Reference or target catch level (eg MSY). Vector of length nsim. Positive real numbers 
#' @slot CV_Cref Log-normal CV for reference or target catch level. Vector of length nsim. Positive real numbers 
#' @slot Bref Reference or target biomass level (eg BMSY). Vector of length nsim. Positive real numbers 
#' @slot CV_Bref Log-normal CV for reference or target biomass level. Vector of length nsim. Positive real numbers 
#' @slot Iref Reference or target relative abundance index level (eg BMSY / B0). Vector of length nsim. Positive real numbers 
#' @slot CV_Iref Log-normalCV for reference or target relative abundance index level. Vector of length nsim. Positive real numbers
#' @slot t The number of years corresponding to AvC and Dt. Single value. Positive integer  
#' @slot AvC Average catch over time t. Vector nsim long. Positive real numbers  
#' @slot CV_AvC Coefficient of variation in average catches over time t. Vector nsim long. Positive real numbers
#' @slot Dt Depletion over time t SSB(now)/SSB(now-t+1). Vector nsim long. Fraction  
#' @slot CV_Dt Coefficient of variation in depletion over time t. Vector nsim long. Positive real numbers  
#' @slot Ref A reference management level (eg a catch limit). Single value. Positive real number  
#' @slot Ref_type Type of reference management level (eg 2009 catch limit). Single value. Character string 
#' @slot Log A record of events. Single value. Character string 
#' @slot params A place to store estimated parameters. An object. R list 
#' @slot PosMPs The methods that can be applied to these data. Vector. Character strings 
#' @slot TAC The calculated catch limits (function TAC). An array with dimensions PosMPs x replicate TAC samples x nsim. Positive real numbers  
#' @slot Sense The results of the sensitivity analysis (function Sense). An array with dimensions PosMPs x sensitivity increments. Positive real numbers  
#' @slot MPs The methods that were applied to these data. Vector. Character strings  
#' @slot OM A table of operating model conditions. R table object of nsim rows. Real numbers  
#' @slot Obs A table of observation model conditions. R table object of nsim rows. Real numbers 
#' @slot Misc Other information for MPs. An object. R list   
#' 
#' @author T. Carruthers and A. Hordyk
#' @export
#' @keywords classes
#' @examples
#' 
#' newdata<-new('Data')
#' 
setClass("Data", 
         representation(Name = "character", Common_Name='character', 
                        Species='character', Region='character',
                        LHYear = "numeric", MPrec = "vector", 
                        Units = "character", MPeff = "vector", 
                        nareas = "numeric",
                        
                        MaxAge = "vector", Mort = "vector", CV_Mort = "vector",
                        vbLinf = "vector", CV_vbLinf = "vector",
                        vbK = "vector", CV_vbK = "vector",
                        vbt0 = "vector", CV_vbt0 = "vector",
                        wla = "vector", CV_wla = "vector", 
                        wlb = "vector", CV_wlb = "vector",
                        steep = "vector", CV_steep = "vector", 
                        sigmaR='vector', CV_sigmaR='vector',
                        L50 = "vector",  CV_L50 = "vector",
                        L95 = "vector",
                        LenCV="vector",
                        
                        LFC = "vector", CV_LFC = "vector",
                        LFS = "vector", CV_LFS = "vector", 
                        Vmaxlen = 'vector',
                        
                        Year = "vector", 
                        Cat = "matrix", CV_Cat = "matrix",
                        Effort = 'matrix', CV_Effort = 'matrix',
                        Ind = "matrix", CV_Ind = "matrix", 
                        SpInd = "matrix", CV_SpInd = "matrix", 
                        VInd = "matrix", CV_VInd = "matrix", 
                        
                        AddInd = "array", CV_AddInd = "array", AddIndV = "array",
                        AddIunits = 'vector', AddIndType='vector',
                        
                        Rec = "matrix", CV_Rec = "matrix", 
                        ML = "matrix",  Lc = "matrix", Lbar = "matrix", 
                        
                        Vuln_CAA = "matrix", CAA = "array",
                        Vuln_CAL = 'matrix', CAL_bins = "numeric",
                        CAL_mids = "numeric", CAL = "array",
                        
                        Dep = "vector", CV_Dep = "vector",
                        Abun = "vector", CV_Abun = "vector",
                        SpAbun= "vector",  CV_SpAbun = "vector",
                        
                        FMSY_M = "vector", CV_FMSY_M = "vector",
                        BMSY_B0 = "vector", CV_BMSY_B0 = "vector", 
                        
                        Cref = "vector", CV_Cref = "vector",
                        Bref = "vector", CV_Bref = "vector", 
                        Iref = "vector", CV_Iref = "vector", 
                        t = "vector",  AvC = "vector", CV_AvC = "vector",  
                        Dt = "vector",   CV_Dt = "vector",
                        
                        Ref = "numeric", Ref_type = "character", 
                        Log = "list", params = "list", PosMPs = "vector",
                        TAC = "array", Sense = "array",
                        MPs = "vector", OM = "data.frame", Obs = "data.frame", 
                        Misc = "list"))       
                              
                                 
setMethod("initialize", "Data", function(.Object, stock="nada", ...) {
  if (file.exists(stock)) {
    .Object <- XL2Data(stock, ...)  
  } else {
    slots <- slotNames('Data')
    for (x in seq_along(slots)) {
      sl <- slots[x]

      cl <- class(slot(.Object, sl))
      if ("logical" %in% cl) slot(.Object, sl) <- as.numeric(NA)
      if ("character" %in% cl) slot(.Object, sl) <- ''
      if ("matrix" %in% cl) {
        if (length(dim(slot(.Object, sl))) > 2) {
          slot(.Object, sl) <- array(NA, dim=c(1,1,1))
        } else {
          slot(.Object, sl) <- matrix(NA)  
        }
      }
      if ("array" %in% cl) {
        if (length(dim(slot(.Object, sl))) > 2) {
          slot(.Object, sl) <- array(NA, dim=c(1,1,1))
        } else {
          slot(.Object, sl) <- matrix(NA)  
        }
      }
      if ("vector" %in% cl) slot(.Object, sl) <- NA
      if ("numeric" %in% cl) slot(.Object, sl) <- as.numeric(NA)
      if ("list" %in% cl) slot(.Object, sl) <- list()
    }
  }
  
  # Default values
  if (all(is.na(.Object@CV_Cat))) .Object@CV_Cat <- matrix(0.2, nrow=1, ncol=1)
  if (all(is.na(.Object@CV_Ind))) .Object@CV_Ind <- matrix(0.2, nrow=1, ncol=1)
  if (all(is.na(.Object@CV_SpInd))) .Object@CV_SpInd <- matrix(0.2, nrow=1, ncol=1)
  if (all(is.na(.Object@CV_VInd))) .Object@CV_VInd <- matrix(0.2, nrow=1, ncol=1)
  if (all(is.na(.Object@CV_Effort))) .Object@CV_Effort <- matrix(0.2, nrow=1, ncol=1)
  if (all(is.na(.Object@CV_Rec))) .Object@CV_Rec <- matrix(0.2, nrow=1, ncol=1)
  
  if (NAor0(.Object@LenCV)) .Object@LenCV <- 0.1
  if (NAor0(.Object@CV_Dt)) .Object@CV_Dt <- 0.25
  if (NAor0(.Object@CV_AvC)) .Object@CV_AvC <- 0.2
  if (NAor0(.Object@CV_Mort)) .Object@CV_Mort <- 0.2
  if (NAor0(.Object@CV_FMSY_M)) .Object@CV_FMSY_M <- 0.2
  if (NAor0(.Object@CV_BMSY_B0)) .Object@CV_BMSY_B0 <- 0.045
  if (NAor0(.Object@CV_Cref)) .Object@CV_Cref <- 0.2
  if (NAor0(.Object@CV_Bref)) .Object@CV_Bref <- 0.2
  if (NAor0(.Object@CV_Iref)) .Object@CV_Iref <- 0.2
  if (NAor0(.Object@CV_Dep)) .Object@CV_Dep <- 0.25
  if (NAor0(.Object@CV_Abun)) .Object@CV_Abun <- 0.25
  if (NAor0(.Object@CV_vbK)) .Object@CV_vbK <- 0.1
  if (NAor0(.Object@CV_vbLinf)) .Object@CV_vbLinf <- 0.1
  if (NAor0(.Object@CV_vbt0)) .Object@CV_vbt0 <- 0.1
  if (NAor0(.Object@CV_L50))  .Object@CV_L50 <- 0.1
  if (NAor0(.Object@CV_LFC))  .Object@CV_LFC <- 0.2
  if (NAor0(.Object@CV_LFS))  .Object@CV_LFS <- 0.2
  if (NAor0(.Object@CV_wla))  .Object@CV_wla <- 0.1
  if (NAor0(.Object@CV_wlb))  .Object@CV_wlb <- 0.1
  if (NAor0(.Object@CV_steep)) .Object@CV_steep <- 0.2
  if (NAor0(.Object@nareas)) .Object@nareas <- 2
  
  if (length(.Object@CAA) == 0) .Object@CAA <- array(NA, c(1, 1, 1))
  if (length(.Object@CAL) == 0) .Object@CAL <- array(NA, c(1, 1, 1))
  if (length(.Object@CAL_bins) == 0) .Object@CAL_bins <- 1
  if (length(.Object@TAC) == 0) .Object@TAC <- array(1, c(1, 1))
  # if (length(.Object@TACbias) == 0) .Object@TACbias <- array(1, c(1, 1))
  if (length(.Object@Sense) == 0) .Object@Sense <- array(1, c(1, 1))
  if (length(.Object@ML) == 0)  .Object@ML <- array(NA, c(1, 1))
  if (length(.Object@Lbar) == 0) .Object@Lbar <- array(NA, c(1, 1))
  if (length(.Object@Lc) == 0) .Object@Lc <- array(NA, c(1, 1))
  
  return(.Object)
})

# setMethod("initialize", "Data", function(.Object, stock = "nada", dec=c(".", ","), silent=TRUE) {
#   if (file.exists(stock)) {
#     dec <- match.arg(dec)
#     Ncol <- max(unlist(lapply(strsplit(readLines(stock), ","), length)))
#     col.names <- paste0("V", 1:Ncol)
# 
#     if (dec == ".")
#       dat <- read.csv(stock, header = F, colClasses = "character", col.names=col.names)  # read 1st sheet
#     if (dec == ",")
#       dat <- read.csv2(stock, header = F, colClasses = "character", col.names=col.names)  # read 1st sheet
#     dname <- dat[, 1]
#     dat <- dat[, 2:ncol(dat)]
# 
#     .Object@Name <- dat[match("Name", dname), 1]
#     .Object@Common_Name <- dat[match("Common Name", dname), 1]
#     .Object@Species <- dat[match("Species", dname), 1]
#     .Object@Region <- dat[match("Region", dname), 1]
#     .Object@Year <- as.numeric(dat[match("Year", dname), dat[match("Year", dname), ] != ""])
#     # .Object@Cat <- matrix(as.numeric(dat[match("Catch", dname), dat[match("Catch", dname), ] != ""]), nrow = 1)
#     .Object@Cat <- matrix(as.numeric(dat[match("Catch", dname), 1:length(.Object@Year)]), nrow = 1)
#     .Object@Ind <- matrix(as.numeric(dat[match("Abundance index", dname), 1:length(.Object@Year)]), nrow = 1)
# 
#     # .Object@Type <- dat[match("Index type", dname), ] %>% as.character()
#     # .Object@Type <- .Object@Type[!is.na(.Object@Type)]
#     # .Object@Type <- .Object@Type[nchar(.Object@Type)>0]
#     # if (length(.Object@Type)>0) {
#     #   .Object@Type <- .Object@Type[nchar(.Object@Type)>0] %>% as.character()
#     #   n.ind <- length(.Object@Type)
#     #   r.ind <- match("Real indices", dname)
#     #   RInd.dat <- dat[r.ind:(r.ind+n.ind-1),] %>% data.matrix()
#     #   .Object@RInd <- array(RInd.dat, dim=c(1, n.ind, ncol(RInd.dat)))
#     # }
# 
#     .Object@Rec <- matrix(as.numeric(dat[match("Recruitment", dname), 1:length(.Object@Year)]), nrow = 1)
#     .Object@t <- as.numeric(dat[match("Duration t", dname), 1])
#     .Object@AvC <- as.numeric(dat[match("Average catch over time t", dname), 1])
#     .Object@Dt <- as.numeric(dat[match("Depletion over time t", dname), 1])
#     .Object@Mort <- as.numeric(dat[match("M", dname), 1])
#     .Object@FMSY_M <- as.numeric(dat[match("FMSY/M", dname), 1])
#     .Object@BMSY_B0 <- as.numeric(dat[match("BMSY/B0", dname), 1])
#     .Object@L50 <- as.numeric(dat[match("Length at 50% maturity", dname), 1])
#     .Object@L95 <- as.numeric(dat[match("Length at 95% maturity", dname), 1])
#     .Object@ML <- matrix(as.numeric(dat[match("Mean length", dname), 1:length(.Object@Year)]), nrow = 1)
#     .Object@Lbar <- matrix(as.numeric(dat[match("Mean length Lc", dname), 1:length(.Object@Year)]), nrow = 1)
#     .Object@Lc <- matrix(as.numeric(dat[match("Modal length", dname), 1:length(.Object@Year)]), nrow = 1)
#     .Object@LFC <- as.numeric(dat[match("Length at first capture",  dname), 1])
#     .Object@LFS <- as.numeric(dat[match("Length at full selection", dname), 1])
# 
#     CAAy <- grep("CAA", dname)[1:length(grep("CAA", dname))]
#     CAAa <- sum(dat[CAAy[1], ] != "")
#     if (!is.na(CAAa)) {
#       .Object@CAA <- array(as.numeric(as.matrix(dat[CAAy, 1:CAAa])),  dim = c(1, length(CAAy), CAAa))
#     }
#     .Object@Dep <- as.numeric(dat[match("Current stock depletion",  dname), 1])
#     .Object@Abun <- as.numeric(dat[match("Current stock abundance",  dname), 1])
#     .Object@SpAbun <- as.numeric(dat[match("Current spawning stock abundance",  dname), 1])
#     .Object@vbK <- as.numeric(dat[match("Von Bertalanffy K parameter", dname), 1])
#     .Object@vbLinf <- as.numeric(dat[match("Von Bertalanffy Linf parameter", dname), 1])
#     .Object@vbt0 <- as.numeric(dat[match("Von Bertalanffy t0 parameter", dname), 1])
#     .Object@LenCV <- as.numeric(dat[match("CV of length-at-age", dname), 1])
#     .Object@wla <- as.numeric(dat[match("Length-weight parameter a", dname), 1])
#     .Object@wlb <- as.numeric(dat[match("Length-weight parameter b", dname), 1])
#     .Object@steep <- as.numeric(dat[match("Steepness", dname), 1])
#     .Object@sigmaR <- as.numeric(dat[match("sigmaR", dname), 1])
# 
#     .Object@CV_Cat <- as.numeric(dat[match("CV Catch", dname), 1])
#     .Object@CV_Dt <- as.numeric(dat[match("CV Depletion over time t", dname), 1])
#     .Object@CV_AvC <- as.numeric(dat[match("CV Average catch over time t", dname), 1])
#     .Object@CV_Ind <- as.numeric(dat[match("CV Abundance index", dname), 1])
#     .Object@CV_Mort <- as.numeric(dat[match("CV M", dname), 1])
#     .Object@CV_FMSY_M <- as.numeric(dat[match("CV FMSY/M", dname),  1])
#     .Object@CV_BMSY_B0 <- as.numeric(dat[match("CV BMSY/B0", dname), 1])
#     .Object@CV_Dep <- as.numeric(dat[match("CV current stock depletion", dname), 1])
#     .Object@CV_Abun <- as.numeric(dat[match("CV current stock abundance", dname), 1])
#     .Object@CV_vbK <- as.numeric(dat[match("CV von B. K parameter", dname), 1])
#     .Object@CV_vbLinf <- as.numeric(dat[match("CV von B. Linf parameter", dname), 1])
#     .Object@CV_vbt0 <- as.numeric(dat[match("CV von B. t0 parameter", dname), 1])
#     .Object@CV_L50 <- as.numeric(dat[match("CV Length at 50% maturity", dname), 1])
#     .Object@CV_LFC <- as.numeric(dat[match("CV Length at first capture", dname), 1])
#     .Object@CV_LFS <- as.numeric(dat[match("CV Length at full selection", dname), 1])
#     .Object@CV_wla <- as.numeric(dat[match("CV Length-weight parameter a", dname), 1])
#     .Object@CV_wlb <- as.numeric(dat[match("CV Length-weight parameter b", dname), 1])
#     .Object@CV_steep <- as.numeric(dat[match("CV Steepness", dname),  1])
#     .Object@sigmaL <- as.numeric(dat[match("Sigma length composition", dname), 1])
# 
#     .Object@MaxAge <- as.numeric(dat[match("Maximum age", dname), 1])
# 
#     if (length(grep("CAL", dname)) > 1) {
#       CAL_bins <- as.numeric(dat[match("CAL_bins", dname), dat[match("CAL_bins", dname), ] != ""])
#       nCAL <- length(CAL_bins) - 1
#       .Object@CAL_bins <- CAL_bins
#       CALdat <- grep("CAL ", dname)
#       if (length(CALdat) > 0) .Object@CAL <- array(as.numeric(as.matrix(dat[CALdat, 1:nCAL])),dim = c(1, length(CALdat), nCAL))
#     }
# 
#     .Object@Units <- dat[match("Units", dname), 1]
#     .Object@Ref <- as.numeric(dat[match("Reference OFL", dname), 1])
#     .Object@Ref_type <- dat[match("Reference OFL type", dname), 1]
# 
#     .Object@Cref <- as.numeric(dat[match("Cref", dname), 1])
#     .Object@Iref <- as.numeric(dat[match("Iref", dname), 1])
#     .Object@Bref <- as.numeric(dat[match("Bref", dname), 1])
# 
#     .Object@CV_Cref <- as.numeric(dat[match("CV Cref", dname), 1])
#     .Object@CV_Iref <- as.numeric(dat[match("CV Iref", dname), 1])
#     .Object@CV_Bref <- as.numeric(dat[match("CV Bref", dname), 1])
#     .Object@CV_Rec <- as.numeric(dat[match("CV Rec", dname), 1])
# 
#     .Object@MPrec <- as.numeric(dat[match("MPrec", dname), 1])
#     .Object@MPeff <- as.numeric(dat[match("MPeff", dname), 1])
# 
#     .Object@LHYear <- as.numeric(dat[match("LHYear", dname), 1])
#     .Object@nareas <- as.numeric(dat[match("nareas", dname), 1])
# 
#     .Object@Log[[1]] <- paste("Created:", Sys.time())
#     .Object@params <- new("list")
#     .Object@OM <- data.frame(NA)
#     .Object@Obs <- data.frame(NA)
#     .Object@TAC <- array(NA, dim = c(1, 1, 1))
#     .Object@Sense <- array(NA, dim = c(1, 1, 1))
#     .Object@PosMPs <- NA
#     .Object@MPs <- NA
# 
#   } else {
#     if (stock != "MSE") {
#       if (!is.na(stock) && !silent) print("Couldn't find specified csv file, blank DLM object created")
#     }
#   }
# 
#   if (is.na(.Object@MPeff) || length(.Object@MPeff)==0) .Object@MPeff <- 1
# 
#   # Standardise Index if not already
#   .Object@Ind <- .Object@Ind/mean(.Object@Ind, na.rm=TRUE)
# 
#   # Default value
#   if (NAor0(.Object@LenCV)) .Object@LenCV <- 0.1
#   if (NAor0(.Object@CV_Cat)) .Object@CV_Cat <- 0.2
#   if (NAor0(.Object@CV_Dt)) .Object@CV_Dt <- 0.25
#   if (NAor0(.Object@CV_AvC)) .Object@CV_AvC <- 0.2
#   if (NAor0(.Object@CV_Ind)) .Object@CV_Ind <- 0.2
#   if (NAor0(.Object@CV_Mort)) .Object@CV_Mort <- 0.2
#   if (NAor0(.Object@CV_FMSY_M)) .Object@CV_FMSY_M <- 0.2
#   if (NAor0(.Object@CV_BMSY_B0)) .Object@CV_BMSY_B0 <- 0.045
#   if (NAor0(.Object@CV_Cref)) .Object@CV_Cref <- 0.2
#   if (NAor0(.Object@CV_Bref)) .Object@CV_Bref <- 0.2
#   if (NAor0(.Object@CV_Iref)) .Object@CV_Iref <- 0.2
#   if (NAor0(.Object@CV_Rec)) .Object@CV_Rec <- 0.2
#   if (NAor0(.Object@CV_Dep)) .Object@CV_Dep <- 0.25
#   if (NAor0(.Object@CV_Abun)) .Object@CV_Abun <- 0.25
#   if (NAor0(.Object@CV_vbK)) .Object@CV_vbK <- 0.1
#   if (NAor0(.Object@CV_vbLinf)) .Object@CV_vbLinf <- 0.1
#   if (NAor0(.Object@CV_vbt0)) .Object@CV_vbt0 <- 0.1
#   if (NAor0(.Object@CV_L50))  .Object@CV_L50 <- 0.1
#   if (NAor0(.Object@CV_LFC))  .Object@CV_LFC <- 0.2
#   if (NAor0(.Object@CV_LFS))  .Object@CV_LFS <- 0.2
#   if (NAor0(.Object@CV_wla))  .Object@CV_wla <- 0.1
#   if (NAor0(.Object@CV_wlb))  .Object@CV_wlb <- 0.1
#   if (NAor0(.Object@CV_steep)) .Object@CV_steep <- 0.2
#   if (NAor0(.Object@nareas)) .Object@nareas <- 2
# 
#   if (length(.Object@sigmaL) == 0) .Object@sigmaL <- 0.2
#   if (length(.Object@CAA) == 0) .Object@CAA <- array(NA, c(1, 1, 1))
#   if (length(.Object@CAL) == 0) .Object@CAL <- array(NA, c(1, 1, 1))
#   if (length(.Object@CAL_bins) == 0) .Object@CAL_bins <- 1
#   if (length(.Object@TAC) == 0) .Object@TAC <- array(1, c(1, 1))
#   # if (length(.Object@TACbias) == 0) .Object@TACbias <- array(1, c(1, 1))
#   if (length(.Object@Sense) == 0) .Object@Sense <- array(1, c(1, 1))
#   if (length(.Object@ML) == 0)  .Object@ML <- array(NA, c(1, 1))
#   if (length(.Object@Lbar) == 0) .Object@Lbar <- array(NA, c(1, 1))
#   if (length(.Object@Lc) == 0) .Object@Lc <- array(NA, c(1, 1))
# 
#   if (length(.Object@Type) == 0) .Object@Type <- NA
#   if (length(.Object@RInd) == 0) .Object@RInd <- array(NA, c(1,1,1))
# 
#   .Object
# })



# # ---- Fease Class ----
# #' Class \code{'Fease'}
# #' 
# #' An object for storing information about what data are available or might be
# #' available
# #' 
# #' @name Fease-class
# #' @docType class
# #' @section Objects from the Class: Objects can be created by calls of the form
# #' \code{new('Fease', stock)}
# #'
# #' @slot Name The name of the data feasibility object
# #' @slot Case The names of the data feasibility cases
# #' @slot Catch Total annual catches
# #' @slot Index An index of relative abundance, catch per unit effort data or of fishing mortality rate (effort)
# #' @slot Natural_mortality_rate From Maximum age, Tagging data, early fishery catch composition data
# #' @slot Maturity_at_length From gonadal analysis, growth and natural mortality rate estimates
# #' @slot Growth Paired length and age observations, maximum length and an estimate of natural mortality rate
# #' @slot Length_weight_conversion Paired weight and length observations, equivalent data from a similar species
# #' @slot Fleet_selectivity Length composition of catches with growth curve and natural mortality rate, estimates from a similar fleet type targetting a similar species
# #' @slot Catch_at_length Length composition of catches (length samples)
# #' @slot Catch_at_age Age composition of catches (age samples)
# #' @slot Recruitment_index Spawn survey, estimates from a stock assessment, VPA analysis of catch composition data
# #' @slot Stock_recruitment_relationship Stock assessment, a stock assessment of a similar species
# #' @slot Target_catch An agreed annual catch target, MSY proxy
# #' @slot Target_biomass An agreed absolute biomass target, mean historical biomass estimate
# #' @slot Target_index An agreed catch rate target
# #' @slot Abundance Fishery independent survey, current fishing mortality rate from recent length composition, natural mortality rate, maturity at age, growth and stock recruitment relationship, habitat and relative density extrapolation
# #'
# #' @author T. Carruthers and A. Hordyk
# #' @keywords classes
# #' @examples
# #' 
# #' newdata<-new('Fease')
# #' 
# setClass("Fease", representation(Name = "character", Case = "character", 
#                                  Catch = "numeric", Index = "numeric", Natural_mortality_rate = "numeric", 
#                                  Maturity_at_length = "numeric", Growth = "numeric", Length_weight_conversion = "numeric", 
#                                  Fleet_selectivity = "numeric", Catch_at_length = "numeric", Catch_at_age = "numeric", 
#                                  Recruitment_index = "numeric", Stock_recruitment_relationship = "numeric", 
#                                  Target_catch = "numeric", Target_biomass = "numeric", Target_index = "numeric", 
#                                  Abundance = "numeric"))
# 
# # initialize Fease
# setMethod("initialize", "Fease", function(.Object, file = "nada", ncases = 1, dec=c(".", ",")) {
#   # run an error check here
#   if (file.exists(file)) {
#     Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
#     dec <- match.arg(dec)
#     if (dec == ".") dat <- read.csv(file, header = F, colClasses = "character", col.names = paste0("V", 
#                                                                                    1:Ncol))  # read 1st sheet
#     if (dec == ",") dat <- read.csv2(file, header = F, colClasses = "character", col.names = paste0("V", 
#                                                                                                    1:Ncol))  # read 1st sheet
#     nr <- nrow(dat)
#     ncases = ncol(dat) - 1
#     dname <- dat[, 1]
#     if (ncases == 1) dat <- array(dat[, 2:ncol(dat)], dim = c(nr, ncases))
#     if (ncases > 1) dat <- dat[, 2:ncol(dat)]
#     .Object@Name <- dat[match("Name", dname), 1]
#     .Object@Case <- as.character(dat[match("Case", dname), 1:ncases])
#     .Object@Catch <- as.numeric(dat[match("Catch", dname), 1:ncases])
#     .Object@Index <- as.numeric(dat[match("Index", dname), 1:ncases])
#     .Object@Natural_mortality_rate <- as.numeric(dat[match("Natural_mortality_rate",  dname), 1:ncases])
#     .Object@Maturity_at_length <- as.numeric(dat[match("Maturity_at_length",  dname), 1:ncases])
#     .Object@Growth <- as.numeric(dat[match("Growth", dname), 1:ncases])
#     .Object@Length_weight_conversion <- as.numeric(dat[match("Length_weight_conversion", dname), 1:ncases])
#     .Object@Fleet_selectivity <- as.numeric(dat[match("Fleet_selectivity", dname), 1:ncases])
#     .Object@Catch_at_length <- as.numeric(dat[match("Catch_at_length", dname), 1:ncases])
#     .Object@Catch_at_age <- as.numeric(dat[match("Catch_at_age", dname), 1:ncases])
#     .Object@Recruitment_index <- as.numeric(dat[match("Recruitment_index", dname), 1:ncases])
#     .Object@Stock_recruitment_relationship <- as.numeric(dat[match("Stock_recruitment_relationship", dname), 1:ncases])
#     .Object@Target_catch <- as.numeric(dat[match("Target_catch", dname), 1:ncases])
#     .Object@Target_biomass <- as.numeric(dat[match("Target_biomass", dname), 1:ncases])
#     .Object@Target_index <- as.numeric(dat[match("Target_index", dname), 1:ncases])
#     .Object@Abundance <- as.numeric(dat[match("Abundance", dname), 1:ncases])
#   } else {
#     .Object@Name <- "Blank DLM_Fease"
#     .Object@Case <- "Case 1"
#     .Object@Catch <- 1
#     .Object@Index <- 1
#     .Object@Natural_mortality_rate <- 1
#     .Object@Maturity_at_length <- 1
#     .Object@Growth <- 1
#     .Object@Length_weight_conversion <- 1
#     .Object@Fleet_selectivity <- 1
#     .Object@Catch_at_length <- 1
#     .Object@Catch_at_age <- 1
#     .Object@Recruitment_index <- 1
#     .Object@Stock_recruitment_relationship <- 1
#     .Object@Target_catch <- 1
#     .Object@Target_biomass <- 1
#     .Object@Target_index <- 1
#     .Object@Abundance <- 1
#   }
#   .Object
#   
# })
# 

# ---- Stock Class ----
#' Class \code{'Stock'}
#' 
#' An operating model component that specifies the parameters of the population
#' dynamics model
#' 

#' @name Stock-class
#' @docType class
#' 
#' @slot Name The name of the Stock object. Single value. Character string 
#' @template Stock_template
#' 
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Stock')}
#' 
#' @author T. Carruthers and A. Hordyk
#' @export
#' @keywords classes
#' @examples
#' 
#' showClass('Stock')
#' 
setClass("Stock", representation(Name = "character", Common_Name='character', Species="character",
                                 maxage = "numeric", 
                                 R0 = "numeric", M = "numeric", M2 = "numeric", 
                                 Mexp="numeric",  Msd = "numeric", Mgrad = "numeric",  
                                 h = "numeric", SRrel = "numeric", Perr = "numeric", AC = "numeric",
                                 Period = "numeric", Amplitude = "numeric",
                                 Linf = "numeric", K = "numeric", t0 = "numeric", LenCV="numeric", 
                                 Ksd = "numeric", Kgrad = "numeric", Linfsd = "numeric", Linfgrad = "numeric",
                                 L50 = "numeric", L50_95 = "numeric", 
                                 D = "numeric", 
                                 a = "numeric", b = "numeric",
                                 Size_area_1 = "numeric", Frac_area_1 = "numeric", Prob_staying = "numeric", 
                                 Fdisc="numeric", Source = "character"))

# initialize Stock
setMethod("initialize", "Stock", function(.Object, file = NA, dec=c(".", ",")) {
  
  if (!is.na(file)) {
    if (file.exists(file)) {
      dec <- match.arg(dec)
      Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
      if (dec == ".") dat <- read.csv(file, header = F, colClasses = "character", col.names = paste0("V", 1:Ncol))  # read 1st sheet
      if (dec == ",") dat <- read.csv2(file, header = F, colClasses = "character", col.names = paste0("V", 1:Ncol))  # read 1st sheet
      dname <- dat[, 1]
      dat <- dat[, 2:ncol(dat)]
      Ncol <- ncol(dat)
      .Object@Name <- dat[match("Name", dname), 1]
      commonname <- dat[match("Common_Name", dname), 1]
      commonname1 <- dat[match("Common Name", dname), 1]
      if (length(commonname)>0) {
        .Object@Common_Name <-commonname
      } else if (length(commonname1)>0) {
        .Object@Common_Name <-commonname1
      } else {
        .Object@Common_Name <- "Not specified"
      }
      .Object@Species <- dat[match("Species", dname), 1]
      .Object@maxage <- as.numeric(dat[match("maxage", dname), 1])
      .Object@R0 <- as.numeric(dat[match("R0", dname), 1])
      options(warn=-1)
      temp <- as.numeric(dat[match("M", dname), 1:Ncol])
      .Object@M <- temp
      nas <- which(is.na(temp))
      if (length(nas) > 0) {
        .Object@M <- temp[1:(nas[1]-1)]
      }
      temp <- as.numeric(dat[match("M2", dname), 1:Ncol])
      .Object@M2 <- temp
      nas <- which(is.na(temp))
      if (length(nas) > 0) {
        .Object@M2 <- temp[1:(nas[1]-1)]
      }
      options(warn=1)
      .Object@Msd <- as.numeric(dat[match("Msd", dname), 1:2])
      .Object@Mgrad <- as.numeric(dat[match("Mgrad", dname), 1:2])
      .Object@Mexp <- as.numeric(dat[match("Mexp", dname), 1:2])
      .Object@Fdisc <- as.numeric(dat[match("Fdisc", dname), 1:2])
      .Object@h <- as.numeric(dat[match("h", dname), 1:2])
      .Object@SRrel <- as.numeric(dat[match("SRrel", dname), 1])
      .Object@Linf <- as.numeric(dat[match("Linf", dname), 1:2])
      .Object@K <- as.numeric(dat[match("K", dname), 1:2])
      .Object@t0 <- as.numeric(dat[match("t0", dname), 1:2])
      .Object@LenCV <- as.numeric(dat[match("LenCV", dname), 1:2])
      .Object@Ksd <- as.numeric(dat[match("Ksd", dname), 1:2])
      .Object@Kgrad <- as.numeric(dat[match("Kgrad", dname), 1:2])
      .Object@Linfsd <- as.numeric(dat[match("Linfsd", dname), 1:2])
      .Object@Linfgrad <- as.numeric(dat[match("Linfgrad", dname),  1:2])
      #.Object@recgrad <- as.numeric(dat[match("recgrad", dname), 1:2])
      .Object@a <- as.numeric(dat[match("a", dname), 1])
      .Object@b <- as.numeric(dat[match("b", dname), 1])
      .Object@D <- as.numeric(dat[match("D", dname), 1:2])
      .Object@Perr <- as.numeric(dat[match("Perr", dname), 1:2])
      .Object@Period <- as.numeric(dat[match("Period", dname), 1:2])
      .Object@Amplitude <- as.numeric(dat[match("Amplitude", dname), 1:2])
      .Object@AC <- as.numeric(dat[match("AC", dname), 1:2])
      .Object@Size_area_1 <- as.numeric(dat[match("Size_area_1", dname), 1:2])
      .Object@Frac_area_1 <- as.numeric(dat[match("Frac_area_1", dname), 1:2])
      .Object@Prob_staying <- as.numeric(dat[match("Prob_staying", dname), 1:2])
      .Object@L50 <- as.numeric(dat[match("L50", dname), 1:2])
      .Object@L50_95 <- as.numeric(dat[match("L50_95", dname), 1:2])
      
      # .Object@FecB <- as.numeric(dat[match("FecB", dname), 1:2])
      
      .Object@Source <- dat[match("Source", dname), 1]
    } else {
      message("File doesn't exist")
    }
  }
  if (all(is.na(.Object@LenCV))) .Object@LenCV <- c(0.08, 0.12)
  # if (all(is.na(.Object@recgrad))) .Object@recgrad <- c(0, 0) # recgrad not currently used
  if (all(is.na(.Object@Size_area_1))) .Object@Size_area_1 <- .Object@Frac_area_1
  if (length(.Object@Size_area_1) == 0) .Object@Size_area_1 <- .Object@Frac_area_1
  
  .Object
  
})

#' Class union for isRel slot
#' 
#' @description Used internally.
#' @keywords internal
#' @export
setClassUnion(name="char.log", members=c("character", "logical"))


# ---- Fleet Class -----
#' Class \code{'Fleet'}
#' 
#' The component of the operating model that controls fishing dynamics
#' 
#' @name Fleet-class
#' @docType class
#' @slot Name Name of the Fleet object. Single value. Character string.
#' @template Fleet_template
#' 
#' @section Creating Object: 
#' Objects can be created by calls of the form \code{new('Fleet')}
#' 
#' @author T. Carruthers and A. Hordyk
#' @export
#' @keywords classes
#' @examples
#' 
#' showClass('Fleet')
#' 
setClass("Fleet", slots = c(Name = "character", nyears = "numeric", Spat_targ = "numeric", 
                            EffYears = "numeric", EffLower = "numeric", EffUpper = "numeric", 
                            Esd = "numeric", qinc = "numeric", qcv = "numeric",   
                            L5 = "numeric", LFS = "numeric", Vmaxlen = "numeric", isRel = "char.log",
                            LR5 = "numeric", LFR = "numeric", Rmaxlen = "numeric", DR = "numeric",
                            SelYears = "numeric", AbsSelYears = "numeric",
                            L5Lower = "numeric", L5Upper = "numeric", LFSLower = "numeric", 
                            LFSUpper = "numeric", VmaxLower = "numeric", 
                            VmaxUpper = "numeric", CurrentYr="numeric", MPA='matrix'))

# initialize Fleet
setMethod("initialize", "Fleet", function(.Object, file = NA, dec=c(".", ",")) {
  if (!is.na(file)) {
    if (file.exists(file)) {
      dec <- match.arg(dec)
      Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
      if (dec == ".") dat <- read.csv(file, header = F, colClasses = "character", 
                                      col.names = paste0("V", 1:Ncol))  # read 1st sheet
      if (dec == ",") dat <- read.csv2(file, header = F, colClasses = "character", 
                                       col.names = paste0("V", 1:Ncol))  # read 1st sheet
      dname <- dat[, 1]
      dat <- dat[, 2:ncol(dat)]
      
      .Object@Name <- dat[match("Name", dname), 1]
      .Object@nyears <- as.numeric(dat[match("nyears", dname), 1])
      
      .Object@CurrentYr <- as.numeric(dat[match("CurrentYr", dname), 1])
      if(is.na(.Object@CurrentYr)).Object@CurrentYr<-.Object@nyears
      
      .Object@Spat_targ <- as.numeric(dat[match("Spat_targ", dname),  1:2])
      if (!is.na(match("Esd", dname))) {
        .Object@Esd <- as.numeric(dat[match("Esd", dname), 1:2])  
      } else {
        .Object@Esd <- as.numeric(dat[match("Fsd", dname), 1:2])  
      }
      
      # .Object@Fgrad<-as.numeric(dat[match('Fgrad',dname),1:2])
      nEffYears <- ncol(dat[match("EffYears", dname), ])
      oldw <- getOption("warn")
      options(warn = -1)
      chk <- as.numeric(dat[match("EffYears", dname), 1:nEffYears])
      options(warn = oldw)
      ind <- which(!is.na(chk))
      nEffYears <- length(ind)
      .Object@EffYears <- as.numeric(dat[match("EffYears", dname),  1:nEffYears])
      .Object@EffLower <- as.numeric(dat[match("EffLower", dname),  1:nEffYears])
      .Object@EffUpper <- as.numeric(dat[match("EffUpper", dname),  1:nEffYears])
      .Object@qinc <- as.numeric(dat[match("qinc", dname), 1:2])
      .Object@qcv <- as.numeric(dat[match("qcv", dname), 1:2])
      
      chkName <- match("SelYears", dname)  # Check if vector of selectivity years exists
      if (is.finite(chkName)) {
        nSelYears <- ncol(dat[match("SelYears", dname), ])
        oldw <- getOption("warn")
        options(warn = -1)
        chk <- as.numeric(dat[match("SelYears", dname), 1:nSelYears])
        options(warn = oldw)
        ind <- which(is.finite(chk))
        nSelYears <- length(ind)
        chk <- length(ind)
        if (is.finite(chk) & chk > 0) {
          # parameters for selectivity years exists
          .Object@SelYears <- as.numeric(dat[match("SelYears", dname), 1:nSelYears])
          .Object@L5Lower <- as.numeric(dat[match("L5Lower", dname), 1:nSelYears])
          .Object@L5Upper <- as.numeric(dat[match("L5Upper", dname), 1:nSelYears])
          .Object@LFSLower <- as.numeric(dat[match("LFSLower", dname), 1:nSelYears])
          .Object@LFSUpper <- as.numeric(dat[match("LFSUpper", dname), 1:nSelYears])
          .Object@VmaxLower <- as.numeric(dat[match("VmaxLower", dname), 1:nSelYears])
          .Object@VmaxUpper <- as.numeric(dat[match("VmaxUpper", dname), 1:nSelYears])
        }
      }
      # These are ignored in MSE if L5Lower etc are set
      .Object@L5 <- as.numeric(dat[match("L5", dname), 1:2])
      .Object@LFS <- as.numeric(dat[match("LFS", dname), 1:2])
      .Object@Vmaxlen <- as.numeric(dat[match("Vmaxlen", dname), 1:2])
      
      # Retention curve parameters 
      .Object@LR5 <- as.numeric(dat[match("LR5", dname), 1:2])
      .Object@LFR <- as.numeric(dat[match("LFR", dname), 1:2])
      .Object@Rmaxlen <- as.numeric(dat[match("Rmaxlen", dname), 1:2])
      .Object@DR <- as.numeric(dat[match("DR", dname), 1:2])
      
      .Object@isRel <- dat[match("isRel", dname), 1]  # Are selecivity parameters relative to maturity?
      if (NAor0(.Object@isRel)) .Object@isRel <- "TRUE"
      .Object@isRel <- as.character(.Object@isRel)
      
      isMPA <- grep('MPA', dname)
      if (length(isMPA)<1) isMPA <- NA
      if (!is.na(isMPA)) {
        MPAdat <- dat[isMPA:nrow(dat),]
        nrow <- min(which(is.na(as.numeric(MPAdat[1,]))))
        MPAdat <- MPAdat[,1:(nrow-1)]
        MPAdat <- as.matrix(sapply(MPAdat, as.numeric))  
        .Object@MPA <- MPAdat
      }
      
    } else {
      message("File doesn't exist")
    }
  }
  .Object
})





#' ~~ Methods for Function \code{initialize} ~~
#' 
#' ~~ Methods for function \code{initialize} ~~
#' 
#' 
#' @name initialize-methods
#' @aliases initialize-methods initialize,Data-method
#' initialize,Fleet-method initialize,MSE-method initialize,Obs-method
#' initialize,OM-method initialize,Stock-method initialize
#' initialize,Fease-method initialize,DLM_general-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list('signature(.Object = \'DLM\')')}{ %% ~~describe this method
#' here~~ }
#' 
#' \item{list('signature(.Object = \'Fleet\')')}{ %% ~~describe this method
#' here~~ }
#' 
#' \item{list('signature(.Object = \'MSE\')')}{ %% ~~describe this method
#' here~~ }
#' 
#' \item{list('signature(.Object = \'Obs\')')}{ %% ~~describe this
#' method here~~ }
#' 
#' \item{list('signature(.Object = \'OM\')')}{ %% ~~describe this method here~~
#' }
#' 
#' \item{list('signature(.Object = \'Stock\')')}{ %% ~~describe this method
#' here~~ }
#' 
#' 
#' \item{list('signature(.Object = \'Fease\')')}{ %% ~~describe this method
#' here~~ } \item{list('signature(.Object = \'DLM_general\')')}{ %% ~~describe
#' this method here~~ }
#' 
#' }
#' @keywords methods ~~ other possible keyword(s) ~~
NULL


# ---- Obs Class ----
#' Class \code{'Obs'}
#' 
#' An operating model component that controls the observation model
#' 
#' 
#' @name Obs-class
#' @docType class
#' @note Its questionable whether the hyperstability/hyperdepletion should be
#' categorised as an observation model characteristic as it is most often
#' driven by fleet dynamics (and therefore should be in the fleet object). Oh
#' well its here and you might want to make it hyperstable beta < 1 or
#' hyperdeplete beta > 1, only.
#' 
#' @slot Name The name of the observation model object. Single value. Character string. 
#' 
#' @template Obs_template
#' 
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Obs')} 
#' 
#' @author T. Carruthers and A. Hordyk
#' @export
#' @keywords classes
#' @examples
#' 
#' showClass('Obs')
#' 
setClass("Obs", representation(Name = "character", 
                               Cobs = "numeric", Cbiascv = "numeric", CAA_nsamp = "numeric", CAA_ESS = "numeric", 
                               CAL_nsamp = "numeric", CAL_ESS = "numeric", 
                               Iobs = "numeric",  Ibiascv = "numeric", Btobs = "numeric", Btbiascv = "numeric", beta = "numeric",
                               LenMbiascv = "numeric", Mbiascv = "numeric", Kbiascv = "numeric",t0biascv = "numeric", Linfbiascv = "numeric",
                               LFCbiascv = "numeric", LFSbiascv = "numeric",
                               FMSYbiascv = "numeric", FMSY_Mbiascv = "numeric", BMSY_B0biascv = "numeric",
                               Irefbiascv = "numeric", Brefbiascv = "numeric", Crefbiascv = "numeric", 
                               Dbiascv = "numeric", Dobs = "numeric",
                               hbiascv = "numeric", Recbiascv = "numeric"))

# initialize Obs
setMethod("initialize", "Obs", function(.Object, file = NA, dec=c(".", ",")) {
  if (!is.na(file)) {
    if (file.exists(file)) {
      dec <- match.arg(dec)
      Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
      if (dec ==".") dat <- read.csv(file, header = F, colClasses = "character", 
                      col.names = paste0("V", 1:Ncol))  # read 1st sheet
      if (dec ==",") dat <- read.csv2(file, header = F, colClasses = "character", 
                                     col.names = paste0("V", 1:Ncol))  # read 1st sheet
      dname <- dat[, 1]
      dat <- dat[, 2:ncol(dat)]
      .Object@Name <- dat[match("Name", dname), 1]
      .Object@LenMbiascv <- as.numeric(dat[match("LenMbiascv", dname), 1])
      .Object@Cobs <- as.numeric(dat[match("Cobs", dname), 1:2])
      .Object@Cbiascv <- as.numeric(dat[match("Cbiascv", dname), 1])
      .Object@CAA_nsamp <- as.numeric(dat[match("CAA_nsamp", dname), 1:2])
      .Object@CAA_ESS <- as.numeric(dat[match("CAA_ESS", dname), 1:2])
      .Object@CAL_nsamp <- as.numeric(dat[match("CAL_nsamp", dname), 1:2])
      .Object@CAL_ESS <- as.numeric(dat[match("CAL_ESS", dname), 1:2])
      # .Object@CALcv <- as.numeric(dat[match("CALcv", dname), 1:2])
      .Object@Iobs <- as.numeric(dat[match("Iobs", dname), 1:2])
      .Object@Mbiascv <- as.numeric(dat[match("Mbiascv", dname), 1])
      .Object@Kbiascv <- as.numeric(dat[match("Kbiascv", dname), 1])
      .Object@t0biascv <- as.numeric(dat[match("t0biascv", dname), 1])
      .Object@Linfbiascv <- as.numeric(dat[match("Linfbiascv", dname), 1])
      .Object@LFCbiascv <- as.numeric(dat[match("LFCbiascv", dname), 1])
      .Object@LFSbiascv <- as.numeric(dat[match("LFSbiascv", dname), 1])
     # .Object@B0cv <- as.numeric(dat[match("B0cv", dname), 1])
      .Object@FMSYbiascv <- as.numeric(dat[match("FMSYbiascv", dname), 1])
      .Object@FMSY_Mbiascv <- as.numeric(dat[match("FMSY_Mbiascv", dname), 1])
      .Object@BMSY_B0biascv <- as.numeric(dat[match("BMSY_B0biascv", dname), 1])
     # .Object@rcv <- as.numeric(dat[match("rcv", dname), 1])
      .Object@Dbiascv <- as.numeric(dat[match("Dbiascv", dname), 1])
      .Object@Dobs <- as.numeric(dat[match("Dobs", dname), 1:2])
      .Object@Btbiascv <- as.numeric(dat[match("Btbiascv", dname), 1:2])
      .Object@Btobs <- as.numeric(dat[match("Btobs", dname), 1:2])
      # .Object@Fcurbiascv <- as.numeric(dat[match("Fcurbiascv", dname), 1])
      # .Object@Fcurcv <- as.numeric(dat[match("Fcurcv", dname), 1:2])
      .Object@hbiascv <- as.numeric(dat[match("hbiascv", dname), 1])
      .Object@Ibiascv <- as.numeric(dat[match("Ibiascv", dname), 1])
      # .Object@maxagecv <- as.numeric(dat[match("maxagecv", dname), 1])
      .Object@Recbiascv <- as.numeric(dat[match("Recbiascv", dname), 1:2])
      .Object@Irefbiascv <- as.numeric(dat[match("Irefbiascv", dname), 1])
      .Object@Crefbiascv <- as.numeric(dat[match("Crefbiascv", dname), 1])
      .Object@Brefbiascv <- as.numeric(dat[match("Brefbiascv", dname), 1])
      .Object@beta <- as.numeric(dat[match("beta", dname), 1:2])
    } else {
      message("File doesn't exist")
    }
  }
  .Object
  
})




# ---- Imp Class ----
#' Class \code{'Imp'}
#' 
#' An operating model component that specifies the degree of adherence to management recommendations (Implementation error)
#'  
#' @name Imp-class
#' @docType class
#' @slot Name The name of the Implementation error object. Single value. Character string.
#' 
#' @template Imp_template
#' 
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Imp')}#' 
#'      
#' @author T. Carruthers and A. Hordyk
#' @export
#' @keywords classes
#' @examples
#' 
#' showClass('Imp')
#' 
setClass("Imp", representation(Name = "character", TACFrac = "numeric", TACSD = "numeric", 
                               TAEFrac = "numeric", TAESD = "numeric", 
                               SizeLimFrac="numeric", SizeLimSD = "numeric"))

# initialize Imp
setMethod("initialize", "Imp", function(.Object, file = NA, dec=c(".", ",")) {
  
  .Object@Name <- "Perfect implementation"
  .Object@TACSD <- c(0,0)
  .Object@TACFrac <- c(1,1)
  .Object@TAESD <- c(0,0)
  .Object@TAEFrac <-c(1,1)
  .Object@SizeLimSD <- c(0,0)
  .Object@SizeLimFrac<-c(1,1)
  # .Object@Source <-"DLMtool generated"
  
  if (!is.na(file)) {
    if (file.exists(file)) {
      dec <- match.arg(dec)
      Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
      if (dec ==".") dat <- read.csv(file, header = F, colClasses = "character", 
                      col.names = paste0("V", 1:Ncol))  # read 1st sheet
      if (dec ==",") dat <- read.csv2(file, header = F, colClasses = "character", 
                                     col.names = paste0("V", 1:Ncol))  # read 1st sheet
      dname <- dat[, 1]
      dat <- dat[, 2:ncol(dat)]
      
      .Object@Name <- dat[match("Name", dname), 1]
      .Object@TACSD <- as.numeric(dat[match("TACSD", dname), 1:2])
      .Object@TACFrac <- as.numeric(dat[match("TACFrac", dname), 1:2])
      .Object@TAESD <- as.numeric(dat[match("TAESD", dname), 1:2])
      .Object@TAEFrac <- as.numeric(dat[match("TAEFrac", dname), 1:2])
      .Object@SizeLimSD <- as.numeric(dat[match("SizeLimSD", dname), 1:2])
      .Object@SizeLimFrac <- as.numeric(dat[match("SizeLimFrac", dname), 1:2])
      # .Object@Source <- dat[match("Source", dname), 1]
      
    } else {
      
      message("File doesn't exist")
      
    }
    
    
  }
  .Object
  
})



# ---- OM Class ----
#' Class \code{'OM'}
#' 
#' An object containing all the parameters needed to control the MSE which can
#' be build from component Stock, Fleet, Obs, and Imp objects. 
#' 
#' Almost all of these inputs are a vector of length 2 which describes the upper and lower
#' bounds of a uniform distribution from which to sample the parameter.
#' 
#' @name OM-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('OM', Stock, Fleet, Obs, Imp)}. 

#' @slot Name Name of the operating model
#' @slot Agency Name of the agency responsible for the management of the fishery. Character string
#' @slot Region Name of the general geographic region of the fishery. Character string
#' @slot Sponsor Name of the organization who sponsored the OM. Character string
#' @slot Latitude Latitude (decimal degrees). Negative values represent the South of the Equator. Numeric. Single value 
#' @slot Longitude Longitude (decimal degrees). Negative values represent the West of the Prime Meridian. Numeric. Single value 

#' @slot nsim The number of simulations
#' @slot proyears The number of projected years
#' @slot interval The assessment interval - how often would you like to update the management system?
#' @slot pstar The percentile of the sample of the management recommendation for each method
#' @slot maxF Maximum instantaneous fishing mortality rate that may be simulated for any given age class
#' @slot reps Number of samples of the management recommendation for each method. Note that when this is set to 1, the mean value of 
#' the data inputs is used. 
#' @slot cpars A list of custom parameters. Time series are a matrix nsim rows by nyears columns. Single parameters are a vector nsim long
#' @slot seed A random seed to ensure users can reproduce results exactly
#' @slot Source A reference to a website or article from which parameters were taken to define the operating model 
#' 
#' @template Stock_template
#' @template Fleet_template
#' @template Obs_template
#' @template Imp_template
#' 
#' @author T. Carruthers and A. Hordyk
#' @export
#' @keywords classes
#' 
setClass("OM", representation(Name = "character", Agency="character",
                              Region="character", Sponsor="character",
                              Latitude="numeric", Longitude="numeric",
                              nsim="numeric", proyears="numeric", 
                              interval='numeric', pstar='numeric', maxF='numeric', reps='numeric',
                              cpars="list",seed="numeric", Source="character"), 
         contains=c("Stock", "Fleet", "Obs", "Imp"))
# initialize OM
setMethod("initialize", "OM", function(.Object, Stock=NULL, Fleet=DLMtool::Generic_Fleet, 
                                       Obs=DLMtool::Generic_Obs, Imp=DLMtool::Perfect_Imp,
                                       interval=4, pstar=0.5, maxF=0.8, reps=1, nsim=48, proyears=50) {
  if (is.null(Stock)) {
    message("No Stock object found. Returning a blank OM object") 
    .Object@seed <- 1
    return(.Object)
  }
  
  if (class(Stock) != "Stock") 
    print(paste("Could not build operating model:", deparse(substitute(Stock)), "not of class Stock"))
  if (class(Fleet) != "Fleet") 
    print(paste("Could not build operating model:", deparse(substitute(Fleet)), "not of class Fleet"))
  if (class(Obs) != "Obs") 
    print(paste("Could not build operating model:", deparse(substitute(Obs)), "not of class Obs"))
  if (class(Imp) != "Imp") 
    print(paste("Could not build operating model:", deparse(substitute(Imp)), "not of class Imp"))
  if (class(Stock) != "Stock" | class(Fleet) != "Fleet" | class(Obs) != "Obs"  | class(Imp) != "Imp") stop()
  .Object@Name <- paste("Stock:", Stock@Name, "  Fleet:", Fleet@Name, "  Obs model:", 
                        Obs@Name, "  Imp model:", Imp@Name, sep = "")
  
  # Now copy the values for stock, fleet and observation slots to same
  # slots in the Sim object
  Sslots <- slotNames(Stock)
  for (i in 2:length(Sslots)) {
    tt <- .hasSlot(Stock, Sslots[i])  # For back-compatibility
    if (tt) slot(.Object, Sslots[i]) <- slot(Stock, Sslots[i])
  }
  Fslots <- slotNames(Fleet)
  for (i in 2:length(Fslots)) {
    tt <- .hasSlot(Fleet, Fslots[i])
    if (tt) slot(.Object, Fslots[i]) <- slot(Fleet, Fslots[i])
  }
  Oslots <- slotNames(Obs)
  for (i in 2:length(Oslots)) {
    tt <- .hasSlot(Obs, Oslots[i])
    if (tt) slot(.Object, Oslots[i]) <- slot(Obs, Oslots[i])
  }
  Islots <- slotNames(Imp)
  for (i in 2:length(Islots)) {
    tt <- .hasSlot(Imp, Islots[i])
    if (tt) slot(.Object, Islots[i]) <- slot(Imp, Islots[i])
  }
  
  # 
  # source <- paste("Stock:", Stock@Source, "Fleet:", Fleet@Source, "Obs:", Obs@Source, "Imp:",Imp@Source)
  slot(.Object, "Source") <- Stock@Source
  
  # Default MSE parameters
  if (.hasSlot(.Object, "nsim")) .Object@nsim <- nsim
  if (.hasSlot(.Object, "proyears")) .Object@proyears <- proyears
  
  # interval, pstar, maxF, reps
  if (.hasSlot(.Object, "interval")) .Object@interval <- interval
  if (.hasSlot(.Object, "pstar")) .Object@pstar <- pstar
  if (.hasSlot(.Object, "maxF")) .Object@maxF <- maxF
  if (.hasSlot(.Object, "reps")) .Object@reps <- reps
  
  if(length(.Object@Mexp) < 2) .Object@Mexp <- c(0,0)
  if(length(.Object@LenCV) < 2) .Object@LenCV <- c(0.08,0.15)
  if(length(.Object@CurrentYr)==0).Object@CurrentYr=.Object@nyears
  
  # if(length(.Object@FecB) < 2) .Object@FecB <- c(3,3)
  # if(all(is.na(.Object@FecB))) .Object@FecB <- c(3,3)  
  if(all(is.na(.Object@Mexp))) .Object@Mexp <- c(0,0)
  if(all(is.na(.Object@LenCV))) .Object@LenCV <- c(0.08,0.15)
  if(all(is.na(.Object@CurrentYr))) .Object@CurrentYr=.Object@nyears
  
  if(length(.Object@LR5) < 2) .Object@LR5 <- c(0,0)
  if(length(.Object@LFR) < 2) .Object@LFR <- c(0,0)
  if(length(.Object@Rmaxlen) < 2) .Object@Rmaxlen <- c(1,1)
  if(length(.Object@Fdisc) < 2) .Object@Fdisc <- c(0,0)
  
  if(all(is.na(.Object@LR5))) .Object@LR5 <- c(0,0)  
  if(all(is.na(.Object@LFR))) .Object@LFR <- c(0,0)  
  if(all(is.na(.Object@Rmaxlen))) .Object@Rmaxlen <- c(1,1)
  if(all(is.na(.Object@Fdisc))) .Object@Fdisc <- c(0,0)  
  
  if (.hasSlot(.Object, "Size_area_1")) {
    if (length(.Object@Size_area_1)==0) .Object@Size_area_1 <- .Object@Frac_area_1
    if (all(is.na(.Object@Size_area_1))) .Object@Size_area_1 <- .Object@Frac_area_1
  } 
  
  .Object@seed=1
  .Object
})


# ---- MSE Class ----
#' Class \code{'MSE'}
#' 
#' A Management Strategy Evaluation object that contains information about
#' simulation conditions and performance of data-limited methods
#' 
#' 
#' @name MSE-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('MSE', Name, nyears, proyears, nMPs, MPs, nsim, OMtable, Obs,
#' B_BMSYa, F_FMSYa, Ba, FMa, Ca, OFLa, Effort, PAA, CAA, CAL, CALbins)} 
#'
#' @slot Name Name of the MSE object. Single value. Character string
#' @slot nyears The number of years for the historical simulation. Single value. Positive integer
#' @slot proyears The number of years for the projections - closed loop simulations. Single value. Positive integer
#' @slot nMPs Number of management procedures simulation tested. Single value. Positive integer. 
#' @slot MPs The names of the MPs that were tested. Vector of length nMPs. Character strings. 
#' @slot nsim Number of simulations. Single value. Positive integer
#' 
#' @template OM_desc 
#' @template Obs_desc 
#' 
#' @slot B_BMSY Simulated spawning biomass relative to spawning BMSY over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers 
#' @slot F_FMSY Simulated fishing mortality rate relative to FMSY over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers
#' @slot B Simulated stock biomass over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers 
#' @slot SSB Simulated spawning stock biomass over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers
#' @slot VB Simulated vulnerable biomass over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers
#' @slot FM Simulated fishing mortality rate over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers
#' @slot C Simulated catches (taken) over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers
#' @slot TAC Simulated Total Allowable Catch (prescribed) over the projection (this is NA for input controls). An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers 
#' @slot SSB_hist Simulated historical spawning stock biomass. An array with dimensions: nsim, nages, nyears, nareas. Non-negative real numbers
#' @slot CB_hist Simulated historical catches in weight. An array with dimensions: nsim, nages, nyears, nareas. Non-negative real numbers
#' @slot FM_hist Simulated historical fishing mortality rate. An array with dimensions: nsim, nages, nyears, nareas. Non-negative real numbers
#' @slot Effort Simulated relative fishing effort in the projection years. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers
#' @slot PAA Population at age in last projection year. An array with dimensions: nsim, nMPs, nages. Non-negative real numbers
#' @slot CAA Catch at age in last projection year. An array with dimensions: nsim, nMPs, nages. Non-negative real numbers
#' @slot CAL Catch at length in last projection year. An array with dimensions: nsim, nMPs, nCALbins. Non-negative real numbers
#' @slot CALbins Mid-points of the catch-at-length bins. Vector of length nCALbins. Positive real numbers. 
#' @slot Misc Miscellanenous output such as posterior predictive data
#'
#' @author T. Carruthers
#' @keywords classes
setClass("MSE", representation(Name = "character", nyears = "numeric", 
                               proyears = "numeric", nMPs = "numeric", MPs = "character", nsim = "numeric", 
                               OM = "data.frame", Obs = "data.frame", B_BMSY = "array", F_FMSY = "array", 
                               B = "array", SSB="array", VB="array", FM = "array", C = "array", 
                               TAC = "array", SSB_hist = "array", 
                               CB_hist = "array", FM_hist = "array", Effort = "array", PAA= "array", CAA= "array", 
                               CAL= "array", CALbins="numeric", Misc="list"))


setMethod("initialize", "MSE", function(.Object, Name, nyears, proyears, 
                                        nMPs, MPs, nsim, OM, Obs, B_BMSY, F_FMSY, B, SSB, VB, FM, C, TAC, 
                                        SSB_hist, CB_hist, FM_hist, Effort = array(), PAA,  CAA, CAL, CALbins, Misc) {
  .Object@Name <- Name
  .Object@nyears <- nyears
  .Object@proyears <- proyears
  .Object@nMPs <- nMPs
  .Object@MPs <- MPs
  .Object@nsim <- nsim
  .Object@OM <- OM
  .Object@Obs <- Obs
  .Object@B_BMSY <- B_BMSY
  .Object@F_FMSY <- F_FMSY
  .Object@B <- B
  .Object@SSB <- SSB
  .Object@VB <- VB
  .Object@FM <- FM
  .Object@C <- C  
  .Object@TAC <- TAC
  .Object@SSB_hist <- SSB_hist
  .Object@CB_hist <- CB_hist
  .Object@FM_hist <- FM_hist
  .Object@Effort <- Effort
  .Object@PAA <- PAA
  .Object@CAA <- CAA
  .Object@CAL <- CAL
  .Object@CALbins <- CALbins
  .Object@Misc <- Misc
  
  .Object
})

# setMethod("initialize", "MSE", function(.Object, Name, nyears, proyears, 
# nMPs, MPs, nsim, OMtable, Obs, B_BMSYa, F_FMSYa, Ba, FMa, Ca, TACa, 
# SSB_hist, CB_hist, FM_hist, Effort = array()) {
# .Object@Name <- Name
# .Object@nyears <- nyears
# .Object@proyears <- proyears
# .Object@nMPs <- nMPs
# .Object@MPs <- MPs
# .Object@nsim <- nsim
# .Object@OM <- OMtable
# .Object@Obs <- Obs
# .Object@B_BMSY <- B_BMSYa
# .Object@F_FMSY <- F_FMSYa
# .Object@B <- Ba
# .Object@FM <- FMa
# .Object@C <- Ca
# .Object@TAC <- TACa
# .Object@SSB_hist <- SSB_hist
# .Object@CB_hist <- CB_hist
# .Object@FM_hist <- FM_hist
# .Object@Effort <- Effort
# .Object
# })


# #' Example data object
# #' 
# #' Example data object with a number of output control MPs run on it, and
# #' includes resulting distributions of TACs
# #' 
# #' 
# #' @name ourReefFish
# #' @docType data
# #' @usage data('ourReefFish')
# #' @keywords datasets
# #' @examples
# #' \dontrun{ 
# #' data(ourReefFish)
# #' str(ourReefFish)  
# #' plot(ourReefFish) 
# #' }
# #'
# NULL


# ---- PMobj Class ----

#' An object for storing data for analysis using data-limited methods
#' 
#' Used interally
#' 
#' @name PMobj-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('PMobj')} 
#' @slot Name Name of the Performance Metric. Character 
#' @slot Caption A caption to be used in plots. Character, call, or function.
#' @slot Stat Statistic of interest for the PM. Dimensions: nsim, nMP, yrs. Array 
#' @slot Ref Reference value to calculate probability for statistic. Numeric.
#' @slot Prob Probability (mean over years) Dimensions: nsim by MP.  Matrix, numeric or data.frame  
#' @slot Mean Mean probability (mean over years and simulations). Numeric. Length nMPs 
#' @slot MPs Name of MPs. Single value. Character string  
#' @author  A. Hordyk
#' @keywords classes

setClass("PMobj", representation(Name = "character",  Caption='label.class', 
                                 Stat='array', Ref='numeric', Prob='prob.class', Mean='numeric',
                                 MPs="character"))



show <- function(object) methods::show(object)


#' Show the output of a PM
#'
#' @param object object of class MSE
#' @rdname show-MSE
#' @export
setMethod("show", signature = (object="PMobj"), function(object) {
  
  sls <- slotNames(object)
  df <- data.frame(Slot=sls, Value=NA, stringsAsFactors = FALSE)
  if (any(sapply(sls, function(sl) length(slot(object, sl))) == 0)) {
    cat("Incomplete PMobj\n")
    for (sl in sls) {
      r <- match(sl, sls)
      slval <- slot(object, sl)
      if (class(slval) == "array" & length(slval)>0) {
        df[r,2] <- 'array'
      } else if (class(slval) == "matrix" & length(slval)>0) {
        df[r,2] <- 'matrix'
      } else if (length(slval) > 0 & class(slval) != "call") {
        df[r,2] <- slval
      } else {
        df[r, 2] <- 'not defined'
      }
      
    }
    print(df)
  } else {
    cat(object@Name)
    cat("\n", object@Caption)
    cat("\n")
    
    nMP <- length(object@MPs)
    if (nMP > 1) nsim <- dim(object@Prob)[1]
    if (nMP == 1) nsim <- length(object@Prob)
    
    nprint <- min(nsim, 10)
    if (nMP > 1) df <- data.frame(object@Prob[1:nprint,])
    if (nMP == 1) df <- data.frame(object@Prob[1:nprint])
    if (nsim > (nprint+1)) {
      if (nMP > 1) lst <- object@Prob[nprint+1,]
      if (nMP == 1) lst <- object@Prob[nprint+1]
    } else {
      if (nMP > 1) lst <- object@Prob[nprint,]
      if (nMP == 1) lst <- object@Prob[nprint]
    }
    df <- round(df,2)
    lst <- round(lst,2)
    colnames(df) <- object@MPs
    names(lst) <- object@MPs
    if (nsim > (nprint+1)) {
      df <- rbind(df,
                  rep(".", nMP),
                  rep(".", nMP),
                  rep(".", nMP),
                  lst)
      rownames(df) <- c(1:(nprint+3), nsim)
    }
    print(df)
    
    cat("\nMean\n")
    print(round(object@Mean,2))  
  
  }
})


# ---- Summary of MSE object ----

#' Summary of MSE object
#'
#' @param object object of class MSE
#' @param ... a list of names of PM methods
#' @param silent Should summary be printed to console? Logical.
#' @param Refs An optional named list (matching the PM names) with numeric values to override the default `Ref` values. See examples.
#' @rdname summary-MSE
#' @export
setMethod('summary', signature="MSE", function(object, ..., silent=FALSE, Refs=NULL) {
  PMlist <- unlist(list(...))
  
  if(length(PMlist) == 0) PMlist <- c("PNOF", "P50", "AAVY", "LTY")
  if (class(PMlist) != 'character') stop("Must provide names of PM methods")
  # check
  for (X in seq_along(PMlist)) 
    if (class(get(PMlist[X])) != "PM") stop(PMlist[X], " is not a valid PM method")
    
  if (!silent) message("Calculating Performance Metrics")
  storeMean <- vector('list', length(PMlist))
  storeName <- vector('list', length(PMlist))
  storeCap <- vector('list', length(PMlist))
  storeHeading <- vector('list', length(PMlist))
  storeMP <- vector('list', length(PMlist))
  for (X in 1:length(PMlist)) {
    ref <- Refs[[PMlist[X]]]
    if (is.null(ref)) {
      runPM <- eval(call(PMlist[X], object))  
    } else {
      runPM <- eval(call(PMlist[X], object, Ref=ref))
    }
    storeMean[[X]] <- runPM@Mean
    storeName[[X]] <- runPM@Name
    storeCap[[X]] <- runPM@Caption
    storeMP[[X]] <- runPM@MPs
  }
  
  df <- data.frame('MP'=storeMP[[1]], signif(do.call('cbind', storeMean),2), stringsAsFactors = FALSE)
  # heading <- do.call('rbind', storeHeading)
  colnames(df)[2:(length(PMlist)+1)] <- PMlist #caps # gsub(" ", "", caps)
  if (!silent) {
    dfprint <- data.frame('Performance Metrics' = do.call('rbind', storeName), gap="", do.call('rbind', storeCap))
    names(dfprint)[2:3] <- ''
    print(dfprint)
    cat("\n")
    cat("\nPerformance Statistics:\n")
    print(df)  
  }
  
  invisible(df)
  
})


# # ---- Plot Data Object -----
# #' Plot Data object
# #'
# #' @rdname plot-Data 
# #' @param x object of class Data
# #' @param upq Upper quantile of TACs for max ylim
# #' @param lwq Lower quantile of TACs for min ylim
# #' @param outline Logical. Include outliers in plot?
# #' @param ...  Optional additional arguments passed to \code{boxplot}
# #' @export
# setMethod("plot",
#           signature(x = "Data"),
#           function(x, upq=0.9, lwq=0.1, outline = FALSE, ...){
#             
#             old_par <- par(no.readonly = TRUE)
#             on.exit(par(list = old_par), add = TRUE)
#             boxplot.Data(x, upq, lwq, outline, ...)
#           })
#             # Data<-x
            # if (class(Data) != "Data") stop("Must supply object of class Data")
            # if (all(is.na(Data@TAC))) stop("No TAC data found")
            # cols<-rep(c('black','red','green','blue','orange','brown','purple','dark grey','violet','dark red','pink','dark blue','grey'),4)
            # ltys<-rep(1:4,each=13)
            # 
            # if(is.na(funcs[1]))funcs<-Data@MPs
            # 
            # nMPs<-length(funcs)
            # nplots<-ceiling(nMPs/maxlines)
            # maxl<-ceiling(nMPs/nplots)
            # mbyp <- split(1:nMPs, ceiling(1:nMPs/maxl))   # assign methods to plots
            # 
            # if(is.na(xlims[1])|length(xlims)!=2){
            #   xlims<-quantile(Data@TAC,c(0.005,0.95),na.rm=T)
            #   if(xlims[1]<0)xlims[1]<-0
            # }
            # if(!NAor0(Data@Ref)){
            #   if(xlims[1]>Data@Ref)xlims[1]<-max(0,0.98*Data@Ref)
            #   if(xlims[2]<Data@Ref)xlims[2]<-1.02*Data@Ref
            # }
            # ylims<-c(0,1)
            # 
            # #for(m in 1:nMPs){
            # # if(sum(!is.na(Data@TAC[m,,1]))>2){
            # # dens<-density(Data@TAC[m,,1],na.rm=T)
            # #print(quantile(dens$y,0.99,na.rm=T))
            # #  if(quantile(dens$y,0.9,na.rm=T)>ylims[2])ylims[2]<-quantile(dens$y,0.90,na.rm=T)
            # #}
            # #}
            # 
            # #dev.new2(width=10,height=0.5+7*nplots)
            # par(mfrow=c(ceiling(nplots/2),2),mai=c(0.4,0.4,0.01,0.01),omi=c(0.35,0.35,0.35,0.05))
            # 
            # for(p in 1:nplots){
            #   m<-mbyp[[p]][1]
            #   plot(NA,NA,xlim=xlims,ylim=ylims,main="",xlab="",ylab="",col="white",lwd=3,type="l")
            #   abline(h=0)
            #   if(!NAor0(Data@Ref)){
            #     abline(v=Data@Ref,col="light grey",lwd=2)
            #     if(!NAor0(Data@Ref_type[1]))legend('right',Data@Ref_type,text.col="grey",bty='n')
            #   }
            #   #plot(density(DLM@TAC[m,,1],from=0,na.rm=T),xlim=xlims,ylim=ylims,main="",xlab="",ylab="",col=coly[m],lty=ltyy[m],type="l")
            #   
            #   if(!is.na(perc[1]))abline(v=quantile(Data@TAC[m,,1],p=perc,na.rm=T),col=cols[m],lty=ltys[m])
            #   #if(length(mbyp[[p]])>0){
            #   for(ll in 1:length(mbyp[[p]])){
            #     m<-mbyp[[p]][ll]
            #     if(sum(!is.na(Data@TAC[m,,1]))>10){  # only plot if there are sufficient non-NA TAC samples
            #       x<-density(Data@TAC[m,,1],from=0,na.rm=T)$x
            #       y<-density(Data@TAC[m,,1],from=0,na.rm=T)$y
            #       y<-y/max(y)
            #       lines(x,y,col=cols[ll])
            #     }else{
            #       print(paste("Method ",funcs[m]," produced too many NA TAC values for plotting densities",sep=""))
            #     }
            #     if(!is.na(perc[1]))abline(v=quantile(Data@TAC[m,,1],p=perc,na.rm=T),col=cols[ll],lty=2)
            #   }
            #   #}
            #   cind<-1:length(mbyp[[p]])
            #   legend('topright',funcs[mbyp[[p]]],text.col=cols[cind],col=cols[cind],lty=1,bty='n',cex=0.75)
            # }
            # 
            # mtext(paste("TAC (",Data@Units,")",sep=""),1,outer=T,line=0.5)
            # mtext(paste("Standardized relative frequency",sep=""),2,outer=T,line=0.5)
            # mtext(paste("TAC calculation for ",Data@Name,sep=""),3,outer=T,line=0.5)
          # })

# # ---- Plot MSE object ----
# #' Plot MSE object
# #'
# #' @rdname plot-MSE
# #' @param x object of class MSE
# #' @export
# setMethod("plot",
#           signature(x = "MSE"),
#           function(x){
#             MSEobj<-x
#             Pplot(MSEobj)
#             Kplot(MSEobj)
#             Tplot(MSEobj)
#           })





# ---- Summary of Data Object ----
#' Summary of Data object
#'
#' @rdname summary-Data
#' @param object An object of class Data
#' @param wait Logical. Wait for key press before next plot?
#' @param x iteration number for the Data object.
#' @param plots Character. What plots to show? `all`, `TS`, `CAA`, `CAL`, `PD` 
#' for all plots, time-series, catch-at-age, catch-at-length, and 
#' probability distributions respectively
#' @param rmd Logical. Used in a rmd file?
#' @param head Character. Heading for rmd file. Default is '##' (second level heading)
#' @param tplot Integer. Number of plots per page. Default 25
#' @export
setMethod("summary",
          signature(object = "Data"),
          function(object, wait=TRUE, x=1, plots='all', rmd=FALSE, head="##", tplot=25){
            plots <- match.arg(plots, c('all', 'TS', 'CAA', 'CAL', 'PD'), several.ok = TRUE)
            if ('all' %in% plots) plots <- c('TS', 'CAA', 'CAL', 'PD')
            
            Freq <- n <- Var2 <- NULL # cran check
            if (class(object) != "Data") stop("Object must be class `Data`", call.=FALSE)
            
            # Time-Series
            Year <- object@Year
            l <- length(Year)
            
            slts <- c("Cat", 'Ind', 'SpInd', 'VInd', 'Rec', 'ML', 'Lc')
            Cat <- Ind <- SpInd <- VInd <- Rec <- ML <- Lc <- NULL # cran checks
            for (sl in slts) {
              tt <- slot(object, sl)[x,]
              if (length(tt)!=l) tt <- c(tt, rep(NA,l-length(tt)))
              assign(sl, tt)
            }
            Val <- c(Cat, Ind, SpInd, VInd, Rec, ML, Lc)
            Var <- rep(c("Catch", "Total Index", "Spawning Index",
                         "Vuln. Index", "Recruitment", "Mean Length", "Mean Length above Lc"), each=length(Year))
            ts.df <- data.frame(Year=Year, Val=Val, Var=Var, stringsAsFactors = TRUE)
            # ts.df$Year <- as.factor(ts.df$Year)
            ts.df$Var <- factor(ts.df$Var, levels=
                                  c("Catch", "Total Index", "Spawning Index",
                                    "Vuln. Index", "Recruitment", "Mean Length", "Mean Length above Lc"))
            ts.df <- subset(ts.df, !is.na(Val))
            if (nrow(ts.df)>0 && 'TS' %in% plots) {
              P1 <- ggplot2::ggplot(ts.df, ggplot2::aes(x=Year, y=Val, group = Var)) +
                ggplot2::facet_wrap(~Var, scales='free_y') + ggplot2::geom_line(size=1.25) +
                ggplot2::theme_classic() + 
                ggplot2::expand_limits(y=0) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
                ggplot2::scale_x_continuous(breaks=pretty(rev(Year), length(Year)/5)) +
                ggplot2::labs(y="")
              
            } else {
              P1 <- NULL
            }
            if (!is.null(P1)) {
              if (rmd) {
                cat('\n')
                cat('\n')
                cat(head, 'Time-Series')
                cat('\n')
              } else { 
                message('Plotting Time-Series')
              }
              print(P1)
            }
            
            if (interactive() & wait & !is.null(P1)) 
              invisible(readline(prompt="Press [enter] to continue..."))
            
            # CAA 
            if (all(is.na(object@CAA))) {
              P2 <- FALSE
            } else {
              P2 <- TRUE
            }
          
            if (P2 == TRUE) {
              CAA <- object@CAA[x,,]
              nyrs <- nrow(CAA); maxage <- ncol(CAA)
    
              dimnames(CAA) <- list(1:nyrs, 1:maxage)
              
              df1 <- as.data.frame.table(CAA, stringsAsFactors = FALSE)
              colnames(df1) <- c("Year", "Val", "Freq")
              df1$Val <- as.numeric(df1$Val)
              
              df1$Year <- as.numeric(df1$Year)

              yr.n <- df1 %>% dplyr::group_by(Year) %>% dplyr::summarise(n=sum(Freq))
              yr.ind <- yr.n %>% dplyr::filter(n>0) %>% dplyr::select(Year)
              
              Years <- object@Year
              nyears <- length(unique(df1$Year))
              df1$Year_val <- (Years[(length(Years)-nyears+1):length(Years)])
              df1 <- df1[!is.na(df1$Freq),]  # Fliter out NA values, so we don't try to plot missing years
            
              if (nrow(df1)>0 && 'CAA' %in% plots) {
                
                if (rmd) {
                  cat('\n')
                  cat('\n')
                  cat(paste(head, 'Catch-at-Age'))
                  cat('\n')
                } else {
                  message('Plotting Catch-at-Age')
                }
                
                nyears <- length(unique(df1$Year))
                nbins <- length(unique(df1$Val))
                if (nyears > tplot) {
                  npages <- ceiling(nyears/tplot) 
                  ncol <- 5 
                  nrow <- 5
                  nplot <- ncol * nrow
                } else {
                  npages <- 1
                  nrow <- ceiling(sqrt(nyears))
                  ncol <- ceiling(nyears/nrow)
                  nplot <- nyears
                }
                pmat <- matrix(1:(nrow*ncol), nrow=nrow, ncol=ncol, byrow=TRUE)
                pmat[pmat >nplot] <- NA
                
                op <- par(mfrow=c(nrow, ncol), no.readonly = TRUE, mar=c(2,2,2,1), oma=c(4,4,2,0))
                on.exit(par(op))
                
                yr1 <- 1 
                col <- "grey"
                for (pg in 1:npages) {
                  if(npages>1)message('Plot ', pg, ' of ', npages)
                  yrind <- unique(df1$Year)[yr1:(yr1+nplot-1)]
                  yr1 <- max(yrind) + 1
                  dat <- df1 %>% dplyr::filter(Year %in% yrind)
                  un.yrs_val <- as.numeric(unique(dat$Year_val))
                  un.yrs <- as.numeric(unique(dat$Year))
                  
                  if (pg >1 && pg == npages) {
                    nplot <- nyears - (npages-1) * tplot
                    ncol <- ceiling(sqrt(nplot))
                    nrow <- ceiling(nplot/ncol)
                    op <- par(mfrow=c(nrow, ncol), no.readonly = TRUE, mar=c(2,2,2,1), oma=c(4,4,2,0))
                    on.exit(par(op))
                    pmat <- matrix(1:nplot, nrow=nrow, ncol=ncol, byrow=TRUE)
                  }
                  for (p in 1:nplot) {
                    pdat <- dat %>% dplyr::filter(Year==un.yrs[p])
                    if (nrow(pdat) > 0) {
                      if (p %in% pmat[nrow,]) {
                        barplot(pdat$Freq, names=round(pdat$Val, 2), axes=FALSE, col=col)  
                      } else {
                        barplot(pdat$Freq, names=FALSE, axes=FALSE, col=col)
                      } 
                      if (p %in% pmat[,1]) axis(side=2)
                      if (!p %in% pmat[,1]) axis(side=2, labels=TRUE)
                      ncount <- round(sum(pdat$Freq),0)
                      title(un.yrs_val[p])
                      text(max(pdat$Val), max(pdat$Freq), paste('n = ', ncount),
                           xpd=NA)
                    }
                    
                  } 
                  mtext(side=1, outer=TRUE, "Age", line=2, cex=1.5)
                  mtext(side=2, outer=TRUE, "Frequency", line=2, cex=1.5)
                  
                }
              }
            }
            if (interactive() & wait & !is.null(P2)) 
              invisible(readline(prompt="Press [enter] to continue..."))
            
            
            # CAL 
            if (all(is.na(object@CAL))) {
              P3 <- FALSE
            } else {
              P3 <- TRUE
            }
            
            if (P3 == TRUE) {
              CAL <- object@CAL[x,,]
              nyrs <- nrow(CAL); nbins <- length(object@CAL_bins) - 1
              By <- object@CAL_bins[2] - object@CAL_bins[1]
              BinsMid <- seq(object@CAL_bins[1] + 0.5*By, by=By,length.out = nbins)
              dimnames(CAL) <- list(1:nyrs, BinsMid)
              
              df1 <- as.data.frame.table(CAL, stringsAsFactors = FALSE)
              colnames(df1) <- c("Year", "Val", "Freq")
              df1$Val <- as.numeric(df1$Val)
              
              df1$Year <- as.numeric(df1$Year)
  
              yr.n <- df1 %>% dplyr::group_by(Year) %>% dplyr::summarise(n=sum(Freq))
              yr.ind <- yr.n %>% dplyr::filter(n>0) %>% dplyr::select(Year)
              
              Years <- object@Year
              nyears <- length(unique(df1$Year))
              df1$Year_val <- (Years[(length(Years)-nyears+1):length(Years)])
  
              if (nrow(df1) > 0  && 'CAL' %in% plots) {
                
                if (rmd) {
                  cat('\n')
                  cat('\n')
                  cat(paste(head, 'Catch-at-Length'))
                  cat('\n')
                } else {
                  message('Plotting Catch-at-Length')
                }
                
                nyears <- length(unique(df1$Year))
                nayears <- df1 %>% group_by(Year) %>% summarize(isna=all(is.na(Freq)))
                nyears <- sum(!nayears$isna)
                
                nbins <- length(unique(df1$Val))
                if (nyears > tplot) {
                  npages <- ceiling(nyears/tplot) 
                  ncol <- 5 
                  nrow <- 5
                  nplot <- ncol * nrow
                } else {
                  npages <- 1
                  nrow <- ceiling(sqrt(nyears))
                  ncol <- ceiling(nyears/nrow)
                  nplot <- nyears
                }
                pmat <- matrix(1:(ncol*nrow), nrow=nrow, ncol=ncol, byrow=TRUE)
                pmat[pmat>nplot] <- NA
                
                op <- par(mfrow=c(nrow, ncol), no.readonly = TRUE, mar=c(2,2,2,1), oma=c(4,4,2,0))
                on.exit(par(op))
                
                yr1 <- 1
                col <- "grey"
                for (pg in 1:npages) {
                  if(npages>1)message('Plot ', pg, ' of ', npages)
                  if (pg ==1) {
                    yrind <- yr.ind$Year[1:(nplot)]  
                  } else {
                    t1 <- nplot * (pg-1) +1 
                    t2 <- t1+nplot
                    yrind <- yr.ind$Year[t1:t2]
                  }
                  
                  yr1 <- max(yrind) + 1
                  dat <- df1 %>% dplyr::filter(Year %in% yrind)
                  un.yrs_val <- as.numeric(unique(dat$Year_val))
                  un.yrs <- as.numeric(unique(dat$Year))
                  
                  if (pg >1 && pg == npages) {
                    nplot <- nyears - (npages-1) * tplot
                    ncol <- ceiling(sqrt(nplot))
                    nrow <- ceiling(nplot/ncol)
                    op <- par(mfrow=c(nrow, ncol), no.readonly = TRUE, mar=c(2,2,2,1), oma=c(4,4,2,0))
                    on.exit(par(op))
                    pmat <- matrix(1:(ncol*nrow), nrow=nrow, ncol=ncol, byrow=TRUE)
                    pmat[pmat>nplot] <- NA
                  }
                  for (p in 1:nplot) {
                    pdat <- dat %>% dplyr::filter(Year==un.yrs[p])
                    if (nrow(pdat) > 0) {
                      if (all(is.na(pdat$Freq))) {
                       
                      } else{
                        if (p %in% pmat[nrow,]) {
                          barplot(pdat$Freq, names.arg=round(pdat$Val, 2), axes=FALSE, col=col, las=2)  
                        } else {
                          barplot(pdat$Freq, names.arg=FALSE, axes=FALSE, col=col)
                        } 
                        if (p %in% pmat[,1]) axis(side=2)
                        if (!p %in% pmat[,1]) axis(side=2, labels=TRUE)
                        ncount <- round(sum(pdat$Freq),0)
                        title(un.yrs_val[p])
                        text(length(unique(df1$Val)), max(pdat$Freq), paste('n = ', ncount), xpd=NA)
                      }
                      
                    }
                    
                  } 
                  mtext(side=1, outer=TRUE, "Length", line=2, cex=1.5)
                  mtext(side=2, outer=TRUE, "Frequency", line=2, cex=1.5)
                  
                }
              }
            }
            
            if (interactive() & wait & !is.null(P3))
              invisible(readline(prompt="Press [enter] to continue..."))
            
            # Biology & Depletion
            slots<-c("Dep","Mort","FMSY_M","Dt","BMSY_B0","vbK", "vbLinf")
            namey<-c("Stock depletion", "Natural Mortality rate","Ratio of FMSY to M",
                     "Depletion over time t","BMSY relative to unfished",
                     "von B. k parameter", "von B. Linf parameter")
            slotsCV<-c("CV_Dep","CV_Mort","CV_FMSY_M","CV_Dt","CV_BMSY_B0","CV_vbK", "CV_vbLinf")
            
            reps <- 5000
            val <- list()
            for (i in seq_along(slots)) {
              mu <- slot(object, slots[i])[x]
              cv <- slot(object, slotsCV[i])[x]
              val[[i]] <- trlnorm(reps, mu,cv)
            }
            vals <- do.call("cbind", val)
            
            colnames(vals) <- namey
            
            df1 <- as.data.frame.table(vals, stringsAsFactors = TRUE)
            df1 <- df1 %>% dplyr::filter(is.na(Freq) == FALSE)
            if (nrow(df1) > 0  && 'PD' %in% plots) {
              P4 <-   ggplot2::ggplot(df1, ggplot2::aes(x=Freq, group=Var2)) +
                ggplot2::facet_wrap(~Var2, scales="free") + ggplot2::geom_histogram(bins=30) +
                ggplot2::labs(y="Frequency", x="Parameter Value") +
                ggplot2::theme_classic()
            } else {
              P4 <- NULL
            }
            
            if (!is.null(P4)) {
              if (rmd) {
                cat('\n')
                cat('\n')
                cat(paste(head, 'Parameter Distributions'))
                cat('\n')
              } else {
                message('Plotting Parameter Distributions')
              }
              print(P4)
            }
            
          })
            
           

# Summary of MSE object
#
# @param object object of class MSE
# @rdname summary-MSE
# @export
# # setMethod("summary",
# #           signature(object = "MSE"),
# summaryold <- function(object){            
#             
#             MSEobj<-object      
#             nm<-MSEobj@nMPs
#             nsim<-MSEobj@nsim
#             proyears<-MSEobj@proyears
#             
#             Yd<-P10<-P50<-P100<-POF<-LTY<-STY<-VY<-array(NA,c(nm,nsim))
#             
#             yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
#             RefYd<-MSEobj@OM$RefY
#             yend<-max(MSEobj@proyears-9,1):MSEobj@proyears
#             ystart<-1:10
#             y1<-1:(MSEobj@proyears-1)
#             y2<-2:MSEobj@proyears
#             
#             for(m in 1:nm){
#               Yd[m,]<-round(apply(MSEobj@C[,m,yind],1,mean,na.rm=T)/RefYd*100,1)
#               POF[m,]<-round(apply(MSEobj@F_FMSY[,m,]>1,1,sum,na.rm=T)/proyears*100,1)
#               P10[m,]<-round(apply(MSEobj@B_BMSY[,m,]<0.1,1,sum,na.rm=T)/proyears*100,1)
#               P50[m,]<-round(apply(MSEobj@B_BMSY[,m,]<0.5,1,sum,na.rm=T)/proyears*100,1)
#               P100[m,]<-round(apply(MSEobj@B_BMSY[,m,]<1,1,sum,na.rm=T)/proyears*100,1)
#               LTY[m]<-round(sum(MSEobj@C[,m,yend]/RefYd>0.5)/(MSEobj@nsim*length(yend))*100,1)
#               STY[m]<-round(sum(MSEobj@C[,m,ystart]/RefYd>0.5)/(MSEobj@nsim*length(ystart))*100,1)
#               AAVY<-apply(((MSEobj@C[,m,y1]-MSEobj@C[,m,y2])^2)^0.5,1,mean)/apply(MSEobj@C[,m,y2],1,mean)
#               VY[m]<-round(sum(AAVY<0.1)/MSEobj@nsim*100,1)
#             }
#             nr<-2
#             out<-cbind(MSEobj@MPs,round(apply(Yd,1,mean,na.rm=T),nr),round(apply(Yd,1,sd,na.rm=T),nr),
#                        round(apply(POF,1,mean,na.rm=T),nr),round(apply(POF,1,sd,na.rm=T),nr),
#                        round(apply(P10,1,mean,na.rm=T),nr),round(apply(P10,1,sd,na.rm=T),nr),
#                        round(apply(P50,1,mean,na.rm=T),nr),round(apply(P50,1,sd,na.rm=T),nr),
#                        round(apply(P100,1,mean,na.rm=T),nr),round(apply(P100,1,sd,na.rm=T),nr),
#                        round(apply(LTY,1,mean,na.rm=T),nr),
#                        round(apply(STY,1,mean,na.rm=T),nr),
#                        round(apply(VY,1,mean,na.rm=T),nr))
#             out<-as.data.frame(out)
#             names(out)<-c("MP","Yield","stdev","POF","stdev ","P10","stdev",
#                           "P50","stdev","P100","stdev","LTY","STY","VY")
#             out[,1]<-as.character(out[,1])
#             for(i in 2:ncol(out))out[,i]<-as.numeric(as.character(out[,i]))
#             out
#           }



# # -- Input Control Recommendation Class -
# #' Class \code{'InputRec'}
# #' 
# #' An object for storing the recommendation for an input control MP 
# #' 
# #' @name InputRec-class
# #' @docType class
# #' @section Objects from the Class: Objects can be created by calls of the form
# #' \code{new('InputRec')} 
# 
# #' @slot Effort A numeric value with the effort recommendation as a fraction of current (nyear) fishing effort
# #' @slot Spatial A boolean vector of length 'nareas' specifying if area is open (1) or closed (0) to fishing 
# #' @slot Allocate A boolean value describing if effort should be re-allocated from close to open areas
# #' @slot LR5 smallest length at 5 per cent retention
# #' @slot LFR smallest length at full retention 
# #' @slot HS upper harvest slot (no retention above this)
# #' @slot Rmaxlen retention of the largest size class
# #' @slot Misc An empty list that can be used to store information and pass on to MPs in future 
# #' @author A. Hordyk
# #' @keywords classes
# 
# setClass("InputRec", representation(Effort = "numeric", 
#                                     Spatial="numeric", Allocate = "numeric", LR5 = "numeric",
#                                     LFR = "numeric", HS="numeric", Rmaxlen="numeric", Misc="list"))
# setMethod("initialize", "InputRec", function(.Object){
#      .Object@Effort<-1
#      .Object@Allocate<-1
#      .Object@Spatial<-c(1,1)
#      .Object
#    })


# -- Management Recommendation Class ----
#' Class \code{'Rec'}
#' 
#' An object for storing the MP recommendations 
#' 
#' @name Rec-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Rec')} 
#' @slot TAC A numeric value with the TAC recommendation 
#' @slot Effort A numeric value with the effort recommendation as a fraction of current (nyear) fishing effort
#' @slot Spatial A boolean vector of length 'nareas' specifying if area is open (1) or closed (0) to fishing 
#' @slot Allocate A boolean value describing if effort should be re-allocated from close to open areas
#' @slot LR5 smallest length at 5 per cent retention - in absolute units - i.e same units as Linf and L50
#' @slot LFR smallest length at full retention  - in absolute units - i.e same units as Linf and L50
#' @slot HS upper harvest slot (no retention above this)  - in absolute units - i.e same units as Linf and L50
#' @slot Rmaxlen retention of the largest size class - fraction between 0 and 1
#' @slot L5 smallest length at 5 per cent selection - in absolute units - i.e same units as Linf and L50
#' @slot LFS smallest length at full selection  - in absolute units - i.e same units as Linf and L50
#' @slot Vmaxlen selection of the largest size class - fraction between 0 and 1
#' @slot Fdisc fraction of discarded fish that die - fraction between 0 and 1
#' @slot DR Discard rate - the fraction of caught fish that are discarded
#' @slot Misc An empty list that can be used to store information and pass on to MPs in future 
#' @author A. Hordyk
#' @keywords classes
setClass("Rec", representation(
  TAC = "numeric",
  Effort = "numeric", 
  Spatial="numeric", Allocate = "numeric", 
  LR5 = "numeric", LFR = "numeric", HS="numeric", Rmaxlen="numeric", 
  L5 = "numeric", LFS = "numeric", Vmaxlen="numeric", 
  Fdisc = "numeric",
  DR='numeric',
  Misc="list"))


#' Show the output of a single MP recommendation
#'
#' @param object object of class Rec
#' @rdname show-Rec
#' @export
setMethod("show", signature = (object="Rec"), function(object) {

 Rec <- object
 slots <- slotNames(Rec)
 recList <- list()
 perc <- 0.5
 for (X in slots) { # sequence along recommendation slots 
   if (X == "Misc") { # convert to a list nsim by nareas
     rec <- slot(Rec, X)
   } else {
     rec <- slot(Rec,X) # unlist(lapply(temp, slot, name=X))
   }
   if (X == "Spatial") { # convert to a matrix nsim by nareas
     nsims <- 1
     nareas <- max(2,length(rec))
     rec <- matrix(rec, nareas, nsims, byrow=FALSE)   
   }
   recList[[X]] <- rec
 }
 
 names <- c("TAC", "Effort", "LR5", "LFR", "HS", "Rmaxlen",
            "L5", "LFS", 'Vmaxlen', 'Fdisc', 'DR', 'Spatial')
            
 mat <- matrix(0, nrow=1, ncol=length(names)+nareas-1)
 count <- 0 
 for (x in names) {
   temp <- recList[[x]]
   count <- count + 1 
   if (x!="Spatial") {
     mat[,count] <- quantile(temp, probs=perc, na.rm=TRUE)
   } else {
     mat[,count:ncol(mat)] <- t(matrix(unlist(temp), nrow=nareas, ncol=1))
   }
 }
 names[length(names)] <- "Area 1"
 names <- c(names, paste('Area', 2:nareas))
 if (perc !=0.5) names[1] <- paste0(names[1], ' (', perc, ' perc')
 if (perc ==0.5) names[1] <- paste0(names[1], ' (median)')
 colnames(mat) <- names
 
 if (nrow(mat) == 1) {
   mat <- as.data.frame(mat)
   matout <- mat[!is.na(mat)]
   names(matout) <- names[!is.na(mat)]
   print(matout)
 }
})


# -- Hist Object Class ----
#' Class \code{'Hist'}
#' 
#' An object for storing information generated by the end of the historical simulations  
#' 
#' @name Hist-class
#' @docType class
#' 
#' @slot Data The Data object at the end of the historical period
#' 
#' @template Obs_desc
#' @slot OM A numeric data.frame with nsim rows with sampled Stock & Fleet 
#' parameters
#' @slot AtAge A named list with arrays (dim nsim, maxage, nyears+proyears):
#'  \itemize{
#'  \item Length: Length-at-age for each simulation, age, and year
#'  \item Weight: Weight-at-age for each simulation, age, and year
#'  \item Select: Selectivity-at-age for each simulation, age, and year 
#'  \item Retention: Retention-at-age for each simulation, age, and year
#'  \item Maturity: Maturity-at-age for each simulation, age, and year
#'  \item N.Mortality: Natural mortality-at-age for each simulation, age, and year
#'  \item Nage: Total numbers by simulation, age, and year
#'  \item SSBage: Spawning stock biomass by simulation, age, and year
#'  \item FM: Fishing mortality by simulation, age, year, and area
#'  }
#' @slot TSdata A named list with population dynamics by simulation and year :
#'  \itemize{
#'  \item VB: Vulnerable biomass
#'  \item SSB: Spawning stock biomass  
#'  \item B: Total biomass
#'  \item Removals: Removals 
#'  \item Catch: Retained catch (will be same as removals unless there is discard mortality)
#'  \item Rec: Recruitment 
#'  \item N: Total numbers 
#'  \item Find: Historical fishing effort 
#'  \item Marray: Average adult natural mortality (historical & projection)
#'  \item RecDev: Recruitment deviations (historical & projection)
#' } 
#' 
#' @slot Ref A numeric data.frame with nsim rows containing biological 
#' reference points:
#'  \itemize{
#'   \item B0: Average unfished total biomass
#'   \item Blow: Spawning stock biomass where it takes MGThorizon x MGT to 
#'   reach Bfrac of BMSY
#'   \item BMSY: Average total biomass corresponding with MSY
#'   \item BMSY_B0: Ratio of BMSY to B0
#'   \item FMSY: Fishing mortality rate corresponding with MSY 
#'   \item FMSY_M: Ratio of FMSY to (adult) M 
#'   \item MGT: Mean generation time 
#'   \item MSY: Maximum sustainable yield
#'   \item N0: Average unfished numbers
#'   \item R0: Average unfished recruitment
#'   \item RefY: Maximum yield obtained in forward projections with a fixed F
#'   \item SSB0: Average unfished spawning biomass 
#'   \item SSBMSY: Average spawning biomass corresponding with MSY 
#'   \item SSBMSY_SSB: Ratio of SSBMSY to SSB0 
#'   \item UMSY: Exploitation rate corresponding with MSY 
#'   \item VBMSY: Average vulnerable biomass corresponding with MSY
#' }
#' 
#' @slot SampPars All sampled Stock, Fleet, Obs, and Imp parameters
#' @slot Misc A list of additional information
#' @author A. Hordyk
#' @keywords classes
setClass("Hist", representation(
  Data = 'Data',
  Obs = 'data.frame',
  OM = 'data.frame',
  AtAge = 'list',
  TSdata = 'list',
  Ref = "data.frame",
  SampPars='list',
  Misc = 'list'
  ))




# #' Class \code{'lmmodel'}
# #'
# #' An object for storing fitted linear model objects in this case the
# #' relationship between M, age-at-maturity and the von B. K parameter.
# #'
# #'
# #' @name lmmodel-class
# #' @docType class
# #' @section Objects from the Class: Objects can be created by calls of the form
# #' \code{new('lmmodel', stock)}. %% ~~ describe objects here ~~
# #' @author T. Carruthers
# #' @keywords classes
# #' @examples
# #'
# #' newdata<-new('lmmodel','Name',new('list'))
# #'
# setClass("lmmodel",representation(Name="character",models="list"))
# # initialize lmmodel
# setMethod("initialize", "lmmodel", function(.Object,Name,models){
#   .Object@Name<-Name
#   .Object@models<-models
#   .Object
# })
