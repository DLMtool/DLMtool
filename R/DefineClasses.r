#' Label class union for performance metric objects
#' 
#' @description Used internally. Nothing to see here!
#'  
#' @export
setClassUnion(name="label.class", members=c("call", "character", "function"))

#' Prob class union for performance metric objects
#' 
#' @description Used internally. Nothing to see here!
#'  
#' @export
setClassUnion(name="prob.class", members=c("matrix", "numeric", "data.frame"))


# ---- Data Class ----

#' Class \code{'Data'}
#' 
#' An object for storing data for analysis using data-limited methods
#' 
#' 
#' @name Data-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Data', stock)} 
#' @slot Name The name of the Data object. Single value. Character string  
#' @slot Year Years that corresponding to catch and relative abundance data. Vector nyears long. Positive integer
#' @slot Cat Total annual catches. Matrix of nsim rows and nyears columns. Non-negative real numbers 
#' @slot Ind Relative abundance index. Matrix of nsim rows and nyears columns. Non-negative real numbers
#' @slot Rec Recent recruitment strength. Vector of length nsim. Positive real numbers 
#' @slot t The number of years corresponding to AvC and Dt. Single value. Positive integer  
#' @slot AvC Average catch over time t. Vector nsim long. Positive real numbers  
#' @slot Dt Depletion over time t SSB(now)/SSB(now-t+1). Vector nsim long. Fraction  
#' @slot Mort Natural mortality rate. Vector nsim long. Positive real numbers 
#' @slot FMSY_M An assumed ratio of FMSY to M. Vector nsim long. Positive real numbers  
#' @slot BMSY_B0 The most productive stock size relative to unfished. Vector nsim long. Fraction  
#' @slot L50 Length at 50 percent maturity. Vector nsim long. Positive real numbers 
#' @slot L95 Length at 95 percent maturity. Vector nsim long. Positive real numbers 
#' @slot ML Mean length time series. Matrix of nsim rows and nyears columns. Non-negative real numbers
#' @slot Lbar Mean length of catches over Lc. Matrix of nsim rows and nyears columns. Positive real numbers  
#' @slot Lc Modal length of catches. Matrix of nsim rows and nyears columns. Positive real numbers  
#' @slot LFC Length at first capture. Matrix of nsim rows and nyears columns. Positive real numbers 
#' @slot LFS Shortest length at full selection. Matrix of nsim rows and nyears columns. Positive real numbers  
#' @slot CAA Catch at Age data. Array of dimensions nsim x nyears x MaxAge. Non-negative integers
#' @slot Dep Stock depletion SSB(current)/SSB(unfished). Vector nsim long. Fraction.  
#' @slot Abun An estimate of absolute current vulnerable abundance. Vector nsim long. Positive real numbers 
#' @slot SpAbun An estimate of absolute current spawning stock abundance. Vector nsim long. Positive real numbers 
#' @slot vbK The von Bertalanffy growth coefficient K. Vector nsim long. Positive real numbers  
#' @slot vbLinf Maximum length. Vector nsim long. Positive real numbers
#' @slot vbt0 Theoretical age at length zero. Vector nsim long. Non-positive real numbers
#' @slot LenCV Coefficient of variation of length-at-age (assumed constant for all age classes). Vector nsim long. Positive real numbers 
#' @slot wla Weight-Length parameter alpha. Vector nsim long. Positive real numbers 
#' @slot wlb Weight-Length parameter beta. Vector nsim long. Positive real numbers  
#' @slot steep Steepness of stock-recruitment relationship. Vector nsim long. Value in the range of one-fifth to 1 
#' @slot CV_Cat Coefficient of variation in annual catches. Vector nsim long. Positive real numbers 
#' @slot CV_Dt Coefficient of variation in depletion over time t. Vector nsim long. Positive real numbers  
#' @slot CV_AvC Coefficient of variation in average catches over time t. Vector nsim long. Positive real numbers
#' @slot CV_Ind Coefficient of variation in the relative abundance index. Vector nsim long. Positive real numbers 
#' @slot CV_Mort Coefficient of variation in natural mortality rate. Vector nsim long. Positive real numbers 
#' @slot CV_FMSY_M Coefficient of variation in the ratio in FMSY/M. Vector nsim long. Positive real numbers 
#' @slot CV_BMSY_B0 Coefficient of variation in the position of the most productive stock size relative to unfished. Vector nsim long. Positive real numbers 
#' @slot CV_Dep Coefficient of variation in current stock depletion. Vector nsim long. Positive real numbers 
#' @slot CV_Abun Coefficient of variation in estimate of absolute current stock size. Vector nsim long. Positive real numbers 
#' @slot CV_vbK Coefficient of variation in the von Bertalanffy K parameter. Vector nsim long. Positive real numbers 
#' @slot CV_vbLinf Coefficient of variation in maximum length. Vector nsim long. Positive real numbers 
#' @slot CV_vbt0 Coefficient of variation in age at length zero. Vector nsim long. Positive real numbers 
#' @slot CV_L50 Coefficient of variation in length at 50 per cent maturity. Vector nsim long. Positive real numbers 
#' @slot CV_LFC Coefficient of variation in length at first capture. Vector nsim long. Positive real numbers 
#' @slot CV_LFS Coefficient of variation in length at full selection. Vector nsim long. Positive real numbers 
#' @slot CV_wla Coefficient of variation in weight-length parameter a. Vector nsim long. Positive real numbers 
#' @slot CV_wlb Coefficient of variation in weight-length parameter b. Vector nsim long. Positive real numbers 
#' @slot CV_steep Coefficient of variation in steepness. Vector nsim long. Positive real numbers   
#' @slot sigmaL Assumed observaton error of the length composition data. Vector nsim long. Positive real numbers 
#' @slot MaxAge Maximum age. Vector nsim long. Positive integer 
#' @slot CAL_bins The values delimiting the length bins for the catch-at-length data. Vector. Non-negative real numbers
#' @slot CAL Catch-at-length data. An array with dimensions nsim x nyears x length(CAL_bins). Non-negative integers 
#' @slot TAC The calculated catch limits (function TAC). An array with dimensions PosMPs x replicate TAC samples x nsim. Positive real numbers  
#' @slot Sense The results of the sensitivity analysis (function Sense). An array with dimensions PosMPs x sensitivity increments. Positive real numbers  
#' @slot Units Units of the catch/absolute abundance estimates. Single value. Character string
#' @slot Ref A reference management level (eg a catch limit). Single value. Positive real number  
#' @slot Ref_type Type of reference management level (eg 2009 catch limit). Single value. Character string 
#' @slot Log A record of events. Single value. Character string 
#' @slot params A place to store estimated parameters. An object. R list    
#' @slot PosMPs The methods that can be applied to these data. Vector. Character strings  
#' @slot MPs The methods that were applied to these data. Vector. Character strings  
#' @slot OM A table of operating model conditions. R table object of nsim rows. Real numbers  
#' @slot Obs A table of observation model conditions. R table object of nsim rows. Real numbers 
#' @slot Cref Reference or target catch level (eg MSY). Vector of length nsim. Positive real numbers 
#' @slot Iref Reference or target relative abundance index level (eg BMSY / B0). Vector of length nsim. Positive real numbers 
#' @slot Bref Reference or target biomass level (eg BMSY). Vector of length nsim. Positive real numbers 
#' @slot CV_Cref Log-normal CV for reference or target catch level. Vector of length nsim. Positive real numbers  
#' @slot CV_Iref Log-normalCV for reference or target relative abundance index level. Vector of length nsim. Positive real numbers
#' @slot CV_Bref Log-normal CV for reference or target biomass level. Vector of length nsim. Positive real numbers 
#' @slot CV_Rec Log-normal CV for recent recruitment strength. Vector of length nsim. Positive real numbers 
#' @slot MPrec The previous recommendation of a management procedure. Vector of length nsim. Positive real numbers   
#' @slot MPeff The current level of effort. Vector of length nsim. Positive real numbers  
#' @slot LHYear The last historical year of the simulation (before projection). Single value. Positive integer  
#' @slot nareas Number of fishing areas. Vector of length nsim. Non-negative integer 
#' @slot Misc Other information for MPs. An object. R list   
#' 
#' @author T. Carruthers and A. Hordyk
#' @keywords classes
#' @examples
#' 
#' newdata<-new('Data')
#' 
setClass("Data", representation(Name = "character", Year = "vector", 
                                Cat = "matrix", Ind = "matrix", Rec = "matrix", t = "vector",
                                AvC = "vector", Dt = "vector", Mort = "vector", FMSY_M = "vector", 
                                BMSY_B0 = "vector", L50 = "vector", L95 = "vector", 
                                ML = "array", Lbar = "array", Lc = "array",
                                LFC = "vector", LFS = "vector", CAA = "array", Dep = "vector", 
                                Abun = "vector", SpAbun="vector", vbK = "vector", vbLinf = "vector", vbt0 = "vector", 
                                LenCV="vector", wla = "vector", wlb = "vector",  steep = "vector", 
                                CV_Cat = "vector", CV_Dt = "vector", CV_AvC = "vector", 
                                CV_Ind = "vector", CV_Mort = "vector", CV_FMSY_M = "vector",
                                CV_BMSY_B0 = "vector", CV_Dep = "vector", CV_Abun = "vector",
                                CV_vbK = "vector", CV_vbLinf = "vector", CV_vbt0 = "vector",
                                CV_L50 = "vector", CV_LFC = "vector", CV_LFS = "vector", 
                                CV_wla = "vector", CV_wlb = "vector", CV_steep = "vector", 
                                sigmaL = "vector", MaxAge = "vector",
                                CAL_bins = "numeric", CAL = "array", 
                                TAC = "array", Sense = "array",
                                Units = "character", Ref = "numeric", Ref_type = "character", 
                                Log = "list", params = "list", PosMPs = "vector",
                                MPs = "vector", OM = "data.frame", Obs = "data.frame", 
                                Cref = "vector", Iref = "vector",  Bref = "vector", 
                                CV_Cref = "vector", CV_Iref = "vector", CV_Bref = "vector", 
                                CV_Rec = "vector",  
                                MPrec = "vector", MPeff = "vector", LHYear = "numeric", 
                                nareas = "numeric", Misc = "list"))

# initialize Data
setMethod("initialize", "Data", function(.Object, stock = "nada", dec=c(".", ",")) {
  # .Object }) .Object<-new('Data') run an error check here
  if (file.exists(stock)) {
    dec <- match.arg(dec)
    if (dec == ".") dat <- read.csv(stock, header = F, colClasses = "character")  # read 1st sheet
    if (dec == ",") dat <- read.csv2(stock, header = F, colClasses = "character")  # read 1st sheet
    dname <- dat[, 1]
    dat <- dat[, 2:ncol(dat)]
    
    .Object@Name <- dat[match("Name", dname), 1]
    .Object@Year <- as.numeric(dat[match("Year", dname), dat[match("Year", dname), ] != ""])
    .Object@Cat <- matrix(as.numeric(dat[match("Catch", dname), dat[match("Catch", dname), ] != ""]), nrow = 1)
    .Object@Ind <- matrix(as.numeric(dat[match("Abundance index", dname), 1:length(.Object@Year)]), nrow = 1)
    .Object@Rec <- matrix(as.numeric(dat[match("Recruitment", dname), 1:length(.Object@Year)]), nrow = 1)
    .Object@t <- as.numeric(dat[match("Duration t", dname), 1])
    .Object@AvC <- as.numeric(dat[match("Average catch over time t", dname), 1])
    .Object@Dt <- as.numeric(dat[match("Depletion over time t", dname), 1])
    .Object@Mort <- as.numeric(dat[match("M", dname), 1])
    .Object@FMSY_M <- as.numeric(dat[match("FMSY/M", dname), 1])
    .Object@BMSY_B0 <- as.numeric(dat[match("BMSY/B0", dname), 1])
    .Object@L50 <- as.numeric(dat[match("Length at 50% maturity", dname), 1])
    .Object@L95 <- as.numeric(dat[match("Length at 95% maturity", dname), 1])
    .Object@ML <- matrix(as.numeric(dat[match("Mean length", dname), 1:length(.Object@Year)]), nrow = 1)
    .Object@Lbar <- matrix(as.numeric(dat[match("Mean length Lc", dname), 1:length(.Object@Year)]), nrow = 1)
    .Object@Lc <- matrix(as.numeric(dat[match("Modal length", dname), 1:length(.Object@Year)]), nrow = 1)
    .Object@LFC <- as.numeric(dat[match("Length at first capture",  dname), 1])
    .Object@LFS <- as.numeric(dat[match("Length at full selection", dname), 1])
  
    CAAy <- grep("CAA", dname)[1:length(grep("CAA", dname))]
    CAAa <- sum(dat[CAAy[1], ] != "")
    if (!is.na(CAAa)) {
      .Object@CAA <- array(as.numeric(as.matrix(dat[CAAy, 1:CAAa])),  dim = c(1, length(CAAy), CAAa))
    }
    .Object@Dep <- as.numeric(dat[match("Current stock depletion",  dname), 1])
    .Object@Abun <- as.numeric(dat[match("Current stock abundance",  dname), 1]) 
    .Object@SpAbun <- as.numeric(dat[match("Current spawning stock abundance",  dname), 1])
    .Object@vbK <- as.numeric(dat[match("Von Bertalanffy K parameter", dname), 1])
    .Object@vbLinf <- as.numeric(dat[match("Von Bertalanffy Linf parameter", dname), 1])
    .Object@vbt0 <- as.numeric(dat[match("Von Bertalanffy t0 parameter", dname), 1])
    .Object@LenCV <- as.numeric(dat[match("CV of length-at-age", dname), 1])
    .Object@wla <- as.numeric(dat[match("Length-weight parameter a", dname), 1])
    .Object@wlb <- as.numeric(dat[match("Length-weight parameter b", dname), 1])
    .Object@steep <- as.numeric(dat[match("Steepness", dname), 1])
    
    .Object@CV_Cat <- as.numeric(dat[match("CV Catch", dname), 1])
    .Object@CV_Dt <- as.numeric(dat[match("CV Depletion over time t", dname), 1])
    .Object@CV_AvC <- as.numeric(dat[match("CV Average catch over time t", dname), 1])
    .Object@CV_Ind <- as.numeric(dat[match("CV Abundance index", dname), 1])
    .Object@CV_Mort <- as.numeric(dat[match("CV M", dname), 1])
    .Object@CV_FMSY_M <- as.numeric(dat[match("CV FMSY/M", dname),  1])
    .Object@CV_BMSY_B0 <- as.numeric(dat[match("CV BMSY/B0", dname), 1])
    .Object@CV_Dep <- as.numeric(dat[match("CV current stock depletion", dname), 1])
    .Object@CV_Abun <- as.numeric(dat[match("CV current stock abundance", dname), 1])
    .Object@CV_vbK <- as.numeric(dat[match("CV von B. K parameter", dname), 1])
    .Object@CV_vbLinf <- as.numeric(dat[match("CV von B. Linf parameter", dname), 1])
    .Object@CV_vbt0 <- as.numeric(dat[match("CV von B. t0 parameter", dname), 1])
    .Object@CV_L50 <- as.numeric(dat[match("CV Length at 50% maturity", dname), 1])
    .Object@CV_LFC <- as.numeric(dat[match("CV Length at first capture", dname), 1])
    .Object@CV_LFS <- as.numeric(dat[match("CV Length at full selection", dname), 1])
    .Object@CV_wla <- as.numeric(dat[match("CV Length-weight parameter a", dname), 1])
    .Object@CV_wlb <- as.numeric(dat[match("CV Length-weight parameter b", dname), 1])
    .Object@CV_steep <- as.numeric(dat[match("CV Steepness", dname),  1])
    .Object@sigmaL <- as.numeric(dat[match("Sigma length composition", dname), 1])
    
    .Object@MaxAge <- as.numeric(dat[match("Maximum age", dname), 1])
    
    if (length(grep("CAL", dname)) > 1) {
      CAL_bins <- as.numeric(dat[match("CAL_bins", dname), dat[match("CAL_bins", dname), ] != ""])
      nCAL <- length(CAL_bins) - 1
      .Object@CAL_bins <- CAL_bins
      CALdat <- grep("CAL ", dname)
      if (length(CALdat) > 0) .Object@CAL <- array(as.numeric(as.matrix(dat[CALdat, 1:nCAL])),dim = c(1, length(CALdat), nCAL))
    }
    
    .Object@Units <- dat[match("Units", dname), 1]
    .Object@Ref <- as.numeric(dat[match("Reference OFL", dname), 1])
    .Object@Ref_type <- dat[match("Reference OFL type", dname), 1]
    
    .Object@Cref <- as.numeric(dat[match("Cref", dname), 1])
    .Object@Iref <- as.numeric(dat[match("Iref", dname), 1])
    .Object@Bref <- as.numeric(dat[match("Bref", dname), 1])
    
    .Object@CV_Cref <- as.numeric(dat[match("CV Cref", dname), 1])
    .Object@CV_Iref <- as.numeric(dat[match("CV Iref", dname), 1])
    .Object@CV_Bref <- as.numeric(dat[match("CV Bref", dname), 1])
    .Object@CV_Rec <- as.numeric(dat[match("CV Rec", dname), 1])

    .Object@MPrec <- as.numeric(dat[match("MPrec", dname), 1])
    .Object@MPeff <- as.numeric(dat[match("MPeff", dname), 1])
    
    .Object@LHYear <- as.numeric(dat[match("LHYear", dname), 1])
    .Object@nareas <- as.numeric(dat[match("nareas", dname), 1])

    .Object@Log[[1]] <- paste("Created:", Sys.time())
    .Object@params <- new("list")
    .Object@OM <- data.frame(NA)
    .Object@Obs <- data.frame(NA)
    .Object@TAC <- array(NA, dim = c(1, 1, 1))
    .Object@Sense <- array(NA, dim = c(1, 1, 1))
    .Object@PosMPs <- NA
    .Object@MPs <- NA
    
  } else {
    if (stock != "MSE") {
      if (!is.na(stock)) print("Couldn't find specified csv file, blank DLM object created")
    }
  }
  
  if (is.na(.Object@MPeff) || length(.Object@MPeff)==0) .Object@MPeff <- 1 
  
  # Standardise Index if not already 
  .Object@Ind <- .Object@Ind/mean(.Object@Ind, na.rm=TRUE)
  
  # Default value
  if (NAor0(.Object@LenCV)) .Object@LenCV <- 0.1
  if (NAor0(.Object@CV_Cat)) .Object@CV_Cat <- 0.2
  if (NAor0(.Object@CV_Dt)) .Object@CV_Dt <- 0.25
  if (NAor0(.Object@CV_AvC)) .Object@CV_AvC <- 0.2
  if (NAor0(.Object@CV_Ind)) .Object@CV_Ind <- 0.2
  if (NAor0(.Object@CV_Mort)) .Object@CV_Mort <- 0.2
  if (NAor0(.Object@CV_FMSY_M)) .Object@CV_FMSY_M <- 0.2
  if (NAor0(.Object@CV_BMSY_B0)) .Object@CV_BMSY_B0 <- 0.045
  if (NAor0(.Object@CV_Cref)) .Object@CV_Cref <- 0.2
  if (NAor0(.Object@CV_Bref)) .Object@CV_Bref <- 0.2
  if (NAor0(.Object@CV_Iref)) .Object@CV_Iref <- 0.2
  if (NAor0(.Object@CV_Rec)) .Object@CV_Rec <- 0.2
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
  
  if (length(.Object@sigmaL) == 0) .Object@sigmaL <- 0.2
  if (length(.Object@CAA) == 0) .Object@CAA <- array(NA, c(1, 1, 1))
  if (length(.Object@CAL) == 0) .Object@CAL <- array(NA, c(1, 1, 1))
  if (length(.Object@CAL_bins) == 0) .Object@CAL_bins <- 1
  if (length(.Object@TAC) == 0) .Object@TAC <- array(1, c(1, 1))
  # if (length(.Object@TACbias) == 0) .Object@TACbias <- array(1, c(1, 1))
  if (length(.Object@Sense) == 0) .Object@Sense <- array(1, c(1, 1))
  if (length(.Object@ML) == 0)  .Object@ML <- array(NA, c(1, 1))
  if (length(.Object@Lbar) == 0) .Object@Lbar <- array(NA, c(1, 1))
  if (length(.Object@Lc) == 0) .Object@Lc <- array(NA, c(1, 1))
  
  .Object
})



# ---- Fease Class ----
#' Class \code{'Fease'}
#' 
#' An object for storing information about what data are available or might be
#' available
#' 
#' @name Fease-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Fease', stock)}
#'
#' @slot Name The name of the data feasibility object
#' @slot Case The names of the data feasibility cases
#' @slot Catch Total annual catches
#' @slot Index An index of relative abundance, catch per unit effort data or of fishing mortality rate (effort)
#' @slot Natural_mortality_rate From Maximum age, Tagging data, early fishery catch composition data
#' @slot Maturity_at_length From gonadal analysis, growth and natural mortality rate estimates
#' @slot Growth Paired length and age observations, maximum length and an estimate of natural mortality rate
#' @slot Length_weight_conversion Paired weight and length observations, equivalent data from a similar species
#' @slot Fleet_selectivity Length composition of catches with growth curve and natural mortality rate, estimates from a similar fleet type targetting a similar species
#' @slot Catch_at_length Length composition of catches (length samples)
#' @slot Catch_at_age Age composition of catches (age samples)
#' @slot Recruitment_index Spawn survey, estimates from a stock assessment, VPA analysis of catch composition data
#' @slot Stock_recruitment_relationship Stock assessment, a stock assessment of a similar species
#' @slot Target_catch An agreed annual catch target, MSY proxy
#' @slot Target_biomass An agreed absolute biomass target, mean historical biomass estimate
#' @slot Target_index An agreed catch rate target
#' @slot Abundance Fishery independent survey, current fishing mortality rate from recent length composition, natural mortality rate, maturity at age, growth and stock recruitment relationship, habitat and relative density extrapolation
#'
#' @author T. Carruthers and A. Hordyk
#' @keywords classes
#' @examples
#' 
#' newdata<-new('Fease')
#' 
setClass("Fease", representation(Name = "character", Case = "character", 
                                 Catch = "numeric", Index = "numeric", Natural_mortality_rate = "numeric", 
                                 Maturity_at_length = "numeric", Growth = "numeric", Length_weight_conversion = "numeric", 
                                 Fleet_selectivity = "numeric", Catch_at_length = "numeric", Catch_at_age = "numeric", 
                                 Recruitment_index = "numeric", Stock_recruitment_relationship = "numeric", 
                                 Target_catch = "numeric", Target_biomass = "numeric", Target_index = "numeric", 
                                 Abundance = "numeric"))

# initialize Fease
setMethod("initialize", "Fease", function(.Object, file = "nada", ncases = 1, dec=c(".", ",")) {
  # run an error check here
  if (file.exists(file)) {
    Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
    dec <- match.arg(dec)
    if (dec == ".") dat <- read.csv(file, header = F, colClasses = "character", col.names = paste0("V", 
                                                                                   1:Ncol))  # read 1st sheet
    if (dec == ",") dat <- read.csv2(file, header = F, colClasses = "character", col.names = paste0("V", 
                                                                                                   1:Ncol))  # read 1st sheet
    nr <- nrow(dat)
    ncases = ncol(dat) - 1
    dname <- dat[, 1]
    if (ncases == 1) dat <- array(dat[, 2:ncol(dat)], dim = c(nr, ncases))
    if (ncases > 1) dat <- dat[, 2:ncol(dat)]
    .Object@Name <- dat[match("Name", dname), 1]
    .Object@Case <- as.character(dat[match("Case", dname), 1:ncases])
    .Object@Catch <- as.numeric(dat[match("Catch", dname), 1:ncases])
    .Object@Index <- as.numeric(dat[match("Index", dname), 1:ncases])
    .Object@Natural_mortality_rate <- as.numeric(dat[match("Natural_mortality_rate",  dname), 1:ncases])
    .Object@Maturity_at_length <- as.numeric(dat[match("Maturity_at_length",  dname), 1:ncases])
    .Object@Growth <- as.numeric(dat[match("Growth", dname), 1:ncases])
    .Object@Length_weight_conversion <- as.numeric(dat[match("Length_weight_conversion", dname), 1:ncases])
    .Object@Fleet_selectivity <- as.numeric(dat[match("Fleet_selectivity", dname), 1:ncases])
    .Object@Catch_at_length <- as.numeric(dat[match("Catch_at_length", dname), 1:ncases])
    .Object@Catch_at_age <- as.numeric(dat[match("Catch_at_age", dname), 1:ncases])
    .Object@Recruitment_index <- as.numeric(dat[match("Recruitment_index", dname), 1:ncases])
    .Object@Stock_recruitment_relationship <- as.numeric(dat[match("Stock_recruitment_relationship", dname), 1:ncases])
    .Object@Target_catch <- as.numeric(dat[match("Target_catch", dname), 1:ncases])
    .Object@Target_biomass <- as.numeric(dat[match("Target_biomass", dname), 1:ncases])
    .Object@Target_index <- as.numeric(dat[match("Target_index", dname), 1:ncases])
    .Object@Abundance <- as.numeric(dat[match("Abundance", dname), 1:ncases])
  } else {
    .Object@Name <- "Blank DLM_Fease"
    .Object@Case <- "Case 1"
    .Object@Catch <- 1
    .Object@Index <- 1
    .Object@Natural_mortality_rate <- 1
    .Object@Maturity_at_length <- 1
    .Object@Growth <- 1
    .Object@Length_weight_conversion <- 1
    .Object@Fleet_selectivity <- 1
    .Object@Catch_at_length <- 1
    .Object@Catch_at_age <- 1
    .Object@Recruitment_index <- 1
    .Object@Stock_recruitment_relationship <- 1
    .Object@Target_catch <- 1
    .Object@Target_biomass <- 1
    .Object@Target_index <- 1
    .Object@Abundance <- 1
  }
  .Object
  
})

# ---- Stock Class ----
#' Class \code{'Stock'}
#' 
#' An operating model component that specifies the parameters of the population
#' dynamics model
#' 
#' 
#' @name Stock-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Stock')}
#' @slot Name The name of the Stock object. Single value. Character string 
#' @slot Common_Name Common name of the species. Character string
#' @slot Species Scientific name of the species. Genus and species name. Character string
#' @slot maxage The maximum age of individuals that is simulated (there is no 'plus group'). Single value. Positive integer
#' @slot R0 The magnitude of unfished recruitment. Single value. Positive real number
#' @slot M Natural mortality rate. Uniform distribution lower and upper bounds. Positive real number 
#' @slot M2 (Optional) Natural mortality rate at age. Vector of length 'maxage'. Positive real number
#' @slot Mexp Exponent of the Lorenzen function assuming an inverse relationship between M and weight. Uniform distribution lower and upper bounds. Real numbers <= 0.
#' @slot Msd Inter-annual variability in natural mortality rate expressed as a coefficient of variation. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Mgrad Mean temporal trend in natural mortality rate, expressed as a percentage change in M per year. Uniform distribution lower and upper bounds. Real numbers 
#' @slot h Steepness of the stock recruit relationship. Uniform distribution lower and upper bounds. Values from 1/5 to 1 
#' @slot SRrel Type of stock-recruit relationship. Single value, switch (1) Beverton-Holt (2) Ricker. Integer 
#' @slot Perr Process error, the CV of lognormal recruitment deviations. Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot AC Autocorrelation in recruitment deviations rec(t)=AC*rec(t-1)+(1-AC)*sigma(t). Uniform distribution lower and upper bounds. Non-negative real numbers 
# #' @slot recgrad Mean temporal trend in log-normal recruitment deviations, expressed as a percentage change per year. Uniform distribution lower and upper bounds. Real numbers 
#' @slot Period (Optional) Period for cyclical recruitment pattern in years. Uniform distribution lower and upper bounds. Non-negative real numbers  
#' @slot Amplitude (Optional) Amplitude in deviation from long-term average recruitment during recruitment cycle (eg a range from 0 to 1 means recruitment decreases or increases by up to 100\% each cycle). Uniform distribution lower and upper bounds. 0 < Amplitude < 1 
#' @slot Linf Maximum length. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot K von Bertalanffy growth parameter k. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot t0 von Bertalanffy theoretical age at length zero. Uniform distribution lower and upper bounds. Non-positive real numbers
#' @slot LenCV Coefficient of variation of length-at-age (assumed constant for all age classes). Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Ksd Inter-annual variability in growth parameter k. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Kgrad Mean temporal trend in growth parameter k, expressed as a percentage change in k per year. Uniform distribution lower and upper bounds. Real numbers 
#' @slot Linfsd Inter-annual variability in maximum length. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Linfgrad Mean temporal trend in maximum length, expressed as a percentage change in Linf per year. Uniform distribution lower and upper bounds. Real numbers 
#' @slot L50 Length at 50 percent maturity. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot L50_95 Length increment from 50 percent to 95 percent maturity. Uniform distribution lower and upper bounds. Positive real numbers 
# @slot FecB Exponent of the length-fecundity relationship, ie, (relative) fecundity-at-length is proportional to length^FecB (uniform distribution)
#' @slot D Current level of stock depletion SSB(current)/SSB(unfished). Uniform distribution lower and upper bounds. Fraction
#' @slot a Length-weight parameter alpha. Single value. Positive real number 
#' @slot b Length-weight parameter beta. Single value. Positive real number
#' @slot Size_area_1 The size of area 1 relative to area 2. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot Frac_area_1 The fraction of the unfished biomass in stock 1. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot Prob_staying The probability of inviduals in area 1 remaining in area 1 over the course of one year. Uniform distribution lower and upper bounds. Positive fraction.
#' @slot Fdisc Fraction of discarded fish that die. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Source A reference to a website or article from which parameters were taken to define the stock object. Single value. Character string. 

#' @author T. Carruthers and A. Hordyk
#' @keywords classes
#' @examples
#' 
#' showClass('Stock')
#' 
# FecB = "numeric"
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

# ---- Fleet Class -----
#' Class \code{'Fleet'}
#' 
#' The component of the operating model that controls fishing dynamics
#' 
#' 
#' @name Fleet-class
#' @docType class
#' @section Creating Object: 
#' Objects can be created by calls of the form \code{new('Fleet')}
#' 
#' @section MPA slot: 
#' Each row should contain year index (e.g 10 for 10th historical year)
#' followed by fraction of area closed to fishing for each area. i.e. each row represents a change and the number of columns is nareas + 1. 
#' The spatial closures are assumed to remain in place for the future projections unless changed by a MP. 
#' Default (if left blank) is all areas are open to fishing in historical period.
#'
#' @slot Name Name of the Fleet object. Single value. Character string. 
#' @slot nyears The number of years for the historical 'spool-up' simulation. Single value. Positive integer 
#' @slot Spat_targ Distribution of fishing in relation to spatial biomass: fishing distribution is proportional to B^Spat_targ. Uniform distribution lower and upper bounds. Real numbers   
#' @slot EffYears Years representing join-points (vertices) of time-varying effort. Vector. Non-negative real numbers 
#' @slot EffLower Lower bound on relative effort corresponding to EffYears. Vector. Non-negative real numbers
#' @slot EffUpper Upper bound on relative effort corresponding to EffYears. Vector. Non-negative real numbers 
#' @slot Esd Additional inter-annual variability in fishing mortality rate. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot qinc Average percentage change in fishing efficiency (applicable only to forward projection and input controls). Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot qcv Inter-annual variability in fishing efficiency (applicable only to forward projection and input controls). Uniform distribution lower and upper bounds. Non-negative real numbers

#' @slot L5 Shortest length corresponding to 5 percent vulnerability. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LFS Shortest length that is fully vulnerable to fishing. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot Vmaxlen The vulnerability of fish at Stock@Linf. Uniform distribution lower and upper bounds. Fraction 
#' @slot isRel Selectivity parameters in units of size-of-maturity (or absolute eg cm). Single value. Boolean.
#' @slot LR5 Shortest length corresponding ot 5 percent retention. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot LFR Shortest length that is fully retained. Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot Rmaxlen The retention of fish at Stock@Linf. Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot DR Discard rate - the fraction of caught fish that are discarded. Uniform distribution lower and upper bounds. Fraction
#' 
#' @slot SelYears (Optional) Years representing join-points (vertices) at which historical selectivity pattern changes. Vector. Positive real numbers
#' @slot AbsSelYears (Optional) Calendar years corresponding with SelYears (eg 1951, rather than 1), used for plotting only. Vector (of same length as SelYears). Positive real numbers  
#' @slot L5Lower (Optional) Lower bound of L5 (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot L5Upper (Optional) Upper bound of L5 (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot LFSLower (Optional) Lower bound of LFS (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot LFSUpper (Optional) Upper bound of LFS (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot VmaxLower (Optional) Lower bound of Vmaxlen (use \code{ChooseSelect} function to set these). Vector. Fraction 
#' @slot VmaxUpper (Optional) Upper bound of Vmaxlen (use \code{ChooseSelect} function to set these). Vector. Fraction
#' @slot CurrentYr The current calendar year (final year) of the historical simulations (eg 2011). Single value. Positive integer. 
#' 
#' @slot MPA (Optional) Matrix specifying spatial closures for historical years. 
#' 
#' @author T. Carruthers and A. Hordyk
#' @keywords classes
#' @examples
#' 
#' showClass('Fleet')
#' 
setClass("Fleet", slots = c(Name = "character", nyears = "numeric", Spat_targ = "numeric", 
                            EffYears = "numeric", EffLower = "numeric", EffUpper = "numeric", Esd = "numeric", 
                            qinc = "numeric", qcv = "numeric",   
                            L5 = "numeric", LFS = "numeric", Vmaxlen = "numeric", isRel = "character",
                            LR5 = "numeric", LFR = "numeric", Rmaxlen = "numeric", DR = "numeric",
                            SelYears = "numeric", AbsSelYears = "numeric",
                            L5Lower = "numeric", L5Upper = "numeric", LFSLower = "numeric", LFSUpper = "numeric", VmaxLower = "numeric", 
                            VmaxUpper = "numeric", CurrentYr="numeric", MPA='matrix'))

# initialize Fleet
setMethod("initialize", "Fleet", function(.Object, file = NA, dec=c(".", ",")) {
  if (!is.na(file)) {
    if (file.exists(file)) {
      dec <- match.arg(dec)
      Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
      if (dec == ".") dat <- read.csv(file, header = F, colClasses = "character", col.names = paste0("V", 1:Ncol))  # read 1st sheet
      if (dec == ",") dat <- read.csv2(file, header = F, colClasses = "character", col.names = paste0("V", 1:Ncol))  # read 1st sheet
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
      
      isMPA <- grep('MPA', dname)
      if (length(isMPA)<1) isMPA <- NA
      if (!is.na(isMPA)) {
        MPA <- temp <- data.matrix(dat[isMPA:nrow(dat),])
        valCols <- !is.na(colSums(MPA))
        MPA <- MPA[,valCols, drop=FALSE]
        valRows <- !is.na(rowSums(MPA))
        MPA <- MPA[valRows, drop=FALSE]
        MPA <- matrix(MPA, nrow=nrow(temp))
        .Object@MPA <- MPA
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
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Obs')} 
#' @slot Name The name of the observation model object. Single value. Character string. 
#' @slot Cobs Log-normal catch observation error expressed as a coefficient of variation. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Cbiascv Log-normal coefficient of variation controlling the sampling of bias in catch observations for each simulation. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot CAA_nsamp Number of catch-at-age observation per time step. Uniform distribution lower and upper bounds. Positive real numbers   
#' @slot CAA_ESS Effective sample size (independent age draws) of the multinomial catch-at-age observation error model. Uniform distribution lower and upper bounds. Positive integers
#' @slot CAL_nsamp Number of catch-at-length observation per time step. Uniform distribution lower and upper bounds. Positive integers
#' @slot CAL_ESS Effective sample size (independent length draws) of the multinomial catch-at-length observation error model. Uniform distribution lower and upper bounds. Positive integers
# #' @slot CALcv Log-normal, CV of length-at-age. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Iobs Observation error in the relative abundance indices expressed as a coefficient of variation. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Ibiascv Log-normal coefficient of variation controlling error in observations of relative abundance index. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Btobs Log-normal coefficient of variation controlling error in observations of current stock biomass among years. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Btbiascv Uniform-log bounds for sampling persistent bias in current stock biomass. Uniform-log distribution lower and upper bounds. Positive real numbers 
#' @slot beta A parameter controlling hyperstability/hyperdepletion where values below 1 lead to hyperstability (an index that decreases slower than true abundance) and values above 1 lead to hyperdepletion (an index that decreases more rapidly than true abundance). Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LenMbiascv Log-normal coefficient of variation for sampling persistent bias in length at 50 percent maturity. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Mbiascv Log-normal coefficient of variation for sampling persistent bias in observed natural mortality rate. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Kbiascv Log-normal coefficient of variation for sampling persistent bias in observed growth parameter K. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot t0biascv Log-normal coefficient of variation for sampling persistent bias in observed t0. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Linfbiascv Log-normal coefficient of variation for sampling persistent bias in observed maximum length. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LFCbiascv Log-normal coefficient of variation for sampling persistent bias in observed length at first capture. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LFSbiascv Log-normal coefficient of variation for sampling persistent bias in length-at-full selection. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot FMSYbiascv Log-normal coefficient of variation for sampling persistent bias in FMSY. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot FMSY_Mbiascv Log-normal coefficient of variation for sampling persistent bias in FMSY/M. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot BMSY_B0biascv Log-normal coefficient of variation for sampling persistent bias in BMSY relative to unfished. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Irefbiascv Log-normal coefficient of variation for sampling persistent bias in relative abundance index at BMSY. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Brefbiascv Log-normal coefficient of variation for sampling persistent bias in BMSY. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Crefbiascv Log-normal coefficient of variation for sampling persistent bias in MSY. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Dbiascv Log-normal coefficient of variation for sampling persistent bias in stock depletion. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Dobs Log-normal coefficient of variation controlling error in observations of stock depletion among years. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot hbiascv Log-normal coefficient of variation for sampling persistent bias in steepness. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Recbiascv Log-normal coefficient of variation for sampling persistent bias in recent recruitment strength. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot B0cv Log-normal coefficient of variation for sampling persistent bias in unfished biomass. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot rcv Log-normal coefficient of variation for sampling persistent bias in intrinsic rate of increase. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot Fcurbiascv Log-normal coefficient of variation for sampling persistent bias in current fishing mortality rate. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot Fcurcv Log-normal coefficient of variation controlling error in observations of current fishing mortality rate among years. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot maxagecv Log-normal coefficient of variation for sampling persistent bias in observation of maximum age. Uniform distribution lower and upper bounds. Positive real numbers  
#'     
#' @author T. Carruthers and A. Hordyk
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
      .Object@CAL_nsamp <- as.numeric(dat[match("CAA_nsamp", dname), 1:2])
      .Object@CAL_ESS <- as.numeric(dat[match("CAA_ESS", dname), 1:2])
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
#' 
#' @name Imp-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Imp')}
#' @slot Name The name of the Implementation error object. Single value. Character string.  
#' @slot TACFrac Mean fraction of TAC taken. Uniform distribution lower and upper bounds. Positive real number. 
#' @slot TACSD Log-normal coefficient of variation in the fraction of Total Allowable Catch (TAC) taken. Uniform distribution lower and upper bounds. Non-negative real numbers. 
#' @slot TAEFrac Mean fraction of TAE taken. Uniform distribution lower and upper bounds. Positive real number. 
#' @slot TAESD Log-normal coefficient of variation in the fraction of Total Allowable Effort (TAE) taken. Uniform distribution lower and upper bounds. Non-negative real numbers.
#' @slot SizeLimFrac The real minimum size that is retained expressed as a fraction of the size. Uniform distribution lower and upper bounds. Positive real number.
#' @slot SizeLimSD Log-normal coefficient of variation controlling mismatch between a minimum size limit and the real minimum size retained. Uniform distribution lower and upper bounds. Non-negative real numbers.
#' @author T. Carruthers and A. Hordyk
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

# Stock slots
#' @slot Common_Name Common name of the species. Character string
#' @slot Species Scientific name of the species. Genus and species name. Character string
#' @slot maxage The maximum age of individuals that is simulated (there is no 'plus group'). Single value. Positive integer
#' @slot R0 The magnitude of unfished recruitment. Single value. Positive real number
#' @slot M Natural mortality rate. Uniform distribution lower and upper bounds. Positive real number 
#' @slot M2 (Optional) Natural mortality rate at age. Vector of length 'maxage'. Positive real number
#' @slot Mexp Exponent of the Lorenzen function assuming an inverse relationship between M and weight. Uniform distribution lower and upper bounds. Real numbers <= 0.
#' @slot Msd Inter-annual variability in natural mortality rate expressed as a coefficient of variation. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Mgrad Mean temporal trend in natural mortality rate, expressed as a percentage change in M per year. Uniform distribution lower and upper bounds. Real numbers 
#' @slot h Steepness of the stock recruit relationship. Uniform distribution lower and upper bounds. Values from 1/5 to 1 
#' @slot SRrel Type of stock-recruit relationship. Single value, switch (1) Beverton-Holt (2) Ricker. Integer 
#' @slot Perr Process error, the CV of lognormal recruitment deviations. Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot AC Autocorrelation in recruitment deviations rec(t)=AC*rec(t-1)+(1-AC)*sigma(t). Uniform distribution lower and upper bounds. Non-negative real numbers 
# #' @slot recgrad Mean temporal trend in log-normal recruitment deviations, expressed as a percentage change per year. Uniform distribution lower and upper bounds. Real numbers 
#' @slot Period (Optional) Period for cyclical recruitment pattern in years. Uniform distribution lower and upper bounds. Non-negative real numbers  
#' @slot Amplitude (Optional) Amplitude in deviation from long-term average recruitment during recruitment cycle (eg a range from 0 to 1 means recruitment decreases or increases by up to 100\% each cycle). Uniform distribution lower and upper bounds. 0 < Amplitude < 1 
#' @slot Linf Maximum length. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot K von Bertalanffy growth parameter k. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot t0 von Bertalanffy theoretical age at length zero. Uniform distribution lower and upper bounds. Non-positive real numbers
#' @slot LenCV Coefficient of variation of length-at-age (assumed constant for all age classes). Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Ksd Inter-annual variability in growth parameter k. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Kgrad Mean temporal trend in growth parameter k, expressed as a percentage change in k per year. Uniform distribution lower and upper bounds. Real numbers 
#' @slot Linfsd Inter-annual variability in maximum length. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Linfgrad Mean temporal trend in maximum length, expressed as a percentage change in Linf per year. Uniform distribution lower and upper bounds. Real numbers 
#' @slot L50 Length at 50 percent maturity. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot L50_95 Length increment from 50 percent to 95 percent maturity. Uniform distribution lower and upper bounds. Positive real numbers 
# @slot FecB Exponent of the length-fecundity relationship, ie, (relative) fecundity-at-length is proportional to length^FecB (uniform distribution)
#' @slot D Current level of stock depletion SSB(current)/SSB(unfished). Uniform distribution lower and upper bounds. Fraction
#' @slot a Length-weight parameter alpha. Single value. Positive real number 
#' @slot b Length-weight parameter beta. Single value. Positive real number
#' @slot Size_area_1 The size of area 1 relative to area 2. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot Frac_area_1 The fraction of the unfished biomass in stock 1. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot Prob_staying The probability of inviduals in area 1 remaining in area 1 over the course of one year. Uniform distribution lower and upper bounds. Positive fraction.
#' @slot Fdisc Fraction of discarded fish that die. Uniform distribution lower and upper bounds. Non-negative real numbers 

# Fleet slots
#' @slot nyears The number of years for the historical 'spool-up' simulation. Single value. Positive integer 
#' @slot Spat_targ Distribution of fishing in relation to spatial biomass: fishing distribution is proportional to B^Spat_targ. Uniform distribution lower and upper bounds. Real numbers   
#' @slot EffYears Years representing join-points (vertices) of time-varying effort. Vector. Non-negative real numbers 
#' @slot EffLower Lower bound on relative effort corresponding to EffYears. Vector. Non-negative real numbers
#' @slot EffUpper Upper bound on relative effort corresponding to EffYears. Vector. Non-negative real numbers 
#' @slot Esd Additional inter-annual variability in fishing mortality rate. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot qinc Average percentage change in fishing efficiency (applicable only to forward projection and input controls). Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot qcv Inter-annual variability in fishing efficiency (applicable only to forward projection and input controls). Uniform distribution lower and upper bounds. Non-negative real numbers

#' @slot L5 Shortest length corresponding to 5 percent vulnerability. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LFS Shortest length that is fully vulnerable to fishing. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot Vmaxlen The vulnerability of fish at Stock@Linf. Uniform distribution lower and upper bounds. Fraction 
#' @slot isRel Selectivity parameters in units of size-of-maturity (or absolute eg cm). Single value. Boolean.
#' @slot LR5 Shortest length corresponding ot 5 percent retention. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot LFR Shortest length that is fully retained. Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot Rmaxlen The retention of fish at Stock@Linf. Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot DR Discard rate - the fraction of caught fish that are discarded. Uniform distribution lower and upper bounds. Fraction
#' 
#' @slot SelYears (Optional) Years representing join-points (vertices) at which historical selectivity pattern changes. Vector. Positive real numbers
#' @slot AbsSelYears (Optional) Calendar years corresponding with SelYears (eg 1951, rather than 1), used for plotting only. Vector (of same length as SelYears). Positive real numbers  
#' @slot L5Lower (Optional) Lower bound of L5 (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot L5Upper (Optional) Upper bound of L5 (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot LFSLower (Optional) Lower bound of LFS (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot LFSUpper (Optional) Upper bound of LFS (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot VmaxLower (Optional) Lower bound of Vmaxlen (use \code{ChooseSelect} function to set these). Vector. Fraction 
#' @slot VmaxUpper (Optional) Upper bound of Vmaxlen (use \code{ChooseSelect} function to set these). Vector. Fraction
#' @slot CurrentYr The current calendar year (final year) of the historical simulations (eg 2011). Single value. Positive integer. .
#' @slot MPA (Optional) Matrix specifying spatial closures for historical years. Each row should contain year index (e.g 10 for 10th historical year)
#' followed by fraction of area closed to fishing for each area. i.e. each row represents a change and the number of columns is nareas + 1. 
#' The spatial closures are assumed to remain in place for the future projections unless changed by a MP. 
#' Default (if left blank) is all areas are open to fishing in historical period.


# Obs slots
#' @slot Cobs Log-normal catch observation error expressed as a coefficient of variation. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Cbiascv Log-normal coefficient of variation controlling the sampling of bias in catch observations for each simulation. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot CAA_nsamp Number of catch-at-age observation per time step. Uniform distribution lower and upper bounds. Positive real numbers   
#' @slot CAA_ESS Effective sample size (independent age draws) of the multinomial catch-at-age observation error model. Uniform distribution lower and upper bounds. Positive integers
#' @slot CAL_nsamp Number of catch-at-length observation per time step. Uniform distribution lower and upper bounds. Positive integers
#' @slot CAL_ESS Effective sample size (independent length draws) of the multinomial catch-at-length observation error model. Uniform distribution lower and upper bounds. Positive integers
# #' @slot CALcv Log-normal, CV of length-at-age. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Iobs Observation error in the relative abundance indices expressed as a coefficient of variation. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Ibiascv Log-normal coefficient of variation controlling error in observations of relative abundance index. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Btobs Log-normal coefficient of variation controlling error in observations of current stock biomass among years. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Btbiascv Uniform-log bounds for sampling persistent bias in current stock biomass. Uniform-log distribution lower and upper bounds. Positive real numbers 
#' @slot beta A parameter controlling hyperstability/hyperdepletion where values below 1 lead to hyperstability (an index that decreases slower than true abundance) and values above 1 lead to hyperdepletion (an index that decreases more rapidly than true abundance). Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LenMbiascv Log-normal coefficient of variation for sampling persistent bias in length at 50 percent maturity. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Mbiascv Log-normal coefficient of variation for sampling persistent bias in observed natural mortality rate. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Kbiascv Log-normal coefficient of variation for sampling persistent bias in observed growth parameter K. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot t0biascv Log-normal coefficient of variation for sampling persistent bias in observed t0. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Linfbiascv Log-normal coefficient of variation for sampling persistent bias in observed maximum length. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LFCbiascv Log-normal coefficient of variation for sampling persistent bias in observed length at first capture. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LFSbiascv Log-normal coefficient of variation for sampling persistent bias in length-at-full selection. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot FMSYbiascv Log-normal coefficient of variation for sampling persistent bias in FMSY. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot FMSY_Mbiascv Log-normal coefficient of variation for sampling persistent bias in FMSY/M. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot BMSY_B0biascv Log-normal coefficient of variation for sampling persistent bias in BMSY relative to unfished. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Irefbiascv Log-normal coefficient of variation for sampling persistent bias in relative abundance index at BMSY. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Brefbiascv Log-normal coefficient of variation for sampling persistent bias in BMSY. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Crefbiascv Log-normal coefficient of variation for sampling persistent bias in MSY. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Dbiascv Log-normal coefficient of variation for sampling persistent bias in stock depletion. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Dobs Log-normal coefficient of variation controlling error in observations of stock depletion among years. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot hbiascv Log-normal coefficient of variation for sampling persistent bias in steepness. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Recbiascv Log-normal coefficient of variation for sampling persistent bias in recent recruitment strength. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot B0cv Log-normal coefficient of variation for sampling persistent bias in unfished biomass. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot rcv Log-normal coefficient of variation for sampling persistent bias in intrinsic rate of increase. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot Fcurbiascv Log-normal coefficient of variation for sampling persistent bias in current fishing mortality rate. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot Fcurcv Log-normal coefficient of variation controlling error in observations of current fishing mortality rate among years. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot maxagecv Log-normal coefficient of variation for sampling persistent bias in observation of maximum age. Uniform distribution lower and upper bounds. Positive real numbers  

# Imp slots
#' @slot TACFrac Mean fraction of TAC taken. Uniform distribution lower and upper bounds. Positive real number. 
#' @slot TACSD Log-normal coefficient of variation in the fraction of Total Allowable Catch (TAC) taken. Uniform distribution lower and upper bounds. Non-negative real numbers. 
#' @slot TAEFrac Mean fraction of TAE taken. Uniform distribution lower and upper bounds. Positive real number. 
#' @slot TAESD Log-normal coefficient of variation in the fraction of Total Allowable Effort (TAE) taken. Uniform distribution lower and upper bounds. Non-negative real numbers.
#' @slot SizeLimFrac The real minimum size that is retained expressed as a fraction of the size. Uniform distribution lower and upper bounds. Positive real number.
#' @slot SizeLimSD Log-normal coefficient of variation controlling mismatch between a minimum size limit and the real minimum size retained. Uniform distribution lower and upper bounds. Non-negative real numbers.

#' @author T. Carruthers and A. Hordyk
#' @keywords classes
#' 
setClass("OM", representation(Name = "character", Agency="character",
                              Region="character", Sponsor="character",
                              Latitude="numeric", Longitude="numeric",
                              nsim="numeric", proyears="numeric", 
                              interval='numeric', pstar='numeric', maxF='numeric', reps='numeric',
                              cpars="list",seed="numeric", Source="character"), contains=c("Stock", "Fleet", "Obs", "Imp"))


# initialize OM
setMethod("initialize", "OM", function(.Object, Stock=NULL, Fleet=DLMtool::Generic_Fleet, 
                                       Obs=DLMtool::Generic_Obs, Imp=DLMtool::Perfect_Imp, 
                                       interval=4, pstar=0.5, maxF=0.8, reps=1, nsim=48, proyears=50) {
  if (is.null(Stock)) {
    message("No Stock object found. Returning a blank OM object") 
    return(makePerf(.Object))
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
#' @slot OM A table of sampled parameter of the operating model. Table object of nsim rows. Real numbers\cr
#'   \itemize{
#'   \item RefY: reference yield, the highest long-term yield (mean over last five years of projection) obtained from a fixed F strategy. This is a useful reference point for framing performance of MPs because it standardizes for starting point and future productivity. 
#'   \item M: instantaneous natural mortality rate
#'   \item Depletion: stock depletion (biomass / unfished biomass) in the final historical year (prior to projection)
#'   \item A: abundance (biomass) updated in each management update of projection
#'   \item BMSY_B0: most productive stock size relative to unfished
#'   \item FMSY_M: fishing mortality rate divided by natural mortality rate
#'   \item Mgrad: mean average percentage gradient in natural mortality rate (percentage per time step)
#'   \item Msd: interannual variability in natural mortality rate (lognormal CV)
#'   \item procsd: process error - CV in log-normal recruitment deviations
#'   \item Esd: interannual variability in historical effort (fishing mortality rate)
#'   \item dFfinal: gradient in fishing mortality rate over final five years of the historical simulation
#'   \item MSY: Maximum Sustainable Yield
#'   \item qinc: mean percentage increase in fishing efficiency (catchability) in projected years (input controls only)
#'   \item qcv: interannual variability in future fishing efficiency (catchability) in projected years (input controls only)
# #'   \item CALcv: variability in lengths at age around the growth curve (normal CV)
#'   \item FMSY: Fishing mortality rate at Maximum Sustainable Yield
#'   \item Linf: maximum length (von Bertalanffy Linf parameter)
#'   \item K: maximum growth rate (von Bertalanffy K parameter)
#'   \item t0: theoretical length at age zero (von Bertalanffy t0 parameter)
#'   \item hs: steepness of the stock recruitment relationship (the fraction of unfished recruitment at a fifth of unfished stock levels)
#'   \item Linfgrad: mean gradient in maximum length (per cent per time step)
#'   \item Kgrad: mean gradient in maximum growth rate (per cent per time step)
#'   \item Linfsd: interannual variability in maximum length (log normal CV)
# #'   \item recgrad: gradient in recruitment strength (age 1 population numbers) over last 10 years of historical simulations
#'   \item Ksd: interannual variability in maximum growth rate (log normal CV)
#'   \item ageM: age at 50 per cent maturity
#'   \item LFS: length at full selection (the shortest length class where fishery selectivity is 100 per cent)
#'   \item age05: the age at 5 percent selectivity (ascending limb of selectivity curve)
#'   \item Vmaxage: the selectivity of the oldest age class (controls dome shape of selectivity curve)
#'   \item LFC: length at first capture, the smallest length that can be caught by the gear
#'   \item OFLreal: the true simulated Over Fishing Limit (FMSY x biomass) updated in each management update of the projection
#'   \item Spat_targ: spatial targetting parameter, fishing mortality rate across areas is proportional to vulnerable biomass raised to the power of this number. 
#'   \item Size_area_1: The size of area 1 relative to area 2
#'   \item Frac_area_1: the fraction of unfished biomass inhabiting area 1 (can be seen as fraction of habitat in area 1 or relative size of area 1)
#'   \item Prob_staying: the probability that individuals in area 1 remain there between time-steps
#'   \item AC: autocorrelation in recruitment
#'  }
#' @slot Obs A table of sampled parameters of the observation model. Table of nsim rows. Real numbers\cr
#'   \itemize{
#'   \item Cbias: bias in observed catches
#'   \item Csd: observation error in observed catches (lognormal CV)
#'   \item CAA_nsamp: the number of catch-at-age observations per time step
#'   \item CAA_ESS: the effective sample size of multinomial catch-at-age observation model (number of independent draws)
#'   \item CAL_nsamp: the number of catch-at-length observations per time step
#'   \item CAL_ESS: the effective sample size of multinomial catch-at-length observation model (number of independent draws)
#'   \item Isd: observation error in relative abundance index (lognormal CV)
#'   \item Dbias: bias in observed stock depletion (also applies to depletion Dt for DCAC)
#'   \item Mbias: bias in observed natural mortality rate
#'   \item FMSY_Mbias: bias in ratio of FMSY to natural mortality rate
#'   \item BMSY_B0bias: bias in ratio of most productive stock size relative to unfished
#'   \item AMbias: bias in age at 50 per cent maturity
#'   \item LFCbias: bias in length at first capture
#'   \item LFSbias: bias in length at full selection
#'   \item Abias: bias in observed current absolute stock biomass
#'   \item Kbias: bias in maximum growth rate (von Bertalanffy K parameter)
#'   \item t0bias: bias in theoretical length at age zero (von Bertalanffy t0 parameter)
#'   \item Linfbias: bias in maximum length (von Bertalanffy Linf parameter)
#'   \item hbias: bias in observed steepness of the stock recruitment relationship
#'   \item Irefbias: bias in abundance index corresponding to BMSY stock levels
#'   \item Crefbias: bias in MSY prediction (target or reference catch)
#'   \item Brefbias: bias in BMSY stock levels (target or reference biomass levels)}
#' @slot B_BMSY Simulated biomass relative to BMSY over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers 
#' @slot F_FMSY Simulated fishing mortality rate relative to FMSY over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers
#' @slot B Simulated stock biomass over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers 
#' @slot SSB Simulated spawning stock biomass over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers
#' @slot VB Simulated vulnerable biomass over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers
#' @slot FM Simulated fishing mortality rate over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers
#' @slot C Simulated catches (taken) over the projection. An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers
#' @slot TAC Simulated Total Allowable Catch (prescribed) over the projection (this is NA for input controls). An array with dimensions: nsim, nMPs, proyears. Non-negative real numbers 
#' @slot SSB_hist Simulated historical spawning stock biomass. An array with dimensions: nsim, nages, nMPs, proyears. Non-negative real numbers
#' @slot CB_hist Simulated historical catches in weight. An array with dimensions: nsim, nages, nMPs, proyears. Non-negative real numbers
#' @slot FM_hist Simulated historical fishing mortality rate. An array with dimensions: nsim, nages, nMPs, proyears. Non-negative real numbers
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


#' Example data object
#' 
#' Example data object with a number of output control MPs run on it, and
#' includes resulting distributions of TACs
#' 
#' 
#' @name ourReefFish
#' @docType data
#' @usage data('ourReefFish')
#' @keywords datasets
#' @examples
#' \dontrun{ 
#' data(ourReefFish)
#' str(ourReefFish)  
#' plot(ourReefFish) 
#' }
#'
NULL


# ---- PMobj Class ----

#' An object for storing data for analysis using data-limited methods
#' 
#' Used interally
#' 
#' @name PMobj-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('PMobj')} 
#' @slot name Name of the Performance Metric. Character 
#' @slot caption A caption to be used in plots. Character, call, or function.
#' @slot Stat Statistic of interest for the PM. Dimensions: nsim, nMP, yrs. Array 
#' @slot Prob Probability (mean over years) Dimensions: nsim by MP.  Matrix, numeric or data.frame  
#' @slot Mean Mean probability (mean over years and simulations). Numeric. Length nMPs 
#' @slot MPs Name of MPs. Single value. Character string  
#' @author  A. Hordyk
#' @importFrom methods show
#' @keywords classes

setClass("PMobj", representation(name = "character",  caption='label.class', 
                                 Stat='array', Prob='prob.class', Mean='numeric',
                                 MPs="character"))


#' Calculate Probabilty
#' 
#' @param PM A PM method 
#' @param MSEobj An object of class MSE
#'
#' @export
#'
calcProb <- function(PM,  MSEobj) {
  mar <- ifelse(MSEobj@nMPs>1, 2, 1) # set margins for apply
  mar <- 1:mar
  apply(PM, mar, mean)
}

#' Calculate Mean Probabilty
#' 
#' @param Prob Prob slot from an object of class PMobj 
#' @param MSEobj An object of class MSE
#'
#' @export
#'
calcMean <- function(Prob, MSEobj) {
  if (class(Prob) == 'matrix') return( apply(Prob , 2, mean))
  if (class(Prob) == 'numeric') return(mean(Prob))
}

show <- function(object) methods::show(object)


#' Show the output of a PM
#'
#' @param object object of class MSE
#' @rdname show-MSE
#' @export
setMethod("show", signature = (object="PMobj"), function(object) {
  cat(object@name)
  cat("\n", object@caption)
  cat("\n")
  
  nMP <- length(object@MPs)
  if (nMP > 1) nsim <- dim(object@Prob)[1]
  if (nMP == 1) nsim <- length(object@Prob)
  
  nprint <- min(nsim, 10)
  if (nMP > 1) df <- data.frame(object@Prob[1:nprint,])
  if (nMP == 1) df <- data.frame(object@Prob[1:nprint])
  if (nMP > 1) lst <- object@Prob[nprint+1,]
  if (nMP == 1) lst <- object@Prob[nprint+1]
  df <- signif(df,2)
  lst <- signif(lst,2)
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
  print(signif(object@Mean,2))
})


#' Summary of MSE object
#'
#' @param object object of class MSE
#' @param ... a list of names of PM methods
#' @param silent Should summary be printed to console? Logical.
#' @rdname summary-MSE
#' @export
setMethod('summary', signature="MSE", function(object, ..., silent=FALSE) {
  PMlist <- unlist(list(...))
  
  if(length(PMlist) == 0) PMlist <- avail("PM")
  if (class(PMlist) != 'character') stop("Must provide names of PM methods")
  # check
  for (X in seq_along(PMlist)) 
    if (!PMlist[X] %in% avail("PM")) stop(PMlist[X], " is not a valid PM method")
  
  if (!silent) message("Calculating Performance Metrics")
  storeMean <- vector('list', length(PMlist))
  storeName <- vector('list', length(PMlist))
  storeHeading <- vector('list', length(PMlist))
  storeMP <- vector('list', length(PMlist))
  for (X in 1:length(PMlist)) {
    runPM <- eval(call(PMlist[[X]],object))
    storeMean[[X]] <- runPM@Mean
    storeName[[X]] <- runPM@name
    # storeHeading[[X]] <- runPM@call
    storeMP[[X]] <- runPM@MPs
  }
  
  df <- data.frame('MP'=storeMP[[1]], signif(do.call('cbind', storeMean),2))
  # heading <- do.call('rbind', storeHeading)
  colnames(df)[2:(length(PMlist)+1)] <- PMlist #caps # gsub(" ", "", caps)
  if (!silent) {
    print(data.frame('Performance Metrics' = do.call('rbind', storeName)))
    cat("\n")
    cat("\nProbability:\n")
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
#' @param object object of class Data
#' @export
setMethod("summary",
          signature(object = "Data"),
          function(object){
            
            old_par <- par(no.readonly = TRUE)
            on.exit(par(list = old_par), add = TRUE)
            
            scols<-c('red','green','blue','orange','brown','purple','dark grey','violet','dark red','pink','dark blue','grey')
            
            #dev.new2(width=8,height=4.5)
            par(mai=c(0.35,0.9,0.2,0.01),c(0.3,0,0,0))
            layout(matrix(c(1,2,1,2,1,2,3,3,3,3),nrow=2))
            plot(object@Year,object@Cat[1,],col="blue",type="l",xlab="Year",ylab=paste("Catch (",object@Units,")",sep=""),ylim=c(0,max(object@Cat[1,],na.rm=T)))
            plot(object@Year,object@Ind[1,],col="orange",type="l",xlab="Year",ylab="Relative abundance",ylim=c(0,max(object@Ind[1,],na.rm=T)))
            
            slots<-c("Dep","Mort","FMSY_M","Dt","BMSY_B0","vbK")
            namey<-c("Stock depletion", "Natural Mortality rate","Ratio of FMSY to M","Depletion over time t","BMSY relative to unfished","Von B. k parameter")
            slotsCV<-c("CV_Dep","CV_Mort","CV_FMSY_M","CV_Dt","CV_BMSY_B0","CV_vbK")
            
            ind<-rep(TRUE,length(slotsCV))
            for(i in 1:length(slotsCV))if(NAor0(attr(object,slots[i]))|NAor0(attr(object,slotsCV[i])))ind[i]<-FALSE
            slots<-slots[ind]
            slotsCV<-slotsCV[ind]
            nrep<-150
            xstore<-array(NA,c(length(slots),nrep))
            ystore<-array(NA,c(length(slots),nrep))
            
            
            for(i in 1:length(slots)){
              mu<-attr(object,slots[i])
              cv<-attr(object,slotsCV[i])
              xstore[i,]<-qlnorm(seq(0,1,length.out=nrep),mconv(mu,cv),sdconv(mu,cv))
              ystore[i,]<-dlnorm(xstore[i,],mconv(mu,cv),sdconv(mu,cv))
            }
            
            plot(xstore[1,],ystore[1,],type="l",xlim=c(0,1.2),ylim=c(0,quantile(ystore,0.97)),xlab="",ylab="Relative frequency",col=scols[1])
            if(length(slots)>1){
              for(i in 2:length(slots)) lines(xstore[i,],ystore[i,],col=scols[i])
            }
            legend('topright',legend=namey[ind],text.col=scols[1:length(slots)],bty='n')
            mtext(paste("Data summary for",deparse(substitute(Data)),sep=" "),3,font=2,line=0.25,outer=T)
            
          })

# ---- Summary of MSE object ----
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
  Misc="list"))

setMethod("initialize", "Rec", function(.Object){
  # .Object@TAC <- as.numeric(NA)
  # .Object@Effort<-1
  # .Object@Allocate<-1
  # .Object@Spatial<-c(1,1)
  .Object
})

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
            "L5", "LFS", 'Vmaxlen', 'Spatial')
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
