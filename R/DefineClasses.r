

#' Class \code{'Data'}
#' 
#' An object for storing data for analysis using data-limited methods
#' 
#' 
#' @name Data-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Data', stock)} 
#' @slot Name The name of the case-study
#' @slot Year A vector of years that correspond to catch and relative abundance data
#' @slot Cat Total annual catches
#' @slot Ind Relative abundance index 
#' @slot t The number of years corresponding to AvC and Dt 
#' @slot AvC Average catch over time t 
#' @slot Dt Depletion over time t e.g. Bnow/Bthen 
#' @slot ML Mean length time series 
#' @slot Mort Natural mortality rate 
#' @slot FMSY_M An assumed ratio of FMSY to M 
#' @slot BMSY_B0 The most productive stock size relative to unfished 
#' @slot L50 Length at 50 percent maturity 
#' @slot L95 Length at 95 percent maturity 
#' @slot Lbar Mean length of catches over Lc (modal length) 
#' @slot Lc Modal length 
#' @slot LFC Length at first capture 
#' @slot LFS smallest Length at full selection 
#' @slot CAA Catch at Age data 
#' @slot Dep Stock depletion Bnow/Bunfished (total stock) 
#' @slot Abun An estimate of absolute current vulnerable abundance
#' @slot SpAbun An estimate of absolute current spawning stock abundance
#' @slot vbK The von Bertalanffy growth coefficient 
#' @slot vbLinf Maximum length 
#' @slot vbt0 Theoretical age at length zero 
#' @slot wla Weight-Length parameter alpha 
#' @slot wlb Weight-Length parameter beta 
#' @slot steep Steepness of the Beverton Holt stock-recruitment relationship 
#' @slot CV_Cat Coefficient of variation in annual catches 
#' @slot CV_Dt Coefficient of variation in depletion over time t 
#' @slot CV_AvC Coefficient of variation in average catches over time t 
#' @slot CV_Ind Coefficient of variation in the relative abundance index 
#' @slot CV_Mort Coefficient of variation in natural mortality rate 
#' @slot CV_FMSY_M Coefficient of variation in the ratio in FMSY/M 
#' @slot CV_BMSY_B0 Coefficient of variation in the position of the most productive stock size relative to unfished 
#' @slot CV_Dep Coefficient of variation in current stock depletion 
#' @slot CV_Abun Coefficient of variation in estimate of absolute current stock size 
#' @slot CV_vbK Coefficient of variation in the von Bert. k parameter 
#' @slot CV_vbLinf Coefficient of variation in maximum length 
#' @slot CV_vbt0 Coefficient of variation in age at length zero 
#' @slot CV_L50 Coefficient of variation in length at 50 per cent maturity 
#' @slot CV_LFC Coefficient of variation in length at first capture 
#' @slot CV_LFS Coefficient of variation in length at full selection 
#' @slot CV_wla Coefficient of variation in weight-length parameter a 
#' @slot CV_wlb Coefficient of variation in weight-length parameter b 
#' @slot CV_steep Coefficient of variation in steepness   
#' @slot sigmaL Assumed observaton error of the length composition data 
#' @slot MaxAge Maximum age 
#' @slot Units Units of the catch/absolute abundance estimates 
#' @slot Ref A reference quota level 
#' @slot Ref_type Its type 
#' @slot Log A log of events 
#' @slot params A place to store estimated parameters 
#' @slot PosMPs The methods that can be applied to these data 
#' @slot MPs The methods that were applied to these data 
#' @slot OM A table of operating model conditions 
#' @slot Obs A table of observation model conditions 
#' @slot TAC The calculated TAC 
#' @slot TACbias The known bias in the calculated TAC 
#' @slot Sense The results of the sensitivity analysis 
#' @slot CAL_bins The length bins for the catch-at-length data 
#' @slot CAL Catch-at-length data 
#' @slot Cref Reference or target catch level 
#' @slot Iref Reference or target relative abundance index level 
#' @slot Bref Reference or target biomass level 
#' @slot CV_Cref CV for reference or target catch level 
#' @slot CV_Iref CV for reference or target relative abundance index level 
#' @slot CV_Bref CV for reference or target biomass level 
#' @slot CV_Rec CV for recent recruitment strength 
#' @slot Rec Recent recruitment strength 
#' @slot MPrec The previous recommendation of a management proceedure 
#' @slot MPeff The current level of effort 
#' @slot LHYear The last historical year of the simulation (before projection) 
#' @slot Misc Optional list which is passed to MPs

#' @author T. Carruthers
#' @keywords classes
#' @examples
#' 
#' newdata<-new('Data')
#' 
setClass("Data", representation(Name = "character", Year = "vector", 
  Cat = "matrix", Ind = "matrix", Rec = "matrix", t = "vector", AvC = "vector", 
  Dt = "vector", Mort = "vector", FMSY_M = "vector", BMSY_B0 = "vector", 
  Cref = "vector", Bref = "vector", Iref = "vector", L50 = "vector", 
  L95 = "vector", LFC = "vector", LFS = "vector", CAA = "array", Dep = "vector", 
  Abun = "vector", SpAbun="vector", vbK = "vector", vbLinf = "vector", vbt0 = "vector", 
  wla = "vector", wlb = "vector", steep = "vector", CV_Cat = "vector", 
  CV_Dt = "vector", CV_AvC = "vector", CV_Ind = "vector", CV_Mort = "vector", 
  CV_FMSY_M = "vector", CV_BMSY_B0 = "vector", CV_Cref = "vector", CV_Bref = "vector", 
  CV_Iref = "vector", CV_Rec = "vector", CV_Dep = "vector", CV_Abun = "vector", 
  CV_vbK = "vector", CV_vbLinf = "vector", CV_vbt0 = "vector", CV_L50 = "vector", 
  CV_LFC = "vector", CV_LFS = "vector", CV_wla = "vector", CV_wlb = "vector", 
  CV_steep = "vector", sigmaL = "vector", MaxAge = "vector", Units = "character", 
  Ref = "numeric", Ref_type = "character", Log = "list", params = "list", 
  PosMPs = "vector", MPs = "vector", OM = "data.frame", Obs = "data.frame", 
  TAC = "array", TACbias = "array", Sense = "array", CAL_bins = "numeric", 
  CAL = "array", MPrec = "vector", MPeff = "vector", ML = "array", Lbar = "array", 
  Lc = "array", LHYear = "numeric", Misc = "list"))

# initialize Data
setMethod("initialize", "Data", function(.Object, stock = "nada") {
  # .Object }) .Object<-new('Data') run an error check here
  if (file.exists(stock)) {
    dat <- read.csv(stock, header = F, colClasses = "character")  # read 1st sheet
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
    .Object@Cref <- as.numeric(dat[match("Cref", dname), 1])
    .Object@Bref <- as.numeric(dat[match("Bref", dname), 1])
    .Object@Iref <- as.numeric(dat[match("Iref", dname), 1])
    .Object@L50 <- as.numeric(dat[match("Length at 50% maturity", dname), 1])
    .Object@L95 <- as.numeric(dat[match("Length at 95% maturity", dname), 1])
    .Object@LFC <- as.numeric(dat[match("Length at first capture",  dname), 1])
    .Object@LFS <- as.numeric(dat[match("Length at full selection", dname), 1])
    .Object@Dep <- as.numeric(dat[match("Current stock depletion",  dname), 1])
    .Object@Abun <- as.numeric(dat[match("Current stock abundance",  dname), 1])
	  .Object@SpAbun <- as.numeric(dat[match("Current spawning stock abundance",  dname), 1])
    .Object@vbK <- as.numeric(dat[match("Von Bertalanffy K parameter", dname), 1])
    .Object@vbLinf <- as.numeric(dat[match("Von Bertalanffy Linf parameter", dname), 1])
    .Object@vbt0 <- as.numeric(dat[match("Von Bertalanffy t0 parameter", dname), 1])
    .Object@wla <- as.numeric(dat[match("Length-weight parameter a", dname), 1])
    .Object@wlb <- as.numeric(dat[match("Length-weight parameter b", dname), 1])
    .Object@steep <- as.numeric(dat[match("Steepness", dname), 1])
    .Object@sigmaL <- as.numeric(dat[match("Sigma length composition", dname), 1])
    .Object@CV_Cat <- as.numeric(dat[match("CV Catch", dname), 1])
    .Object@CV_Dt <- as.numeric(dat[match("CV Depletion over time t", dname), 1])
    .Object@CV_AvC <- as.numeric(dat[match("CV Average catch over time t", dname), 1])
    .Object@CV_Ind <- as.numeric(dat[match("CV Abundance index", dname), 1])
    .Object@CV_Mort <- as.numeric(dat[match("CV M", dname), 1])
    .Object@CV_Rec <- as.numeric(dat[match("CV Rec", dname), 1])
    .Object@CV_FMSY_M <- as.numeric(dat[match("CV FMSY/M", dname),  1])
    .Object@CV_BMSY_B0 <- as.numeric(dat[match("CV BMSY/B0", dname), 1])
    .Object@CV_Cref <- as.numeric(dat[match("CV Cref", dname), 1])
    .Object@CV_Bref <- as.numeric(dat[match("CV Bref", dname), 1])
    .Object@CV_Iref <- as.numeric(dat[match("CV Iref", dname), 1])
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
    .Object@MaxAge <- as.numeric(dat[match("Maximum age", dname), 1])
    .Object@MPrec <- as.numeric(dat[match("MPrec", dname), 1])
    .Object@MPeff <- as.numeric(dat[match("MPeff", dname), 1])
    
    if (length(grep("CAL", dname)) > 1) {
      CAL_bins <- as.numeric(dat[match("CAL_bins", dname), dat[match("CAL_bins", dname), ] != ""])
      nCAL <- length(CAL_bins) - 1
      .Object@CAL_bins <- CAL_bins
      CALdat <- grep("CAL ", dname)
      .Object@CAL <- array(as.numeric(as.matrix(dat[CALdat, 1:nCAL])),dim = c(1, length(CALdat), nCAL))
    }
    
    CAAy <- grep("CAA", dname)[1:length(grep("CAA", dname))]
    CAAa <- sum(dat[CAAy[1], ] != "")
    if (!is.na(CAAa)) {
      .Object@CAA <- array(as.numeric(as.matrix(dat[CAAy, 1:CAAa])),  dim = c(1, length(CAAy), CAAa))
    }
    
    .Object@ML <- matrix(as.numeric(dat[match("Mean length", dname), 1:length(.Object@Year)]), nrow = 1)
    .Object@Lbar <- matrix(as.numeric(dat[match("Mean length Lc", dname), 1:length(.Object@Year)]), nrow = 1)
    .Object@Lc <- matrix(as.numeric(dat[match("Modal length", dname), 1:length(.Object@Year)]), nrow = 1)
    
    .Object@LHYear <- as.numeric(dat[match("LHYear", dname), 1])
    .Object@Units <- dat[match("Units", dname), 1]
    .Object@Ref <- as.numeric(dat[match("Reference OFL", dname), 1])
    .Object@Ref_type <- dat[match("Reference OFL type", dname), 1]
    .Object@Log[[1]] <- paste("Created:", Sys.time())
    .Object@params <- new("list")
    .Object@OM <- data.frame(NA)
    .Object@Obs <- data.frame(NA)
    .Object@TAC <- array(NA, dim = c(1, 1, 1))
    .Object@TACbias <- array(NA, dim = c(1, 1, 1))
    .Object@Sense <- array(NA, dim = c(1, 1, 1))
    .Object@PosMPs <- NA
    .Object@MPs <- NA
    
  } else {
    if (stock != "MSE") {
      if (!is.na(stock)) print("Couldn't find specified csv file, blank DLM object created")
    }
  }
  # Default values
  # -------------------------------------------------------------
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
  if (length(.Object@sigmaL) == 0) .Object@sigmaL <- 0.2
  if (length(.Object@CAA) == 0) .Object@CAA <- array(NA, c(1, 1, 1))
  if (length(.Object@CAL) == 0) .Object@CAL <- array(NA, c(1, 1, 1))
  if (length(.Object@CAL_bins) == 0) .Object@CAL_bins <- 1
  if (length(.Object@TAC) == 0) .Object@TAC <- array(1, c(1, 1))
  if (length(.Object@TACbias) == 0) .Object@TACbias <- array(1, c(1, 1))
  if (length(.Object@Sense) == 0) .Object@Sense <- array(1, c(1, 1))
  if (length(.Object@ML) == 0)  .Object@ML <- array(NA, c(1, 1))
  if (length(.Object@Lbar) == 0) .Object@Lbar <- array(NA, c(1, 1))
  if (length(.Object@Lc) == 0) .Object@Lc <- array(NA, c(1, 1))
  .Object
})




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
#' @author T. Carruthers
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
setMethod("initialize", "Fease", function(.Object, file = "nada", ncases = 1) {
  # run an error check here
  if (file.exists(file)) {
    Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
    dat <- read.csv(file, header = F, colClasses = "character", col.names = paste0("V", 
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


#' Class \code{'Fleet'}
#' 
#' The component of the operating model that controls fishing dynamics
#' 
#' 
#' @name Fleet-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Fleet')}
#'
#' @slot Name Name of the Fleet object
#' @slot nyears The number of years for the historical simulation
#' @slot Spat_targ Distribution of fishing in relation to spatial biomass: F is proportional to B^Spat_targ (uniform distribution) 
#' @slot Esd Inter-annual variability in fishing mortality rate
#' @slot EffYears Vector of verticies, years at which to simulate varying relative effort
#' @slot EffLower Lower bound on relative effort corresponding to EffYears (uniform distribution)
#' @slot EffUpper Uppper bound on relative effort corresponding to EffYears (uniform distribution)
#' @slot LFS Shortest length that is fully vulnerable to fishing (uniform distribution)
#' @slot L5 Shortest length corresponding ot 5 percent vulnerability (uniform distribution)
#' @slot Vmaxlen The vulnerability of the longest (oldest) fish (uniform distribution)
#' @slot SelYears Vector of verticies, index for years at which historical selectivity pattern changed. If left empty, historical selectivity is constant
#' @slot AbsSelYears Optional values for SelYears, used for plotting only. Must be of same length as SelYears
#' @slot L5Lower Optional vector of values of length SelYears, specifiying lower limits of L5 (use \code{ChooseSelect} function to set these)
#' @slot L5Upper Optional vector of values of length SelYears, specifiying upper limits of L5 (use \code{ChooseSelect} function to set these)
#' @slot LFSLower Optional vector of values of length SelYears, specifiying lower limits of LFS (use \code{ChooseSelect} function to set these)
#' @slot LFSUpper Optional vector of values of length SelYears, specifiying upper limits of LFS (use \code{ChooseSelect} function to set these)
#' @slot VmaxLower Optional vector of values of length SelYears, specifiying lower limits of Vmaxlen (use \code{ChooseSelect} function to set these)
#' @slot VmaxUpper Optional vector of values of length SelYears, specifiying upper limits of Vmaxlen (use \code{ChooseSelect} function to set these)
#' @slot qinc Average percentage change in fishing efficiency (uniform distribution)(applicable only to forward projection and input controls)
#' @slot qcv Inter-annual variability in fishing efficiency (uniform distribution)(applicable only to forward projection and input controls)
#' @slot isRel Are the selectivity parameters relative to size-of-maturity? TRUE or FALSE
#' @slot CurrentYr The current calendar year (final year) of the historical simulations (e.g. 2011)
#'
#' @author T. Carruthers
#' @keywords classes
#' @examples
#' 
#' showClass('Fleet')
#' 
setClass("Fleet", slots = c(Name = "character", nyears = "numeric", Spat_targ = "numeric", 
  Esd = "numeric", qinc = "numeric", qcv = "numeric", EffYears = "numeric", 
  EffLower = "numeric", EffUpper = "numeric", SelYears = "numeric", AbsSelYears = "numeric", 
  L5 = "numeric", LFS = "numeric", Vmaxlen = "numeric", L5Lower = "numeric", 
  L5Upper = "numeric", LFSLower = "numeric", LFSUpper = "numeric", VmaxLower = "numeric", 
  VmaxUpper = "numeric", isRel = "character",CurrentYr="numeric"))

# initialize Fleet
setMethod("initialize", "Fleet", function(.Object, file = NA) {
  if (!is.na(file)) {
    if (file.exists(file)) {
      Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
      dat <- read.csv(file, header = F, colClasses = "character", col.names = paste0("V", 1:Ncol))  # read 1st sheet
      dname <- dat[, 1]
      dat <- dat[, 2:ncol(dat)]
      
      .Object@Name <- dat[match("Name", dname), 1]
      .Object@nyears <- as.numeric(dat[match("nyears", dname), 1])
      
      .Object@CurrentYr <- as.numeric(dat[match("CurrentYr", dname), 1])
      if(is.na(.Object@CurrentYr)).Object@CurrentYr<-.Object@nyears
      
      .Object@Spat_targ <- as.numeric(dat[match("Spat_targ", dname),  1:2])
      .Object@Esd <- as.numeric(dat[match("Esd", dname), 1:2])
	    .Object@Esd <- as.numeric(dat[match("Fsd", dname), 1:2])
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
      .Object@isRel <- dat[match("isRel", dname), 1]  # Are selecivity parameters relative to maturity?
      if (NAor0(.Object@isRel)) .Object@isRel <- "TRUE"
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
#' @slot Name Name of the MSE object
#' @slot nyears The number of years for the historical simulation
#' @slot proyears The number of years for the projections - closed loop simulations
#' @slot nMPs Number of management procedures simulation tested
#' @slot MPs The names of the MPs that were tested
#' @slot nsim Number of simulations
#' @slot OM A table of nsim rows with a column for each sampled parameter of the operating model\cr
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
#'   \item CALcv: variability in lengths at age around the growth curve (normal CV)
#'   \item FMSY: Fishing mortality rate at Maximum Sustainable Yield
#'   \item Linf: maximum length (von Bertalanffy Linf parameter)
#'   \item K: maximum growth rate (von Bertalanffy K parameter)
#'   \item t0: theoretical length at age zero (von Bertalanffy t0 parameter)
#'   \item hs: steepness of the stock recruitment relationship (the fraction of unfished recruitment at a fifth of unfished stock levels)
#'   \item Linfgrad: mean gradient in maximum length (per cent per time step)
#'   \item Kgrad: mean gradient in maximum growth rate (per cent per time step)
#'   \item Linfsd: interannual variability in maximum length (log normal CV)
#'   \item recgrad: gradient in recruitment strength (age 1 population numbers) over last 10 years of historical simulations
#'   \item Ksd: interannual variability in maximum growth rate (log normal CV)
#'   \item ageM: age at 50 per cent maturity
#'   \item LFS: length at full selection (the shortest length class where fishery selectivity is 100 per cent)
#'   \item age05: the age at 5 percent selectivity (ascending limb of selectivity curve)
#'   \item Vmaxage: the selectivity of the oldest age class (controls dome shape of selectivity curve)
#'   \item LFC: length at first capture, the smallest length that can be caught by the gear
#'   \item OFLreal: the true simulated Over Fishing Limit (FMSY x biomass) updated in each management update of the projection
#'   \item Spat_targ: spatial targetting parameter, fishing mortality rate across areas is proportional to vulnerable biomass raised to the power of this number. 
#'   \item Frac_area_1: the fraction of unfished biomass inhabiting area 1 (can be seen as fraction of habitat in area 1 or relative size of area 1)
#'   \item Prob_staying: the probability that individuals in area 1 remain there between time-steps
#'   \item AC: autocorrelation in recruitment
#'  }
#' @slot Obs A table of nsim rows with a column for each sampled parameter of the observation model\cr
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
#' @slot B_BMSY Stored biomass relative to BMSY over the projection (an array with dimensions nsim, nMPs, proyears) 
#' @slot F_FMSY Stored fishing mortality rate relative to FMSY over the projection (an array with dimensions nsim, nMPs, proyears) 
#' @slot B Stored stock biomass over the projection (an array with dimensions nsim, nMPs, proyears)
#' @slot SSB Stored spawning stock biomass over the projection (an array with dimensions nsim, nMPs, proyears)
#' @slot VB Stored vulnerable biomass over the projection (an array with dimensions nsim, nMPs, proyears) 
#' @slot FM Stored fishing mortality rate over the projection (an array with dimensions nsim, nMPs, proyears)
#' @slot C Stored catches (taken) over the projection (an array with dimensions nsim, nMPs, proyears) 
#' @slot TAC Stored Total Allowable Catch (prescribed) over the projection (an array with dimensions nsim, nMPs, proyears)(note that this is NA for input controls) 
#' @slot SSB_hist Stored historical spawning stock biomass (historical simulations - an array with dimensions nsim, nages, nyears, nareas)
#' @slot CB_hist Stored historical catches in weight (historical simulations - an array with dimensions nsim, nages, nyears, nareas)
#' @slot FM_hist Stored historical fishing mortality rate (historical simulations - an array with dimensions nsim, nages, nyears, nareas)
#' @slot Effort Stored relative fishing effort in the projection years
#' @slot PAA Population at age in last projection year (an array with dimensions nsim, nMPs, nages)
#' @slot CAA Catch at age in last projection year (an array with dimensions nsim, nMPs, nages)
#' @slot CAL Catch at length in last projection year (an array with dimensions nsim, nMPs, nCALbins)
#' @slot CALbins Mid-points of the catch-at-length bins
#'
#' @author T. Carruthers
#' @keywords classes
setClass("MSE", representation(Name = "character", nyears = "numeric", 
  proyears = "numeric", nMPs = "numeric", MPs = "character", nsim = "numeric", 
  OM = "data.frame", Obs = "data.frame", B_BMSY = "array", F_FMSY = "array", 
  B = "array", SSB="array", VB="array", FM = "array", C = "array", 
  TAC = "array", SSB_hist = "array", 
  CB_hist = "array", FM_hist = "array", Effort = "array", PAA= "array", CAA= "array", 
  CAL= "array", CALbins="numeric"))

  
setMethod("initialize", "MSE", function(.Object, Name, nyears, proyears, 
  nMPs, MPs, nsim, OM, Obs, B_BMSY, F_FMSY, B, SSB, VB, FM, C, TAC, 
  SSB_hist, CB_hist, FM_hist, Effort = array(), PAA,  CAA, CAL, CALbins) {
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
#' @slot Name The name of the observation model object 
#' @slot Cobs Log-normal catch observation error expressed as a coefficient of variation (uniform distribution) 
#' @slot Cbiascv A coefficient of variation controlling the sampling of bias in catch observations for each simulation (uniform distribution) 
#' @slot CAA_nsamp Number of catch-at-age observation per time step (uniform distribution) 
#' @slot CAA_ESS Effective sample size (independent age draws) of the multinomial catch-at-age observation error model (uniform distribution) 
#' @slot CAL_nsamp Number of catch-at-length observation per time step (uniform distribution) 
#' @slot CAL_ESS Effective sample size (independent length draws) of the multinomial catch-at-length observation error model (uniform distribution) 
#' @slot CALcv Lognormal, variability in the length at age (uniform distribution) 
#' @slot Iobs Observation error in the relative abundance indices expressed as a coefficient of variation (uniform distribution) 
#' @slot Mcv Persistent bias in the prescription of natural mortality rate sampled from a log-normal distribution with coefficient of variation (Mcv)(uniform distribution) 
#' @slot Kcv Persistent bias in the prescription of growth parameter k sampled from a log-normal distribution with coefficient of variation (Kcv)(uniform distribution) 
#' @slot t0cv Persistent bias in the prescription of t0 sampled from a log-normal distribution with coefficient of variation (t0cv)(uniform distribution) 
#' @slot Linfcv Persistent bias in the prescription of maximum length sampled from a log-normal distribution with coefficient of variation (Linfcv)(uniform distribution) 
#' @slot LFCcv Persistent bias in the prescription of lenght at first capture sampled from a log-normal distribution with coefficient of variation (LFCcv)(uniform distribution) 
#' @slot LFScv Persistent bias in the prescription of length-at-fully selection sampled from a log-normal distribution with coefficient of variation (LFScv)(uniform distribution) 
#' @slot B0cv Persistent bias in the prescription of maximum lengthunfished biomass sampled from a log-normal distribution with coefficient of variation (B0cv)(uniform distribution) 
#' @slot FMSYcv Persistent bias in the prescription of FMSY sampled from a log-normal distribution with coefficient of variation (FMSYcv)(uniform distribution) 
#' @slot FMSY_Mcv Persistent bias in the prescription of FMSY/M sampled from a log-normal distribution with coefficient of variation (FMSY_cv)(uniform distribution) 
#' @slot BMSY_B0cv Persistent bias in the prescription of BMsY relative to unfished sampled from a log-normal distribution with coefficient of variation (BMSY_B0cv)(uniform distribution) 
#' @slot rcv Persistent bias in the prescription of intrinsic rate of increase sampled from a log-normal distribution with coefficient of variation (rcv)(uniform distribution) 
#' @slot LenMcv Persistent bias in the prescription of length at 50 percent maturity sampled from a log-normal distribution with coefficient of variation (A50cv)(uniform distribution) 
#' @slot Dbiascv Persistent bias in the prescription of stock depletion sampled from a log-normal distribution with coefficient of variation (Linfcv)(uniform distribution) 
#' @slot Dcv Imprecision in the prescription of stock depletion among years, expressed as a coefficient of variation (uniform distribution) 
#' @slot Btbias Persistent bias in the prescription of current stock biomass sampled from a uniform-log distribution with range (Btbias)(uniform distribution) 
#' @slot Btcv Imprecision in the prescription of current stock biomass among years expressed as a coefficient of variation (uniform distribution) 
#' @slot Fcurbiascv Persistent bias in the prescription of current fishing mortality rate sampled from a log-normal distribution with coefficient of variation (Fcurcv)(uniform distribution) 
#' @slot Fcurcv Imprecision in the prescription of current fishing mortality rate among years expressed as a coefficient of variation (uniform distribution) 
#' @slot hcv Persistent bias in steepness (uniform distribution) 
#' @slot Icv Observation error in realtive abundance index expressed as a coefficient of variation (uniform distirbution) 
#' @slot maxagecv Bias in the prescription of maximum age (uniform distribution) 
#' @slot beta A parameter controlling hyperstability/hyperdepletion. I^beta therefore values below 1 lead to hyperstability (an index that decreases slower than true abundance) and values above 1 lead to hyperdepletion (an index that decreases more rapidly than true abundance)(uniform distribution) 
#' @slot Reccv Bias in the knowledge of recent recruitment strength (uniform distribution) 
#' @slot Irefcv Bias in the knowledge of the relative abundance index at BMSY (uniform distribution) 
#' @slot Brefcv Bias in the knowledge of BMSY (uniform distribution) 
#' @slot Crefcv Bias in the knowledge of MSY(uniform distribution) 
#'     
#' @author T. Carruthers
#' @keywords classes
#' @examples
#' 
#' showClass('Obs')
#' 
setClass("Obs", representation(Name = "character", LenMcv = "numeric", 
  Cobs = "numeric", Cbiascv = "numeric", CAA_nsamp = "numeric", CAA_ESS = "numeric", 
  CAL_nsamp = "numeric", CAL_ESS = "numeric", CALcv = "numeric", Iobs = "numeric", 
  Mcv = "numeric", Kcv = "numeric", t0cv = "numeric", Linfcv = "numeric", 
  LFCcv = "numeric", LFScv = "numeric", B0cv = "numeric", FMSYcv = "numeric", 
  FMSY_Mcv = "numeric", BMSY_B0cv = "numeric", rcv = "numeric", Dbiascv = "numeric", 
  Dcv = "numeric", Btbias = "numeric", Btcv = "numeric", Fcurbiascv = "numeric", 
  Fcurcv = "numeric", hcv = "numeric", Icv = "numeric", maxagecv = "numeric", 
  Reccv = "numeric", Irefcv = "numeric", Crefcv = "numeric", Brefcv = "numeric", 
  beta = "numeric"))

# initialize Obs
setMethod("initialize", "Obs", function(.Object, file = NA) {
  if (!is.na(file)) {
    if (file.exists(file)) {
      Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
      dat <- read.csv(file, header = F, colClasses = "character", 
        col.names = paste0("V", 1:Ncol))  # read 1st sheet
      dname <- dat[, 1]
      dat <- dat[, 2:ncol(dat)]
      .Object@Name <- dat[match("Name", dname), 1]
      .Object@LenMcv <- as.numeric(dat[match("LenMcv", dname), 1])
      .Object@Cobs <- as.numeric(dat[match("Cobs", dname), 1:2])
      .Object@Cbiascv <- as.numeric(dat[match("Cbiascv", dname), 1])
      .Object@CAA_nsamp <- as.numeric(dat[match("CAA_nsamp", dname), 1:2])
      .Object@CAA_ESS <- as.numeric(dat[match("CAA_ESS", dname), 1:2])
      .Object@CAL_nsamp <- as.numeric(dat[match("CAA_nsamp", dname), 1:2])
      .Object@CAL_ESS <- as.numeric(dat[match("CAA_ESS", dname), 1:2])
      .Object@CALcv <- as.numeric(dat[match("CALcv", dname), 1:2])
      .Object@Iobs <- as.numeric(dat[match("Iobs", dname), 1:2])
      .Object@Mcv <- as.numeric(dat[match("Mcv", dname), 1])
      .Object@Kcv <- as.numeric(dat[match("Kcv", dname), 1])
      .Object@t0cv <- as.numeric(dat[match("t0cv", dname), 1])
      .Object@Linfcv <- as.numeric(dat[match("Linfcv", dname), 1])
      .Object@LFCcv <- as.numeric(dat[match("LFCcv", dname), 1])
      .Object@LFScv <- as.numeric(dat[match("LFScv", dname), 1])
      .Object@B0cv <- as.numeric(dat[match("B0cv", dname), 1])
      .Object@FMSYcv <- as.numeric(dat[match("FMSYcv", dname), 1])
      .Object@FMSY_Mcv <- as.numeric(dat[match("FMSY_Mcv", dname), 1])
      .Object@BMSY_B0cv <- as.numeric(dat[match("BMSY_B0cv", dname), 1])
      .Object@rcv <- as.numeric(dat[match("rcv", dname), 1])
      .Object@Dbiascv <- as.numeric(dat[match("Dbiascv", dname), 1])
      .Object@Dcv <- as.numeric(dat[match("Dcv", dname), 1:2])
      .Object@Btbias <- as.numeric(dat[match("Btbias", dname), 1:2])
      .Object@Btcv <- as.numeric(dat[match("Btcv", dname), 1:2])
      .Object@Fcurbiascv <- as.numeric(dat[match("Fcurbiascv", dname), 1])
      .Object@Fcurcv <- as.numeric(dat[match("Fcurcv", dname), 1:2])
      .Object@hcv <- as.numeric(dat[match("hcv", dname), 1])
      .Object@Icv <- as.numeric(dat[match("Icv", dname), 1])
      .Object@maxagecv <- as.numeric(dat[match("maxagecv", dname), 1])
      .Object@Reccv <- as.numeric(dat[match("Reccv", dname), 1:2])
      .Object@Irefcv <- as.numeric(dat[match("Irefcv", dname), 1])
      .Object@Crefcv <- as.numeric(dat[match("Crefcv", dname), 1])
      .Object@Brefcv <- as.numeric(dat[match("Brefcv", dname), 1])
      .Object@beta <- as.numeric(dat[match("beta", dname), 1:2])
    } else {
      message("File doesn't exist")
    }
  }
  .Object
  
})





#' Class \code{'OM'}
#' 
#' An object containing all the parameters needed to control the MSE which can
#' be build from component Stock, Fleet and Obs objects. Almost all of
#' these inputs are a vector of length 2 which describes the upper and lower
#' bounds of a uniform distribution from which to sample the parameter.
#' 
#' 
#' @name OM-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('OM', Stock, Fleet, Obs, Imp)}. 

#' @slot Name Name of the operating model
#' @slot nsim The number of simulations
#' @slot proyears The number of projected years

#' @slot nyears The number of years for the historical simulation 
#' @slot maxage The maximum age of individuals that is simulated (there is no 'plus group': individuals die off beyone the maximum age so there isn't a huge cost to simulating more older age classes) 
#' @slot R0 The magnitude of unfished recruitment. This is normally fixed to some arbitrary value since it simply scales the simulated numbers) 
#' @slot M Natural mortality rate (uniform distribution) 
#' @slot Msd Inter-annual variability in natural mortality rate expressed as a coefficient of variation (uniform distribution) 
#' @slot Mgrad Mean temporal trend in natural mortality rate, expressed as a percentage change in M per year (uniform distribution) 
#' @slot h Steepness of the stock recruit relationship (uniform distribution) 
#' @slot SRrel Type of stock-recruit relationship (1)Beverton-Holt (2) Ricker 
#' @slot Linf Maximum length (uniform distribution) 
#' @slot K von B. growth parameter k (uniform distribution) 
#' @slot t0 von B. theoretical age at length zero (uniform distribution) 
#' @slot Ksd Inter-annual variability in growth parameter k (uniform distribution) 
#' @slot Kgrad Mean temporal trend in growth parameter k, expressed as a percentage change in k per year (uniform distribution) 
#' @slot Linfsd Inter-annual variability in maximum length - uniform distribution 
#' @slot Linfgrad Mean temporal trend in maximum length, expressed as a percentage change in Linf per year (uniform distribution) 
#' @slot recgrad Mean temporal trend in log-normal recruitment deviations (uniform distribution) 
#' @slot AC Autocorrelation in recruitment deviations rec(t)=AC*rec(t-1)+(1-AC)*sigma(t) (uniform distribution) 
#' @slot a Length-weight parameter alpha (uniform distribution) 
#' @slot b Length-weight parameter beta (uniform distribution) 
#' @slot D Current level of stock depletion (Bcurrent/Bunfished) (uniform distribution) 
#' @slot Size_area_1 The size of area 1 relative to area 2 (uniform distribution) 
#' @slot Frac_area_1 The fraction of the unfished biomass in stock 1 (uniform distribution) 
#' @slot Prob_staying The probability of inviduals in area 1 remaining in area 1 over the course of one year 

#' @slot beta A parameter controlling hyperstability/hyperdepletion. I^beta therefore values below 1 lead to hyperstability (an index that decreases slower than true abundance) and values above 1 lead to hyperdepletion (an index that decreases more rapidly than true abundance)(uniform distribution) 
#' @slot Spat_targ Distribution of fishing in relation to spatial biomass: F is proportional to B^Spat_targ (uniform distribution)  
#' @slot LFS Shortest length that is fully vulnerable to fishing (uniform distribution) 
#' @slot L5 Shortest length at 5 percent vulnerability (uniform distribution) 
#' @slot Vmaxlen The vulnerability of the longest (oldest) fish (uniform distribution) 
#' @slot SelYears Vector of verticies that index years where historical selectivity pattern changed. Leave empty to ignore 
#' @slot AbsSelYears vector of absolute year values that correspond to year indices in SelYears. Used only for plotting 
#' @slot L5Lower Optional vector of values of length SelYears, specifiying lower limits of L5 (use \code{ChooseSelect}  function to set these. Overrides L5 above) 
#' @slot L5Upper Optional vector of values of length SelYears, specifiying upper limits of L5 (use \code{ChooseSelect}  function to set these. Overrides L5 above) 
#' @slot LFSLower Optional vector of values of length SelYears, specifiying lower limits of LFS (use \code{ChooseSelect}  function to set these. Overrides LFS above) 
#' @slot LFSUpper Optional vector of values of length SelYears, specifiying upper limits of LFS (use \code{ChooseSelect}  function to set these. Overrides LFS above) 
#' @slot VmaxLower Optional vector of values of length SelYears, specifiying lower limits of Vmaxlen (use \code{ChooseSelect}  function to set these. Overrides Vmaxlen above) 
#' @slot VmaxUpper Optional vector of values of length SelYears, specifiying upper limits of Vmaxlen (use \code{ChooseSelect}  function to set these. Overrides Vmaxlen above) 
#' @slot isRel Are the selectivity parameters relative to size-of-maturity? TRUE or FALSE 
#' @slot L50 Length at 50 percent maturity (uniform distribution) 
#' @slot L50_95 Length increment from 50 to 95 percent maturity (uniform distribution) 
#' @slot Esd Inter-annual variability in fishing mortality rate 
#' @slot EffYears Vector of verticies, years at which to simulate varying relative effort 
#' @slot EffLower Lower bound on relative effort corresponding to EffYears (uniform distribution) 
#' @slot EffUpper Uppper bound on relative effort corresponding to EffYears (uniform distribution) 
#' @slot qinc Average percentage change in fishing efficiency (uniform distribution)(applicable only to forward projection and input controls) 
#' @slot qcv Inter-annual variability in fishing efficiency (uniform distribution)(applicable only to forward projection and input controls) 
#' @slot CurrentYr The current calendar year (final year) of the historical simulations (e.g. 2011)

#' @slot Cobs Log-normal catch observation error expressed as a coefficient of variation (uniform distribution) 
#' @slot Cbiascv A coefficient of variation controlling the sampling of bias in catch observations for each simulation (uniform distribution) 
#' @slot CAA_nsamp Number of catch-at-age observation per time step (uniform distribution) 
#' @slot CAA_ESS Effective sample size (independent age draws) of the multinomial catch-at-age observation error model (uniform distribution) 
#' @slot CAL_nsamp Number of catch-at-length observation per time step (uniform distribution) 
#' @slot CAL_ESS Effective sample size (independent length draws) of the multinomial catch-at-length observation error model (uniform distribution) 
#' @slot CALcv Lognormal, variability in the length at age (uniform distribution) 
#' @slot Iobs Observation error in the relative abundance indices expressed as a coefficient of variation (uniform distribution) 
#' @slot Perr The extent of inter-annual log-normal recruitment variability (sigma R)(uniform distribution) 
#' @slot Period Period for cylical recruitment pattern in years (uniform distribution). Leave empty to ignore  
#' @slot Amplitude Amplitude in deviation from long-term average recruitment during recruitment cycle, both positive and negative (uniform distribution). E.g., a range from 0 to 0.5 means recruitment decreases or increases by up to 50\% each cycle. Leave empty to ignore 
#' @slot Mcv Persistent bias in the prescription of natural mortality rate sampled from a log-normal distribution with coefficient of variation (Mcv)(uniform distribution) 
#' @slot Kcv Persistent bias in the prescription of growth parameter k sampled from a log-normal distribution with coefficient of variation (Kcv)(uniform distribution) 
#' @slot t0cv Persistent bias in the prescription of t0 sampled from a log-normal distribution with coefficient of variation (t0cv)(uniform distribution) 
#' @slot Linfcv Persistent bias in the prescription of maximum length sampled from a log-normal distribution with coefficient of variation (Linfcv)(uniform distribution) 
#' @slot LFCcv Persistent bias in the prescription of lenght at first capture sampled from a log-normal distribution with coefficient of variation (LFCcv)(uniform distribution) 
#' @slot LFScv Persistent bias in the prescription of length-at-fully selection sampled from a log-normal distribution with coefficient of variation (LFScv)(uniform distribution) 
#' @slot B0cv Persistent bias in the prescription of maximum lengthunfished biomass sampled from a log-normal distribution with coefficient of variation (B0cv)(uniform distribution) 
#' @slot FMSYcv Persistent bias in the prescription of FMSY sampled from a log-normal distribution with coefficient of variation (FMSYcv)(uniform distribution) 
#' @slot FMSY_Mcv Persistent bias in the prescription of FMSY/M sampled from a log-normal distribution with coefficient of variation (FMSY_cv)(uniform distribution) 
#' @slot BMSY_B0cv Persistent bias in the prescription of BMsY relative to unfished sampled from a log-normal distribution with coefficient of variation (BMSY_B0cv)(uniform distribution) 
#' @slot rcv Persistent bias in the prescription of intrinsic rate of increase sampled from a log-normal distribution with coefficient of variation (rcv)(uniform distribution) 
#' @slot LenMcv Persistent bias in the prescription of length at 50 percent maturity sampled from a log-normal distribution with coefficient of variation (A50cv)(uniform distribution) 
#' @slot Dbiascv Persistent bias in the prescription of stock depletion sampled from a log-normal distribution with coefficient of variation (Linfcv)(uniform distribution) 
#' @slot Dcv Imprecision in the prescription of stock depletion among years, expressed as a coefficient of variation (uniform distribution) 
#' @slot Btbias Persistent bias in the prescription of current stock biomass sampled from a uniform-log distribution with range (Btbias)(uniform distribution) 
#' @slot Btcv Imprecision in the prescription of current stock biomass among years expressed as a coefficient of variation (uniform distribution) 
#' @slot Fcurbiascv Persistent bias in the prescription of current fishing mortality rate sampled from a log-normal distribution with coefficient of variation (Fcurcv)(uniform distribution) 
#' @slot Fcurcv Imprecision in the prescription of current fishing mortality rate among years expressed as a coefficient of variation (uniform distribution) 
#' @slot hcv Persistent bias in steepness (uniform distribution) 
#' @slot Icv Observation error in realtive abundance index expressed as a coefficient of variation (uniform distirbution) 
#' @slot maxagecv Bias in the prescription of maximum age (uniform distribution) 
#' @slot Reccv Bias in the knowledge of recent recruitment strength (uniform distribution) 
#' @slot Irefcv Bias in the knowledge of the relative abundance index at BMSY (uniform distribution) 
#' @slot Brefcv Bias in the knowledge of BMSY (uniform distribution) 
#' @slot Crefcv Bias in the knowledge of MSY(uniform distribution) 

#' @slot cpars A list of custom parameters (single parameters are a vector nsim long, time series are a matrix nsim rows by nyears columns)
#' @slot seed A random seed to ensure users can reproduce results exactly
#' @slot Source A reference to a website or article form which parameters were taken to define the operating model 

#' @slot TACSD lognormal standard deviation in fraction of TAC taken (uniform distribution) 
#' @slot TACFrac Mean fraction of TAC taken (uniform distribution) (can be an improper fraction greater than 1)
#' @slot ESD lognormal standard deviation in fraction of TAE taken(uniform distribution)
#' @slot EFrac Mean fraction of recommended effort taken (uniform distribution)
#' @slot SizeLimSD lognormal error in size limit implementation (uniform distribution)
#' @slot SizeLimFrac Mean fraction of the size limit (uniform distribution) (can be an improper fraction greater than 1)
#' @slot DiscMort Discard mortality rate (uniform distribution) (can be an improper fraction greater than 1)
#' 
#' @author T. Carruthers
#' @keywords classes
#' @examples
#' 
#' showClass('OM')
#' 
setClass("OM", representation(Name = "character", nsim="numeric",proyears="numeric",
  nyears = "numeric", maxage = "numeric", 
  R0 = "numeric", M = "numeric", Msd = "numeric", Mgrad = "numeric", 
  h = "numeric", SRrel = "numeric", Linf = "numeric", K = "numeric", 
  t0 = "numeric", Ksd = "numeric", Kgrad = "numeric", Linfsd = "numeric", 
  Linfgrad = "numeric", recgrad = "numeric", a = "numeric", b = "numeric", 
  D = "numeric", Size_area_1 = "numeric", Frac_area_1 = "numeric", Prob_staying = "numeric", 
  Source = "character", L50 = "numeric", L50_95 = "numeric", SelYears = "numeric", 
  AbsSelYears = "numeric", L5 = "numeric", LFS = "numeric", Vmaxlen = "numeric", 
  L5Lower = "numeric", L5Upper = "numeric", LFSLower = "numeric", LFSUpper = "numeric", 
  VmaxLower = "numeric", VmaxUpper = "numeric", isRel = "character", 
  beta = "numeric", Spat_targ = "numeric", Esd = "numeric",
  Period = "numeric", Amplitude = "numeric", EffYears = "numeric", EffLower = "numeric", 
  EffUpper = "numeric", qinc = "numeric", qcv = "numeric", AC = "numeric", 
  
  Cobs = "numeric", Cbiascv = "numeric", CAA_nsamp = "numeric", CAA_ESS = "numeric", 
  CAL_nsamp = "numeric", CAL_ESS = "numeric", CALcv = "numeric", Iobs = "numeric", 
  Perr = "numeric", Mcv = "numeric", Kcv = "numeric", t0cv = "numeric", 
  Linfcv = "numeric", LFCcv = "numeric", LFScv = "numeric", B0cv = "numeric", 
  FMSYcv = "numeric", FMSY_Mcv = "numeric", BMSY_B0cv = "numeric", LenMcv = "numeric", 
  rcv = "numeric", Dbiascv = "numeric", Dcv = "numeric", Btbias = "numeric", 
  Btcv = "numeric", Fcurbiascv = "numeric", Fcurcv = "numeric", hcv = "numeric", 
  Icv = "numeric", maxagecv = "numeric", Reccv = "numeric", Irefcv = "numeric", 
  Crefcv = "numeric", Brefcv = "numeric",
  
  TACSD = "numeric", TACFrac = "numeric", 
  ESD = "numeric", EFrac = "numeric",
  SizeLimSD = "numeric", SizeLimFrac="numeric",
  DiscMort = "numeric", 
  
  cpars="list",seed="numeric",CurrentYr="numeric"))


# This hack is neccessary to get around CRAN checks but causes problems that 
# Generic_fleet etc objects are considered class NULL ONLY when accessed by internal functions
# in the package - i.e., they are ok in interactive mode
# Commented out for now but need a fix before we release the next version to CRAN. 

# Generic_fleet<-NULL
# Generic_obs<-NULL
# Perfect_Imp<-NULL
# testOM<-NULL


# initialize OM
setMethod("initialize", "OM", function(.Object, Stock=NULL, Fleet=DLMtool::Generic_fleet, 
                                       Obs=DLMtool::Generic_obs, Imp=DLMtool::Perfect_Imp, 
                                       nsim=48, proyears=50) {
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
  
  # Default MSE parameters
  .Object@nsim <- nsim       
  .Object@proyears <- proyears
  
  if(length(.Object@CurrentYr)==0).Object@CurrentYr=.Object@nyears
  
  .Object@seed=1
  .Object
})



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




#' Plot Data object
#'
#' @rdname plot-Data 
#' @param x object of class Data
#' @param funcs MPs 
#' @param maxlines maximum number of lines
#' @param perc percentile of TAC recommendation
#' @param xlims limits of x-axis
#' @export
setMethod("plot",
  signature(x = "Data"),
  function(x,funcs=NA,maxlines=6,perc=0.5,xlims=NA){
   
    old_par <- par(no.readonly = TRUE)
    on.exit(par(list = old_par), add = TRUE)
	  
    Data<-x
	if (class(Data) != "Data") stop("Must supply object of class Data")
	if (all(is.na(Data@TAC))) stop("No TAC data found")
    cols<-rep(c('black','red','green','blue','orange','brown','purple','dark grey','violet','dark red','pink','dark blue','grey'),4)
    ltys<-rep(1:4,each=13)
    
    if(is.na(funcs[1]))funcs<-Data@MPs

    nMPs<-length(funcs)
    nplots<-ceiling(nMPs/maxlines)
    maxl<-ceiling(nMPs/nplots)
    mbyp <- split(1:nMPs, ceiling(1:nMPs/maxl))   # assign methods to plots

    if(is.na(xlims[1])|length(xlims)!=2){
      xlims<-quantile(Data@TAC,c(0.005,0.95),na.rm=T)
      if(xlims[1]<0)xlims[1]<-0
    }
    if(!NAor0(Data@Ref)){
      if(xlims[1]>Data@Ref)xlims[1]<-max(0,0.98*Data@Ref)
      if(xlims[2]<Data@Ref)xlims[2]<-1.02*Data@Ref
    }
    ylims<-c(0,1)

    #for(m in 1:nMPs){
     # if(sum(!is.na(Data@TAC[m,,1]))>2){
       # dens<-density(Data@TAC[m,,1],na.rm=T)
        #print(quantile(dens$y,0.99,na.rm=T))
      #  if(quantile(dens$y,0.9,na.rm=T)>ylims[2])ylims[2]<-quantile(dens$y,0.90,na.rm=T)
      #}
    #}

    #dev.new2(width=10,height=0.5+7*nplots)
    par(mfrow=c(ceiling(nplots/2),2),mai=c(0.4,0.4,0.01,0.01),omi=c(0.35,0.35,0.35,0.05))

    for(p in 1:nplots){
      m<-mbyp[[p]][1]
      plot(NA,NA,xlim=xlims,ylim=ylims,main="",xlab="",ylab="",col="white",lwd=3,type="l")
      abline(h=0)
      if(!NAor0(Data@Ref)){
        abline(v=Data@Ref,col="light grey",lwd=2)
        if(!NAor0(Data@Ref_type[1]))legend('right',Data@Ref_type,text.col="grey",bty='n')
      }
      #plot(density(DLM@TAC[m,,1],from=0,na.rm=T),xlim=xlims,ylim=ylims,main="",xlab="",ylab="",col=coly[m],lty=ltyy[m],type="l")

      if(!is.na(perc[1]))abline(v=quantile(Data@TAC[m,,1],p=perc,na.rm=T),col=cols[m],lty=ltys[m])
      #if(length(mbyp[[p]])>0){
        for(ll in 1:length(mbyp[[p]])){
          m<-mbyp[[p]][ll]
          if(sum(!is.na(Data@TAC[m,,1]))>10){  # only plot if there are sufficient non-NA TAC samples
            x<-density(Data@TAC[m,,1],from=0,na.rm=T)$x
            y<-density(Data@TAC[m,,1],from=0,na.rm=T)$y
            y<-y/max(y)
            lines(x,y,col=cols[ll])
          }else{
            print(paste("Method ",funcs[m]," produced too many NA TAC values for plotting densities",sep=""))
          }
          if(!is.na(perc[1]))abline(v=quantile(Data@TAC[m,,1],p=perc,na.rm=T),col=cols[ll],lty=2)
        }
      #}
      cind<-1:length(mbyp[[p]])
      legend('topright',funcs[mbyp[[p]]],text.col=cols[cind],col=cols[cind],lty=1,bty='n',cex=0.75)
    }

    mtext(paste("TAC (",Data@Units,")",sep=""),1,outer=T,line=0.5)
    mtext(paste("Standardized relative frequency",sep=""),2,outer=T,line=0.5)
    mtext(paste("TAC calculation for ",Data@Name,sep=""),3,outer=T,line=0.5)
})

#' Plot MSE object
#'
#' @rdname plot-MSE
#' @param x object of class MSE
#' @export
setMethod("plot",
  signature(x = "MSE"),
  function(x){
  MSEobj<-x
  Pplot(MSEobj)
  Kplot(MSEobj)
  Tplot(MSEobj)
})


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
#' @slot Name The name of the Stock object 
#' @slot maxage The maximum age of individuals that is simulated (there is no 'plus group': individuals die off beyone the maximum age so there isn't a huge cost to simulating more older age classes) 
#' @slot R0 The magnitude of unfished recruitment. This is normally fixed to some arbitrary value since it simply scales the simulated numbers) 
#' @slot M Natural mortality rate (uniform distribution) 
#' @slot Msd Inter-annual variability in natural mortality rate expressed as a coefficient of variation (uniform distribution) 
#' @slot Mgrad Mean temporal trend in natural mortality rate, expressed as a percentage change in M per year (uniform distribution) 
#' @slot h Steepness of the stock recruit relationship (uniform distribution) 
#' @slot SRrel Type of stock-recruit relationship (1)Beverton-Holt (2) Ricker 
#' @slot Linf Maximum length (uniform distribution) 
#' @slot K von B. growth parameter k (uniform distribution) 
#' @slot t0 von B. theoretical age at length zero (uniform distribution) 
#' @slot Ksd Inter-annual variability in growth parameter k (uniform distribution) 
#' @slot Kgrad Mean temporal trend in growth parameter k, expressed as a percentage change in k per year (uniform distribution) 
#' @slot Linfsd Inter-annual variability in maximum length - uniform distribution 
#' @slot Linfgrad Mean temporal trend in maximum length, expressed as a percentage change in Linf per year (uniform distribution) 
#' @slot recgrad Mean temporal trend in log-normal recruitment deviations (uniform distribution) 
#' @slot AC Autocorrelation in recruitment deviations rec(t)=AC*rec(t-1)+(1-AC)*sigma(t) (uniform distribution) 
#' @slot a Length-weight parameter alpha (uniform distribution) 
#' @slot b Length-weight parameter beta (uniform distribution) 
#' @slot L50 Length-at- 50 percent maturity (uniform distribution) 
#' @slot L50_95 Length increment from 50 percent to 95 percent maturity 
#' @slot D Current level of stock depletion (Bcurrent/Bunfished) (uniform distribution) 
#' @slot Perr Process error, the CV of lognormal recruitment deviations  (uniform distribution) 
#' @slot Period Period for cylical recruitment pattern in years (uniform distribution). Leave empty to ignore  
#' @slot Amplitude Amplitude in deviation from long-term average recruitment during recruitment cycle, both positive and negative (uniform distribution). E.g., a range from 0 to 0.5 means recruitment decreases or increases by up to 50\% each cycle. Leave empty to ignore 
#' @slot Size_area_1 The size of area 1 relative to area 2 (uniform distribution) 
#' @slot Frac_area_1 The fraction of the unfished biomass in stock 1 (uniform distribution) 
#' @slot Prob_staying The probability of inviduals in area 1 remaining in area 1 over the course of one year 
#' @slot Source A reference to a website or article form which parameters were taken to define the operating model

#' @author T. Carruthers
#' @keywords classes
#' @examples
#' 
#' showClass('Stock')
#' 
setClass("Stock", representation(Name = "character", maxage = "numeric", 
  R0 = "numeric", M = "numeric", Msd = "numeric", Mgrad = "numeric", 
  h = "numeric", SRrel = "numeric", Linf = "numeric", K = "numeric", 
  t0 = "numeric", Ksd = "numeric", Kgrad = "numeric", Linfsd = "numeric", 
  Linfgrad = "numeric", recgrad = "numeric", a = "numeric", b = "numeric", 
  D = "numeric", Perr = "numeric", Period = "numeric", Amplitude = "numeric", 
  Size_area_1 = "numeric", Frac_area_1 = "numeric", Prob_staying = "numeric", 
  AC = "numeric", L50 = "numeric", L50_95 = "numeric", Source = "character"))

# initialize Stock
setMethod("initialize", "Stock", function(.Object, file = NA) {
  
  if (!is.na(file)) {
    if (file.exists(file)) {
      Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
      dat <- read.csv(file, header = F, colClasses = "character", col.names = paste0("V", 1:Ncol))  # read 1st sheet
      dname <- dat[, 1]
      dat <- dat[, 2:ncol(dat)]
      
      .Object@Name <- dat[match("Name", dname), 1]
      .Object@maxage <- as.numeric(dat[match("maxage", dname), 1])
      .Object@R0 <- as.numeric(dat[match("R0", dname), 1])
      .Object@M <- as.numeric(dat[match("M", dname), 1:2])
      .Object@Msd <- as.numeric(dat[match("Msd", dname), 1:2])
      .Object@Mgrad <- as.numeric(dat[match("Mgrad", dname), 1:2])
      .Object@h <- as.numeric(dat[match("h", dname), 1:2])
      .Object@SRrel <- as.numeric(dat[match("SRrel", dname), 1])
      .Object@Linf <- as.numeric(dat[match("Linf", dname), 1:2])
      .Object@K <- as.numeric(dat[match("K", dname), 1:2])
      .Object@t0 <- as.numeric(dat[match("t0", dname), 1:2])
      .Object@Ksd <- as.numeric(dat[match("Ksd", dname), 1:2])
      .Object@Kgrad <- as.numeric(dat[match("Kgrad", dname), 1:2])
      .Object@Linfsd <- as.numeric(dat[match("Linfsd", dname), 1:2])
      .Object@Linfgrad <- as.numeric(dat[match("Linfgrad", dname),  1:2])
      .Object@recgrad <- as.numeric(dat[match("recgrad", dname), 1:2])
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
      .Object@Source <- dat[match("Source", dname), 1]
    } else {
      message("File doesn't exist")
    }
  }
  .Object
  
})



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

#' Summary of MSE object
#'
#' @param object object of class MSE
#' @rdname summary-MSE
#' @export
setMethod("summary",
          signature(object = "MSE"),
          function(object){            

    MSEobj<-object      
    nm<-MSEobj@nMPs
    nsim<-MSEobj@nsim
    proyears<-MSEobj@proyears
    
    Yd<-P10<-P50<-P100<-POF<-LTY<-STY<-VY<-array(NA,c(nm,nsim))
    
    yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
    RefYd<-MSEobj@OM$RefY
    yend<-max(MSEobj@proyears-9,1):MSEobj@proyears
    ystart<-1:10
    y1<-1:(MSEobj@proyears-1)
    y2<-2:MSEobj@proyears
    
    for(m in 1:nm){
      Yd[m,]<-round(apply(MSEobj@C[,m,yind],1,mean,na.rm=T)/RefYd*100,1)
      POF[m,]<-round(apply(MSEobj@F_FMSY[,m,]>1,1,sum,na.rm=T)/proyears*100,1)
      P10[m,]<-round(apply(MSEobj@B_BMSY[,m,]<0.1,1,sum,na.rm=T)/proyears*100,1)
      P50[m,]<-round(apply(MSEobj@B_BMSY[,m,]<0.5,1,sum,na.rm=T)/proyears*100,1)
      P100[m,]<-round(apply(MSEobj@B_BMSY[,m,]<1,1,sum,na.rm=T)/proyears*100,1)
      LTY[m]<-round(sum(MSEobj@C[,m,yend]/RefYd>0.5)/(MSEobj@nsim*length(yend))*100,1)
      STY[m]<-round(sum(MSEobj@C[,m,ystart]/RefYd>0.5)/(MSEobj@nsim*length(ystart))*100,1)
      AAVY<-apply(((MSEobj@C[,m,y1]-MSEobj@C[,m,y2])^2)^0.5,1,mean)/apply(MSEobj@C[,m,y2],1,mean)
      VY[m]<-round(sum(AAVY<0.1)/MSEobj@nsim*100,1)
    }
    nr<-2
    out<-cbind(MSEobj@MPs,round(apply(Yd,1,mean,na.rm=T),nr),round(apply(Yd,1,sd,na.rm=T),nr),
                             round(apply(POF,1,mean,na.rm=T),nr),round(apply(POF,1,sd,na.rm=T),nr),
                             round(apply(P10,1,mean,na.rm=T),nr),round(apply(P10,1,sd,na.rm=T),nr),
                             round(apply(P50,1,mean,na.rm=T),nr),round(apply(P50,1,sd,na.rm=T),nr),
                             round(apply(P100,1,mean,na.rm=T),nr),round(apply(P100,1,sd,na.rm=T),nr),
                             round(apply(LTY,1,mean,na.rm=T),nr),
                             round(apply(STY,1,mean,na.rm=T),nr),
                             round(apply(VY,1,mean,na.rm=T),nr))
    out<-as.data.frame(out)
    names(out)<-c("MP","Yield","stdev","POF","stdev ","P10","stdev",
                  "P50","stdev","P100","stdev","LTY","STY","VY")
    out[,1]<-as.character(out[,1])
    for(i in 2:ncol(out))out[,i]<-as.numeric(as.character(out[,i]))
    out
  })



#' Class \code{'Imp'}
#' 
#' An operating model component that specifies the degree of adherence to management recommendations (Implementation error)
#' 
#' 
#' @name Imp-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('Imp')}
#' @slot Name The name of the Implementation error object 
#' @slot TACSD lognormal standard deviation in fraction of TAC taken (uniform distribution) 
#' @slot TACFrac Mean fraction of TAC taken (uniform distribution) (can be an improper fraction greater than 1)
#' @slot ESD lognormal standard deviation in fraction of TAE taken(uniform distribution)
#' @slot EFrac Mean fraction of recommended effort taken (uniform distribution)
#' @slot SizeLimSD lognormal error in size limit implementation (uniform distribution)
#' @slot SizeLimFrac Mean fraction of the size limit (uniform distribution) (can be an improper fraction greater than 1)
#' @slot DiscMort Discard mortality rate (uniform distribution) (can be an improper fraction greater than 1)
#' @slot Source A reference to a website or article form which parameters were taken to define the operating model
#' @author T. Carruthers
#' @keywords classes
#' @examples
#' 
#' showClass('Imp')
#' 
setClass("Imp", representation(Name = "character", TACSD = "numeric", TACFrac = "numeric", 
                               ESD = "numeric", EFrac = "numeric",
                               SizeLimSD = "numeric", SizeLimFrac="numeric",
                               DiscMort = "numeric", 
                               Source = "character"))

# initialize Imp
setMethod("initialize", "Imp", function(.Object, file = NA) {
  
  .Object@Name <- "Perfect implementation"
  .Object@TACSD <- c(0,0)
  .Object@TACFrac <- c(1,1)
  .Object@ESD <- c(0,0)
  .Object@EFrac <-c(1,1)
  .Object@SizeLimSD <- c(0,0)
  .Object@SizeLimFrac<-c(1,1)
  .Object@DiscMort <- c(0,0)
  .Object@Source <-"DLMtool generated"
  
  if (!is.na(file)) {
    if (file.exists(file)) {
      Ncol <- max(unlist(lapply(strsplit(readLines(file), ","), length)))
      dat <- read.csv(file, header = F, colClasses = "character", 
                      col.names = paste0("V", 1:Ncol))  # read 1st sheet
      dname <- dat[, 1]
      dat <- dat[, 2:ncol(dat)]
      
      .Object@Name <- dat[match("Name", dname), 1]
      .Object@TACSD <- as.numeric(dat[match("TACSD", dname), 1:2])
      .Object@TACFrac <- as.numeric(dat[match("TACFrac", dname), 1:2])
      .Object@ESD <- as.numeric(dat[match("ESD", dname), 1:2])
      .Object@EFrac <- as.numeric(dat[match("EFrac", dname), 1:2])
      .Object@SizeLimSD <- as.numeric(dat[match("SizeLimSD", dname), 1:2])
      .Object@SizeLimFrac <- as.numeric(dat[match("SizeLimFrac", dname), 1:2])
      .Object@DiscMort <- as.numeric(dat[match("DiscMort", dname), 1:2])
      .Object@Source <- dat[match("Source", dname), 1]
      
    } else {
      
      message("File doesn't exist")
      
    }
    
    
  }
  .Object
  
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
