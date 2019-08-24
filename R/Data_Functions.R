# Functions that operate on Data Object


#' Initialize an empty Data workbook or CSV
#'
#' @param name Name of the data file. Default is Data.xlsx
#' @param ext Optional file extension. 'xlsx' (default) or 'csv'
#' @param overwrite Logical. Overwrite existing files?
#' @param dir Optional directory path to create the Data file. Default is `getwd()``
#'
#' @return Nothing. Creates a data file in the working directory.
#' @export
#' 
#' @author A. Hordyk
#' @examples
#' \dontrun{
#' DataInit("MyData")
#' }
DataInit <- function(name="Data", ext=c("xlsx", "csv"), overwrite=FALSE, dir=NULL) {
  ext <- match.arg(ext)
  if (is.null(dir)) dir <- getwd()
  name <- paste(name, ext, sep=".")
  # Copy xlsx file over to working directory 
 
  message("Creating ", name, " in ", dir)
  
  if (ext == "xlsx") {
    path <- system.file("Data.xlsx", package = "DLMtool")
    pathout <- gsub("Data.xlsx", name, path)
    pathout <- gsub(dirname(pathout), dir, pathout)
  } else {
    path <- system.file("Data.csv", package = "DLMtool")
    pathout <- gsub("Data.csv", name, path)
    pathout <- gsub(dirname(pathout), dir, pathout) 
  }

  # Check if file exists 
  exist <- file.exists(pathout)
  if (exist & !overwrite) stop(name, " already exists in ", dir, ". Use 'overwrite=TRUE' to overwrite", call.=FALSE)
  copy <- file.copy(path, pathout, overwrite = overwrite)
  if (!copy) stop("Excel file not copied from ", path)

}



#' Import a Data object from Excel file
#'
#' @param name Name of the data file, with or without file extension. 
#' Include full file path if not in working directory
#' @param dec the character used in the file for decimal points.
#' @param sheet Sheet number if importing Data from XL file
#' @param silent Logical. Hide messages?
#' @return An object of class 'Data'
#' @export
#' @author A. Hordyk
#' @examples 
#' \dontrun{
#' MyData <- XL2Data("MyData.xlsx")
#' }
XL2Data <- function(name, dec=c(".", ","), sheet=1, silent=TRUE) {
  dec <- match.arg(dec)
  dir <- dirname(name)
  if (dir ==".") {
    dir <- NULL
  } else {
    name <- basename(name)
  }
  if (is.null(dir)) dir <- getwd()
  if (nchar(tools::file_ext(name)) == 0) {
    xl.fname1 <- paste0(name, ".xlsx")
    xl.fname2 <- paste0(name, ".csv")
    fls <- file.exists(c(file.path(dir, xl.fname1), file.path(dir,xl.fname2)))
    if (sum(fls) == 0) stop(xl.fname1, " or ", xl.fname2, " not found in ", dir)
    if (sum(fls) > 1) stop(name, " found with multiple extensions. Specify file extension.", call.=FALSE)
    name <- c(xl.fname1, xl.fname2)[fls]
  }
  if (tools::file_ext(name) == "csv") {
    Ncol <- max(unlist(lapply(strsplit(readLines(file.path(dir,name)), ","), length)))
    col.names <- paste0("V", 1:Ncol)
    if (dec == ".") {
      datasheet <- read.csv(file.path(dir,name), header = F, 
                            colClasses = "character", col.names=col.names,
                            stringsAsFactors = FALSE)  
    } else {
      datasheet <- read.csv2(file.path(dir,name), header = F, 
                             colClasses = "character", col.names=col.names,
                             stringsAsFactors = FALSE)  
    }
  } else if(tools::file_ext(name) %in% c("xls", "xlsx")) {
    datasheet <- readxl::read_excel(file.path(dir,name), sheet = sheet, 
                                    col_names = TRUE, .name_repair = "minimal")
  } else {
    stop("File extension must be .csv, .xls, or .xlsx")
  }
  colnames(datasheet) <- c("Name", "Data", colnames(datasheet)[3:length(colnames(datasheet))])
  
  # check names in Column 1 
  input <- file.path(system.file(package = 'DLMtool'), "Data.csv")
  valnames <- read.csv(input, header=FALSE, stringsAsFactors = FALSE)
  valnames <- valnames[,1]
  valnames <- c(valnames, "LHYear", "MPrec", 'MPeff')
  InNames <- datasheet$Name
  inval <- InNames[!InNames %in% valnames]
  inval <- inval[!grepl("CAL", inval)]
  inval <- inval[!grepl("CAA", inval)]
  inval <- inval[!grepl("Index", inval)]
  inval <- inval[!grepl("CV Index", inval)]
  inval <- inval[!grepl("Index", inval)]
  inval <- inval[!grepl("Modal length", inval)]
  inval <- inval[!grepl("Mean length Lc", inval)]
  inval <- inval[!grepl("Current spawning stock abundance", inval)]
  inval <- inval[!grepl("Cref", inval)]
  inval <- inval[!grepl("Iref", inval)]
  inval <- inval[!grepl("Bref", inval)]
  inval <- inval[!grepl("CV Cref", inval)]
  inval <- inval[!grepl("CV Iref", inval)]
  inval <- inval[!grepl("CV Bref", inval)]
  inval <- inval[!grepl("CV Rec", inval)]
  inval <- inval[!grepl("Recruitment", inval)]
  inval <- inval[!grepl("LenCV", inval)]
  
  inval <- inval[!is.na(inval)]
  if (length(inval)>0)
    warning("These rows in the Data file are not valid names and were not imported:\n", paste(inval, collapse=', '))
  
  datasheet$Name[datasheet$Name == "MPrec"] <- 'Previous TAC'
  datasheet$Name[datasheet$Name == "MPeff"] <- 'Previous TAE'
  
  datasheet$Name[datasheet$Name == "Modal length"] <- 'Modal length (Lc)'
  datasheet$Name[datasheet$Name == "Mean length Lc"] <- 'Mean length above Lc'
  datasheet$Name[datasheet$Name == "Current spawning stock abundance"] <- 'Current spawning abundance'
  datasheet$Name[datasheet$Name == "Cref"] <- 'Catch Reference'
  datasheet$Name[datasheet$Name == "Iref"] <- 'Index Reference'
  datasheet$Name[datasheet$Name == "Bref"] <- 'Biomass Reference'
  datasheet$Name[datasheet$Name == "CV Cref"] <- 'CV Catch Reference'
  datasheet$Name[datasheet$Name == "CV Iref"] <- 'CV Index Reference'
  datasheet$Name[datasheet$Name == "CV Bref"] <- 'CV Biomass Reference'
  datasheet$Name[datasheet$Name == "Recruitment"] <- 'Recruitment Index'
  datasheet$Name[datasheet$Name == "LenCV"] <- 'CV of length-at-age'
  
  
  Data <- new("Data")
  
  # ---- Main ----
  Data@Name <- datasheet$Data[which(datasheet$Name=="Name")]
  Data@Common_Name <- datasheet$Data[which(datasheet$Name=="Common Name")]
  Data@Species <- datasheet$Data[which(datasheet$Name=="Species")]
  Data@Region <- datasheet$Data[which(datasheet$Name=="Region")]
  tryLHyear <- datasheet$Data[which(datasheet$Name=="LHYear")]
  if (length(tryLHyear)<1) {
    tryLHyear <- datasheet$Data[which(datasheet$Name=="Last Historical Year")]
  }
  tryLHyear <- as.numeric(tryLHyear)
  if (length(tryLHyear)<1 | is.na(tryLHyear))
    stop("Last Historical Year must be specified (single numeric value)")
  
  Data@LHYear <- tryLHyear
  Data@MPrec <- datasheet$Data[which(datasheet$Name=="Previous TAC")] %>% as.numeric()
  Data@Units <- datasheet$Data[which(datasheet$Name=="Units")]
  Data@MPeff <- datasheet$Data[which(datasheet$Name=="Previous TAE")]  %>% as.numeric()
  Data@nareas <- datasheet$Data[which(datasheet$Name=="nareas")] %>% as.numeric()
  
  
  # ---- Biology ----
  Data@MaxAge <- datasheet$Data[which(datasheet$Name=="Maximum age")] %>% as.numeric()
  Data@Mort <- datasheet$Data[which(datasheet$Name=="M")] %>% as.numeric()
  Data@CV_Mort  <- datasheet$Data[which(datasheet$Name=="CV M")] %>% as.numeric()
  Data@vbLinf <- datasheet$Data[which(datasheet$Name=="Von Bertalanffy Linf parameter")] %>% as.numeric()
  Data@CV_vbLinf  <- datasheet$Data[which(datasheet$Name=="CV von B. Linf parameter")] %>% as.numeric()
  Data@vbK <- datasheet$Data[which(datasheet$Name=="Von Bertalanffy K parameter")] %>% as.numeric()
  Data@CV_vbK <- datasheet$Data[which(datasheet$Name=="CV von B. K parameter")] %>% as.numeric()
  Data@vbt0 <- datasheet$Data[which(datasheet$Name=="Von Bertalanffy t0 parameter")] %>% as.numeric() 
  Data@CV_vbt0 <- datasheet$Data[which(datasheet$Name=="CV von B. t0 parameter")] %>% as.numeric() 
  Data@wla <- datasheet$Data[which(datasheet$Name=="Length-weight parameter a")] %>% as.numeric()
  Data@CV_wla <- datasheet$Data[which(datasheet$Name=="CV Length-weight parameter a")] %>% as.numeric()
  Data@wlb <- datasheet$Data[which(datasheet$Name=="Length-weight parameter b")] %>% as.numeric()
  Data@CV_wlb <- datasheet$Data[which(datasheet$Name=="CV Length-weight parameter b")] %>% as.numeric()
  Data@steep <- datasheet$Data[which(datasheet$Name=="Steepness")] %>% as.numeric()
  Data@CV_steep <- datasheet$Data[which(datasheet$Name=="CV Steepness")] %>% as.numeric()
  Data@sigmaR <- datasheet$Data[which(datasheet$Name=="sigmaR")] %>% as.numeric()
  Data@CV_sigmaR <- datasheet$Data[which(datasheet$Name=="CV sigmaR")] %>% as.numeric()
  Data@L50  <- datasheet$Data[which(datasheet$Name=="Length at 50% maturity")] %>% as.numeric() 
  Data@CV_L50 <- datasheet$Data[which(datasheet$Name=="CV Length at 50% maturity")] %>% as.numeric() 
  Data@L95 <- datasheet$Data[which(datasheet$Name=="Length at 95% maturity")] %>% as.numeric() 
  Data@LenCV <- datasheet$Data[which(datasheet$Name=="CV of length-at-age")] %>% as.numeric()
  
  
  # ---- Selectivity ----
  Data@LFC <- datasheet$Data[which(datasheet$Name=="Length at first capture")] %>% as.numeric() 
  Data@CV_LFC <- datasheet$Data[which(datasheet$Name=="CV Length at first capture")] %>% as.numeric() 
  Data@LFS <- datasheet$Data[which(datasheet$Name=="Length at full selection")] %>% as.numeric() 
  Data@CV_LFS <- datasheet$Data[which(datasheet$Name=="CV Length at full selection")] %>% as.numeric() 
  Data@Vmaxlen <-datasheet$Data[which(datasheet$Name=="Vulnerability at asymptotic length")] %>% as.numeric()  
  
  # ---- Time-Series ----
  row <- which(datasheet$Name == "Year")
  Year <- datasheet[row,2:ncol(datasheet)] %>% as.numeric()
  if (!Data@LHYear %in% Year) 
    stop("`Year` must include Last Historical Year", call. = FALSE)
  Year <- Year[!is.na(Year)]
  if (!all(seq(Year[1], Year[length(Year)], 1) == Year))
    stop("`Year` must be sequential and include all years")
  Data@Year <- Year   
  Nyears <- length(Year)
  # Catch time-series
  Data@Cat <- datasheet[which(datasheet$Name=="Catch"), 2:(Nyears+1)] %>% 
    as.numeric() %>% matrix(., nrow=1)
  CV_Cat <- datasheet[which(datasheet$Name=="CV Catch"), 2:(Nyears+1)]
  if (!is.na(CV_Cat[1]) & all(is.na(CV_Cat[2:length(CV_Cat)])))
    CV_Cat <- rep(CV_Cat[1], Nyears)
  Data@CV_Cat <- CV_Cat %>% as.numeric %>% matrix(., nrow=1)
  
  # Effort time-series
  Data@Effort <- suppressWarnings(datasheet[which(datasheet$Name=="Effort"), 2:(Nyears+1)] %>% 
                                    as.numeric() %>% matrix(., nrow=1))
  CV_Effort <- datasheet[which(datasheet$Name=="CV Effort"), 2:(Nyears+1)]
  if (nrow(CV_Effort)>0) {
    if (!is.na(CV_Effort[1]) & all(is.na(CV_Effort[2:length(CV_Effort)])))
      CV_Effort <- rep(CV_Effort[1], Nyears)
    Data@CV_Effort <- CV_Effort %>% as.numeric %>% matrix(., nrow=1)
  }
  
  # Vulnerable abundance index - fishery dependant
  Data@Ind <- suppressWarnings(datasheet[which(datasheet$Name=="Abundance index"), 2:(Nyears+1)] %>%
                                 as.numeric() %>% matrix(., nrow=1) )
  CV_Ind <- datasheet[which(datasheet$Name=="CV Abundance index"), 2:(Nyears+1)]
  if (!is.na(CV_Ind[1]) & all(is.na(CV_Ind[2:length(CV_Ind)])))
    CV_Ind <- rep(CV_Ind[1], Nyears)
  Data@CV_Ind <- CV_Ind %>% as.numeric %>% matrix(., nrow=1)
  
  # # Spawning abundance index - subject to hyper-stability beta
  # Data@SpInd <- suppressWarnings(datasheet[which(datasheet$Name=="SpAbun index"), 2:(Nyears+1)] %>%
  #                                  as.numeric() %>% matrix(., nrow=1) )
  # CV_SpInd <- datasheet[which(datasheet$Name=="CV SpAbun index"), 2:(Nyears+1)]
  # if (nrow(CV_SpInd)>0) {
  #   if (!is.na(CV_SpInd[1]) & all(is.na(CV_SpInd[2:length(CV_SpInd)])))
  #     CV_SpInd <- rep(CV_SpInd[1], Nyears)
  #   Data@CV_SpInd <- CV_SpInd %>% as.numeric %>% matrix(., nrow=1)
  # }
  
  # Additional indices
  ind <- grepl("Index", datasheet$Name)
  index_text <- datasheet$Name[ind]
  index_text <- index_text[! index_text %in% c("Index Reference", "CV Index Reference")]
  n_indices <- sum(sub(" .*$", "", index_text) == "Index")
  n_cv <- sum(sub(" .*$", "", index_text) == "CV")
  n_vuln <- sum(sub(" .*$", "", index_text) == "Vuln")
  if (n_cv != n_indices) stop("CV missing for some or all additional indices", call. = FALSE)
  if (n_vuln != n_indices) stop("Vulnerability-at-age schedule missing for some or all additional indices", call. = FALSE)
  
  if (!all(is.na(datasheet[which(datasheet$Name == "Index 1"),2:(Nyears+1)]))) {
    indexexist <- TRUE
  } else {
    indexexist <- FALSE
  }
  if (indexexist) {
    if (is.na(Data@MaxAge))
      stop("Require `Maximum age` to use Additional Indices", call. = FALSE)
    Data@AddInd <- Data@CV_AddInd <- array(NA, dim=c(1, n_indices, Nyears))
    Data@AddIndV <- array(NA, dim=c(1, n_indices, Data@MaxAge))
  } else {
    Data@AddInd <- Data@CV_AddInd <- array(NA, dim=c(1,1,1))
    Data@AddIndV <- array(NA, dim=c(1,1,1))
  }
  
  # Loop over additional indices 
  if (indexexist) {
    for (x in 1:n_indices) {
      ind <- which(datasheet$Name == paste("Index", x))
      Data@AddInd[1,x,] <- datasheet[ind, 2:(Nyears+1)] %>% as.numeric()
      ind <- which(datasheet$Name == paste("CV Index", x))
      Data@CV_AddInd[1,x,] <- datasheet[ind, 2:(Nyears+1)] %>% as.numeric()
      ind <- which(datasheet$Name == paste("Vuln Index", x))
      Data@AddIndV[1,x,] <- datasheet[ind, 2:(Data@MaxAge+1)] %>% as.numeric() 
      if (any(is.na(Data@AddIndV[1,x,])))
        warning("Vuln Index must be length `Maximum age` and contain only numeric values (no NA)")
    }
  }
  
  # Standarise indices 
  slot(Data, 'Ind') <- slot(Data, 'Ind')/mean(slot(Data,'Ind'), na.rm=TRUE)
  
  if (indexexist) {
    dims <- dim(Data@AddInd)
    nind <- dims[2]
    for (i in 1:nind) {
     Data@AddInd[1,i,] <-   Data@AddInd[1,i,]/mean(Data@AddInd[1,i,], na.rm=TRUE)
    }
  } 
 
  # Recruitment index 
  ind <- grepl('Recruitment', datasheet$Name) %>% which()
  Data@Rec <- suppressWarnings(datasheet[ind[1], 2:(Nyears+1)] %>% as.numeric() %>%
                                 matrix(., nrow=1))
  Data@CV_Rec <- suppressWarnings(datasheet[ind[2], 2:(Nyears+1)] %>% as.numeric() %>%
                                    matrix(., nrow=1))
  
  # Mean length
  ind <- which(datasheet$Name == "Mean length")
  Data@ML <- suppressWarnings(datasheet[ind[1], 2:(Nyears+1)] %>% as.numeric() %>%
                                matrix(., nrow=1))
  ind <- which(datasheet$Name == "Modal length (Lc)")
  Data@Lc <- suppressWarnings(datasheet[ind[1], 2:(Nyears+1)] %>% as.numeric() %>%
                                matrix(., nrow=1))
  ind <- which(datasheet$Name == "Mean length above Lc")
  Data@Lbar <- suppressWarnings(datasheet[ind[1], 2:(Nyears+1)] %>% as.numeric() %>%
                                  matrix(., nrow=1))
  
  # ---- Catch-at-Age ----
  ind <- which(datasheet$Name == "Vuln CAA")
  if (length(ind)>0 && (!is.na(datasheet[ind,2]))) {
    VulnCAA <- datasheet[ind, 2:(Data@MaxAge+1)]
    nonNA <- as.vector(!is.na(VulnCAA[1,]))
    if (!all(is.na(VulnCAA)) & !all(nchar(VulnCAA[nonNA])==0)) {
      if (any(is.na(VulnCAA))) 
        warning("Vuln CAA must be length `Maximum age` and contain only numeric values (no NA)")
    }
    Data@Vuln_CAA <- suppressWarnings(matrix(as.numeric(VulnCAA), nrow=1))
  }
  
  ind <- which(grepl('CAA', datasheet$Name))
  ind2 <-  which(datasheet$Name == "Vuln CAA")
  if (length(ind2)>0) 
    ind <- ind[!ind ==ind2]
  CAAexists <- datasheet$Name[ind] == "CAA"
  if (length(CAAexists) < 1) CAAexists <- FALSE
  if (length(ind) <=1 || !CAAexists) {
    CAA_Yrs <- numeric(0)
  } else {
    CAA_Yrs <- sapply(strsplit(datasheet$Name[ind], " "), function(x) unlist(strsplit(x[2], " ")))
    if (any(is.na(CAA_Yrs))) 
      stop("Missing years for CAA data. Must be of the form `CAA YEAR` e.g., `CAA 2019`", call. = FALSE)
    
    if(!all(CAA_Yrs %in% Data@Year)) stop("All CAA Years must be included in `Year`")
  }
 
  
  if (length(CAA_Yrs)>0) {
    Data@CAA <- array(NA, dim=c(1, Nyears, Data@MaxAge))
    CAAdat <- datasheet[ind, 2:(Data@MaxAge+1)] %>% as.matrix() %>% as.numeric() %>% matrix(., nrow=length(ind), ncol=Data@MaxAge)
    if (any(is.na(CAAdat))) 
      stop("NAs in CAA data. Is each row of length `Maximum age`?", call. = FALSE)
    yrind <- match(CAA_Yrs, Data@Year)
    Data@CAA[1, yrind,] <- CAAdat
  } else{
    Data@CAA <- array(NA, dim=c(1,1,1))
  }
  
  # ---- Catch-at-Length ----
  ind <- which(datasheet$Name == "Vuln CAL")
  if (length(ind)>0 && !is.na(datasheet$Data[ind])) {
    VulnCAL <- datasheet[ind, 2:(Data@MaxAge+1)]
    nonNA <- as.vector(!is.na(VulnCAL[1,]))
    if (!all(is.na(VulnCAL)) & !all(nchar(VulnCAL[nonNA])==0)) {
      if (any(is.na(VulnCAL))) 
        warning("Vuln CAL must be length `Maximum age` and contain only numeric values (no NA)")
    }
    Data@Vuln_CAL <- suppressWarnings(matrix(as.numeric(VulnCAL), nrow=1))
  }
  
  CAL_bins <- suppressWarnings(datasheet[which(datasheet$Name == "CAL_bins"),] %>% as.numeric())
  CAL_bins <- CAL_bins[!is.na(CAL_bins)]
  
  CAL_mids <- suppressWarnings(datasheet[which(datasheet$Name == "CAL_mids"),] %>% as.numeric())
  CAL_mids <- CAL_mids[!is.na(CAL_mids)]
  
  ind <- which(grepl('CAL', datasheet$Name) & !grepl('CAL_bins', datasheet$Name) &
                 !grepl('Vuln CAL', datasheet$Name) &
                 !grepl('CAL_mids', datasheet$Name))
  
  CALexists <- grepl("CAL", datasheet$Name[ind])
  if (length(CALexists) < 1) CALexists <- FALSE
  if (length(ind) <=1 || !any(CALexists)) {
    ind <- numeric(0)
    CAL_Yrs <- numeric(0)
  } 
  if (length(ind)>0) {
    # CAL data exists
    if (length(CAL_bins)<1 & length(CAL_mids) < 1)
      stop("Require either CAL_mids or CAL_bins", call. = FALSE)
    if (length(CAL_bins)>1  & length(CAL_mids) > 1) {
      if (!length(CAL_mids) == length(CAL_bins)-1)
        stop("CAL_mids should be of length `length(CAL_bins)-1`", call. = FALSE)
      by <- CAL_mids[2] - CAL_mids[1]
      if (all(CAL_mids != seq(CAL_bins[1]+0.5*by, by=by, length.out = length(CAL_bins)-1)))
        stop("CAL_mids should be mid-points of CAL_bins", call. = FALSE)
    } 
    if (length(CAL_bins)>1) {
      if (!all(diff(CAL_bins[2:length(CAL_bins)]) == diff(CAL_bins)[2])) 
        stop("Inconsistent bin width in CAL_bins", call. = FALSE)
    }
    if (length(CAL_mids)>1) {
      if(!all(diff(CAL_mids) == mean(diff(CAL_mids)))) 
        stop("Inconsistent bin width in CAL_mids", call. = FALSE)
    }
    if (length(CAL_bins)>1 & length(CAL_mids)<1) {
      by <- CAL_bins[2] - CAL_bins[1]
      CAL_mids <- seq(CAL_bins[1]+0.5*by, by=by, length.out = length(CAL_bins)-1)
    }
    
    CAL_Yrs <- sapply(strsplit(datasheet$Name[ind], " "), function(x) unlist(strsplit(x[2], " ")))
    
  }
  
  ind2 <- which(datasheet$Name == "Vuln CAL")
  if (length(ind2)>01)
    ind <- ind[!ind == ind2]
  
  CAL_Yrs <- sapply(strsplit(datasheet$Name[ind], " "), function(x) unlist(strsplit(x[2], " ")))
  if(!all(CAL_Yrs %in% Data@Year)) stop("All CAL Years must be included in `Year`")
  
  NMids <- length(CAL_mids)
  Data@CAL <- array(NA, dim=c(1, Nyears, NMids))
  
  CALdat <- datasheet[ind, 2:(NMids+1)] %>% as.matrix() %>% as.numeric() %>% matrix(., nrow=length(ind), ncol=NMids)
  if (any(is.na(CALdat))) 
    stop("NAs in CAL data. Is each row same length as `CAL_mids`? (or `length(CAL_bins)-1`)", call. = FALSE)
  yrind <- match(CAL_Yrs, Data@Year)
  Data@CAL_bins <- CAL_bins
  Data@CAL_mids <- CAL_mids
  Data@CAL[1, yrind,] <- CALdat
  
  # ---- Reference ----
  Data@Dep <- datasheet[which(datasheet$Name=="Current stock depletion"),2] %>% as.numeric()
  Data@CV_Dep <- datasheet[which(datasheet$Name=="CV current stock depletion"),2] %>% as.numeric()
  Data@Abun <- datasheet[which(datasheet$Name=="Current stock abundance"),2] %>% as.numeric()
  Data@CV_Abun <- datasheet[which(datasheet$Name=="CV current stock abundance"),2] %>% as.numeric()
  Data@SpAbun <- datasheet[which(datasheet$Name=="Current spawning abundance"),2] %>% as.numeric()
  Data@CV_SpAbun <- datasheet[which(datasheet$Name=="CV current spawning abundance"),2] %>% as.numeric()
  
  Data@FMSY_M <- datasheet[which(datasheet$Name=="FMSY/M"),2] %>% as.numeric()
  Data@CV_FMSY_M <- datasheet[which(datasheet$Name=="CV FMSY/M"),2] %>% as.numeric()
  Data@BMSY_B0 <- datasheet[which(datasheet$Name=="BMSY/B0"),2] %>% as.numeric()
  Data@CV_BMSY_B0 <- datasheet[which(datasheet$Name=="CV BMSY/B0"),2] %>% as.numeric()
  Data@Cref <- datasheet[which(datasheet$Name=="Catch Reference"),2] %>% as.numeric()
  Data@CV_Cref <- datasheet[which(datasheet$Name=="CV Catch Reference"),2] %>% as.numeric()
  Data@Bref <- datasheet[which(datasheet$Name=="Biomass Reference"),2] %>% as.numeric()
  Data@CV_Bref <- datasheet[which(datasheet$Name=="CV Biomass Reference"),2] %>% as.numeric()
  
  Data@Iref <- datasheet[which(datasheet$Name=="Index Reference"),2] %>% as.numeric()
  Data@CV_Iref <- datasheet[which(datasheet$Name=="CV Index Reference"),2] %>% as.numeric()
  Data@t <- datasheet[which(datasheet$Name=="Duration t"),2] %>% as.numeric()
  Data@AvC <- datasheet[which(datasheet$Name=="Average catch over time t"),2] %>% as.numeric()
  Data@CV_AvC <- datasheet[which(datasheet$Name=="CV Average catch over time t"),2] %>% as.numeric()
  Data@Dt <- datasheet[which(datasheet$Name=="Depletion over time t"),2] %>% as.numeric()
  Data@CV_Dt <- datasheet[which(datasheet$Name=="CV Depletion over time t"),2] %>% as.numeric()
  Data@Ref <- datasheet[which(datasheet$Name=="Reference OFL"),2] %>% as.numeric()
  Data@Ref_type <- datasheet[which(datasheet$Name=="Reference OFL type"),2] %>% as.character() 
  
  # Default values
  if (all(is.na(Data@CV_Cat))) Data@CV_Cat <- matrix(0.2, nrow=1, ncol=Nyears)
  if (all(is.na(Data@CV_Ind))) Data@CV_Ind <- matrix(0.2, nrow=1, ncol=Nyears)
  # if (all(is.na(Data@CV_SpInd))) Data@CV_SpInd <- matrix(0.2, nrow=1, ncol=Nyears)
  if (all(is.na(Data@CV_Effort))) Data@CV_Effort <- matrix(0.2, nrow=1, ncol=Nyears)
  if (all(is.na(Data@CV_Rec))) Data@CV_Rec <- matrix(0.2, nrow=1, ncol=Nyears)
  
  if (NAor0(Data@LenCV)) Data@LenCV <- 0.1
  if (NAor0(Data@CV_Dt)) Data@CV_Dt <- 0.25
  if (NAor0(Data@CV_AvC)) Data@CV_AvC <- 0.2
  if (NAor0(Data@CV_Mort)) Data@CV_Mort <- 0.2
  if (NAor0(Data@CV_FMSY_M)) Data@CV_FMSY_M <- 0.2
  if (NAor0(Data@CV_BMSY_B0)) Data@CV_BMSY_B0 <- 0.045
  if (NAor0(Data@CV_Cref)) Data@CV_Cref <- 0.2
  if (NAor0(Data@CV_Bref)) Data@CV_Bref <- 0.2
  if (NAor0(Data@CV_Iref)) Data@CV_Iref <- 0.2
  if (NAor0(Data@CV_Dep)) Data@CV_Dep <- 0.25
  if (NAor0(Data@CV_Abun)) Data@CV_Abun <- 0.25
  if (NAor0(Data@CV_vbK)) Data@CV_vbK <- 0.1
  if (NAor0(Data@CV_vbLinf)) Data@CV_vbLinf <- 0.1
  if (NAor0(Data@CV_vbt0)) Data@CV_vbt0 <- 0.1
  if (NAor0(Data@CV_L50))  Data@CV_L50 <- 0.1
  if (NAor0(Data@CV_LFC))  Data@CV_LFC <- 0.2
  if (NAor0(Data@CV_LFS))  Data@CV_LFS <- 0.2
  if (NAor0(Data@CV_wla))  Data@CV_wla <- 0.1
  if (NAor0(Data@CV_wlb))  Data@CV_wlb <- 0.1
  if (NAor0(Data@CV_steep)) Data@CV_steep <- 0.2
  if (NAor0(Data@nareas)) Data@nareas <- 2
  
  if (length(Data@CAA) == 0) Data@CAA <- array(NA, c(1, 1, 1))
  if (length(Data@CAL) == 0) Data@CAL <- array(NA, c(1, 1, 1))
  if (length(Data@CAL_bins) == 0) Data@CAL_bins <- 1
  if (length(Data@TAC) == 0) Data@TAC <- array(1, c(1, 1))
  # if (length(Data@TACbias) == 0) Data@TACbias <- array(1, c(1, 1))
  if (length(Data@Sense) == 0) Data@Sense <- array(1, c(1, 1))
  if (length(Data@ML) == 0)  Data@ML <- array(NA, c(1, 1))
  if (length(Data@Lbar) == 0) Data@Lbar <- array(NA, c(1, 1))
  if (length(Data@Lc) == 0) Data@Lc <- array(NA, c(1, 1))
  
  if (is.na(Data@MPeff) || length(Data@MPeff)==0) Data@MPeff <- 1
  
  Data@Log[[1]] <- paste("Created:", Sys.time())
  Data@params <- new("list")
  Data@OM <- data.frame(NA)
  Data@Obs <- data.frame(NA)
  Data@PosMPs <- NA
  Data@MPs <- NA
  
  Data  
}

# XL2Data <- function(name="Data") {
#   if (class(name) != 'character') stop("file name must be provided", call.=FALSE)
#   
#   dir <- dirname(name)
#   if (dir ==".") {
#     dir <- NULL
#   } else {
#     name <- basename(name)
#   }
#   
#   if (is.null(dir)) dir <- getwd()
#   if (nchar(tools::file_ext(name)) == 0) {
#     xl.fname1 <- paste0(name, ".xlsx")
#     xl.fname2 <- paste0(name, ".csv")
#     fls <- file.exists(c(file.path(dir, xl.fname1), file.path(dir,xl.fname2)))
#     if (sum(fls) == 0) stop(xl.fname1, " or ", xl.fname2, " not found in ", dir)
#     if (sum(fls) > 1) stop(name, " found with multiple extensions. Specify file extension.", call.=FALSE)
#     name <- c(xl.fname1, xl.fname2)[fls]
#   }
#   
#   if (!file.exists(file.path(dir, name))) stop(file.path(dir, name), " not found", call.=FALSE) 
#   
#   isCSV <- grepl('.csv', name)
#   message("Reading ", name)
#   if (isCSV) {
#     Data <- new("Data", file.path(dir,name))
#   } else {
#     sheetnames <- readxl::excel_sheets(file.path(dir,name))  # names of the sheets
#     
#     # DataXLSlot <- DLMtool:::DataXLSlot
#     NewSheetNames <- names(DataXLSlot)
#     if (all(NewSheetNames %in% sheetnames)) {
#       Data <- importnewXLData(dir,name, NewSheetNames)
#     } else {
#       datasheet <- suppressMessages(as.data.frame(readxl::read_excel(file.path(dir,name), sheet = 1, col_names = FALSE)))
#       if (datasheet[1,1]== "Slot") 
#         datasheet <- suppressMessages(as.data.frame(readxl::read_excel(file.path(dir,name), sheet = 1, col_names = FALSE, skip=1)))
#       
#       if (all(dim(datasheet) == 0)) stop("Nothing found in first sheet", call.=FALSE)
#       tmpfile <- tempfile(fileext=".csv")
#       writeCSV2(inobj = datasheet, tmpfile, objtype = "Data")
#       
#       if (ncol(datasheet)<2) {
#         unlink(tmpfile)
#         stop("No parameter values found in first worksheet ", call.=FALSE)
#       } else {
#         Data <- new("Data", tmpfile)
#         unlink(tmpfile)
#       }
#     }
#     return(Data)
#   }
# }





# importnewXLData_old <- function(dir,name, NewSheetNames) {
#   Data <- new("Data", silent=TRUE)
#   BlankDat <-new("Data", silent=TRUE)
#   ignoreSheet <- NULL
#   
#   # sh <- NewSheetNames[1]
#   for (sh in NewSheetNames) { # import from individual worksheets
#     datasheet <- suppressMessages(as.data.frame(readxl::read_excel(file.path(dir,name), 
#                                                   sheet = sh, col_names = FALSE)))
#     if (dim(datasheet)[2] <= 1) {
#       message('No data found in sheet: ', sh)
#     } else {
#       message("Importing from: ", sh)
#       dname <- datasheet[, 1]
#       dat  <- datasheet[, 2:ncol(datasheet), drop=FALSE]
#       df <- data.frame(XLRow=DataXLSlot[[sh]]$XLRow, 
#                        Slot=DataXLSlot[[sh]]$Slot, 
#                        Class=DataXLSlot[[sh]]$Class, 
#                        Ignore=DataXLSlot[[sh]]$Ignore,
#                        stringsAsFactors = FALSE)
#       df$Ignore[is.na(df$Slot)] <- TRUE
#       df$Ignore[df$XLRow=="Fleet Type"] <- TRUE # temporary until new build
#       df$Ignore <- as.logical(df$Ignore)
#       df <- df[!df$Ignore,]
#       
#       # if (sh %in% c("Main", "Biology", "Reference")) {
#       if (sh %in% c("Main", "Biology", "Selectivity", "Reference")) {
#         for (sl in 1:nrow(df)) {
#           temp <- dat[match(df$XLRow[sl], dname),1]
#           if (substr(df$Class[sl],start=1, stop=1) == "c") 
#             slot(Data, df$Slot[sl]) <- temp
#           if (substr(df$Class[sl],start=1, stop=1) == "n") 
#             slot(Data, df$Slot[sl]) <- as.numeric(temp)
#         }
#       } else if (sh == "Time-Series") {
#         YearInd <- match("Year", dname)
#         Year <- dat[YearInd,]
#         Year <- Year[!is.na(Year)]
#         Data@Year <- as.numeric(Year)
#         if (!is.finite(Data@LHYear)) 
#           stop("'Current Year' in 'Main' sheet is missing", call.=FALSE)
#         if (max(Data@Year) != Data@LHYear) 
#           stop("Last year should be equal to 'Current Year' in 'Main' sheet", call.=FALSE)
#         ncol_ts <- length(Year)
#         ncol_cv <- 1
#         for (sl in 1:nrow(df)) {
#           ncol <- ifelse(grepl("CV_", df$Slot[sl]), ncol_cv, ncol_ts)
#           temp <- as.numeric(dat[match(df$XLRow[sl], dname),1:ncol])
#           
#           if (substr(df$Class[sl],start=1, stop=1) == "n")
#             slot(Data, df$Slot[sl]) <- temp
#           if (substr(df$Class[sl],start=1, stop=1) == "m")
#             slot(Data, df$Slot[sl]) <- matrix(temp, nrow=1)
#           if(all(is.na(slot(Data, df$Slot[sl]))))  
#             slot(Data, df$Slot[sl]) <- slot(BlankDat, df$Slot[sl]) 
#         }
#         
#       } else if (sh == "CAA") {
#         if (!is.finite(Data@MaxAge)) 
#           stop("'Maximum age' in 'Biology' sheet is missing", call.=FALSE)
#         if (any(!is.finite(Data@Year)))
#           stop("'Year' in 'Time-Series' sheet is missing", call.=FALSE)
#         CAAMat <- array(NA, dim=c(1,length(Data@Year),Data@MaxAge))
#         Year <- Data@Year
#         
#         CAAYr <- as.numeric(dname[dname !="Year"])
#         if (length(CAAYr)<1) {
#           message("No catch-at-age data found")
#         } else {
#           YrInd <- match(CAAYr, Year) # match years 
#           if (max(CAAYr) > max(Year))
#             stop("More years in CAA than Time-Series Year", call.=FALSE)
#           CAAdat <- data.matrix(dat[2:nrow(dat),])
#           
#           if (ncol(CAAdat)> Data@MaxAge) 
#             stop("Number of age-classes in CAA data > MaxAge", call.=FALSE)
#           if (ncol(CAAdat) < Data@MaxAge) {
#             message("Number of age-classes in CAA data < MaxAge. Filling with 0s")
#             zeromat <- matrix(0, nrow=nrow(CAAdat), ncol=Data@MaxAge-ncol(CAAdat))
#             CAAdat <- cbind(CAAdat, zeromat)
#           }
#           CAAMat[1, YrInd, 1:ncol(CAAdat)] <- CAAdat
#           Data@CAA <- CAAMat
#         }
#       } else if (sh == "CAL") {
#         CAL_bins <- as.numeric(dat[1,])
#         if (!all(diff(CAL_bins) == diff(CAL_bins))) 
#           stop('Length bins do not have equal intervals' , call.=FALSE)
#         if (any(!is.finite(Data@Year)))
#           stop("'Year' in 'Time-Series' sheet is missing", call.=FALSE)
#         CALMat <- array(NA, dim=c(1,length(Data@Year), length(CAL_bins)-1))
#         Year <- Data@Year
#         CALYr <- as.numeric(dname[dname !="Year" & dname !="CAL_bins"])
#         if (length(CALYr)<1) {
#           message("No catch-at-length data found")
#         } else {
#           YrInd <- match(CALYr, Year) # match years
#           if (max(CALYr) > max(Year))
#             stop("More years in CAL than Time-Series Year", call.=FALSE)
#           ncol <- ncol(dat)
#           CALdat <- data.matrix(dat[3:nrow(dat),1:ncol])
#           if (!all(is.na(CALdat[,ncol]))) {
#             stop("Number of Catch-at-Length bins (CAL_bins) should \nbe 1 greater than number columns of CAL data", call.=FALSE)
#           }
#           CALdat <- CALdat[,1:(ncol-1)]
#           CALMat[1, YrInd, ] <- CALdat
#         }
#         Data@CAL_bins <- CAL_bins
#         Data@CAL <- CALMat
#       }
#     }
#     
#   } # end sheet loop
#   Data
# }




#' Run a Management Procedure
#'
#' @param Data A DLMtool Data object
#' @param MPs The name of the MP to run (or a vector or names)
#' @param reps Number of repititions
#' @param perc Percentile to summarize reps (default is median)
#' @param chkMPs Logical. Should the MPs be checked before attempting to run them?
#' @param silent Logical. Should messages by suppressed?
#'
#' @export
#' @examples
#' Data_TAc <- runMP(DLMtool::Cobia)
#' @return invisibly returns the Data object
#'
runMP <- function(Data, MPs = NA, reps = 100, perc=0.5, chkMPs=TRUE, silent=FALSE) {
  if (class(MPs) != 'character' && !all(is.na(MPs))) stop('MPs must be character string', call.=FALSE)
  if (class(Data) != 'Data') stop("Data must be class 'Data'", call.=FALSE)
  if (all(is.na(MPs))) {
    MPs <- avail("MP")
    if (!silent) message("running all available MPs")
  }
  if (chkMPs) {
    cans <- Can(Data, MPs=MPs)
    MPs <- MPs[MPs %in% cans]
  }
  if (length(MPs) <1) stop("No MPs possible")
  
  MPrecs <- applyMP(Data, MPs, reps, nsims=1, silent=silent)
  
  names <- c("TAC", "Effort", "LR5", "LFR", "HS", "Rmaxlen",
             "L5", "LFS", 'Vmaxlen', 'Spatial')
  mat <- matrix(0, nrow=length(MPs), ncol=length(names)+Data@nareas-1)
  for (x in seq_along(names)) {
    temp <- lapply(MPrecs[[1]], '[[', names[x])
    if (names[x]!="Spatial") {
      mat[,x] <- unlist(lapply(temp, quantile, probs=perc, na.rm=TRUE))
    } else {
      mat[,x:ncol(mat)] <- t(matrix(unlist(temp), nrow=Data@nareas, ncol=length(MPs)))
    }
  }
  rownames(mat) <- MPs 
  names[length(names)] <- "Area 1"
  if (Data@nareas > 1) names <- c(names, paste('Area', 2:Data@nareas))
  colnames(mat) <- names
  
  if (nrow(mat) > 1) {
    allNA <- colSums(apply(mat, 2, is.na)) == length(MPs)
    matout <- data.frame(round(mat[,!allNA], 2), stringsAsFactors = FALSE)
    names(matout) <- names[!allNA]
    if (!silent) print(as.matrix(matout), na.print="")
  }
  if (nrow(mat) == 1) {
    mat <- data.frame(mat)
    matout <- mat[!is.na(mat)]
    matout <- matrix(matout, nrow=nrow(mat))
    colnames(matout) <- names[!is.na(mat)]
    rownames(matout) <- MPs
    if (!silent) print(as.matrix(round(matout,2)), na.print="")
  }
  
  invisible(MPrecs[[2]])
}

#' Apply Management Procedures to an object of class Data
#'
#' @param Data An object of class Data
#' @param MPs Name(s) of the MPs to run
#' @param reps Number of samples
#' @param nsims Optional. Number of simulations. 
#' @param silent Logical. Should messages be suppressed?
#'
#' @return A list with the first element a list of management recommendations,
#' and the second the updated Data object
#' @export
#'
applyMP <- function(Data, MPs = NA, reps = 100, nsims=NA, silent=FALSE) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  Dataout <- Data
  if (is.na(nsims)) nsims <- nrow(Data@Cat)
  nMPs <- length(MPs)
  
  if (.hasSlot(Data, "nareas")) {
    nareas <- Data@nareas   
  } else {
    nareas <- 2 
  }
  returnList <- list() # a list nMPs long containing MPs recommendations
  recList <- list() # a list containing nsim recommendations from a single MP 
  TACout <- array(NA, dim=c(nMPs, reps, nsims))
  # if (!sfIsRunning() | (nMPs < 8 & nsims < 8)) {
  for (mp in 1:nMPs) {
    temp <- lapply(1:nsims, MPs[mp], Data = Data, reps = reps)  
    slots <- slotNames(temp[[1]])
    for (X in slots) { # sequence along recommendation slots 
      if (X == "Misc") { # convert to a list nsim by nareas
        rec <- lapply(temp, slot, name=X)
      } else {
        rec <- do.call("cbind", lapply(temp, slot, name=X)) # unlist(lapply(temp, slot, name=X))
      }
      if (X == "Spatial") { # convert to a matrix nsim by nareas
        rec <- matrix(rec, nareas, nsims, byrow=FALSE)   
      }
      recList[[X]] <- rec
      for (x in 1:nsims) Dataout@Misc[[x]] <- recList$Misc[[x]]
      recList$Misc <- NULL
    }
    if (length(recList$TAC)>0)  TACout[mp,,] <- recList$TAC 
    returnList[[mp]] <- recList
    if (!silent && any(apply(is.na(recList$TAC), 2, sum) > rep(0.5 * reps, nsims)))
      message("Method ", MPs[mp], " produced greater than 50% NA values")
  }

  Dataout@TAC <- TACout
  Dataout@MPs <- MPs
  
  list(returnList, Dataout)
}

# applyMP2 <- function(Data, MPs = NA, reps = 100, nsims=NA, silent=FALSE) {
#   if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
#   Data <- updateMSE(Data)
#   if (is.na(nsims)) nsims <- length(Data@Mort)
#   nMPs <- length(MPs)
#   
#   if (.hasSlot(Data, "nareas")) {
#     nareas <- Data@nareas   
#   } else {
#     nareas <- 2 
#   }
#   returnList <- list() # a list nMPs long containing MPs recommendations
#   recList <- list() # a list containing nsim recommendations from a single MP 
#   TACout <- array(NA, dim=c(nMPs, reps, nsims))
#   # if (!sfIsRunning() | (nMPs < 8 & nsims < 8)) {
#     for (mp in 1:nMPs) {
#       temp <- sapply(1:nsims, MPs[mp], Data = Data, reps = reps)  
#       slots <- slotNames(temp[[1]])
#       for (X in slots) { # sequence along recommendation slots 
#         if (X == "Misc") { # convert to a list nsim by nareas
#           rec <- lapply(temp, slot, name=X)
#         } else {
#           rec <- do.call("cbind", lapply(temp, slot, name=X)) # unlist(lapply(temp, slot, name=X))
#         }
#         if (X == "Spatial") { # convert to a matrix nsim by nareas
#           rec <- matrix(rec, nareas, nsims, byrow=FALSE)   
#         }
#         recList[[X]] <- rec
#         for (x in 1:nsims) Data@Misc[[x]] <- recList$Misc[[x]]
#         recList$Misc <- NULL
#       }
#       if (length(recList$TAC)>0)  TACout[mp,,] <- recList$TAC 
#       returnList[[mp]] <- recList
#       if (!silent && sum(is.na(recList$TAC)) > 0.5 * reps)
#         message("Method ", MPs[mp], " produced greater than 50% NA values")
#     }
#   # } else {
#   #   for (mp in 1:nMPs) {
#   #     temp <- sfSapply(1:nsims, MPs[mp], Data = Data, reps = reps)  
#   #     slots <- slotNames(temp[[1]])
#   #     for (X in slots) { # sequence along recommendation slots 
#   #       if (X == "Misc") { # convert to a list nsim by nareas
#   #         rec <- lapply(temp, slot, name=X)
#   #       } else {
#   #         rec <- do.call("cbind", lapply(temp, slot, name=X)) # unlist(lapply(temp, slot, name=X))
#   #       }
#   #       if (X == "Spatial") { # convert to a matrix nsim by nareas
#   #         rec <- matrix(rec, nareas, nsims, byrow=FALSE)  
#   #       }
#   #       recList[[X]] <- rec
#   #       for (x in 1:nsims) Data@Misc[[x]] <- recList$Misc[[x]]
#   #       recList$Misc <- NULL
#   #     }
#   #     if (length(recList$TAC)>0) TACout[mp,,] <- recList$TAC
#   #     returnList[[mp]] <- recList
#   #     
#   #     if (!silent && sum(is.na(recList$TAC)) > 0.5 * reps)
#   #       message("Method ", MPs[mp], " produced greater than 50% NA values")
#   #   }
#   #   
#   # }
#   
#   Data@TAC <- TACout
#   Data@MPs <- MPs
#   
#   list(returnList, Data)
# }


#' Identify management procedures (MPs) based on data availability
#' 
#' Diagnostic tools that look up the slot requirements of each MP and
#' compares to the data available in the Data object.
#'  
#' @param Data A data-limited methods data object (class Data)
#' @param timelimit The maximum time (seconds) taken for an MP to undertake
#' 5 reps (this filters out methods that are too slow)
#' @param MPs Optional list of MP names
#' @param dev Logical. Run in development mode?
#' @seealso \link{avail} \linkS4class{Data}
#' @examples 
#' CanMPs <- Can(DLMtool::Cobia)
#' CantMPs <- Cant(DLMtool::Cobia)
#' Needs <- Needed(DLMtool::Cobia)
#' @describeIn Can Identifies MPs that have the correct data, do not produce errors,
#' and run within the time limit.
#' @export 
Can <- function(Data, timelimit = 1, MPs=NA, dev=FALSE) {
  DLMdiag(Data, "available",  timelimit = timelimit, funcs1=MPs, dev=dev)
}


#' @describeIn Can Identifies MPs that don't have sufficient data, lead to errors, or don't run in
#' time along with a list of their data requirements.
#' @export Cant
Cant <- function(Data, timelimit = 1) {
  DLMdiag(Data, "not available", timelimit = timelimit)
}

#' @describeIn Can Internal function called by `Can` and `Cant`
#' @param command What to calculate? Character. Options = c("available", "not available", "needed")
#' @param reps The number of replicates for the MP
#' @param funcs1 A character vector of the MP names (optional)
#' @export
DLMdiag <- function(Data, command = c("available", "not available", "needed"), reps = 5, 
                    timelimit = 1, funcs1=NA, dev=FALSE) {
  command <- match.arg(command)
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  set.seed(101)
  Data <- updateMSE(Data)
  if (all(is.na(funcs1))) funcs1 <- avail("MP")
  isMP <- vapply(funcs1, function(x) inherits(get(x), "MP"), logical(1))
  if (any(!isMP)) stop(paste0("Not an MP: ", paste(funcs1[!isMP], collapse = ", ")))
  
  good <- rep(TRUE, length(funcs1))
  report <- rep("Passed test.", length(funcs1))
  test <- list()
  timey <- numeric(length(funcs1))
  
  rr <- try(slot(Data, "Misc"), silent = TRUE)
  if (class(rr) == "try-error") Data@Misc <- list()
  if (!dev) {
    ReqData <- DLMtool::ReqData
    builtin <- funcs1[funcs1 %in% ReqData$MP]
    custom <- funcs1[!funcs1 %in% ReqData$MP]
    inMPind <- which(funcs1 %in% builtin)
    cMPind <- which(funcs1 %in% custom)
    repp <- rep('', length(funcs1))
    # built in MPs (doesn't require 'dependencies' line)
    reqdata <-  ReqData[match(builtin, ReqData$MP),2]
    
    repp[inMPind] <- vapply(reqdata, match_slots, character(1), Data = Data, internal=TRUE)
    
    # custom MPs 
    temp <- lapply(custom, function(x) paste(format(match.fun(x)), collapse = " "))
    repp[cMPind] <- vapply(temp, match_slots, character(1), Data = Data)
    
    chk_needed <- nzchar(repp) # TRUE = has missing data
  } else {
    repp <- rep('', length(funcs1))
    chk_needed <- rep(FALSE, length(funcs1))
  }

  
  if (command == "needed") {
    # repp[!chk_needed] <- "All required data are present. Test MP with Can()"
    repp[chk_needed] <- paste0("Missing data: ", repp[chk_needed]) 
    ind <- nchar(repp)>0
    output <- matrix(repp[ind], ncol = 1, dimnames = list(funcs1[ind]))
    return(output)
  }
  
  report[chk_needed] <- paste0("Missing data: ", repp[chk_needed]) 
  good[chk_needed] <- FALSE
  
  for (y in 1:length(funcs1)) {
    if(!chk_needed[y]) {
      setTimeLimit(timelimit * 1.5)
      time1 <- Sys.time()
      test[[y]] <- tryCatch(do.call(funcs1[y], list(x = 1, Data = Data, reps = reps)), 
                            error = function(e) as.character(e))
      time2 <- Sys.time()
      setTimeLimit(Inf)
      timey[[y]] <- time2 - time1
      
      if (timey[[y]] > timelimit) {
        report[y] <- "Exceeded the user-specified time limit."
        good[y] <- FALSE
      # } else if (is.character(test[[y]]) | is.na(test[[y]])) { # Error message 
      } else if (is.character(test[[y]])) { # Error message 
        report[y] <- "MP returned an error. Check MP function and/or Data object."
        good[y] <- FALSE
      } else if (inherits(test[[y]], "Rec")) {
        # Rec_test <- vapply(slotNames("Rec"), function(x) NAor0(slot(test[[y]], x)), logical(1))
        Rec_test <- vapply(slotNames("Rec"), function(x) all(is.na(slot(test[[y]], x))), logical(1))
        if(all(Rec_test)) { # If all NAor0
          report[y] <- "Produced all NA scores. Check MP function and/or Data object."
          good[y] <- FALSE
        }
      }
    }
  }
  
  if (command == "available") {
    # return(matrix(report[good], ncol = 1, dimnames = list(funcs1[good])))
    return(funcs1[good])
  }
  if (command == "not available") {
    # return(matrix(report[!good], ncol = 1, dimnames = list(funcs1[!good])))
    return(matrix(c(funcs1[!good], report[!good]), ncol = 2))
  }
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
#' @examples 
#' Input(DLMtool::Cobia)
#' @export 
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
      if (length(runIn$Effort) > 0) Out[mm, 1] <- runIn$Effort[1]
      if (length(runIn$LR5) > 0) Out[mm, 2] <- runIn$LR5[1]
      if (length(runIn$LFR) > 0) Out[mm, 3] <- runIn$LFR[1]
      if (length(runIn$HS) > 0) Out[mm, 4] <- runIn$HS[1]
      Out[mm, 5:ncol(Out)] <- runIn$Spatial[1,]
    }
  }
  Out <- round(Out,2)
  print(Out, na.print='')
  invisible(Out)
  
}

needed <- function(Data, funcs) {
  rr <- try(slot(Data, "Misc"), silent = TRUE)
  if (class(rr) == "try-error") Data@Misc <- list()
  Data@Misc <- list()
  
  temp <- lapply(funcs, function(x) paste(format(match.fun(x)), collapse = " "))
  repp <- vapply(temp, match_slots, character(1), Data = Data)
  repp[!nzchar(repp)] <- "All data available. MP returns NA."
  matrix(repp, ncol = 1, dimnames = list(funcs))
}


# Internal function for:
## Required(): 
## needed(): matches slotnames of Data class to those that are required in an MP func
##           but also return NAor0 = TRUE
match_slots <- function(func, slotnams = paste0("Data@", slotNames("Data")), 
                        slots = slotNames("Data"), Data = NULL, internal=FALSE) {
  # check if each slotname in Data class is required in an MP (T/F)
  
  if(internal) ind_MP <- vapply(slots, grepl, logical(1), x = func)
  if(!internal) ind_MP <- vapply(slotnams, grepl, logical(1), x = func)
  if(!is.null(Data) && inherits(Data, "Data")) { # check if Data slots return NA or zero
    ind_NAor0 <- vapply(slots, function(x) all(is.na(slot(Data, x))), logical(1))
    repp <- slots[ind_MP & ind_NAor0] # returns slots where both tests are true
  } else {
    repp <- slots[ind_MP]
  }
  return(paste(repp, collapse = ", "))
}


#' @describeIn Can Identifies what data are needed to run
#' the MPs that are currently not able to run given a Data
#' object
#' @export Needed
Needed <- function(Data, timelimit = 1) {
  DLMdiag(Data, "needed", timelimit = timelimit)
}


#' Make stochastic variables certain for only one rep
#' 
#' As title.
#'
#' @param Data An object of class Data that has been run though TAC()
#' @author T. Carruthers
#' @export 
#' @keywords internal 
OneRep <- function(Data) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  Data@CV_Cat =  Data@CV_Ind = Data@CV_Rec = matrix(tiny, nrow=1, ncol=1)
  Data@CV_Dt = Data@CV_AvC = Data@CV_Mort = Data@CV_FMSY_M = Data@CV_BMSY_B0 = 
    Data@CV_Cref = Data@CV_Bref = Data@CV_Iref = Data@CV_Dep = Data@CV_Abun = 
    Data@CV_L50 = Data@CV_vbK = Data@CV_vbLinf = Data@CV_vbt0 = Data@CV_LFC = 
    Data@CV_LFS = Data@CV_wla = Data@CV_wlb = Data@CV_steep = tiny
  Data
}



#' A generic OFL plot for NOAA use
#' 
#' As title.
#' 
#' 
#' @param Data An object of class Data that has been run though TAC()
#' @param xlims x axis limits
#' @param perc The percentile of the OFL distribution to be plotted
#' @return A table of performance metrics.
#' @author T. Carruthers
#' @export 
plotOFL <- function(Data, xlims = NA, perc = 0.5) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  cols <- rep(c("black", "red", "green", "blue", "orange", "brown", "purple", 
                "dark grey", "violet", "dark red", "pink", "dark blue", "grey"),  4)
  ltys <- rep(1:4, each = 13)
  
  funcs <- Data@MPs
  nMPs <- length(funcs)
  
  if (is.na(xlims[1]) | length(xlims) != 2) {
    xlims <- quantile(Data@TAC, c(0.005, 0.9), na.rm = T)
    if (xlims[1] < 0) xlims[1] <- 0
  }
  if (!NAor0(Data@Ref)) {
    if (xlims[1] > Data@Ref) xlims[1] <- max(0, 0.98 * Data@Ref)
    if (xlims[2] < Data@Ref) xlims[2] <- 1.02 * Data@Ref
    if (xlims[2] > Data@Ref * 2) xlims[2] <- 2 * Data@Ref
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
  legend("topright", funcs, text.col = cols[cind], col = cols[cind], lty = 1, bty = "n", cex = 0.75)
  
  mtext(paste("OFL (", Data@Units, ")", sep = ""), 1, outer = F, line = 2.6)
  mtext(paste("Standardized relative frequency", sep = ""), 2, outer = F, line = 2.6)
  # mtext(paste('OFL calculation for
  # ',Data@Name,sep=''),3,outer=F,line=1)
  
}

#' Enlarge (replicate) a DLM data object to create an additional dimension for
#' simulation / sensitivity testing
#' 
#' Replicates position 1 data to multiple positions for sensitivity testing etc
#' 
#' 
#' @param Data A data-limited methods data object
#' @param nrep The number of positions to expand the DLM object to
#' @author T. Carruthers
#' @export
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


#' Sensitivity analysis
#' 
#' A function that determines the inputs for a given data-limited method of
#' class Output and then analyses the sensitivity of TAC estimates to 
#' marginal differences in each input. The range used for sensitivity is based 
#' on the user-specified CV for that input (e.g. CV_Mort, Mort)
#' 
#' @param Data A data-limited methods data object
#' @param MP A character string representing an MP applied in calculating the TAC recommendations in the DLM object
#' @param nsense The number of points over which to calculate the TAC (resolution)
#' @param reps The number of samples of the quota taken for the calculation of the TAC
#' @param perc The percentile of the sample TAC
#' @param ploty A logical switch, (T/F, should a plot be drawn?)
#'
#' @author T. Carruthers
#' @export 
#' @examples 
#' \dontrun{
#' Data <- Sense(DLMtool::Cobia, "AvC")
#' }
Sense <- function(Data, MP, nsense = 6, reps = 100, perc = c(0.05, 0.5, 0.95), ploty = T) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  
  DLM_data2 <- Data
  nm <- deparse(substitute(DLM_data2))
  # refTAC <- quantile(getTAC(DLM_data2, MP, reps)[[1]], perc, na.rm = T)
  refTAC <- quantile(applyMP(DLM_data2, MPs=MP, reps=reps)[[1]][[1]]$TAC, perc, na.rm = T)
  
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
      attr(Data, slots[i])[ind, ] <- attr(Data, slots[i])[ind, ] * qlnorm(pss, mconv(1, cv), sdconv(1, cv))
      vals[i, ] <- qlnorm(pss, mconv(1, cv), sdconv(1, cv))
    }
  }
  
  # TACa <- getTAC(Data, MPs = MP, reps = reps)[[1]]
  TACa <- applyMP(Data, MPs=MP, reps=reps)[[1]][[1]]$TAC
  TACa <- apply(TACa, 2, quantile, p = perc, na.rm = T)
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



#' Calculate TAC recommendations for more than one MP
#' 
#' A function that returns the stochastic TAC recommendations from a vector of
#' data-limited MPs (Output) given a data-limited data object Data
#' 
#' 
#' @param Data A data-limited methods data object
#' @param MPs optional vector of MP names
#' @param reps Number of repititions
#' @param timelimit The maximum time (seconds) taken to complete 10 reps
#' @author T. Carruthers
#' @examples 
#' \dontrun{
#' Data <- TAC(DLMtool::Cobia)
#' plot(Data)
#' }
#' @export 
TAC <- function(Data, MPs = NA, reps = 100, timelimit = 1) {
  if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
  Data <- updateMSE(Data)
  nm <- deparse(substitute(Data))
  PosMPs <- Can(Data, timelimit = timelimit)
  PosMPs <- PosMPs[PosMPs %in% avail("Output")]
  Data@PosMPs <- PosMPs
  if (!is.na(MPs[1])) Data@MPs <- MPs[MPs %in% PosMPs]
  if (is.na(MPs[1])) Data@MPs <- PosMPs
  funcs <- Data@MPs
  
  if (length(funcs) == 0) {
    stop("None of the methods 'MPs' are possible given the data available")
  } else {
    Data <- applyMP(Data, MPs = funcs, reps)[[2]]
    return(Data)
    # assign(nm,DLM,envir=.GlobalEnv)
  }
  
}


#' Join Data objects present in a list 
#' 
#' A function that combined a list of data objects into a single data object (same dimensions but can have different numbers of simulations)
#' 
#' 
#' @param DataList A list of data objects of identical dimension (except for simulation)
#' @author T. Carruthers
#' @export 
joinData<-function(DataList){
  
  if (class(DataList) != "list") stop("DataList must be a list")
  if (length(DataList) < 2) stop("DataList list doesn't contain multiple MSE objects")
  
  Data<-DataList[[1]]
  nD<-length(DataList)
 
  slots<-slotNames(Data)
  slots_identical <- c("Name", "Common_Name", "Species", "Region", "Year", "MaxAge", "Units", "Ref_type", "PosMPs", "MPs", "nareas", "LHYear")
  #slots <- slots[!slots%in%c("Name","Ref","OM","MaxAge","CAL_bins","Year","Units","Ref","Ref_type","Log","params","PosMPs","MPs","Obs","Misc","nareas","LHYear")]
  slots <- slots[!slots %in% slots_identical]
  
  nslots<-length(slots)
  getslot<-function(obj,name)slot(obj,name) # weird issue with namespace conflict and the @Cat slot
  getslotclass<-function(obj,name)class(slot(obj,name))
  sclass<-sapply(1:nslots,function(x,obj,slots)getslotclass(obj,slots[x]),obj=DataList[[1]],slots=slots)
  #getdim<-function(x){
  #  dim<-dim(x)
  #  if(is.null(dim))dim=length(x)
  #  dim
  #}
  #sdims<-sapply(1:nslots,function(x,obj,slots)getdim(getslot(obj,slots[x])),obj=DataList[[1]],slots=slots)
  #nsims<-sapply(1:nD,function(x,DataList)length(DataList[[x]]@Dt),DataList)

  for (sn in 1:nslots){
    templist<-lapply(DataList,getslot,name=slots[sn])
     
    if (sclass[sn] == "numeric"|sclass[sn]=="integer") {
      if (slots[sn] == "CAL_bins") {
        nbin <- vapply(templist, length, numeric(1))
        attr(Data, slots[sn]) <- templist[[which.max(nbin)]]
      } else {
        attr(Data, slots[sn]) <- unlist(templist)
      }
    } else if (sclass[sn]== "matrix"|sclass[sn]=="array") {
      
      if(slots[sn] == "CAL") {
        nbin <- vapply(templist, function(x) dim(x)[3], numeric(1))
        templist2 <- vector("list", nD)
        for (i in 1:nD) {
          templist2[[i]] <- array(0, dim = c(dim(templist[[i]])[1:2], max(nbin)))
          templist2[[i]][ , , 1:nbin[i]] <- templist[[i]]
        }
        if (all(diff(do.call("rbind", lapply(templist2, dim))[,2]) == 0)) {
          # arrays may be different dimensions if MPs fail
          attr(Data, slots[sn]) <- abind(templist2, along=1)
        }
      } else {
        if (all(diff(do.call("rbind", lapply(templist, dim))[,2]) == 0)) {
          # arrays may be different dimensions if MPs fail
          attr(Data, slots[sn]) <- abind(templist, along=1)
        }
      }
      
    } else if (sclass[sn] == "list") {
      attr(Data, slots[sn]) <- do.call(c, templist)
    } else if (sclass[sn] == "data.frame") {
      attr(Data, slots[sn]) <- do.call(rbind, templist)
    }
  }
  
  for (sn in 1:length(slots_identical)) {
    templist <- lapply(DataList, getslot, name = slots_identical[sn])
    attr(Data, slots_identical[sn]) <- unique(do.call(c, templist))
  }
    
  #checkdims<-sapply(1:nslots,function(x,obj,slots)getdim(getslot(obj,slots[x])),obj=Data,slots=slots)
  Data  

}


#' Find the Management Procedures that use a particular data slot
#'
#' @param slot A slot from an object of class `Data`. Character string.
#' @param silent Logical. Should messages be printed? 
#'
#' @return A character string of MPs that use the slot.
#' @author A. Hordyk
#' @export
#'
#' @examples
#' Uses("Mort")
Uses <- function(slot, silent=FALSE) {
  if (class(slot) !="character") stop("Slot must be character", call. = FALSE)
  if(length(slot)>1) stop("Slot must be length 1", call. = FALSE)
  if (!slot %in% slotNames('Data')) stop("Slot is not a valid slot in Data object. Use slotNames('Data')", call.=FALSE)
  MPs <- avail("MP")
  List <- lapply(seq_along(MPs), function(x) Required(MPs[x]))
  df <- data.frame(matrix(unlist(List), nrow=length(MPs), byrow=T), stringsAsFactors = FALSE)
  mps <- df[grepl(slot, df[,2]),1]
  if (length(mps) >0) {
    if(!silent)
      message("MPs that require Data slot '" , slot, "' are:")
    return(mps)
  } else {
    if(!silent)
      message("No MPs used Data slot '" , slot, "'.")
    return(NULL)
  }
}








# parallelMPs <- function(x, Data, reps, MPs, ss) sapply(ss, MPs[x], Data, reps = reps)

# getTAC <- function(Data, MPs = NA, reps = 100) {
#   if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
#   Data <- updateMSE(Data)
#   nsims <- length(Data@Mort)
#   nMPs <- length(MPs)
#   TACa <- array(NA, dim = c(nMPs, reps, nsims))
#   
#   if (!sfIsRunning()) {
#     for (ff in 1:nMPs) {
#       temp <- sapply(1:nsims, MPs[ff], Data = Data, reps = reps)
#       tt <- sapply(1:1, AvC,  Data=Data, reps=reps)
#       
#       tt@TAC
#       
#       
#       if (mode(temp) == "numeric") 
#         TACa[ff, , ] <- temp
#       if (mode(temp) == "list") {
#         TACa[ff, , ] <- unlist(temp[1, ])
#         for (x in 1:nsims) Data@Misc[[x]] <- temp[2, x][[1]]
#       }
#     }
#   } else {
#     sfExport("Data") 
#     if (nsims < 8) {
#       sfExport(list = c("MPs", "reps"))
#       for (ss in 1:nsims) {
#         temp <- (sfSapply(1:length(MPs), parallelMPs, Data = Data, 
#           reps = reps, MPs = MPs, ss = ss))
#         if (mode(temp) == "list") {
#           Lens <- unlist(lapply(temp, length))
#           for (X in 1:length(Lens)) {
#           Classes <- unlist(lapply(temp[, X][[1]], class))
#           if (length(unique(Classes)) == 1) {
#             # no Misc object
#             TACa[X, , ss] <- unlist(temp[, X])
#           } else {
#             # a Misc object is include
#             ind <- which(Classes == "list")
#             TACa[X, , ss] <- unlist(temp[, X][[1]][1:(ind - 1), 
#             ])
#             Data@Misc[[ss]] <- temp[, X][[1]][ind, ]
#           }
#           }
#         } else {
#           temp <- matrix(temp, nrow = nMPs, ncol = reps, byrow = TRUE)
#           TACa[, , ss] <- temp
#         }
#       }
#     } else {
#       for (ff in 1:nMPs) {
#         temp <- sfSapply(1:nsims, MPs[ff], Data = Data, 
#           reps = reps)
#         if (mode(temp) == "numeric") 
#           TACa[ff, , ] <- temp
#         if (mode(temp) == "list") {
#           TACa[ff, , ] <- unlist(temp[1, ])
#           for (x in 1:nsims) Data@Misc[[x]] <- temp[2, x][[1]]
#         }
#       }
#     }
#   }
#   for (ff in 1:nMPs) {
#     if (sum(is.na(TACa[ff, , ])) > sum(!is.na(TACa[ff, , ]))) {
#       # only plot if there are sufficient non-NA TAC samples
#       print(paste("Method ", MPs[ff], " produced greater than 50% NA values", 
#         sep = ""))
#     }
#   }
#   out <- list(TACa, Data)
#   return(out)
# }
# 
# 

# #' Conduct stock assessment
# #' 
# #' A wrapper function that gets the OFL recommendation in cases where a method
# #' of DLM quota has been specified
# #' 
# #' 
# #' @usage Sam(Data, MPs = NA, reps = 100, perc = 0.5)
# #' @param Data A data-limited methods data object
# #' @param MPs A character vector of methods of DLM quota, DLM space or DLM size
# #' @param reps The number of samples of quota recommendations by method
# #' @param perc quantile of the recommendation to use
# #' @author T. Carruthers
# #' @export Sam
# Sam <- function(Data, MPs = NA, reps = 100, perc = 0.5) {
#   if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
#   Data <- updateMSE(Data)
#   nm <- deparse(substitute(DLM))
#   Data@PosMPs <- MPs
#   funcs <- Data@PosMPs
#   nMPs <- length(funcs)
#   Data@MPs <- funcs
#   temp <- getTAC(Data, MPs = funcs, reps)
#   TACa <- temp[[1]]
#   Data <- temp[[2]]
#   nsim <- length(Data@Mort)
#   ref <- array(rep(Data@Ref, nMPs), c(nsim, nMPs))
#   TACm <- apply(TACa, c(3, 1), quantile, p = perc, na.rm = T)
#   TACbias <- (TACm - ref)/ref * 100
#   POF <- round(apply(TACbias > 0, 2, sum)/length(Data@Mort) * 100, 
#     1)
#   Data@TAC <- TACa
#   Data@TACbias <- TACbias
#   Data
# }
# 
