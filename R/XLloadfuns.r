
#' Read in operating model parameters from Excel spreadsheet
#' 
#' A function to read in operating model parameters from an Excel spreadsheet
#' with tabs named following specific convention.
#' 
#' The Excel spreadsheet must have tabs named with the following convention.
#' For example if \code{stkname} is 'myFish', the Stock parameters are in a tab
#' named 'myFishStock', Fleet parameters in a tab named 'myFishFleet', 
#' Observation parameters in a tab named 'myFishObs', and Implementation in 'myFishImp'. 
#' All tabs (Stock, Fleet, Obs, and Imp) must be present for a single stock. You can have multiple
#' stocks in a single spreadsheet, provided that the stock names are different.
#' 
#' @usage OM_xl(fname, stkname, fpath = '', saveCSV = FALSE)
#' @param fname Name of the Excel spreadsheet file. Must include file
#' extension.
#' @param stkname Name of the Stock.
#' @param fpath Full file path, if file is not in current working directory
#' @param saveCSV Do you also want to save the Stock, Fleet and Observation
#' parameters to CSV files?
#' @return A object of class OM
#' @author A. Hordyk
#' @examples
#' 
#' \dontrun{
#' OM <- OM_xl(fname='OMTables.xlsx', stkname='myFish')
#' }
#' 
#' @export 
OM_xl <- function(fname, stkname, fpath = "", saveCSV = FALSE) {
  infile <- paste0(fpath, fname)  # full path and name 
  shtname <- readxl::excel_sheets(infile)  # names of the sheets 
  # Stock
  stock <- readxl::read_excel(infile, sheet = grep(paste0(stkname, "Stock"), shtname), col_names = FALSE)
  stock <- as.data.frame(stock)
  tmpfile <- paste0(fpath, stkname, "Stock.csv")
  if (file.exists(tmpfile)) unlink(tmpfile)
  writeCSV(inobj = stock, tmpfile, objtype = "Stock")
  tmpstock <- new("Stock", tmpfile)
  if (!saveCSV) unlink(tmpfile)
  
  # Fleet
  index <- which(pmatch(shtname, paste0(stkname, "Fleet")) == 1)
  if (length(index) > 1) stop("More than one match")
  fleet <- readxl::read_excel(infile, sheet = index, col_names = FALSE)
  fleet <- as.data.frame(fleet)
  tmpfile <- paste0(fpath, stkname, "Fleet.csv")
  if (file.exists(tmpfile)) unlink(tmpfile)
  writeCSV(inobj = fleet, tmpfile, objtype = "Fleet")
  tmpfleet <- new("Fleet", tmpfile)
  if (!saveCSV) unlink(tmpfile)
  
  # Observation
  index <- which(pmatch(shtname, paste0(stkname, "Obs")) == 1)
  if (length(index) > 1) stop("More than one match")
  obs <- readxl::read_excel(infile, sheet = index, col_names = FALSE)
  obs <- as.data.frame(obs)
  tmpfile <- paste0(fpath, stkname, "Obs.csv")
  if (file.exists(tmpfile)) unlink(tmpfile)
  writeCSV(inobj = obs, tmpfile, objtype = "Obs")
  tmpobs <- new("Obs", tmpfile)
  if (!saveCSV) unlink(tmpfile)
 
  # Implementation
  index <- which(pmatch(shtname, paste0(stkname, "Imp")) == 1)
  if (length(index) > 1)  stop("More than one match")
  obs <- readxl::read_excel(infile, sheet = index, col_names = FALSE)
  obs <- as.data.frame(obs)
  tmpfile <- paste0(fpath, stkname, "Imp.csv")
  if (file.exists(tmpfile)) unlink(tmpfile)
  writeCSV(inobj = obs, tmpfile, objtype = "Imp")
  tmpimp <- new("Imp", tmpfile)
  if (!saveCSV) unlink(tmpfile)
  
  # Operating Model
  OM <- new("OM", Stock = tmpstock, Fleet = tmpfleet, Obs = tmpobs, Imp=tmpimp)
  OM
}


#' Read in Data object from Excel spreadsheet
#' 
#' A function to read in Data object from an Excel spreadsheet
#' with tabs named following specific convention.
#' 
#' The Excel spreadsheet must have tabs named with the following convention.
#' For example if \code{stkname} is 'myFish', the Data parameters are in a tab
#' named 'myFishData'.
#' 
#' @param fname Name of the Excel spreadsheet file. Must include file
#' extension.
#' @param stkname Name of the Stock.
#' @param fpath Full file path, if file is not in current working directory
#' @param saveCSV Do you also want to the Data parameters to a CSV file?
#' @return A object of class Data
#' @author A. Hordyk
#' @examples
#' 
#' \dontrun{
#' OM <- OM_xl(fname='OMTables.xlsx', stkname='myFish')
#' }
#' 
#' @export 
Data_xl <- function(fname, stkname, fpath = "", saveCSV = FALSE) {
  infile <- paste0(fpath, fname)  # full path and name 
  shtname <- readxl::excel_sheets(infile)  # names of the sheets 
  # Data
  index <- which(pmatch(shtname, paste0(stkname, "Data")) == 1)
  if (length(index) > 1)  stop("More than one match")
  data <- readxl::read_excel(infile, sheet = index, col_names = FALSE)
  data <- as.data.frame(data)
  tmpfile <- paste0(fpath, stkname, "Data.csv")
  if (file.exists(tmpfile)) unlink(tmpfile)
  writeCSV(inobj = data, tmpfile, objtype = "Data")
  tmpimp <- new("Data", tmpfile)
  data <- new("Data",tmpfile)
  if (!saveCSV) unlink(tmpfile)
  return(data)
  
}

#' Read in feasibility parameters from Excel spreadsheet
#' 
#' A function to read in feasibility parameters from an Excel spreadsheet with
#' tabs named following specific convention
#' 
#' The Excel spreadsheet must have tabs named with the following convention.
#' For example if \code{stkname} is 'myFish', the tab must be named
#' 'myFishFease,
#' 
#' @usage Fease_xl(fname, stkname, fpath = '', saveCSV = FALSE)
#' @param fname Name of the Excel spreadsheet file. Must include file
#' extension.
#' @param stkname Name of the Stock.
#' @param fpath Full file path, if file is not in current working directory
#' @param saveCSV Do you also want to save the Stock, Fleet and Observation
#' parameters to CSV files?
#' @return A object of class Fease
#' @author A. Hordyk
#' @examples
#' 
#'  \dontrun{
#'  myFease <- Fease_xl(fname='FeaseTables.xlsx', stkname='myFish')
#' }
#' 
#' @export Fease_xl
Fease_xl <- function(fname, stkname, fpath = "", saveCSV = FALSE) {
  infile <- paste0(fpath, fname)  # full path and name 
  shtname <- readxl::excel_sheets(infile)  # names of the sheets 
  # Fease
  feasedat <- readxl::read_excel(infile, sheet = grep(paste0(stkname, "Fease"), 
    shtname), col_names = FALSE)
  feasedat <- feasedat[, 1:2]
  tmpfile <- paste0(fpath, stkname, "Fease.csv")
  if (file.exists(tmpfile)) 
    unlink(tmpfile)
  writeCSV(inobj = feasedat, tmpfile, objtype = "Fease")
  fease <- new("Fease", tmpfile)
  if (!saveCSV) 
    unlink(tmpfile)
  
  fease
}



#' Internal function to write CSVs for objects
#' 
#' Used internally in the DLMtool package to write CSV files from an existing
#' DLMtool object
#' 
#' 
#' @param inobj A object of class Stock, Fleet, Obs, Imp, Data, OM, or
#' Fease
#' @param tmpfile The full file path and name for the saved CSV file
#' @param objtype The class corresonding to the \code{inobj}
#' @author A. Hordyk
#' @export writeCSV
writeCSV <- function(inobj, tmpfile = NULL, objtype = c("Stock", "Fleet", 
  "Obs", "Imp", "Data", "OM", "Fease")) {
  objtype <- match.arg(objtype)
  
  for (X in 1:nrow(inobj)) {
    indat <- inobj[X, ]
    index <- which(!is.na(indat))
    index <- 2:max(index)
    if (X == 1) 
      write(do.call(paste, c(indat[1], as.list(indat[index]), sep = ",")), tmpfile, 1)
    if (X > 1) 
      write(do.call(paste, c(indat[1], as.list(indat[index]), sep = ",")), tmpfile, 1, append = TRUE)
  }
  
  
  # tmpobj <- new(objtype)
  # sn <- slotNames(tmpobj)
  # ind <- which(inobj[, 1] %in% sn == FALSE)
  # if (length(ind) > 0) {
  #   message("Input file names don't match slot names for ", objtype, " object")
  #   message("Unknown input name:", inobj[ind, 1])
  #   stop("Check the input file row names")
  # }
  # for (X in seq_along(sn)) {
  #   ind <- match(sn[X], inobj[, 1])
  #   if (!is.na(ind)) {
  #     indat <- inobj[ind, ]
  #     index <- which(!is.na(indat))
  #     index <- 2:max(index)
  #     if (X == 1) 
  #       write(do.call(paste, c(sn[X], as.list(indat[index]), sep = ",")), 
  #         tmpfile, 1)
  #     if (X > 1) 
  #       write(do.call(paste, c(sn[X], as.list(indat[index]), sep = ",")), 
  #         tmpfile, 1, append = TRUE)
  #   }
  # }
}


