#' Initialize Operating Model
#'
#' Generates an Excel spreadsheet and a text file in the current working directory for 
#' specifying and documenting a DLMtool Operating Model.
#' 
#' @param name The name of the Excel and text file to be created in the working directory (character)
#' @param templates An optional named list of existing DLMtool Stock, Fleet, Obs, or Imp objects
#' to use as templates for the Operating Model.
#' @param overwrite Logical. Should files be overwritten if they already exist?
#'
#' @return An xlsx and text file are created in the working directory.  
#' @export
#' @author A. Hordyk
#'
#' @examples
#' \dontrun{
#' # Create an Excel OM template and text file called 'myOM.xlxs' and 'myOM.txt': 
#' OMinit('myOM')
#' 
#' # Create an Excel OM template and text file called 'myOM.xlxs' and 'myOM.txt', using
#' the Stock object 'Herring' as a template: 
#' OMinit('myOM', templates=list(Stock='Herring'))
#' 
#' # Create an Excel OM template and text file called 'myOM.xlxs' and 'myOM.txt', using
#' the Stock object 'Herring', and Obs object 'Generic_obs' as templates: 
#' OMinit('myOM', templates=list(Stock='Herring', Obs='Generic_obs'), overwrite=TRUE)
#' }
#' 
OMinit <- function(name="Example", templates=NULL, overwrite=FALSE) {
  if (class(name) != 'character') stop("name must be text")
  ## Write Excel skeleton ####
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, sheetName = "OM", gridLines = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "Stock", gridLines = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "Fleet", gridLines = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "Obs", gridLines = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "Imp", gridLines = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "Data", gridLines = TRUE)
  
  if (name == "Example") {
    message("Creating Example.xlsx")
    message("Using example templates")
    templates <- list(Stock=DLMtool::Albacore, Fleet=DLMtool::DecE_Dom, 
                      Obs=DLMtool::Generic_obs, Imp=DLMtool::Perfect_Imp)
    overwrite <- TRUE
  }
  ObTemplates <- ObjTemps(templates)
  if (!is.null(ObTemplates) && class(ObTemplates) == 'list') {
    nm <- names(ObTemplates)  
    message("\n\nUsing Object Templates:")
    for (X in nm) {
      message(ObTemplates[[X]]@Name)
    }
  }
  
  # Load the Excel File ####
  if (nchar(tools::file_ext(name)) == 0) {
    nameNoExt <- name
    name <- paste0(name, ".xlsx")
  } else {
    ext <- tools::file_ext(name)
    if (!ext %in% c("xlsx", "xls")) stop("File extension must be 'xlsx' or 'xls'")
    nameNoExt <- tools::file_path_sans_ext(name)
  }

  message("Creating ", name, " in ", getwd())
  
  # OM sheet####
  df <- data.frame(Slot=c("Name", "nsim", "proyears", "interval", "pstar", "maxF", "reps"))
  
  # write slots 
  openxlsx::writeData(wb, sheet = "OM", x = df, 
                      startCol = 1, startRow = 1,
                      colNames = TRUE, rowNames = FALSE, 
                      withFilter = FALSE,
                      keepNA = FALSE) 
  
  
  
  # write slots 
  df <- data.frame(Defaults=nameNoExt)
  openxlsx::writeData(wb, sheet = "OM", x = df, 
                      startCol = 2, startRow = 1,
                      colNames = TRUE, rowNames = FALSE, 
                      withFilter = FALSE,
                      keepNA = FALSE) 
  
  df <- data.frame(Values=c( 48, 50, 4, 0.5, 0.8, 1))
  openxlsx::writeData(wb, sheet = "OM", x = df, 
                      startCol = 2, startRow = 3,
                      colNames = FALSE, rowNames = FALSE, 
                      withFilter = FALSE,
                      keepNA = FALSE) 
  
  
  # Stock sheet ####
  df <- data.frame(Slot=slotNames("Stock"))
  
  # write slots 
  openxlsx::writeDataTable(wb, sheet = "Stock", x = df, 
                           startCol = 1, startRow = 1,
                           tableStyle = "none",
                           colNames = TRUE, rowNames = FALSE, 
                           withFilter = FALSE,
                           bandedRows = FALSE,
                           bandedCols = FALSE,
                           keepNA = FALSE,
                           firstColumn =TRUE)
  
  openxlsx::setColWidths(wb, sheet = "Stock", cols = 1, widths = 'auto')
  
  # loop through slot values if Obj template provided
  if (!is.null(ObTemplates$StockTemp)) {
    obj <- ObTemplates$StockTemp
    slots <- slotNames(obj)
    
    for (sl in seq_along(slots)) {
      val <- slot(obj, slotNames("Stock")[sl])
      ln <- length(val)
      if (ln >0 && !is.na(ln)) {
        df <- data.frame(t(val))
        openxlsx::writeData(wb, sheet = "Stock", x = df, 
                            startCol = 2, startRow = sl+1,
                            colNames = FALSE, rowNames = FALSE, 
                            withFilter = FALSE,
                            keepNA = FALSE)         
      }
    }
  }
  
  
  # Fleet sheet ####
  df <- data.frame(Slot=slotNames("Fleet"))
  
  # write slots 
  openxlsx::writeDataTable(wb, sheet = "Fleet", x = df, 
                           startCol = 1, startRow = 1,
                           tableStyle = "none",
                           colNames = TRUE, rowNames = FALSE, 
                           withFilter = FALSE,
                           bandedRows = FALSE,
                           bandedCols = FALSE,
                           keepNA = FALSE,
                           firstColumn =TRUE)
  
  openxlsx::setColWidths(wb, sheet = "Fleet", cols = 1, widths = 'auto')
  
  # loop through slot values if Obj template provided
  if (!is.null(ObTemplates$FleetTemp)) {
    obj <- ObTemplates$FleetTemp
    slots <- slotNames(obj)
    
    for (sl in seq_along(slots)) {
      val <- slot(obj, slotNames("Fleet")[sl])
      ln <- length(val)
      if (ln >0 && !is.na(ln)) {
        df <- data.frame(t(val))
        openxlsx::writeData(wb, sheet = "Fleet", x = df, 
                            startCol = 2, startRow = sl+1,
                            colNames = FALSE, rowNames = FALSE, 
                            withFilter = FALSE,
                            keepNA = FALSE)         
      }
      
    }
  }
  
  # Obs sheet ####
  df <- data.frame(Slot=slotNames("Obs"))
  
  # write slots 
  openxlsx::writeDataTable(wb, sheet = "Obs", x = df, 
                           startCol = 1, startRow = 1,
                           tableStyle = "none",
                           colNames = TRUE, rowNames = FALSE, 
                           withFilter = FALSE,
                           bandedRows = FALSE,
                           bandedCols = FALSE,
                           keepNA = FALSE,
                           firstColumn =TRUE)
  
  openxlsx::setColWidths(wb, sheet = "Obs", cols = 1, widths = 'auto')
  
  # loop through slot values if Obj template provided
  if (!is.null(ObTemplates$ObsTemp)) {
    obj <- ObTemplates$ObsTemp
    slots <- slotNames(obj)
    
    for (sl in seq_along(slots)) {
      val <- slot(obj, slotNames("Obs")[sl])
      ln <- length(val)
      if (ln >0 && !is.na(ln)) {
        df <- data.frame(t(val))
        openxlsx::writeData(wb, sheet = "Obs", x = df, 
                            startCol = 2, startRow = sl+1,
                            colNames = FALSE, rowNames = FALSE, 
                            withFilter = FALSE,
                            keepNA = FALSE)         
      }
    }
  }
  
  # Imp sheet ####
  df <- data.frame(Slot=slotNames("Imp"))
  
  # write slots 
  openxlsx::writeDataTable(wb, sheet = "Imp", x = df, 
                           startCol = 1, startRow = 1,
                           tableStyle = "none",
                           colNames = TRUE, rowNames = FALSE, 
                           withFilter = FALSE,
                           bandedRows = FALSE,
                           bandedCols = FALSE,
                           keepNA = FALSE,
                           firstColumn =TRUE)
  
  openxlsx::setColWidths(wb, sheet = "Imp", cols = 1, widths = 'auto')
  
  # loop through slot values if Obj template provided
  if (!is.null(ObTemplates$ImpTemp)) {
    obj <- ObTemplates$ImpTemp
    slots <- slotNames(obj)
    
    for (sl in seq_along(slots)) {
      val <- slot(obj, slotNames("Imp")[sl])
      ln <- length(val)
      if (ln >0 && !is.na(ln)) {
        df <- data.frame(t(val))
        openxlsx::writeData(wb, sheet = "Imp", x = df, 
                            startCol = 2, startRow = sl+1,
                            colNames = FALSE, rowNames = FALSE, 
                            withFilter = FALSE,
                            keepNA = FALSE)         
      }
    }
  }
  
  # Data sheet ####
  df <- data.frame(Slot=slotNames("Data"))
  
  # write slots 
  openxlsx::writeDataTable(wb, sheet = "Data", x = df, 
                           startCol = 1, startRow = 1,
                           tableStyle = "none",
                           colNames = TRUE, rowNames = FALSE, 
                           withFilter = FALSE,
                           bandedRows = FALSE,
                           bandedCols = FALSE,
                           keepNA = FALSE,
                           firstColumn =TRUE)
  
  openxlsx::setColWidths(wb, sheet = "Data", cols = 1, widths = 'auto')
  
  
  
  
  # Write Excel file ####
  if (file.exists(name) & !overwrite) {
    stop(name, " already exists.\n Use 'overwrite=TRUE'. \nNOTE: this will overwrite both .xlsx and txt files if they exist.")
  } else {
    options(warn=2)
    tryWrite <- try(openxlsx::saveWorkbook(wb, name, overwrite = overwrite), 
                    silent=TRUE) ## save to working directory
    options(warn=0)
    if (tryWrite != 1) stop("Can't write to ", name, ". If the file open?")
  }
  
  
  
  ## Write Text skeleton ####
  message("Creating ", nameNoExt, ".txt in ", getwd())
  Textfile <- paste0(nameNoExt, ".txt")
  if (file.exists(Textfile) & !overwrite) {
    stop(Textfile, " already exists.\n Use 'overwrite=TRUE'.")
  } else {
    tt <- file.create(Textfile)
  }
  
  # Title
  cat("\n# Title\n", sep="", append=TRUE, file=Textfile) 
  cat("Title. One line only.\n", sep="", append=TRUE, file=Textfile) 
  
  # Subtitle - optional 
  cat("\n# Subtitle\n", sep="", append=TRUE, file=Textfile) 
  cat("Optional. Subtitle. One line only. Delete text and heading if not required.\n", 
      sep="", append=TRUE, file=Textfile) 
  
  
  # Author(s) 
  cat("\n# Author(s)\n", sep="", append=TRUE, file=Textfile) 
  cat("Name and contact details (e.g email, affiliation) for each author.\n", 
      sep="", append=TRUE, file=Textfile) 
  cat("One line per author.\n", 
      sep="", append=TRUE, file=Textfile) 
  
  # # Affiliation/Email - optional 
  # cat("\n# Affiliation/Email\n", sep="", append=TRUE, file=Textfile) 
  # cat("Affiliation and/or email for each author. One line for each author.\n", 
  #     sep="", append=TRUE, file=Textfile) 
  # cat("Will be recycled for each author if more authors than affiliations.\n", 
  #     sep="", append=TRUE, file=Textfile) 
  
  # Date - optional 
  cat("\n# Date\n", sep="", append=TRUE, file=Textfile) 
  cat("Optional. Date that the operating model was created. If none provided, today's date will be used.\n", 
      sep="", append=TRUE, file=Textfile) 
 
  # Introduction ####
  cat("\n# Introduction\n\n", sep="", append=TRUE, file=Textfile)
  
  cat("## Completing the OM Documentation\n", sep="", append=TRUE, file=Textfile)
  cat("This document is used to generate a HTML OM report document.\n\n", sep="", append=TRUE, file=Textfile)
  cat("The document is separated into 7 sections:\n", sep="", append=TRUE, file=Textfile)
  cat("1. Introduction (this section)\n", sep="", append=TRUE, file=Textfile)
  cat("2. Custom Parameters (optional)\n", sep="", append=TRUE, file=Textfile)
  cat("3. Stock Parameters\n", sep="", append=TRUE, file=Textfile)
  cat("4. Fleet Parameters\n", sep="", append=TRUE, file=Textfile)
  cat("5. Obs (Observation) Parameters\n", sep="", append=TRUE, file=Textfile)
  cat("6. Imp (Implementation) Parameters\n", sep="", append=TRUE, file=Textfile)
  cat("7. References\n\n", sep="", append=TRUE, file=Textfile)
  
  cat("The Introduction section is used to briefly describe the fishery and the details of the Operating Model.\n", sep="", append=TRUE, file=Textfile)
  cat("It should include an explanation for the OM parameters:\n ", sep="", append=TRUE, file=Textfile)
  cat("* nsim: the number of simulations.\n ", sep="", append=TRUE, file=Textfile)
  cat("* proyears: the number of projectio years.\n ", sep="", append=TRUE, file=Textfile)
  cat("* interval: the management interval.\n ", sep="", append=TRUE, file=Textfile)
  cat("* pstar: the percentile of the sample of the management recommendation for each method.\n ", sep="", append=TRUE, file=Textfile)
  cat("* maxF: the maximum instantaneous fishing mortality rate that may be simulated for any given age class.\n ", sep="", append=TRUE, file=Textfile)
  cat("* reps: the number of samples of the management recommendation for each method.\n\n", sep="", append=TRUE, file=Textfile)
  
  cat("The Custom Parameters section is optional and should only be included if the cpars feature of DLMtool is used in the OM.\n", sep="", append=TRUE, file=Textfile)
  cat("Delete both the heading and the text in this section if `cpars` are not used.\n\n", sep="", append=TRUE, file=Textfile)
  
  cat("The Stock, Fleet, Obs, and Imp sections include each slot in these components of the OM object.\n", sep="", append=TRUE, file=Textfile)
  cat("Provide details (including references where applicable) for the choice of values for each slot below the corresponding slot name (e.g., ## M).\n", sep="", append=TRUE, file=Textfile)
  cat("For example: \n", sep="", append=TRUE, file=Textfile)
  cat("**M**\n", sep="", append=TRUE, file=Textfile)
  cat("An explanation for the values of the natural mortality rate in the OM (Smith et al. 1999).\n\n", sep="", append=TRUE, file=Textfile)
  cat("You do not need to include the actual values from the OM. These will be included automatically in the final compiled document.\n\n", sep="", append=TRUE, file=Textfile)
  cat("References should be included in the 'References' section at the end of the document.\n\n", sep="", append=TRUE, file=Textfile)
  
  cat("Once complete, this text file will be compiled into an OM Report Document.\n", sep="", append=TRUE, file=Textfile)
  cat("This text file is linked to the Excel spreadsheet that was generated with the same file name.\n", sep="", append=TRUE, file=Textfile)
  cat("It serves as a single documentation source for a DLMtool OM, and should be updated whenever parameter values in the OM spreadsheet are updated.\n", sep="", append=TRUE, file=Textfile)
  
  cat("\n## Tips on filling this Document\n\n", sep="", append=TRUE, file=Textfile)
  cat("This document is uses the Markdown format. All first and second level headings have been provided, and in general you\n", sep="", append=TRUE, file=Textfile)
  cat("should only need to enter plain text.\n\n", sep="", append=TRUE, file=Textfile)
  
  cat("You can have multiple paragraphs throughout the document.\n\n", sep="", append=TRUE, file=Textfile)
  
  cat("The Introduction and Custom Parameters sections also support second and third level headings.\n\n", sep="", append=TRUE, file=Textfile)
  
  cat("## An example Second level heading\n\n", sep="", append=TRUE, file=Textfile)
  
  cat("### An example third level heading\n\n", sep="", append=TRUE, file=Textfile)
  
  cat("Delete all text below 'Introduction' and replace with a description of the OM.\n\n", sep="", append=TRUE, file=Textfile)
  
  
  # Cpars ####
  cat("\n\n# Custom Parameters\n", sep="", append=TRUE, file=Textfile)  
  cat("Optional. Only required if the `cpars` feature is used in the OM.\n\n", sep="", append=TRUE, file=Textfile)
  
  cat("Provide details for the parameters included in 'cpars' here instead of in the corresponding slot sections below.\n", sep="", append=TRUE, file=Textfile)
  cat("Text in the slot section below will be ignored if a parameter is included in 'cpars'.\n", sep="", append=TRUE, file=Textfile)
  
  cat("Delete this section (including heading) if the `cpars` feature is not used in the OM.\n", sep="", append=TRUE, file=Textfile)
  
  # Stock Parameters ####
  cat("\n\n# Stock Parameters\n\n", sep="", append=TRUE, file=Textfile) 
  slots <- slotNames("Stock")
  for (X in slots) {
    cat("## ", X, "\n", sep="", append=TRUE, file=Textfile)
    if (is.null(ObTemplates$StockTemp)) {
      cat("Justification/explanation for parameter values \n\n", sep="", append=TRUE, file=Textfile)  
    } else {
      cat("Borrowed from ", ObTemplates$StockTemp@Name, "\n\n", sep="", append=TRUE, file=Textfile)
    }
  }
  
  # Fleet Parameters ####
  cat("\n\n# Fleet Parameters\n\n", sep="", append=TRUE, file=Textfile) 
  slots <- slotNames("Fleet")
  for (X in slots) {
    cat("## ", X, "\n", sep="", append=TRUE, file=Textfile)
    if (is.null(ObTemplates$FleetTemp)) {
      cat("Justification/explanation for parameter values \n\n", sep="", append=TRUE, file=Textfile)  
    } else {
      cat("Borrowed from ", ObTemplates$FleetTemp@Name, "\n\n", sep="", append=TRUE, file=Textfile)
    }
  }
  
  
  # Obs Parameters ####
  cat("\n\n# Obs Parameters\n\n", sep="", append=TRUE, file=Textfile) 
  slots <- slotNames("Obs")
  for (X in slots) {
    cat("## ", X, "\n", sep="", append=TRUE, file=Textfile)
    if (is.null(ObTemplates$ObsTemp)) {
      cat("Justification/explanation for parameter values \n\n", sep="", append=TRUE, file=Textfile)  
    } else {
      cat("Borrowed from ", ObTemplates$ObsTemp@Name, "\n\n", sep="", append=TRUE, file=Textfile)
    }
  }
  
  # Imp Parameters ####
  cat("\n\n# Imp Parameters\n\n", sep="", append=TRUE, file=Textfile) 
  slots <- slotNames("Imp")
  for (X in slots) {
    cat("## ", X, "\n", sep="", append=TRUE, file=Textfile)
    if (is.null(ObTemplates$ImpTemp)) {
      cat("Justification/explanation for parameter values \n\n", sep="", append=TRUE, file=Textfile)  
    } else {
      cat("Borrowed from ", ObTemplates$ImpTemp@Name, "\n\n", sep="", append=TRUE, file=Textfile)
    }
  }
  
  # References ####
  cat("\n\n# References\n\n", sep="", append=TRUE, file=Textfile) 
  
  
  message("\nOM spreadsheet and text file initialized\n")
  message("Populate OM parameters in ", name)
  message("Document OM parameters in ", Textfile)
  
 
  
}



ObjTemps <- function(templates=NULL) {
  StockTemp <- NULL; FleetTemp <- NULL; ObsTemp <- NULL; ImpTemp <- NULL
  
  if (class(templates) != "NULL" && class(templates) != "list") stop("'templates' must be a named list")
  if (class(templates) == 'list') {
    if (is.null(names(templates))) stop("'templates' must be a named list")
    if (any(!names(templates) %in% c("Stock", "Fleet", "Obs", "Imp"))) 
      stop("invalid names in list 'templates'. Must be one or all of: Stock, Fleet, Obs, Imp")
    ind <- match("Stock", names(templates))
    if (!is.na(ind)) {
      StockTemp <- templates$Stock
      if (class(StockTemp)=='character' && !StockTemp %in% avail("Stock")) stop(StockTemp, " is not an available object of class Stock")
      if (class(StockTemp)!='character' && class(StockTemp) != 'Stock') stop(StockTemp, " is not an available object of class Stock")
      if (class(StockTemp) == 'character') StockTemp <- get(StockTemp)
    }
    ind <- match("Fleet", names(templates))
    if (!is.na(ind)) {
      FleetTemp <- templates$Fleet
      if (class(FleetTemp)=='character' && !FleetTemp %in% avail("Fleet")) stop(FleetTemp, " is not an available object of class Fleet")
      if (class(FleetTemp)!='character' && class(FleetTemp) != 'Fleet') stop(FleetTemp, " is not an available object of class Fleet")
      if (class(FleetTemp) == 'character') FleetTemp <- get(FleetTemp)
    }
    ind <- match("Obs", names(templates))
    if (!is.na(ind)) {
      ObsTemp <- templates$Obs
      if (class(ObsTemp)=='character' && !ObsTemp %in% avail("Obs")) stop(ObsTemp, " is not an available object of class Obs")
      if (class(ObsTemp)!='character' && class(ObsTemp) != 'Obs') stop(ObsTemp, " is not an available object of class Obs")
      if (class(ObsTemp) == 'character') ObsTemp <- get(ObsTemp)
    }
    ind <- match("Imp", names(templates))
    if (!is.na(ind)) {
      ImpTemp <- templates$Imp
      if (class(ImpTemp)=='character' && !ImpTemp %in% avail("Imp")) stop(ImpTemp, " is not an available object of class Imp")
      if (class(ImpTemp)!='character' && class(ImpTemp) != 'Imp') stop(ImpTemp, " is not an available object of class Imp")
      if (class(ImpTemp) == 'character') ImpTemp <- get(ImpTemp)
    }
  }
  return(list=c(StockTemp=StockTemp, FleetTemp=FleetTemp, ObsTemp=ObsTemp, ImpTemp=ImpTemp))
}



#' Load OM from Excel file
#' 
#' Imports an OM from a correctly formatted Excel file. Create the Excel spreadsheet template
#' using `OMinit` and document each slot in the corresponding text file.
#' 
#' An error message will alert if any slots are missing values, or if the Excel file is missing
#' the required tabs.
#'
#' @param name Name of the OM Excel file in the current working directory.
#' @param cpars An optional list of custom parameters (single parameters are a vector nsim 
#' long, time series are a matrix nsim rows by nyears columns)
#'
#' @return An object of class OM
#' @export
#' @author A. Hordyk
#'
#' @examples 
#' \dontrun{
#' OMinit('myOM', templates=list(Stock='Herring', Fleet='Generic_fleet', Obs='Generic_obs',
#' Imp='Perfect_Imp'), overwrite=TRUE)
#' myOM <- XL2OM('myOM.xlsx')
#' 
#' }
XL2OM <- function(name=NULL, cpars=NULL) {
  
  if (class(name) != 'character') stop("file name must be provided")
  # Load the Excel File ####
  if (nchar(tools::file_ext(name)) == 0) {
    xl.fname1 <- paste0(name, ".xlsx")
    xl.fname2 <- paste0(name, ".xls")
    fls <- file.exists(c(xl.fname1, xl.fname2))
    if (sum(fls) == 0) stop('file ',  xl.fname1, " or ", xl.fname2, " not found")
    if (sum(fls) > 1) stop(name, " found with multiple extensions. Specify file extension.")
    name <- c(xl.fname1, xl.fname2)[fls]
  }
  if (!file.exists(name)) stop('file ', name, " not found") 
  message("Reading ", name)
  sheetnames <- readxl::excel_sheets(name)  # names of the sheets 
  reqnames <- c("OM", "Stock", "Fleet", "Obs", "Imp") 
  ind <- which(!reqnames%in% sheetnames)
  if (length(ind)>0) {
    message("Sheets: ", paste(reqnames[ind], ""), "not found in ", name)
  }
  count <- 1
  tempObj <- vector("list", 4)
  for (obj in c("Stock", "Fleet", "Obs", "Imp")) {
    sht <- as.data.frame(readxl::read_excel(name, sheet = obj, col_names = FALSE))
    rows <- sht[,1] 
    rows <- rows[!rows == "Slot"]
    ind <- which(!rows %in% slotNames(obj))
    if (length(ind)>0) {
      warning(paste(rows[ind], ""), "are not valid slots in object class ", obj)
    }
    
    if (all(dim(sht) == 0)) stop("Nothing found in sheet: ", obj)
    tmpfile <- tempfile(fileext=".csv")
    writeCSV2(inobj = sht, tmpfile, objtype = obj)
    if (ncol(sht)<2) {
      unlink(tmpfile)
      stop("No parameter values found in Sheet ", obj)
    } else {
      tempObj[[count]] <- new(obj, tmpfile)  
    }
    unlink(tmpfile)
    count <- count + 1
  }
  
  # Operating Model
  OM <- new("OM", Stock = tempObj[[1]], Fleet = tempObj[[2]], 
            Obs = tempObj[[3]], Imp=tempObj[[4]])
  
  # Read in the OM sheet
  sht <- as.data.frame(readxl::read_excel(name, sheet = "OM", col_names = FALSE))
  dat <- sht[,1:2] 
  dat <- dat[which(dat[,1] != "Slot"),]
  if (ncol(sht)>2) warning("More than two columns found in Sheet OM. Values in columns C+ are ignored")
  if (ncol(sht)<2) {
    message("No values found for OM slots in Sheet OM. Using defaults")
  } else {
    for (xx in 1:nrow(dat)) {
      val <- dat[xx, 2]
      if (length(val)) {
        if (dat[xx,1] != "Name") slot(OM, dat[xx, 1]) <- as.numeric(val)
        if (dat[xx,1] == "Name") slot(OM, dat[xx, 1]) <- val
      } else{
        message("No value found for OM slot ", dat[xx,1], ". Using default: ", slot(OM, dat[xx, 1]))
      }
    }
  }
  
  if (!is.null(cpars)) {
    if (class(cpars) == "list") {
      OM@cpars <- cpars
    } else {
      stop("'cpars' must be a list")
    }
  }
  ChkObj(OM)
  message('OM successfully imported\n')
  message("Document OM slots in text file (probably ", tools::file_path_sans_ext(name), ".txt),
  and run 'OMdoc' if OM parameter values have changed." )
  OM
}

writeCSV2 <- function(inobj, tmpfile = NULL, objtype = c("Stock", "Fleet", 
                                                         "Obs", "Imp", "Data", "OM", "Fease")) {
  objtype <- match.arg(objtype)
  
  for (X in 1:nrow(inobj)) {
    indat <- inobj[X, ]
    index <- which(!is.na(indat))
    if (length(index) >1) {
      index <- 2:max(index)
      if (X == 1) 
        write(do.call(paste, c(indat[1], as.list(indat[index]), sep = ",")), tmpfile, 1)
      if (X > 1) 
        write(do.call(paste, c(indat[1], as.list(indat[index]), sep = ",")), tmpfile, 1, append = TRUE)    
    } else if (indat[1] != "Slot") {
      write(unlist(indat[1]), tmpfile, 1, append = TRUE)  
    }
    
  }
}




#' Generate OM Documentation Report
#'
#' @param OM An object of class 'OM'
#' @param txt.file Optional. Name of the text file corresponding to the 'OM'. Default assumption
#' is that the txt file is 'OM@Name.txt'
#' @param overwrite Logical. Should existing files be overwritten?
#' @param toc Logical. Should a Table of Contents be included?
#' @param color Character. Color of the OM values in the text. 
#' @param out.file Optional. Character. Name of the output file. Default is the same as the text file.
#' @param output Character. Format of the compiled output document. Default is 'html_document'. Other 
#' options may require additional software and settings.
#' @param theme Optional. Theme for the html output file.
#' @param inc.plot Logical. Should the plots be included?
#' @param render Logical. Should the document be compiled? May be useful if there are problems 
#' with compililing the Rmd file.
#'
#' @return Creates a Rmarkdown file and compiles a HTML report file in the working directory.
#' @export
#' @importFrom methods getSlots
#' @importFrom knitr kable
#' @author A. Hordyk
#' @examples 
#' \dontrun{
#' OMinit('myOM', templates=list(Stock='Herring', Fleet='Generic_fleet', Obs='Generic_obs',
#' Imp='Perfect_Imp'), overwrite=TRUE)
#' myOM <- XL2OM('myOM.xlsx')
#' OMdoc(myOM)
#' }
OMdoc <- function(OM, txt.file=NULL, overwrite=FALSE,
                  toc=TRUE, color="blue",
                  out.file=NULL, output="html_document", theme="flatly", 
                  inc.plot=TRUE, render=TRUE) {
  
  if (class(OM) != "OM") stop('OM must be class "OM"')
  
  ## Read in text file ####
  if (is.null(txt.file)) {
    txt.file <- list.files(pattern=".txt")
    if (length(txt.file) == 0) stop("txt.file' not specified and no .txt files found in working directory")
    if (length(txt.file) == 1) {
      message("txt.file not specified. Reading ", txt.file, " found in working directory")
      textIn <- readLines(txt.file)
    } else {
      NoExt <- tools::file_path_sans_ext(txt.file)
      ind <- which(tolower(NoExt) == tolower(OM@Name))
      if (length(ind) > 0) {
        txt.file <- txt.file[ind]
        message("Reading ", txt.file)
        textIn <- readLines(txt.file)
      } else {
        stop("'txt.file' not specified and multiple .txt files found in working directory")
      }
    }
  } else {
    if (nchar(tools::file_ext(txt.file)) == 0) {
      txt.file <- paste0(txt.file, ".txt")
    } else if (tools::file_ext(txt.file) != "txt") stop("txt.file extension must be txt")
    
    if (!file.exists(txt.file)) stop(txt.file, " not found in working directory")
    message("Reading ", txt.file)
    textIn <- readLines(txt.file)
  }
  
  ## Create Markdown file ####
  if (is.null(out.file)) out.file <- tools::file_path_sans_ext(txt.file)
  RMDfile <- paste0(out.file, ".Rmd")
  if (file.exists(RMDfile) & !overwrite) {
    stop(RMDfile, " already exists.\n Provide a different output file name ('out.file') or use 'overwrite=TRUE'")
  } else {
    message('Writing ', RMDfile)
    tt <- file.create(RMDfile)
  }
  
  ## Write YAML ####
  ind <- grep("^# Title", textIn)
  if (length(ind)>0) {
    title <- trimws(textIn[ind+1])
    if (nchar(title) == 0) title <- paste("Operating Model:", OM@Name)
  } else {
    title <- paste("Operating Model:", OM@Name)
  }
  
  ind <- grep("^# Subtitle", textIn)
  if (length(ind)>0) {
    subtitle <- trimws(textIn[ind+1])
    if (nchar(subtitle) == 0) subtitle <- NULL
  } else {
    subtitle <- NULL
  }
  
  ind <- grep("# Author", textIn)
  if (length(ind)>0) {
    temp <- min(which(nchar(textIn[(ind+1):length(textIn)]) == 0))
    if (temp > 1) {
      authors <- trimws(textIn[(ind+1):(ind+temp-1)])
    } else {
      authors <- "No author provided"
    }
  } else {
    authors <- "No author provided"
  }
  
  # ind <- grep("# Affiliation/Email", textIn)
  # if (length(ind)>0) {
  #   temp <- min(which(nchar(textIn[(ind+1):length(textIn)]) == 0))
  #   if (temp > 1) {
  #     affiliation <- trimws(textIn[(ind+1):(ind+temp-1)])
  #   } else {
  #     affiliation <- NULL
  #   }
  # } else {
  #   affiliation <- NULL
  # }
  # if (length(affiliation) >0) {
  #   affiliation <- rep(affiliation, length(authors))[1:length(authors)]
  #   affiliation[is.na(affiliation)] <- NULL
  # }
  # if (length(affiliation) > length(authors)) {
  #   warning("Should not be more lines in affiliation than Authors")
  #   affiliation <- affiliation[1:length(authors)]
  # }
  
  ind <- grep("# Date", textIn)
  if (length(ind)>0) {
    date <- trimws(textIn[(ind+1)])
    if (grepl("Optional. Date", date)) date <- NULL
  } else {
    date <- NULL
  }
  
  
  cat("---\n", sep="", append=TRUE, file=RMDfile)
  cat("title: '", title, "'\n", append=TRUE, file=RMDfile, sep="") 
  if (!is.null(subtitle)) cat("subtitle: '", subtitle, "'\n", append=TRUE, file=RMDfile, sep="")
  if (length(authors) > 1) {
    cat("author:", "\n", append=TRUE, file=RMDfile, sep="")
    for (xx in 1:length(authors)) {
      # if (!any(is.null(affiliation)) && affiliation[xx] != 'NA') {
      #   cat("- ", authors[xx], "^[", affiliation[xx], "]\n", append=TRUE, file=RMDfile, sep="")  
      # } else {
      #   cat("- ", authors[xx], "\n", append=TRUE, file=RMDfile, sep="")  
      # }
      cat("- ", authors[xx], "\n", append=TRUE, file=RMDfile, sep="")
      
    }
  } else {
    cat("author: ", authors, "\n", append=TRUE, file=RMDfile, sep="")
    # if (is.null(affiliation)) {
    #   cat("author: ", authors, "\n", append=TRUE, file=RMDfile, sep="")
    # } else {
    #   cat("author: ", authors, "^[", affiliation, "]\n", append=TRUE, file=RMDfile, sep="")
    # }
    
  }
  if (is.null(date)) date <- format(Sys.time(), '%d %B %Y')
  cat("date: ", date, "\n", append=TRUE, file=RMDfile, sep="")
  
  if (toc) {
    cat("output: ", "\n", append=TRUE, file=RMDfile, sep="")  
    cat("   ", output, ":", "\n", append=TRUE, file=RMDfile, sep="")  
    cat("     toc: true\n", append=TRUE, file=RMDfile, sep="")  
    cat("     toc_depth: 3\n", append=TRUE, file=RMDfile, sep="")  
    cat("     toc_float: true\n", append=TRUE, file=RMDfile, sep="")
    cat("     theme: ", theme, "\n", append=TRUE, file=RMDfile, sep="")
    
  } else {
    cat("output: ", output, "\n", append=TRUE, file=RMDfile, sep="")
  }
  cat("---\n\n", sep="", append=TRUE, file=RMDfile) 
  
  
  ## OM Details ####
  cat("# Operating Model Details \n", append=TRUE, file=RMDfile, sep="")
  cat("**Name**: '", OM@Name, "'\n\n", append=TRUE, file=RMDfile, sep="")
  cat("**nsim**: ", OM@nsim, " | ", append=TRUE, file=RMDfile, sep="")
  cat(" **proyears**: ", OM@proyears, " | ", append=TRUE, file=RMDfile, sep="")
  cat(" **interval**: ", OM@interval, "\n\n", append=TRUE, file=RMDfile, sep="")
  cat("**pstar**: ", OM@pstar, " | ", append=TRUE, file=RMDfile, sep="")
  cat(" **maxF**: ", OM@maxF, " | ", append=TRUE, file=RMDfile, sep="")
  cat(" **reps**: ", OM@reps, " \n\n", append=TRUE, file=RMDfile, sep="")
  cat("**Source**: '",  OM@Source, "' \n\n", append=TRUE, file=RMDfile, sep="")
  useCpars <- length(OM@cpars) > 0
  if (useCpars) cat("**Uses cpars**: ", useCpars, " \n\n", append=TRUE, file=RMDfile, sep="")
  
  ## Generate Sampled Parameters ####
  
  Pars <- NULL
  if (inc.plot) {
    
    # --- Generate Historical Samples ----
    message("\n\nRunning Historical Simulations\n\n")
    OM <- updateMSE(OM) # update and add missing slots with default values
    if (OM@nsim > 48) setup()
    out<- runMSE(OM,Hist=T)
    
    Pars <- c(out$SampPars, out$TSdata, out$MSYs)
  
  }
  
  
  ## Input text ####
  ind <- which(unlist(lapply(textIn, nchar)) == 0) # delete empty lines 
  if (length(ind) > 0) textIn <- textIn[-ind]
  ind <- grep("# Introduction", textIn)
  if (length(ind)>1) stop("# Introduction heading found more than once")
  if (length(ind)>0) {
    textIn <- textIn[ind:length(textIn)]
  } else {
    ind <- grep("# Stock Parameters", textIn)
    if (length(ind)>1) stop("# Stock Parameters heading found more than once")
    if (length(ind) == 0) stop("# Stock Parameters not found")
    textIn <- textIn[ind:length(textIn)]
  }
  
  
  ## Introduction ####
  writeSection(class="Intro", OM, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## Cpars ####
  if (useCpars) {
    writeSection(class="cpars", OM, textIn, RMDfile, color=color, inc.plot=inc.plot)
  }
  
  ## Stock Parameters ####
  writeSection(class="Stock", OM, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## Fleet Parameters ####
  writeSection(class="Fleet", OM, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## Observation Parameters ####
  writeSection(class="Obs", OM, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## Implementation Parameters ####
  writeSection(class="Imp", OM, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## OM Plots ####
  if (inc.plot) {
    cat("# OM Plots\n\n", sep="", append=TRUE, file=RMDfile) # write heading
    cat("```{r, echo=FALSE, fig.asp=2}\n", append=TRUE, file=RMDfile, sep="")
    cat("plot.OM(out)\n", append=TRUE, file=RMDfile, sep="")
    cat("```\n\n\n", append=TRUE, file=RMDfile, sep="")

  }
  
  
  
  
  
  ## References ####
  writeSection(class="References", OM, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## Render Markdown ####
  if (render) {
    message("\n\nRendering markdown document as ", RMDfile)
    EffYears <- seq(from=(OM@CurrentYr -  OM@nyears + 1), to=OM@CurrentYr, length.out=length(OM@EffYears))
    EffYears <- round(EffYears,0)
    Effvals <- data.frame(EffYears=EffYears, EffLower=OM@EffLower, EffUpper=OM@EffUpper)
    params <- list(OM=OM, Pars=Pars, Effvals=Effvals, out=out)
    rmarkdown::render(input=RMDfile, output_format=output, param=params)
    
    utils::browseURL(file.path(getwd(), paste0(out.file, ".html")))
    
  } else {
    
  }
  
}




# Text templates for the OM documentation ####
Template <- function(type=c("Stock", "Fleet", "Obs", "Imp")) {
  type <- match.arg(type)
  if (type == "Stock") mat <- 
      matrix(c("Mortality and age:  maxage, R0, M, M2, Mexp, Msd, Mgrad",
               "Recruitment: h, SRrel, Perr, AC, recgrad",
               "Non-stationarity in stock productivity: Period, Amplitude",
               "Growth: Linf, K, t0, LenCV, Ksd, Kgrad, Linfsd, Linfgrad",
               "Maturity: L50, L50_95",
               "Stock depletion: D",
               "Length-weight conversion parameters: a, b",
               "Spatial distribution and movement: Frac_area_1, Prob_staying",
               "Discard Mortality: Fdisc "), ncol=1)
  if (type == "Fleet") mat <- 
      matrix(c(
        "Historical years of fishing, spatial targeting: nyears, Spat_targ",
        "Trend in historical fishing effort (exploitation rate), interannual variability in fishing effort: EffYears, EffLower, EffUpper, Esd",
        "Annual increase in catchability, interannual variability in catchability: qinc, qcv",
        "Fishery gear length selectivity: L5, LFS, Vmaxlen, isRel",
        "Fishery length retention: LR5, LFR, Rmaxlen, DR",
        "Time-varying selectivity: SelYears, AbsSelYears, L5Lower, L5Upper, LFSLower, LFSUpper, VmaxLower, VmaxUpper",
        "Current Year: CurrentYr"), ncol=1)
  if (type == "Obs") mat <- 
      matrix(c(
        "Catch statistics: Cobs, Cbiascv, CAA_nsamp, CAA_ESS, CAL_nsamp, CAL_ESS, CALcv",
        "Index imprecision, bias and hyperstability: Iobs, Icv, Btcv, Btbias, beta",
        "Bias in maturity, natural mortality rate and growth parameters: LenMcv, Mcv, Kcv, t0cv, Linfcv",
        "Bias in length at first capture, length at full selection: LFCcv, LFScv",
        "Bias in fishery reference points, unfished biomass, FMSY, FMSY/M ratio, biomass at MSY relative to unfished: FMSYcv, FMSY_Mcv, BMSY_B0cv",
        "Management targets in terms of the index (i.e., model free), the total annual catches and absolute biomass levels: Irefcv, Crefcv, Brefcv",
        "Depletion bias and imprecision: Dbiascv, Dcv",
        "Recruitment compensation and trend: hcv, Reccv",
        "Currently unused observation processes - bias in unfished biomass, intrinsic rate of increase, annual increase in fishing efficiency and age at 50% vulnerability, bias and imprecision in current fishing rate, bias in maximum age: B0cv, rcv, Fcurbiascv, Fcurcv, maxagecv"), ncol=1)
  if (type == "Imp") mat <-
      matrix(c(
        "Output Control Implementation Error: TACFrac, TACSD",
        "Effort Control Implementation Error: EFrac, ESD", 
        "Size Limit Control Implementation Error: SizeLimFrac, SizeLimSD"), ncol=1)
  
  # Check slots 
  Slots <- names(methods::getSlots(type))
  for (X in Slots) {
    tt <- grep(paste0("\\<", X, "\\>"), mat) 
    if (X != "Name" && X != "Source") {
      if (length(tt) < 1) stop("slot ", X, " not found in ", type, " template")
      if (length(tt) > 1) stop("slot ", X, " found multiple times in ", type, " template")    
    }
  }
  return(mat)
}


writeSection <- function(class=c("Intro", "Stock", "Fleet", "Obs", "Imp", "References", "cpars"), OM,
                         textIn, RMDfile, color, inc.descript=TRUE, inc.plot=TRUE) {
  class <- match.arg(class)
  
  useCpars <- length(OM@cpars) > 0
  
  
  fLH <- grep("^#[^##]", textIn)
  fstH <- trimws(gsub("#", "", textIn[fLH])) # first level headings
  fstHmd <- trimws(textIn[fLH]) # first level headings
  
  if (class == "Intro") {
    intro <- which(trimws(gsub("#", "", fstH)) == "Introduction")
    if (intro == 0) stop("# Introduction heading not found")
    if (length(intro) == 1) {
      cat(fstHmd[intro], "\n\n", sep="", append=TRUE, file=RMDfile) # write first heading
      for (ll in (intro+1):(fLH[intro+1] - 1)) {
        cat(textIn[ll], "\n\n", sep="", append=TRUE, file=RMDfile) # write intro paragraphs
      }
    } else {
      stop("More than one section # Introduction")
    }
    
  } else if (class == "References") {
    refText <- which(trimws(gsub("#", "", fstH)) == "References")
    if (length(refText) == 1) {
      cat("# References\n\n", sep="", append=TRUE, file=RMDfile) # write heading
      st <- which(textIn == '# References')
      end <- length(textIn)
      if (st+1 < end) {
        for (ll in (st+1):(end)) {
          cat(textIn[ll], "\n\n", sep="", append=TRUE, file=RMDfile) # 
        }   
      }
      
    } else stop("More than one section # References")
  } else if (class == "cpars") {
    cparstext <- which(trimws(gsub("#", "", fstH)) == "Custom Parameters")
    if (length(cparstext) == 1) {
      # get cpars text 
      cat("# Custom Parameters\n\n", sep="", append=TRUE, file=RMDfile) # write heading
      
      st <- which(textIn == '# Custom Parameters')
      end <- textIn %in% paste0('# ', fstH)
      temp <- which(textIn[end] != "# Custom Parameters") 
      temp <- min(temp[temp>cparstext])
      end <- which(textIn == paste0("# ", fstH[temp]))
      for (ll in (st+1):(end - 1)) {
        cat(textIn[ll], "\n\n", sep="", append=TRUE, file=RMDfile) # write cpars paragraphs
      }
    }  else stop("More than one section # Custom Parameters")
    
  } else {
    # Write class heading 
    st <- which(trimws(gsub("#", "", textIn)) == paste(class, "Parameters"))
    sta <- which(fstH == paste(class, "Parameters"))
    if (length(st) > 1) stop("Multiple '# ", class, " Parameters' headings in document.")
    if (length(st) < 1) stop("'# ", class, " Parameters' heading not found in document.")
    cat("# ", class, " Parameters \n\n", append=TRUE, file=RMDfile, sep="")
    
    # Find second level headings and check that they match slots in class
    bg <- st+1
    end <- fLH[sta+1]-1
    if (is.na(end)) end <- length(textIn)
    Text <- textIn[bg:end]
    scLHloc <- grep("^##[^###]", Text) # second level headings
    scLHmd <- Text[scLHloc]
    scLH <- trimws(gsub("##", "",scLHmd))
    
    Slots <- slotNames(class)
    if (any(!scLH %in% Slots)) {
      invalid <- scLH[!scLH %in% Slots]
      stop("Invalid second level headings (must match slots in class ", class, "): ", paste(invalid, ""))
    }
    
    # Get template for class section 
    ClTemp <- Template(class)
    
    # loop through template and write section 
    for (rr in 1:nrow(ClTemp)) {
      if(grepl("^Currently unused", ClTemp[rr,1])) {
        temptext <- trimws(unlist(strsplit(ClTemp[rr,], "-")))
        cat("### ", temptext[1], "\n\n", append=TRUE, file=RMDfile, sep="")
        cat("*", temptext[2], "*\n\n", append=TRUE, file=RMDfile, sep="")
      } else {
        slots <- trimws(unlist(strsplit(strsplit(ClTemp[rr,], ":")[[1]][2], ",")))
        cat("### ", ClTemp[rr,], "\n\n", append=TRUE, file=RMDfile, sep="")
        for (sl in slots) {
          # get slot value if not in cpars 
          if (useCpars && sl %in% names(OM@cpars)) {
            val <- range(OM@cpars[[sl]])
            used <- TRUE
            val <- gsub('"', "", paste(val, collapse="\", \""))
            valtext <- paste0("(", trimws(val), ")")
            # currently not used. Manually describe in the custom parameters section
            
          } else {
            val <- slot(OM, sl)
            if (is.numeric(val)) val <- round(val,2)
            used <- length(val)>0 && !is.na(val) && !is.null(val) # is the slot used?
            if (used) {
              val <- gsub('"', "", paste(val, collapse="\", \""))
              valtext <- paste0("(", trimws(val), ")")
            } else {
              valtext <- val
            }           
          }
          
          valtext <- paste0("<span style='color:", color, "'>", " ", valtext, "</span>\n\n")
          
          loc <- which(scLH == sl)
          if (length(loc) > 0) {
            bg <- scLHloc[loc]+1
            end <- scLHloc[loc+1]-1
            if (is.na(end)) end <- length(Text)
            description <- Text[bg:end]
            if (!used) description <- c("Slot not used.")
            if (used && ! sl%in% c("EffYears", "EffLower", "EffUpper")) description <- c(description, valtext)

          } else {
            if (used & sl != "Source") description <- paste("No justification provided.", valtext)
            if (!used) description <- "Slot not used. "
          }

          
          if (inc.descript) {
            des <- switch(class, 
                          "Stock" = DLMtool::StockDescription,
                          "Fleet" = DLMtool::FleetDescription,
                          "Obs" = DLMtool::ObsDescription,
                          "Imp" = DLMtool::ImpDescription)
            
            rownames(des) <- des[,1]
            cat("**", sl, "**: ", des[sl, 2], append=TRUE, file=RMDfile, sep="")
            cat("\n\n", description,"\n\n", append=TRUE, file=RMDfile, sep="")
            
          } else {
            cat("**", sl, "**:\n",  description,"\n\n", append=TRUE, file=RMDfile, sep="")
          }
          
          if (used && sl == "EffUpper") {
            cat("<style type='text/css'>\n", append=TRUE, file=RMDfile, sep="")
            cat(".table {\n", append=TRUE, file=RMDfile, sep="")
            cat("    width: 75%; \n", append=TRUE, file=RMDfile, sep="")
            cat("}\n", append=TRUE, file=RMDfile, sep="")
            cat("</style>\n", append=TRUE, file=RMDfile, sep="")
            
            cat("```{r, echo=FALSE, results='asis'}\n", append=TRUE, file=RMDfile, sep="")
            cat("knitr::kable(round(Effvals,2), format='markdown', caption='')\n", append=TRUE, file=RMDfile, sep="")
            cat("```\n\n", append=TRUE, file=RMDfile, sep="")
            
          }
          
          
          # Plots ####
          if (inc.plot) {
            if (class %in% c("Stock", "Fleet")) {
              if (sl == slots[length(slots)]) plotText(OM, slots, RMDfile)
            } 
        
          }
        }
      }
    }
    if (class == "Obs" | class =='Imp') {
      plotText(OM, slots=class, RMDfile)
    }
  }
}


