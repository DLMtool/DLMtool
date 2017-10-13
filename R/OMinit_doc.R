#' Copy example OM XL and OM Documentation 
#'
#' @export
#'
#' @examples
#' \dontrun{
#' OMexample()
#' }
OMexample <- function() {
  fromRMD <- system.file("Example_Chile_hake_source.Rmd", package="DLMtool")
  tt <- file.copy(fromRMD, getwd(), overwrite = TRUE)
  fromXL <- system.file("Example_Chile_hake.xlsx", package="DLMtool")
  tt <- file.copy(fromXL, getwd(), overwrite = TRUE)
}

#' Initialize Operating Model
#'
#' Generates an Excel spreadsheet and a source.rmd file in the current working directory for 
#' specifying and documenting a DLMtool Operating Model.
#' 
#' @param name The name of the Excel and source.rmd file to be created in the working directory (character). 
#' Use 'example' for a populated example OM XL and documentation file.
#' @param files What files should be created: 'xlsx', 'rmd', or c('xlsx', 'rmd') (default: both)
#' @param templates An optional named list of existing DLMtool Stock, Fleet, Obs, or Imp objects
#' to use as templates for the Operating Model.
#' @param overwrite Logical. Should files be overwritten if they already exist?
#'
#' @return OM xlsx and source.rmd files are created in the working directory.  
#' @export
#' @author A. Hordyk
#'
#' @examples
#' \dontrun{
#' # Create an Excel OM template and rmd file called 'myOM.xlsx' and 'myOM_source.rmd': 
#' OMinit('myOM')
#' 
#' # Create an Excel OM template and text file called 'myOM.rmd' and 'myOM_source.rmd', using
#' the Stock object 'Herring' as a template: 
#' OMinit('myOM', templates=list(Stock='Herring'))
#' 
#' # Create an Excel OM template and text file called 'myOM.rmd' and 'myOM_source.rmd', using
#' the Stock object 'Herring', and Obs object 'Generic_obs' as templates: 
#' OMinit('myOM', templates=list(Stock='Herring', Obs='Generic_obs'), overwrite=TRUE)
#' }
#' 
OMinit <- function(name=NULL, files=c('xlsx', 'rmd'), templates=NULL, overwrite=FALSE) {
  files <- match.arg(files, several.ok = TRUE)
  
  if (is.null(name)) stop("Require OM name", call.=FALSE)
  
  if (tolower(name) == 'example') {
    OMexample()
    return(message("Creating Example Files in ", getwd()))
  }
  if (class(name) != 'character') stop("name must be text", call.=FALSE)
 
  ## Write Excel skeleton ####
  if (nchar(tools::file_ext(name)) == 0) {
    nameNoExt <- name
    name <- paste0(name, ".xlsx")
  } else {
    ext <- tools::file_ext(name)
    if (!ext %in% c("xlsx", "xls")) stop("File extension must be 'xlsx' or 'xls'", call.=FALSE)
    nameNoExt <- tools::file_path_sans_ext(name)
  }
  
  
  ObTemplates <- ObjTemps(templates)
  if (!is.null(ObTemplates) && class(ObTemplates) == 'list') {
    nm <- names(ObTemplates)  
    message("\n\nUsing Object Templates:")
    for (X in nm) {
      message(ObTemplates[[X]]@Name)
    }
  }
  
  
  
  if ('xlsx' %in% files) {
    wb <- openxlsx::createWorkbook()
    
    openxlsx::addWorksheet(wb, sheetName = "Stock", gridLines = TRUE)
    openxlsx::addWorksheet(wb, sheetName = "Fleet", gridLines = TRUE)
    openxlsx::addWorksheet(wb, sheetName = "Obs", gridLines = TRUE)
    openxlsx::addWorksheet(wb, sheetName = "Imp", gridLines = TRUE)
    openxlsx::addWorksheet(wb, sheetName = "OM", gridLines = TRUE)
    openxlsx::addWorksheet(wb, sheetName = "Data", gridLines = TRUE)
    
    

    # Load the Excel File ####
    
    
    message("Creating ", name, " in ", getwd())
    
    
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
      stop(name, " already exists.\n Use 'overwrite=TRUE'. \nNOTE: this will overwrite both .xlsx and source.rmd files if they exist.", call.=FALSE)
    } else {
      options(warn=2)
      tryWrite <- try(openxlsx::saveWorkbook(wb, name, overwrite = overwrite), 
                      silent=TRUE) ## save to working directory
      options(warn=0)
      if (tryWrite != 1) stop("Can't write to ", name, ". If the file open?", call.=FALSE)
    }
  
 
  }
 
  if ('rmd' %in% files) {
    ## Write Rmd source skeleton ####
    message("Creating ", nameNoExt, "_source.rmd in ", getwd())
    RmdSource <- paste0(nameNoExt, "_source.rmd")
    if (file.exists(RmdSource) & !overwrite) {
      stop(RmdSource, " already exists.\n Use 'overwrite=TRUE'.", call.=FALSE)
    } else {
      tt <- file.create(RmdSource)
    }
    
    # Title
    cat("\n# Title\n", sep="", append=TRUE, file=RmdSource) 
    cat("Title. One line only.\n", sep="", append=TRUE, file=RmdSource) 
    
    # Subtitle - optional 
    cat("\n# Subtitle\n", sep="", append=TRUE, file=RmdSource) 
    cat("Optional. Subtitle. One line only. Delete text and heading if not required.\n", 
        sep="", append=TRUE, file=RmdSource) 
    
    
    # Author(s) 
    cat("\n# Author(s)\n", sep="", append=TRUE, file=RmdSource) 
    cat("Name and contact details (e.g email, affiliation) for each author.\n", 
        sep="", append=TRUE, file=RmdSource) 
    cat("One line per author.\n", 
        sep="", append=TRUE, file=RmdSource) 
    
    # # Affiliation/Email - optional 
    # cat("\n# Affiliation/Email\n", sep="", append=TRUE, file=RmdSource) 
    # cat("Affiliation and/or email for each author. One line for each author.\n", 
    #     sep="", append=TRUE, file=RmdSource) 
    # cat("Will be recycled for each author if more authors than affiliations.\n", 
    #     sep="", append=TRUE, file=RmdSource) 
    
    # Date - optional 
    cat("\n# Date\n", sep="", append=TRUE, file=RmdSource) 
    cat("Optional. Date that the operating model was created. If none provided, today's date will be used.\n", 
        sep="", append=TRUE, file=RmdSource) 
    
    # Introduction ####
    cat("\n# Introduction\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("## Completing the OM Documentation\n", sep="", append=TRUE, file=RmdSource)
    cat("This document is used to generate a HTML OM report document.\n\n", sep="", append=TRUE, file=RmdSource)
    cat("The document is separated into 7 sections:\n", sep="", append=TRUE, file=RmdSource)
    cat("1. Introduction (this section)\n", sep="", append=TRUE, file=RmdSource)
    cat("2. Custom Parameters (optional)\n", sep="", append=TRUE, file=RmdSource)
    cat("3. Stock Parameters\n", sep="", append=TRUE, file=RmdSource)
    cat("4. Fleet Parameters\n", sep="", append=TRUE, file=RmdSource)
    cat("5. Obs (Observation) Parameters\n", sep="", append=TRUE, file=RmdSource)
    cat("6. Imp (Implementation) Parameters\n", sep="", append=TRUE, file=RmdSource)
    cat("7. References\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("The Introduction section is used to briefly describe the fishery and the details of the Operating Model.\n", sep="", append=TRUE, file=RmdSource)
    cat("It should include an explanation for the OM parameters:\n ", sep="", append=TRUE, file=RmdSource)
    cat("* nsim: the number of simulations.\n ", sep="", append=TRUE, file=RmdSource)
    cat("* proyears: the number of projectio years.\n ", sep="", append=TRUE, file=RmdSource)
    cat("* interval: the management interval.\n ", sep="", append=TRUE, file=RmdSource)
    cat("* pstar: the percentile of the sample of the management recommendation for each method.\n ", sep="", append=TRUE, file=RmdSource)
    cat("* maxF: the maximum instantaneous fishing mortality rate that may be simulated for any given age class.\n ", sep="", append=TRUE, file=RmdSource)
    cat("* reps: the number of samples of the management recommendation for each method.\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("The Custom Parameters section is optional and should only be included if the cpars feature of DLMtool is used in the OM.\n", sep="", append=TRUE, file=RmdSource)
    cat("Delete both the heading and the text in this section if `cpars` are not used.\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("The Stock, Fleet, Obs, and Imp sections include each slot in these components of the OM object.\n", sep="", append=TRUE, file=RmdSource)
    cat("Provide details (including references where applicable) for the choice of values for each slot below the corresponding slot name (e.g., ## M).\n", sep="", append=TRUE, file=RmdSource)
    cat("For example: \n", sep="", append=TRUE, file=RmdSource)
    cat("**M**\n", sep="", append=TRUE, file=RmdSource)
    cat("An explanation for the values of the natural mortality rate in the OM (Smith et al. 1999).\n\n", sep="", append=TRUE, file=RmdSource)
    cat("You do not need to include the actual values from the OM. These will be included automatically in the final compiled document.\n\n", sep="", append=TRUE, file=RmdSource)
    cat("References should be included in the 'References' section at the end of the document.\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("Once complete, this text file will be compiled into an OM Report Document.\n", sep="", append=TRUE, file=RmdSource)
    cat("This text file is linked to the Excel spreadsheet that was generated with the same file name.\n", sep="", append=TRUE, file=RmdSource)
    cat("It serves as a single documentation source for a DLMtool OM, and should be updated whenever parameter values in the OM spreadsheet are updated.\n", sep="", append=TRUE, file=RmdSource)
    
    cat("\n## Tips on filling this Document\n\n", sep="", append=TRUE, file=RmdSource)
    cat("This document is uses the Markdown format. All first and second level headings have been provided, and in general you\n", sep="", append=TRUE, file=RmdSource)
    cat("should only need to enter plain text.\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("You can have multiple paragraphs throughout the document.\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("The Introduction and Custom Parameters sections also support second and third level headings.\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("## An example Second level heading\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("### An example third level heading\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("### Technical Tip\n\n", sep="", append=TRUE, file=RmdSource)
    cat("This document will be compiled into Rmarkdown, and then a HTML document using Pandoc. Equations can be included by\n\n", sep="", append=TRUE, file=RmdSource)
    cat("using Latex (see [here](https://www.sharelatex.com/learn/Mathematical_expressions) for some examples).\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("## Delete all text below 'Introduction' and replace with a description of the OM.\n\n", sep="", append=TRUE, file=RmdSource)
    
    
    # Cpars ####
    cat("\n\n# Custom Parameters\n", sep="", append=TRUE, file=RmdSource)  
    cat("Optional. Only required if the `cpars` feature is used in the OM.\n\n", sep="", append=TRUE, file=RmdSource)
    
    cat("Provide details for the parameters included in 'cpars' here instead of in the corresponding slot sections below.\n", sep="", append=TRUE, file=RmdSource)
    cat("Text in the slot section below will be ignored if a parameter is included in 'cpars'.\n", sep="", append=TRUE, file=RmdSource)
    
    cat("Delete this section (including heading) if the `cpars` feature is not used in the OM.\n", sep="", append=TRUE, file=RmdSource)
    
    # Stock Parameters ####
    cat("\n\n# Stock Parameters\n\n", sep="", append=TRUE, file=RmdSource) 
    slots <- slotNames("Stock")
    for (X in slots) {
      cat("## ", X, "\n", sep="", append=TRUE, file=RmdSource)
      if (is.null(ObTemplates$StockTemp)) {
        cat("No justification provided. \n\n", sep="", append=TRUE, file=RmdSource)  
      } else {
        cat("Borrowed from ", ObTemplates$StockTemp@Name, "\n\n", sep="", append=TRUE, file=RmdSource)
      }
    }
    
    # Fleet Parameters ####
    cat("\n\n# Fleet Parameters\n\n", sep="", append=TRUE, file=RmdSource) 
    slots <- slotNames("Fleet")
    for (X in slots) {
      cat("## ", X, "\n", sep="", append=TRUE, file=RmdSource)
      if (is.null(ObTemplates$FleetTemp)) {
        cat("No justification provided. \n\n", sep="", append=TRUE, file=RmdSource)  
      } else {
        cat("Borrowed from ", ObTemplates$FleetTemp@Name, "\n\n", sep="", append=TRUE, file=RmdSource)
      }
    }
    
    
    # Obs Parameters ####
    cat("\n\n# Obs Parameters\n\n", sep="", append=TRUE, file=RmdSource) 
    slots <- slotNames("Obs")
    for (X in slots) {
      cat("## ", X, "\n", sep="", append=TRUE, file=RmdSource)
      if (is.null(ObTemplates$ObsTemp)) {
        cat("No justification provided. \n\n", sep="", append=TRUE, file=RmdSource)  
      } else {
        cat("Borrowed from ", ObTemplates$ObsTemp@Name, "\n\n", sep="", append=TRUE, file=RmdSource)
      }
    }
    
    # Imp Parameters ####
    cat("\n\n# Imp Parameters\n\n", sep="", append=TRUE, file=RmdSource) 
    slots <- slotNames("Imp")
    for (X in slots) {
      cat("## ", X, "\n", sep="", append=TRUE, file=RmdSource)
      if (is.null(ObTemplates$ImpTemp)) {
        cat("No justification provided. \n\n", sep="", append=TRUE, file=RmdSource)  
      } else {
        cat("Borrowed from ", ObTemplates$ImpTemp@Name, "\n\n", sep="", append=TRUE, file=RmdSource)
      }
    }
    
    # References ####
    cat("\n\n# References\n\n", sep="", append=TRUE, file=RmdSource) 
    
    
    
  } 

  if ('xlsx' %in% files) message("Populate OM parameters in ", name)  
  if ('rmd' %in% files) message("Document OM parameters in ", RmdSource)
}



ObjTemps <- function(templates=NULL) {
  StockTemp <- NULL; FleetTemp <- NULL; ObsTemp <- NULL; ImpTemp <- NULL
  
  if (class(templates) != "NULL" && class(templates) != "list") stop("'templates' must be a named list", call.=FALSE)
  if (class(templates) == 'list') {
    if (is.null(names(templates))) stop("'templates' must be a named list", call.=FALSE)
    if (any(!names(templates) %in% c("Stock", "Fleet", "Obs", "Imp"))) 
      stop("invalid names in list 'templates'. Must be one or all of: Stock, Fleet, Obs, Imp", call.=FALSE)
    ind <- match("Stock", names(templates))
    if (!is.na(ind)) {
      StockTemp <- templates$Stock
      if (class(StockTemp)=='character' && !StockTemp %in% avail("Stock")) stop(StockTemp, " is not an available object of class Stock", call.=FALSE)
      if (class(StockTemp)!='character' && class(StockTemp) != 'Stock') stop(StockTemp, " is not an available object of class Stock", call.=FALSE)
      if (class(StockTemp) == 'character') StockTemp <- get(StockTemp)
    }
    ind <- match("Fleet", names(templates))
    if (!is.na(ind)) {
      FleetTemp <- templates$Fleet
      if (class(FleetTemp)=='character' && !FleetTemp %in% avail("Fleet")) stop(FleetTemp, " is not an available object of class Fleet", call.=FALSE)
      if (class(FleetTemp)!='character' && class(FleetTemp) != 'Fleet') stop(FleetTemp, " is not an available object of class Fleet", call.=FALSE)
      if (class(FleetTemp) == 'character') FleetTemp <- get(FleetTemp)
    }
    ind <- match("Obs", names(templates))
    if (!is.na(ind)) {
      ObsTemp <- templates$Obs
      if (class(ObsTemp)=='character' && !ObsTemp %in% avail("Obs")) stop(ObsTemp, " is not an available object of class Obs", call.=FALSE)
      if (class(ObsTemp)!='character' && class(ObsTemp) != 'Obs') stop(ObsTemp, " is not an available object of class Obs", call.=FALSE)
      if (class(ObsTemp) == 'character') ObsTemp <- get(ObsTemp)
    }
    ind <- match("Imp", names(templates))
    if (!is.na(ind)) {
      ImpTemp <- templates$Imp
      if (class(ImpTemp)=='character' && !ImpTemp %in% avail("Imp")) stop(ImpTemp, " is not an available object of class Imp", call.=FALSE)
      if (class(ImpTemp)!='character' && class(ImpTemp) != 'Imp') stop(ImpTemp, " is not an available object of class Imp", call.=FALSE)
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
#' @param msg Should messages be printed?
#'
#' @return An object of class OM
#' @export
#' @author A. Hordyk
#'
#' @examples 
#' \dontrun{
#' OMinit('myOM', templates=list(Stock='Herring', Fleet='Generic_Fleet', Obs='Generic_Obs',
#' Imp='Perfect_Imp'), overwrite=TRUE)
#' myOM <- XL2OM('myOM.xlsx')
#' 
#' }
XL2OM <- function(name=NULL, cpars=NULL, msg=TRUE) {
  

  # Load the Excel File ####
  if (is.null(name)) {
    fls <- list.files(pattern=".xlsx", ignore.case = TRUE)
    if (length(fls) == 0) stop('Name not provided and no .xlsx files found.', call.=FALSE)
    if (length(fls) > 1) stop("Name not provided and multiple .xlsx files found", call.=FALSE)
    name <- fls
  }
  
  if (class(name) != 'character') stop("file name must be provided", call.=FALSE)

  if (nchar(tools::file_ext(name)) == 0) {
    xl.fname1 <- paste0(name, ".xlsx")
    xl.fname2 <- paste0(name, ".xls")
    fls <- file.exists(c(xl.fname1, xl.fname2))
    if (sum(fls) == 0) stop(xl.fname1, " or ", xl.fname2, " not found")
    if (sum(fls) > 1) stop(name, " found with multiple extensions. Specify file extension.", call.=FALSE)
    name <- c(xl.fname1, xl.fname2)[fls]
  }
  if (!file.exists(name)) stop(name, " not found", call.=FALSE) 
  message("Reading ", name)
  sheetnames <- readxl::excel_sheets(name)  # names of the sheets 
  reqnames <- c("OM", "Stock", "Fleet", "Obs", "Imp") 
  ind <- which(!reqnames%in% sheetnames)
  if (length(ind)>0) stop("Sheets: ", paste(reqnames[ind], ""), "not found in ", name, call.=FALSE)
  
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
    
    if (all(dim(sht) == 0)) stop("Nothing found in sheet: ", obj, call.=FALSE)
    tmpfile <- tempfile(fileext=".csv")
    writeCSV2(inobj = sht, tmpfile, objtype = obj)
    if (ncol(sht)<2) {
      unlink(tmpfile)
      stop("No parameter values found in Sheet ", obj, call.=FALSE)
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
  # if (ncol(sht)>2) warning("More than two columns found in Sheet OM. Values in columns C+ are ignored")
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
      stop("'cpars' must be a list", call.=FALSE)
    }
  }
  ChkObj(OM)
  if (msg) {
    message('OM successfully imported\n')
    message("Document OM slots in source.rmd file (probably ", tools::file_path_sans_ext(name), "_source.rmd),
  and run 'OMdoc' if OM parameter values have changed." )
  }

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
#' @param OM An object of class 'OM' or the name of an OM xlsx file 
#' @param rmd.source Optional. Name of the source.rmd file corresponding to the 'OM'. Default assumption
#' is that the file is 'OM@Name_source.Rmd'
#' @param overwrite Logical. Should existing files be overwritten?
#' @param out.file Optional. Character. Name of the output file. Default is the same as the text file.
#' @param inc.plot Logical. Should the plots be included?
#' @param render Logical. Should the document be compiled? May be useful to turn off if 
#' there are problems with compililing the Rmd file.
#' @param output Character. Output file type. Default is 'html_document'. 'pdf_document' is available
#' but may require additional software and have some formatting issues.
#'
#' @return Creates a Rmarkdown file and compiles a HTML report file in the working directory.
#' @export
#' @importFrom methods getSlots
#' @importFrom knitr kable
#' @author A. Hordyk
#' @examples 
#' \dontrun{
#' OMinit('myOM', templates=list(Stock='Herring', Fleet='Generic_Fleet', Obs='Generic_Obs',
#' Imp='Perfect_Imp'), overwrite=TRUE)
#' myOM <- XL2OM('myOM.xlsx')
#' OMdoc(myOM)
#' }
OMdoc <- function(OM=NULL, rmd.source=NULL, overwrite=FALSE, out.file=NULL,  
                  inc.plot=TRUE, render=TRUE, output="html_document") {
  # markdown compile options
  toc=TRUE; color="blue";  theme="flatly"
  OMXLname <- NULL
  if (class(OM) == "OM") {
    # nothing
  } else if (class(OM) == 'character') {
    OMXLname <- OM
    OM <- XL2OM(OM, msg=FALSE)
  } else if (is.null(OM)) {
    fls <- list.files(pattern=".xlsx", ignore.case=TRUE)
    if (length(fls)==1) OM <- XL2OM(fls, msg=FALSE)
    if (length(fls)>1) stop('argument "OM" not provided and multiple .xlsx files in working directory', call.=FALSE)
  } else stop('OM must be class "OM" or name of OM xlsx file.', call.=FALSE)
  
  if (is.null(OM)) stop('OM not imported. Is the name correct?', call.=FALSE)
  ## Read in Rmd.source file ####
  if (is.null(rmd.source)) {
    rmd.source <- list.files(pattern=".rmd", ignore.case=TRUE)
    if (length(rmd.source) == 0) stop("rmd.source' not specified and no .rmd files found in working directory", call.=FALSE)
    if (length(rmd.source) == 1) {
      # message("rmd.source not specified. Reading ", rmd.source, " found in working directory")
      textIn <- readLines(rmd.source)
    } else {
      NoExt <- tools::file_path_sans_ext(rmd.source)
      if (!is.null(OMXLname)) ind <- which(tolower(NoExt) == tolower(paste0(OMXLname, "_source")))
      if (is.null(OMXLname)) ind <- which(tolower(NoExt) == tolower(paste0(OM@Name, "_source")))
      if (length(ind) > 0) {
        rmd.source <- rmd.source[ind]
        message("Reading ", rmd.source)
        textIn <- readLines(rmd.source)
      } else {
        stop("'rmd.source' not specified and multiple .rmd files found in working directory", call.=FALSE)
      }
    }
  } else {
    if (nchar(tools::file_ext(rmd.source)) == 0) {
      rmd.source <- paste0(rmd.source, ".rmd")
    } else if (tools::file_ext(rmd.source) != "rmd") stop("rmd.source extension must be rmd", call.=FALSE)
    
    if (!file.exists(rmd.source)) stop(rmd.source, " not found in working directory", call.=FALSE)
    message("Reading ", rmd.source)
    textIn <- readLines(rmd.source)
  }
  
  ## Create Markdown file ####
  if (!dir.exists('build')) {
    dir.create('build')
    tt <- file.create('build/readme.txt')
    cat("This directory was created by DLMtool function OMdoc\n\n", sep="", append=TRUE, file='build/readme.txt') 
    cat("Files in this directory are used to generate the OM report.\n\n", sep="", append=TRUE, file='build/readme.txt') 
  } 

  if (is.null(out.file)) out.file <- tools::file_path_sans_ext(rmd.source)
  out.file <- gsub("_source", "_compiled", out.file)
  
  RMDfile <- paste0("build/", out.file, ".Rmd")
  # if (file.exists(RMDfile) & !overwrite) {
  #   stop(RMDfile, " already exists.\n Provide a different output file name ('out.file') or use 'overwrite=TRUE'")
  # } else {
  #   message('Writing ', RMDfile)
  tt <- file.create(RMDfile)
  # }

  
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
    if (output == "html_document") {
      cat("     toc_float: true\n", append=TRUE, file=RMDfile, sep="")
      cat("     theme: ", theme, "\n", append=TRUE, file=RMDfile, sep="")
    }

    
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
  
  ## knitr options ####
  # cat("```{r setup, include=FALSE}\n", append=TRUE, file=RMDfile, sep="")
  # cat("library(knitr)\n", append=TRUE, file=RMDfile, sep="")
  # cat("opts_chunk$set(dev = 'pdf')\n", append=TRUE, file=RMDfile, sep="")
  # cat("```\n", append=TRUE, file=RMDfile, sep="")

  ## Generate Sampled Parameters ####
  
  Pars <- NULL
  if (inc.plot) {
    
    # --- Generate Historical Samples ----
    # Only required if the OM has changed #
    runSims <- FALSE
    fileName <- OM@Name
    fileName <- gsub(" ", "", gsub("[[:punct:]]", "", fileName))
    if (nchar(fileName)>15) fileName <-  substr(fileName, 1, 15)
      
    
    if (file.exists(paste0('build/', fileName, '.dat'))) {
      # OM has been documented before - check if it has changed
      testOM <- readRDS(paste0("build/", fileName, '.dat'))
      if (class(testOM) == 'OM') {
        # check if OM has changed 
        changed <- rep(FALSE, length(slotNames(OM)))
        for (sl in seq_along(slotNames(OM))) {
          oldOM <- slot(OM, slotNames(OM)[sl])
          newOM <- slot(testOM, slotNames(OM)[sl])
          if (class(oldOM) !='character') {
            if (class(oldOM) != 'list') {
              if (length(oldOM)<1 || !is.finite(oldOM)) oldOM <- 0
              if (length(newOM)<1 || !is.finite(newOM)) newOM <- 0
              if (any(oldOM != newOM)) changed[sl] <- TRUE
            } else {
              if (length(oldOM) != length(newOM)) changed[sl] <- TRUE
              if (length(setdiff(oldOM, newOM)) > 0) changed[sl] <- TRUE
            }
          }
        }
        if (sum(changed)>0) runSims <- TRUE 
        if (sum(changed) == 0) {
          out <-  readRDS(paste0('build/', fileName, 'Hist.dat'))
          Pars <- c(out$SampPars, out$TSdata, out$MSYs)  
        }
      } else {
        file.remove(paste0('build/',fileName, '.dat'))
        runSims <- TRUE
      }
     
    } else{
      saveRDS(OM, file=paste0('build/', fileName, '.dat'))
      runSims <- TRUE
    }
    
    if (runSims) {
      message("\n\nRunning Historical Simulations\n\n")
      OM <- updateMSE(OM) # update and add missing slots with default values
      if (OM@nsim > 48) setup()
      out<- runMSE(OM,Hist=T)
      Pars <- c(out$SampPars, out$TSdata, out$MSYs)  
      saveRDS(out, file=paste0('build/', fileName, 'Hist.dat'))
    }
  }
  
  
  ## Input text ####
  # ind <- which(unlist(lapply(textIn, nchar)) == 0) # delete empty lines 
  # if (length(ind) > 0) textIn <- textIn[-ind]
  ind <- grep("# Introduction", textIn)
  if (length(ind)>1) stop("# Introduction heading found more than once", call.=FALSE)
  if (length(ind)>0) {
    textIn <- textIn[ind:length(textIn)]
  } else {
    ind <- grep("# Stock Parameters", textIn)
    if (length(ind)>1) stop("# Stock Parameters heading found more than once", call.=FALSE)
    if (length(ind) == 0) stop("# Stock Parameters not found", call.=FALSE)
    textIn <- textIn[ind:length(textIn)]
  }
  
  
  ## Introduction ####
  writeSection(class="Intro", OM, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## Cpars ####
  if (useCpars) writeSection(class="cpars", OM, textIn, RMDfile, color=color, 
                             inc.plot=inc.plot)

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
    cat("```{r plot.OM, echo=FALSE, fig.asp=2}\n", append=TRUE, file=RMDfile, sep="")
    cat("plot.OM(out)\n", append=TRUE, file=RMDfile, sep="")
    cat("```\n\n\n", append=TRUE, file=RMDfile, sep="")

  }
  
  
  ## References ####
  writeSection(class="References", OM, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## Render Markdown ####
  if (render) {
    RMDfileout <- gsub("_compiled", "", tools::file_path_sans_ext(RMDfile))
    if (output == "html_document") RMDfileout <- paste0(unlist(strsplit(RMDfileout, "/"))[2], ".html")
    if (output == "pdf_document") RMDfileout <- paste0(unlist(strsplit(RMDfileout, "/"))[2], ".pdf")

    message("\n\nRendering markdown document as ", RMDfileout)
    EffYears <- seq(from=(OM@CurrentYr -  OM@nyears + 1), to=OM@CurrentYr, length.out=length(OM@EffYears))
    EffYears <- round(EffYears,0)
    Effvals <- data.frame(EffYears=EffYears, EffLower=OM@EffLower, EffUpper=OM@EffUpper)
    params <- list(OM=OM, Pars=Pars, Effvals=Effvals, out=out)
    rmarkdown::render(input=RMDfile, output_file=RMDfileout, output_format=output, output_dir=getwd(), param=params)
    
    utils::browseURL(file.path(getwd(), RMDfileout))
    
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
               "Spatial distribution and movement: Size_area_1, Frac_area_1, Prob_staying",
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
      if (length(tt) < 1) stop("slot ", X, " not found in ", type, " template", call.=FALSE)
      if (length(tt) > 1) stop("slot ", X, " found multiple times in ", type, " template", call.=FALSE)    
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
    if (intro == 0) stop("# Introduction heading not found", call.=FALSE)
    if (length(intro) == 1) {
      cat(fstHmd[intro], "\n\n", sep="", append=TRUE, file=RMDfile) # write first heading
      for (ll in (intro+1):(fLH[intro+1] - 1)) {
        cat(textIn[ll], "\n\n", sep="", append=TRUE, file=RMDfile) # write intro paragraphs
      }
    } else {
      stop("More than one section # Introduction", call.=FALSE)
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
      
    } else stop("More than one section # References", call.=FALSE)
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
    }  else stop("More than one section # Custom Parameters", call.=FALSE)
    
  } else {
    # Write class heading 
    st <- which(trimws(gsub("#", "", textIn)) == paste(class, "Parameters"))
    sta <- which(fstH == paste(class, "Parameters"))
    if (length(st) > 1) stop("Multiple '# ", class, " Parameters' headings in document.", call.=FALSE)
    if (length(st) < 1) stop("'# ", class, " Parameters' heading not found in document.", call.=FALSE)
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
      stop("Invalid second level headings (must match slots in class ", class, "): ", paste(invalid, ""), call.=FALSE)
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
            if (used & sl != "Source" & sl != "Name") description <- paste("No justification provided.", valtext)
            if (!used) description <- "Slot not used. "
          }

          description[nchar(description) == 0] <- "\n\n"
          description[length(description)-1] <- "  "
          
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


