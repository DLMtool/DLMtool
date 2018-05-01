#' Copy example OM XL and OM Documentation 
#'
#' @param dir the file path to copy the files to 
#' @export
#'
#' @examples
#' \dontrun{
#' OMexample()
#' }
OMexample <- function(dir) {
  fromRMD <- system.file("Example_Chile_Hake.Rmd", package="DLMtool")
  tt <- file.copy(fromRMD, dir, overwrite = TRUE)
  fromXL <- system.file("Example_Chile_hake.xlsx", package="DLMtool")
  tt <- file.copy(fromXL, dir, overwrite = TRUE)
}

#' Initialize Operating Model
#'
#' Generates an Excel spreadsheet and a source.rmd file in the current working directory for 
#' specifying and documenting a DLMtool Operating Model.
#' 
#' @param name The name of the Excel and source.rmd file to be created in the working directory (character). 
#' Use 'example' for a populated example OM XL and documentation file.
#' @param ... Optional DLMtool objects to use as templates: OM, Stock, Fleet, Obs, or Imp objects
#' @param files What files should be created: 'xlsx', 'rmd', or c('xlsx', 'rmd') (default: both)
#' to use as templates for the Operating Model.
#' @param dir Optional file path to create the xlsx and rmd files. Default is `getwd()`
#' @param overwrite Logical. Should files be overwritten if they already exist?
#'
#' @return name.xlsx and name.rmd files are created in the working directory.  
#' @export
#' @author A. Hordyk
#'
#' @examples
#' \dontrun{
#' # Create an Excel OM template and rmd file called 'myOM.xlsx' and 'myOM.rmd': 
#' OMinit('myOM')
#' 
#' # Create an Excel OM template and text file called 'myOM.rmd' and 'myOM.rmd', using
#' another OM as a template: 
#' OMinit('myOM', myOM)
#' 
#' # Create an Excel OM template and text file called 'myOM.rmd' and 'myOM.rmd', using
#' the Stock object 'Herring' as a template: 
#' OMinit('myOM', Herring)
#' 
#' # Create an Excel OM template and text file called 'myOM.rmd' and 'myOM.rmd', using
#' the Stock object 'Herring', and Obs object 'Generic_obs' as templates: 
#' OMinit('myOM', Herring, Generic_obs)
#' }
#' 
OMinit <- function(name=NULL, ..., files=c('xlsx', 'rmd'), dir=NULL, overwrite=FALSE) {
  files <- match.arg(files, several.ok = TRUE)
  
  if(is.null(dir)) dir <- getwd()
  
  if (is.null(name)) stop("Require OM name", call.=FALSE)
  
  if (tolower(name) == 'example') {
    OMexample(dir)
    return(message("Creating Example Files in ", dir))
  }
  if (class(name) != 'character') stop("name must be text", call.=FALSE)
 
  ## Create Folders ####
  if(!dir.exists(file.path(dir,'data'))) dir.create(file.path(dir,'data'))
  if(!dir.exists(file.path(dir,'docs'))) dir.create(file.path(dir,'docs'))
  if(!dir.exists(file.path(dir,'images'))) dir.create(file.path(dir,'images'))
  if(!dir.exists(file.path(dir,'robustness'))) dir.create(file.path(dir,'robustness'))

  ## Write Excel skeleton ####
  if (nchar(tools::file_ext(name)) == 0) {
    nameNoExt <- name
    name <- paste0(name, ".xlsx")
  } else {
    ext <- tools::file_ext(name)
    if (!ext %in% c("xlsx", "xls")) stop("File extension must be 'xlsx' or 'xls'", call.=FALSE)
    nameNoExt <- tools::file_path_sans_ext(name)
  }
  
  InTemplates <- list(...)
  ObTemplates <- list()
  useTemp <- FALSE
  if (length(InTemplates) >0) {
     inclasses <- unlist(lapply(InTemplates, class))
    if (!is.null(inclasses)) {
      # check if zip application exists
      chck <- Sys.which("zip") # requires 'zip.exe' on file path
      if (nchar(chck) <1) {
        message('zip application is required for templates. If a zip application is installed on your machine you may need to add it to the path. Try:')
        message('path <- Sys.getenv("PATH")')
        message('Sys.setenv("PATH" = paste(path, "path_to_zip.exe", sep = ";"))')
        stop("Can't use templates without zip application. You may need to install Rtools to use templates", call.=FALSE)
      }
    }
   
    for (x in seq_along(inclasses)) {
      if (!inclasses[x] %in% c("Stock", "Fleet", "Obs", "Imp", "OM")) stop(InTemplates[[x]], " is not a valid DLMtool object")
    }
    isOM <- which(inclasses == "OM")
    if (length(isOM)>0) {
      message("\nUsing OM Template")
      ObTemplates$Stock <- SubOM(InTemplates[[isOM]], "Stock")
      if (is.na(ObTemplates$Stock@Name) || nchar(ObTemplates$Stock@Name)==0) ObTemplates$Stock@Name <- InTemplates[[isOM]]@Name
      ObTemplates$Fleet <- SubOM(InTemplates[[isOM]], "Fleet")
      if (is.na(ObTemplates$Fleet@Name) || nchar(ObTemplates$Fleet@Name)==0) ObTemplates$Fleet@Name <- InTemplates[[isOM]]@Name
      ObTemplates$Obs <- SubOM(InTemplates[[isOM]], "Obs")
      if (is.na(ObTemplates$Obs@Name) || nchar(ObTemplates$Obs@Name)==0) ObTemplates$Obs@Name <- InTemplates[[isOM]]@Name
      ObTemplates$Imp <- SubOM(InTemplates[[isOM]], "Imp")
      if (is.na(ObTemplates$Imp@Name) || nchar(ObTemplates$Imp@Name)==0) ObTemplates$Imp@Name <- InTemplates[[isOM]]@Name
      useTemp <- TRUE
    } else {
      for (x in seq_along(inclasses)) {
        if (inclasses[x] == 'Stock') ObTemplates$Stock <- InTemplates[[x]]
        if (inclasses[x] == 'Fleet') ObTemplates$Fleet <- InTemplates[[x]]
        if (inclasses[x] == 'Obs') ObTemplates$Obs <- InTemplates[[x]]
        if (inclasses[x] == 'Imp') ObTemplates$Imp <- InTemplates[[x]]
      }
      nm <- names(ObTemplates)  
      message("\n\nUsing Object Templates:")
      useTemp <- TRUE
      for (X in nm) {
        message(ObTemplates[[X]]@Name)
      }
    }
  }

  if ('xlsx' %in% files) {
   
    # Copy xlsx file over to working directory 
    # Copy the Excel File ####
    message("Creating ", name, " in ", dir)
    path <- system.file("OM.xlsx", package = "DLMtool")
    pathout <- gsub("OM.xlsx", name, path)
    pathout <- gsub(dirname(pathout), dir, pathout)
    
    # Check if file exists 
    exist <- file.exists(pathout)
    if (exist & !overwrite) stop(name, " already exists in working directory. Use 'overwrite=TRUE' to overwrite", call.=FALSE)
    copy <- file.copy(path, pathout, overwrite = overwrite)
    if (!copy) stop("Excel file not copied from ", path)
    
    # loop through slot values if Obj template provided
    if (useTemp) {
      wb <- openxlsx::loadWorkbook(file.path(dir, name))
      names <- c("Stock", "Fleet", "Obs", "Imp")
      for (objname in names) {
        if (!is.null(ObTemplates[objname])) {
          obj <- ObTemplates[objname][[1]]
          slots <- slotNames(obj)
          
          for (sl in seq_along(slots)) {
            val <- slot(obj, slotNames(objname)[sl])
            ln <- length(val)
            if (ln >0 && !is.na(ln)) {
              df <- data.frame(t(val))
              openxlsx::writeData(wb, sheet = objname, x = df, 
                                  startCol = 2, startRow = sl+1,
                                  colNames = FALSE, rowNames = FALSE, 
                                  withFilter = FALSE,
                                  keepNA = FALSE)         
            }
            
          }
          openxlsx::setColWidths(wb, sheet = objname, cols = 1, widths = 'auto')
        }
      }
      
      # OM tab not currently updated
      openxlsx::saveWorkbook(wb, file.path(dir,name), overwrite = TRUE)
    
      
    }
    
      
    }
    

  if ('rmd' %in% files) { 
    # RMD File ####
    rmdname <- paste0(nameNoExt, '.rmd')
    message("Creating ", rmdname, " in ", dir)
    path <- system.file("OM.rmd", package = "DLMtool")
    if (nchar(path) <1) stop("OM.rmd not found in DLMtool package")
    pathout <- gsub("OM.rmd", rmdname, path)
    pathout <- gsub(dirname(pathout), dir, pathout)
    
    # Check if file exists 
    exist <- file.exists(pathout)
    if (exist & !overwrite) stop(rmdname, " alread exists in ", dir, ". Use 'overwrite=TRUE' to overwrite", call.=FALSE)
    copy <- file.copy(path, pathout, overwrite = overwrite)
    if (!copy) stop("Rmd file not copied from ", path)
    
    # Copy over templates - if used ####
    if (length(ObTemplates)>0) {
      names <- c("Stock", "Fleet", "Obs", "Imp")
      textIn <- readLines(file.path(dir,rmdname))
      for (objname in names) {
        if (!is.null(ObTemplates[objname])) {
          obj <- ObTemplates[objname][[1]]
          slots <- slotNames(obj)
          
          for (sl in slots) {
            if (!sl %in% c("Name", "Source")) {
              lineno <- grep(paste0("^## ", sl, "$"), textIn)
              textIn[lineno+1] <- paste("Borrowed from:", obj@Name)
            }
            
          }
        }
      }
      writeLines(textIn, con = file.path(dir, rmdname), sep = "\n", useBytes = FALSE)
    }
  } 
  
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
    fls <- fls[!grepl('~', fls)]
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
      warning("No parameter values found in Sheet: ", obj, ". Using defaults", call.=FALSE)
      tempObj[[count]] <- new(obj)
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
  dat <- sht # sht[,1:2] 
  dat <- dat[which(dat[,1] != "Slot"),]
  # if (ncol(sht)>2) warning("More than two columns found in Sheet OM. Values in columns C+ are ignored")
  if (ncol(sht)<2) {
    message("No values found for OM slots in Sheet OM. Using defaults")
  } else {
    for (xx in 1:nrow(dat)) {
      val <- dat[xx, 2:ncol(dat)]
      if (length(val)) {
        if (!dat[xx,1] %in% c("Name", "Agency", "Region", "Sponsor")) {
          options(warn=-1)
          val <- as.numeric(val)
          options(warn=1)
          val <- val[!is.na(val)]
          if (.hasSlot(OM, dat[xx,1])) slot(OM, dat[xx, 1]) <- val
        } else  {
          val <- val[!is.na(val)]
          if (.hasSlot(OM, dat[xx,1])) slot(OM, dat[xx, 1]) <- val
        }
        
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
  ChkObj(OM, FALSE)
  if (msg) {
    message('OM successfully imported\n')
    message("Document OM slots in .rmd file (probably ", tools::file_path_sans_ext(name), ".rmd),
  and run 'OMdoc' if OM parameter values have changed." )
  }

  OM
}





#' Generate OM Documentation Report
#'
#' @param OM An object of class 'OM' or the name of an OM xlsx file 
#' @param rmd.source Optional. Name of the source.rmd file corresponding to the 'OM'. Default assumption
#' is that the file is 'OM@Name.Rmd'
#' @param overwrite Logical. Should existing files be overwritten?
#' @param out.file Optional. Character. Name of the output file. Default is the same as the text file.
#' @param inc.plot Logical. Should the plots be included?
#' @param render Logical. Should the document be compiled? May be useful to turn off if 
#' there are problems with compililing the Rmd file.
#' @param output Character. Output file type. Default is 'html_document'. 'pdf_document' is available
#' but may require additional software and have some formatting issues.
#' @param openFile Logical. Should the compiled file be opened in web browser?
#' @param quiet TRUE to supress printing of the pandoc command line.
#' @param dir Optional file path to read the xlsx and rmd files. Default is `getwd()`
#' @param ... Optional additional named arguments provided to `runMSE`
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
                  inc.plot=TRUE, render=TRUE, output="html_document", openFile=TRUE, quiet=FALSE,
                  dir=NULL, ...) {
  # markdown compile options
  toc=TRUE; color="blue";  theme="flatly"
  if (is.null(dir)) dir <- getwd()
  OMXLname <- NULL
  
  if (class(OM) == "OM") {
    # nothing
  } else if (class(OM) == 'character') {
    OMXLname <- OM
    OM <- XL2OM(file.path(dir,OM), msg=FALSE)
  } else if (is.null(OM)) {
    fls <- list.files(path=dir, pattern=".xlsx", ignore.case=TRUE)
    fls <- fls[!grepl('~', fls)]
    if (length(fls)==1) OM <- XL2OM(file.path(dir,fls), msg=FALSE)
    if (length(fls)>1) stop('argument "OM" not provided and multiple .xlsx files in ', dir, call.=FALSE)
  } else stop('OM must be class "OM" or name of OM xlsx file.', call.=FALSE)
  
  if (is.null(OM)) stop('OM not imported. Is the name correct?', call.=FALSE)
  ## Read in Rmd.source file ####
  if (is.null(rmd.source)) {
    rmd.source <- list.files(path=dir, pattern=".rmd", ignore.case=TRUE)
    if (length(rmd.source) == 0) stop("rmd.source' not specified and no .rmd files found in ", dir, call.=FALSE)
    if (length(rmd.source) == 1) {
      textIn <- readLines(file.path(dir,rmd.source))
    } else {
      NoExt <- tools::file_path_sans_ext(rmd.source)
      if (!is.null(OMXLname)) ind <- which(tolower(NoExt) == tolower(paste0(OMXLname, "_source")))
      if (is.null(OMXLname)) ind <- which(tolower(NoExt) == tolower(paste0(OM@Name, "_source")))
      if (length(ind) > 0) {
        rmd.source <- rmd.source[ind]
        message("Reading ", file.path(dir,rmd.source))
        textIn <- readLines(file.path(dir,rmd.source))
      } else {
        stop("'rmd.source' not specified and multiple .rmd files found in ", dir, call.=FALSE)
      }
    }
  } else {
    if (nchar(tools::file_ext(rmd.source)) == 0) {
      rmd.source <- paste0(rmd.source, ".rmd")
    } else if (tools::file_ext(rmd.source) != "rmd") stop("rmd.source extension must be rmd", call.=FALSE)
    
    if (!file.exists(rmd.source)) stop(rmd.source, " not found in ", dir, call.=FALSE)
    message("Reading ", file.path(dir,rmd.source))
    textIn <- readLines(file.path(dir,rmd.source))
  }
  
  ## Create Markdown file ####
  if (!dir.exists(file.path(dir, 'build'))) {
    dir.create(file.path(dir, 'build'))
    tt <- file.create(file.path(dir,'build/readme.txt'))
    cat("This directory was created by DLMtool function OMdoc\n\n", sep="", append=TRUE, file=file.path(dir,'build/readme.txt')) 
    cat("Files in this directory are used to generate the OM report.\n\n", sep="", append=TRUE, file=file.path(dir,'build/readme.txt')) 
  } 
  
  if(dir.exists(file.path(dir,"images"))) {
    dir.create(file.path(dir,'build/images'), showWarnings = FALSE)
    cpy <- file.copy(file.path(dir,'images'), file.path(dir,'build'), overwrite=TRUE, recursive = TRUE)
  }

  if (is.null(out.file)) out.file <- tools::file_path_sans_ext(rmd.source)
  # out.file <- gsub("_source", "_compiled", out.file)
  
  RMDfile <- file.path(dir, paste0("build/", out.file, ".Rmd"))
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
  

  ## knitr options ####
  # cat("```{r setup, include=FALSE}\n", append=TRUE, file=RMDfile, sep="")
  # cat("library(knitr)\n", append=TRUE, file=RMDfile, sep="")
  # cat("opts_chunk$set(dev = 'pdf')\n", append=TRUE, file=RMDfile, sep="")
  # cat("```\n", append=TRUE, file=RMDfile, sep="")

  ## Generate Sampled Parameters ####
  
  Pars <- NULL
  out <- NULL
  if (inc.plot) {
    
    # --- Generate Historical Samples ----
    # Only required if the OM has changed #
    runSims <- FALSE
    fileName <- OM@Name
    fileName <- gsub(" ", "", gsub("[[:punct:]]", "", fileName))
    if (nchar(fileName)>15) fileName <-  substr(fileName, 1, 15)
    
    if (file.exists(paste0(file.path(dir, 'build/', paste0(fileName, '.dat'))))) {
      # OM has been documented before - check if it has changed
      testOM <- readRDS(file.path(dir,paste0("build/", fileName, '.dat')))
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
              if (length(oldOM) != length(newOM)) {
                changed[sl] <- TRUE
              } else if (length(oldOM)>0){
                for (xx in 1:length(oldOM)) {
                  if(any(oldOM[[xx]] != newOM[[xx]]))changed[sl] <- TRUE
                  
                }
              }
            }
          }
        }
        if (sum(changed)>0) runSims <- TRUE 
        if (sum(changed) == 0) {
          out <-  readRDS(file.path(dir,paste0('build/', fileName, 'Hist.dat')))
          Pars <- c(out$SampPars, out$TSdata, out$MSYs)  
        }
      } else {
        file.remove(file.path(dir,paste0('build/',fileName, '.dat')))
        runSims <- TRUE
      }
      
    } else{
      saveRDS(OM, file=file.path(dir,paste0('build/', fileName, '.dat')))
      runSims <- TRUE
    }
    
    if (runSims) {
      message("\nRunning Historical Simulations")
      OM2 <- updateMSE(OM) # update and add missing slots with default values
      if (OM2@nsim > 48) {
        message("nsim too high for documentation purposes. Running here with nsim=48")
        OM2@nsim <- 48
      }
      out<- runMSE(OM2,Hist=T, parallel = FALSE, silent=TRUE, ...)
      Pars <- c(out$SampPars, out$TSdata, out$MSYs)  
      saveRDS(out, file=file.path(dir,paste0('build/', fileName, 'Hist.dat')))
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
  writeSection(class="Intro", OM, Pars, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## OM Details ####
  OMdesc <- DLMtool::OMDescription
  cat("# Operating Model \n", append=TRUE, file=RMDfile, sep="")
  
  ## Link to OM object ####
  chkFile <- file.exists("OM.rdata")
  if (chkFile) {
    cat("The OM rdata file can be downloaded from [here](OM.rdata)\n\n", append=TRUE, file=RMDfile, sep="")
    cat("Download and import into R using `myOM <- readRDS('OM.rdata')` \n\n", append=TRUE, file=RMDfile, sep="")
  } 
  
  
  # Taxonomic Info and Location ####
  if (.hasSlot(OM, 'Species') && !is.na(OM@Species)) {
    cat("## Species Information \n\n", append=TRUE, file=RMDfile, sep="")
    cat("**Species**: *", OM@Species, "*\n\n", append=TRUE, file=RMDfile, sep="")
    cat("**Common Name**: *", OM@Common_Name, "*\n\n", append=TRUE, file=RMDfile, sep="")
    cat("**Management Agency**: ", OM@Agency, "\n\n", append=TRUE, file=RMDfile, sep="")
    cat("**Region**: ", OM@Region, "\n\n", append=TRUE, file=RMDfile, sep="")
    if (length(OM@Sponsor)>0) cat("**Sponsor**: ", OM@Sponsor, "\n\n", append=TRUE, file=RMDfile, sep="")
    if (length(OM@Latitude)>0) {
      lat <- paste0(OM@Latitude, sep="", collapse=", ")
      cat("**Latitude**: ", lat, "\n\n", append=TRUE, file=RMDfile, sep="")
    }
    if (length(OM@Longitude)>0) {
      long <- paste0(OM@Longitude, sep="", collapse=", ")
      cat("**Longitude**: ", long, "\n\n", append=TRUE, file=RMDfile, sep="")
    }
  }
  
  cat("## OM Parameters \n", append=TRUE, file=RMDfile, sep="")
  cat("**OM Name**: ", OMdesc$Description[OMdesc$Slot =='Name'], ": ", "<span style='color:", color, "'>", " ", OM@Name, "</span>", "\n\n", append=TRUE, file=RMDfile, sep="")
  
  cat("**nsim**: ", OMdesc$Description[OMdesc$Slot =='nsim'], ": ", "<span style='color:", color, "'>", " ", OM@nsim, "</span>", "\n\n", "\n\n", append=TRUE, file=RMDfile, sep="")
  
  cat("**proyears**: ", OMdesc$Description[OMdesc$Slot =='proyears'], ": ", "<span style='color:", color, "'>", " ", OM@proyears, "</span>", "\n\n", "\n\n", "\n\n", append=TRUE, file=RMDfile, sep="")

  cat("**interval**: ", OMdesc$Description[OMdesc$Slot =='interval'], " ", "<span style='color:", color, "'>", " ", OM@interval, "</span>", "\n\n",append=TRUE, file=RMDfile, sep="")
  
  cat("**pstar**: ", OMdesc$Description[OMdesc$Slot =='pstar'], ": ", "<span style='color:", color, "'>", " ", OM@pstar, "</span>", "\n\n",append=TRUE, file=RMDfile, sep="")
  
  cat("**maxF**: ", OMdesc$Description[OMdesc$Slot =='maxF'], ": ", "<span style='color:", color, "'>", " ", OM@maxF, "</span>", "\n\n", append=TRUE, file=RMDfile, sep="")
  
  cat("**reps**: ", OMdesc$Description[OMdesc$Slot =='reps'], " ", "<span style='color:", color, "'>", " ", OM@reps, "</span>", "\n\n", append=TRUE, file=RMDfile, sep="")
 
  cat("**Source**: ", OMdesc$Description[OMdesc$Slot =='Source'], "\n\n", "<span style='color:", color, "'>", " ", OM@Source, "</span>", "\n\n", append=TRUE, file=RMDfile, sep="")
  
  ## 
  

  
  useCpars <- length(OM@cpars) >0 
  ## Cpars ####
  if (useCpars) {
    message("writing cpars section")
    writeSection(class="cpars", OM, Pars, textIn, RMDfile, color=color, 
                             inc.plot=inc.plot)
  }
  ## Stock Parameters ####
  message("writing Stock section")
  writeSection(class="Stock", OM, Pars, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## Fleet Parameters ####
  message("writing Fleet section")
  writeSection(class="Fleet", OM, Pars, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## Observation Parameters ####
  message("writing Obs section")
  writeSection(class="Obs", OM, Pars, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## Implementation Parameters ####
  message("writing Imp section")
  writeSection(class="Imp", OM, Pars, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## OM Plots ####
  if (inc.plot) {
    cat("# OM Plots\n\n", sep="", append=TRUE, file=RMDfile) # write heading
    cat("```{r plotOM, echo=FALSE, fig.asp=2}\n", append=TRUE, file=RMDfile, sep="")
    cat("plot.OM(out)\n", append=TRUE, file=RMDfile, sep="")
    cat("```\n\n\n", append=TRUE, file=RMDfile, sep="")

  }
  
  
  ## References ####
  message("writing Reference section")
  writeSection(class="References", OM, Pars, textIn, RMDfile, color=color, inc.plot=inc.plot)
  
  ## Render Markdown ####
  if (render) {
    RMDfileout <- gsub("_compiled", "", tools::file_path_sans_ext(RMDfile))
  
    if (output == "html_document") RMDfileout <- paste0(basename(RMDfileout), ".html")
    if (output == "pdf_document") RMDfileout <- paste0(basename(RMDfileout), ".pdf")

    message("\n\nRendering markdown document as ", RMDfileout)
    
    EffYears <- seq(from=(OM@CurrentYr -  OM@nyears + 1), to=OM@CurrentYr, length.out=length(OM@EffYears))
    EffYears <- round(EffYears,0)
    if (length(OM@cpars$Find)>0) {
      lower <- as.numeric(signif(apply(OM@cpars$Find, 2, min),3))
      upper <- as.numeric(signif(apply(OM@cpars$Find, 2, max),3))
      Effvals <- data.frame(EffYears=EffYears, EffLower=lower, EffUpper=upper)
    } else {
      Effvals <- data.frame(EffYears=EffYears, EffLower=signif(OM@EffLower,3), EffUpper=signif(OM@EffUpper,3))
    }
  
    params <- list(OM=OM, Pars=Pars, Effvals=Effvals, out=out)
    knitr::knit_meta(class=NULL, clean = TRUE)
    rmarkdown::render(input=RMDfile, output_file=RMDfileout, output_format=output, 
                      output_dir=dir, param=params, quiet=quiet)
    
    if (openFile) utils::browseURL(RMDfileout)
    
  } else {
    
  }
  
}




# Text templates for the OM documentation ####
Template <- function(type=c("Stock", "Fleet", "Obs", "Imp")) {
  type <- match.arg(type)
  if (type == "Stock") mat <- 
      matrix(c("Mortality and age:  maxage, R0, M, M2, Mexp, Msd, Mgrad",
               "Recruitment: h, SRrel, Perr, AC",
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
        "Current Year: CurrentYr",
        "Existing Spatial Closures: MPA"), ncol=1)
  if (type == "Obs") mat <- 
      matrix(c(
        "Catch statistics: Cobs, Cbiascv, CAA_nsamp, CAA_ESS, CAL_nsamp, CAL_ESS",
        "Index imprecision, bias and hyperstability: Iobs, Ibiascv, Btobs, Btbiascv, beta",
        "Bias in maturity, natural mortality rate and growth parameters: LenMbiascv, Mbiascv, Kbiascv,t0biascv, Linfbiascv",
        "Bias in length at first capture, length at full selection: LFCbiascv, LFSbiascv",
        "Bias in fishery reference points, unfished biomass, FMSY, FMSY/M ratio, biomass at MSY relative to unfished: FMSYbiascv, FMSY_Mbiascv, BMSY_B0biascv",
        "Management targets in terms of the index (i.e., model free), the total annual catches and absolute biomass levels: Irefbiascv, Crefbiascv, Brefbiascv",
        "Depletion bias and imprecision: Dbiascv, Dobs",
        "Recruitment compensation and trend: hbiascv, Recbiascv"), ncol=1)
        # "Currently unused observation processes - bias in unfished biomass, intrinsic rate of increase, annual increase in fishing efficiency and age at 50% vulnerability, bias and imprecision in current fishing rate, bias in maximum age: B0cv, rcv, Fcurbiascv, Fcurcv, maxagecv"), ncol=1)
  if (type == "Imp") mat <-
      matrix(c(
        "Output Control Implementation Error: TACFrac, TACSD",
        "Effort Control Implementation Error: TAEFrac, TAESD", 
        "Size Limit Control Implementation Error: SizeLimFrac, SizeLimSD"), ncol=1)
 
  # Check slots 
  Slots <- names(methods::getSlots(type))
  for (X in Slots) {
    tt <- grep(paste0("\\<", X, "\\>"), mat) 
    if (X != "Name" && X != "Source" && X!="Species" && X!="Common_Name" && X!="Latitude" && X!='Longitude') {
      if (length(tt) < 1) stop("slot ", X, " not found in ", type, " template", call.=FALSE)
      if (length(tt) > 1) stop("slot ", X, " found multiple times in ", type, " template", call.=FALSE)    
    }
  }
  return(mat)
}


writeSection <- function(class=c("Intro", "Stock", "Fleet", "Obs", "Imp", "References", "cpars"), OM, Pars,
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
    } 
    if (length(cparstext) > 1) stop("More than one section # Custom Parameters", call.=FALSE)
    
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
    
    # check for text after class heading 
    preText <- textIn[bg:(bg+scLHloc[1]-2)]
    if (any(nchar(preText)>0)) {
      cat(preText, "\n\n", append=TRUE, file=RMDfile, sep="")
    }
  
    
    # Get template for class section 
    ClTemp <- Template(class)
    
    # loop through template and write section 
    for (rr in 1:nrow(ClTemp)) {
      if (grepl("^Currently unused", ClTemp[rr,1])) {
        temptext <- trimws(unlist(strsplit(ClTemp[rr,], "-")))
        cat("### ", temptext[1], "\n\n", append=TRUE, file=RMDfile, sep="")
        cat("*", temptext[2], "*\n\n", append=TRUE, file=RMDfile, sep="")
      } else {
        slots <- trimws(unlist(strsplit(strsplit(ClTemp[rr,], ":")[[1]][2], ",")))
        cat("### ", ClTemp[rr,], "\n\n", append=TRUE, file=RMDfile, sep="")
        for (sl in slots) {
          # get slot value if not in cpars 
          if (useCpars && sl %in% names(OM@cpars)) {
            if (length(Pars[[sl]])>0) {
              val <- range(Pars[[sl]])
            } else {
              val <- range(OM@cpars[[sl]])
            }
            val <- round(val,2)
            used <- TRUE
            val <- gsub('"', "", paste(val, collapse="\", \""))
            valtext <- paste0("Specified in cpars: ", "<span style='color:", color, "'>", " ", trimws(val), "</span>", "\n\n")
          } else {
            val <- slot(OM, sl) 
            # if (length(Pars[[sl]])>0) {
            #   if (length(Pars[[sl]])==1) val <- (Pars[[sl]])
            #   if (length(Pars[[sl]])>1) {
            #     if (all(Pars[[sl]] == mean(Pars[[sl]]))) {
            #       val <- mean(Pars[[sl]])
            #     } else  val <- range(Pars[[sl]])
            #   }
            # } else {
            #   val <- slot(OM, sl)  
            # }
            if (is.numeric(val)) val <- round(val,2)
            used <- length(val)>0 && !is.na(val) && !is.null(val) # is the slot used?
            if (used) {
              val <- gsub('"', "", paste(val, collapse="\", \""))
              valtext <- paste0("Specified Value(s): ", "<span style='color:", color, "'>", " ", trimws(val), "</span>", "\n\n")
            } else {
              valtext <- val
            }           
          }
          
          loc <- which(scLH == sl)
          if (length(loc) > 0) {
            bg <- scLHloc[loc]+1
            end <- scLHloc[loc+1]-1
            if (is.na(end)) end <- length(Text)
            description <- Text[bg:end]
            description <- paste(description, "\n")
            if (!used) description <- c("Slot not used.")
            if (used && ! sl%in% c("EffYears", "EffLower", "EffUpper")) description <- c(valtext, description)

          } else {
            if (used & sl != "Source" & sl != "Name") description <- paste(valtext, "No justification provided.")
            if (!used) description <- "Slot not used. "
          }

          description[description == " \n"] <- "\n\n"
          # description[length(description)-1] <- "  "
          
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
            cat("knitr::kable(Effvals, format='markdown', caption='')\n", append=TRUE, file=RMDfile, sep="")
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
    
    if ((class == "Obs" | class =='Imp') & inc.plot) {
      plotText(OM, slots=class, RMDfile)
    }
  }
}



writeCSV2 <- function(inobj, tmpfile = NULL, objtype = c("Stock", "Fleet", 
                                                         "Obs", "Imp", "Data", "OM")) {
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





plotText <- function(OM, slots, RMDfile) {
  if (any(c("M", "h", "Linf", "L50", "D", "EffUpper", "qcv", "Vmaxlen", "DR") %in% slots)) {
    # slotstext <- paste("c(", paste(slots, sep=",", collapse = ","), ")")
    slotstext <- slots[slots %in% c("M", "h", "Linf", "L50", "D", "EffUpper", "qcv", "Vmaxlen", 
                                    "DR")]
    fig.asp <- switch(slotstext,
                      "M" = 1.5,
                      "h" = 1,
                      "Linf" =1.5,
                      "L50" = 1,
                      "D" = 0.5,
                      "EffUpper" = 1/2,
                      "qcv" = 1,
                      "Vmaxlen"=1,
                      "DR" = 0.75)
    cat("```{r plot.", slotstext, ", echo=FALSE, fig.asp=", fig.asp, "}\n", append=TRUE, file=RMDfile, sep="")
    cat("plotSlot(OM, Pars, slot='", slotstext, "')\n", append=TRUE, file=RMDfile, sep="")
    cat("```\n\n\n", append=TRUE, file=RMDfile, sep="")   
    
  } else if ('Obs' %in% slots) {
    cat("\n### Obs Plots\n", append=TRUE, file=RMDfile, sep="")
    cat("```{r plot.Obs, echo=FALSE, fig.asp=1}\n", append=TRUE, file=RMDfile, sep="")
    cat("plotObs(OM)\n", append=TRUE, file=RMDfile, sep="")
    cat("```\n\n", append=TRUE, file=RMDfile, sep="")   
    
  } else if ("Imp" %in% slots) {
    cat("\n### Imp Plots\n", append=TRUE, file=RMDfile, sep="")
    cat("```{r plot.Imp, echo=FALSE, fig.asp=1}\n", append=TRUE, file=RMDfile, sep="")
    cat("plotImp(OM)\n", append=TRUE, file=RMDfile, sep="")
    cat("```\n\n", append=TRUE, file=RMDfile, sep="") 
  } 
}


plotSlot <- function(OM, Pars, slot) {
  if (OM@nsim > 48) OM@nsim <- 48
  if (slot == 'M') plotM2(OM, Pars) 
  if (slot == "h") plotRec(OM, Pars) 
  if (slot == "Linf") plotGrowth(OM, Pars) 
  if (slot == "L50") plotMat(OM, Pars) 
  if (slot == "D") plotDep(OM, Pars) 
  if (slot == "EffUpper") plotEff(OM, Pars) 
  if (slot == "qcv") plotQcv(OM, Pars)
  if (slot == "Vmaxlen") plotSelHists(OM, Pars)
  if (slot == "DR") plotSelect(OM, Pars)
  
}

plotSelHists <- function(OM, Pars, nsamp=3, col="darkgray", 
                         breaks=10, lwd=2) {
  
  
  ncol <- 3
  m <- layout(matrix(c(1,2,3,
                       4,5,6,
                       7,0,0), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 2, 2, 1), oma=c(3,3,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  its <- sample(1:OM@nsim, nsamp)
  
  hist2(Pars$L5, col=col, axes=FALSE, main="L5", breaks=breaks)
  abline(v=Pars$L5[OM@nyears,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$LFS, col=col, axes=FALSE, main="LFS", breaks=breaks)
  abline(v=Pars$LFS[OM@nyears,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$Vmaxlen, col=col, axes=FALSE, main="Vmaxlen", breaks=breaks)
  abline(v=Pars$Vmaxlen[OM@nyears,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  
  hist2(Pars$LR5, col=col, axes=FALSE, main="LR5", breaks=breaks)
  abline(v=Pars$LR5[OM@nyears,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$LFR, col=col, axes=FALSE, main="LFR", breaks=breaks)
  abline(v=Pars$LFR[OM@nyears,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$Rmaxlen, col=col, axes=FALSE, main="Rmaxlen", breaks=breaks)
  abline(v=Pars$Rmaxlen[OM@nyears,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$DR, col=col, axes=FALSE, main="DR", breaks=breaks)
  abline(v=Pars$DR[OM@nyears,its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
}
plotQcv <- function(OM, Pars, nsamp=3, col="darkgray", 
                    breaks=10, lwd=2) {
  
  ncol <- 2
  m <- layout(matrix(c(1,2,
                       3,3), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 2, 2, 1), oma=c(3,3,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  its <- sample(1:OM@nsim, nsamp)
  
  hist2(Pars$qinc, col=col, axes=FALSE, main="qinc", breaks=breaks)
  abline(v=Pars$qinc[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$qcv, col=col, axes=FALSE, main="qcv", breaks=breaks)
  abline(v=Pars$qcv[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  matplot(t(Pars$qvar), type="l", bty="l", xlab="Projection Years",
          ylab="Catchability", las=1, ylim=c(0, max(Pars$qvar)), xpd=NA, cex.lab=1.5)
  matplot(t(Pars$qvar)[,its], type="l", bty="l", xlab="", add=TRUE,
          ylab="", las=1, ylim=c(0, max(Pars$qvar)), lwd=3, lty=1)
  
}

plotEff <- function(OM, Pars, nsamp=3, col="darkgray", 
                    breaks=10, lwd=2) {
  
  ncol <- 3
  m <- layout(matrix(c(1,1,2), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 2, 2, 1), oma=c(3,3,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  its <- sample(1:OM@nsim, nsamp)
  
  yrs <- (OM@CurrentYr -  OM@nyears + 1) : OM@CurrentYr
  matplot(yrs, t(Pars$Find), type="l", bty="l", xlab="Historical Years",
          ylab="Fishing Effort", las=1, ylim=c(0, max(Pars$Find)), cex.lab=1.5, xpd=NA)
  matplot(t(Pars$Find)[,its], type="l", bty="l", xlab="", add=TRUE,
          ylab="", las=1, ylim=c(0, max(Pars$Find)), lwd=3, lty=1)
  
  hist2(Pars$Esd, col=col, axes=FALSE, main="Esd", breaks=breaks)
  abline(v=Pars$Esd[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
}





plotDep <- function(OM, Pars=NULL, nsim=48, nyears=50, proyears=50, nsamp=3, col="darkgray", 
                    breaks=10, lwd=2) {
  if (class(OM) != "OM") stop("Must supply object of class 'OM'")
  
  if (is.finite(OM@nyears)) nyears <- OM@nyears
  if (is.finite(OM@proyears)) proyears <- OM@proyears
  if (is.finite(OM@nsim)) nsim <- OM@nsim	
  
  if (is.null(Pars)) {
    OM <- updateMSE(OM) # update and add missing slots with default values
    out<- runMSE(OM,Hist=T)
    Pars <- c(out$SampPars, out$TSdata, out$MSYs)
  }
  
  its <- sample(1:nsim, nsamp)
  
  ncol <- 2
  m <- layout(matrix(c(1,2), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 2, 2, 1), oma=c(2,3,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  ssb0 <- matrix(rep(Pars$SSB0, nyears), nrow=nyears, byrow=TRUE)
  dep <- Pars$SSB/ssb0
  ylim <- c(0, max(dep))
  matplot(dep,  type="l", bty="l", ylab="SB/SB0", xlab="Historical Years", xpd=NA, ylim=ylim)
  matplot(dep[, its],  type="l", bty="l", ylab="", xlab="", add=TRUE, lwd=4, col=1:nsamp, 
          lty=1, ylim=ylim)
  
  hist2(Pars$D, col=col, axes=FALSE, main="Depletion (SB/SB0)", breaks=breaks)
  abline(v=Pars$D[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  
  
  
}



plotMat <- function(OM, Pars=NULL, nsim=48, nyears=50, proyears=50, nsamp=3, col="darkgray", breaks=10, lwd=2) {
  if (class(OM) != "OM") stop("Must supply object of class 'OM'")
  
  if (is.finite(OM@nyears)) nyears <- OM@nyears
  if (is.finite(OM@proyears)) proyears <- OM@proyears
  if (is.finite(OM@nsim)) nsim <- OM@nsim	
  
  if (is.null(Pars)) {
    OM <- updateMSE(OM) # update and add missing slots with default values
    out<- runMSE(OM,Hist=T)
    Pars <- c(out$SampPars, out$TSdata, out$MSYs)
  }
  
  its <- sample(1:nsim, nsamp)
  
  ncol <- 2
  m <- layout(matrix(c(1,2,3,4), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 2, 2, 1), oma=c(2,3,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  
  # Maturity 
  hist2(Pars$L50, col=col, axes=FALSE, main="L50", breaks=breaks)
  abline(v=Pars$L50[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  hist2(Pars$L95, col=col, axes=FALSE, main="L95", breaks=breaks)
  abline(v=Pars$L95[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  
  slope <- log(19)/(Pars$L95-Pars$L50)
  Ls <- seq(0, to=max(Pars$Linf), length.out=200)
  
  Mat_len <- sapply(its, function(X)  plogis(Ls, Pars$L50[X], 1/slope[X]))
  matplot(Ls, Mat_len, type="l", bty="l", main="Maturity-at-length", lwd=lwd, lty=1, 
          ylab="Probability", xlab="Length", ylim=c(0,1), xpd=NA)
  
  matplot(t(Pars$Mat_age[its,,nyears]), type="l", bty="l", main="Maturity-at-age", lwd=lwd, 
          lty=1, axes=FALSE, xlim=c(0, Pars$maxage), ylab="", xlab="Age", ylim=c(0,1), xpd=NA)
  axis(side=1)
  axis(side=2, labels=FALSE)
  
}



plotGrowth <- function(OM, Pars=NULL, nsim=48, nyears=50, proyears=50, nsamp=3, col="darkgray", 
                       breaks=10, lwd=2) {
  if (class(OM) != "OM") stop("Must supply object of class 'OM'")
  
  if (is.finite(OM@nyears)) nyears <- OM@nyears
  if (is.finite(OM@proyears)) proyears <- OM@proyears
  if (is.finite(OM@nsim)) nsim <- OM@nsim	
  
  if (is.null(Pars)) {
    OM <- updateMSE(OM) # update and add missing slots with default values
    out<- runMSE(OM,Hist=T)
    Pars <- c(out$SampPars, out$TSdata, out$MSYs)
  }
  
  its <- sample(1:nsim, nsamp)
  
  ncol <- 9
  m <- layout(matrix(c(1,1, 2,2, 0, 3,3, 4,4,
                       5,5, 6,6, 0, 7,7, 8,8,
                       9,9,9,9,  0, 10,10,10,10,
                       11,11,11, 12,12,12, 13,13,13,
                       14,14,14, 15,15,15, 16,16,16), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(5, 3, 3, 1), oma=c(4,6,2,1), las=0, no.readonly = TRUE)
  on.exit(par(op))
  
  # Histograms #### 
  # Linf
  hist2(Pars$Linf, col=col, axes=FALSE, main="Linf", breaks=breaks)
  abline(v=Pars$Linf[its], col=1:nsamp, lwd=lwd)
  axis(side=1)   
  
  # K 
  hist2(Pars$K, col=col, axes=FALSE, main="K", breaks=breaks)
  abline(v=Pars$K[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # t0 
  hist2(Pars$t0, col=col, axes=FALSE, main="t0", breaks=breaks)
  abline(v=Pars$t0[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # LenCV 
  hist2(Pars$LenCV, col=col, axes=FALSE, main="LenCV", breaks=breaks)
  abline(v=Pars$LenCV[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # Linfsd 
  hist2(Pars$Linfsd, col=col, axes=FALSE, main="Linfsd", breaks=breaks)
  abline(v=Pars$Linfsd[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # Linfgrad 
  hist2(Pars$Linfgrad, col=col, axes=FALSE, main="Linfgrad", breaks=breaks)
  abline(v=Pars$Linfgrad[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  
  # Ksd 
  hist2(Pars$Ksd, col=col, axes=FALSE, main="Ksd", breaks=breaks)
  abline(v=Pars$Ksd[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # Kgrad 
  hist2(Pars$Kgrad, col=col, axes=FALSE, main="Kgrad", breaks=breaks)
  abline(v=Pars$Kgrad[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  
  
  # By year 
  # Linf 
  matplot(t(Pars$Linfarray[its,]), type="l", bty="l", main="Linf by Year", lwd=lwd, lty=1, ylab="", xpd=NA) 
  # K 
  matplot(t(Pars$Karray[its,]), type="l", bty="l", main="K by Year", lwd=lwd, lty=1, ylab="", xpd=NA) 
  
  
  # Growth curves
  Len_age <- Pars$Len_age
  Wt_age <- Pars$Wt_age
  cex.lab <- 1.25
  fstYr <- Len_age[its,,1]
  curYr <- Len_age[its,,nyears]
  lstYr <- Len_age[its,,proyears+nyears]
  MaxL <- max(Len_age)
  matplot(t(fstYr), type="l", bty="l", main="First historical year", ylim=c(0, MaxL), 
          xlab="Age", ylab="Length", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)
  matplot(t(curYr), type="l", bty="l", main="Last historical year", ylim=c(0, MaxL),  
          axes=FALSE, xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)
  axis(side=1)
  axis(side=2, labels=FALSE)  
  matplot(t(lstYr), type="l", bty="l", main="Last projected year", ylim=c(0, MaxL), axes=FALSE, 
          xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)	
  axis(side=1)
  axis(side=2, labels=FALSE)  
  title("Sampled length-at-age curves", outer=TRUE, cex.main=2)
  
  fstYr <- Wt_age[its,,1]
  curYr <- Wt_age[its,,nyears]
  lstYr <- Wt_age[its,,proyears+nyears]
  MaxL <- max(Wt_age)
  matplot(t(fstYr), type="l", bty="l", main="First historical year", ylim=c(0, MaxL), 
          xlab="Age", ylab="Weight", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)
  matplot(t(curYr), type="l", bty="l", main="Last historical year", ylim=c(0, MaxL), 
          axes=FALSE, xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)
  axis(side=1)
  axis(side=2, labels=FALSE)  
  matplot(t(lstYr), type="l", bty="l", main="Last projected year", ylim=c(0, MaxL), 
          axes=FALSE, xlab="Age", ylab="", cex.lab=cex.lab, lwd=lwd, lty=1, xpd=NA)	
  axis(side=1)
  axis(side=2, labels=FALSE)  
  title("Sampled length-at-age curves", outer=TRUE, cex.main=2)
  
}


plotRec <- function(OM, Pars=NULL, nsim=48, nyears=50, proyears=50, nsamp=3, col="darkgray", breaks=10, lwd=2) {
  if (class(OM) != "OM") stop("Must supply object of class 'OM'")
  
  if (is.finite(OM@nyears)) nyears <- OM@nyears
  if (is.finite(OM@proyears)) proyears <- OM@proyears
  if (is.finite(OM@nsim)) nsim <- OM@nsim	
  
  if (is.null(Pars)) {
    OM <- updateMSE(OM) # update and add missing slots with default values
    out<- runMSE(OM,Hist=T)
    Pars <- c(out$SampPars, out$TSdata, out$MSYs)
  }
  
  its <- sample(1:nsim, nsamp)
  
  ncol <- 3
  m <- layout(matrix(c(1, 2, 3, 
                       4, 5, 0,
                       6, 6, 6), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 1, 3, 1), oma=c(3,4,2,1), las=0, no.readonly = TRUE)
  on.exit(par(op))
  
  ## Histograms ####
  
  # h 
  hist2(Pars$hs, col=col, axes=FALSE, main="h", breaks=breaks)
  abline(v=Pars$hs[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # Perr 
  hist2(Pars$procsd, col=col, axes=FALSE, main="Perr", breaks=breaks)
  abline(v=Pars$procsd[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # AC 
  hist2(Pars$AC, col=col, axes=FALSE, main="AC", breaks=breaks)
  abline(v=Pars$AC[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  # recgrad - not used 
  # hist2(Pars$recgrad, col=col, axes=FALSE, main="recgrad- currently not used", breaks=breaks)
  # abline(v=Pars$recgrad[its], col=1:nsamp, lwd=lwd)
  # axis(side=1)
  
  # Period - if used
  if (!is.null(Pars$Period)) {
    hist2(Pars$Period, col=col, axes=FALSE, main="Period", breaks=breaks)
    abline(v=Pars$Period[its], col=1:nsamp, lwd=lwd)
    axis(side=1)  
  } else {
    hist2(0, col=col, axes=FALSE, main="Period", breaks=breaks)
  }
  
  # Amplitude
  if (!is.null(Pars$Amplitude)) {
    hist2(Pars$Amplitude, col=col, axes=FALSE, main="Amplitude", breaks=breaks)
    abline(v=Pars$Amplitude[its], col=1:nsamp, lwd=lwd)
    axis(side=1)  
  } else {
    hist2(0, col=col, axes=FALSE, main="Amplitude", breaks=breaks)
  }
  
  
  # Recruitment
  matplot(t(Pars$Perr[its,]), type="l", bty="l", main="Rec Devs by Year", lwd=lwd, lty=1, ylab="")
  
  
}





# Plot Natural Mortality 
plotM2 <- function(OM, Pars=NULL, nsim=48, nyears=50, proyears=50, nsamp=3, col="darkgray", breaks=10, lwd=2) {
  if (class(OM) != "OM") stop("Must supply object of class 'OM'")
  
  if (is.finite(OM@nyears)) nyears <- OM@nyears
  if (is.finite(OM@proyears)) proyears <- OM@proyears
  if (is.finite(OM@nsim)) nsim <- OM@nsim	
  
  if (is.null(Pars)) {
    OM <- updateMSE(OM) # update and add missing slots with default values
    out<- runMSE(OM,Hist=T)
    Pars <- c(out$SampPars, out$TSdata, out$MSYs)
  }
  
  its <- sample(1:nsim, nsamp)
  
  
  ncol <- 4
  m <- layout(matrix(c(1, 2, 3, 4,
                       rep(5, 4),
                       6, 7, 8, 0,
                       9, 10, 11, 0,
                       12, 13, 14, 0), ncol=ncol, byrow = TRUE))
  
  op <- par(mar = c(3, 1, 3, 1), oma=c(3,4,2,1), las=1, no.readonly = TRUE)
  on.exit(par(op))
  
  ## Plot histograms of M parameters ####
  hist2(Pars$M, col=col, axes=FALSE, main="M", breaks=breaks)
  abline(v=Pars$M[its], col=1:nsamp, lwd=lwd)
  axis(side=1)  
  
  hist2(Pars$Mexp, col=col, axes=FALSE, main="Mexp", breaks=breaks)
  abline(v=Pars$Mexp[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  hist2(Pars$Msd, col=col, axes=FALSE, main="Msd", breaks=breaks)
  abline(v=Pars$Msd[its], col=1:nsamp, lwd=lwd)
  axis(side=1) 
  
  hist2(Pars$Mgrad, col=col, axes=FALSE, main="Mgrad", breaks=breaks)
  abline(v=Pars$Mgrad[its], col=1:nsamp, lwd=lwd)
  axis(side=1)
  
  # M by year 
  ylims <- range(Pars$M_ageArray[its,, ]) * c(0.95, 1.05)
  ylims[1] <- min(0, ylims[1] )
  matplot(t(Pars$Marray[its,]), type="l", lty=1, bty="l", main="average adult M by Year", lwd=lwd, ylab="M", ylim=ylims)
  abline(v=nyears, col="gray", lty=2)
  text(nyears, min(Pars$Marray[its,]), "Last historical year", pos=4, col="gray")
  
  # M at age 
  M_ageArray <- Pars$M_ageArray
  Len_age <- Pars$Len_age
  Wt_at_age <- Pars$Wt_age
  
  
  matplot(t(M_ageArray[its,,1]), type="l", lty=1, bty="l", lwd=lwd, ylim=ylims, ylab="M")
  mtext(side=3, "First historical year", cex=0.8, line=-1)
  mtext(side=1, "Age", line=2, cex=0.7)
  
  matplot(t(M_ageArray[its,,nyears]), type="l", lty=1, bty="l", main="M-at-age", lwd=lwd, ylim=ylims, axes=FALSE, ylab="")
  mtext(side=3, "Last historical year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Age", line=2, cex=0.7)
  
  matplot(t(M_ageArray[its,,nyears+proyears]), type="l", lty=1, bty="l", lwd=lwd, ylim=ylims, axes=FALSE, ylab="")
  mtext(side=3, "Last projected year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Age", line=2, cex=0.7)
  
  
  # M at length 
  xlims <- range(Len_age[its,, c(1, nyears, nyears+proyears)]) * c(0.95, 1.05)
  
  matplot(t(Len_age[its,,1]), t(M_ageArray[its,,1]), type="l", lty=1, bty="l", lwd=lwd, 
          ylim=ylims, xlim=xlims, ylab="M", xlab="")
  mtext(side=3, "First historical year", cex=0.8, line=-1)
  mtext(side=1, "Length", line=2, cex=0.7)
  
  matplot(t(Len_age[its,,nyears]), t(M_ageArray[its,,nyears]), type="l", lty=1, bty="l", 
          main="M-at-length", lwd=lwd, ylim=ylims, xlim=xlims, axes=FALSE, ylab="", xlab="")
  mtext(side=3, "Last historical year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Length", line=2, cex=0.7)
  
  matplot(t(Len_age[its,,nyears+proyears]), t(M_ageArray[its,,nyears+proyears]), type="l", 
          lty=1, bty="l", lwd=lwd, ylim=ylims, axes=FALSE, xlim=xlims, ylab="", xlab="")
  mtext(side=3, "Last projected year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Length", line=2, cex=0.7)
  
  
  # M at weight
  xlims <- range(Wt_at_age[its,, c(1, nyears, nyears+proyears)]) * c(0.95, 1.05)
  
  matplot(t(Wt_at_age[its,,1]), t(M_ageArray[its,,1]), type="l", lty=1, bty="l", lwd=lwd, 
          ylim=ylims, xlim=xlims, ylab="M", xlab="")
  mtext(side=3, "First historical year", cex=0.8, line=-1)
  mtext(side=1, "Weight", line=2, cex=0.7)
  
  matplot(t(Wt_at_age[its,,nyears]), t(M_ageArray[its,,nyears]), type="l", lty=1, bty="l", 
          main="M-at-weight", lwd=lwd, ylim=ylims, xlim=xlims, axes=FALSE, ylab="", xlab="")
  mtext(side=3, "Last historical year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Weight", line=2, cex=0.7)
  
  matplot(t(Wt_at_age[its,,nyears+proyears]), t(M_ageArray[its,,nyears+proyears]), type="l", 
          lty=1, bty="l", lwd=lwd, ylim=ylims, axes=FALSE, xlim=xlims, ylab="", xlab="")
  mtext(side=3, "Last projected year", cex=0.8, line=-1)
  axis(side=1)
  mtext(side=1, "Weight", line=2, cex=0.7)
  
  
}



## Old functions to be removed at some stage if not used/needed ####

# #' Read in operating model parameters from Excel spreadsheet
# #' 
# #' A function to read in operating model parameters from an Excel spreadsheet
# #' with tabs named following specific convention. Since DLMtool 4.5 this function is no
# #' longer recommended. Use 'OMinit' instead.
# #' 
# #' The Excel spreadsheet must have tabs named with the following convention.
# #' For example if \code{stkname} is 'myFish', the Stock parameters are in a tab
# #' named 'myFishStock', Fleet parameters in a tab named 'myFishFleet', 
# #' Observation parameters in a tab named 'myFishObs', and Implementation in 'myFishImp'. 
# #' All tabs (Stock, Fleet, Obs, and Imp) must be present for a single stock. You can have multiple
# #' stocks in a single spreadsheet, provided that the stock names are different.
# #' 
# #' @usage OM_xl(fname, stkname, fpath = '', saveCSV = FALSE)
# #' @param fname Name of the Excel spreadsheet file. Must include file
# #' extension.
# #' @param stkname Name of the Stock. Only required if more than one Stock in the Excel file.
# #' @param fpath Full file path, if file is not in current working directory
# #' @param saveCSV Do you also want to save the Stock, Fleet and Observation
# #' parameters to CSV files?
# #' @return A object of class OM
# #' @author A. Hordyk
# #' @examples
# #' 
# #' \dontrun{
# #' OM <- OM_xl(fname='OMTables.xlsx', stkname='myFish')
# #' }
# #' 
# #' @export 
# OM_xl <- function(fname, stkname=NULL, fpath = "", saveCSV = FALSE) {
#   .Deprecated('XL2OM')
#   infile <- paste0(fpath, fname)  # full path and name 
#   shtname <- readxl::excel_sheets(infile)  # names of the sheets 
#  
#   if (is.null(stkname)) {
#     names <- c(unlist(strsplit(shtname[grep('Stock', shtname)], "Stock")),
#                unlist(strsplit(shtname[grep('Fleet', shtname)], "Fleet")),
#                unlist(strsplit(shtname[grep('Obs', shtname)], "Obs")))
#     if (length(unique(names)) == 1) {
#       stkname <- unique(names)
#     } else stop("stkname not provided and multiple stocks found in spreadsheet.")
#   }
# 
#   
#   # Stock
#   sheet <- grep(paste0(stkname, "Stock"), shtname)
#   if(length(sheet)<1) stop("No Stock sheet found. Looking for: ", paste0(stkname, "Stock"))
#   stock <- readxl::read_excel(infile, sheet = sheet, col_names = FALSE)
#   stock <- as.data.frame(stock) 
#   if (all(dim(stock) == 0)) stop("No data found in Stock tab")
#   tmpfile <- paste0(fpath, stkname, "Stock.csv")
#   if (file.exists(tmpfile)) unlink(tmpfile)
#   writeCSV(inobj = stock, tmpfile, objtype = "Stock")
#   tmpstock <- new("Stock", tmpfile)
#   if (!saveCSV) unlink(tmpfile)
#   
#   # Fleet
#   index <- which(pmatch(shtname, paste0(stkname, "Fleet")) == 1)
#   if (length(index) > 1) stop("More than one match")
#   
#   if(length(index)<1) stop("No Fleet sheet found. Looking for: ", paste0(stkname, "Fleet"))
#   
#   fleet <- readxl::read_excel(infile, sheet = index, col_names = FALSE)
#   fleet <- as.data.frame(fleet)
#   if (all(dim(fleet) == 0)) stop("No data found in Fleet tab")
#   tmpfile <- paste0(fpath, stkname, "Fleet.csv")
#   if (file.exists(tmpfile)) unlink(tmpfile)
#   writeCSV(inobj = fleet, tmpfile, objtype = "Fleet")
#   tmpfleet <- new("Fleet", tmpfile)
#   if (!saveCSV) unlink(tmpfile)
#   
#   # Observation
#   index <- which(pmatch(shtname, paste0(stkname, "Obs")) == 1)
#   if (length(index) > 1) stop("More than one match")
#   if(length(index)<1) stop("No Obs sheet found. Looking for: ", paste0(stkname, "Obs")) 
#   obs <- readxl::read_excel(infile, sheet = index, col_names = FALSE)
#   obs <- as.data.frame(obs)
#   if (all(dim(obs) == 0)) stop("No data found in Obs tab")
#   tmpfile <- paste0(fpath, stkname, "Obs.csv")
#   if (file.exists(tmpfile)) unlink(tmpfile)
#   writeCSV(inobj = obs, tmpfile, objtype = "Obs")
#   tmpobs <- new("Obs", tmpfile)
#   if (!saveCSV) unlink(tmpfile)
#  
#   # Implementation
#   index <- which(pmatch(shtname, paste0(stkname, "Imp")) == 1)
#   if (length(index) > 1)  stop("More than one match")
#   if(length(index)<1) {
#     warning("No Imp sheet found. Looking for: ", paste0(stkname, "Imp"), ". Defaulting to 'Perfect_Imp'", call.=FALSE) 
#     tmpimp <- DLMtool::Perfect_Imp
#   } else {
#     imp <- readxl::read_excel(infile, sheet = index, col_names = FALSE)
#     imp <- as.data.frame(imp)
#     if (all(dim(imp) == 0)) stop("No data found in Imp tab")
#     tmpfile <- paste0(fpath, stkname, "Imp.csv")
#     if (file.exists(tmpfile)) unlink(tmpfile)
#     writeCSV(inobj = imp, tmpfile, objtype = "Imp")
#     tmpimp <- new("Imp", tmpfile)
#     if (!saveCSV) unlink(tmpfile)
#   }
# 
#   
#   # Operating Model
#   OM <- new("OM", Stock = tmpstock, Fleet = tmpfleet, Obs = tmpobs, Imp=tmpimp)
#   OM
# }
# 

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

# #' Read in feasibility parameters from Excel spreadsheet
# #' 
# #' A function to read in feasibility parameters from an Excel spreadsheet with
# #' tabs named following specific convention
# #' 
# #' The Excel spreadsheet must have tabs named with the following convention.
# #' For example if \code{stkname} is 'myFish', the tab must be named
# #' 'myFishFease,
# #' 
# #' @usage Fease_xl(fname, stkname, fpath = '', saveCSV = FALSE)
# #' @param fname Name of the Excel spreadsheet file. Must include file
# #' extension.
# #' @param stkname Name of the Stock.
# #' @param fpath Full file path, if file is not in current working directory
# #' @param saveCSV Do you also want to save the Stock, Fleet and Observation
# #' parameters to CSV files?
# #' @return A object of class Fease
# #' @author A. Hordyk
# #' @examples
# #' 
# #'  \dontrun{
# #'  myFease <- Fease_xl(fname='FeaseTables.xlsx', stkname='myFish')
# #' }
# #' 
# #' @export Fease_xl
# Fease_xl <- function(fname, stkname, fpath = "", saveCSV = FALSE) {
#   infile <- paste0(fpath, fname)  # full path and name 
#   shtname <- readxl::excel_sheets(infile)  # names of the sheets 
#   # Fease
#   feasedat <- readxl::read_excel(infile, sheet = grep(paste0(stkname, "Fease"), 
#                                                       shtname), col_names = FALSE)
#   feasedat <- feasedat[, 1:2]
#   tmpfile <- paste0(fpath, stkname, "Fease.csv")
#   if (file.exists(tmpfile)) 
#     unlink(tmpfile)
#   writeCSV(inobj = feasedat, tmpfile, objtype = "Fease")
#   fease <- new("Fease", tmpfile)
#   if (!saveCSV) 
#     unlink(tmpfile)
#   
#   fease
# }
# 


#' Internal function to write CSVs for objects
#' 
#' Used internally in the DLMtool package to write CSV files from an existing
#' DLMtool object
#' 
#' 
#' @param inobj A object of class Stock, Fleet, Obs, Imp, Data, or OM
#' 
#' @param tmpfile The full file path and name for the saved CSV file
#' @param objtype The class corresonding to the \code{inobj}
#' @author A. Hordyk
writeCSV <- function(inobj, tmpfile = NULL, objtype = c("Stock", "Fleet", 
                                                        "Obs", "Imp", "Data", "OM")) {
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


