
#' Create a table of Performance Limits and Performance Objectives
#'
#' @param MSE An object of class 'MSE'
#' @param ... PM objects to be used as performance limits. Characters (i.e names of PM objects)
#' @param Prob Minimum probability threshold
#' @param Labels Optional named list specifying new labels for MPs. For example: `Labels = list(AvC="Average Catch", CC1="Constant Catch")`
#' @param FeaseMPs Optional. Character vector of MP names that are considered feasible. e.g. the output from `Fease()`
#' @param out.file Name of the output file. If none provided, output file will be named 'PerfLimTable' 
#' @param output_format Output file format. Currently only 'html_document' is supported
#' @param openFile Logical. Should the file be opened in browser?
#' @param quiet Logical. An option to suppress printing of the pandoc command line.
#' @param dir Optional. Directory for output file. Default is working directory.
#' @param RMDfile Optional. RMD template file
#' @param font_size Numeric. Font size for text in the table
#' @param auto_width Logical. Should table be width be automatic?
#' @param enableSearch Currently disabled. Logical. Should search be enabled in the html table? 
#' @param PMlist Optional. List of PM names.
#' @param build Logical. Build the html table?
#' @describeIn PMLimit Create a table of Performance Limits 
#' @author A. Hordyk
#' @return `PMLimit` invisibly returns names of MPs that pass all performance limits 
#' @examples
#' \dontrun{
#' MSE <- runMSE()
#' PMLimit(MSE, "P50", "PNOF", Prob=0.9)
#' PMObj(MSE, "P100", "LTY")
#' }
#' 
#' @export
#'
PMLimit <- function(MSE, ..., Prob=NULL, Labels=NULL, FeaseMPs=NULL,
                    out.file=NULL,
                      output_format="html_document", openFile=TRUE,
                      quiet=TRUE, dir=NULL, RMDfile=NULL, font_size=14,
                      auto_width=FALSE, enableSearch=TRUE, PMlist=NULL, build=TRUE) {

  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("DT is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("kableExtra is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package \"knitr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package \"rmarkdown\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package \"tidyr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  MP <- prob <- Feasible <- NULL # hacks for cran checks
  if (class(MSE) != 'MSE') stop("Object is not class 'MSE'", call. = FALSE)
  nMPs <- MSE@nMPs
  if (is.null(Prob)) stop("Must specify 'Prob'", call. = FALSE)

  # Calculate prob of performace limits
  if (is.null(PMlist)) {
    PMlist <- unlist(list(...))
  } else {
    PMlist <- unlist(PMlist)
  }
  if (length(PMlist)<1) stop("Must provided at least one PM", call.=FALSE)
  
  nPM <- length(PMlist)
  for (X in seq_along(PMlist))
    if (!PMlist[X] %in% avail("PM")) stop(PMlist[X], " is not a valid PM method")
  runPM <- vector("list", length(PMlist))
  for (X in 1:length(PMlist)) {
    runPM[[X]] <- eval(call(PMlist[X], MSE))
  }

  # Create data frame of probs
  df <- data.frame(MP=lapply(runPM, function(x) x@MPs) %>% unlist(),
             prob=lapply(runPM, function(x) x@Mean) %>% unlist(),
             PM=rep(1:nPM, each=nMPs))
  df$prob <- round(df$prob,2)
  temp <- df %>% dplyr::group_by(MP) %>% dplyr::summarize(min=min(prob))
  df <- dplyr::left_join(df, temp, by='MP') %>% dplyr::arrange(MP)
  df$MP <- as.character(df$MP)
  df$url <- sapply(df$MP, MPurl) %>% unlist()
  types <- MPtype(df$MP)
  df$Type <- NA
  ind <- match(df$MP, types[,1])
  df$Type <- types[ind,2]
  
  labels <- MSE@MPs
  if (class(Labels) == "list") {
    repnames <- names(Labels)
    invalid <- repnames[!repnames %in% labels]
    if (length(invalid >0)) {
      warning("Labels: ", paste(invalid, collapse=", "), " are not MPs in MSE")
      Labels[invalid] <- NULL
      repnames <- names(Labels)
    }
    labels[labels %in% repnames] <- Labels %>% unlist()
  }
  
  
  df$Feasible <- NA
  if(!is.null(FeaseMPs)) {
    if(!class(FeaseMPs) == "character") stop("FeaseMPs must be character vector", call.=FALSE)
    # df$Feasible <- rep(MSE@MPs %in% FeaseMPs, each=2)
    df$Feasible <- df$MP %in% FeaseMPs
    df$Feasible[is.na(df$Feasible)] <- FALSE
  }
  df$Feasible[df$Feasible != TRUE] <- "No"
  df$Feasible[df$Feasible == TRUE] <- "Yes"
  df$MP <- labels[match(df$MP,MSE@MPs)]
 
  if (output_format == "html_document") {
    ext <- '.html'
  } else if (output_format == "pdf_document") {
    ext <- '.pdf'
  # } else if (output_format == "word_document") {
  #   ext <- '.docx'
  # } else if (output_format == "odt_document") {
  #   ext <- '.odt'
  } else {
    stop("output_format '", output_format, "' is not valid. Available options are: 'html_document', 'pdf_document'", call.=FALSE)
  }
  if (!is.null(out.file)) {
    out.file <- tools::file_path_sans_ext(out.file)
    out.file <- paste0(out.file, ext)
  }
    
  if (is.null(out.file)) out.file <- paste0('PerfLimTable', ext)
  if (is.null(dir)) dir <- getwd()
  RMDfileout <- file.path(dir, out.file)

  params <- list(df=df, runPM=runPM, Name=MSE@Name, font_size=font_size,
                 auto_width=auto_width, output_format=output_format,
                 enableSearch=enableSearch, Prob=Prob)
  knitr::knit_meta(class=NULL, clean = TRUE)

  if(is.null(RMDfile))
    RMDfile <- file.path(system.file(package = "DLMtool"), 'PLimitTable.Rmd')

  out <- df %>% dplyr::filter(min >= Prob & Feasible=="Yes") %>% dplyr::select(MP) %>% unique()
  if (build) {
    if (file.exists(RMDfileout)) unlink(RMDfileout)
    rmarkdown::render(input=RMDfile, output_file=RMDfileout, output_format=output_format,
                      output_dir=dir, param=params, quiet=quiet)
    if (file.exists(RMDfileout)) message("Table successfully built: ", RMDfileout)
    if (!file.exists(RMDfileout)) {
      warning("Rmarkdown file could not be rendered. Run with quiet=FALSE for error messages")
      return(invisible(out))
    }
    if (openFile) utils::browseURL(RMDfileout)
  }  
  return(invisible(out))

}

#' @describeIn PMLimit Create a table of Performance Objectives. 
#' @param use.colors Logical. Color scale the probability text?
#' @param cols Optional character vector of colors for probability text
#' @param show.legend Logical. Show the legend??
#' @param cex.tex Size of legend text
#' @param inc.title Logical. Include title for legend?
#' @param title Title for the legend
#' @export
#'
PMObj <- function(MSE, ..., Labels=NULL, out.file=NULL,
                  output_format="html_document", openFile=TRUE,
                  quiet=TRUE, dir=NULL, RMDfile=NULL, font_size=14,
                  use.colors=TRUE,
                  cols=NULL, show.legend=TRUE,
                  auto_width=FALSE, enableSearch=TRUE, PMlist=NULL, build=TRUE,
                  cex.tex=0.75, inc.title=TRUE, title="Legend") {
  
  MP <- prob <- NULL # hacks for cran checks
  
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("DT is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop("kableExtra is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package \"knitr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package \"rmarkdown\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package \"tidyr\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (class(MSE) != 'MSE') stop("Object is not class 'MSE'", call. = FALSE)
  nMPs <- MSE@nMPs

  # Calculate prob of performace limits
  if (is.null(PMlist)) {
    PMlist <- unlist(list(...))
  } else {
    PMlist <- unlist(PMlist)
  }
  if(length(PMlist) == 0) PMlist <- c("STY", "LTY", "P10", "AAVY")

  nPM <- length(PMlist)
  for (X in seq_along(PMlist))
    if (!PMlist[X] %in% avail("PM")) stop(PMlist[X], " is not a valid PM method")
  runPM <- vector("list", length(PMlist))
  for (X in 1:length(PMlist)) {
    runPM[[X]] <- eval(call(PMlist[X], MSE))
  }
  
  # Create data frame of probs
  df <- data.frame(MP=lapply(runPM, function(x) x@MPs) %>% unlist(),
                   prob=lapply(runPM, function(x) x@Mean) %>% unlist(),
                   PM=rep(1:nPM, each=nMPs))
  df$prob <- round(df$prob,2)
  temp <- df %>% dplyr::group_by(MP) %>% dplyr::summarize(min=min(prob))
  df <- dplyr::left_join(df, temp, by='MP') %>% dplyr::arrange(MP)
  df$MP <- as.character(df$MP)
  df$url <- sapply(df$MP, MPurl) %>% unlist()
  types <- MPtype(df$MP)
  df$Type <- NA
  ind <- match(df$MP, types[,1])
  df$Type <- types[ind,2]
  
  labels <- MSE@MPs
  if (class(Labels) == "list") {
    repnames <- names(Labels)
    invalid <- repnames[!repnames %in% labels]
    if (length(invalid >0)) {
      warning("Labels: ", paste(invalid, collapse=", "), " are not MPs in MSE")
      Labels[invalid] <- NULL
      repnames <- names(Labels)
    }
    labels[labels %in% repnames] <- Labels %>% unlist()
  }
  
  df$MP <- labels[match(df$MP,MSE@MPs)]

  if (output_format == "html_document") {
    ext <- '.html'
  } else if (output_format == "pdf_document") {
    ext <- '.pdf'
    # } else if (output_format == "word_document") {
    #   ext <- '.docx'
    # } else if (output_format == "odt_document") {
    #   ext <- '.odt'
  } else {
    stop("output_format '", output_format, "' is not valid. Available options are: 'html_document', 'pdf_document'", call.=FALSE)
  }
  if (!is.null(out.file)) {
    out.file <- tools::file_path_sans_ext(out.file)
    out.file <- paste0(out.file, ext)
  }
  
  if (is.null(out.file)) out.file <- paste0('PerfObjTable', ext)
  if (is.null(dir)) dir <- getwd()
  RMDfileout <- file.path(dir, out.file)
  
  if (is.null(cols)) {
    colsfun <- colorRampPalette(c("forestgreen", "orange", "red"))
    cols <- rev(colsfun(5))
  }
  params <- list(df=df, runPM=runPM, Name=MSE@Name, font_size=font_size,
                 auto_width=auto_width, output_format=output_format,
                 enableSearch=enableSearch, cols=cols, use.colors=use.colors,
                 show.legend=show.legend)
  knitr::knit_meta(class=NULL, clean = TRUE)
  
  if(is.null(RMDfile))
    RMDfile <- file.path(system.file(package = "DLMtool"), 'PObjTable.Rmd')
  
  if (build) {
    if (file.exists(RMDfileout)) unlink(RMDfileout)
    rmarkdown::render(input=RMDfile, output_file=RMDfileout, output_format=output_format,
                      output_dir=dir, param=params, quiet=quiet)
    if (file.exists(RMDfileout)) message("Table successfully built: ", RMDfileout)
    if (!file.exists(RMDfileout)) {
      warning("Rmarkdown file could not be rendered. Run with quiet=FALSE for error messages")
    }
    if (openFile) utils::browseURL(RMDfileout)
  }  

}




