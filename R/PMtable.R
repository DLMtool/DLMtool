
#' Create a table of Performance Limits
#'
#' @param MSE An object of class 'MSE'
#' @param ... PM objects to be used as performance limits. Characters (i.e names of PM objects)
#' @param Pthresh Minimum probability threshold
#' @param out.file Name of the output file. If none provided, output file will be named 'PerfLimTable' 
#' @param output_format Output file format. Currently only 'html_document' is supported
#' @param openFile Logical. Should the file be opened in browser?
#' @param quiet Logical. An option to suppress printing of the pandoc command line.
#' @param dir Optional. Directory for output file. Default is working directory.
#' @param RMDfile Optional. RMD template file
#' @param font_size Numeric. Font size for text in the table
#' @param full_width Logical. Should table be full width?
#' @param enableSearch Logical. Should search be enabled in the html table?
#' @param PMlist Optional. List of PM names.
#'
#' @return Invisibly returns names of MPs that pass all performance limits
#' @export
#'
PMLimit <- function(MSE, ..., Pthresh=NULL, out.file=NULL,
                      output_format="html_document", openFile=TRUE,
                      quiet=TRUE, dir=NULL, RMDfile=NULL, font_size=14,
                      full_width=TRUE, enableSearch=TRUE, PMlist=NULL) {

  if (class(MSE) != 'MSE') stop("Object is not class 'MSE'", call. = FALSE)
  nMPs <- MSE@nMPs
  if (is.null(Pthresh)) stop("Must specify 'Pthresh'", call. = FALSE)

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
             Prob=lapply(runPM, function(x) x@Mean) %>% unlist(),
             PM=rep(1:nPM, each=nMPs))
  df$Prob <- round(df$Prob,2)
  temp <- df %>% group_by(MP) %>% summarize(min=min(Prob))
  df <- dplyr::left_join(df, temp, by='MP') %>% arrange(MP)

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
  if (!is.null(out.file))
    out.file <- tools::file_path_sans_ext(out.file)

  if (is.null(out.file)) out.file <- paste0('PerfLimTable', ext)
  if (is.null(dir)) dir <- getwd()
  RMDfileout <- file.path(dir, out.file)

  params <- list(df=df, runPM=runPM, Name=MSE@Name, font_size=font_size,
                 full_width=full_width, output_format=output_format,
                 enableSearch=enableSearch, Pthresh=Pthresh)
  knitr::knit_meta(class=NULL, clean = TRUE)

  if(is.null(RMDfile))
    RMDfile <- file.path(system.file(package = "DLMtool"), 'PLimitTable.Rmd')

  rmarkdown::render(input=RMDfile, output_file=RMDfileout, output_format=output_format,
                    output_dir=dir, param=params, quiet=quiet)
  if (openFile) utils::browseURL(RMDfileout)

  invisible(df %>% filter(min >= Pthresh) %>% select(MP) %>% unique())
}




