
# Still to add Help Manual for these 
OM_xl <- function(fname, stkname, fpath="", saveCSV=FALSE) {
  infile <- paste0(fpath, fname)	# full path and name 
  shtname <- excel_sheets(infile) # names of the sheets 
  # Stock 
  stock <- read_excel(infile, sheet=grep(paste0(stkname, "Stock"), shtname),
			col_names=FALSE)
  stock <- as.data.frame(stock)
  tmpfile <- paste0(fpath, stkname, "Stock.csv")
  if (file.exists(tmpfile)) unlink(tmpfile) 
  writeCSV(inobj=stock, tmpfile, objtype="Stock")
  tmpstock <- new("Stock", tmpfile)
  if (!saveCSV) unlink(tmpfile)

  # Fleet 
  fleet <- read_excel(infile, sheet=grep(paste0(stkname, "Fleet"), shtname),
			col_names=FALSE)	
  fleet <- as.data.frame(fleet)			
  tmpfile <- paste0(fpath, stkname, "Fleet.csv")
  if (file.exists(tmpfile)) unlink(tmpfile) 
  writeCSV(inobj=fleet, tmpfile, objtype="Fleet")
  tmpfleet <- new("Fleet", tmpfile)
  if (!saveCSV) unlink(tmpfile)  

  # Observation 
  index <- which(pmatch(shtname, paste0(stkname, "Observation")) == 1)
  if (length(index) > 1) stop("More than one match")
  obs <- read_excel(infile, sheet=index, col_names=FALSE)
  obs <- as.data.frame(obs)		
  tmpfile <- paste0(fpath, stkname, "Observation.csv")
  if (file.exists(tmpfile)) unlink(tmpfile) 
  writeCSV(inobj=obs, tmpfile, objtype="Observation")
  tmpobs <- new("Observation", tmpfile)
  if (!saveCSV) unlink(tmpfile)  
  
  # Operating Model 
  OM <- new("OM",  Stock=tmpstock,  Fleet=tmpfleet,  Observation=tmpobs)
  OM 
}

Fease_xl <- function(fname, stkname, fpath="", saveCSV=FALSE) {
  infile <- paste0(fpath, fname)	# full path and name 
  shtname <- excel_sheets(infile) # names of the sheets 
  # Fease
  feasedat <- read_excel(infile, sheet=grep(paste0(stkname, "Fease"), shtname),
			col_names=FALSE)
  feasedat <- feasedat[,1:2]			
  tmpfile <- paste0(fpath, stkname, "Fease.csv")
  if (file.exists(tmpfile)) unlink(tmpfile) 
  writeCSV(inobj=feasedat, tmpfile, objtype="DLM_fease")
  fease <- new("DLM_fease", tmpfile)
  if (!saveCSV) unlink(tmpfile)

  fease
}

writeCSV <- function(inobj, tmpfile=NULL,
  objtype=c("Stock", "Fleet", "Observation",  "DLM_data", "OM", "DLM_fease")) {
  objtype <- match.arg(objtype)
  tmpobj <- new(objtype)
  sn <- slotNames(tmpobj)
  ind <- which(inobj[,1] %in% sn == FALSE)
  if (length(ind) > 0) {
    message("Input file names don't match slot names for ", objtype, " object")
	message("Unknown input name:", inobj[ind,1])
    stop("Check the input file row names")
  }
  for (X in seq_along(sn)) {
    ind <- match(sn[X], inobj[,1])
	if (!is.na(ind)) {
	  indat <- inobj[ind,]
	  index <- which(!is.na(indat))
	  index <- 2:max(index)
      if (X == 1) write(do.call(paste, c(sn[X], as.list(indat[index]), sep=","))
	    ,tmpfile,1)
	  if (X > 1)  write(do.call(paste, c(sn[X], as.list(indat[index]), sep=",")),
	    tmpfile,1, append=TRUE) 
	} 
  }
}


