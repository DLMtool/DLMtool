# 
# MSEobj <- runMSE()
# 
# dom <- DOM(MSEobj)
# 
# PMlist <- avail("PM")

Dom <- function(MSEobj, ..., PMlist=NULL, Refs=NULL, Yrs=NULL) {
  if (class(MSEobj) != 'MSE') stop("Object must be class `MSE`", call.=FALSE)
  if (is.null(PMlist)) {
    PMlist <- unlist(list(...))
  } else {
    PMlist <- unlist(PMlist)
  }
  
  if (class(PMlist) != 'character') stop("Must provide names of PM methods")
  
  runPM <- vector("list", length(PMlist))
  for (X in 1:length(PMlist)) {
    ref <- Refs[[PMlist[X]]]
    yrs <- Yrs[[PMlist[X]]]
    if (is.null(ref)) {
      if (is.null(yrs)) {
        runPM[[X]] <- eval(call(PMlist[X], MSEobj))    
      } else {
        runPM[[X]] <- eval(call(PMlist[X], MSEobj, Yrs=yrs))  
      }
    } else {
      if (is.null(yrs)) {
        runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref=ref))    
      } else {
        runPM[[X]] <- eval(call(PMlist[X], MSEobj, Ref=ref, Yrs=yrs))  
      }
    }
  }

  pm <- 5
  low <- grepl("<", runPM[[pm]]@Caption)
  
  mean(runPM[[pm]]@Prob[,1] >= runPM[[pm]]@Prob[,2])
  
  plot(runPM[[pm]]@Prob[,1], runPM[[pm]]@Prob[,2], xlim=c(0,1), 
       ylim=c(0,1), pch=16)
  abline(c(0,0), c(1,1), lty=3)
  
  
  runPM[[pm]]@Mean
  
}