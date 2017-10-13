
#' Valid custom parameters (cpars)
#'
#' @param print Print the valid names for cpars?
#'
#' @return invisibly returns vector of valid cpars names
#' @export
#'
validcpars <- function(print=TRUE) {
  vnames <- sort(c("dep","Esd","Find","procsd","AC","M","Msd", 
                   "Mgrad","hs","Linf","Linfsd","Linfgrad",
                   "K","Ksd","Kgrad","t0","L50", "L95", "L50_95","Spat_targ",
                   "Frac_area_1","Prob_staying","Size_area_1", 
                   "Csd","Cbias","CAA_nsamp","CAA_ESS","CAL_nsamp",
                   "CAL_ESS","betas","Isd","Derr","Dbias", 
                   "Mbias","FMSY_Mbias","lenMbias","LFCbias",
                   "LFSbias","Aerr","Abias","Kbias","t0bias", 
                   "Linfbias","Irefbias","Crefbias","Brefbias",
                   "Recsd","qinc","qcv","L5","LFS","Vmaxlen","Perr","R0","Mat_age", 
                   "Mrand","Linfrand","Krand","maxage","V",  
                   "ageM", "age95", "EffYears", "EffLower", "EffUpper",
                   "Wt_age", "Len_age", "Marray", "M_at_Length", "LenCV", 
                   "CAL_binsmid", "CAL_bins", "LatASD", "dFfinal",
                   "LR5", "LFR", "Rmaxlen", "DR", "Fdisc")) 
  
  if (print) {
    n <- length(vnames)
    vec <- 3:7
    nc <- vec[which.min(n  %% vec)]
    
    options(warn=-1)
    temp <- matrix(vnames, ncol=nc, byrow=TRUE)
    options(warn=1)
    temp[duplicated(as.vector(temp))] <- ""
    print(temp)
  }
  invisible(vnames)
}





#' Sample custom pars
#'
#' @param cpars A named list containing custom parameters for the OM
#' @param nsim number of simulations
#' @param msg logical - print the names of the cpars? Turn off when using the function in a loop
#' @return A named list of sampled custom parameters
#' @export
#'
SampleCpars <- function(cpars, nsim=48, msg=TRUE) {
 
  # Vector of valid names for custompars list or data.frame. Names not in this list will be printed out in warning and ignored #	
  ParsNames <- validcpars(FALSE)
  
  sampCpars <- list()
  ncparsim<-cparscheck(cpars)
  Names <- names(cpars)
  # report invalid names 
  invalid <- which(!Names %in% ParsNames)
  if (length(invalid) > 0) {
    outNames <- paste(Names[invalid], "")
    for (i in seq(5, by=5, length.out=floor(length(outNames)/5))) outNames <- gsub(outNames[i], paste0(outNames[i], "\n"), outNames)
    if(msg) message("ignoring invalid names found in custom parameters (OM@cpars) \n", outNames)	
  }
  # report found names
  valid <- which(Names %in% ParsNames)
  cpars <- cpars[valid]
  if (length(valid) == 0) stop("No valid names found in custompars (OM@cpars)", call.=FALSE)
  Names <- names(cpars)
  outNames <- paste(Names, "")
  for (i in seq(5, by=5, length.out=floor(length(outNames)/5)))
    outNames <- gsub(outNames[i], paste0(outNames[i], "\n"), outNames)
  if(msg) message("valid custom parameters (OM@cpars) found: \n", outNames)
  
  # Sample custom pars 
  if (ncparsim < nsim) ind <- sample(1:ncparsim, nsim, replace=TRUE)
  if (!ncparsim < nsim) ind <- sample(1:ncparsim, nsim, replace=FALSE)

  for (i in 1:length(cpars)) {
    samps <- cpars[[i]]
    name <- names(cpars)[i]
    if (any(c("EffUpper", "EffLower", "EffYears", "maxage", "M_at_Length", "CAL_binsmid", "CAL_bins") %in% name)) {
      sampCpars[[name]] <- samps
    } else {
      if (class(samps) == "numeric" | class(samps) == "integer") sampCpars[[name]] <- samps[ind]
      
      if (class(samps) == "matrix") sampCpars[[name]] <- samps[ind,, drop=FALSE] 
      
      if (class(samps) == "array") {
        if (length(dim(samps)) == 3)  sampCpars[[name]] <- samps[ind, , ,drop=FALSE]
      }
      if (class(samps) == "data.frame")   sampCpars[[name]] <- samps 
    }
  }
  

  
  sampCpars
}