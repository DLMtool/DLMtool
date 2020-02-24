

#' Converts a Data object into a .csv data file
#'
#' @description A function that writes a correctly formatted .csv file from a DLMtool / MSEtool Data object
#' @param Data An object of class 'Data'. 
#' @param file Character string. The name of the location and file you wish to create (e.g. "C:/temp/mydata.csv")
#' @param simno Integer. An optional argument to specify the simulation number if writing simulated data
#' @param overwrite Boolean. Should existing data files be automatically overwritten. 
#' @param keepNAs Boolean. Should slots with NAs still be written to the data file. 
#' @author T. Carruthers
#' @export
Data2csv<-function(Data, file=NULL, simno = 1, overwrite=F, keepNAs=T) {
  
  if(!inherits(Data,'Data')) stop("First argument 'Data' not an object of class 'Data' ")
  
  if(is.null(file)){
    file=paste0(getwd(),"/",deparse(substitute(Data)),".csv")
    message(paste0("Second argument 'file' not specified - data to be written to ",file))
  }
  
  if(substr(file,nchar(file),nchar(file))=="/"){
    file=paste0(file,deparse(substitute(Data)),".csv")
    message(paste0("Second argument 'file' appears to be only a folder, no file name specified - data to be written to ",file))
  }  
  
  if(substr(file,nchar(file)-3,nchar(file))!=".csv") file=paste0(file,".csv")
  
  dirsplit<-unlist(strsplit(file,split="/"))
  folder<-paste(dirsplit[1:(length(dirsplit)-1)],collapse="/") 
  if (folder != file) {
    if(!dir.exists(folder)){
      file=paste0(getwd(),"/",deparse(substitute(Data)),".csv")
      message(paste0("Folder < ",folder," > does not exist -  data to be written to ",file))
    }  
  } else {
    message("Data to be written to ", getwd(), "/", file)
  }

  if(file.exists(file)){
    if(!overwrite){
      ANSWER<- readline("This file already exists. Is it ok to continue and overwrite it (y/n) ?")
      ## a better version would check the answer less cursorily, and
      ## perhaps re-prompt
      if (substr(ANSWER, 1, 1) == "n"){
        stop("Data file writing aborted")
      }
    }
    message(paste0("File < ",file," > is being over-written"))
  }
  
  slots<-slotNames(Data)
  ns<-length(slots)
  
  appendy<-c(F,rep(T,ns-1))
  
  yrs<-Data@Year
  
  slottest<-function(obj){
    #if(sum(!is.na(obj))==0){
    #  return(TRUE)
    if(inherits(obj,"character")){
      return(FALSE)
    }else if(sum(obj,na.rm=T)==0){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  lexicon<-t(matrix(
    c("Name","Name",
      "Common_Name", "Common Name",
      "Species", "Species",
      "Region",  "Region",
      "LHYear","Last Historical Year",
      "MPrec","Previous TAC",
      "Units", "Units",
      "MPeff","Previous TAE",
      "nareas","nareas",
      
      "MaxAge","Maximum age",
      "Mort","M",
      "CV_Mort","CV M",
      "vbLinf","Von Bertalanffy Linf parameter",
      "CV_vbLinf","CV von B. Linf parameter",
      "vbK","Von Bertalanffy K parameter",
      "CV_vbK","CV von B. K parameter",
      "vbt0","Von Bertalanffy t0 parameter",
      "CV_vbt0","CV von B. t0 parameter",
      "wla","Length-weight parameter a",
      "CV_wla","CV Length-weight parameter a",
      "wlb","Length-weight parameter b",
      "CV_wlb","CV Length-weight parameter b",
      
      "steep","Steepness",
      "CV_steep","CV Steepness",
      "sigmaR","sigmaR",
      "CV_sigmaR", "CV sigmaR",
      "L50","Length at 50% maturity",
      "CV_L50","CV Length at 50% maturity",
      "L95","Length at 95% maturity",
      "LenCV", "CV of length-at-age",
      
      "LFC","Length at first capture",
      "CV_LFC","CV Length at first capture",
      "LFS","Length at full selection",
      "CV_LFS","CV Length at full selection",
      "Vmaxlen", "Vulnerability at asymptotic length",
      
      "Year","Year", 
      "Cat", "Catch",
      "CV_Cat","CV Catch",
      "Effort", "Effort",
      "CV_Effort", "CV Effort",
      "Ind", "Abundance index",
      "CV_Ind","CV Abundance index",
      
      'SpInd', 'Spawning Abundance index',
      'CV_SpInd', 'CV Spawning Abundance index',
      
      'VInd', 'Vulnerable Abundance index',
      'CV_VInd', 'CV Vulnerable Abundance index',
      
      "Rec", "Recruitment index",
      "CV_Rec","CV Recruitment index",
      
      "ML", "Mean length",
      "Lc", "Modal length (Lc)",
      "Lbar", "Mean length above Lc",
  
      "Vuln_CAA", "Vuln CAA",
      "CAA", "CAA",
      
      'Vuln_CAL', "Vuln CAL",
      "CAL_bins","CAL_bins",
      "CAL_mids","CAL_mids",
      "CAL","CAL",
      
      "Dep","Current stock depletion",
      "CV_Dep","CV current stock depletion",
      "Abun","Current stock abundance",
      "CV_Abun","CV current stock abundance",
      "SpAbun","Current spawning abundance",
      "CV_SpAbun","CV current spawning abundance",
      
      "FMSY_M","FMSY/M",
      "CV_FMSY_M","CV FMSY/M",
      "BMSY_B0","BMSY/B0",
      "CV_BMSY_B0","CV BMSY/B0",
      
      "Cref","Catch Reference",
      "CV_Cref","CV Catch Reference",
      "Bref","Biomass Reference",
      "CV_Bref","CV Biomass Reference",
      "Iref","Index Reference",
      "CV_Iref","CV Index Reference",
      "t","Duration t",
      "AvC","Average catch over time t",
      "CV_AvC","CV Average catch over time t",
      "Dt","Depletion over time t",
      "CV_Dt","CV Depletion over time t",
      "Ref","Reference OFL",
      "Ref_type", "Reference OFL type"),
      nrow=2))
  
  # TO DO - add AddInd & CV Index 1 etc
  # if(!all(is.na(Data@AddInd))) {
  #   nind <- dim(Data@AddInd)[2]
  #   paste("Index", 1:nind)
  #   
  # }
 
  nsim <- length(Data@Mort)
  for(i in 1:ns){
    obj<-slot(Data,slots[i])
    lex<-lexicon[match(slots[i],lexicon[,1]),2]
    if(!is.na(lex)){
      if(!inherits(obj,"list") & !inherits(obj,"data.frame")){
        allNA<-sum(!is.na(obj))==0
        if(!(allNA & !keepNAs)){ # write NA values
          if(slottest(obj)){ # is the slot empty or all NAs?
            write(paste(c(lex,"NA"),collapse=","),file,1,append=appendy[i])
          }else{
            if(is.null(dim(obj))){  # vector or single value
              if (length(obj)==nsim) {
                write(paste(c(lex,obj[simno]),collapse=","),file,1,append=appendy[i])  
              } else {
                write(paste(c(lex,obj),collapse=","),file,1,append=appendy[i])
              }
            }else if(length(dim(obj))==2){ # a matrix (time series)
              write(paste(c(lex,obj[simno,]),collapse=","),file,1,append=appendy[i])
             
            }else if(length(dim(obj))==3){ # 3d array of composition data (CAL, CAA)
               ny<-dim(obj)[2]
               if (slots[i] == "RInd") {
                 for(yy in 1:ny){
                   if (yy ==1) write(paste(c(lex,obj[simno,yy,]),collapse=","),file,1,append=appendy[i])
                   if (yy !=1) write(paste(c('',obj[simno,yy,]),collapse=","),file,1,append=appendy[i])
                 }
               } else {
                 for(yy in 1:ny){
                   write(paste(c(paste(lex,yrs[yy]),obj[simno,yy,]),collapse=","),file,1,append=appendy[i])
                 }
               }
              
            } 
          } # end of slottest
        } # end of removeNAs
      }
    }#end of na lex  
  }
}


#' Simplifies the CAL slot of data object
#'
#' @description A function that condenses the number of catch-at-length bins in a data object
#' @param Data An object of class 'Data'. 
#' @param nbins Integer. The target number of catch at length bins
#' @param simno Integer. An optional argument to specify the simulation number if writing simulated data
#' @author T. Carruthers
#' @export
CALsimp<-function(Data,nbins=10,simno=1){
  
  oldbins<-Data@CAL_bins
  nold<-length(oldbins)
  ind<-rep((1:nold),each=floor(nold/nbins))[1:nold]
  maxbin<-max(ind)
  newCAL_bins<-c(Data@CAL_bins[match(1:maxbin,ind)],Data@CAL_bins[nold])
  ny<-dim(Data@CAL)[2]
  newCAL<-array(0,c(1,ny,maxbin))
  for(b in 1:(nold-1)) newCAL[1,,ind[b]]<-newCAL[1,,ind[b]]+Data@CAL[simno,,b]
    
  Data@CAL_bins<-newCAL_bins
  Data@CAL<-newCAL
  Data  
}

  