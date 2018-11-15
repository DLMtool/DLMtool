

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
Data2csv <- function(Data, file=NULL, simno = 1, overwrite=F, keepNAs=T) {
  
  if(class(Data)!='Data') stop("First argument 'Data' not an object of class 'Data' ")
  
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
  if(!dir.exists(folder)){
    file=paste0(getwd(),"/",deparse(substitute(Data)),".csv")
    message(paste0("Folder < ",folder," > does not exist -  data to be written to ",file))
  }   
  
  if(file.exists(file)){
    if(!overwrite){
      ANSWER <- readline("This file already exists. Is it ok to continue and overwrite it (y/n) ?")
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
    if(class(obj)=="character"){
      return(FALSE)
    }else if(sum(obj,na.rm=T)==0){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  for(i in 1:ns){
    obj<-slot(Data,slots[i])
    if(class(obj)!="list" & class(obj)!="data.frame"){
      allNA<-sum(!is.na(obj))==0
      if(!(allNA & !keepNAs)){ # write NA values
        if(slottest(obj)){ # is the slot empty or all NAs?
          write(paste(c(slots[i],"NA"),collapse=","),file,1,append=appendy[i])
        }else{
          if(is.null(dim(obj))){  # vector or single value
            write(paste(c(slots[i],obj),collapse=","),file,1,append=appendy[i])
          
          }else if(length(dim(obj))==2){ # a matrix (time series)
            write(paste(c(slots[i],obj[simno,]),collapse=","),file,1,append=appendy[i])
           
          }else if(length(dim(obj))==3){ # 3d array of composition data (CAL, CAA)
             ny<-dim(obj)[2]
            for(yy in 1:ny){
              write(paste(c(paste(slots[i],yrs[yy]),obj[simno,yy,]),collapse=","),file,1,append=appendy[i])
            }
          } 
        } # end of slottest
      } # end of removeNAs
    }
  }
  
}




  