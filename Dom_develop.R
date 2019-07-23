
MSEobj <- runMSE(MPs=NA)

dom <- DOM(MSEobj)
library(dplyr)

PMlist <- c("PNOF", "LTY", "P50")
Refs=NULL; Yrs=NULL

Dom(MSEobj, PMlist=c("PNOF", "LTY", "P10", "AAVY"))

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

  # Create Table
  DF <- Required(MSEobj@MPs, TRUE) %>% data.frame(., stringsAsFactors = FALSE)
  DF <- dplyr::left_join(DF, MPtype(MSEobj@MPs), by="MP")

  DF3 <- lapply(runPM, slot, name='Mean') %>% do.call("cbind", .)
  colnames(DF3) <- PMlist
  DF3 <- data.frame(MP=runPM[[1]]@MPs, DF3, stringsAsFactors = FALSE)
  DF <- dplyr::left_join(DF, DF3, by="MP")

  DomMPs <- Top <- list()
  for (i in 1:MSEobj@nMPs) {
    # ind <- which(DF$DataClass == DF$DataClass[i] & DF$Rec == DF$Rec[i] & DF$Type ==DF$Type[i])
    ind <- which(DF$Rec == DF$Rec[i] & DF$Type ==DF$Type[i])
    ind <- ind[!ind==i]
    if (length(ind)>0) {
      df1 <- DF[i,] %>% dplyr::select(-MP, -Data, -DataClass, -Type, -Recs) %>% unlist()
      m1 <- matrix(df1, nrow=length(ind), ncol=length(df1), byrow=TRUE) %>% round(2)
      df2 <- DF[ind,] %>% dplyr::select(-MP, -Data, -DataClass, -Type, -Recs) %>% as.matrix() %>% round(2)
      ind2 <- which(rowSums(m1 < df2) ==0) 
      if (length(ind2)>0) {
        # Top[[i]] <-  data.frame(MP=DF$MP[i],Rec=DF$Recs[i], DataClass=DF$DataClass[i])
        DomMPs <- append(DomMPs, DF$MP[ind[ind2]])
      } else {
    
      }
    } else {
   
    }
  }
  DomMPs <- unlist(DomMPs) %>% unique()
  NonDom <- MSEobj@MPs[!MSEobj@MPs %in% DomMPs]

  TradePlot(Sub(MSEobj, MPs=NonDom), PMlist=PMlist)
  

}
