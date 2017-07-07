# Input Control Functions Wrapper function for input control methods


#' Runs input control MPs on a Data object.
#' 
#' Function runs a MP (or MPs) of class 'Input' and returns a list: input
#' control recommendation(s) in element 1 and Data object in element 2.
#' 
#' 
#' @usage runInMP(Data, MPs = NA, reps = 100)
#' @param Data A object of class Data
#' @param MPs A vector of MPs of class 'Input'
#' @param reps Number of stochastic repititions - often not used in input
#' control MPs.
#' @author A. Hordyk
#' @export runInMP
runInMP <- function(Data, MPs = NA, reps = 100) {
  
  nsims <- length(Data@Mort)
  nareas <- Data@nareas 
  nMPs <- length(MPs)
  
  returnList <- list() # a list nMPs long containing MPs recommendations
  recList <- list() # a list containing nsim recommendations from a single MP 
  
  if (!sfIsRunning() | (nMPs < 8 & nsims < 8)) {
    for (ff in 1:nMPs) {
      temp <- sapply(1:nsims, MPs[ff], Data = Data, reps = reps)
      slots <- slotNames(temp[[1]])
      for (X in slots) { # sequence along recommendation slots 
        if (X == "Misc") { # convert to a list nsim by nareas
          rec <- lapply(temp, slot, name=X)
        } else {
          rec <- unlist(lapply(temp, slot, name=X))
        }
        if (X == "Spatial") { # convert to a matrix nsim by nareas
          rec <- matrix(rec, nsims, nareas, byrow=TRUE)  
        }
    
        recList[[X]] <- rec
        for (x in 1:nsims) Data@Misc[[x]] <- recList$Misc[[x]]
        recList$Misc <- NULL
      }
      returnList[[ff]] <- recList
    }
  } else {
    sfExport(list = c("Data"))
    for (ff in 1:nMPs) {
      temp <- sfSapply(1:nsims, MPs[ff], Data = Data, reps = reps)
      slots <- slotNames(temp[[1]])
      for (X in slots) { # sequence along recommendation slots 
        if (X == "Misc") { # convert to a list nsim by nareas
          rec <- lapply(temp, slot, name=X)
        } else {
          rec <- unlist(lapply(temp, slot, name=X))
        }
        if (X == "Spatial") { # convert to a matrix nsim by nareas
          rec <- matrix(rec, nsims, nareas, byrow=TRUE)  
        }
        
        recList[[X]] <- rec
        for (x in 1:nsims) Data@Misc[[x]] <- recList$Misc[[x]]
        recList$Misc <- NULL
      }
      returnList[[ff]] <- recList
    }   
  }
  
  return(list(returnList, Data))
}