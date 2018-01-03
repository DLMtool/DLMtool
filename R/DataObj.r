
# parallelMPs <- function(x, Data, reps, MPs, ss) sapply(ss, MPs[x], Data, reps = reps)

# getTAC <- function(Data, MPs = NA, reps = 100) {
#   if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
#   Data <- updateMSE(Data)
#   nsims <- length(Data@Mort)
#   nMPs <- length(MPs)
#   TACa <- array(NA, dim = c(nMPs, reps, nsims))
#   
#   if (!sfIsRunning()) {
#     for (ff in 1:nMPs) {
#       temp <- sapply(1:nsims, MPs[ff], Data = Data, reps = reps)
#       tt <- sapply(1:1, AvC,  Data=Data, reps=reps)
#       
#       tt@TAC
#       
#       
#       if (mode(temp) == "numeric") 
#         TACa[ff, , ] <- temp
#       if (mode(temp) == "list") {
#         TACa[ff, , ] <- unlist(temp[1, ])
#         for (x in 1:nsims) Data@Misc[[x]] <- temp[2, x][[1]]
#       }
#     }
#   } else {
#     sfExport("Data") 
#     if (nsims < 8) {
#       sfExport(list = c("MPs", "reps"))
#       for (ss in 1:nsims) {
#         temp <- (sfSapply(1:length(MPs), parallelMPs, Data = Data, 
#           reps = reps, MPs = MPs, ss = ss))
#         if (mode(temp) == "list") {
#           Lens <- unlist(lapply(temp, length))
#           for (X in 1:length(Lens)) {
#           Classes <- unlist(lapply(temp[, X][[1]], class))
#           if (length(unique(Classes)) == 1) {
#             # no Misc object
#             TACa[X, , ss] <- unlist(temp[, X])
#           } else {
#             # a Misc object is include
#             ind <- which(Classes == "list")
#             TACa[X, , ss] <- unlist(temp[, X][[1]][1:(ind - 1), 
#             ])
#             Data@Misc[[ss]] <- temp[, X][[1]][ind, ]
#           }
#           }
#         } else {
#           temp <- matrix(temp, nrow = nMPs, ncol = reps, byrow = TRUE)
#           TACa[, , ss] <- temp
#         }
#       }
#     } else {
#       for (ff in 1:nMPs) {
#         temp <- sfSapply(1:nsims, MPs[ff], Data = Data, 
#           reps = reps)
#         if (mode(temp) == "numeric") 
#           TACa[ff, , ] <- temp
#         if (mode(temp) == "list") {
#           TACa[ff, , ] <- unlist(temp[1, ])
#           for (x in 1:nsims) Data@Misc[[x]] <- temp[2, x][[1]]
#         }
#       }
#     }
#   }
#   for (ff in 1:nMPs) {
#     if (sum(is.na(TACa[ff, , ])) > sum(!is.na(TACa[ff, , ]))) {
#       # only plot if there are sufficient non-NA TAC samples
#       print(paste("Method ", MPs[ff], " produced greater than 50% NA values", 
#         sep = ""))
#     }
#   }
#   out <- list(TACa, Data)
#   return(out)
# }
# 
# 

# #' Conduct stock assessment
# #' 
# #' A wrapper function that gets the OFL recommendation in cases where a method
# #' of DLM quota has been specified
# #' 
# #' 
# #' @usage Sam(Data, MPs = NA, reps = 100, perc = 0.5)
# #' @param Data A data-limited methods data object
# #' @param MPs A character vector of methods of DLM quota, DLM space or DLM size
# #' @param reps The number of samples of quota recommendations by method
# #' @param perc quantile of the recommendation to use
# #' @author T. Carruthers
# #' @export Sam
# Sam <- function(Data, MPs = NA, reps = 100, perc = 0.5) {
#   if (class(Data) != "Data") stop("First argument must be object of class 'Data'", call.=FALSE)
#   Data <- updateMSE(Data)
#   nm <- deparse(substitute(DLM))
#   Data@PosMPs <- MPs
#   funcs <- Data@PosMPs
#   nMPs <- length(funcs)
#   Data@MPs <- funcs
#   temp <- getTAC(Data, MPs = funcs, reps)
#   TACa <- temp[[1]]
#   Data <- temp[[2]]
#   nsim <- length(Data@Mort)
#   ref <- array(rep(Data@Ref, nMPs), c(nsim, nMPs))
#   TACm <- apply(TACa, c(3, 1), quantile, p = perc, na.rm = T)
#   TACbias <- (TACm - ref)/ref * 100
#   POF <- round(apply(TACbias > 0, 2, sum)/length(Data@Mort) * 100, 
#     1)
#   Data@TAC <- TACa
#   Data@TACbias <- TACbias
#   Data
# }
# 






