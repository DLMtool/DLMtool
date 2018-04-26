#' MP feasibility diagnostic
#' 
#' What MPs may be run (best case scenario) for various data-availability
#' scenarios?
#' 
#' 
#' @param feaseobj An object of class 'Fease'
#' @param outy Determines whether you would like a full table or some column of
#' the table for a specific case of the feasibility object. When set equal to
#' table, the full table is produced. When set equal to an integer number the
#' names of MPs that are feasible for that case are returned.
#' @author T. Carruthers
#' @export Fease



# Fease <- function(feaseobj, outy = "table") {
#   
#   if (class(feaseobj) != "Fease") 
#     stop("Incorrect format: you need an object of class Fease")
#   
#   sloty <- c("Cat", "Ind", "AvC", "Dt", "Rec", "CAA", "CAL", "Mort", 
#              "L50", "L95", "vbK", "vbLinf", "vbt0", "wla", "wlb", "steep", "LFC", 
#              "LFS", "Cref", "Bref", "Iref", "Dep", "Abun", "ML")
#   
#   type <- c("Catch", "Index", "Catch", "Index", "Recruitment_index", 
#             "Catch_at_age", "Catch_at_length", "Natural_mortality_rate", "Maturity_at_length", 
#             "Maturity_at_length", "Growth", "Growth", "Growth", "Length_weight_conversion", 
#             "Length_weight_conversion", "Stock_recruitment_relationship", "Fleet_selectivity", 
#             "Fleet_selectivity", "Target_catch", "Target_biomass", "Target_index", 
#             "Index", "Abundance")
#   
#   ncases <- length(feaseobj@Case)
#   slots <- slotNames(feaseobj)
#   ns <- length(slots)
#   ftab <- array(TRUE, c(ns - 2, ncases))
#   for (j in 3:ns) ftab[j - 2, ] <- as.logical(as.numeric(slot(feaseobj, 
#                                                               slots[j])))
#   
#   req <- Required()
#   nMPs <- nrow(req)
#   gridy <- array("", c(nMPs, ncases))
#   for (i in 1:ncases) {
#     types <- slotNames(feaseobj)[3:17][ftab[, i]]
#     slots <- sloty[type %in% types]
#     for (m in 1:nMPs) {
#       brec <- unlist(strsplit(req[m, 2], ", "))
#       brec <- brec[grep("CV_", brec, invert = T)]  #remove CV dependencies (we think we can guess these...)
#       brec <- brec[brec != "Year" & brec != "MaxAge" & brec != "FMSY_M" & 
#                      brec != "BMSY_B0" & brec != "t" & brec != "OM" & brec != 
#                      "MPrec" & brec != "CAL_bins" & brec != "MPeff" & brec != 
#                      "LHYear"]
#       nr <- length(brec)
#       if (nr == 0) {
#         gridy[m, i] <- "Yes"
#       } else {
#         cc <- 0
#         for (r in 1:nr) {
#           # loop over requirements
#           if (brec[r] %in% slots) 
#             cc <- cc + 1
#         }
#         if (cc == nr) 
#           gridy[m, i] <- "Yes"
#       }
#     }
#   }
#   gridy <- as.data.frame(gridy)
#   row.names(gridy) = req[, 1]
#   names(gridy) = feaseobj@Case
#   if (outy == "table") 
#     return(gridy)
#   if (outy != "table" & class(outy) != "numeric") 
#     return(req[, 1][gridy[, 1] == "Yes"])
#   if (class(outy) == "numeric") {
#     if (outy < (ncases + 1)) {
#       return(req[, 1][gridy[, as.integer(outy)] == "Yes"])
#     } else {
#       return(req[, 1][gridy[, 1] == "Yes"])
#     }
#   }
#   
# }