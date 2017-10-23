

#' Double-normal selectivity curve
#'
#' @param lens Vector of lengths 
#' @param lfs Length at full selection
#' @param sl Sigma of ascending limb
#' @param sr Sigma of descending limb
#'
#' @export
#'
dnormal<-function(lens,lfs,sl,sr){
  cond<-lens<=lfs
  sel<-rep(NA,length(lens))
  sel[cond]<-2.0^-((lens[cond]-lfs)/sl*(lens[cond]-lfs)/sl)
  sel[!cond]<-2.0^-((lens[!cond]-lfs)/sr*(lens[!cond]-lfs)/sr)
  sel
}


#' Calculate selectivity curve
#'
#' @param x Simulation number
#' @param lens Matrix of lengths (nsim by nlengths)
#' @param lfs Vector of length at full selection (nsim long)
#' @param sls Vector of sigmas of ascending limb (nsim long)
#' @param srs Vector of sigmas of descending limb (nsim long)
#'
#' @export
#'
getsel <- function(x, lens, lfs, sls, srs) {
  if (is.null(ncol(lens))) return(dnormal(lens, lfs[x], sls[x], srs[x]))
  dnormal(lens[x,], lfs[x], sls[x], srs[x])
}


# 
# 
# #' Internal Two sided selectivity curve
# #'
# #' @param LFS first length at fulll selection 
# #' @param s1 ascending slope
# #' @param s2 descending slope 
# #' @param lens vector of lengths 
# #' @keywords internal
# #' @export TwoSidedFun
# TwoSidedFun <- function(LFS, s1, s2, lens) {
#   Sl <- rep(0, length(lens))
#   Sl[lens < LFS] <- exp(-((lens[lens < LFS] - LFS)^2)/(2 * s1^2))
#   Sl[lens >= LFS] <- exp(-((lens[lens >= LFS] - LFS)^2)/(2 * s2^2))
#   return(Sl)
# }
# 
# #' Internal function to calculate ascending slope of selectivity curve 
# #'
# #' @param s1 slope 1 
# #' @param LFS first length at fulll selection 
# #' @param L0.05 length at 5 percent selection 
# #' @keywords internal
# #' @export getSlope1
# getSlope1 <- function(s1, LFS, L0.05) 
#   (0.05 - TwoSidedFun(LFS, s1 = s1, s2 = 1E5, lens=L0.05))^2
# 
# 
# #' Internal function to calculate slope
# #'
# #' @param s2 desceding slope 
# #' @param LFS length one
# #' @param s1 ascending slope 
# #' @param maxlen length of oldest age class
# #' @param MaxSel selectivity of maxlen 
# #' @keywords internal
# #' @export getSlope2
# getSlope2 <- function(s2, LFS, s1, maxlen, MaxSel) 
#   (MaxSel - TwoSidedFun(LFS, s1, s2, maxlen))^2
# 
# 
# 
# #' Selectivity at length function
# #'
# #' @param i index 
# #' @param SL0.05 length at 5 percent selection
# #' @param SL1 length at full selection
# #' @param MaxSel Maximum selectivity
# #' @param maxlens maximum length 
# #' @param Lens vector of lengths 
# #' @keywords internal
# #' @export SelectFun
# SelectFun <- function(i, SL0.05, SL1, MaxSel, maxlens, Lens) {
#   s1 <- optimise(getSlope1, interval = c(0, 1e+06), LFS = SL1[i], L0.05 = SL0.05[i])$minimum
#   s2 <- optimise(getSlope2, interval = c(0, 1e+06), LFS = SL1[i], s1 = s1, maxlen = maxlens[i], 
#                  MaxSel = MaxSel[i])$minimum
#   if (is.vector(Lens)) TwoSidedFun(LFS = SL1[i], s1 = s1, s2 = s2, lens = Lens)  #nsim = 1
#   else TwoSidedFun(LFS = SL1[i], s1 = s1, s2 = s2, lens = Lens[i, ])  #nsim > 1
# }






