#' @slot nyears The number of years for the historical 'spool-up' simulation. Single value. Positive integer 
#' @slot Spat_targ Distribution of fishing in relation to spatial biomass: fishing distribution is proportional to B^Spat_targ. Uniform distribution lower and upper bounds. Real numbers   
#' @slot EffYears Years representing join-points (vertices) of time-varying effort. Vector. Non-negative real numbers 
#' @slot EffLower Lower bound on relative effort corresponding to EffYears. Vector. Non-negative real numbers
#' @slot EffUpper Upper bound on relative effort corresponding to EffYears. Vector. Non-negative real numbers 
#' @slot Esd Additional inter-annual variability in fishing mortality rate. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot qinc Average percentage change in fishing efficiency (applicable only to forward projection and input controls). Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot qcv Inter-annual variability in fishing efficiency (applicable only to forward projection and input controls). Uniform distribution lower and upper bounds. Non-negative real numbers

#' @slot L5 Shortest length corresponding to 5 percent vulnerability. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LFS Shortest length that is fully vulnerable to fishing. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot Vmaxlen The vulnerability of fish at \code{Stock@Linf}. Uniform distribution lower and upper bounds. Fraction 
#' @slot isRel Selectivity parameters in units of size-of-maturity (or absolute eg cm). Single value. Boolean.
#' @slot LR5 Shortest length corresponding ot 5 percent retention. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot LFR Shortest length that is fully retained. Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot Rmaxlen The retention of fish at \code{Stock@Linf}. Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot DR Discard rate - the fraction of caught fish that are discarded. Uniform distribution lower and upper bounds. Fraction
#' 
#' @slot SelYears (Optional) Years representing join-points (vertices) at which historical selectivity pattern changes. Vector. Positive real numbers
#' @slot AbsSelYears (Optional) Calendar years corresponding with SelYears (eg 1951, rather than 1), used for plotting only. Vector (of same length as SelYears). Positive real numbers  
#' @slot L5Lower (Optional) Lower bound of L5 (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot L5Upper (Optional) Upper bound of L5 (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot LFSLower (Optional) Lower bound of LFS (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot LFSUpper (Optional) Upper bound of LFS (use \code{ChooseSelect} function to set these). Vector. Non-negative real numbers 
#' @slot VmaxLower (Optional) Lower bound of Vmaxlen (use \code{ChooseSelect} function to set these). Vector. Fraction 
#' @slot VmaxUpper (Optional) Upper bound of Vmaxlen (use \code{ChooseSelect} function to set these). Vector. Fraction
#' @slot CurrentYr The current calendar year (final year) of the historical simulations (eg 2011). Single value. Positive integer. 
#' 
#' @slot MPA (Optional) Matrix specifying spatial closures for historical years. 
#' 
#' @section MPA slot: 
#' Each row should contain year index (e.g 10 for 10th historical year)
#' followed by fraction of area closed to fishing for each area. i.e. each row represents a change and the number of columns is nareas + 1. 
#' The spatial closures are assumed to remain in place for the future projections unless changed by a MP. 
#' Default (if left blank) is all areas are open to fishing in historical period.
