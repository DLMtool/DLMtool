#' @slot Cobs Log-normal catch observation error expressed as a coefficient of variation. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Cbiascv Log-normal coefficient of variation controlling the sampling of bias in catch observations for each simulation. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot CAA_nsamp Number of catch-at-age observation per time step. Uniform distribution lower and upper bounds. Positive real numbers   
#' @slot CAA_ESS Effective sample size (independent age draws) of the multinomial catch-at-age observation error model. Uniform distribution lower and upper bounds. Positive integers
#' @slot CAL_nsamp Number of catch-at-length observation per time step. Uniform distribution lower and upper bounds. Positive integers
#' @slot CAL_ESS Effective sample size (independent length draws) of the multinomial catch-at-length observation error model. Uniform distribution lower and upper bounds. Positive integers
# #' @slot CALcv Log-normal, CV of length-at-age. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Iobs Observation error in the relative abundance indices expressed as a coefficient of variation. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Ibiascv Not Used. Log-normal coefficient of variation controlling error in observations of relative abundance index. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Btobs Log-normal coefficient of variation controlling error in observations of current stock biomass among years. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Btbiascv Uniform-log bounds for sampling persistent bias in current stock biomass. Uniform-log distribution lower and upper bounds. Positive real numbers 
#' @slot beta A parameter controlling hyperstability/hyperdepletion where values below 1 lead to hyperstability (an index that decreases slower than true abundance) and values above 1 lead to hyperdepletion (an index that decreases more rapidly than true abundance). Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LenMbiascv Log-normal coefficient of variation for sampling persistent bias in length at 50 percent maturity. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Mbiascv Log-normal coefficient of variation for sampling persistent bias in observed natural mortality rate. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Kbiascv Log-normal coefficient of variation for sampling persistent bias in observed growth parameter K. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot t0biascv Log-normal coefficient of variation for sampling persistent bias in observed t0. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Linfbiascv Log-normal coefficient of variation for sampling persistent bias in observed maximum length. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LFCbiascv Log-normal coefficient of variation for sampling persistent bias in observed length at first capture. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot LFSbiascv Log-normal coefficient of variation for sampling persistent bias in length-at-full selection. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot FMSYbiascv Not used. Log-normal coefficient of variation for sampling persistent bias in FMSY. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot FMSY_Mbiascv Log-normal coefficient of variation for sampling persistent bias in FMSY/M. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot BMSY_B0biascv Log-normal coefficient of variation for sampling persistent bias in BMSY relative to unfished. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Irefbiascv Log-normal coefficient of variation for sampling persistent bias in relative abundance index at BMSY. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Brefbiascv Log-normal coefficient of variation for sampling persistent bias in BMSY. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Crefbiascv Log-normal coefficient of variation for sampling persistent bias in MSY. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Dbiascv Log-normal coefficient of variation for sampling persistent bias in stock depletion. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Dobs Log-normal coefficient of variation controlling error in observations of stock depletion among years. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot hbiascv Log-normal coefficient of variation for sampling persistent bias in steepness. Uniform distribution lower and upper bounds. Positive real numbers  
#' @slot Recbiascv Log-normal coefficient of variation for sampling persistent bias in recent recruitment strength. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot B0cv Log-normal coefficient of variation for sampling persistent bias in unfished biomass. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot rcv Log-normal coefficient of variation for sampling persistent bias in intrinsic rate of increase. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot Fcurbiascv Log-normal coefficient of variation for sampling persistent bias in current fishing mortality rate. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot Fcurcv Log-normal coefficient of variation controlling error in observations of current fishing mortality rate among years. Uniform distribution lower and upper bounds. Positive real numbers 
# #' @slot maxagecv Log-normal coefficient of variation for sampling persistent bias in observation of maximum age. Uniform distribution lower and upper bounds. Positive real numbers  
