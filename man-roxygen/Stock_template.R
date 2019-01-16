#' @slot Common_Name Common name of the species. Character string
#' @slot Species Scientific name of the species. Genus and species name. Character string
#' @slot maxage The maximum age of individuals that is simulated (there is no 'plus group'). Single value. Positive integer
#' @slot R0 The magnitude of unfished recruitment. Single value. Positive real number
#' @slot M Natural mortality rate. Uniform distribution lower and upper bounds. Positive real number 
#' @slot M2 (Optional) Natural mortality rate at age. Vector of length 'maxage'. Positive real number
#' @slot Mexp Exponent of the Lorenzen function assuming an inverse relationship between M and weight. Uniform distribution lower and upper bounds. Real numbers <= 0.
#' @slot Msd Inter-annual variability in natural mortality rate expressed as a coefficient of variation. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Mgrad No longer used. Previously mean temporal trend in natural mortality rate, expressed as a percentage change in M per year. 
#' @slot h Steepness of the stock recruit relationship. Uniform distribution lower and upper bounds. Values from 1/5 to 1 
#' @slot SRrel Type of stock-recruit relationship. Single value, switch (1) Beverton-Holt (2) Ricker. Integer 
#' @slot Perr Process error, the CV of lognormal recruitment deviations. Uniform distribution lower and upper bounds. Non-negative real numbers
#' @slot AC Autocorrelation in recruitment deviations rec(t)=AC*rec(t-1)+(1-AC)*sigma(t). Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Period (Optional) Period for cyclical recruitment pattern in years. Uniform distribution lower and upper bounds. Non-negative real numbers  
#' @slot Amplitude (Optional) Amplitude in deviation from long-term average recruitment during recruitment cycle (eg a range from 0 to 1 means recruitment decreases or increases by up to 100\% each cycle). Uniform distribution lower and upper bounds. 0 < Amplitude < 1 
#' @slot Linf Maximum length. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot K von Bertalanffy growth parameter k. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot t0 von Bertalanffy theoretical age at length zero. Uniform distribution lower and upper bounds. Non-positive real numbers
#' @slot LenCV Coefficient of variation of length-at-age (assumed constant for all age classes). Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot Ksd Inter-annual variability in growth parameter k expressed as a coefficient of variation. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Kgrad No longer used. Previously mean temporal trend in growth parameter k, expressed as a percentage change in k per year. 
#' @slot Linfsd Inter-annual variability in maximum length expressed as a coefficient of variation. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Linfgrad No longer used. Previously mean temporal trend in maximum length, expressed as a percentage change in Linf per year. 
#' @slot L50 Length at 50 percent maturity. Uniform distribution lower and upper bounds. Positive real numbers 
#' @slot L50_95 Length increment from 50 percent to 95 percent maturity. Uniform distribution lower and upper bounds. Positive real numbers 
# @slot FecB Exponent of the length-fecundity relationship, ie, (relative) fecundity-at-length is proportional to length^FecB (uniform distribution)
#' @slot D Current level of stock depletion SSB(current)/SSB(unfished). Uniform distribution lower and upper bounds. Fraction
#' @slot a Length-weight parameter alpha. Single value. Positive real number 
#' @slot b Length-weight parameter beta. Single value. Positive real number
#' @slot Size_area_1 The size of area 1 relative to area 2. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot Frac_area_1 The fraction of the unfished biomass in stock 1. Uniform distribution lower and upper bounds. Positive real numbers
#' @slot Prob_staying The probability of inviduals in area 1 remaining in area 1 over the course of one year. Uniform distribution lower and upper bounds. Positive fraction.
#' @slot Fdisc Fraction of discarded fish that die. Uniform distribution lower and upper bounds. Non-negative real numbers 
#' @slot Source A reference to a website or article from which parameters were taken to define the stock object. Single value. Character string. 
