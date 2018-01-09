
# Functions for sampling steepness parameter (h) from the Data object
# 
# Code from Quang Huynh that fixes the bug where h is sometimes sampled > 1 or < 0.2

#' Sample steepness given mean and cv
#'
#' @param n number of samples
#' @param mu mean h  
#' @param cv cv of h
#'
#' @author Q. Huynh
#' @export
#'
sample_steepness2 <- function(n, mu, cv) {
  if(n == 1) h <- mu
  else {
    sigma <- mu * cv
    mu.beta.dist <- (mu - 0.2)/0.8
    sigma.beta.dist <- sigma/0.8
    beta.par <- derive_beta_par(mu.beta.dist, sigma.beta.dist)
    h.transformed <- rbeta(n, beta.par[1], beta.par[2])
    h <- 0.8 * h.transformed + 0.2
    h[h > 0.99] <- 0.99
    h[h < 0.2] <- 0.2
  }
  return(h)
}



#' This function reduces the CV by 5 per cent until steepness values can be sampled without error
#'
#'
#' @param mu mean h 
#' @param sigma sd of h 
#'
#' @author Q. Huynh
#' @export
#'
derive_beta_par <- function(mu, sigma) {
  
  a <- alphaconv(mu, sigma)
  b <- betaconv(mu, sigma)
  
  if(a <= 0 || b <= 0) {
    sigma <- 0.95 * sigma
    Recall(mu, sigma)
  }
  else return(c(a, b))
  
}
