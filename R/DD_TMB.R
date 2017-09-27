## ADD AUTHORSHIP DETAILS 

## IN DEVELOPMENT


#' Delay - Difference Stock Assessment with UMSY and MSY leading using TMB
#' 
#' A simple delay-difference assessment that estimates the TAC using a
#' time-series of catches and a relative abundance index and coded with TMB
#' 
#' 
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @return A numeric vector of TAC recommendations
#' @note This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Similar to many other assessment
#' models it depends on a whole host of dubious assumptions such as temporally
#' stationary productivity and proportionality between the abundance index and
#' real abundance. Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers & Z. Siders. Zach Siders coded the TMB function 
#' @references Method based on equations of Carl Walters (bug him with
#' questions and expect colourful responses)
#' @export
#' @importFrom TMB MakeADFun sdreport
#' @useDynLib DD_tmb_sim
DD_TMB <- function(x, Data, reps = 100) {
  dependencies = "Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@vbt0, Data@CV_vbt0, Data@Mort, Data@CV_Mort, Data@wla, Data@wlb"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  a <- Data@wla[x]
  b <- Data@wlb[x]
  
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])
  a50V <- max(a50V, 1)
  yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x,   ])]
  C_hist <- Data@Cat[x, yind]
  E_hist <- C_hist/Data@Ind[x, yind]
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYprior <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3)
  
  data <- list(So_DD = So_DD, Alpha_DD = Alpha_DD, Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, C_hist = C_hist, UMSYprior = UMSYprior)
  params <- list(log_UMSY_DD = log(Data@Mort[x]), log_MSY_DD = log(mean(C_hist, na.rm = T)), log_q_DD = log(Data@Mort[x]))
  info <- list(data = data, params = params)
  
  # Fit model 
  Obj <- TMB::MakeADFun(data=info$data, parameters=info$params, DLL="DD_tmb_sim", silent=TRUE) #makes TMB function
  Opt <- stats::nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr) #optimizes
  SD <- TMB::sdreport(Obj) #gets parameter estimates and std. error using ADREPORT
  
  SDval <- SD$value
  SDsd <- sapply(SDval,function(x) (x^2)^0.5*0.1)
  TAC <- rep(NA, reps)
  samps <- cbind(rnorm(reps, SDval[1], SDsd[1]), rnorm(reps, SDval[2], SDsd[2]), rnorm(reps, SDval[3], SDsd[3]))
  if (reps == 1)  samps <- matrix(c(SDval[1], SDval[2], SDval[3]), nrow = 1)
  for (i in 1:reps) 
    TAC[i] <- DD_R(samps[i, ], opty = 2, So_DD = data$So_DD, Alpha_DD = data$Alpha_DD, 
                   Rho_DD = data$Rho_DD, ny_DD = data$ny_DD, k_DD = data$k_DD, 
                   wa_DD = data$wa_DD, E_hist = data$E_hist, C_hist = data$C_hist, 
                   UMSYprior = data$UMSYprior)
  TACfilter(TAC)
  
}
class(DD_TMB) <- "Output"


