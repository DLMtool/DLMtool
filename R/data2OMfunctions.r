

#' Generate bootstrapped estimates of von Bertalanffy growth parameters from length-at-age data
#' 
#' The von Bertalanffy model is fitted to length-at-age data and bootstrapped to 
#' provide either `OM@nsim` or `nsim` estimates of Linf, K, and t0 parameters. 
#' LenCV is also estimated from mean growth curve and the observed length-at-age data.
#' 
#' The function either returns an OM with the cpars slot updated with estimated values,
#' or a data.frame with the estimated values.
#'
#' @param data A data.frame with columns names 'Age' and 'Length'
#' @param OM Optional. Object of class `OM`. Function will return OM with `cpars`
#' slot populated if an OM is provided, otherwise it returns a data.frame
#' @param nsim Optional numeric. If an OM is not provided, nsim is used.
#' @param seed Optional numierc. If an OM is not provided, seed is used in `set.seed`.
#' @param plot Logical. Plot the data and model fits?
#' @param msg Logical. Display messages?
#'
#' @return An OM with cpars slot populated with Linf, K, t0 and LenCV values, 
#' or a data.frame.
#' @export 
#' 
#' @examples
#' # Simple model to generate length-at-age data
#' OM <- DLMtool::testOM 
#' OM@nsim <- 2
#' Hist <- runMSE(OM, Hist=TRUE)
#' N <- Hist@AtAge$Nage[1,,1] * Hist@AtAge$Select[1,,1]
#' meanL <- Hist@AtAge$Length[1,,1]
#' sdL <- Hist@AtAge$Length[1,,1] * 0.1
#' nsamp <- ceiling(N/sum(N) * 1000)
#' Length <- unlist(sapply(1:length(meanL), function(i) rnorm(nsamp[i], meanL[i], sdL[i])))
#' Ages <- rep(1:length(N), nsamp)
#' data <- data.frame(Age=Ages, Length=Length)
#' 
#' # Return an OM with cpars populated
#' OM@nsim <- 48
#' newOM <- Growth2OM(data, OM)
#' 
#' # Return a data.frame
#' estPars <- Growth2OM(data)
#' 
Growth2OM <- function(data=NULL, OM=NULL, nsim=48, seed=101, plot=TRUE, msg=TRUE) {
  om <- TRUE
  if (all(!inherits(OM,'OM'))) om <- FALSE
  if (om) {
    set.seed(OM@seed)
  } else {
    set.seed(seed)
  }
  if (!inherits(data, 'data.frame')) stop("data must be a data.frame")
  if (!all(c("Age", "Length") %in% colnames(data))) stop("data must have columns: Age, Length")
  
  Length <- av.len <- n <- st.age <- NULL # hacks from cran checks
  
  VB <- function(Linf,k,t0,age) Linf*(1- exp(-k*(age-t0)))
  nsamp <- ifelse(om, OM@nsim, nsim) 
  if (msg) message("Running ", nsamp, " bootstrap samples")
  options(warn=-1)
  boot1 <- boot::boot(data, function(DF, i) {
    tryCatch(coef(nls(Length ~ SSasympOff(Age, Asym, lrc, c0), data = DF[i,])),
             error = function(e) c(Asym = NA, lrc = NA, c0 = NA), silent=TRUE)
  }, R = nsamp) 
  ests <- data.frame(Linf=boot1$t[,1], K=exp(boot1$t[,2]), t0=boot1$t[,3])
  if (sum(is.na(ests)) >0 ) {
    if (msg)  message("Some NAs. Running again with more sims")
    boot1 <- boot::boot(data, function(DF, i) {
      tryCatch(coef(nls(Length ~ SSasympOff(Age, Asym, lrc, c0), data = DF[i,])),
               error = function(e) c(Asym = NA, lrc = NA, c0 = NA), silent=TRUE)
    }, R = nsamp*3) # run three times as many and remove NAs
    ests <- data.frame(Linf=boot1$t[,1], K=exp(boot1$t[,2]), t0=boot1$t[,3])
  }
  options(warn=0)
  ests <- ests[which(apply(!is.na(ests), 1, prod)==1),]
  if (nrow(ests) >= nsamp) {
    ests <- ests[1:nsamp,]  
  } else {
    if (msg) message("Note: Less than ", nsamp, "estimates")
  }
  
  if (plot) {
    op <- par(mfrow=c(1,1))
    on.exit(par(op))
    ages <- 0:ceiling(max(data$Age))
    fitVBs <- sapply(1:nsamp, function(x) VB(ests$Linf[x], ests$K[x], ests$t0[x], ages))
    xlim <- range(ages)
    plot(data$Age, data$Length, xlab="Age", ylab="Length", bty="l", las=1, 
         xlim=xlim, pch=16, ylim=c(0, max(data$Length, na.rm=TRUE)))
    matplot(ages, fitVBs, type="l", col=1:nsamp, add=TRUE)
    
    mod <- nls(Length ~ SSasympOff(Age, Asym, lrc, c0), data = data)
    bestests <- round(c(coef(mod)[1],exp(coef(mod)[2]), coef(mod)[3]),2)

    leg <-c(paste0('Linf = ', bestests[1]), 
            paste0('K = ', bestests[2]), 
            paste0('t0 = ', bestests[3]))
    legend('topleft', legend=leg, bty="n")
    lines(ages, VB(bestests[1], bestests[2], bestests[3], ages), lwd=3)
  }
  
  Ages <- 0:max(data$Age)
  meanCurve <- data.frame(Age=Ages, 
                          av.len=mean(ests$Linf, na.rm=TRUE) * 
                            (1-exp(-mean(ests$K, na.rm=TRUE) * 
                                     (Ages-mean(ests$t0, na.rm=TRUE)))))
  data$st.age <- round(data$Age,0)
  data <- dplyr::left_join(data, meanCurve, by="Age")
  
  estcv <- data %>% dplyr::group_by(st.age) %>% 
    dplyr::summarize(cv=mean(sd(Length, na.rm=TRUE)/av.len, na.rm=TRUE),
                     n=length(Length))
  
  LenCV <- estcv %>% dplyr::filter(n>10) %>% 
    dplyr::summarise(low=quantile(cv, 0.05, na.rm=TRUE),
                     upp=quantile(cv, 0.95, na.rm=TRUE))
  
  if (om) {
    if (msg) message('Returning OM')
    parname <- c("Linf", "K", "t0")
    matchname <- names(OM@cpars) %in% parname
    if (length(matchname)<1) matchname <- FALSE
    matchname <- matchname[matchname]
    if (msg) {
      if (sum(matchname)) {
        for (x in 1:length(matchname)) {
          message('Note: ', parname[x], ' already in cpars. Replacing with new estimates')
        }
      }
    }
    OM@LenCV <- as.numeric(LenCV )
    OM@Linf <- c(0,0)
    OM@K <- c(0,0)
    OM@t0 <- c(0,0)
    OM@cpars$Linf <- ests$Linf
    OM@cpars$K <- ests$K
    OM@cpars$t0 <- ests$t0
    return(OM)
  } else {
    if (msg) message('Returning data.frame')
    ests$LenCV_low <- rep(LenCV[1], nrow(ests))
    ests$LenCV_upp <- rep(LenCV[2], nrow(ests))
    return(ests)
  }
}


#' Estimate length-weight parameters from data
#'
#' Function estimates alpha and beta parameter from length-weight data and populates
#' the relevant slots in the OM
#' 
#' @param data A data frame with columnns 'Length' and 'Weight' with numeric data 
#' @param OM An object of class `OM`
#' @param plot Logical. Show plot of data and fit?
#'
#' @return An object of class `OM` with `OM@a` and `OM@b` slots populated
#' @export
#'
LW2OM <- function(data=NULL, OM=NULL, plot=TRUE) {
  if (!inherits(OM, 'OM')) stop("OM must be class OM")
  if (!inherits(data, 'data.frame')) stop("data must be a data.frame")
  if (!all(c("Weight", "Length") %in% colnames(data))) stop("data must have columns: Length, Weight")
  
  Length <- Weight <- NULL # hacks from cran checks
  
  dat <- dplyr::select(data, Length, Weight)
  options(warn=-1)
  dat$Length <- as.numeric(dat$Length)
  dat$Weight <- as.numeric(dat$Weight)
  options(warn=0)
  dat <- dat[!as.logical(apply(is.na(dat), 1, sum)),]
  mod <- lm(log(Weight)~log(Length), data=dat)
  OM@a <- as.numeric(exp(coef(mod)[1]))
  OM@b <- as.numeric(coef(mod)[2])
  
  if (plot) {
    op <- par(mfrow=c(1,1))
    on.exit(par(op))
    plot(dat$Length, dat$Weight, xlab="Length", ylab="Weight", las=1, pch=16)
    Lengths <- seq(min(dat$Length), max(dat$Length), length.out=100)
    lines(Lengths,  OM@a*Lengths^  OM@b , lwd=2)
    
    leg <-c(paste0('a = ', signif(OM@a,2)), 
            paste0('b = ', signif(OM@b,2)))
    legend('topleft', legend=leg, bty="n")
    
  }
  return(OM)

}

# 
# OM <- DLMtool::testOM 
# OM@nsim <- 2
# Hist <- runMSE(OM, Hist=TRUE)
# N <- Hist@AtAge$Nage[1,,1] * Hist@AtAge$Select[1,,1]
# meanL <- Hist@AtAge$Length[1,,1]
# sdL <- Hist@AtAge$Length[1,,1] * 0.1
# nsamp <- ceiling(N/sum(N) * 1000)
# Length <- unlist(sapply(1:length(meanL), function(i) rnorm(nsamp[i], meanL[i], sdL[i])))
# prob <- 1/(1 + exp(-log(19) * ((Length - OM@L50[1])/(OM@L50_95[1])))) # prob mature
# Mature <- rbinom(n=sum(nsamp), size=1, prob=prob)
# data <- data.frame(Length=Length, Mature=Mature)
#  
# 
# Mat2OM <- function(data=NULL, OM=NULL, plot=TRUE) {
#   if (!inherits(OM, 'OM')) stop("OM must be class OM")
#   if (!inherits(data, 'data.frame')) stop("data must be a data.frame")
#   if (!all(c("Length", "Mature") %in% colnames(data))) stop("data must have columns: Length, Mature")
#   
#   
#   plot(data, alpha=0.1, pch=16)
#   
#   test <- glm(Mature~Length, family=binomial(link="logit"), data=data)
#   coef(test)
#   summary(test)
#   plot(meanL, predict(test, newdata=data.frame(Length=meanL), type="response"))
#   
#  
#   
# }
