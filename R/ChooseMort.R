#' Manually map natural mortality at age or size.
#' 
#' Interactive plot which allows users to specify M by age or size class
#' 
#' @param OM An object of class 'OM'
#' @param type A character string - is M to be mapped by 'age' or 'length'?
#' @param x Optional vector for x-axis
#' @param y Optional vector for y-axis
#' @author A. Hordyk
#' @export ChooseM
ChooseM <- function(OM, type=c("age", "length"), x=NULL, y=NULL) {
 
  type <- match.arg(type)
  if (class(OM) != "OM") stop("First argument must be object of class 'OM'", call.=FALSE)
  
  ylab <- "Natural mortality (M)"
  if (is.null(y)) y <- seq(0, 1, 0.05)
  if (type == "age") {
    if (is.null(x)) x <- 1:OM@maxage
    xlab <- "Age Class"
    runSketch <- SketchM(x, y, xlab, ylab)  
    M1 <- approx(x = runSketch[,1], y = runSketch[,2], xout=x, method = "linear")$y  # linear interpolation
    M2 <- approx(x = runSketch[,1], y = runSketch[,3], xout=x, method = "linear")$y  # linear interpolation
    OM@M <- M1 
    OM@M2 <- M2 
    return(OM)
  }
  if (type == "length") {
    if (is.null(x)) x <- seq(from=0, to=max(OM@Linf), 5)
    xlab <- "Length Class"
    runSketch <- SketchM(x, y, xlab, ylab)  
    M1 <- approx(x = runSketch[,1], y = runSketch[,2], xout=x, method = "linear")$y  # linear interpolation
    M2 <- approx(x = runSketch[,1], y = runSketch[,3], xout=x, method = "linear")$y  # linear interpolation
    M_at_Length <- data.frame(Lens=x, M1, M2)
    OM@cpars[["M_at_Length"]] <- M_at_Length
    return(OM)
  }
}


SketchM <- function(x, y, xlab, ylab) {
  x1 <- x[seq(1, by=2, to=length(x))]
  x2 <- x[seq(2, by=2, to=length(x))]
  grd1 <- expand.grid(x1, y)
  grd2 <- expand.grid(x2, y)
  grd3 <- expand.grid(x, y)
  Xs <- grd3[, 1]
  Ys <- grd3[, 2]
  
  op <- par(mfrow=c(1,1), las=1, oma=c(1,2,1,2))
  plot(range(x), range(y), type="n", xlab = xlab, ylab = ylab, cex.lab = 1.5)
  points(grd1[, 1], grd1[, 2], col = "black", pch = 19, cex = 0.2)
  points(grd2[, 1], grd2[, 2], col = "darkgrey", pch = 19, cex = 0.2)
  axis(side=3)
  axis(side=4)  
  mtext(side = 3, "Click 'Finish' or press 'Escape' to Finish.", xpd = NA, 
        cex = 1.25, line=2.5)
  line1 <- "Use mouse to select points on the grid"
  line2 <- "First and last age/size classes must be selected."
  line3 <- "Select two points for each age/size class to represent range of uncertainty"
  message(line1, "\n", line2, "\n", line3, "\n")
 
  
  out <- NULL
  out <- identifyPch(x = Xs, y = Ys, tolerance = 0.1)
  while (is.null(dim(out))) {
    message("Must choose more than one point")
    out <- identifyPch(x = Xs, y = Ys, tolerance = 0.1)
  }
  while (min(out[, 1]) != x[1]) {
    message("Choose point(s) for first age/size class")
    dat <- rbind(out, identifyPch(x = Xs, y = Ys, tolerance = 0.1))
    out <- dat[order(dat[, 1]), ]
  }
  while (max(out[, 1]) != max(x)) {
    message("Choose point(s) for last age/size class")
    dat <- rbind(out, identifyPch(x = Xs, y = Ys, tolerance = 0.1))
    out <- dat[order(dat[, 1]), ]
  }
  ord <- order(out[, 1], out[, 2])
  out <- out[ord, ]
  yrs <- unique(out[, 1])
  mat <- matrix(NA, nrow = length(yrs), ncol = 3)
  mat[, 1] <- yrs
  ind <- which(!duplicated(out[, 1]))
  mat[, 2:3] <- out[ind, 2]
  for (X in seq_along(yrs)) {
    chk <- out[, 1] %in% yrs[X]
    ind <- range(which(chk))
    if (sum(chk) > 1) {
      mat[X, 2:3] <- out[ind, 2]
    }
  }
  
  lines(mat[, 1], mat[, 2])
  lines(mat[, 1], mat[, 3])
  
  par(op)
  return(mat)
}

#' Plot M-at-Age and Size
#'
#' @param Stock An object of class 'Stock' or 'OM' 
#' @param nsim The number of simulations to plot
#'
#' @author A. Hordyk
#' @export
#'
#' @examples plotM(Albacore)
plotM <- function(Stock, nsim=5) {
  if (class(Stock) != "Stock" && class(Stock) != "OM") stop("Must supply object of class 'Stock' or 'OM'")
  
  nyears <- 30
  proyears <- 30
  SampCpars <- list() # empty list 
  if (class(Stock) == "OM") {
    # custom parameters exist - sample and write to list
    if(length(Stock@cpars)>0){
      ncparsim<-cparscheck(Stock@cpars)   # check each list object has the same length and if not stop and error report
      SampCpars <- SampleCpars(Stock@cpars, nsim) 
    }
    nyears <- Stock@nyears 
    proyears <- Stock@proyears
    Stock@nsim <- nsim
  }
  
  StockPars <- SampleStockPars(Stock, nsim, nyears, proyears, SampCpars)
  # Assign Stock pars to function environment
  for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])
  
  M_at_age <- StockPars$M_ageArray
  Len_at_age <- StockPars$Len_age
  Wt_at_age <- StockPars$Wt_age
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  

  par(mfrow=c(3,3), bty="l", las=1, mar=c(3,3,2,1), oma=c(2,2,0,0))
  ylim <- c(0, max(M_at_age))
  lwd <- 2
  matplot(t(M_at_age[,,1]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  mtext(side=2, "M", line=3)
  # mtext(side=1, "Age", line=2.5)
  matplot(t(Len_at_age[,,1]), t(M_at_age[,,1]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  # mtext(side=1, "Length", line=2.5)
  title("First Historical Year")
  matplot(t(Wt_at_age[,,1]), t(M_at_age[,,1]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  # mtext(side=1, "Weight", line=2.5)
  
  matplot(t(M_at_age[,,nyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  mtext(side=2, "M", line=3)
  # mtext(side=1, "Age", line=2.5)
  matplot(t(Len_at_age[,,nyears]), t(M_at_age[,,nyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  # mtext(side=1, "Length", line=2.5)
  title("Last Historical Year")
  matplot(t(Wt_at_age[,,nyears]), t(M_at_age[,,nyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  # mtext(side=1, "Weight", line=2.5)
  
  matplot(t(M_at_age[,,proyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  mtext(side=2, "M", line=3)
  mtext(side=1, "Age", line=2.5)
  matplot(t(Len_at_age[,,proyears]), t(M_at_age[,,proyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  mtext(side=1, "Length", line=2.5)
  title("Last Projection Year")
  matplot(t(Wt_at_age[,,proyears]), t(M_at_age[,,proyears]), type="l", xlab="", ylab="", ylim=ylim, lwd=lwd)
  mtext(side=1, "Weight", line=2.5)
  
  
}
