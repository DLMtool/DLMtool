#' Check Convergence
#' 
#' Have I undertaken enough simulations (nsim)? Has my MSE converged on stable
#' (reliable) peformance metrics?
#' 
#' Performance metrics are plotted against the number of simulations. Convergence diagonostics 
#' are calculated over the last `ref.it` (default = 20) iterations. The convergence diagnostics are:
#' \enumerate{
#'   \item Is the order of the MPs stable over the last `ref.it` iterations? 
#'   \item Is the average difference in performance statistic over the last `ref.it` iterations < `thresh`?
#' } 
#' 
#' By default three commonly used performance metrics are used: 
#' \enumerate{
#'   \item Average Yield Relative to Reference Yield 
#'   \item Probability Spawning Biomass is above 0.1BMSY  
#'   \item Probability Average Annual Variability in Yield is < 20 per cent
#' }
#' Additional or alternative performance metrics objects can be supplied. Advanced users can develop their own performance metrics.
#'  
#' @param MSEobj An MSE object of class \code{'MSE'}
#' @param PMs A list of PM objects
#' @param maxMP Maximum number of MPs to include in a single plot
#' @param thresh The convergence threshold. Maximum root mean square deviation over the last `ref.it` iterations
#' @param ref.it The number of iterations to calculate the convergence statistics. For example,
#' a value of 20 means convergence diagnostics are calculated over last 20 simulations
#' @param inc.leg Logical. Should the legend be displayed?
#' @param all.its Logical. Plot all iterations? Otherwise only (nsim-ref.it):nsim
#' @param nrow Numeric. Optional. Number of rows
#' @param ncol Numeric. Optional. Number of columns
#' 
#' @templateVar url checking-convergence
#' @templateVar ref NULL
#' @template userguide_link
#' 
#' @examples 
#' \dontrun{
#' MSE <- runMSE()
#' Converge(MSE)
#' }
#' 
#' @author A. Hordyk
#' @export 
#' 
Converge <- function(MSEobj, PMs=list(Yield, P10, AAVY), maxMP=15, thresh=0.5, ref.it=20,
                     inc.leg=FALSE, all.its=FALSE, nrow=NULL, ncol=NULL) {
  
  
  if(class(MSEobj) != "MSE") stop("MSEobj must be object of class 'MSE'", call.=FALSE)
  if(MSEobj@nMPs <2) stop("Converge requires more than 1 MP in the MSE object", call.=FALSE)
  
  nPMs <- length(PMs)
  if (is.null(ncol)) ncol <- floor(sqrt(nPMs))
  if (is.null(nrow)) nrow <- ceiling(nPMs)/ncol
  if (ncol * nrow < nPMs) stop("ncol x nrow must be > length(PMs)")
 
  if (MSEobj@nMPs > maxMP) {
    nplots <- ceiling(MSEobj@nMPs /maxMP)
  } else{
    nplots <- 1
  }
  
  nsim <- MSEobj@nsim
  if (nsim-ref.it < 1) {
    ref.it.new <- nsim - 1
    message("nsim (", nsim, ") - ref.it (", ref.it,") < 1. Setting ref.it to ", ref.it.new, "\n")
    ref.it <- ref.it.new 
  }
  
  message("Checking if order of MPs is changing in last ", ref.it, " iterations\n")
  message("Checking average difference in PM over last ", ref.it, " iterations is > ", thresh, "\n")

  SwitchOrd <- vector(mode = "list", length = nPMs)
  NonCon <- vector(mode = "list", length = nPMs)
  PMName <- vector("character", length=nPMs)
  
  getPalette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  st <- 1 
  end <- min(maxMP, MSEobj@nMPs)
  it <- MP <- NULL # cran check hacks 
  
  for (nplot in 1:nplots) {
    
    message('Plotting MPs ', st, ' - ', end)
    subMSE <- Sub(MSEobj, MPs=MSEobj@MPs[st:end])
    
    nMPs <-  subMSE@nMPs
    values <- rep(c("solid", "dashed", "dotted"),  nMPs)
    values <- values[1:nMPs]
    
    plist <- list()
    for (xx in 1:nPMs) {
      PMval <- PMs[[xx]](subMSE)
      PMName[xx] <- PMval@Name
      PMval@Prob[!is.finite(PMval@Prob)] <- 0
      cum_mean <- apply(PMval@Prob, 2, cumsum)/apply(PMval@Prob, 2, seq_along) * 100
      vals <- as.vector(cum_mean) 
      mp <- rep(subMSE@MPs, each=subMSE@nsim)
      
      df <- data.frame(it=1:subMSE@nsim, vals, MP=mp, name=PMval@Name)
      if (!all.its) df <- subset(df, it %in% (nsim-ref.it+1):nsim)
  
      p <- ggplot2::ggplot(df, ggplot2::aes(x=it, y=vals, color=MP, linetype=MP)) + 
        ggplot2::scale_linetype_manual(values = values) +
        ggplot2::scale_color_manual(values=getPalette(nMPs)) + 
        ggplot2::geom_line() +
        ggplot2::theme_classic() +
        ggplot2::labs(x="# Simulations", y=PMval@Caption) +
        ggplot2::ggtitle(PMval@Name) +
        ggplot2::geom_vline(xintercept=MSEobj@nsim-ref.it+1, color="darkgray", linetype="dashed") +
        ggplot2::geom_vline(xintercept=MSEobj@nsim, color="darkgray", linetype="dashed") +
        ggplot2::coord_cartesian(xlim = c(min(df$it), max(df$it))) 
      if (!inc.leg) p <- p + ggplot2::theme(legend.position = "none")
      
      # check positions & convergence
      ords <- apply(cum_mean[(MSEobj@nsim-ref.it+1):MSEobj@nsim,], 1, order, decreasing=FALSE)
      rownames(ords) <- subMSE@MPs
      mat <- matrix(subMSE@MPs[ords], nrow=subMSE@nMPs, ncol=ref.it)
      tab <- table(unlist(apply(mat, 1, unique)))
      SwitchOrd[[xx]] <- append(SwitchOrd[[xx]], rownames(tab)[which(tab > 1)])
      NonCon[[xx]] <- append(NonCon[[xx]], subMSE@MPs[apply(cum_mean, 2, Chk, MSEobj=MSEobj, thresh=thresh, ref.it)])
      
      noncoverg <- unique(c(SwitchOrd[[xx]], NonCon[[xx]]))
      if (length(noncoverg)>0) {
        df2 <- subset(df, MP %in% noncoverg)
        if (dim(df2)[1] > 0) {
          p <- p + 
            ggrepel::geom_text_repel(
              data = subset(df2, it == max(it)),
              ggplot2::aes(label = MP),
              size = 4,
              segment.color = "grey",
              direction = "both",
              show.legend = FALSE)
        }
  
      }
    
      plist[[xx]] <- p

    }
    
    if (inc.leg) join_plots(plist, nrow=nrow, ncol=ncol, position="right")
    if (!inc.leg) gridExtra::grid.arrange(grobs=plist, nrow=nrow, ncol=ncol) 
    
    st <- st + maxMP 
    end <- min(end + maxMP, MSEobj@nMPs)
  
  }
  
  ## Report Performance 
  for (x in 1:nPMs) {
    SwitchOrds <- SwitchOrd[[x]]
    NonCons <- NonCon[[x]]
    if (length(SwitchOrds)>0  | length(NonCons)>0) {
      message("\n", PMName[x], "\n")
      if (length(SwitchOrds)>0) {
        message("Order over last ", ref.it, ' iterations is not consistent for:\n ', paste(SwitchOrds, ""), ' \n')
      }
      if (length(NonCons)>0) {
        message("Mean difference over last ", ref.it, " iterations is > ", thresh, " for:\n",
                paste(NonCons, ""), '\n')
      }
    } 
  }
}




Chk <- function(X, MSEobj, thresh, ref.it) {
  L <- length(X)
  Y <- min(MSEobj@nsim, ref.it) + 1 
  x <- X[(L-Y):L]
  
  # root mean square deviation in last ref.it iterations is greater than thresh
  sqrt((sum((x-mean(x))^2))/(L-1)) > thresh
}

  
#   nm <- MSEobj@nMPs
#   nsim <- MSEobj@nsim
#   proyears <- MSEobj@proyears
#   
#   Yd <- CumlYd <- array(NA, c(nm, nsim))
#   P10 <- CumlP10 <- array(NA, c(nm, nsim))
#   P50 <- CumlP50 <- array(NA, c(nm, nsim))
#   P100 <- CumlP100 <- array(NA, c(nm, nsim))
#   POF <- CumlPOF <- array(NA, c(nm, nsim))
#   yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
#   RefYd <- MSEobj@OM$RefY
#   
#   for (m in 1:nm) {
#     Yd[m, ] <- round(apply(MSEobj@C[, m, yind], 1, mean, na.rm = T)/RefYd *  100, 1)
#     POF[m, ] <- round(apply(MSEobj@F_FMSY[, m, ] >= 1, 1, sum, na.rm = T)/proyears * 100, 1)
#     P10[m, ] <- round(apply(MSEobj@B_BMSY[, m, ] <= 0.1, 1, sum, na.rm = T)/proyears * 100, 1)
#     P50[m, ] <- round(apply(MSEobj@B_BMSY[, m, ] <= 0.5, 1, sum, na.rm = T)/proyears * 100, 1)
#     P100[m, ] <- round(apply(MSEobj@B_BMSY[, m, ] <= 1, 1, sum, na.rm = T)/proyears * 100, 1)
#     CumlYd[m, ] <- cumsum(Yd[m, ])/seq_along(Yd[m, ])  #/ mean(Yd[m,], na.rm=TRUE) 
#     CumlPOF[m, ] <- cumsum(POF[m, ])/seq_along(POF[m, ])  # / mean(POF[m,], na.rm=TRUE)
#     CumlP10[m, ] <- cumsum(P10[m, ])/seq_along(P10[m, ])  # / mean(P10[m,], na.rm=TRUE)
#     CumlP50[m, ] <- cumsum(P50[m, ])/seq_along(P50[m, ])  # / mean(P50[m,], na.rm=TRUE)
#     CumlP100[m, ] <- cumsum(P100[m, ])/seq_along(P100[m, ])  # / mean(P100[m,], na.rm=TRUE)
#   }
#   
#   
#   
#   # CumlYd[is.nan(CumlYd)] <- 1 CumlPOF[is.nan(CumlPOF)] <- 1
#   # CumlP10[is.nan(CumlP10)] <- 1 CumlP50[is.nan(CumlP50)] <- 1
#   # CumlP100[is.nan(CumlP100)] <- 1
#   
#   if (Plot) {	  
#     op <- par(mfrow = c(2, 3), cex.axis = 1.5, cex.lab = 1.7, oma = c(1, 1, 0, 0), mar = c(5, 5, 1, 1), bty = "l")
#     matplot(t(CumlYd), type = "l", xlab = "Iteration", ylab = "Rel. Yield")
#     matplot(t(CumlPOF), type = "l", xlab = "Iteration", ylab = "Prob. F/FMSY > 1")
#     matplot(t(CumlP10), type = "l", xlab = "Iteration", ylab = "Prob. B/BMSY < 0.1")
#     matplot(t(CumlP50), type = "l", xlab = "Iteration", ylab = "Prob. B/BMSY < 0.5")
#     matplot(t(CumlP100), type = "l", xlab = "Iteration", ylab = "Prob. B/BMSY < 1")
#     
#   }
#   
#   Chk <- function(X) {
#     # checks if difference in last 20 iterations is greater than thresh
#     L <- length(X)
#     Y <- 1:min(nsim, 20)
#     
#     # return(mean(abs((X[L-1:10] - X[L]))/X[L], na.rm=TRUE) > thresh)
#     return(mean(abs((X[L - Y] - X[L])), na.rm = TRUE) > thresh)
#   }
# 
#   NonCon <- sort(unique(c(which(apply(CumlYd, 1, Chk)), 
#                           which(apply(CumlPOF, 1, Chk)), 
#                           which(apply(CumlP10, 1, Chk)), 
#                           which(apply(CumlP50, 1, Chk)), 
#                           which(apply(CumlP100, 1, Chk)))))
#   
#   if (length(NonCon) > 0) {
#     if (Plot) {
#       plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", 
#            ylab = "")
#       text(0.5, 0.5, "Some MPs have not converged", cex = 1)
#       # ev.new()
#       par(mfrow = c(2, 3), cex.axis = 1.5, cex.lab = 1.7, oma = c(1, 
#                                                                   1, 0, 0), mar = c(5, 5, 1, 1), bty = "l")
#       if (length(NonCon) > 1) {
#         matplot(t(CumlYd[NonCon, ]), type = "l", xlab = "Iteration",  ylab = "Rel. Yield", lwd = 2, lty=1:length(NonCon))
#       } else matplot((CumlYd[NonCon, ]), type = "l", xlab = "Iteration",  ylab = "Rel. Yield", lwd = 2, lty=1:length(NonCon))
#       if (length(NonCon) > 1) {
#         matplot(t(CumlPOF[NonCon, ]), type = "l", xlab = "Iteration", ylab = "Prob. F/FMSY > 1", lwd = 2, lty=1:length(NonCon))
#       } else matplot((CumlPOF[NonCon, ]), type = "l", xlab = "Iteration", ylab = "Prob. F/FMSY > 1", lwd = 2, lty=1:length(NonCon))  
#       if (length(NonCon) > 1) {
#         matplot(t(CumlP10[NonCon, ]), type = "l", xlab = "Iteration", ylab = "Prob. B/BMSY < 0.1", lwd = 2, lty=1:length(NonCon))
#       } else matplot((CumlP10[NonCon, ]), type = "l", xlab = "Iteration", ylab = "Prob. B/BMSY < 0.1", lwd = 2, lty=1:length(NonCon))
#       if (length(NonCon) > 1) {
#         matplot(t(CumlP50[NonCon, ]), type = "l", xlab = "Iteration", ylab = "Prob. B/BMSY < 0.5", lwd = 2, lty=1:length(NonCon))
#       } else matplot((CumlP50[NonCon, ]), type = "l", xlab = "Iteration", ylab = "Prob. B/BMSY < 0.5", lwd = 2, lty=1:length(NonCon))
#       if (length(NonCon) > 1) {
#         matplot(t(CumlP100[NonCon, ]), type = "l", xlab = "Iteration", ylab = "Prob. B/BMSY < 1", lwd = 2, lty=1:length(NonCon))
#       } else matplot((CumlP100[NonCon, ]), type = "l", xlab = "Iteration", ylab = "Prob. B/BMSY < 1", lwd = 2, lty=1:length(NonCon))
#       
#       legend(nsim * 1.25, 50, legend = MSEobj@MPs[NonCon], col = 1:length(NonCon), 
#              bty = "n", xpd = NA, lty = 1:length(NonCon), lwd = 2, cex = 1.25)
#     }
#     
#     message("Some MPs may not have converged in ", nsim, " iterations (threshold = ", thresh, "%)")
#     message("MPs are: ", paste(MSEobj@MPs[NonCon], " "))
#     message("MPs #: ", paste(NonCon, " "))
#     return(data.frame(Num = NonCon, MP = MSEobj@MPs[NonCon]))
#   }
#   if (length(NonCon) == 0) {
#     if (Plot) {
#       plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", 
#            ylab = "")
#       text(0.5, 0.5, "All MPs converged", cex = 1)
#     }
#     message("All MPs appear to have converged in ", nsim, " iterations (threshold = ", 
#             thresh, "%)")
#   }
#   par(op)
# }







# Kobe plot Kplot<-function(MSEobj,maxsim=60,nam=NA){
# nr<-floor((MSEobj@nMPs)^0.5) nc<-ceiling((MSEobj@nMPs)/nr)

# FMSYr<-quantile(MSEobj@F_FMSY,c(0.001,0.90),na.rm=T)
# BMSYr<-quantile(MSEobj@B_BMSY,c(0.001,0.975),na.rm=T)

# #dev.new2(width=nc*3,height=nr*3.6)
# par(mfrow=c(nr,nc),mai=c(0.45,0.45,0.45,0.01),omi=c(0.45,0.3,0.35,0.01))

# colsse<-rainbow(MSEobj@proyears,start=0.63,end=0.95)[1:MSEobj@proyears]
# colsse<-makeTransparent(colsse,95)

# for(mm in 1:MSEobj@nMPs){
# plot(c(MSEobj@B_BMSY[1,mm,1],MSEobj@B_BMSY[1,mm,2]),
# c(MSEobj@F_FMSY[1,mm,1],MSEobj@F_FMSY[1,mm,2]),xlim=BMSYr,ylim=FMSYr,
# col=colsse[1],type='l')

# OO<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]<1&MSEobj@F_FMSY[,mm,MSEobj@proyears]>1,na.rm=T)/MSEobj@nsim*100,1)
# OU<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]>1&MSEobj@F_FMSY[,mm,MSEobj@proyears]>1,na.rm=T)/MSEobj@nsim*100,1)
# UO<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]<1&MSEobj@F_FMSY[,mm,MSEobj@proyears]<1,na.rm=T)/MSEobj@nsim*100,1)
# UU<-round(sum(MSEobj@B_BMSY[,mm,MSEobj@proyears]>1&MSEobj@F_FMSY[,mm,MSEobj@proyears]<1,na.rm=T)/MSEobj@nsim*100,1)

# #alp<-80
# #polygon(c(1,-1000,-1000,1),c(1,1,1000,1000),col=makeTransparent('orange',alp),border=makeTransparent('orange',alp))
# #polygon(c(1,1000,1000,1),c(1,1,1000,1000),col=makeTransparent('yellow',alp),border=makeTransparent('yellow',alp))
# #polygon(c(1,-1000,-1000,1),c(1,1,-1000,-1000),col=makeTransparent('yellow',alp),border=makeTransparent('yellow',alp))
# #polygon(c(1,1000,1000,1),c(1,1,-1000,-1000),col=makeTransparent('green',alp),border=makeTransparent('yellow',alp))


# abline(h=1,col='grey',lwd=3) abline(v=1,col='grey',lwd=3)
# #abline(v=c(0.1,0.5),col='grey',lwd=2) rng<-1:min(maxsim,MSEobj@nsim)
# for(i in rng){ for(y in 1:(MSEobj@proyears-1)){
# lines(c(MSEobj@B_BMSY[i,mm,y],MSEobj@B_BMSY[i,mm,y+1]),
# c(MSEobj@F_FMSY[i,mm,y],MSEobj@F_FMSY[i,mm,y+1]),
# col=colsse[y],lwd=1.6) } }

# points(MSEobj@B_BMSY[rng,mm,1],MSEobj@F_FMSY[rng,mm,1],pch=19,cex=0.8,col=colsse[1])
# points(MSEobj@B_BMSY[rng,mm,MSEobj@proyears],MSEobj@F_FMSY[rng,mm,MSEobj@proyears],pch=19,cex=0.8,col=colsse[MSEobj@proyears])

# if(mm==1)legend('right',c('Start','End'),bty='n',text.col=c(colsse[1],colsse[MSEobj@proyears]),pch=19,col=c(colsse[1],colsse[MSEobj@proyears]))
# legend('topleft',paste(OO,'%',sep=''),bty='n',text.font=2)
# legend('topright',paste(OU,'%',sep=''),bty='n',text.font=2)
# legend('bottomleft',paste(UO,'%',sep=''),bty='n',text.font=2)
# legend('bottomright',paste(UU,'%',sep=''),bty='n',text.font=2)

# mtext(MSEobj@MPs[mm],3,line=0.45) }
# mtext('B/BMSY',1,outer=T,line=1.4) mtext('F/FMSY',2,outer=T,line=0.2)
# if(is.na(nam))mtext(deparse(substitute(MSEobj)),3,outer=T,line=0.25,font=2)
# if(!is.na(nam))mtext(MSEobj@Name,3,outer=T,line=0.25,font=2) }





# Value of information analysis
# Value of information




# Manipulation of MSE Object
# --------------------------------------------------- Subset the MSE
# object by particular MPs (either MP number or name), or particular
# simulations

#' Utility functions for MSE objects
#' 
#' @param MSEobj A MSE object. For `updateMSE`, a MSE object from a previous version of 
#' DLMtool. Also works with Stock, Fleet, Obs, Imp, and Data objects.
#' @param MSEobjs A list of MSE objects. Must all have identical operating
#' model and MPs. MPs which don't appear in all MSE objects will be dropped.
#' @return An object of class \code{MSE}
#' @examples 
#' # An example of joinMSE
#' \dontrun{
#' OM1 <- DLMtool::testOM
#' MSE1 <- runMSE(OM1) 
#' OM2 <- OM1 
#' OM2@seed <- OM1@seed + 1
#' MSE2 <- runMSE(OM2)
#' MSE <- joinMSE(list(MSE1, MSE2))
#' MSE@nsim
#' }
#' @author A. Hordyk
#' @describeIn checkMSE Check that an MSE object includes all slots in the latest version of DLMtool
#' Use `updateMSE` to update the MSE object
#' @export checkMSE
checkMSE <- function(MSEobj) {
  nms <- slotNames(MSEobj)
  errs <- NULL
  for (x in seq_along(nms)) {
    chk <- try(slot(MSEobj, nms[x]), silent=TRUE)
    if (class(chk) == "try-error") errs <- c(errs, x)
  }
  if (length(errs) > 0) {
    message("MSE object slots not found: ", paste(nms[errs], ""))
    stop("slot names of MSEobj don't match MSE object class. Try use `updateMSE`", call.=FALSE)
    return(FALSE)
  }
  return(TRUE)
}



#' Subset MSE object by management procedure (MP) or simulation.
#' 
#' Subset the MSE object by particular MPs (either MP number or name), or
#' particular simulations, or a subset of the projection years (e.g., 1: <
#' projection years).
#' 
#' @param MSEobj A MSE object.
#' @param MPs A vector MPs names or MP numbers to subset the MSE object.
#' Defaults to all MPs.
#' @param sims A vector of simulation numbers to subset the MSE object. Can
#' also be a logical vector. Defaults to all simulations.
#' @param years A numeric vector of projection years. Should start at 1 and
#' increase by one to some value equal or less than the total number of
#' projection years.
#' 
#' @templateVar url subsetting-the-mse-object
#' @templateVar ref NULL
#' @template userguide_link
#' 
#' @author A. Hordyk
#' @examples
#' \dontrun{
#' MSE <- runMSE() 
#' MSE_1 <- Sub(MSE, MPs=1:2)
#' MSE_1@MPs
#' MSE_2 <- Sub(MSE, sims=1:10)
#' MSE_2@nsim
#' }
#' @export Sub
Sub <- function(MSEobj, MPs = NULL, sims = NULL, years = NULL) {
  
  checkMSE(MSEobj) # check that MSE object contains all slots 
  
  Class <- class(MPs) 
  if (Class == "NULL") subMPs <- MSEobj@MPs
  if (Class == "integer" | Class == "numeric") subMPs <- MSEobj@MPs[as.integer(MPs)]
  if (Class == "character") subMPs <- MPs
  if (Class == "factor") subMPs <- as.character(MPs)
  subMPs <- subMPs[!is.na(subMPs)]
  
  SubMPs <- match(subMPs, MSEobj@MPs)  #  which(MSEobj@MPs %in% subMPs)
  if (any(is.na(SubMPs))) {
    missing <- subMPs[is.na(SubMPs)]
    stop(paste0(missing, collapse=','), ' not found in MSE object', call.=FALSE)
  }
  
  not <- (subMPs %in% MSEobj@MPs)  # Check for MPs misspelled
  ind <- which(not == FALSE)
  newMPs <- MSEobj@MPs[SubMPs]
  if (length(SubMPs) < 1) stop("MPs not found")
  if (length(ind > 0)) {
    message("Warning: MPs not found - ", paste0(subMPs[ind], " "))
    message("Subsetting by MPs: ", paste0(newMPs, " "))
  }
  
  
  ClassSims <- class(sims)
  if (ClassSims == "NULL")  SubIts <- 1:MSEobj@nsim
  if (ClassSims == "integer" | ClassSims == "numeric") {
    # sims <- 1:min(MSEobj@nsim, max(sims))
    SubIts <- as.integer(sims)
  }
  if (ClassSims == "logical")  SubIts <- which(sims)
  nsim <- length(SubIts)
  
  ClassYrs <- class(years)
  AllNYears <- MSEobj@proyears
  if (ClassYrs == "NULL") 
    Years <- 1:AllNYears
  if (ClassYrs == "integer" | ClassYrs == "numeric") 
    Years <- years
  if (max(Years) > AllNYears) 
    stop("years exceeds number of years in MSE")
  if (min(Years) <= 0) 
    stop("years must be positive")
  if (min(Years) != 1) {
    message("Not starting from first year. Are you sure you want to do this?")
    message("Probably a bad idea!")
  }
  if (!all(diff(Years) == 1)) 
    stop("years are not consecutive")
  if (length(Years) <= 1) 
    stop("You are going to want more than 1 projection year")
  MSEobj@proyears <- max(Years)
  
  SubF <- MSEobj@F_FMSY[SubIts, SubMPs, Years, drop = FALSE]
  SubB <- MSEobj@B_BMSY[SubIts, SubMPs, Years, drop = FALSE]
  SubC <- MSEobj@C[SubIts, SubMPs, Years, drop = FALSE]
  SubBa <- MSEobj@B[SubIts, SubMPs, Years, drop = FALSE]
  SubFMa <- MSEobj@FM[SubIts, SubMPs, Years, drop = FALSE]
  SubTACa <- MSEobj@TAC[SubIts, SubMPs, Years, drop = FALSE]
  
  OutOM <- MSEobj@OM[SubIts, ]
  # check if slot exists
  tt <- try(slot(MSEobj, "Effort"), silent = TRUE)
  if (class(tt) == "try-error")  slot(MSEobj, "Effort") <- array(NA)
  if (all(is.na(tt)) || all(tt == 0)) slot(MSEobj, "Effort") <- array(NA)
  if (all(is.na(MSEobj@Effort))) {
    SubEffort <- array(NA)
  } else {
    SubEffort <- MSEobj@Effort[SubIts, SubMPs, Years, drop = FALSE]
  }
  
  # check if slot exists
  tt <- try(slot(MSEobj, "SSB"), silent = TRUE)
  if (class(tt) == "try-error") slot(MSEobj, "SSB") <- array(NA)
  if (all(is.na(tt)) || all(tt == 0))slot(MSEobj, "SSB") <- array(NA)
  if (all(is.na(MSEobj@SSB))) {
    SubSSB <- array(NA)
  } else {
    SubSSB <- MSEobj@SSB[SubIts, SubMPs, Years, drop = FALSE]
  }
  
  # check if slot exists
  tt <- try(slot(MSEobj, "VB"), silent = TRUE)
  if (class(tt) == "try-error") slot(MSEobj, "VB") <- array(NA)
  if (all(is.na(tt)) || all(tt == 0)) slot(MSEobj, "VB") <- array(NA)
  if (all(is.na(MSEobj@VB))) {
    SubVB <- array(NA)
  } else {
    SubVB <- MSEobj@VB[SubIts, SubMPs, Years, drop = FALSE]
  }
  
  # check if slot exists
  tt <- try(slot(MSEobj, "PAA"), silent = TRUE)
  if (class(tt) == "try-error") slot(MSEobj, "PAA") <- array(NA)
  if (all(is.na(tt)) || all(tt == 0))slot(MSEobj, "PAA") <- array(NA)
  if (all(is.na(MSEobj@PAA))) {
    SubPAA <- array(NA)
  } else {
    SubPAA <- MSEobj@PAA[SubIts, SubMPs, , drop = FALSE]
  }  
  
  # check if slot exists
  tt <- try(slot(MSEobj, "CAL"), silent = TRUE)
  if (class(tt) == "try-error") slot(MSEobj, "CAL") <- array(NA)
  if (all(is.na(tt)) || all(tt == 0)) slot(MSEobj, "CAL") <- array(NA)
  if (all(is.na(MSEobj@CAL))) {
    SubCAL <- array(NA)
  } else {
    SubCAL <- MSEobj@CAL[SubIts, SubMPs, , drop = FALSE]
  } 
  
  # check if slot exists
  tt <- try(slot(MSEobj, "CAA"), silent = TRUE)
  if (class(tt) == "try-error") slot(MSEobj, "CAA") <- array(NA)
  if (all(is.na(tt)) || all(tt == 0)) slot(MSEobj, "CAA") <- array(NA)
  if (all(is.na(MSEobj@CAA))) {
    SubCAA <- array(NA)
  } else {
    SubCAA <- MSEobj@CAA[SubIts, SubMPs, , drop = FALSE]
  } 
  
  CALbins <- MSEobj@CALbins 
  
  SubResults <- new("MSE", Name = MSEobj@Name, nyears = MSEobj@nyears, 
                    proyears = MSEobj@proyears, nMPs = length(SubMPs), MPs = newMPs, 
                    nsim = length(SubIts), OM = OutOM, Obs = MSEobj@Obs[SubIts, , drop = FALSE],
                    B_BMSY = SubB, F_FMSY = SubF, B = SubBa, SSB=SubSSB, VB=SubVB, 
                    FM = SubFMa,  SubC, 
                    TAC = SubTACa, SSB_hist = MSEobj@SSB_hist[SubIts, , , , drop = FALSE], 
                    CB_hist = MSEobj@CB_hist[SubIts, , , , drop = FALSE], 
                    FM_hist = MSEobj@FM_hist[SubIts, , , , drop = FALSE], 
                    Effort = SubEffort, PAA=SubPAA, CAL=SubCAL, CAA=SubCAA , CALbins=CALbins,
                    Misc=list())
  
  return(SubResults)
}

#' @describeIn checkMSE Joins two or more MSE objects together. MSE objects must have identical
#' number of historical years, and projection years. Also works for Hist objects returned
#' by `runMSE(Hist=TRUE)`
#' @export
joinMSE <- function(MSEobjs = NULL) {
  # join two or more MSE objects
  if (class(MSEobjs) != "list") stop("MSEobjs must be a list")
  if (length(MSEobjs) < 2) stop("MSEobjs list doesn't contain multiple MSE objects")
  
  lapply(MSEobjs, checkMSE) # check that MSE objects contains all slots 
  
  ishist <- all(lapply(MSEobjs, slotNames) %>% unlist() %>% unique() %in% slotNames('Hist'))
  
  if (ishist) {
    out <- new("Hist")
    sls <- slotNames('Hist')
    nsim <- MSEobjs[[1]]@Ref$B0 %>% length()
    for (sl in sls) {
      obj <-lapply(MSEobjs, slot, name=sl)
      if (sl == "Data") {
        out@Data <- joinData(obj)
      } else {
        if (class(obj[[1]]) == "data.frame") {
          slot(out, sl) <- do.call('rbind', obj)
        }
        if (class(obj[[1]]) == "list") {
          out.list <- list()
          for (nm in names(obj[[1]])) {
            obj2 <- lapply(obj, '[[', nm)
            ind <- which(dim(obj2[[1]]) == nsim)
            if (length(ind) >0) {
              if (class(obj2[[1]]) == "array") {
                tempVal <- lapply(obj2, dim)
                # check all dimensions the same (hack for different CAL bins)
                tdf <- lapply(obj2, dim) %>% unlist() %>% 
                  matrix(nrow=length(obj2), ncol=length(dim(obj2[[1]])),byrow=TRUE)
                nBins <- tdf[,2]
                Max <- max(nBins)
                nyrs <- max(tdf[,3])
                nsims <- sapply(tempVal, function(x) x[1])
                if (!mean(nBins) == max(nBins)) { # not all same size
                  index <- which(nBins < Max)
                  for (kk in index) {
                    dif <- Max - dim(obj2[[kk]])[2]
                    obj2[[kk]] <- abind::abind(obj2[[kk]], array(0, dim=c(nsims[kk], dif, nyrs)), along=2)
                  }
                }
              }
              out.list[[nm]] <- abind::abind(obj2, along=ind)  
            } else {
              out.list[[nm]] <- unlist(obj2) %>% unique()
            }
          }
          slot(out, sl) <- out.list
        }
      }
    }
    return(out)
  }  
  
  MPNames <- lapply(MSEobjs, getElement, name = "MPs")  # MPs in each object 
  allsame <- length(unique(lapply(MPNames, unique))) == 1
  
  if (!allsame) {
    # drop the MPs that don't appear in all MSEobjs
    mpnames <- unlist(MPNames)
    npack <- length(MSEobjs)
    tab <- table(mpnames)
    ind <- tab == npack
    commonMPs <- names(tab)[ind]
    if (length(commonMPs)<1) stop("No common MPs in MSE objects", call.=FALSE)
    MSEobjs <- lapply(MSEobjs, Sub, MPs = commonMPs)
    message("MPs not in all MSE objects:")
    message(paste(names(tab)[!ind], ""))
    message("Dropped from final MSE object.")
  }
  
  Nobjs <- length(MSEobjs)
  for (X in 1:Nobjs) {
    tt <- MSEobjs[[X]]
    assign(paste0("obj", X), tt)
    if (X > 1) {
      tt <- MSEobjs[[X]]
      tt2 <- MSEobjs[[X - 1]]
      if (!all(slotNames(tt) == slotNames(tt2))) 
        stop("The MSE objects don't have the same slots")
      if (any(tt@MPs != tt2@MPs)) 
        stop("MPs must be the same for all MSE objects")
    }
  }
  
  # Check that nyears and proyears are the same for all
  chkmat <- matrix(NA, nrow = Nobjs, ncol = 2)
  nms <- NULL
  for (X in 1:Nobjs) {
    tt <- get(paste0("obj", X))
    chkmat[X, ] <- c(tt@nyears, tt@proyears)
    if (X > 1) 
      if (!any(grepl(tt@Name, nms))) 
        stop("MSE objects have different names")
    nms <- append(nms, tt@Name)
  }
  chk <- all(colSums(chkmat) == chkmat[1, ] * Nobjs)
  if (!chk) stop("The MSE objects have different number of nyears or proyears")
  
  # Join them together
  Allobjs <- mget(paste0("obj", 1:Nobjs))
  sns <- slotNames(Allobjs[[1]])
  sns<-sns[sns!="Misc"] # ignore the Misc slot
  outlist <- vector("list", length(sns))
  for (sn in 1:length(sns)) {
    templs <- lapply(Allobjs, slot, name = sns[sn])
    if (class(templs[[1]]) == "character") {
      outlist[[sn]] <- templs[[1]]
    }
    if (class(templs[[1]]) == "numeric" | class(templs[[1]]) == "integer") {
      if (sns[sn] == "CALbins") {
        tempInd <- which.max(unlist(lapply(templs, length)))
        CALbins <- templs[[tempInd]]
      } else {
        outlist[[sn]] <- do.call(c, templs)
      }
    }
    if (class(templs[[1]]) == "matrix" | class(templs[[1]]) == "data.frame") {
      outlist[[sn]] <- do.call(rbind, templs)
    }
    if (class(templs[[1]]) == "array") {
      if (sns[sn] == "CAL") { # hack for different sized CAL arrays 
        tempVal <- lapply(templs, dim)
        if (all(unlist(lapply(tempVal, length)) == 3)) {
          nBins <- sapply(tempVal, function(x) x[3])
          nsims <- sapply(tempVal, function(x) x[1])
          nMPs <- sapply(tempVal, function(x) x[2])
          if (!mean(nBins) == max(nBins)) { # not all same size 
            Max <- max(nBins)
            index <- which(nBins < Max)
            for (kk in index) {
              dif <- Max - dim(templs[[kk]])[3]
              templs[[kk]] <- abind::abind(templs[[kk]], array(0, dim=c(nsims[kk], nMPs[kk], dif)), along=3)
            } 
          }      
          outlist[[sn]] <- abind::abind(templs, along = 1)
        } else {
          outlist[[sn]] <- templs[[1]]
        }
      } else {
        outlist[[sn]] <- abind::abind(templs, along = 1)
      }
      
    }
  }
  
  names(outlist) <- sns
  
  Misc<-list()
  if (length(MSEobjs[[1]]@Misc)>0) {
    if (!is.null(MSEobjs[[1]]@Misc$Data)) {
      Misc$Data <- list()
      # Posterior predicted data joining
      for(i in 1:length(MSEobjs[[1]]@Misc$Data)) Misc$Data[[i]]<-joinData(lapply(MSEobjs,function(x)slot(x,"Misc")$Data[[i]]))
    }
    
    if (!is.null(MSEobjs[[1]]@Misc$RInd.stats)) {
      # Error from real indices
      nms <- unique(MSEobjs[[1]]@Misc$RInd.stats$Index) %>% as.character()
      temp <- list()
      for (nm in seq_along(nms)) {
        temp1 <- list()
        for(i in 1:length(MSEobjs)) {
          temp1[[i]] <- MSEobjs[[i]]@Misc$RInd.stats %>% dplyr::filter(Index==nms[nm])
        }
        temp[[nm]] <- do.call('rbind', temp1)
      }
      Misc$RInd.stats <- do.call('rbind', temp)
    }
  }
  
  newMSE <- new("MSE", Name = outlist$Name, nyears = unique(outlist$nyears), 
                proyears = unique(outlist$proyears), nMP = unique(outlist$nMP), 
                MPs = unique(outlist$MPs), nsim = sum(outlist$nsim), OM = outlist$OM, 
                Obs = outlist$Obs, B_BMSY = outlist$B_BMSY, F_FMSY = outlist$F_FMSY, 
                outlist$B, outlist$SSB, outlist$VB,
                outlist$FM, outlist$C, outlist$TAC, outlist$SSB_hist, 
                outlist$CB_hist, outlist$FM_hist, outlist$Effort, outlist$PAA,
                outlist$CAA, outlist$CAL, CALbins, Misc=Misc)
  
  newMSE
}

# Evaluate Peformance of MPs
# --------------------------------------------------- Function examines
# how consistently an MP outperforms another.


#' How dominant is an MP?
#' 
#' The DOM function examines how consistently an MP outperforms another. For
#' example DCAC might provide higher yield than AvC on average but outperforms
#' AvC in less than half of simulations.
#' 
#' 
#' @param MSEobj An object of class 'MSE'
#' @param MPtg A character vector of management procedures for cross
#' examination
#' @return A matrix of performance comparisons length(MPtg) rows by MSE@nMPs
#' columns
#' @author A. Hordyk
#' @export DOM
DOM <- function(MSEobj, MPtg = NA) {
  if (any(is.na(MPtg))) 
    MPtg <- MSEobj@MPs
  proyears <- MSEobj@proyears
  nMP <- MSEobj@nMPs
  nsim <- MSEobj@nsim
  ind <- which(MSEobj@MPs %in% MPtg)
  MPr <- which(!(MSEobj@MPs %in% MPtg))
  yind <- max(MSEobj@proyears - 4, 1):MSEobj@proyears
  y1 <- 1:(MSEobj@proyears - 1)
  y2 <- 2:MSEobj@proyears
  Mat <- matrix(0, nrow = length(MPtg), ncol = nMP)
  rownames(Mat) <- MPtg
  colnames(Mat) <- MSEobj@MPs
  POF <- P100 <- YieldMat <- IAVmat <- Mat
  for (X in 1:length(MPtg)) {
    # Overfishing (F > FMSY)
    ind1 <- as.matrix(expand.grid(1:nsim, ind[X], 1:proyears))
    ind2 <- as.matrix(expand.grid(1:nsim, 1:nMP, 1:proyears))
    t1 <- apply(array(MSEobj@F_FMSY[ind1] > 1, dim = c(nsim, 1, proyears)), 
                c(1, 2), sum, na.rm = TRUE)
    t2 <- apply(array(MSEobj@F_FMSY[ind2] > 1, dim = c(nsim, nMP, proyears)), 
                c(1, 2), sum, na.rm = TRUE)
    POF[X, ] <- round(apply(matrix(rep(t1, nMP), nrow = nsim) < t2, 
                            2, sum)/nsim * 100, 0)
    # B < BMSY
    t1 <- apply(array(MSEobj@B_BMSY[ind1] < 1, dim = c(nsim, 1, proyears)), 
                c(1, 2), sum, na.rm = TRUE)
    t2 <- apply(array(MSEobj@B_BMSY[ind2] < 1, dim = c(nsim, nMP, proyears)), 
                c(1, 2), sum, na.rm = TRUE)
    P100[X, ] <- round(apply(matrix(rep(t1, nMP), nrow = nsim) < t2, 
                             2, sum, na.rm = TRUE)/nsim * 100, 0)
    # Relative yield in last 5 years
    ind1 <- as.matrix(expand.grid(1:nsim, ind[X], yind))
    ind2 <- as.matrix(expand.grid(1:nsim, 1:nMP, yind))
    t1 <- apply(array(MSEobj@C[ind1], dim = c(nsim, 1, length(yind))), 
                c(1, 2), sum, na.rm = TRUE)
    t2 <- apply(array(MSEobj@C[ind2], dim = c(nsim, nMP, length(yind))), 
                c(1, 2), sum, na.rm = TRUE)
    YieldMat[X, ] <- round(apply(matrix(rep(t1, nMP), nrow = nsim) > 
                                   t2, 2, sum, na.rm = TRUE)/nsim * 100, 0)
    # interannual variation in catch
    ind1 <- as.matrix(expand.grid(1:nsim, ind[X], y1))
    ind2 <- as.matrix(expand.grid(1:nsim, ind[X], y2))
    AAVY1 <- apply(array(((MSEobj@C[ind1] - MSEobj@C[ind2])^2)^0.5, 
                         dim = c(nsim, 1, length(y1))), 1, mean, na.rm = T)/apply(array(MSEobj@C[ind2], 
                                                                                        dim = c(nsim, 1, length(y1))), 1, mean, na.rm = T)
    ind1 <- as.matrix(expand.grid(1:nsim, 1:nMP, y1))
    ind2 <- as.matrix(expand.grid(1:nsim, 1:nMP, y2))
    AAVY2 <- apply(array(((MSEobj@C[ind1] - MSEobj@C[ind2])^2)^0.5, 
                         dim = c(nsim, nMP, length(y1))), c(1, 2), mean, na.rm = T)/apply(array(MSEobj@C[ind2], 
                                                                                                dim = c(nsim, nMP, length(y1))), c(1, 2), mean, na.rm = T)
    IAVmat[X, ] <- round(apply(matrix(rep(AAVY1, nMP), nrow = nsim) < 
                                 AAVY2, 2, sum, na.rm = TRUE)/nsim * 100, 0)
  }
  out <- list()
  out$POF <- POF
  out$P100 <- P100
  out$Yd <- YieldMat
  out$AAVY <- IAVmat
  return(out)
}






