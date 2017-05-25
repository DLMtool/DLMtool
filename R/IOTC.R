
#' Indian Ocean Tuna Commission trade-off plot
#'
#' A one-panel trade-off plot showing the probability of exceeding a biomass reference level
#' and a yield reference level
#'
#' @param MSEobj An object of class MSE created by the function runMSE()
#' @param Bref A biomass reference level (an improper fraction of BMSY)
#' @param Yref A yield reference level (an improper fraction of yield given FMSY management)
#' @param Bsat The satisficing level for biomass (required fraction of simulations exceeding Bref)
#' @param Ysat The satisficing level for yield (required fraction of simulations exceeding Yref)
#' @param xlim The limits of the x axis plotting
#' @param ylim The limits of the y axis plotting 
#' @author T. Carruthers
#' @export IOTC_plot
IOTC_plot = function( MSEobj, Bref = 0.75, Yref = 0.75,    # Reference levels for fractions of BMSY and yield at FMSY
                      Bsat = 0.8, Ysat = 0.8,      # Satisficing lines plotted at these levels
                      xlim=c(0,1.1), ylim=c(0,1.1) ){    # The lower bound on the axes
  
  cols = rep(makeTransparent(c("Black","Red","blue","green","orange")),99)
  yend = max(MSEobj@proyears) - (4:0)
  PB = apply(MSEobj@B_BMSY > Bref, 2, mean)
  LTY = apply(MSEobj@C[,,yend] / MSEobj@OM$RefY > Yref,2,mean)
 
   plot(PB, LTY, col = "white", xlim=xlim, ylim=ylim, 
       xlab = paste0("Prob. biomass over ", round(Bref*100,1),"% of BMSY"),
       ylab = paste0("Prob. yield over ", round(Yref*100,1),"% of FMSY"),
       main = paste("MSE tradeoffs for:", MSEobj@Name))
  abline( v = c(0,1), col='#99999950', lwd = 2 )
  abline( h = c(0,1), col='#99999950', lwd = 2 )
  abline( v = Bsat, col = '#99999950', lwd = 2, lty = 2 )
  abline( h = Ysat, col = '#99999950', lwd = 2, lty = 2 )
  text( PB, LTY, MSEobj@MPs, col=cols, font=2)
  data.frame(PB,LTY,row.names=MSEobj@MPs)
  
}


