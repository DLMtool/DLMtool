#' Data-Limited Methods Toolkit
#' 
#' Simulation testing and implementation of data-limited fishery stock
#' assessment methods
#' 
#' @name DLMtool-package
#' @aliases DLMtool-package DLMtool
#' @docType package
#' @author Tom Carruthers \email{t.carruthers@@fisheries.ubc.ca} 
#' @author Adrian Hordyk \email{a.hordyk@oceans.ubc.ca}
#' @section Additional Information:
#' See the \href{https://dlmtool.github.io/DLMtool/userguide/index.html}{DLMtool User Guide} for
#' a detailed description of how to use the DLMtool package. 
#' 
#' The help documentation for the DLMtool package can also be accessed \href{https://dlmtool.github.io/DLMtool/}{here}.
#' 
#' See the \href{http://www.datalimitedtoolkit.org/}{Data-Limited Toolkit Website} for more information on the DLMtool,
#' including an interactive demo of the main features of the toolkit, information on case studies where the toolkit has 
#' been applied, and more about the history and development of the DLMtool.
#' @references Carruthers, T.R., Punt, A.E., Walters, C.J., MacCall, A.,
#' McAllister, M.K., Dick, E.J., Cope, J. 2014. Evaluating methods for setting
#' catch limits in data-limited fisheries. Fisheries Research. 153: 48-68.
#' 
#' Carruthers, T.R., Kell, L.T., Butterworth, D.S., Maunder, M.N., Geromont,
#' H.F., Walters, C., McAllister, M.K., Hillary, R., Levontin, P., Kitakado,
#' T., Davies, C.R. Performance review of simple management procedures. ICES
#' Journal of Marine Science.
#' @keywords data-limited fisheries stock assessment data-poor management
#' strategy evaluation
#' 
NULL
#' @examples
#' 
#' \dontrun{
#' # --- Application to real fishery data ---
#' 
#' library(DLMtool)
#' setup()                        # setup parallel processing
#' mydata<-new('Data')            # create a new DLM data object and define:
#' mydata@Year<-2001:2010         # years
#' mydata@Cat<-matrix((11:20)*10*runif(10,0.5,1.5),nrow=1) # make up some annual catches
#' mydata@Ind<-matrix(seq(1.1,0.9,length.out=10)*runif(10,0.5,1.5),nrow=1) 
#' mydata@Mort<-0.2               # instantaneous natural mortality rate
#' mydata@Abun<-1000              # current abundance estimate (biomass)
#' mydata@FMSY_M<-0.5             # guess of the ratio of FMSY to natural mortality rate
#' mydata@vbLinf<-200             # maximum length
#' mydata@vbK<-0.2                # von B growth coefficient k
#' mydata@LFC<-50                 # length at first capture
#' mydata<-TAC(mydata)            # calculate quotas
#' plot(mydata)                   # plot them
#' mydata<-Sense(mydata,'Fratio') # conduct a sensitivity analysis for one of the methods
#' }