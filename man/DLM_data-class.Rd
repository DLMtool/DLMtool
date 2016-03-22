\name{DLM_data-class}
\Rdversion{1.1}
\docType{class}
\alias{DLM_data-class}

\title{Class \code{"DLM_data"}}
\description{
An object for storing data for analysis using data-limited methods
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DLM_data", stock)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{The name of the case-study}
    \item{\code{Year}:}{A vector of years that correspond to catch and relative abundance data}
    \item{\code{Cat}:}{Total annual catches}
    \item{\code{Ind}:}{Relative abundance index}
    \item{\code{t}:}{The number of years corresponding to AvC and Dt}
    \item{\code{AvC}:}{Average catch over time t}
    \item{\code{Dt}:}{Depletion over time t e.g. Bnow/Bthen}
    \item{\code{ML}:}{Mean length time series}
    \item{\code{Mort}:}{Natural mortality rate}
    \item{\code{FMSY_M}:}{An assumed ratio of FMSY to M}
    \item{\code{BMSY_B0}:}{The most productive stock size relative to unfished}
    \item{\code{L50}:}{Length at 50 percent maturity}
    \item{\code{L95}:}{Length at 95 percent maturity}
    \item{\code{Lbar}:}{Mean length of catches over Lc (modal length)}
    \item{\code{Lc}:}{Modal length}
    \item{\code{LFC}:}{Length at first capture}
    \item{\code{LFS}:}{smallest Length at full selection}
    \item{\code{CAA}:}{Catch at Age data}
    \item{\code{Dep}:}{Stock depletion Bnow/Bunfished (total stock)}
    \item{\code{Abun}:}{An estimate of absolute current vulnerable abundance}
    \item{\code{vbK}:}{The von Bertalanffy growth coefficient}
    \item{\code{vbLinf}:}{Maximum length}
    \item{\code{vbt0}:}{Theoretical age at length zero}
    \item{\code{wla}:}{Weight-Length parameter alpha}
    \item{\code{wlb}:}{Weight-Length parameter beta}
    \item{\code{steep}:}{Steepness of the Beverton Holt stock-recruitment relationship}
    \item{\code{CV_Cat}:}{Coefficient of variation in annual catches}
    \item{\code{CV_Dt}:}{Coefficient of variation in depletion over time t}
    \item{\code{CV_AvC}:}{Coefficient of variation in average catches over time t}
    \item{\code{CV_Ind}:}{Coefficient of variation in the relative abundance index}
    \item{\code{CV_Mort}:}{Coefficient of variation in natural mortality rate}
    \item{\code{CV_FMSY_M}:}{Coefficient of variation in the ratio in FMSY/M}
    \item{\code{CV_BMSY_B0}:}{Coefficient of variation in the position of the most productive stock size relative to unfished}
    \item{\code{CV_Dep}:}{Coefficient of variation in current stock depletion}
    \item{\code{CV_Abun}:}{Coefficient of variation in estimate of absolute current stock size}
    \item{\code{CV_vbK}:}{Coefficient of variation in the von Bert. k parameter}
    \item{\code{CV_vbLinf}:}{Coefficient of variation in maximum length}
    \item{\code{CV_vbt0}:}{Coefficient of variation in age at length zero}
    \item{\code{CV_L50}:}{Coefficient of variation in length at 50 per cent maturity}
    \item{\code{CV_LFC}:}{Coefficient of variation in length at first capture}
    \item{\code{CV_LFS}:}{Coefficient of variation in length at full selection}
    \item{\code{CV_wla}:}{Coefficient of variation in weight-length parameter a}
    \item{\code{CV_wlb}:}{Coefficient of variation in weight-length parameter b}
    \item{\code{CV_steep}:}{Coefficient of variation in steepness  }
    \item{\code{sigmaL}:}{Assumed observaton error of the length composition data}
    \item{\code{MaxAge}:}{Maximum age}
    \item{\code{Units}:}{Units of the catch/absolute abundance estimates}
    \item{\code{Ref}:}{A reference quota level}
    \item{\code{Ref_type}:}{Its type}
    \item{\code{Log}:}{A log of events}
    \item{\code{params}:}{A place to store estimated parameters}
    \item{\code{PosMPs}:}{The methods that can be applied to these data}
    \item{\code{MPs}:}{The methods that were applied to these data}
    \item{\code{OM}:}{A table of operating model conditions}
    \item{\code{Obs}:}{A table of observation model conditions}
    \item{\code{TAC}:}{The calculated TAC}
    \item{\code{TACbias}:}{The known bias in the calculated TAC}
    \item{\code{Sense}:}{The results of the sensitivity analysis}
    \item{\code{CAL_bins}:}{The length bins for the catch-at-length data}
    \item{\code{CAL}:}{Catch-at-length data}
    \item{\code{Cref}:}{Reference or target catch level}
    \item{\code{Iref}:}{Reference or target relative abundance index level}
    \item{\code{Bref}:}{Reference or target biomass level}
    \item{\code{CV_Cref}:}{CV for reference or target catch level}
    \item{\code{CV_Iref}:}{CV for reference or target relative abundance index level}
    \item{\code{CV_Bref}:}{CV for reference or target biomass level}
    \item{\code{CV_Rec}:}{CV for recent recruitment strength}
    \item{\code{Rec}:}{Recent recruitment strength}
    \item{\code{MPrec}:}{The previous recommendation of a management proceedure}
    \item{\code{LHYear}:}{The last historical year of the simulation (before projection)}
    \item{\code{Misc}:}{Optional list which is passed to MPs}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "DLM_data")}: ... }
    \item{plot}{\code{signature(x = "DLM_data")}: ... }
    \item{summary}{\code{signature(object = "DLM_data")}: ... }
	 }
}

\author{
T. Carruthers
}
\examples{
newdata<-new('DLM_data')
}
\keyword{classes}
