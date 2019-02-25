#include <Rcpp.h>
using namespace Rcpp;

//' Internal estimation function for LBSPR MP
//'
//' @param SL50 Length at 50 percent selectivity
//' @param SL95 Length at 95 percent selectivity
//' @param FM Ratio of apical fishing mortality to natural mortality
//' @param nage Number of pseudo age-classes
//' @param nlen Number of length bins
//' @param CVLinf CV of length-at-age
//' @param LenBins Vector of length bins
//' @param LenMids Vector of mid-points of length bins
//' @param MK Ratio of M/K
//' @param Linf Asymptotic length
//' @param rLens Vector of relative length at ate
//' @param Prob ALK 
//' @param Ml Maturity at age vector
//' @param L50 Length at 50 percent maturity
//' @param L95 Length at 95 percent maturity
//' @param Beta Exponent of the length-weight relationship
//' @author A. Hordyk
//' @useDynLib DLMtool
//' @keywords internal
//' @export
// [[Rcpp::export]]
List LBSPRgen(double SL50, double SL95, double FM, int nage, int nlen, double CVLinf, 
                       NumericVector LenBins, NumericVector LenMids, double MK, double Linf,  
                       NumericVector rLens, NumericMatrix Prob, NumericVector Ml,
                       double L50, double L95, double Beta) {
  
  NumericMatrix Cx(nage, nlen);
  NumericVector SL = 1/(1+exp(-log(19.0)*(LenMids-SL50)/(SL95-SL50)));
  NumericVector Sx(nage);
  NumericVector MSX(nage);
  NumericVector Ns(nage);
  NumericVector Ma(nage);
  for (int a=0; a<nage; a++) {
    Sx(a) = sum(SL * Prob(a,_));
    MSX(a) = sum(Sx)/(a+1);
    Cx(a,_) = Prob(a,_) * SL;
    Ns(a) = pow((1-rLens(a)), (MK+(MK*FM)*MSX(a)));
    Ma(a) = sum(Prob(a,_) * Ml);
  }
  
  NumericVector Nc(nlen);
  for (int l=0; l<nlen; l++) {
    Nc(l) = sum(Ns * Cx(_,l));
  }
  Nc = Nc/sum(Nc);
  
  // SPR 
  NumericVector N0 = pow((1-rLens), MK); // unfished numbers-at-age
  double FishedEgg = sum(Ma * Ns * pow(rLens, Beta));
  double UnfishedEgg = sum(Ma * N0 * pow(rLens, Beta));
  
  double SPR = FishedEgg/UnfishedEgg;
  
  List out(2);
  out(0) = Nc;
  out(1) = SPR;
  return(out);
}


// [[Rcpp::export]]
double LBSPRopt(NumericVector pars, NumericVector CAL, int nage, int nlen, double CVLinf, 
                NumericVector LenBins, NumericVector LenMids,  
                double MK, double Linf, NumericVector rLens, NumericMatrix Prob, 
                NumericVector Ml, double L50, double L95, double Beta) {
  
  double SL50 = exp(pars(0)) * Linf;
  double dSL50 = exp(pars(1));
  double SL95 = SL50 + dSL50 * SL50; 
  double FM = exp(pars(2));
  
  NumericVector predLen = LBSPRgen(SL50, SL95, FM, nage, nlen, CVLinf, 
                                   LenBins, LenMids, MK, Linf, rLens, Prob, Ml,
                                   L50, L95, Beta)(0);
  NumericVector CAL_st = CAL/sum(CAL);
  
  
  
  double NLL = 0;
  LogicalVector above0 = (CAL_st > 0) & (predLen > 0);
  for (int l=0; l<nlen; l++) {
    if (above0(l)) {
      NLL +=  (CAL(l) * log(predLen(l)/CAL_st(l)));

    }
  }
  

  // add penalty for selectivity
  double Pen=0;
  double PenVal=NLL;
  Pen=R::dbeta(exp(pars(0)), 5.0, 0.1,0) * PenVal;
  if (exp(pars(0)) >= 1) Pen=PenVal*exp(pars(0));
  NLL = NLL + Pen;
  
  return(-NLL);
}




