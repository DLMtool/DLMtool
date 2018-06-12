#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List LBSPRgen(double SL50, double SL95, double FM, int nage, int nlen, double CVLinf, 
                       NumericVector LenBins, NumericVector LenMids, NumericVector x, double MK, double Linf, double P,
                       double L50, double L95, double Beta) {
  
  NumericVector rLens(nage);
  NumericMatrix Prob (nage, nlen);
  NumericMatrix Cx(nage, nlen);
  for (int a=0; a<nage; a++) {
    rLens(a) = (1-pow(P, x[a]/MK));
  }
  NumericVector EL = rLens * Linf;
  NumericVector SDL = EL * CVLinf;
  
  for (int l=0; l<nlen; l++) {
    NumericVector temp1 = (LenBins(l+1) - EL)/SDL;
    NumericVector temp2 = (LenBins(l) - EL)/SDL;
    if(l==0) {
      Prob(_, l) = pnorm(temp1, 0.0, 1.0);
    } else {
      Prob(_, l) = pnorm(temp1, 0.0, 1.0) - pnorm(temp2, 0.0, 1.0);
    }
  }
  NumericVector temp2 = (LenBins(nlen-1) - EL)/SDL;
  Prob(_, nlen-1) = 1 - pnorm(temp2, 0.0, 1.0);
  
  NumericVector SL = 1/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50)));
  NumericVector Sx(nage);
  NumericVector MSX(nage);
  NumericVector Ns(nage);
  NumericVector Ml = 1/(1+exp(-log(19)* (LenMids-L50)/(L95-L50)));
  NumericVector Ma(nage);
  for (int a=0; a<nage; a++) {
    Sx(a) = sum(SL * Prob(a,_));
    MSX(a) = sum(Sx)/(a+1);
    Cx(a,_) = Prob(a,_) * SL;
    Ns(a) = pow((1-rLens(a)), (MK+(MK*FM)*MSX(a)));\
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
                NumericVector LenBins, NumericVector LenMids, NumericVector x, 
                double MK, double Linf, double P, double L50, double L95, double Beta) {
  
  double SL50 = exp(pars(0)) * Linf;
  double dSL50 = exp(pars(1));
  double SL95 = SL50 + dSL50 * SL50; 
  double FM = exp(pars(2));
  
  NumericVector predLen = LBSPRgen(SL50, SL95, FM, nage, nlen, CVLinf, 
                                   LenBins, LenMids, x, MK, Linf, P, L50, L95, Beta)(0);
  NumericVector CAL_st = CAL/sum(CAL);
  
  double NLL = 0;
  LogicalVector above0 = CAL_st > 0;
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




