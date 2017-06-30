#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//[[Rcpp::export]]

Rcpp::NumericVector popdynOneTScpp(double nareas, double maxage, Rcpp::NumericVector SSBcurr,
                                Rcpp::NumericMatrix Ncurr,  Rcpp::NumericMatrix Zcurr, double PerrYr,
                                double hs,  Rcpp::NumericVector R0a,  Rcpp::NumericVector SSBpR,
                                Rcpp::NumericVector aR,  Rcpp::NumericVector bR,  Rcpp::NumericMatrix mov,
                                double SRrel) {
  
  Rcpp::NumericMatrix Nnext(maxage, nareas);
  Rcpp::NumericMatrix tempMat2(nareas, nareas);	
  Rcpp::NumericMatrix tempMat3(nareas, nareas);
  Rcpp::NumericMatrix Nstore(maxage, nareas); 
  
  // Recruitment assuming regional R0 and stock wide steepness
  for (int A=0; A < nareas; A++) {
    if (SRrel == 1) {
      // BH SRR
      Nnext(0, A) = PerrYr * (4*R0a(A) * hs * SSBcurr(A))/(SSBpR(A) * R0a(A) * (1-hs) + (5*hs-1) * SSBcurr(A));
    }	
    if (SRrel == 2) {
      // most transparent form of the Ricker uses alpha and beta params
      Nnext(0, A) = PerrYr * aR(A) * SSBcurr(A) * exp(-bR(A) * SSBcurr(A));
    }
    
    // Mortality
    for (int age=1; age<maxage; age++) {
      Nnext(age, A) = Ncurr(age-1, A) * exp(-Zcurr(age-1, A)); // Total mortality
    }
  }
  
  // Move stock
  for (int age=0; age<maxage; age++) {
    for (int AA = 0; AA < nareas; AA++) {
      for (int BB = 0; BB < nareas; BB++) {
        if (AA != BB) tempMat2(AA, BB) = (Nnext(age, AA) * mov(AA, AA)) +   (Nnext(age, BB) * mov(BB, AA));
        tempMat3(AA, BB) = tempMat2(AA, BB);				 
      }
      Nstore(age, AA) = sum(tempMat2.row(AA));
    }
  }   
  
  return Nstore;
} 




  

