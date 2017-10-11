#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Population dynamics model for one annual time-step
//'
//' Project population forward one time-step given current numbers-at-age and total mortality
//'
//' @param nareas The number of spatial areas
//' @param maxage The maximum age 
//' @param SSBcurr A numeric vector of length nareas with the current spawning biomass in each area
//' @param Ncurr A numeric matrix (maxage, nareas) with current numbers-at-age in each area
//' @param Zcurr A numeric matrix (maxage, nareas) with total mortality-at-age in each area
//' @param PerrYr A numeric value with recruitment deviation for current year 
//' @param hs Steepness of SRR
//' @param R0a Numeric vector with unfished recruitment by area
//' @param SSBpR Numeric vector with unfished spawning stock per recruit by area 
//' @param aR Numeric vector with Ricker SRR a parameter by area
//' @param bR Numeric vector with Ricker SRR b parameter by area
//' @param mov Numeric matrix (nareas by nareas) with the movement matrix
//' @param SRrel Integer indicating the stock-recruitment relationship to use (1 for Beverton-Holt, 2 for Ricker)
//' 
//' @author A. Hordyk
//' 
//' @export
//' @keywords internal
//[[Rcpp::export]]
Rcpp::NumericMatrix popdynOneTScpp(double nareas, double maxage, Rcpp::NumericVector SSBcurr,
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



//' Population dynamics model in CPP
//'
//' Project population forward pyears given current numbers-at-age and total mortality, etc 
//' for the future years
//'
//' @param nareas The number of spatial areas
//' @param maxage The maximum age 
//' @param SSBcurr A numeric vector of length nareas with the current spawning biomass in each area
//' @param Ncurr A numeric matrix (maxage, nareas) with current numbers-at-age in each area
//' @param pyears The number of years to project the population forward
//' @param M_age Numeric matrix (maxage, pyears) with natural mortality by age and year
//' @param Asize_c Numeric vector (length nareas) with size of each area
//' @param MatAge Numeric vector with proportion mature by age
//' @param WtAge Numeric matrix (maxage, pyears) with weight by age and year
//' @param Vuln Numeric matrix (maxage, pyears) with vulnerability by age and year
//' @param Retc Numeric matrix (maxage, pyears) with retention by age and year
//' @param Prec Numeric vector (pyears) with recruitment error
//' @param mov Numeric matrix (nareas by nareas) with the movement matrix
//' @param SRrelc Integer indicating the stock-recruitment relationship to use (1 for Beverton-Holt, 2 for Ricker)
//' @param Effind Numeric vector (length pyears) with the fishing effort by year
//' @param Spat_targc Integer. Spatial targetting
//' @param hc Numeric. Steepness of stock-recruit relationship
//' @param R0c Numeric vector of length nareas with unfished recruitment by area
//' @param SSBpRc Numeric vector of length nareas with unfished spawning per recruit by area
//' @param aRc Numeric. Ricker SRR a value by area
//' @param bRc Numeric. Ricker SRR b value by area
//' @param Qc Numeric. Catchability coefficient
//' @param Fapic Numeric. Apical F value
//' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
//' @param control Integer. 1 to use q and effort to calculate F, 2 to use Fapic (apical F) and 
//' vulnerablity to calculate F.
//' 
//' @author A. Hordyk
//' 
//' @export
//' @keywords internal
//[[Rcpp::export]]
List popdynCPP(double nareas, double maxage, arma::mat Ncurr, double pyears,
               arma::mat M_age, arma::vec Asize_c, arma::mat MatAge, arma::mat WtAge,
               arma::mat Vuln, arma::mat Retc, arma::vec Prec,
               NumericMatrix movc, double SRrelc, arma::vec Effind,
               double Spat_targc, double hc, NumericVector R0c, NumericVector SSBpRc,
               NumericVector aRc, NumericVector bRc, double Qc, double Fapic, double maxF, int control) {
  
  arma::cube Narray(maxage, pyears, nareas, arma::fill::zeros);
  arma::cube Barray(maxage, pyears, nareas, arma::fill::zeros);
  arma::cube SSNarray(maxage, pyears, nareas, arma::fill::zeros);
  arma::cube SBarray(maxage, pyears, nareas, arma::fill::zeros);
  arma::cube VBarray(maxage, pyears, nareas, arma::fill::zeros);
  arma::cube Marray(maxage, pyears, nareas, arma::fill::zeros);
  arma::cube FMarray(maxage, pyears, nareas, arma::fill::zeros);
  arma::cube FMretarray(maxage, pyears, nareas, arma::fill::zeros);
  arma::cube Zarray(maxage, pyears, nareas, arma::fill::zeros);
  
  NumericVector tempVec(nareas);
  arma::vec fishdist(nareas);
  
  // Initial year
  Narray.subcube(0, 0, 0, maxage-1, 0, nareas-1) = Ncurr;
  for (int A=0; A<nareas; A++) {
    Barray.subcube(0, 0, A, maxage-1, 0, A) = Ncurr.col(A) % WtAge.col(0);
    SSNarray.subcube(0, 0, A, maxage-1, 0, A) = Ncurr.col(A) % MatAge.col(0); 
    SBarray.subcube(0, 0, A, maxage-1, 0, A) = Ncurr.col(A) % WtAge.col(0) % MatAge.col(0);
    VBarray.subcube(0, 0, A, maxage-1, 0, A) = Ncurr.col(A) % WtAge.col(0) % Vuln.col(0);
    Marray.subcube(0, 0, A, maxage-1, 0, A) = M_age.col(0);
    tempVec(A) = accu(VBarray.slice(A));
  }
  
  fishdist = (pow(tempVec, Spat_targc))/mean((pow(tempVec, Spat_targc)));

  // calculate F at age for first year
  if (control == 1) {
    for (int A=0; A<nareas; A++) {
      FMarray.subcube(0,0, A, maxage-1, 0, A) =  (Effind(0) * Qc * fishdist(A) * Vuln.col(0))/Asize_c(A);
      FMretarray.subcube(0,0, A, maxage-1, 0, A) =  (Effind(0) * Qc * fishdist(A) * Retc.col(0))/Asize_c(A);
    }
  }
  if (control == 2) {
    for (int A=0; A<nareas; A++) {
      FMarray.subcube(0,0, A, maxage-1, 0, A) =  (Fapic * fishdist(A) * Vuln.col(0))/Asize_c(A);
      FMretarray.subcube(0,0, A, maxage-1, 0, A) =  (Fapic * fishdist(A) * Retc.col(0))/Asize_c(A);
    }
  }
  
  // apply Fmax condition 
  arma::uvec tempvals = arma::find(FMarray > (1-exp(-maxF)));
  FMarray.elem(tempvals).fill(1-exp(-maxF));
  arma::uvec tempvals2 = arma::find(FMretarray > (1-exp(-maxF)));
  FMretarray.elem(tempvals2).fill(1-exp(-maxF));
  
  Zarray.subcube(0,0, 0, maxage-1, 0, nareas-1) = Marray.subcube(0,0, 0, maxage-1, 0, nareas-1) + FMarray.subcube(0,0, 0, maxage-1, 0, nareas-1);
  
  
  for (int yr=0; yr<(pyears-1); yr++) {
    
    arma::vec SB(nareas);
    
    for (int A=0; A<nareas; A++) SB(A) = accu(SBarray.subcube(0, yr, A, maxage-1, yr, A));
    arma::mat Ncurr2 = Narray.subcube(0, yr, 0, maxage-1, yr, nareas-1);
    arma::mat Zcurr = Zarray.subcube(0, yr, 0, maxage-1, yr, nareas-1);
    NumericMatrix NextYrNa = popdynOneTScpp(nareas, maxage, wrap(SB), wrap(Ncurr2), wrap(Zcurr), 
                                       Prec(yr+1+maxage), hc, R0c, SSBpRc, aRc, bRc, movc, SRrelc); 
    arma::mat NextYrN = as<arma::mat>(NextYrNa);

    Narray.subcube(0, yr+1, 0, maxage-1, yr+1, nareas-1) = NextYrN;
    for (int A=0; A<nareas; A++) {
      Barray.subcube(0, yr+1, A, maxage-1, yr+1, A) = NextYrN.col(A) % WtAge.col(yr+1);
      SSNarray.subcube(0, yr+1, A, maxage-1, yr+1, A) = NextYrN.col(A) % MatAge.col(yr+1);
      SBarray.subcube(0, yr+1, A, maxage-1, yr+1, A) = NextYrN.col(A) % WtAge.col(yr+1) % MatAge.col(yr+1);
      VBarray.subcube(0, yr+1, A, maxage-1, yr+1, A) = NextYrN.col(A) % WtAge.col(yr+1) % Vuln.col(yr+1);
      Marray.subcube(0, yr+1, A, maxage-1, yr+1, A) = M_age.col(yr+1);
      tempVec(A) = accu(VBarray.subcube(0, yr+1, A, maxage-1, yr+1, A));
    }

    fishdist = (pow(tempVec, Spat_targc))/mean((pow(tempVec, Spat_targc)));

    // calculate F at age for next year
    if (control == 1) {
      for (int A=0; A<nareas; A++) {
        FMarray.subcube(0,yr+1, A, maxage-1, yr+1, A) =  (Effind(yr+1) * Qc * fishdist(A) * Vuln.col(yr+1))/Asize_c(A);
        FMretarray.subcube(0,yr+1, A, maxage-1, yr+1, A) =  (Effind(yr+1) * Qc * fishdist(A) * Retc.col(yr+1))/Asize_c(A);
      }
    }
    if (control == 2) {
      for (int A=0; A<nareas; A++) {
        FMarray.subcube(0,yr+1, A, maxage-1, yr+1, A) =  (Fapic * fishdist(A) * Vuln.col(yr+1))/Asize_c(A);
        FMretarray.subcube(0,yr+1, A, maxage-1, yr+1, A) =  (Fapic * fishdist(A) * Retc.col(yr+1))/Asize_c(A);
      }
    }
    // apply Fmax condition 
    arma::uvec tempvals3 = arma::find(FMarray > (1-exp(-maxF)));
    FMarray.elem(tempvals3).fill(1-exp(-maxF));
    arma::uvec tempvals4 = arma::find(FMretarray > (1-exp(-maxF)));
    FMretarray.elem(tempvals4).fill(1-exp(-maxF));
    
    Zarray.subcube(0,yr+1, 0, maxage-1, yr+1, nareas-1) = Marray.subcube(0,yr+1, 0, maxage-1, yr+1, nareas-1) + FMarray.subcube(0,yr+1, 0, maxage-1, yr+1, nareas-1);

  }

  List out(8);
  out(0) = Narray;
  out(1) = Barray;
  out(2) = SSNarray;
  out(3) = SBarray;
  out(4) = VBarray;
  out(5) = FMarray;
  out(6) = FMretarray;
  out(7) = Zarray;
  
  return out;
}

  

