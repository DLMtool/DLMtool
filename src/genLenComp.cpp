#include <Rcpp.h>
using namespace Rcpp;

//' Generate length composition of catch 
//'
//' Generate size composition of catch given sample of catch-at-age,
//' expected length-at-age, and standard deviation of length-at-age.
//' Model assumes length-at-age is normally distributed, and that
//' selectivity is size-dependant
//'
//' @param CAL_bins vector of catch-at-length size bins
//' @param CAL_binsmid vector (nbins = length(CAL_bins) - 1) of mid-points for catch-at-length size bins 
//' @param SL matrix (nbins, nyears) of selectivity-at-length class for each year
//' @param CAL_ESS effective sample size of catch-at-length data 
//' @param CAL_nsamp sample size of catch-at-length data 
//' @param CN matrix (nyears, maxage) of catch-at-age for each year
//' @param LaA matrix (maxage, nyears) of expected length-at-age for each year
//' @param LaASD matrix (maxage, nyears) of standard deviation of length-at-age for each year
//' @param truncSD optional argument to truncate the length-at-age distribution at `truncSD` standard deviations 
//' e.g., a value of 2 truncates the length-at-age distribution at two standard deviations (set to 0 to ignore (default))
//' 
//' @export
// [[Rcpp::export]]

NumericMatrix  genLenComp(NumericVector CAL_bins, NumericVector CAL_binsmid, NumericMatrix SL,
  double CAL_ESS, double CAL_nsamp, NumericMatrix CN, NumericMatrix LaA, NumericMatrix LaASD,
  double truncSD) {
  
  int nbins = CAL_binsmid.size();
  int nyears = CN.nrow();
  int k = CN.ncol();
  
  NumericVector probs(k);
  NumericMatrix CAL(nyears, nbins);
  IntegerVector ans(k);
  NumericVector ans2(k);
  NumericVector Sx(k);
  NumericMatrix Cx(k, nbins);
  NumericMatrix Cx2(k, nbins);
  
  NumericMatrix Prob(k, nbins); // prob dist length at age 
  
  for (int yr=0; yr < nyears; yr++) {		
	// sample ages from catch 
	NumericVector tempCN(k);
	probs = CN.row(yr)/sum(CN.row(yr));
  rmultinom(CAL_ESS, probs.begin(), k, ans.begin()); // multinom age sample with ess 
	ans2 = ans; // convert to numeric vector
  tempCN = ans2 * (CAL_nsamp/CAL_ESS); // scale up from ess to sample size
  ;   
	NumericVector EL(k); 
	NumericVector SDL(k); 
	
	// probablility of length-at-age - assume normally distributed
	EL = LaA.column(yr);  // expected length-at-age 
	SDL = LaASD.column(yr); // sd for length-at-age
	for (int age=0; age < k; age++) {
		Prob(age, 0) = R::pnorm((CAL_bins(1) - EL(age))/SDL(age), 0, 1, 1, 0); 
		for (int L=1; L < (nbins-1); L++) {
			Prob(age, L) = R::pnorm((CAL_bins(L+1) - EL(age))/SDL(age), 0, 1, 1, 0) - 
			               R::pnorm((CAL_bins(L) - EL(age))/SDL(age), 0, 1, 1, 0);
		}
		Prob(age, nbins-1) = 1 -  R::pnorm((CAL_bins(nbins-1) - EL(age))/SDL(age), 0, 1, 1,0); 
		
		if (truncSD > 0 ) {
		  // Truncate at truncSD standard deviations
		  for (int L=0; L < nbins; L++) {
		  	double temp = 0;
		  	temp = (CAL_binsmid(L) - EL[age])/SDL(age);
		  	if (temp >= truncSD) Prob(age,L) = 0;
		  	if (temp <= -truncSD) Prob(age,L) = 0;
		  }	
		}  
	  Sx(age) = sum(SL.column(yr) * Prob.row(age));	// Calculate selectivity at age given length-at-age dist 	
	  Cx.row(age) = Prob.row(age) * SL.column(yr);  // Length-at-age of catch conditional on selectivity
	  Cx2.row(age) = Cx.row(age)/sum(Cx.row(age));  // standardise so prob sum to 1 for each age 
	  }
	
    for (int L=0; L < (nbins); L++) {	
	    for (int age=0; age < k; age++) {
	      double tempVal = 0;
	      tempVal = (tempCN(age) * Cx2(age, L));
		    if (tempVal > 0) CAL(yr, L) += tempVal;
	    }
	  }
	
  }
  
  return(CAL);
}

