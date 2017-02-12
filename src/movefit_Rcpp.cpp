#include <Rcpp.h>
using namespace Rcpp;

//' Rcpp version of the Optimization function that returns the squared difference between user
//' specified and calculated movement parameters. 
//'
//' The user specifies the probability of staying in the same area and spatial
//' heterogeneity (both in the unfished state). This function returns the
//' squared difference between these values and those produced by the three
//' logit movement model.
//'
//' This is paired with getmov to find the correct movement model. 
//' 
//' @param par Three parameters in the logit space that control the four
//' probabilities of moving between 2 areas
//' @param prb User specified probability that individuals in area 1 remain in
//' that area (unfished conditions)
//' @param frac User specified fraction of individuals found in area 1 (unfished
//' conditions)
//' 
//' @author T. Carruthers with an amateur attempt at converting to Rcpp by A. Hordyk (but it works!)
//' @useDynLib DLMtool
//' @export
// [[Rcpp::export]]
double  movfit_Rcpp(NumericVector par, double prb, double frac) {
  double nrow = 2;
  double ncol = 2; 
  double NLL = 0;
  NumericMatrix temp(2,2);  
  NumericMatrix mov(2,2);
  NumericMatrix mov2(2,2);
  NumericMatrix mov3(2,2);
  NumericVector dis(2);
  dis(0) = frac; 
  dis(1) = 1 - frac;
  
  mov(0,0) = exp(par(0));
  mov(0,1) = exp(0.0);
  mov(1,0) = exp(par(1));
  mov(1,1) = exp(par(2));
  for (int i = 0; i < nrow; i++) {
	double total = 0;
	for (int j = 0; j < ncol; j++) {
  	  total +=mov(i, j);	
	}	 
    mov2(i, 0) = mov(i, 0)/total;
	mov2(i, 1) = mov(i, 1)/total;
  }
  
   for (int X = 0; X < 100; X++) {
    for (int i = 0; i < nrow; i++) {	  
	  for (int j = 0; j < ncol; j++) {
  	    temp(i, j) = dis(i) * mov2(i, j); 	
	  }		  
	}
 	dis(0) = temp(0, 0) + temp(1, 0); 
	dis(1) = temp(0, 1) + temp(1, 1);
   }  
  NLL = pow(log(mov2(0,0)) - log(prb),2) +
    pow(log(frac) - log(dis(0)),2);
 return NLL;
}



