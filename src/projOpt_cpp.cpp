#include <Rcpp.h>
using namespace Rcpp;

//' Rcpp version of the Projection Optimizer
//'
//' Optimize for MSY and calculate MSY reference points 
//'
//' @param lnIn internal
//' @param Mc internal
//' @param hc internal
//' @param Mac internal
//' @param Wac internal
//' @param R0c internal
//' @param Vc internal
//' @param nyears internal
//' @param maxage internal
//' @param movc internal
//' @param Spat_targc internal
//' @param SRrelc internal
//' @param aRc internal
//' @param bRc internal
//' @param movc internal
//' @param SSBpRc internal
//' @param proyears internal
//' @param Control internal
//' 
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector projOpt_cpp(double lnIn,  
                          NumericVector Mc, double hc, NumericVector Mac, NumericVector Wac, 
                          double R0c, NumericVector Vc, double nyears, double maxage, NumericMatrix movc,
                          double Spat_targc, double SRrelc, NumericVector aRc, NumericVector bRc, 
                          double proyears, double Control) {
  
  double B0 = 0; 
  double nareas = movc.nrow();
  double FMSYc = 0 ;
  double NYears = 0;
  // double RetVal = 0;
  NumericVector idist(nareas);
  NumericVector store(nareas);
  NumericMatrix idistA(nareas, nareas);
  NumericMatrix N(maxage, nareas); 
  NumericMatrix CB(maxage, nareas); 
  NumericMatrix CN(maxage, nareas); 
  NumericMatrix Nstore(maxage, nareas); 
  NumericMatrix SSN(maxage, nareas); 
  NumericMatrix Biomass(maxage, nareas); 
  NumericMatrix VB(maxage, nareas);  
  NumericMatrix SSB(maxage, nareas);  
  NumericMatrix FMc(maxage, nareas); 
  NumericMatrix Zc(maxage, nareas); 
  
  NumericMatrix tempMat(maxage, nareas);
  NumericMatrix tempMat3(nareas, nareas);
  NumericVector R0a(nareas);
  NumericVector SSBpR(nareas);
  NumericVector SSBbyA(nareas);
  NumericVector targ(nareas);
  NumericVector SSB0(nareas); 
  NumericVector tempVec(nareas);

  // Controls for switching 
  NYears = nyears + proyears;
  FMSYc = exp(lnIn); 	 
    
  for (int i = 0; i < nareas; i++) idist[i] = 1/nareas;
  for (int X = 0; X < 300; X++) {
	for (int A= 0; A < nareas; A++) {  
      for (int i = 0; i < nareas; i++) {	  
	    for (int j = 0; j < nareas; j++) { 
         idistA(i,j) = idist(i) * movc(i,j);
        }
	  }
      store(A) = sum(idistA.column(A));	  
	}
	for (int A= 0; A < nareas; A++) idist(A) = store(A);
  }
  
  for (int A=0; A < nareas; A++) {
	  for (int age=0; age < maxage; age++) {
	  	if (age==0) N(age, A) = R0c * idist(A);
	  	if (age > 0) N(age, A) = N(age-1, A) * exp(-Mc(age));
	  	SSN(age, A) = Mac(age) * N(age, A); 
	  	Biomass(age, A) = Wac(age) * N(age, A);
	  	SSB(age, A) = SSN(age, A) * Wac(age);	
	  	B0 = sum(Biomass);
	  	R0a(A) = idist(A) * R0c;
	  	SSB0(A) = sum(SSB.column(A));
	  	SSBpR(A) = SSB0(A) / R0a(A);
	  } 
  }	
	  
  for (int yr=0; yr < (NYears-1); yr++) {
	  for (int A=0; A < nareas; A++) {
		for (int age=0; age < maxage; age++) tempMat(age, A) = Vc(age) * Biomass(age, A);
		tempVec(A) = pow(sum(tempMat.column(A)), Spat_targc); 
	  }
	  		
	  for (int A=0; A < nareas; A++) {
	    targ(A) = tempVec(A) / (sum(tempVec)/nareas); // calculate spatial targetting
		SSBbyA(A) = sum(SSB.column(A));

		for (int age=0; age < maxage; age++) {
      FMc(age, A) = FMSYc * Vc(age) * targ(A); // 
		  Zc(age, A) = FMc(age, A) + Mc(age);  
		  if (age > 0) Nstore(age, A) = N(age-1, A) * exp(-Zc(age-1, A));
	      
		  // Recruitment 
		  // no process error when estimating MSY			   
		  if (age == 0) {
		    if (SRrelc == 1) {
 		      Nstore(0, A) = (0.8 * R0a(A) * hc * SSBbyA(A))/
                          (0.2 * SSBpR(A) * R0a(A) * (1-hc) + (hc - 0.2) * SSBbyA(A));     						   
		    }	
		    if (SRrelc == 2) {
 		      Nstore(0, A) = aRc(A) * SSBbyA(A) * exp(-bRc(A) * SSBbyA(A));
		    }		
		  }	
		}
		  for (int age=0; age < maxage; age++) N(age, A) = Nstore(age, A);
	  }	
      
	  for (int age=0; age < maxage; age++) {
		// Movement  				
      NumericMatrix tempMat2(nareas, nareas);		
		  for (int AA = 0; AA < nareas; AA++) {
		    for (int BB = 0; BB < nareas; BB++) {
			    if (AA != BB) tempMat2(AA, BB) = (N(age, AA) * movc(AA, AA)) + (N(age, BB) * movc(BB, AA));
			    tempMat3(AA, BB) = tempMat2(AA, BB);				 
		    }
		    Nstore(age, AA) = sum(tempMat2.row(AA));
		  }	
        // if (yr == 0) Rcpp::Rcout << "val: " << Nstore << std::endl;
		for (int A =0; A < nareas; A++) {
		  N(age, A) = Nstore(age, A);
		  CN(age, A) = FMc(age, A)/Zc(age,A) * N(age,A) * (1-exp(-Zc(age, A)));
		  CB(age, A) = CN(age, A) * Wac(age);	  
		  SSN(age, A) = N(age, A) * Mac(age);
      SSB(age, A) = SSN(age, A) * Wac(age);
      Biomass(age, A) = N(age, A) * Wac(age);
      VB(age, A) = N(age, A) * Wac(age) * Vc(age);		  
		}	  
	  }
     // if (yr == 0) Rcpp::Rcout << "CB: " << CB << std::endl;	  
	 // Rcpp::Rcout << "N: " << sum(N) << std::endl;	 
  }
  
  if (Control == 1) {
	NumericVector out(1);
    out(0) = -sum(CB);
	return out; 
  } else {
	  NumericVector out(5);
    out(0) = sum(SSB);
    out(1) = sum(Biomass);
    out(2) = sum(VB);
	  out(3) = sum(SSB)/sum(SSB0);
	  out(4) = sum(Biomass)/B0;
	return out; 
  } 
}






