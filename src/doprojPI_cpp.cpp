#include <Rcpp.h>
using namespace Rcpp;

//' Rcpp version of the Projection function for calculating Reference Yield
//'
//' 
//' @param lnF internal
//' @param Mmat internal
//' @param Wac internal
//' @param Mac internal
//' @param Pc internal
//' @param N_c internal
//' @param SSN_c internal
//' @param Biomass_c internal
//' @param VBiomass_c internal
//' @param SSB_c internal
//' @param Vc internal
//' @param hc internal
//' @param R0ac internal
//' @param proyears internal
//' @param nareas internal
//' @param maxage internal
//' @param movc internal
//' @param SSBpRc internal
//' @param aRc internal
//' @param bRc internal
//' @param SRrelc internal
//' @param Spat_targc internal
//' 
//' @export
//' @keywords internal
// [[Rcpp::export]]
double doprojPI_cpp(double lnF, NumericMatrix Mmat,  
  NumericVector Wac, NumericVector Mac, NumericVector Pc, NumericMatrix N_c,
  NumericMatrix SSN_c, NumericMatrix Biomass_c, NumericMatrix VBiomass_c, 
  NumericMatrix SSB_c, NumericMatrix Vc, double hc, NumericVector R0ac, double proyears,
  double nareas, double maxage, NumericMatrix movc, double SSBpRc,   
  NumericVector aRc, NumericVector bRc, double SRrelc, double Spat_targc)  {
  
  double FF = exp(lnF);
  
  NumericMatrix N_Pcurr(clone(N_c));
  NumericMatrix N_Pnext(maxage, nareas);  
  NumericMatrix SSN_P(clone(SSN_c));  
  NumericMatrix Biomass_P(clone(Biomass_c));
  NumericMatrix VBiomass_P(clone(VBiomass_c)); 

  NumericMatrix SSB_P(clone(SSB_c)); 
  NumericMatrix Biomass_Pnext(maxage, nareas); 
  NumericMatrix VBiomass_Pnext(maxage, nareas);   
  NumericMatrix SSB_Pnext(maxage, nareas);    
  NumericMatrix FM_P(maxage, nareas); 
  NumericMatrix Z_P(maxage, nareas);  
  NumericMatrix Cyr(maxage, nareas);  
  NumericMatrix Nstore(maxage, nareas); 
  NumericMatrix tempMat(maxage, nareas); 
  NumericMatrix tempMat3(maxage, nareas); 
  NumericVector tempVec(nareas);
  NumericVector C_P(proyears);
  NumericVector fishdist(nareas);
  NumericVector SSBbyA(nareas);
  
  // movement   
  for (int A=0; A < nareas; A++) {
	for (int age=0; age < maxage; age++) tempMat(age, A) = VBiomass_P(age, A);
	tempVec(A) = pow(sum(tempMat.column(A)), Spat_targc); 
  }
  for (int A=0; A < nareas; A++) {
	  fishdist(A) = tempVec(A) / (sum(tempVec)/nareas);
	  for (int age=0; age<maxage; age++) {
		  FM_P(age, A) = FF * Vc(age,0) * fishdist(A);
		  Z_P(age, A) = FM_P(age, A) + Mmat(age,0);  
      }
  }   
  
  for (int yr=1; yr < proyears; yr++) {
	  for (int A=0; A < nareas; A++) {
          SSBbyA(A) = sum(SSB_P.column(A));		  
	  	  for (int age=1; age<maxage; age++) {
			  N_Pnext(age, A) = N_Pcurr(age-1,A) * exp(-Z_P(age-1, A));
		  }
  	      if (SRrelc == 1) {
 		      N_Pnext(0, A) = Pc(yr)* (0.8 * R0ac(A) * hc * SSBbyA(A))/
			    (0.2 * SSBpRc * R0ac(A) * (1-hc) + (hc - 0.2) * SSBbyA(A));     						   
		  }	
		  if (SRrelc == 2) {
 		      N_Pnext(0, A) = Pc(yr) * aRc(A) * SSBbyA(A) * exp(-bRc(A) * SSBbyA(A));
		  }
	  } 
	  
      // move fish 
	  for (int age=0; age < maxage; age++) {				
          NumericMatrix tempMat2(nareas, nareas);		
		  for (int AA = 0; AA < nareas; AA++) {
		    for (int BB = 0; BB < nareas; BB++) {
			  if (AA != BB) tempMat2(AA, BB) = (N_Pnext(age, AA)* movc(AA, AA)) +
			                (N_Pnext(age, BB) * movc(BB, AA));				 
		    }
		    Nstore(age, AA) = sum(tempMat2.row(AA));
			tempMat3 = tempMat2;
		  }
          for (int A=0; A < nareas; A++) {
		    N_Pnext(age, A) = Nstore(age, A);
		    Biomass_P(age, A) = N_Pnext(age, A) * Wac(age,yr); 
            VBiomass_P(age, A) = Biomass_P(age, A) * Vc(age,yr);		
		    SSN_P(age, A) = N_Pnext(age, A) * Mac(age);
		    SSB_P(age, A) = SSN_P(age, A) * Wac(age,yr);	    
		}		  
	  }
	   // if (yr==1) Rcpp::Rcout << Biomass_P << std::endl;
	   
	  // move fishing effort 
      for (int A=0; A < nareas; A++) {
	    for (int age=0; age < maxage; age++) tempMat(age, A) = VBiomass_P(age, A);
	    tempVec(A) = pow(sum(tempMat.column(A)), Spat_targc); 
      }	 
     	  
      for (int A=0; A < nareas; A++) {
	    fishdist(A) = tempVec(A) / (sum(tempVec)/nareas);
	    // if (yr==1) Rcpp::Rcout << "The value of fishdist is " << fishdist << std::endl;
		for (int age=0; age<maxage; age++) {
		  FM_P(age, A) = FF * Vc(age,yr) * fishdist(A);		  
		  Z_P(age, A) = FM_P(age, A) + Mmat(age,yr);  
		  Cyr(age, A) = FM_P(age,A)/Z_P(age,A) * (1-exp(-Z_P(age, A))) * Biomass_P(age,A);
		  // Cyr(age, A) = FM_P(age,A)/Z_P(age,A) * (1-exp(-Z_P(age, A))) * VBiomass_P(age,A);
		  N_Pcurr(age, A) = N_Pnext(age, A);
        }
      }
       // if (yr==1) Rcpp::Rcout << fishdist << std::endl;	  
	  
    C_P(yr) = sum(Cyr); 
  }
  
  double RetVal = 0; 
  double temp = 0;
  double temp2 = 0;
  temp = std::min(4.0, (proyears-1));
  temp2 = proyears - temp - 1; 
  NumericVector Store(temp+1);
  // Rcpp::Rcout << Store <<  std::endl;
  int XX = 0;
  for (int i=temp2; i<proyears; i++) {
	  Store(XX) = C_P(i);
      XX = XX + 1;	  
  }
  RetVal = sum(Store)/(temp+1);
 
  
  return -RetVal;
}

