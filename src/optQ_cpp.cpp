#include <Rcpp.h>
using namespace Rcpp;


//' Rcpp version of the q Optimizer
//'
//' Optimize for catchability coefficient
//'
//' @param lnIn internal
//' @param Fc internal
//' @param Perrc internal
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
//' 
//' @export
//' @keywords internal
// [[Rcpp::export]]
double optQ_cpp(double lnIn, double depc, NumericVector Fc, 
  NumericVector Perrc, NumericMatrix Mc, double hc, NumericVector Mac, NumericMatrix Wac, 
  double R0c, NumericMatrix Vc, double nyears, double maxage, NumericMatrix movc, 
  double Spat_targc, double SRrelc, NumericVector aRc, NumericVector bRc) {
  
  double nareas = movc.nrow();
  double qc = 0;
  double RetVal = 0;

  NumericVector idist(nareas);
  NumericVector store(nareas);
  NumericMatrix idistA(nareas, nareas);
  NumericMatrix Neq(maxage, nareas);   
  NumericMatrix N(maxage, nareas); 
  NumericMatrix CB(maxage, nareas); 
  NumericMatrix CN(maxage, nareas); 
  NumericMatrix Nstore(maxage, nareas); 
  NumericMatrix SSNeq(maxage, nareas);   
  NumericMatrix SSN(maxage, nareas); 
  NumericMatrix Biomass(maxage, nareas); 
  NumericMatrix VB(maxage, nareas);  
  NumericMatrix SSBeq(maxage, nareas);  
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

  qc = exp(lnIn);  
    
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
  
  // Equilibrium initial conditions
  for (int A=0; A < nareas; A++) {
	  for (int age=0; age < maxage; age++) {
	  	if (age == 0) Neq(age, A) = R0c  * idist(A);
	  	if (age > 0) Neq(age, A) = Neq(age-1, A) * exp(-Mc(age-1, 0));
	  	SSNeq(age, A) = Mac(age) * Neq(age, A); 
	  	SSBeq(age, A) = SSNeq(age, A) * Wac(age, 0);
	  }
		// B0 = sum(Biomass);
		R0a(A) = idist(A) * R0c;
		SSB0(A) = sum(SSBeq.column(A));
		SSBpR(A) = SSB0(A) / R0a(A);
  }	
  
  // Rcout << Neq << "\n";
  // Rcout << sum(SSB0) << "\n";
  NumericVector tempVec2(maxage);


  // Non-equilibrium initial conditions
  for (int A=0; A < nareas; A++) {
	  for (int age=0; age < maxage; age++) {
	    if (age == 0) N(age, A) = R0c * idist(A) * Perrc(maxage-age-1);
	    tempVec2 = Mc(_,0);
	   
	    if (age > 0) {
	      NumericVector tempVec3(tempVec2.begin(), tempVec2.begin()+age) ;
	      N(age, A) = N(age, A) = R0c * idist(A) * Perrc(maxage-age-1) * exp(-sum(tempVec3));
	    }
	    
	    // if (age > 0) N(age, A) = R0c * idist(A) * Perrc(maxage-age-1) * exp(-(age-1, 0));
	  	SSN(age, A) = Mac(age) * N(age, A); 
	  	Biomass(age, A) = Wac(age, 0) * N(age, A);
	  	SSB(age, A) = SSN(age, A) * Wac(age, 0);	
	  } 
  }	
  // Rcout << tempVec2 << "\n\n";
  // Rf_PrintValue(N) ;

  // Calculate initial Fs 
  for (int A=0; A < nareas; A++) {
	  for (int age=0; age < maxage; age++) tempMat(age, A) = Vc(age, 0) * Biomass(age, A);
	  tempVec(A) = pow(sum(tempMat.column(A)), Spat_targc); 
  }
  for (int A=0; A < nareas; A++) {
    targ(A) = tempVec(A) / (sum(tempVec)/nareas);
    for (int age=0; age < maxage; age++) {
      FMc(age, A) = qc * Fc(0) * Vc(age, 0) * targ(A);
      Zc(age, A) = FMc(age, A) + Mc(age, 0);
    }	
  }
  for (int yr=0; yr < (nyears-1); yr++) {
	  for (int A=0; A < nareas; A++) {
		for (int age=0; age < maxage; age++) tempMat(age, A) = Vc(age, yr) * Biomass(age, A);
		tempVec(A) = pow(sum(tempMat.column(A)), Spat_targc); 
	  }
	  for (int A=0; A < nareas; A++) {
	    targ(A) = tempVec(A) / (sum(tempVec)/nareas);
		
		SSBbyA(A) = sum(SSB.column(A));
		for (int age=0; age < maxage; age++) {
		  if (age > 0) Nstore(age, A) = N(age-1, A) * exp(-Zc(age-1, A)); 	      
		  // Recruitment 	  			   
		  if (age == 0) {
		    if (SRrelc == 1) {
 		      Nstore(0, A) = Perrc(yr+maxage-1) * (0.8 * R0a(A) * hc * SSBbyA(A))/
                          (0.2 * SSBpR(A) * R0a(A) * (1-hc) + (hc - 0.2) * SSBbyA(A));     						   
		    }	
		    if (SRrelc == 2) {
 		      Nstore(0, A) = Perrc(yr+maxage-1) * aRc(A) * SSBbyA(A) * exp(-bRc(A) * SSBbyA(A));
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
			    if (AA != BB) tempMat2(AA, BB) = (N(age, AA) * movc(AA, AA)) +
			                (N(age, BB) * movc(BB, AA));
			    tempMat3(AA, BB) = tempMat2(AA, BB);				 
		    }
		    Nstore(age, AA) = sum(tempMat2.row(AA));
		  }	
        
		for (int A =0; A < nareas; A++) {
		  N(age, A) = Nstore(age, A);
		  CN(age, A) = FMc(age, A)/Zc(age,A) * N(age,A) * (1-exp(-Zc(age, A)));
		  CB(age, A) = CN(age, A) * Wac(age, yr);
		  SSN(age, A) = N(age, A) * Mac(age);
      SSB(age, A) = SSN(age, A) * Wac(age, yr);
      Biomass(age, A) = N(age, A) * Wac(age, yr);
      VB(age, A) = N(age, A) * Wac(age, yr) * Vc(age, yr);
	    FMc(age, A) = qc * Fc(yr+1) * Vc(age, yr) * targ(A); // 
		  Zc(age, A) = FMc(age, A) + Mc(age, yr); 		  
		}	  
	  }
   // if (yr==1) Rcout << "N = " << N << "\n";
   // Rcout << sum(N) << "\n";	 
  // Rcout << sum(SSB)/sum(SSB0) << "\n";   
  }
  // 
  // Rcout << "SSB = " << sum(SSB) << "\n";
  // Rcout << "SSB0 = " << sum(SSB0) << "\n";
  // Rcout << "SB/SB0 = " << sum(SSB)/sum(SSB0) << "\n";
  // Rcout << "dep = " << depc << "\n";

  
  RetVal = pow(log(depc) - log(sum(SSB)/sum(SSB0)),2); // optimize q
  return RetVal; 
}
