#include <Rcpp.h>
using namespace Rcpp;


//' Internal estimation function for LSRA and LSRA2 functions
//'
//' Rcpp version of R code 
//' @param param a numeric value representing log(R0)
//' @param FF_a numeric value, recent fishign mortality rate (apical F)
//' @param Chist a vector of historical catch observations [nyears]
//' @param M_a numeric value, natural mortality rate
//' @param Mat_age_a a vector of maturity at age [nage]
//' @param Wt_age_a a vector of weight at age [nage]
//' @param sel_a a vector of selectivity at age [nage]
//' @param Recdevs_a a vector of recruitment deviations [nyears]
//' @param h_a a numeric value of steepness values of the Bev-Holt Stock-Recruitment relationship
//' @param Umax maximum harvest rate per year
//' @author T. Carruthers with an amateur attempt at converting to Rcpp by A. Hordyk (but it works!)
//' @useDynLib DLMtool
//' @keywords internal
//' @export
// [[Rcpp::export]]
List  LSRA_opt_cpp(double param,  
                          double FF_a, NumericVector Chist, double M_a, 
                          NumericVector Mat_age_a, NumericVector Wt_age_a, 
                          NumericVector sel_a, NumericVector Recdevs_a, double h_a,
                          double Umax) {
  
  double nyears = Chist.length();
  double maxage = Mat_age_a.length();
  double R0 = exp(param);
  NumericMatrix Nstr(nyears, maxage);
  NumericVector N(maxage);
  NumericVector N2(maxage);
  double SSB0 = 0;
  for (int a = 0; a < maxage; a++) {
    N(a)=R0*exp(-M_a*(a));
    SSB0+=N(a)*Mat_age_a(a)*Wt_age_a(a);
  } 
  double SSBpR = SSB0/R0;
  double pen=0;
  double mupredF=0; 
  
  NumericVector PredF(nyears);
  NumericVector SSB(nyears);
  NumericVector PredN(maxage);
  NumericVector PredVN(maxage);
  NumericVector PredVW(maxage);
  NumericVector Predfrac(maxage);
  NumericVector Cat(maxage);
  NumericVector predU(maxage);
  NumericVector cond(maxage);
  
  NumericVector  pentemp = 0;
  NumericVector Cattemp = 0;
  NumericVector Cattemp2 = 0;
  
  Nstr(0,_) = N;

  for (int y=0; y<nyears; y++) {
    N = Nstr(y,_);
    SSB(y) = sum(N*Mat_age_a*Wt_age_a);
    PredN = N*exp(-M_a/2);
    PredVN = PredN*sel_a;
    PredVW = PredVN*Wt_age_a;
    Predfrac = PredVW/sum(PredVW);
    Cat = Chist(y)*Predfrac;
    
    predU = Cat/(PredN*Wt_age_a);
    cond = predU>Umax;
    if(sum(cond)>0){
      pentemp =  predU[cond>0];
      pen+=sum(pow(abs(pentemp -Umax),2));
      Cattemp = Cat[cond>0];
      Cattemp2 = Cattemp/(pentemp/Umax);
      Cat[cond>0] = Cattemp2;
    }
    PredF(y) = -log(1-max(Cat/(N*Wt_age_a)));
    
    N = N *exp(-M_a-PredF(y)*sel_a);
    for (int a = 1; a < maxage; a++)  N2(a) = N(a-1);
    N2(0) =  Recdevs_a(y)*(0.8*R0*h_a*SSB(y))/(0.2*SSBpR*R0*(1-h_a)+(h_a-0.2)*SSB(y));
    if (y<(nyears-1)) Nstr(y+1,_) = N2;
    
  }

  NumericVector ind(11);
  for (int x=0; x<11; x++)  ind(x) = nyears-(16-x);
  NumericVector Val =  PredF[ind];
  mupredF = sum(Val)/11;

  List out(5);
  out[0]= pen+pow((log(mupredF)-log(FF_a)),2);
  out[1]= PredF;
  out[2]= SSB/SSB0;
  out[3]= param;
  out[4]= SSB;
  return(out);
  
}

//' Internal SRA MCMC CPP code
//'
//' Rcpp version of R code 
//' @param nits number of iterations
//' @param pars vector of parameters
//' @param JumpCV jump cv vector
//' @param adapt adapt vector
//' @param parLB lower bounds
//' @param parUB upper bounds
//' @param R0ind index for R0
//' @param inflind index for inflection
//' @param slpind index for slope
//' @param RDind index for recruitment deviations
//' @param nyears number of projection years
//' @param maxage maximum age
//' @param M Natural mortality
//' @param Mat_age A vector of maturity at age 
//' @param Wt_age A vector of weight at age 
//' @param Chist_a A vector of historical catch observations (nyears long) going back to unfished conditions
//' @param Umax A numeric value representing the maximum harvest rate for any age class (rejection of sims where this occurs)
//' @param h steepness of SRR
//' @param CAA A matrix nyears (rows) by nages (columns) of catch at age (age 1 to maxage in length)
//' @param CAAadj internal parameter
//' @param sigmaR A numeric value representing the prior standard deviation of log space recruitment deviations
//' 
//' @author A. Hordyk
//' @export
// [[Rcpp::export]]
List  LSRA_MCMC_sim(double nits, NumericVector pars,
                NumericVector JumpCV, NumericVector adapt, NumericVector parLB,
                NumericVector parUB, int R0ind, int inflind,
                int slpind, IntegerVector RDind, int nyears, int maxage,
                double M, NumericVector Mat_age, NumericVector Wt_age, 
                NumericVector Chist_a, double Umax, double h, NumericMatrix CAA,
                double CAAadj, double sigmaR) {

  double npars = pars.length();
  NumericVector acceptpars(npars);
  NumericVector pars2(clone(pars));
  NumericVector LHstr(nits);
  NumericVector SSB(nyears);
  NumericVector PredF(nyears);
  NumericVector RD(nyears);
  NumericVector sel(maxage);
  
  NumericMatrix parstr(npars, nits);
  NumericMatrix CAA_pred(nyears, maxage);
  
  double SSB0 = 0; 
  
  for (int i=0; i<nits; i++) {
    
    // bool temp = is_true(any(update == i));
    // if (temp == TRUE)   Rcout << ".";
    NumericVector Reject(1);
    NumericVector Accept(1);
    
    // Sample new parameters 
    NumericVector nupars_t(npars);
    NumericVector nupars(npars);
   
    for (int x=0; x<npars; x++) {
      if (i == 0) {
        nupars(x) = exp(pars2(x));
      } else {
        nupars_t(x) = R::rnorm(acceptpars(x), JumpCV(x)*adapt(i));
        if (nupars_t(x) < parLB(x)) nupars_t(x) = parLB(x);
        if (nupars_t(x) > parUB(x)) nupars_t(x) = parUB(x);
        nupars(x) = exp(nupars_t(x));
      }
      
    }
    
    double R0 = nupars[R0ind];
    double infl = nupars[inflind];
    double slp = nupars[slpind];

    NumericVector recdevs = nupars[RDind];
    RD = recdevs/mean(recdevs);
    
    NumericVector initN(maxage);
    
    for (int a=0; a<maxage; a++) {
      double age = a+1;
      sel(a) = 1/(1+exp((infl-age)/slp));
      initN(a) = R0 * exp(-M*a); // equilibrium numbers at age 
    }
  
    NumericVector N(maxage);
    NumericVector PredN(maxage);
    NumericVector PredVN(maxage);
    NumericVector PredVW(maxage);
    NumericVector Predfrac(maxage);
    NumericVector Cat(maxage);
    NumericVector predU(maxage);
    LogicalVector cond(maxage);
    
    
 
    
    SSB0 = sum(initN*Mat_age*Wt_age);
    double SSBpR = SSB0/R0;
    
    for (int y=0; y<nyears; y++) {
      if(y==0) {
        Range idx(0, maxage-1);
        NumericVector recs2 = RD[rev(idx)];
        N = initN * recs2;
      }
   
      SSB(y) = sum(N*Mat_age*Wt_age);
      PredN = N*exp(-M);
      
      PredVN = PredN * sel;
      CAA_pred(y,_) = PredVN/sum(PredVN);
      PredVW = PredVN * Wt_age;
      Predfrac = PredVW/sum(PredVW);
      Cat = (Chist_a(y) * Predfrac)/Wt_age;
      predU = Cat/PredN;
      cond = predU>Umax; // Check max U
      if (is_true(any(cond > 0))) {

        Reject(0) = 1; // Reject where U > Umax for any age class
        NumericVector tempCat =  Cat[cond>0];
        NumericVector temppredU = predU[cond>0];
        NumericVector tempCat2 = (tempCat)/(temppredU/Umax);
        Cat[cond>0] = tempCat2;
      } else {
        Reject(0) = 0;
      } 
      PredF(y) = -log(1-max(Cat/PredN));
      N = PredN - Cat;
    
      // Aging  
      for (int a = (maxage-1); a >0; a--)  N(a) = N(a-1); // aging
      N(0) =  RD(y+maxage)*(0.8*R0*h*SSB(y))/(0.2*SSBpR*R0*(1-h)+(h-0.2)*SSB(y));
     
    }   
    // Rcout << "here" << "\n\n";
    double CAALH = 0;
    for (int y=0; y<nyears; y++) {
      // Rcout << "y = " << y << "\n\n";
      NumericVector tempVec = CAA_pred(y,_);
      CAALH += sum(log(tempVec)*(CAA(y,_)/CAAadj));
    }

    // Rcout << "here" << "\n\n";
    
    double RDLH = sum(dnorm(log(recdevs),-(pow(sigmaR,2))/2, sigmaR, true));
    // Rcout << RDLH << "\n\n";
    double LH = CAALH + RDLH;
  
    NumericVector rand = runif(1);
    if (i > 0) {
      Accept(0) = rand(0) < exp(LH - LHstr(i-1));
      if (Reject(0) > 0) Accept(0) = 0;
      parstr(_,i) = parstr(_,i-1);
      LHstr(i) = LHstr(i-1);
      if (Accept(0) > 0) {
        acceptpars = log(nupars);
        parstr(_,i) = log(nupars);
        LHstr(i) = LH;
      }
    } else {
      parstr(_,i) = pars2;
      acceptpars = pars2;
      LHstr(i) = LH;
      
    }
  
  }
    
  
  List out(7);
  out[0] = parstr;
  out[1] = CAA_pred;
  out[2] = SSB;
  out[3] = SSB0;
  out[4] = RD;
  out[5] = PredF;
  out[6] = sel;
  return(out);
}

