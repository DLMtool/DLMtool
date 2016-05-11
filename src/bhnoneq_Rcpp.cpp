#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double bhnoneq_LL(NumericVector stpar, NumericVector year, NumericVector Lbar, NumericVector ss,
  double Linf, double K, double Lc, int nbreaks) {
    
  int i;
  int j;
  int m;
  
  int count = year.size();
  int nbr = nbreaks-1;
  
  NumericVector Z(nbreaks+1);
  NumericVector ggyr(nbreaks);
  NumericMatrix dy(nbreaks,count);
  
  NumericMatrix a(nbreaks+1,count);
  NumericMatrix s(nbreaks+1,count);
  NumericMatrix r(nbreaks+1,count);
  NumericVector denom(count);
  NumericVector numsum(count);
  NumericVector num(count);
  NumericVector Lpred(count);

  double sum_square_Lpred = 0.;
  double nyear = 0.;
  double sigma;
  double nLL;
  
  for(i=0;i<=nbr;i++) {
    Z[i] = stpar[i];
    ggyr[i] = stpar[nbr+2+i];
  }
  Z[nbr+1] = stpar[nbr+1];  
    
  for(i=0;i<=nbr;i++) {
    for(m=1;m<=count;m++) {
      if(ggyr[i]>=m) dy(i,m-1) = 0.;
      else dy(i,m-1) = m-ggyr[i];
    }
  }
  if(nbreaks>1) {
    for(i=0;i<nbr;i++){
      for(m=0;m<count;m++) dy(i,m) -= dy(i+1,m);
    }
  }
  
  for(m=0;m<count;m++) {
    denom[m] = 0.;
    numsum[m] = 0.;

    for(i=0;i<=nbr+1;i++) {
      a(i,m) = 1.;
      r(i,m) = 1.;

      if(i<nbr+1) s(i,m) = 1. - exp(-(Z[nbr+1-i]+K) * dy(nbr-i,m));
      if(i==nbr+1) s(i,m) = 1.;

      if(i>0) {
        for(j=0;j<=i-1;j++) {
          a(i,m) *= exp(-Z[nbr+1-j] * dy(nbr-j,m));
          r(i,m) *= exp(-(Z[nbr+1-j] + K) * dy(nbr-j,m));
        }
      }

      if(i<=nbr) denom[m] += a(i,m)*(1. - exp(-Z[nbr+1-i] * dy(nbr-i,m)))/Z[nbr+1-i];
      if(i==nbr+1) denom[m] += a(i,m)/Z[nbr+1-i];

      numsum[m] += r(i,m) * s(i,m) / (Z[nbr+1-i] + K);
    }
    
    num[m] = Linf * (denom[m] - (1. - Lc/Linf) * numsum[m]);
    Lpred[m] = num[m]/denom[m];
    
    if(ss[m]>0) {
      sum_square_Lpred += ss[m] * pow(Lbar[m]-Lpred[m],2.);
      nyear += 1.;
    }
  }

  sigma = sqrt(sum_square_Lpred/nyear);

  nLL = -nyear * log(sigma) - 0.5 * sum_square_Lpred/(sigma*sigma);
  nLL *= -1;
  
  return nLL;
}

