#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



// https://stackoverflow.com/questions/30175104/how-to-effectively-combine-a-list-of-numericvectors-into-one-large-numericvector 
// [[Rcpp::export]]
NumericVector combine(const List& list)
{
  std::size_t n = list.size();
  
  // Figure out the length of the output vector
  std::size_t total_length = 0;
  for (std::size_t i = 0; i < n; ++i)
    total_length += Rf_length(list[i]);
  
  // Allocate the vector
  NumericVector output = no_init(total_length);
  
  // Loop and fill
  std::size_t index = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    NumericVector el = list[i];
    std::copy(el.begin(), el.end(), output.begin() + index);
    
    // Update the index
    index += el.size();
  }
  
  return output;
  
}



// https://stackoverflow.com/questions/13661065/superimpose-histogram-fits-in-one-plot-ggplot

// [[Rcpp::export]]
NumericVector get_freq(NumericVector x, double width, double origin = 0, 
                       int outlen=0) {
  int bin= 0;
  int nmissing = 0;
  std::vector<int> out(outlen);
  
  NumericVector::iterator x_it = x.begin(), x_end;
  for(; x_it != x.end(); ++x_it) {
    double val = *x_it;
    if (ISNAN(val)) {
      ++nmissing;
    } else {
      bin = (val - origin) / width;
      ++out[bin];
    }
  }
  
  return wrap(out);
}


// [[Rcpp::export]]
double which_maxC(NumericVector x){
  int out;
  out = std::distance(x.begin(),std::max_element(x.begin(),x.end()));
  out++;
  
  return out;
}



// https://stackoverflow.com/questions/14034200/efficient-random-number-generation-from-a-truncated-normal-distribution

// [[Rcpp::export]]
NumericVector rnormSelect2(int N, int mi, int ma) {
  RNGScope scope;
  int N2 = N * 1.25;
  NumericVector X = rnorm(N2, 0, 1);
  LogicalVector ind = (X >= mi) & (X <= ma);
  NumericVector Y(N);
  int k=0;
  for (int i=0; (i<N2) & (k<N); i++) {
    if (ind[i]) Y(k++) = X(i);
  }
  return Y;
}

// [[Rcpp::export]]
NumericVector tdnorm(NumericVector x, double mi, double ma) {
  RNGScope scope;
  NumericVector dist = dnorm(x, 0.0, 1.0);
  
  NumericVector cdist = pnorm(x, 0.0, 1.0);
  LogicalVector ind = (cdist < R::pnorm(mi, 0.0, 1.0,1,0)) | (cdist >R::pnorm(ma, 0.0, 1.0,1,0));

  double sz = dist.size();
  int maxind = which_maxC(dist);

  NumericVector Y = dist;
  for (int i=0; (i<sz); i++) { // truncate
    if (ind[i]) Y(i) = 0;
  }
  if (sum(Y)==0) Y[maxind] = 1;
  Y = Y/sum(Y);
  return Y;
}



// [[Rcpp::export]]
NumericMatrix  genSizeComp(NumericMatrix VulnN, NumericVector CAL_binsmid, NumericMatrix selCurve,
                           double CAL_ESS, double CAL_nsamp,
                           NumericVector Linfs, NumericVector Ks, NumericVector t0s,
                           double LenCV, double truncSD) {

  int nyears = VulnN.nrow();
  int k = VulnN.ncol();
  int nbins = CAL_binsmid.size();
  NumericMatrix CAL(nyears, nbins);
  double width = CAL_binsmid(1) - CAL_binsmid(0);
  double origin = CAL_binsmid(0) - 0.5* width;
  NumericVector temp(k);

  NumericVector varAges = NumericVector::create(-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5); // monthly ages

  for (int yr=0; yr < nyears; yr++) {
    NumericVector Nage = (VulnN.row(yr)); // numbers of catch-at-age this year
    double Ncatch = sum(Nage); // total catch this year
    if (Ncatch>0) {
      NumericVector Nage2 = (Nage/Ncatch) * CAL_ESS; // number-at-age effective sample
      List Lens(k*12);
      int count = 0;
      for (int age=1; age <= k; age++) { // loop over 1:maxage
        double Nage3 =  round(Nage2(age-1)); // number at this age
        NumericVector rands = RcppArmadillo::sample(NumericVector::create(0,1,2,3,4,5,6,7,8,9,10,11), Nage3, TRUE, NumericVector::create()) ; //  assume ages are uniformly distributed across months
        NumericVector subAgeVec = get_freq(rands, 1, 0, 12); // distribute n across months
        for (int subage=0; subage<=11; subage++) { // loop over 12 months
          if (subAgeVec(subage) > 0) {
            double sage = varAges(subage) + age;
            double mean = Linfs(yr) * (1-exp(-Ks(yr)* (sage - t0s(yr)))); // calculate mean length at sub-age;
            if (mean < 0) mean = 0.01;
            NumericVector dist = tdnorm((CAL_binsmid-mean)/(LenCV*mean), -truncSD, truncSD); // prob density of lengths for this age
           
            NumericVector newdist = dist * selCurve; // probability = dist * size-selection curve
            newdist = newdist/sum(newdist);
            
            Lens(count) = RcppArmadillo::sample(CAL_binsmid, subAgeVec(subage), TRUE, newdist); // sample lengths for this sub-age class
          } else {
            Lens(count) = NA_INTEGER;
          }
          count += 1;
        }
        // double mean = Linfs(yr) * (1-exp(-Ks(yr)* (age - t0s(yr)))); // calculate mean length at age;
        // if (mean < 0) mean = 0;
        // NumericVector dist = tdnorm((CAL_binsmid-mean)/(LenCV*mean), -truncSD, truncSD); // prob density of lengths for this age
        // NumericVector newdist = dist * selCurve; // probability = dist * size-selection curve
        // newdist = newdist/sum(newdist);
        // double Nage3 =  floor(Nage2(age-1));
        // 
        // if (Nage3 > 0) {
        //   Lens(age-1) = RcppArmadillo::sample(CAL_binsmid, Nage3, TRUE, newdist); // sample lengths for this age class
        // } else {
        //   Lens(age-1) = NA_INTEGER;
        // }
        
      }
      
      NumericVector LenVals = combine(Lens); // unlist
      NumericVector templens = get_freq(LenVals, width, origin, nbins); // calculate frequencies
      double rat = CAL_nsamp/sum(templens);
      templens =  templens * rat; // scale to CAL_nsamp
      CAL(yr,_) = templens;
    } else {
      NumericVector zeros(nbins);
      CAL(yr,_) = zeros;
    }

  }


  return(CAL);
}


// [[Rcpp::export]]
NumericMatrix  genSizeComp2(NumericMatrix VulnN, NumericVector CAL_binsmid, 
                           double CAL_ESS, double CAL_nsamp,
                           NumericVector Linfs, NumericVector Ks, NumericVector t0s,
                           double LenCV, double truncSD) {
  
  int nyears = VulnN.nrow();
  int k = VulnN.ncol();
  int nbins = CAL_binsmid.size();
  NumericMatrix CAL(nyears, nbins);
  double width = CAL_binsmid(1) - CAL_binsmid(0);
  double origin = CAL_binsmid(0) - 0.5* width;
  NumericVector temp(k);
  
  for (int yr=0; yr < nyears; yr++) {
    
    NumericVector probs = VulnN.row(yr)/sum(VulnN.row(yr)); // probability of each age class
    IntegerVector ageSamp1(k);
    rmultinom(CAL_ESS, probs.begin(), k, ageSamp1.begin()); // multinom age sample with ess
    temp = ageSamp1; // convert to numeric vector
    NumericVector ageSamp1a = temp * (CAL_nsamp/CAL_ESS); // scale up to CAL_nsamp
    NumericVector ageSamp1b = round(ageSamp1a,0); // round to integers at age
    List Lens(k);
    
    
    for (int age=1; age <= k; age++) {
      if (ageSamp1b(age-1)>0) {
        int ns = ageSamp1b(age-1);
        NumericVector Age(1);
        Age(0)= age;
        NumericVector ageSamp2 = rep(Age(0), ns) - 0.5 + runif(ageSamp1b(age-1)); // add variability to ages
        NumericVector means = Linfs(yr) * (1-exp(-Ks(yr)* (ageSamp2 - t0s(yr)))); // calculate mean length at age;
        NumericVector rands = rnormSelect2(ageSamp1b(age-1), -truncSD, truncSD); // random normal truncated at truncSD 
        NumericVector tempLens = (rands * means * LenCV) + means; // random sample of length-at-age
        tempLens[tempLens<0] = 0;
        Lens(age-1) = tempLens;
      } else {
        Lens(age-1) = NA_INTEGER; 
      }
      
    }
    
    NumericVector LenVals = combine(Lens); // unlist 
    CAL(yr,_) = get_freq(LenVals, width, origin, nbins); // calculate frequencies
    
  }
  
  return(CAL);
}


// // [[Rcpp::export]]
// NumericVector repcpp(NumericVector x, NumericVector y) {
//   int n = y.size();
//   NumericVector myvector(sum(y));
//   int ind=0;
//   for (int i=0; i < n; ++i) {
//     int p = y[i];
//     std::fill(myvector.begin()+ind, myvector.begin()+ind+p, x[i]);
//     ind += p;
//   }
//   return myvector;
// }

// 
// NumericMatrix  makeLenComp(NumericVector AgeVec, NumericVector SubAgeVec, 
//                            NumericVector Linfarray_c, 
//                            NumericVector Karray_c, NumericVector t0array_c, 
//                            double LenCV_c, NumericVector CAL_bins, 
//                            NumericVector CAL_binsmid, NumericMatrix retLength, 
//                            double CAL_ESS, double CAL_nsamp, 
//                            NumericMatrix VulnN, double truncSD) {
//   
//   int nyears = VulnN.nrow();
//   int k = SubAgeVec.size() - 1;
//   int nbins = CAL_binsmid.size();
//   
//   
//   NumericMatrix CN(nyears,k);
//   
//   NumericVector probs(k);
//   NumericMatrix CAL(nyears, nbins);
//   IntegerVector ans(k);
//   NumericVector ans2(k);
//   NumericVector Sx(k);
//   NumericMatrix Cx(k, nbins);
//   NumericMatrix Cx2(k, nbins);
//   
//   NumericMatrix Prob(k, nbins); // prob dist length at age 
//   double width = SubAgeVec(1) - SubAgeVec(0);
//   double origin = SubAgeVec(0);
//   double outlen = SubAgeVec.size() - 1;
//   for (int yr=0; yr < nyears; yr++) {
//     NumericVector nvec = VulnN(yr,_);
//     NumericVector tempval = repcpp(AgeVec, nvec) - 0.5 + runif(sum(nvec)); // add variability to ages
//     CN(yr,_) =  get_freq(tempval, width, origin, outlen); 
//     
//     // sample ages from catch
//     probs = CN.row(yr)/sum(CN.row(yr));
//     rmultinom(CAL_ESS, probs.begin(), k, ans.begin()); // multinom age sample with ess
//     ans2 = ans; // convert to numeric vector
//     NumericVector tempCN = ans2 * (CAL_nsamp/CAL_ESS); // scale up from ess to sample size
//     
//     // probablility of length-at-age - assume normally distributed
//     
//     
//     for (int age=0; age < k; age++) {
//       
//       double EL = Linfarray_c(yr) * (1-exp(-Karray_c(yr) * (SubAgeVec(age)-t0array_c(yr))));  // expected length-at-age
//       double SDL = EL * LenCV_c; // sd for length-at-age
//       
//       Prob(age, 0) = R::pnorm((CAL_bins(1) - EL)/SDL, 0, 1, 1, 0);
//       for (int L=1; L < (nbins-1); L++) {
//         Prob(age, L) = R::pnorm((CAL_bins(L+1) - EL)/SDL, 0, 1, 1, 0) -
//           R::pnorm((CAL_bins(L) - EL)/SDL, 0, 1, 1, 0);
//       }
//       Prob(age, nbins-1) = 1 -  R::pnorm((CAL_bins(nbins-1) - EL)/SDL, 0, 1, 1,0);
//       
//       if (truncSD > 0 ) {
//         // Truncate at truncSD standard deviations
//         for (int L=0; L < nbins; L++) {
//           double temp = (CAL_binsmid(L) - EL)/SDL;
//           if (temp >= truncSD) Prob(age,L) = 0;
//           if (temp <= -truncSD) Prob(age,L) = 0;
//         }
//       }
//       Sx(age) = sum(retLength.column(yr) * Prob.row(age));	// Calculate selectivity at age given length-at-age dist
//       Cx.row(age) = Prob.row(age) * retLength.column(yr);  // Length-at-age of catch conditional on selectivity
//       Cx2.row(age) = Cx.row(age)/sum(Cx.row(age));  // standardise so prob sum to 1 for each age
//     }
//     
//     for (int L=0; L < (nbins); L++) {
//       for (int age=0; age < k; age++) {
//         double tempVal = 0;
//         tempVal = (tempCN(age) * Cx2(age, L));
//         if (tempVal > 0) CAL(yr, L) += tempVal;
//       }
//     }
//     
//   }
//   
//   return(CAL);
// }
// 
// 
// 
// 
// 
// // NumericVector get_freq(NumericVector x, NumericVector breaks) {
// //   int nbreaks = breaks.size();
// //   NumericVector out(nbreaks-1);
// //   for (int i=0; i<nbreaks-1; i++) {
// //     LogicalVector temp = (x>breaks(i)) & (x<=breaks(i+1));
// //     out[i] = sum(temp);
// //   }
// // 
// //   return(out);
// // }
// 
// 
// 
// 
// 
// 
// 
// 
// //#include <Rcpp.h>
// //using namespace Rcpp;
// //
// //' Generate length composition of catch
// //'
// //' Generate size composition of catch given sample of catch-at-age,
// //' expected length-at-age, and standard deviation of length-at-age.
// //' Model assumes length-at-age is normally distributed, and that
// //' selectivity is size-dependant
// //'
// //' @param CAL_bins vector of catch-at-length size bins
// //' @param CAL_binsmid vector (nbins = length(CAL_bins) - 1) of mid-points for catch-at-length size bins
// //' @param SL matrix (nbins, nyears) of selectivity-at-length class for each year
// //' @param CAL_ESS effective sample size of catch-at-length data
// //' @param CAL_nsamp sample size of catch-at-length data
// //' @param CN matrix (nyears, maxage) of catch-at-age for each year
// //' @param LaA matrix (maxage, nyears) of expected length-at-age for each year
// //' @param LaASD matrix (maxage, nyears) of standard deviation of length-at-age for each year
// //' @param truncSD optional argument to truncate the length-at-age distribution at `truncSD` standard deviations
// //' e.g., a value of 2 truncates the length-at-age distribution at two standard deviations (set to 0 to ignore (default))
// //'
// //' @export
// // [[Rcpp::export]]
// 
// NumericMatrix  genLenComp(NumericVector CAL_bins, NumericVector CAL_binsmid, NumericMatrix SL,
//                           double CAL_ESS, double CAL_nsamp, NumericMatrix CN, NumericMatrix LaA, NumericMatrix LaASD,
//                           double truncSD) {
//   
//   int nbins = CAL_binsmid.size();
//   int nyears = CN.nrow();
//   int k = CN.ncol();
//   
//   NumericVector probs(k);
//   NumericMatrix CAL(nyears, nbins);
//   IntegerVector ans(k);
//   NumericVector ans2(k);
//   NumericVector Sx(k);
//   NumericMatrix Cx(k, nbins);
//   NumericMatrix Cx2(k, nbins);
//   
//   NumericMatrix Prob(k, nbins); // prob dist length at age
//   
//   for (int yr=0; yr < nyears; yr++) {
//     // sample ages from catch
//     NumericVector tempCN(k);
//     probs = CN.row(yr)/sum(CN.row(yr));
//     rmultinom(CAL_ESS, probs.begin(), k, ans.begin()); // multinom age sample with ess
//     ans2 = ans; // convert to numeric vector
//     tempCN = ans2 * (CAL_nsamp/CAL_ESS); // scale up from ess to sample size
//     
//     NumericVector EL(k);
//     NumericVector SDL(k);
//     
//     // probablility of length-at-age - assume normally distributed
//     EL = LaA.column(yr);  // expected length-at-age
//     SDL = LaASD.column(yr); // sd for length-at-age
//     for (int age=0; age < k; age++) {
//       Prob(age, 0) = R::pnorm((CAL_bins(1) - EL(age))/SDL(age), 0, 1, 1, 0);
//       for (int L=1; L < (nbins-1); L++) {
//         Prob(age, L) = R::pnorm((CAL_bins(L+1) - EL(age))/SDL(age), 0, 1, 1, 0) -
//           R::pnorm((CAL_bins(L) - EL(age))/SDL(age), 0, 1, 1, 0);
//       }
//       Prob(age, nbins-1) = 1 -  R::pnorm((CAL_bins(nbins-1) - EL(age))/SDL(age), 0, 1, 1,0);
//       
//       if (truncSD > 0 ) {
//         // Truncate at truncSD standard deviations
//         for (int L=0; L < nbins; L++) {
//           double temp = 0;
//           temp = (CAL_binsmid(L) - EL[age])/SDL(age);
//           if (temp >= truncSD) Prob(age,L) = 0;
//           if (temp <= -truncSD) Prob(age,L) = 0;
//         }
//       }
//       Sx(age) = sum(SL.column(yr) * Prob.row(age));	// Calculate selectivity at age given length-at-age dist
//       Cx.row(age) = Prob.row(age) * SL.column(yr);  // Length-at-age of catch conditional on selectivity
//       Cx2.row(age) = Cx.row(age)/sum(Cx.row(age));  // standardise so prob sum to 1 for each age
//     }
//     
//     for (int L=0; L < (nbins); L++) {
//       for (int age=0; age < k; age++) {
//         double tempVal = 0;
//         tempVal = (tempCN(age) * Cx2(age, L));
//         if (tempVal > 0) CAL(yr, L) += tempVal;
//       }
//     }
//     
//   }
//   
//   return(CAL);
// }
// 
