#include <Rcpp.h>
using namespace Rcpp;

//' @name movEmean
//' @title Basic Exponential Moving Mean 
//' 
//' 
//' @description This function returns a basic exponential moving average (EMA) of the given vector.
//'
//'
//' @param vec A numeric vector.
//' @param n An integer: moving window size, with \eqn{1} as default
//' @param smFac A number: smoothing factor, with default \eqn{2/(n+1)},  see \strong{details} below. 
//' 
//' 
//' @details
//' This function makes fairly efficient the computation of EMA, which dubbed as basic 
//' exponential smoothing, the same section of \url{https://en.wikipedia.org/wiki/Exponential_smoothing}.
//' It provides an access to define \emph{smFac} yourself, i.e the smoothing factor, 
//' whose default is \eqn{2/(n+1)}.
//' 
//' 
//' @return This function returns a vector whose length is the same as that of \emph{vec}.
//' 
//' 
//' @examples
//' movEmean(c(1, 4, 3, 6, 8), 2, smFac = 1/3)
//' movEmean(c(1, 4, 3, 6, 8), 2)
//' @export
// [[Rcpp::export]]
NumericVector 
  movEmean(NumericVector vec, int n = 1L, Nullable<double> smFac = R_NilValue) {
    
    // Declare loop counts, nonNA count, vector sizes
    int i, j, k, L = vec.size();
    // Declare sum and smoothing constant of Exponential Moving Average
    double sum = 0.0, sc;
    
    // Create vector filled with NA(R version)
    NumericVector ma(L, NumericVector::get_na()); //ma: moving average
    NumericVector ema(L, NumericVector::get_na()); //ema: exponential moving average
    
    //Crux of the algorithm
    for (i = n - 1; i < L; i++) {
      sum = 0.0; 
      for (j = i - n - 1; j < i + 1; j++) {
        sum += vec[j];
      }
      ma[i] = sum / n;
    } 
    
    
    if (smFac.isNull()){
      sc = (double) 2 / (n + 1); //smoothing constant
    } else {
      sc = as<double>(smFac);
    }
    
    for (k = n - 1; k < L ; k++) {
      if (k == n - 1) {
        ema[k] = ma[k];
      } else {
        ema[k] = sc*vec[k] + (1-sc)*ema[k-1];
      }
    }
    // return result
    return ema;
  }



/*** R
movEmean(c(1, 4, 3, 6, 8), 2, smFac = 1/3)
movEmean(c(1, 4, 3, 6, 8), 2)
TTR::EMA(c(1, 4, 3, 6, 8), 2)

*/
