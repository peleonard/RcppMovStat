#include <Rcpp.h>
using namespace Rcpp;


//' @name movSum
//' @title Weighted Simple Moving Sum
//' 
//' 
//' @description This function returns a simple moving sum of the given vector. The weight 
//' argument is optional. 
//'
//'
//' @param vec A numeric vector.
//' @param n An integer: moving window size, with 1 as default
//' @param ss An integer: step size, only calculating at points with an equal distance \emph{ss}.
//' Namely, there are \emph{ss-1} number between each two 'consecutive' points    
//' @param w An optional weight vector of length \emph{n}. 
//' @param na_rm logical. Should missing values (including NaN) be removed?
//' @param sizeD logical. Only applied when \emph{ss > 1}, it decides whether to get a result of 
//' smaller size. If \eqn{sizeD = T}, \emph{align} does not affect the output. 
//' @param align A string denotes how to align the moving average, three options: 
//' "left", "middle", "right"
//' 
//' 
//' @details
//' This function can obtain the moving sum efficiently. It serves as somehow a generalized  
//' version of \code{\link{movMean}}. The difference is that it will not automatically 
//' normalized the weights vector, \emph{w} argument. \cr
//' If there is no missing value in \emph{vec}, and \emph{w} is normalized, which means
//' the sum of all elements is \emph{1}, this function will return a moving average.    
//' 
//' 
//' @return This function returns a vector whose length is the same as that of \emph{vec} or is 
//' \code{\link[base]{ceiling}}\eqn{((L - n + 1)/ss)}, (when \eqn{sizeD = T}), where \eqn{L} is the 
//' length of \eqn{vec}.
//' 
//' @examples
//' movSum(c(1, 4, 3, NA, 8), 3, align = "right", na_rm = TRUE)
//' movSum(c(1, 4, 3, NA, 8), 3, w = c(0.5, 0.2, 0.3), na_rm = TRUE, align = "right")
//' movSum(c(1, 4, 3, NA, 8, 4, 5, 9, 6, 0), n = 3, ss = 4, na_rm = TRUE, align = "right")
//' @export
// [[Rcpp::export]]
NumericVector 
  movSum(NumericVector vec, int n = 1L, int ss = 1L, Nullable<NumericVector> w = R_NilValue, 
          bool na_rm = false, bool sizeD = false, std::string align = "left") {
    
    // Declare loop counts, nonNA count, vector sizes, alignment position, sum
    int i, j, k, n1 = 0, L = vec.size(), p;
    double sum = 0.0;
    
    // Create vector filled with NA(R version)
    NumericVector ma(L, NumericVector::get_na()); //ma: moving average
    NumericVector wt(n, 1.0);
    
    if (w.isNotNull()) {
      wt = as<NumericVector>(w);
    }
    
    // assign value to p according to alignment rule
    if (align == "left") {
      p = 0;
    } else if (align == "middle") {
      p = n / 2;
    } else if (align == "right") {
      p = n - 1;
    } else {
      Rcpp::stop("Wrong alignment rule!");
    }
    
    //Crux of the algorithm
    if (na_rm) {
      for (i = p; i < L - n + p + 1; i = i + ss) {
        sum = 0.0; 
        n1 = 0;
        k = 0;
        for (j = i - p; j < i - p + n; j++) {
          if (NumericVector::is_na(vec[j])) {
          } else {
            sum += wt[k] * vec[j];
            n1++; 
          }
          k++;
        }
        if (n1 == 0) {
          ma[i] = NA_REAL;
        } else {
          ma[i] = sum;
        }
      }
      
      
    } else {
      for (i = p; i < L - n + p + 1; i = i + ss) {
        sum = 0.0;
        k = 0;
        for (j = i - p; j < i - p + n; j++) {
          sum += wt[k] * vec[j];
          k++;
        }
        ma[i] = sum;
      }
    }
    
    if (sizeD) {
      int L1;
      L1 = (L - n + 1 + ss - 1)/ss; //This is a round-up method for a ratio of two integers: (L - n + 1)/ss
      NumericVector ma_sD(L1);
      for (i = p; i < L - n + p + 1; i = i + ss){
        ma_sD[(i-p)/ss] = ma[i];
      }
      return ma_sD;
    } else {
      return ma;
    }
    
  }


//' @describeIn movSum An function equivalent to \code{movMean(..., align = "right")}
//' @export
// [[Rcpp::export]]
NumericVector 
  movSumr(NumericVector vec, int n = 1L, int ss = 1L,
           Nullable<NumericVector> w = R_NilValue, bool na_rm = false, bool sizeD = false) {
    
    // Declare loop counts, nonNA count, vector sizes, alignment position, sum
    int i, j, k, n1 = 0, L = vec.size();
    double sum = 0.0;
    
    // Create vector filled with NA(R version)
    NumericVector ma(L, NumericVector::get_na()); //ma: moving average
    NumericVector wt(n, 1.0);
    
    if (w.isNotNull()) {
      wt = as<NumericVector>(w);
    }
    
    
    //Crux of the algorithm
    if (na_rm) {
      for (i = n - 1; i < L; i = i + ss) {
        sum = 0.0; 
        n1 = 0;
        k = 0;
        for (j = i - n + 1; j < i + 1; j++) {
          if (NumericVector::is_na(vec[j])) {
          } else {
            sum += wt[k] * vec[j];
            n1++; 
          }
          k++;
        }
        if (n1 == 0) {
          ma[i] = NA_REAL;
        } else {
          ma[i] = sum;
        }
      } 
      
      
    } else {
      for (i = n - 1; i < L; i = i + ss) {
        sum = 0.0;
        k = 0;
        for (j = i - n + 1; j < i + 1; j++) {
          sum += wt[k] * vec[j];
          k++;
        }
        ma[i] = sum;
      } 
    }
    
    if (sizeD) {
      int L1;
      L1 = (L - n + 1 + ss - 1)/ss; //This is a round-up method for a ratio of two integers: (L - n + 1)/ss
      NumericVector ma_sD(L1);
      for (i = n - 1; i < L; i = i + ss){
        ma_sD[(i- n + 1)/ss] = ma[i];
      }
      return ma_sD;
    } else {
      return ma;
    }
  }


/*** R
movSum(c(1, 4, 3, NA, 8), 2)
movSum(c(1, 4, 3, NA, 8), 2, na_rm = TRUE)
movSum(c(1, 4, 3, NA, NA), 2, na_rm = TRUE)
movSum(c(1, 4, 3, NA, NA), 2, na_rm = TRUE, align = 'right')
movSumr(c(1, 4, 3, NA, NA), 2, na_rm = TRUE)

*/
