#include <Rcpp.h>
using namespace Rcpp;

//' @name movMean
//' @title Weighted Simple Moving Mean
//' 
//' 
//' @description This function returns a simple moving average of the given vector. The weight 
//' argument is optional. 
//'
//'
//' @param vec A numeric vector.
//' @param n An integer: moving window size, with 1 as default
//' @param ss An integer: step size, only calculating at points with an equal distance \emph{ss}.
//' Namely, there are \emph{ss-1} number between each two 'consecutive' points    
//' @param w An optional weight vector of length \emph{n}. It will be automatically normalized
//' (sum to 1). 
//' @param na_rm logical. Should missing values (including NaN) be removed?
//' @param sizeD logical. Only applied when \emph{ss > 1}, it decides whether to get a result of 
//' smaller size. If \eqn{sizeD = T}, \emph{align} does not affect the output. 
//' @param align A string denotes how to align the moving average, three options: 
//' "left", "middle", "right"
//' 
//' 
//' @details
//' Despite of Efficient computation, usually 5~6 times faster than the moving average function in 
//' package \code{\link[RcppRoll]{roll_mean}}, it is able to handle potential missing values
//' (\emph{NA} or \emph{NaN}) in the \emph{vec}. \cr
//' For instace, the output of the second examle is \eqn{NA, NA, 2.200000 3.714286 4.875000}. The 
//' last number \eqn{5.5} is obtained by using renormalized weight, namely omitting \eqn{0.2}.
//' The weight applied would be \eqn{0.5/(0.5+0.3)} and \eqn{0.3/(0.5+0.3)}. Hence, 
//' \deqn{4.875 = 3 * 0.5/(0.5+0.3) + 8 * 0.3/(0.5+0.3)}
//' 
//' 
//' @return This function returns a vector whose length is the same as that of \emph{vec} or is 
//' \code{\link[base]{ceiling}}\eqn{((L - n + 1)/ss)}, (when \eqn{sizeD = T}), where \eqn{L} is the 
//' length of \eqn{vec}.
//' 
//' 
//' @examples
//' movMean(c(1, 4, 3, NA, 8), 3, align = "right", na_rm = TRUE)
//' movMean(c(1, 4, 3, NA, 8), 3, w = c(0.5, 0.2, 0.3), na_rm = TRUE, align = "right")
//' movMean(c(1, 4, 3, NA, 8, 4, 5, 9, 6, 0), n = 3, ss = 4, na_rm = TRUE, align = "right")
//' @export
// [[Rcpp::export]]
NumericVector 
  movMean(NumericVector vec, int n = 1L, int ss = 1L, Nullable<NumericVector> w = R_NilValue, 
          bool na_rm = false, bool sizeD = false, std::string align = "left") {
    
    // Declare loop counts, nonNA count, vector sizes, alignment position, sum
    int i, j, k, n1 = 0, L = vec.size(), p;
    double sum = 0.0, sum_wt, sum_wtNA;
    
    // Create vector filled with NA(R version)
    NumericVector ma(L, NumericVector::get_na()); //ma: moving average
    NumericVector wt(n, 1.0);
    
    if (w.isNotNull()) {
      wt = as<NumericVector>(w);
    } 
    sum_wt = 0.0;
    for (k = 0; k < n; k++) {
      sum_wt += wt[k];
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
        sum_wtNA = 0.0;
        for (j = i - p; j < i - p + n; j++) {
          if (NumericVector::is_na(vec[j])) {
          } else {
            sum += wt[k] * vec[j];
            sum_wtNA += wt[k];
            n1++; 
          }
          k++;
        }
        if (n1 == 0) {
          ma[i] = NA_REAL;
        } else {
          ma[i] = sum/sum_wtNA;
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
        ma[i] = sum / sum_wt;
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




//' @describeIn movMean An function equivalent to \code{movMean(..., align = "right")}
//' @export
// [[Rcpp::export]]
NumericVector 
  movMeanr(NumericVector vec, int n = 1L, int ss = 1L,
           Nullable<NumericVector> w = R_NilValue, bool na_rm = false, bool sizeD = false) {
    
    // Declare loop counts, nonNA count, vector sizes, alignment position, sum
    int i, j, k, n1 = 0, L = vec.size();
    double sum = 0.0, sum_wt, sum_wtNA;
    
    // Create vector filled with NA(R version)
    NumericVector ma(L, NumericVector::get_na()); //ma: moving average
    NumericVector wt(n, 1.0);
    
    if (w.isNotNull()) {
      wt = as<NumericVector>(w);
    } 
    sum_wt = 0.0;
    for (k = 0; k < n; k++) {
      sum_wt += wt[k];
    }

    
    //Crux of the algorithm
    if (na_rm) {
      for (i = n - 1; i < L; i = i + ss) {
        sum = 0.0; 
        n1 = 0;
        k = 0;
        sum_wtNA = 0.0;
        for (j = i - n + 1; j < i + 1; j++) {
          if (NumericVector::is_na(vec[j])) {
          } else {
            sum += wt[k] * vec[j];
            sum_wtNA += wt[k];
            n1++; 
          }
          k++;
        }
        if (n1 == 0) {
          ma[i] = NA_REAL;
        } else {
          ma[i] = sum/sum_wtNA;
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
        ma[i] = sum / sum_wt;
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

