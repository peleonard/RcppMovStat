#include <Rcpp.h>
using namespace Rcpp;

//' @name movSumUE
//' @title Weighted Simple Moving Sum for Unevenly Spaced Time Series 
//' 
//' 
//' @description This function returns A matrix: the first column is the position, the second column 
//' the input vector, and third column moving sum of the given vector. The weight 
//' argument is optional. 
//'
//' @param vec A numeric vector.
//' @param pos A numeric vector with all integers. Its length must be the SAME as \eqn{vec}.
//' N.B. We use integers to represent the (relative) positions of every point.
//' 
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
//' This function can obtain the moving sum efficiently. It serves as somehow a generalized  
//' version of \code{\link{movMeanUE}}. The difference is that it will not automatically 
//' normalized the weights vector, \emph{w} argument. \cr
//' If there is no missing value in \emph{vec}, and \emph{w} is normalized, which means
//' the sum of all elements is \emph{1}, this function will return a moving average.
//' For matrix details, please refer to details of \code{movMeanUE}.
//'  
//' 
//' @return This function returns A MATRIX of size: \eqn{L*3}, where L is the length of vector, or 
//' of size: \eqn{L1*3}, where \eqn{L1 =} \code{\link[base]{ceiling}}\eqn{((nrow - n + 1)/ss)}, 
//' (when \eqn{sizeD = T}). In the matrix, the first column denotes the position, the second column the 
//' original vector, and the third column the moving average. 
//' 
//' 
//' @examples
//' movSumUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), 2)
//' movSumUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, na_rm = TRUE)
//' movSumUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = TRUE, 
//' sizeD = TRUE)
//' movSumUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), w = c(0, 1),  n = 2, 
//' ss = 3, na_rm = TRUE, align = "right")
//' movSumUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, 
//' na_rm = TRUE, sizeD = TRUE, align = "right")
//' movSumUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9),  n = 2, ss = 3, 
//' na_rm = TRUE)
//' movSumUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = TRUE, 
//' sizeD = TRUE)
//' @export
// [[Rcpp::export]]
NumericMatrix
  movSumUE(NumericVector vec, NumericVector pos, int n = 1L, int ss = 1L, Nullable<NumericVector> w = R_NilValue,
            bool na_rm = false, bool sizeD = false, std::string align = "left") {
    
    // Declare loop counts, nonNA count, size of vec(or pos), number of rows of matrix,  alignment position
    int i, j, k, n1 = 0, L = vec.size(), nrow, p, pos_init;
    // Declare three sums: sum, sum w/ weight, sum w/ weight and NA in vec
    double sum = 0.0;
    
    //create weight vector with default 1.0
    NumericVector wt(n, 1.0);
    
    // Initial matrix with NA
    pos_init = pos[0];
    pos = pos - pos_init + 1;
    nrow = pos[L-1] - pos[0] + 1;
    NumericMatrix posVecMs(nrow, 3); //A matrix with 3 columns: pos, vec and moving average
    NumericVector NAvec(nrow, NumericVector::get_na());
    for (j = 0; j < nrow; j++) {
      posVecMs(j, 0) = j + 1;
    }
    posVecMs(_, 1) = NAvec;
    posVecMs(_, 2) = NAvec;
    
    // Create the 2nd col based on the input
    for (i = 0; i < L; i++) {
      posVecMs(pos[i] - 1, 1) = vec[i];
    }
    
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
      for (i = p; i < nrow - n + p + 1; i = i + ss) {
        sum = 0.0;
        n1 = 0;
        k = 0;
        for (j = i - p; j < i - p + n; j++) {
          if (NumericVector::is_na(posVecMs(j, 1))) {
          } else {
            sum += wt[k] * posVecMs(j, 1);
            n1++;
          }
          k++;
        }
        if (n1 == 0) {
          posVecMs(i, 2) = NA_REAL;
        } else {
          posVecMs(i, 2) = sum;
        }
      }
      
      
    } else {
      for (i = p; i < nrow - n + p + 1; i = i + ss) {
        sum = 0.0;
        k = 0;
        for (j = i - p; j < i - p + n; j++) {
          sum += wt[k] * posVecMs(j, 1);
          k++;
        }
        posVecMs(i, 2) = sum;
      }
    }
    posVecMs(_, 0) = posVecMs(_, 0) + pos_init - 1; 
    
    if (sizeD) {
      int nrow1;
      nrow1 = (nrow - n + 1 + ss - 1)/ss;
      NumericMatrix posVecMs_sD(nrow1, 3);
      NumericVector NAvec_sD(nrow1, NumericVector::get_na());
      posVecMs_sD(_, 1) = NAvec_sD;
      posVecMs_sD(_, 2) = NAvec_sD;
      for (i = p; i < nrow - n + p + 1; i = i + ss){
        for (j = 0; j < posVecMs_sD.ncol(); j++) {
          posVecMs_sD((i-p)/ss, j) = posVecMs(i, j);
        }
      }
      return posVecMs_sD;
    } else {
      return posVecMs;
    }
  }




//' @describeIn movSumUE An function equivalent to \code{movSumUE(..., align = "right")}
//' @export
// [[Rcpp::export]]
NumericMatrix
  movSumUEr(NumericVector vec, NumericVector pos, int n = 1L, int ss = 1L,
             Nullable<NumericVector> w = R_NilValue, bool na_rm = false, bool sizeD = false) {
    
    // Declare loop counts, nonNA count, size of vec(or pos), number of rows of matrix
    int i, j, k, n1 = 0, L = vec.size(), nrow, pos_init;
    // Declare three sums: sum, sum w/ weight, sum w/ weight and NA in vec
    double sum = 0.0;
    
    //create weight vector with default 1.0
    NumericVector wt(n, 1.0);
    
    // Create matrix filled with NA
    pos_init = pos[0];
    pos = pos - pos_init + 1;
    nrow = pos[L-1] - pos[0] + 1;
    NumericMatrix posVecMs(nrow, 3); //A matrix with two vectors: moving average and position
    NumericVector NAvec(nrow, NumericVector::get_na());
    for (j = 0; j < nrow; j++) {
      posVecMs(j, 0) = j + 1;
    }
    posVecMs(_, 1) = NAvec;
    posVecMs(_, 2) = NAvec;
    
    
    // Create the 2nd col based on the input
    for (i = 0; i < L; i++) {
      posVecMs(pos[i] - 1, 1) = vec[i];
    }
    
    // std::cout<<posVecMs;
    
    if (w.isNotNull()) {
      wt = as<NumericVector>(w);
    }
    
    //Crux of the algorithm
    if (na_rm) {
      for (i = n - 1; i < nrow; i = i + ss) {
        sum = 0.0;
        n1 = 0;
        k = 0;
        for (j = i - n + 1; j < i + 1; j++) {
          if (NumericVector::is_na(posVecMs(j, 1))) {
          } else {
            sum += wt[k] * posVecMs(j, 1);
            n1++;
          }
          k++;
        }
        if (n1 == 0) {
          posVecMs(i, 2) = NA_REAL;
        } else {
          posVecMs(i, 2) = sum;
        }
      }
      
      
    } else {
      for (i = n - 1; i < nrow; i = i + ss) {
        sum = 0.0;
        k = 0;
        for (j = i - n + 1; j < i + 1; j++) {
          sum += wt[k] * posVecMs(j, 1);
          k++;
        }
        posVecMs(i, 2) = sum;
      }
    }
    posVecMs(_, 0) = posVecMs(_, 0) + pos_init - 1; 
    
    if (sizeD) {
      int nrow1;
      //This is a round-up method for a ratio of two integers: (nrow - n + 1)/ss
      nrow1 = (nrow - n + 1 + ss - 1)/ss; 
      NumericMatrix posVecMs_sD(nrow1, 3);
      NumericVector NAvec_sD(nrow1, NumericVector::get_na());
      posVecMs_sD(_, 1) = NAvec_sD;
      posVecMs_sD(_, 2) = NAvec_sD;
      for (i = n - 1; i < nrow; i = i + ss){
        for (j = 0; j < posVecMs_sD.ncol(); j++) {
          posVecMs_sD((i-n+1)/ss, j) = posVecMs(i, j);
        }
      }
      return posVecMs_sD;
    } else {
      return posVecMs;
    }
  }

/***R
# movSumUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), 2)
# movSumUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, na_rm = T)
movSumUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), w = c(0, 1),  n = 2, ss = 3, na_rm = T, align = "right")
# movSumUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = T, sizeD = T, align = "right")
movSumUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = T)
movSumUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = T, sizeD = T)
*/


