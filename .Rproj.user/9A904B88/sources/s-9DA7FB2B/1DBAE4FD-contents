#include <Rcpp.h>
using namespace Rcpp;

//' @name movCountUE
//' @title Weighted Simple Moving Count for Unevenly Spaced Time Series 
//' 
//' 
//' @description This function returns A matrix: the first column is the position, the second column 
//' the input vectorthe, and third column Moving Count of the given vector. The weight 
//' argument is optional. 
//'
//' @param vec A numeric vector.
//' @param pos A numeric vector with all integers. Its length must be the SAME as \eqn{vec}.
//' N.B. We use integers to represent the (relative) postions of every point.
//' The first element MUST BE 1, which is design to follow THE 1-INDEXED RULE of R.
//' 
//' @param n An integer: moving window size, with 1 as default
//' @param ss An integer: step size, only calculating at points with an equal distance \emph{ss}.
//' Namely, there are \emph{ss-1} number between each two 'consecutive' points    
//' @param na_rm logical. Should missing values (including NaN) be removed?
//' @param sizeD logical. Only applied when \emph{ss > 1}, it decides whether to get a result of 
//' smaller size. If \eqn{sizeD = T}, \emph{align} does not affect the output. 
//' @param align A string denotes how to align the moving average, three options: 
//' "left", "middle", "right"
//' 
//' @details This function counts the number of non-missing values for each moving window. It is 
//' especially designed for \emph{vec} vector with missing values. Otherwise, it will return a trival 
//' vector with all elements \emph{n}. \cr
//' This function is more helpful than \code{movCount}, as we would have missing values for an Unevenly 
//' Spaced Time Series. \cr
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
//' movCountUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), 2)
//' movCountUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, na_rm = TRUE)
//' movCountUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = TRUE, 
//' sizeD = TRUE)
//' movCountUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, 
//' ss = 3, na_rm = TRUE, align = "right")
//' movCountUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, 
//' na_rm = TRUE, sizeD = TRUE, align = "right")
//' movCountUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9),  n = 2, ss = 3, 
//' na_rm = TRUE)
//' movCountUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = TRUE, 
//' sizeD = TRUE)
//' @export
// [[Rcpp::export]]
NumericMatrix
  movCountUE(NumericVector vec, NumericVector pos, int n = 1L, int ss = 1L,
         bool na_rm = false, bool sizeD = false, std::string align = "left") {
    
    // Declare loop counts, nonNA count, size of vec(or pos), number of rows of matrix,  alignment position
    int i, j, n1 = 0, L = vec.size(), nrow, p;
    
    // Initial matrix with NA
    nrow = pos[L-1] - pos[0] + 1;
    NumericMatrix posVecMc(nrow, 3); //A matrix with 3 columns: pos, vec and moving average
    NumericVector NAvec(nrow, NumericVector::get_na());
    for (j = 0; j < nrow; j++) {
      posVecMc(j, 0) = j + 1;
    }
    posVecMc(_, 1) = NAvec;
    posVecMc(_, 2) = NAvec;
    
    // Create the 2nd col based on the input
    for (i = 0; i < L; i++) {
      posVecMc(pos[i] - 1, 1) = vec[i];
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
        n1 = 0;
        for (j = i - p; j < i - p + n; j++) {
          if (NumericVector::is_na(posVecMc(j, 1))) {
          } else {
            n1++;
          }
        }
        if (n1 == 0) {
          posVecMc(i, 2) = NA_REAL;
        } else {
          posVecMc(i, 2) = n1;
        }
      }
      
      
    } else {
      for (i = p; i < nrow - n + p + 1; i = i + ss) {
        bool Ind = false;
        for (j = i - p; j < i - p + n; j++) {
          Ind = Ind|NumericVector::is_na(posVecMc(j, 1));
        }
        if (Ind) {
          posVecMc(i, 2) = NA_REAL;
        } else {
          posVecMc(i, 2) = n;
        }
      }
    }
    
    if (sizeD) {
      int nrow1;
      nrow1 = (nrow - n + 1 + ss - 1)/ss;
      NumericMatrix posVecMc_sD(nrow1, 3);
      NumericVector NAvec_sD(nrow1, NumericVector::get_na());
      posVecMc_sD(_, 1) = NAvec_sD;
      posVecMc_sD(_, 2) = NAvec_sD;
      for (i = p; i < nrow - n + p + 1; i = i + ss){
        for (j = 0; j < posVecMc_sD.ncol(); j++) {
          posVecMc_sD((i-p)/ss, j) = posVecMc(i, j);
        }
      }
      return posVecMc_sD;
    } else {
      return posVecMc;
    }
  }




//' @describeIn movCountUE An function equivalent to \code{movCountUE(..., align = "right")}
//' @export
// [[Rcpp::export]]
NumericMatrix
  movCountUEr(NumericVector vec, NumericVector pos, int n = 1L, int ss = 1L,
              bool na_rm = false, bool sizeD = false) {
    
    // Declare loop counts, nonNA count, size of vec(or pos), number of rows of matrix
    int i, j, n1 = 0, L = vec.size(), nrow;
    
    //create weight vector with default 1.0
    NumericVector wt(n, 1.0);
    
    // Create matrix filled with NA
    nrow = pos[L-1] - pos[0] + 1;
    NumericMatrix posVecMc(nrow, 3); //A matrix with two vectors: moving average and position
    NumericVector NAvec(nrow, NumericVector::get_na());
    for (j = 0; j < nrow; j++) {
      posVecMc(j, 0) = j + 1;
    }
    posVecMc(_, 1) = NAvec;
    posVecMc(_, 2) = NAvec;
    
    
    // Create the 2nd col based on the input
    for (i = 0; i < L; i++) {
      posVecMc(pos[i] - 1, 1) = vec[i];
    }
  
    
    //Crux of the algorithm
    if (na_rm) {
      for (i = n - 1; i < nrow; i = i + ss) {
        n1 = 0;
        for (j = i - n + 1; j < i + 1; j++) {
          if (NumericVector::is_na(posVecMc(j, 1))) {
          } else {
            n1++;
          }
        }
        if (n1 == 0) {
          posVecMc(i, 2) = NA_REAL;
        } else {
          posVecMc(i, 2) = n1;
        }
      }
      
      
    } else {
      for (i = n - 1; i < nrow; i = i + ss) {
        bool Ind = false;
        for (j = i - n + 1; j < i + 1; j++) {
          Ind = Ind|NumericVector::is_na(posVecMc(j, 1));
        }
        if (Ind) {
          posVecMc(i, 2) = NA_REAL;
        } else {
          posVecMc(i, 2) = n;
        }
      }
    }
    
    if (sizeD) {
      int nrow1;
      //This is a round-up method for a ratio of two integers: (nrow - n + 1)/ss
      nrow1 = (nrow - n + 1 + ss - 1)/ss; 
      NumericMatrix posVecMc_sD(nrow1, 3);
      NumericVector NAvec_sD(nrow1, NumericVector::get_na());
      posVecMc_sD(_, 1) = NAvec_sD;
      posVecMc_sD(_, 2) = NAvec_sD;
      for (i = n - 1; i < nrow; i = i + ss){
        for (j = 0; j < posVecMc_sD.ncol(); j++) {
          posVecMc_sD((i-n+1)/ss, j) = posVecMc(i, j);
        }
      }
      return posVecMc_sD;
    } else {
      return posVecMc;
    }
  }

/***R
# movCountUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), 2)
movCountUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, na_rm = T)
movCountUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9),  n = 2, ss = 3, na_rm = T, align = "right")
# movCountUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = T, sizeD = T, align = "right")
  movCountUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = T)
  movCountUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = T, sizeD = T)
  */


