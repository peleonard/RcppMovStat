#include <Rcpp.h>
using namespace Rcpp;


//' @name movCount
//' @title Moving Count
//' 
//' 
//' @description This function returns a moving count of the given vector. 
//'
//'
//' @param vec A numeric vector.
//' @param n An integer: moving window size, with 1 as default
//' @param ss An integer: step size, only calculating at points with an equal distance \emph{ss}.
//' Namely, there are \emph{ss-1} number between each two 'consecutive' points 
//' @param na_rm logical. Should missing values (including NaN) be removed?
//' @param sizeD logical. Only applied when \emph{ss > 1}, it decides whether to get a result of 
//' smaller size. If \eqn{sizeD = T}, \emph{align} does not affect the output. 
//' @param align A string denotes how to align the moving average, three options: 
//' "left", "middle", "right"
//' 
//' 
//' @details This function counts the number of non-missing values for each moving window. It is 
//' especially designed for \emph{vec} vector with missing values. Otherwise, it will return a trivial 
//' vector with all elements \emph{n}.
//' 
//' 
//' @return This function returns a vector whose length is the same as that of \emph{vec} or is 
//'\code{\link[base]{ceiling}}\eqn{((L - n + 1)/ss)}, (when \eqn{sizeD = T}), where \eqn{L} is the 
//' length of \eqn{vec}.
//' 
//' @examples
//' movCount(c(1, 4, 3, NA, 8), 2, na_rm = TRUE)
//' movCount(c(1, 4, 3, NA, 8), 2, na_rm = TRUE, align = 'right')
//' movCountr(c(1, 4, 3, NA, 8), 2, na_rm = TRUE)
//' movCount(c(1, 4, 3, NA, NA), 2, na_rm = TRUE)
//' @export
// [[Rcpp::export]]
NumericVector 
  movCount(NumericVector vec, int n = 1L, int ss = 1L, 
              bool na_rm = false, bool sizeD = false, std::string align = "left") {
    
    // Declare loop counters, nonNA counter, vector sizes and sum
    int i, j, n1, L = vec.size(), p;
    
    // Create vector filled with NA(R version)
    NumericVector mc(L, NumericVector::get_na());
    
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
        n1 = 0;
        for (j = i - p; j < i - p + n; j++) {
          if (NumericVector::is_na(vec[j])) {
          } else {
            n1++;
          }
        }
        mc[i] = n1;
      } 
    } else {
      for (i = p; i < L - n + p + 1; i = i + ss) {
        bool Ind = false;
        for (j = i - p; j < i - p + n; j++) {
          Ind = Ind|NumericVector::is_na(vec[j]);
        }
        if (Ind) {
          mc[i] = NA_REAL;
        } else {
          mc[i] = n;
        }
      } 
    }
    
    if (sizeD) {
      int L1;
      L1 = (L - n + 1 + ss - 1)/ss; //This is a round-up method for a ratio of two integers: (L - n + 1)/ss
      NumericVector mc_sD(L1);
      for (i = p; i < L - n + p + 1; i = i + ss){
        mc_sD[(i-p)/ss] = mc[i];
      }
      return mc_sD;
    } else {
      return mc;
    }
  }


//' @describeIn movCount An function equivalent to \code{movCount(..., align = "right")}
//' @export
// [[Rcpp::export]]
NumericVector 
  movCountr(NumericVector vec, int n = 1L, int ss = 1L,  
           bool na_rm = false, bool sizeD = false) {
    
    // Declare loop counters, nonNA counter, vector sizes and sum
    int i, j, n1, L = vec.size();
    
    // Create vector filled with NA(R version)
    NumericVector mc(L, NumericVector::get_na());
    
    
    //Crux of the algorithm
    if (na_rm) {
      for (i = n - 1; i < L; i++) {
        n1 = 0;
        for (j = i - n + 1; j < i + 1; j++) {
          if (NumericVector::is_na(vec[j])) {
          } else {
            n1++;
          }
        }
        mc[i] = n1;
      } 
    } else {
      for (i = n - 1; i < L; i++) {
        bool Ind = false;
        for (j = i - n + 1; j < i + 1; j++) {
          Ind = Ind|NumericVector::is_na(vec[j]);
        }
        if (Ind) {
          mc[i] = NA_REAL;
        } else {
          mc[i] = n;
        }
      } 
    }
    
    if (sizeD) {
      int L1;
      L1 = (L - n + 1 + ss - 1)/ss; //This is a round-up method for a ratio of two integers: (L - n + 1)/ss
      NumericVector mc_sD(L1);
      for (i = n - 1; i < L; i = i + ss){
        mc_sD[(i-n+1)/ss] = mc[i];
      }
      return mc_sD;
    } else {
      return mc;
    }
  }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
movCount(c(1, 4, 3, NA, 8), 2)

movCount(c(1, 4, 3, NA, 8), 2, na_rm = TRUE)

movCount(c(1, 4, 3, NA, 8), 2, na_rm = TRUE, align = 'right')
movCountr(c(1, 4, 3, NA, 8), 2, na_rm = TRUE)

movCount(c(1, 4, 3, NA, NA), 2, na_rm = TRUE)
*/
