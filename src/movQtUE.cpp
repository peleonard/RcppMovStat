#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// quantile_UE setup
struct quantile_UE {
  int lo, hi, n, total;
  double hlo, hhi;
};

quantile_UE make_quantile_UE(int n, double prob) {
  quantile_UE res; 
  double index = (n - 1) * prob;
  res.lo = floor(index);
  res.hi = res.lo + 1;
  res.hhi = index - res.lo;
  res.hlo = 1 - res.hhi;
  return res;	
}

//' @name movQtUE
//' @title Moving quantile_UE(Moving Median, Moving Minimum, Moving Maximum)
//' for Unevenly Spaced Time Series 
//' 
//' @description This function returns A matrix: the first column is the position, the second column 
//' the input vectorthe, and third column moving quantile_UE of the given vector. 
//'
//' @param vec A numeric vector.
//' @param pos A numeric vector with all integers. Its length must be the SAME as \eqn{vec}.
//' N.B. We use integers to represent the (relative) postions of every point.
//' The first element MUST BE 1, which is design to follow THE 1-INDEXED RULE of R.
//' 
//' @param n An integer: moving window size, with 1 as default
//' @param prob A number: between \emph{0} and \emph{1}, meaning \emph{prob} quantile_UE
//' @param ss An integer: step size, only calculating at points with an equal distance \emph{ss}.
//' Namely, there are \emph{ss-1} number between each two 'consecutive' points    
//' @param na_rm logical. Should missing values (including NaN) be removed?
//' @param sizeD logical. Only applied when \emph{ss > 1}, it decides whether to get a result of 
//' smaller size. If \eqn{sizeD = T}, \emph{align} does not affect the output. 
//' @param align A string denotes how to align the moving average, three options: 
//' "left", "middle", "right"
//' 
//' 
//' @details
//' This function is especially designed for Unevenly Spaced Time Series. It is efficient as it herits the 
//' similiar routine of \code{movQt}. \cr
//' The result is kind of tricky. To make it clear, it is written to return a MATRIX. For instance, the 
//' third column of the output of second example is \eqn{2.5, NA, NA, NA, NA, NA, 3.0, NA, NA}. 
//' 2.5 is the median of 1 and 4, and 4.0 the average of 4. The third column of the output of 
//' third example is the every third element starting from \eqn{n}th number. \cr
//' For how weights, \eqn{w}, work, one can refer to \code{movQt}. 
//' 
//' 
//' @return This function returns A MATRIX of size: \eqn{L*3}, where L is the length of vector, or 
//' of size: \eqn{L1*3}, where \eqn{L1 =} \code{\link[base]{ceiling}}\eqn{((nrow - n + 1)/ss)}, 
//' (when \eqn{sizeD = T}). In the matrix, the first column denotes the position, the second column the 
//' original vector, and the third column the moving average. 
//' 
//' @examples
//' movQtUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = TRUE)
//' movQtUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = TRUE, sizeD = TRUE)
//' movQtUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, na_rm = TRUE, align = "middle")
//' movQtUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = TRUE, sizeD = TRUE,
//'  align = "right")
//' @export
// [[Rcpp::export]]
NumericMatrix
  movQtUE(NumericVector vec, NumericVector pos, int n = 1L, double prob = .5, int ss = 1L,
        bool na_rm = false, bool sizeD = false, std::string align = "left") {
    
    // Declare loop counters, nonNA counter, vector sizes and sum
    int i, k, j, p, nrow;
    int L = vec.size();
    
    // Initial matrix with NA
    nrow = pos[L-1] - pos[0] + 1;
    NumericMatrix posVecRq(nrow, 3); //A matrix with 3 columns: pos, vec and moving average
    NumericVector NAvec(nrow, NumericVector::get_na());
    for (j = 0; j < nrow; j++) {
      posVecRq(j, 0) = j + 1;
    }
    posVecRq(_, 1) = NAvec;
    posVecRq(_, 2) = NAvec;
    
    // Create the 2nd col based on the input
    for (i = 0; i < L; i++) {
      posVecRq(pos[i] - 1, 1) = vec[i];
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
        k = 0;
        std::vector<double> r(n);
        for (j = i - p; j < i - p + n; j++) {
          if (NumericVector::is_na(posVecRq(j, 1))) {
          } else {
            r[k] = posVecRq(j, 1);
            k++;
          }
        }
        if (k == 0) {
          posVecRq(i, 2) = NA_REAL;
        } else {
          r.resize(k); 
          auto q = make_quantile_UE(k, prob);
          vector<double> z(&r[0], &r[0+k]);
          sort(z.begin(), z.end());
          posVecRq(i, 2) = q.hlo * z[q.lo] + q.hhi * z[q.hi]; 
        }
      } 
      
    } else {
      // quantile_UE setup
      auto q = make_quantile_UE(n, prob);
      for (i = p; i < nrow - n + p + 1; i = i + ss) {
        vector<double> z(&posVecRq(_, 1)[i - p], &posVecRq(_, 1)[i - p + n]);
        bool Ind = false;
        for (j = i - p; j < i - p + n; j++) {
          Ind = Ind|NumericVector::is_na(posVecRq(j, 1));
        }
        if (Ind) {
          posVecRq(i, 2) = NA_REAL;
        } else {
          sort(z.begin(), z.end());
          posVecRq(i, 2) = q.hlo * z[q.lo] + q.hhi * z[q.hi];  
        }
      } 
    }
    if (sizeD) {
      int nrow1;
      nrow1 = (nrow - n + 1 + ss - 1)/ss;
      NumericMatrix posVecRq_sD(nrow1, 3);
      NumericVector NAvec_sD(nrow1, NumericVector::get_na());
      posVecRq_sD(_, 1) = NAvec_sD;
      posVecRq_sD(_, 2) = NAvec_sD;
      for (i = p; i < nrow - n + p + 1; i = i + ss){
        for (j = 0; j < posVecRq_sD.ncol(); j++) {
          posVecRq_sD((i-p)/ss, j) = posVecRq(i, j);
        }
      }
      return posVecRq_sD;
    } else {
      return posVecRq;
      }
    
    }

//' @describeIn movQtUE An function equivalent to \code{movQtUE(..., align = "right")}
//' @export
// [[Rcpp::export]]
NumericMatrix
  movQtUEr(NumericVector vec, NumericVector pos, int n = 1L, double prob = .5, int ss = 1L,
           bool na_rm = false, bool sizeD = false) {

    // Declare loop counters, nonNA counter, vector sizes and sum
    int i, k, j, nrow;
    int L = vec.size();

    // Initial matrix with NA
    nrow = pos[L-1] - pos[0] + 1;
    NumericMatrix posVecRq(nrow, 3); //A matrix with 3 columns: pos, vec and moving average
    NumericVector NAvec(nrow, NumericVector::get_na());
    for (j = 0; j < nrow; j++) {
      posVecRq(j, 0) = j + 1;
    }
    posVecRq(_, 1) = NAvec;
    posVecRq(_, 2) = NAvec;
    
    // Create the 2nd col based on the input
    for (i = 0; i < L; i++) {
      posVecRq(pos[i] - 1, 1) = vec[i];
    }

    //Crux of the algorithm
    if (na_rm) {
      for (i = n - 1; i < nrow; i = i + ss) {
        k = 0;
        std::vector<double> r(n);
        for (j = i - n + 1; j < i + 1; j++) {
          if (NumericVector::is_na(posVecRq(j, 1))) {
          } else {
            r[k] = posVecRq(j, 1);
            k++;
          }
        }
        if (k == 0) {
          posVecRq(i, 2) = NA_REAL;
        } else {
          r.resize(k);
          auto q = make_quantile_UE(k, prob);
          vector<double> z(&r[0], &r[0+k]);
          sort(z.begin(), z.end());
          posVecRq(i, 2) = q.hlo * z[q.lo] + q.hhi * z[q.hi];
        }
      }

    } else {
      // quantile_UE setup
      auto q = make_quantile_UE(n, prob);
      for (i = n - 1; i < L; i = i + ss) {
        vector<double> z(&posVecRq(_, 1)[i - n + 1], &posVecRq(_, 1)[i + 1]);
        bool Ind = false;
        for (j = i - n + 1; j < i + 1; j++) {
          Ind = Ind|NumericVector::is_na(posVecRq(j, 1));
        }
        if (Ind) {
          posVecRq(i, 2) = NA_REAL;
        } else {
          sort(z.begin(), z.end());
          posVecRq(i, 2) = q.hlo * z[q.lo] + q.hhi * z[q.hi];
        }
      }
    }

    if (sizeD) {
      int nrow1;
      nrow1 = (nrow - n + 1 + ss - 1)/ss;
      NumericMatrix posVecRq_sD(nrow1, 3);
      NumericVector NAvec_sD(nrow1, NumericVector::get_na());
      posVecRq_sD(_, 1) = NAvec_sD;
      posVecRq_sD(_, 2) = NAvec_sD;
      for (i = n - 1; i < nrow; i = i + ss){
        for (j = 0; j < posVecRq_sD.ncol(); j++) {
          posVecRq_sD((i-n+1)/ss, j) = posVecRq(i, j);
        }
      }
      return posVecRq_sD;
    } else {
      return posVecRq;
      }
    
    }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
# # movQtUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), 2)
# movQtUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = T)
# movQtUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = T, sizeD = T)
# 
# movQtUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, na_rm = T, align = "right")
# movQtUE(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = T, sizeD = T, align = "right")
# # movQtUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = T)
# # movQtUEr(c(1, 4, 3, NA, 8), pos = c(1, 2, 7, 8, 9), n = 2, ss = 3, na_rm = T, sizeD = T)
*/
  