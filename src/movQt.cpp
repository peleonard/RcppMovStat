#include <Rcpp.h>
#include "shared/makeQt.hpp"
using namespace Rcpp;
using namespace std;



//' @name movQt
//' @title Moving Quantile(Moving Median, Moving Minimum, Moving Maximum)
//' 
//' 
//' @description This function returns a moving quantile of the given vector. 
//'
//'
//' @param vec A numeric vector.
//' @param n An integer: moving window size, with 1 as default
//' @param prob A number: between \emph{0} and \emph{1}, meaning \emph{prob} quantile
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
//' Despite of Efficient computation, this function can return differents kinds of moving quantile,
//' e.g. moving median(\eqn{prob = 0.5}), moving minimun(\eqn{prob = 0}), 
//' and moving maximum(\eqn{prob = 1}). It can handle potential 
//' missing values(\emph{NA} or \emph{NaN}) in the \emph{vec}. When we move to one specific fragment, missing 
//' values can be removed by setting \eqn{na_rm = TRUERUE}. If all values of this fragment is missing, 
//' it will return \emph{NA}.  \cr
//' In terms of the quantile algorithm, please consult type \emph{7} in function 
//' \code{\link[stats]{quantile}}.
//' 
//' 
//' @return This function returns a vector whose length is the same as that of \emph{vec} or is 
//' \code{\link[base]{ceiling}}\eqn{((L - n + 1)/ss)}, (when \eqn{sizeD = T}), where \eqn{L} is the 
//' length of \eqn{vec}.
//' 
//' @examples
//' movQt(vec = c(1, 4, 3, NA, 8, 4, 5, 9, 6, 0), n = 3, ss = 4, na_rm = TRUE, align = "right")
//' movQt(vec = c(1, 4, 3, NA, 8, 4, 5, 9, 6, 0), n = 3, na_rm = TRUE, align = "right")
//' movQt(vec = c(1, 4, 3, NA, NA, NA, 5, 9, 6, 0), n = 3, ss = 4, na_rm = TRUE, align = "middle")
//' @export
// [[Rcpp::export]]
NumericVector 
  movQt(NumericVector vec, int n = 1L, double prob = .5, int ss = 1L,
        bool na_rm = false, bool sizeD = false, std::string align = "left") {
    
    // Declare loop counters, nonNA counter, vector sizes and sum
    int i, k, j, p;
    int L = vec.size();
    
    // Create vector filled with NA(R version)
    NumericVector rq(L, NumericVector::get_na());
    
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
        k = 0;
        std::vector<double> r(n);
        for (j = i - p; j < i - p + n; j++) {
          if (NumericVector::is_na(vec[j])) {
          } else {
            r[k] = vec[j];
            k++;
          }
        }
        if (k == 0) {
          rq[i] = NA_REAL;
        } else {
          r.resize(k); 
          quantile q = make_quantile(k, prob);
          vector<double> z(&r[0], &r[0+k]);
          sort(z.begin(), z.end());
          rq[i] = q.hlo * z[q.lo] + q.hhi * z[q.hi]; 
        }
      } 
      
    } else {
      // quantile setup
      quantile q = make_quantile(n, prob);
      for (i = p; i < L - n + p + 1; i = i + ss) {
        vector<double> z(&vec[i - p], &vec[i - p + n]);
        bool Ind = false;
        for (j = i - p; j < i - p + n; j++) {
          Ind = Ind|NumericVector::is_na(vec[j]);
        }
        if (Ind) {
          rq[i] = NA_REAL;
        } else {
          sort(z.begin(), z.end());
          rq[i] = q.hlo * z[q.lo] + q.hhi * z[q.hi];  
        }
      } 
    }
    if (sizeD) {
      int L1;
      L1 = (L - n + 1 + ss - 1)/ss; //This is a round-up method for a ratio of two integers: (L - n + 1)/ss
      NumericVector rq_sD(L1);
      for (i = p ; i < L - n + p + 1; i = i + ss){
        rq_sD[(i- p)/ss] = rq[i];
      }
      return rq_sD;
    } else {
      return rq;
    }
  }

//' @describeIn movQt An function equivalent to \code{movQt(..., align = "right")}
//' @export
// [[Rcpp::export]]
NumericVector 
  movQtr(NumericVector vec, int n = 1L, double prob = .5, int ss = 1L, 
        bool na_rm = false, bool sizeD = false) {
    
    // Declare loop counters, nonNA counter, vector sizes and sum
    int i, k, j;
    int L = vec.size();
    
    // Create vector filled with NA(R version)
    NumericVector rq(L, NumericVector::get_na());
    
    //Crux of the algorithm
    if (na_rm) {
      for (i = n - 1; i < L; i = i + ss) {
        k = 0;
        std::vector<double> r(n);
        for (j = i - n + 1; j < i + 1; j++) {
          if (NumericVector::is_na(vec[j])) {
          } else {
            r[k] = vec[j];
            k++;
          }
        }
        if (k == 0) {
          rq[i] = NA_REAL;
        } else {
          r.resize(k); 
          quantile q = make_quantile(k, prob);
          vector<double> z(&r[0], &r[0+k]);
          sort(z.begin(), z.end());
          rq[i] = q.hlo * z[q.lo] + q.hhi * z[q.hi]; 
        }
      } 
      
    } else {
      // quantile setup
      quantile q = make_quantile(n, prob);
      for (i = n - 1; i < L; i = i + ss) {
        vector<double> z(&vec[i - n + 1], &vec[i + 1]);
        bool Ind = false;
        for (j = i - n + 1; j < i + 1; j++) {
          Ind = Ind|NumericVector::is_na(vec[j]);
        }
        if (Ind) {
          rq[i] = NA_REAL;
        } else {
          sort(z.begin(), z.end());
          rq[i] = q.hlo * z[q.lo] + q.hhi * z[q.hi];  
        }
      } 
    }
    
    if (sizeD) {
      int L1;
      L1 = (L - n + 1 + ss - 1)/ss; //This is a round-up method for a ratio of two integers: (L - n + 1)/ss
      NumericVector rq_sD(L1);
      for (i = n - 1; i < L; i = i + ss){
        rq_sD[(i- n + 1)/ss] = rq[i];
      }
      return rq_sD;
    } else {
      return rq;
    }
  }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
movQt(vec = c(1, 4, 3, NA, 8, 4, 5, 9, 6, 0), n = 3, ss = 4, na_rm = TRUE, align = "right")
movQt(vec = c(1, 4, 3, NA, 8, 4, 5, 9, 6, 0), n = 3, na_rm = TRUE, align = "right")
  movQt(vec = c(1, 4, 3, NA, 8, 4, 5, 9, 6, 0), n = 3, ss = 4, na_rm = TRUE, sizeD = T)
  
  movQtr(vec = c(1, 4, 3, NA, 8, 4, 5, 9, 6, 0), n = 3, ss = 4, na_rm = TRUE)
  movQtr(vec = c(1, 4, 3, NA, 8, 4, 5, 9, 6, 0), n = 3, na_rm = TRUE)
  movQtr(vec = c(1, 4, 3, NA, 8, 4, 5, 9, 6, 0), n = 3, ss = 4, na_rm = TRUE, sizeD = T)
  
  
  movQt(c(1, 4, 3, NA, 8), 3, align = "left")
  movQt(c(1, 4, 3, NA, 8), 3, na_rm = TRUE, align = "right")

# movQt(c(1, 4, 3, NA, 8), 2, prob = .5, align = 'left')
# movQt(c(1, 4, 3, NA, 8), 3, prob = .5, na_rm = TRUE)
# movQt(c(1, 4, 3, NA, 8), 3, prob = .5, na_rm = TRUE, align = 'right')
# movQtr(c(1, 4, 3, NA, 8), 3, prob = .5, na_rm = TRUE)
# 
# 
# movQt(c(1, 5, 3, NA, 3), 2, prob = .5)
# movQt(c(1, 5, 3, NA, NA), 2, prob = .5, na_rm = TRUE)


*/