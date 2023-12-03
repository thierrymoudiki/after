#include <Rcpp.h>
using namespace Rcpp;

// ----- 0 - utils

// calculate the crossproduct of x and y
// [[Rcpp::export]]
double crossprod_cpp(NumericVector x, NumericVector y)
{
  unsigned long int n = x.size();
  if (y.size() != n) {
    ::Rf_error("both input vectors must have the same length");
  }
  double res = 0; // variable containing the result

  for(int i = 0; i < n; i++) {
    res += x(i)*y(i);
  }
  return(res);
}

// calculate the crossproduct of eL's columns
// [[Rcpp::export]]
NumericVector columns_crossprod_cpp(NumericMatrix eL)
{
  unsigned long int m = eL.ncol();
  NumericVector res(m); // variable containing the result
  NumericVector temp(eL.nrow());

  for(long int i = 0; i < m; i++) {
    temp = eL(_, i);
    res(i) = crossprod_cpp(temp, temp);
  }
  return(res);
}

// calculate the squared crossproduct of eL's columns and hL
// [[Rcpp::export]]
NumericVector squared_crossprod_cpp(NumericMatrix eL, NumericVector hL)
{
  unsigned long int m = eL.ncol();
  unsigned long int N = eL.nrow();
  if (hL.size() != N) {
    ::Rf_error("both input vectors must have the same length");
  }
  NumericVector res(m); // variable containing the result

  for(long int i = 0; i < m; i++) {
    res(i) = pow(crossprod_cpp(eL(_, i), hL), 2);
  }

  return(res);
}

// ----- 1 - algo's elements

// compute the regressor at step L
// only with bounded activation functions (sigmoid and tanh here)
// [[Rcpp::export]]
NumericVector calculate_hL(NumericMatrix x, NumericVector w, Rcpp::String activation)
{
  unsigned long int N = x.nrow();
  NumericVector res(N); // variable containing the result

    if (w.size() != x.ncol()) {
      ::Rf_error("incompatible dimensions: requires x.ncol() == w.size()");
    }

    if (activation == "tanh")  {
      for(long int i = 0; i < N; i++) {
        res(i) = std::tanh(crossprod_cpp(x(i, _), w));
      }
    }

    if (activation == "sigmoid") {
      for(long int i = 0; i < N; i++) {
        res(i) = 1/(1 + std::exp(-(crossprod_cpp(x(i, _), w))));
      }
    }

  return(res);
}

// calculate xsi, that serve for determining the condition of convergence
// [[Rcpp::export]]
NumericVector calculate_xsiL(NumericMatrix eL, NumericVector hL, double nu,
                             double r, unsigned long int L)
{
  return(nu*(2-nu)*squared_crossprod_cpp(eL, hL)/crossprod_cpp(hL, hL) - (1 - r - (1 - r)/(L + 1))*columns_crossprod_cpp(eL));
}

// regression of current error eL on hL => obtain the betas
// [[Rcpp::export]]
NumericVector calculate_betasL(NumericMatrix eL, NumericVector hL)
{
  unsigned long int m = eL.ncol();
  if (hL.size() != eL.nrow()) {
    ::Rf_error("incompatible dimensions: requires hL.size() == eL.nrow()");
  }
  NumericVector res(m); // variable containing the result

  for (long int i = 0; i < m; i++)
  {
    res(i) = crossprod_cpp(eL(_, i), hL);
  }
  return(res/crossprod_cpp(hL, hL));
}

// calculate fitted values at step L, with learnung rate = nu
// [[Rcpp::export]]
NumericMatrix calculate_fittedeL(NumericVector betasL, NumericVector hL, double nu)
{
  unsigned long int m = betasL.size();
  unsigned long int N = hL.size();
  NumericMatrix res(N, m); // variable containing the result
  for (long int i = 0; i < N; i++){
    for (long int j = 0; j < m; j++)
    {
      res(i, j) = nu*betasL(j)*hL(i);
    }
  }
  return(res);
}

/*** R
# set.seed(156)
# n <- 20 ; p <- 5
# X <- matrix(rnorm(n * p), n, p) # no intercept!
# y <- matrix(rnorm(4*n), ncol = 4)
#
# fit_obj_scn <- fit_SCN(x = X, y = y, B = 500, lam = 100,
#                                        nu = 0.5, col_sample = 0.8)
#  preds <- predict_SCN(fit_obj_scn, newx = X)
#
# par(mfrow=c(2, 2))
# plot(1:nrow(X), y[,1], type = 'l')
#  lines(1:nrow(X), preds[,1],
#       col = "red")
# plot(1:nrow(X), y[,2], type = 'l')
# lines(1:nrow(X), preds[,2],
#       col = "red")
# plot(1:nrow(X), y[,3], type = 'l')
# lines(1:nrow(X), preds[,3],
#       col = "red")
# plot(1:nrow(X), y[,4], type = 'l')
# lines(1:nrow(X), preds[,4],
#       col = "red")

*/
