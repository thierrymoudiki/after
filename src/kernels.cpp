#include <Rcpp.h>
#include <Math.h>
using namespace Rcpp;

/* utils */

// [[Rcpp::export]]
double sum_cpp(NumericVector x)
{
  unsigned long int n = x.size();
  double res = 0;
  for(int i = 0; i < n; i++) {
    res += x(i);
  }
  return(res);
}

// [[Rcpp::export]]
double l2_norm(NumericVector x)
{
  unsigned long int n = x.size();
  double res = 0;
  for(int i = 0; i < n; i++) {
    res += pow(x(i), 2);
  }
  return(sqrt(res));
}

/* Gaussian kernels */

// [[Rcpp::export]]
NumericMatrix gaussian_kxx_cpp(NumericMatrix x, double sigma, double l)
{
  unsigned long int n = x.nrow();
  NumericMatrix res(n, n);
  double sigma_sq = pow(sigma, 2);

  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      res(i , j) = exp(-0.5*pow(l2_norm(x(i, _) - x(j, _)), 2)/sigma_sq);
      res(j , i) = res(i , j);
    }
  }
  return(pow(l, 2)*res);
}

// [[Rcpp::export]]
NumericVector gaussian_kxy_cpp(NumericMatrix x, NumericVector y, double sigma, double l)
{
  unsigned long int m = x.ncol();
  unsigned long int n = y.size();
  if (m != n) {
    ::Rf_error("you must have x.col() == y.size()");
  }
  unsigned long int k = x.nrow();

  NumericVector res(k);
  double sigma_sq = pow(sigma, 2);

  for(int i = 0; i < k; i++) {
    res(i) = exp(-0.5*pow(l2_norm(x(i, _) - y), 2)/sigma_sq);
  }
  return(pow(l, 2)*res);
}

/* Polynomial kernels */

// [[Rcpp::export]]
NumericMatrix poly_kxx_cpp(NumericMatrix x, double sigma, int d, double l)
{
  unsigned long int n = x.nrow();
  NumericMatrix res(n, n);
  double sigma_sq = pow(sigma, 2);

  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      res(i , j) = pow((1 + sum_cpp(x(i, _)*x(j, _))/sigma_sq), d);
      res(j , i) = res(i , j);
    }
  }
  return(pow(l, 2)*res);
}

// [[Rcpp::export]]
NumericVector poly_kxy_cpp(NumericMatrix x, NumericVector y, double sigma, int d, double l)
{
  unsigned long int m = x.ncol();
  unsigned long int n = y.size();
  if (m != n) {
    ::Rf_error("you must have x.col() == y.size()");
  }
  unsigned long int k = x.nrow();

  NumericVector res(k);
  double sigma_sq = pow(sigma, 2);

  for(int i = 0; i < k; i++) {
    res(i) = pow((1 + sum_cpp(x(i, _)*y)/sigma_sq), d);
  }
  return(pow(l, 2)*res);
}

/* Matern kernels 3/2 */ /* GPML P.103 */

// [[Rcpp::export]]
NumericMatrix matern32_kxx_cpp(NumericMatrix x, double sigma, double l)
{
  unsigned long int n = x.nrow();
  NumericMatrix res(n, n);
  NumericVector r(n);
  double sqrt3 = sqrt(3);
  double temp = 0;

  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      temp = sqrt3*l2_norm(x(i, _) - x(j, _))/sigma;
      res(i , j) = (1 + temp)*exp(-temp);
      res(j , i) = res(i , j);
    }
  }
  return(pow(l, 2)*res);
}

// [[Rcpp::export]]
NumericVector matern32_kxy_cpp(NumericMatrix x, NumericVector y, double sigma, double l)
{
  unsigned long int m = x.ncol();
  unsigned long int n = y.size();
  if (m != n) {
    ::Rf_error("you must have x.col() == y.size()");
  }
  unsigned long int k = x.nrow();

  NumericVector res(k);
  double sqrt3 = sqrt(3);
  double temp = 0;

  for(int i = 0; i < k; i++) {
    temp = sqrt3*l2_norm(x(i, _) - y)/sigma;
    res(i) = (1 + temp)*exp(-temp);
  }
  return(pow(l, 2)*res);
}

/* Matern kernels 5/2 */ /* GPML P.103 */

// [[Rcpp::export]]
NumericMatrix matern52_kxx_cpp(NumericMatrix x, double sigma, double l)
{
  unsigned long int n = x.nrow();
  NumericMatrix res(n, n);
  NumericVector r(n);
  double sqrt5 = sqrt(5);
  double temp = 0;

  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      r = x(i, _) - x(j, _);
      temp = sqrt5*l2_norm(r)/sigma;
      res(i , j) = (1 + temp + pow(temp, 2)/3)*exp(-temp);
      res(j , i) = res(i , j);
    }
  }
  return(pow(l, 2)*res);
}

// [[Rcpp::export]]
NumericVector matern52_kxy_cpp(NumericMatrix x, NumericVector y, double sigma, double l)
{
  unsigned long int m = x.ncol();
  unsigned long int n = y.size();
  if (m != n) {
    ::Rf_error("you must have x.col() == y.size()");
  }
  unsigned long int k = x.nrow();

  NumericVector res(k);
  NumericVector r(k);
  double temp;
  double sqrt5 = sqrt(5);

  for(int i = 0; i < k; i++) {
    r = x(i, _) - y;
    temp = sqrt5*l2_norm(r)/sigma;
    res(i) = (1 + temp + pow(temp, 2)/3)*exp(-temp);
  }
  return(pow(l, 2)*res);
}

