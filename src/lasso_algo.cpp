#include <Rcpp.h>
#include <Math.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sign_cpp(double x)
{
  return (x < 0)?-1:1;
}

// [[Rcpp::export]]
double abs_cpp(double x)
{
  return (x < 0)?-x:x;
}

// [[Rcpp::export]]
double max_cpp(double x, double y)
{
  return (x < y)?y:x;
}

// [[Rcpp::export]]
double norm_cpp(NumericVector x)
{
  unsigned long int n = x.size();
  double res = 0;

    for (unsigned long int i = 0; i < n; i++)
    {
      res += abs_cpp(x(i));
    }

  return res;
}

// [[Rcpp::export]]
double cross_prod_cpp(NumericVector x, NumericVector y)
{
  unsigned long int n = x.size();
    if (y.size() != n) {
      ::Rf_error("you must have x.size() == y.size()");
    }
  double res = 0;

      for (unsigned long int i = 0; i < n; i++)
      {
        res += x(i)*y(i);
      }

      // Rcpp::Rcout << "res is " << std::endl << res << std::endl;
  return res;
}

// [[Rcpp::export]]
double soft_thres_cpp(double x, double y) {
  return sign_cpp(x)*max_cpp(abs_cpp(x) - y, 0);
}

// [[Rcpp::export]]
NumericVector pass_by_value(NumericVector x, NumericVector y) {
  unsigned long int n = x.size();
  for (unsigned long int i = 0; i < n; i++)
  {
    x(i) = y(i);
  }
  return x;
}

// #Lasso via coordinate descent: the "Shooting Algorithm" of Fu (1998). Adapted from FDiTraglia's Github code +
// pseudocode algorithm 13.1 of Murphy (2012) + matlab code LassoShooting.m by Mark Schmidt.
// [[Rcpp::export]]
List lasso_shoot_cpp(NumericVector beta, NumericMatrix XX2, NumericVector Xy2,
                     double lambda, double tol = 1e-05, unsigned long int max_iter = 1e05){

  unsigned long int p = XX2.nrow();
  bool converged = false;
  unsigned long int iteration = 0;
  NumericVector beta_opt(p), beta_prev(p);
  double aj, cj;
  double error;
      // pass_by_value beta_opt with beta to avoid modifying
      // directly beta into the function (find out how to this differently)
      beta_opt = pass_by_value(beta_opt, beta);

    while ((converged == false) && (iteration < max_iter)){

      beta_prev = pass_by_value(beta_prev, beta_opt);

        for (unsigned long int j = 0; j < p; j++){
          aj = XX2(j,j);
          cj = Xy2(j) - cross_prod_cpp(XX2(j, _), beta_opt) + beta_opt(j) * aj;
          beta_opt(j) = soft_thres_cpp(cj/aj, lambda/aj);
           // Rcpp::Rcout << "beta is " << std::endl << beta_opt(j) << std::endl;
        }

      error = norm_cpp(beta_prev - beta_opt);
      iteration = iteration + 1;
      converged =  (error < tol);

    }

  return Rcpp::List::create(Rcpp::Named("beta") = beta_opt,
                    Rcpp::Named("n_iter") = iteration,
                    Rcpp::Named("converged") = converged,
                    Rcpp::Named("error") = error);
}

/*** R
# library(regularizedmts)
# fit_obj <- fit_var_mts(housing, penalization = "l1", lags = 2)
# lasso_shoot_cpp(beta = fit_obj$beta0[[1]], XX2= fit_obj$XX2,
#                 Xy2 = fit_obj$Xy2[[1]], lambda = 0.1,
#                 tol = 1e-05, max_iter = 1e05)
*/
