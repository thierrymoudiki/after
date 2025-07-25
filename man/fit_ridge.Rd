% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_lm_ridge.R
\name{fit_ridge}
\alias{fit_ridge}
\title{Title}
\usage{
fit_ridge(x, y, lambda = 10^seq(-10, 10, length.out = 100))
}
\arguments{
\item{lambda}{}
}
\description{
Title
}
\examples{

n <- 100 ; p <- 10
X <- matrix(rnorm(n * p), n, p) # no intercept!
y <- rnorm(n)

(fit_obj <- after::fit_ridge(X, y))

matplot(fit_obj$lambda, t(fit_obj$coef), type = 'l')

}
