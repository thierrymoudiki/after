% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/99_forecast_mts_interface.R
\name{krlsf}
\alias{krlsf}
\title{Title}
\usage{
krlsf(
  x,
  lags = 1,
  lambda_krls = 0.1,
  l = 0.5,
  sigma = 10,
  d = 1,
  kernel_type = c("gaussian", "polynomial", "matern32", "matern52"),
  inv_method = c("chol", "ginv"),
  h = 5,
  type_ci = "none",
  type_forecast = c("recursive", "direct"),
  level = 95
)
}
\arguments{
\item{level}{}
}
\description{
Title
}
