% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/99_forecast_after_interface.R
\name{arimaf}
\alias{arimaf}
\title{Title}
\usage{
arimaf(
  y,
  h = 5,
  level = c(80, 95),
  ci = c("gaussian", "E", "A", "T", "garch"),
  ...
)
}
\arguments{
\item{y}{lorem ipsum}

\item{h}{lorem ipsum}

\item{...}{additional parameters to be passed to \code{stats::arima0}}
}
\description{
Title
}
\examples{

try(arimaf(AirPassengers, h=7), silent=TRUE)

try(arimaf(Nile, h=6), silent=TRUE)

try(arimaf(WWWusage, h=8), silent=TRUE)

try(arimaf(USAccDeaths, h=10), silent=TRUE)

}
