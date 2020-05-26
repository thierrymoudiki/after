

#' Title
#'
#' @param n_estimators
#' @param lags
#'
#' @return
#' @export
#'
#' @examples
#'
#' par(mfrow=c(3, 2))
#'
#' y <- Nile
#' (res <- after::smoothf(y, h=10))
#' plot(res)
#'
#' y <- AirPassengers
#' (res <- after::smoothf(y, h = 10))
#' plot(res)
#'
#' y <- USAccDeaths
#' (res <- after::smoothf(y, h = 10))
#' plot(res)
#'
#' y <- WWWusage
#' (res <- after::smoothf(y, h = 10))
#' plot(res)
#'
#' y <- WWWusage
#' (res <- after::smoothf(y, h = 10, lags = 2L))
#' plot(res)
#'
#'
#' y <- WWWusage
#' (res <- after::smoothf(y, h = 10, lags = 3L))
#' plot(res)
#'
#'
#' grid_params <- expand.grid(n_estimators=1:10, lags=c(1, 2, 3, 4))
#'
#'
smoothf <- function(y, h=5,
                    level=c(80, 95),
                    n_estimators=5, lags=1,
                    method=c("mean", "median"),
                    ci = c("E", "A", "T", "garch"))
{
  method <- match.arg(method)
  ci <- match.arg(ci)
  freq_y <- frequency(y)
  start_preds <- tsp(y)[2] + 1/freq_y



  # Training set fitting --------------------------------------------------------

  Y <- after::embedc(x = y, lags = lags)
  obs <- Y[, 1]


  W_ <- randtoolbox::sobol(n = n_estimators,
                           dim = lags + 1)
  W <- unique(W_/rowSums(W_))
  fitted <- switch(method,
                   "mean" = rowMeans(tcrossprod(Y, W)),
                   "median" = apply(tcrossprod(Y, W), 1,
                                    median))


  res_obj <- vector("list", 0)
  res_obj$fitted <- ts(fitted, end = end(y), frequency = freq_y)
  res_obj$x <- ts(obs, end = end(y), frequency = freq_y)
  res_obj$residuals <- ts(obs - fitted, end = end(y),
                      frequency = freq_y)



  # Forecasting loop --------------------------------------------------------

  y <- c(y, tail(fitted, 1))
  for (i in 1:(h-1))
  {
    Y <- after::embedc(x = y, lags = lags)
    fitted <- switch(method,
                     "mean" = rowMeans(tcrossprod(Y, W)),
                     "median" = apply(tcrossprod(Y, W), 1,
                                      median))
    y <- c(y, tail(fitted, 1))
  }


  ans_mean <- ts(tail(y, h), start = start_preds,
                 frequency = freq_y)



  # Forecasting residuals --------------------------------------------------------

  if (ci == "E")
  {
    resid_fcast <-
      forecast::forecast(forecast::ets(res_obj$residuals),
                         level=level, h = h)

    mean_residuals <- resid_fcast$mean
    res_obj$mean <- ans_mean + mean_residuals
    res_obj$lower <- ans_mean + resid_fcast$lower
    res_obj$upper <- ans_mean + resid_fcast$upper
  }

  if (ci == "A")
  {
    resid_fcast <-
      forecast::forecast(forecast::auto.arima(res_obj$residuals),
                         level=level, h = h)
    res_obj$mean <- ans_mean + resid_fcast$mean
    res_obj$lower <- ans_mean + resid_fcast$lower
    res_obj$upper <- ans_mean + resid_fcast$upper
  }

  if (ci == "T")
  {
    resid_fcast <- forecast::thetaf(res_obj$residuals,
                                    level=level, h = h)
    res_obj$mean <- ans_mean + resid_fcast$mean
    res_obj$lower <- ans_mean + resid_fcast$lower
    res_obj$upper <- ans_mean + resid_fcast$upper
  }

  if (ci == "garch")
  {
    resid_fcast <- after::garch11f(res_obj$residuals,
                                   level=level, h = h)
    res_obj$mean <- ans_mean + resid_fcast$mean
    res_obj$lower <- ans_mean + resid_fcast$lower
    res_obj$upper <- ans_mean + resid_fcast$upper
  }



  # Return --------------------------------------------------------

  res_obj$model <- list(n_estimators=n_estimators,
                        lags=lags)

  res_obj$level <- level

  res_obj$method <- "smoothf"

  n <- length(res_obj$x)

  loglik <- n*log(mean(res_obj$residuals^2))

  k <- (n_estimators + lags)

  res_obj$aic <- loglik + 2*k

  res_obj$bic <- loglik + log(n)*k

  class(res_obj) <- "forecast"

  return(res_obj)
}
