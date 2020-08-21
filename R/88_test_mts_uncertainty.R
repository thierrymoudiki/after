

x <- fpp::usconsumption[,1]
mean_forecast <- forecast::thetaf(x)$mean
resids <- forecast::thetaf(x)$resid
(res <- get_uncertainty(mean_forecast, resids, h=5, level=95, x))
print(str(res))

get_uncertainty <- function(mean_forecast, resids,
                            h, level, x,
                            ci=c("E", "A",
                                 "T", "garch"))
{
  ci <- match.arg(ci)
  res_obj <- vector("list", 0)
  res_obj$residuals <- resids

  if (ci == "E")
  {
    resid_fcast <-
      forecast::forecast(forecast::ets(res_obj$residuals),
                         level=level, h = h)

    mean_residuals <- resid_fcast$mean
    res_obj$mean <- mean_forecast + mean_residuals
    res_obj$lower <- mean_forecast + resid_fcast$lower
    res_obj$upper <- mean_forecast + resid_fcast$upper
  }

  if (ci == "A")
  {
    resid_fcast <-
      forecast::forecast(forecast::auto.arima(res_obj$residuals),
                         level=level, h = h)
    res_obj$mean <- mean_forecast + resid_fcast$mean
    res_obj$lower <- mean_forecast + resid_fcast$lower
    res_obj$upper <- mean_forecast + resid_fcast$upper
  }

  if (ci == "T")
  {
    resid_fcast <- forecast::thetaf(res_obj$residuals,
                                    level=level, h = h)
    res_obj$mean <- mean_forecast + resid_fcast$mean
    res_obj$lower <- mean_forecast + resid_fcast$lower
    res_obj$upper <- mean_forecast + resid_fcast$upper
  }

  if (ci == "garch")
  {
    resid_fcast <- garch11f(res_obj$residuals,
                            level=level, h = h)
    res_obj$mean <- mean_forecast + resid_fcast$mean
    res_obj$lower <- mean_forecast + resid_fcast$lower
    res_obj$upper <- mean_forecast + resid_fcast$upper
  }

  res_obj$fitted <- 0
  res_obj$x <- x
  res_obj$model <- "None"
  res_obj$level <- level
  res_obj$method <- "None"
  class(res_obj) <- "forecast"

  return(res_obj)
}

