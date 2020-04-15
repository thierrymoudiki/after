

#' GARCH(1, 1) forecasting function
#'
#' @param x a univariate time series
#' @param h number of periods for forecasting
#' @param level confidence levels for prediction intervals
#'
#' @return
#' @export
#'
#' @examples
#'
#' x <- ts(fGarch::garchSim()$garch)
#' # (obj <- after::garch11f(x))
#'
garch11f <- function(x, h=5, level=c(80, 95))
{

  # container for the results
  fit <- fGarch::garchFit(~ garch(1,1), data = x, trace = FALSE)
  preds <- predict(fit, h)
  critical_values <- sapply(level, function(q) qnorm(1 - (1 - q/100)/2))
  freq_x <- frequency(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  res_obj <- vector("list", 0)
  res_obj$model <- fit
  res_obj$x <- x
  res_obj$residuals <- residuals(fit)
  res_obj$fitted <- fitted(fit)
  res_obj$level <- level
  res_obj$mean <- ts(preds$meanForecast, start = start_preds, frequency = freq_x)
  res_obj$upper <- ts(sapply(critical_values,
                             function(q) preds$meanForecast + q*preds$standardDeviation),
                      start=start_preds, frequency = freq_x)
  res_obj$lower <- ts(sapply(critical_values,
                             function(q) preds$meanForecast - q*preds$standardDeviation),
                      start=start_preds, frequency = freq_x)
  res_obj$method <- "garch(1, 1)"

  class(res_obj) <- "forecast"

  return (res_obj)
}
