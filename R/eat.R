#' combined ets-arima-theta forecasts
#'
#' combined ets, arima, and theta (eat) forecasting
#'
#' @export
#' @param y a univariate time series
#' @param h number of periods for forecasting
#' @param level confidence levels for prediction intervals
#' @param method forecasting method: "E": \code{forecast::ets};
#' "A": \code{forecast::auto.arima}; "T": \code{forecast::thetaf};
#'  or "EAT" for the combination of the three methods (with weights)
#' @param weights weights of each method, in method 'EAT'
#' @param ... additional parameters to be passed to \code{forecast::ets},
#' \code{forecast::auto.arima}, \code{forecast::thetaf} and
#' \code{forecast::forecast}
#'
#' @details ensemble forecasts obtained from \code{forecast::ets},
#' \code{forecast::auto.arima} and \code{forecast::theta}
#'
#' @return a list with point forecasts and prediction intervals; returns
#' from \code{forecast::ets}, \code{forecast::auto.arima},
#' \code{forecast::theta}, only point forecasts for eat
#'
#' \code{mean} point forecasts as a time series
#'
#' @export
#'
#' @examples
#'
#'require(forecast)
#'
#'print(after::eatf(WWWusage, method = "EAT",
#'weights = c(0.5, 0, 0.5)))
#'
#'print(after::eatf(WWWusage, method = "EAT"))
#'
#'print(after::eatf(WWWusage, method = "T"))
#'
#'
#'obj <- after::eatf(WWWusage, method = "EAT",
#'weights = c(0, 0.5, 0.5),
#'ci = "T")
#'plot(obj)
#'
#'
#'obj <- after::eatf(WWWusage, method = "EAT",
#'weights = c(0, 0.5, 0.5))
#'plot(obj)
#'
eatf <- function(y, h = 5,
                level = c(80, 95),
                method = c("E", "A", "T", "EAT"),
                weights = rep(1/3, 3),
                ci = c("E", "A", "T"),
                ...) {

  method <- match.arg(method)
  # stopifnot(length(level) == 1)
  ci <- match.arg(ci)
  n_y <- length(y)

  stopifnot(sum(weights) == 1)

  if (method == "E")
  {
    return(forecast::forecast(forecast::ets(y = y, ...),
                              h = h, level = level,...))
  }

  if (method == "A")
  {
    return(forecast::forecast(forecast::auto.arima(y = y, ...),
                              h = h, level = level, ...))
  }

  if (method == "T")
  {
    return(forecast::thetaf(y = y, h = h,
                            level = level, ...))
  }

  if (method == "EAT")
  {
    #stopifnot(sum(weights) == 1)

    if (all(weights != 0))
    {
      obj_ets <- forecast::forecast(forecast::ets(y = y, ...),
                                    h = h, level = level, ...)
      obj_arima <- forecast::forecast(forecast::auto.arima(y = y, ...),
                                      h = h, level = level, ...)
      obj_theta <- forecast::thetaf(y = y, h = h,
                                    level = level, ...)

      fcasts <- rbind(
        E = obj_ets$mean,
        A = obj_arima$mean,
        T = obj_theta$mean)

      resids <- rbind(
        E = obj_ets$residuals,
        A = obj_arima$residuals,
        T = obj_theta$residuals)

      res_obj <- switch(ci,
                        "E" = obj_ets,
                        "A" = obj_arima,
                        "T" = obj_theta)

      res_obj$mean <- ts(drop(crossprod(weights, fcasts)),
                         start = start(res_obj$mean),
                         frequency = frequency(y))
      res_obj$residuals <- ts(drop(crossprod(weights, resids)),
                              start = start(res_obj$mean),
                              frequency = frequency(y))
      res_obj$method <- paste0("EAT(", ci, ")")

    } else { # if (any(weights == 0))

      if (weights[1] == 0){ # model AT

        obj_arima <- forecast::forecast(forecast::auto.arima(y = y, ...),
                                        h = h, level = level, ...)
        obj_theta <- forecast::thetaf(y = y, h = h,
                                      level = level, ...)

        fcasts <- rbind(
          E = rep(0, h),
          A = obj_arima$mean,
          T = obj_theta$mean)

        resids <- rbind(
          E = rep(0, n_y),
          A = obj_arima$residuals,
          T = obj_theta$residuals)

        if((ci %in% c("A", "T") == FALSE))
           ci <- "A"

        res_obj <- switch(ci,
                          "A" = obj_arima,
                          "T" = obj_theta)

        res_obj$mean <- ts(drop(crossprod(weights, fcasts)),
                           start = start(res_obj$mean),
                           frequency = frequency(y))
        res_obj$residuals <- ts(drop(crossprod(weights, resids)),
                                start = start(res_obj$mean),
                                frequency = frequency(y))
        res_obj$method <- paste0("AT(", ci, ")")
      }

      if (weights[2] == 0){ # model ET

        obj_ets <- forecast::forecast(forecast::ets(y = y, ...),
                                      h = h, level = level, ...)
        obj_theta <- forecast::thetaf(y = y, h = h,
                                      level = level, ...)

        fcasts <- rbind(
          E = obj_ets$mean,
          A = rep(0, h),
          T = obj_theta$mean)

        resids <- rbind(
          E = obj_ets$residuals,
          A = rep(0, n_y),
          T = obj_theta$residuals)

        if((ci %in% c("E", "T") == FALSE))
          ci <- "E"

        res_obj <- switch(ci,
                          "E" = obj_ets,
                          "T" = obj_theta)

        res_obj$mean <- ts(drop(crossprod(weights, fcasts)),
                           start = start(res_obj$mean),
                           frequency = frequency(y))
        res_obj$residuals <- ts(drop(crossprod(weights, resids)),
                                start = start(res_obj$mean),
                                frequency = frequency(y))
        res_obj$method <- paste0("ET(", ci, ")")
      }

      if (weights[3] == 0){ # model EA

        obj_ets <- forecast::forecast(forecast::ets(y = y, ...),
                                      h = h, level = level, ...)
        obj_arima <- forecast::forecast(forecast::auto.arima(y = y, ...),
                                        h = h, level = level, ...)

        fcasts <- rbind(
          E = obj_ets$mean,
          A = obj_arima$mean,
          T = rep(0, h))

        resids <- rbind(
          E = obj_ets$residuals,
          A = obj_arima$residuals,
          T = rep(0, n_y))

        if((ci %in% c("E", "A") == FALSE))
          ci <- "E"

        res_obj <- switch(ci,
                          "E" = obj_ets,
                          "A" = obj_arima)

        res_obj$mean <- ts(drop(crossprod(weights, fcasts)),
                           start = start(res_obj$mean),
                           frequency = frequency(y))
        res_obj$residuals <- ts(drop(crossprod(weights, resids)),
                                start = start(res_obj$mean),
                                frequency = frequency(y))
        res_obj$method <- paste0("EA(", ci, ")")
      }

    }

    res_obj$ci <- ci

    res_obj$level <- level

    return(res_obj)

  }

}
