#' Ensemble forecasting from fitted forecasting objects
#'
#' @param fitted_objects
#' @param weights
#' @param level
#' @param ci
#'
#' @return
#' @export
#'
#' @examples
#'
#' x <- ts(rnorm(25))
#'
#' obj1 <- forecast::thetaf(x)
#' obj2 <- forecast::forecast(forecast::ets(x))
#'
#' list_obj <- list(obj1, obj2)
#'
#' # par(mfrow=c(2, 2))
#' # plot(ensemblef(list_obj, ci="E"))
#' # plot(ensemblef(list_obj, ci="A"))
#' # plot(ensemblef(list_obj, ci="T"))
#' # plot(ensemblef(list_obj, ci="garch"))
#'
ensemblef <- function(fitted_objects,
                      weights = rep(1 / length(fitted_objects),
                                    length(fitted_objects)),
                      ci = c("E", "A", "T",
                             "garch", "gaussian"))

{
  n_objects <- length(fitted_objects)
  stopifnot(n_objects == length(weights))
  if (!check_list_obj(fitted_objects))
    stop("objects not fitted on the same input time series")
  obj1 <- fitted_objects[[1]]
  x <- obj1$x
  level <- obj1$level
  nlevels_ <- length(level)
  start_preds <- start(obj1$mean)
  start_x <- start(x)
  freq <- frequency(obj1$x)
  if (!is.null(dim(obj1$lower)))
  {
    h <- length(obj1$lower[, 1])
  } else {
    h <- length(obj1$lower)
  }
  ci <- match.arg(ci)



  fcasts_matrix <- matrix(0, nrow = n_objects, ncol = h)
  resids_matrix <- matrix(0, nrow = n_objects, ncol = length(fitted_objects[[1]]$residuals))

  # container for the results
  #res_obj <- forecast::forecast(forecast::ets(x), h = h,
  #                                            level = level)
  res_obj <- vector("list", 0)
  res_obj$x <- x
  res_obj$model <- fitted_objects

  # indiv. forecasts
  `%op%` <- foreach::`%do%`
  pb <- txtProgressBar(min = 0, max = n_objects, style = 3)
  foreach::foreach (i = 1:n_objects)%op%{
    fcasts_matrix[i, ] <- as.numeric(fitted_objects[[i]]$mean)
    resids_matrix[i, ] <- as.numeric(fitted_objects[[i]]$residuals)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  fcasts_matrix <- as.matrix(fcasts_matrix)
  resids_matrix <- as.matrix(resids_matrix)

  # mean forecasts
  ans_mean <- ts(drop(crossprod(weights, fcasts_matrix)),
                 start = start_preds, frequency = freq)
  ans_residuals <- ts(drop(crossprod(weights, resids_matrix)),
                      start = start_x, frequency = freq)

  res_obj$residuals <- ans_residuals
  res_obj$method <- paste0(c(paste0(sapply(1:n_objects, function(i) fitted_objects[[i]]$method),
                                              collapse="+"), " (", ci, ")"), collapse = "")
  res_obj$ci <- ci

  if (ci == "E")
  {
    resid_fcast <-
      forecast::forecast(forecast::ets(ans_residuals), h = h)

    mean_residuals <- resid_fcast$mean
    res_obj$mean <- ans_mean + mean_residuals
    res_obj$lower <- ans_mean + resid_fcast$lower
    res_obj$upper <- ans_mean + resid_fcast$upper
  }

  if (ci == "A")
  {
    resid_fcast <-
      forecast::forecast(forecast::auto.arima(ans_residuals), h = h)
    res_obj$mean <- ans_mean + resid_fcast$mean
    res_obj$lower <- ans_mean + resid_fcast$lower
    res_obj$upper <- ans_mean + resid_fcast$upper
  }

  if (ci == "T")
  {
    resid_fcast <- forecast::thetaf(ans_residuals, h = h)
    res_obj$mean <- ans_mean + resid_fcast$mean
    res_obj$lower <- ans_mean + resid_fcast$lower
    res_obj$upper <- ans_mean + resid_fcast$upper
  }

  if (ci == "gaussian")
  {
    rep_1_h <- rep(1, h)
    conf_ints <- lapply(level,
                        function (x)
                          tcrossprod(rep_1_h, t.test(ans_residuals,
                                                     conf.level = x /
                                                       100)$conf))
    res_obj$mean <- ts(ans_mean + mean(ans_residuals),
                       start = start_preds)

    residuals_lower <- ts(sapply(1:nlevels_,
                        function(idx)
                        conf_ints[[idx]][, 1]),
                        start = start_preds)
    residuals_upper <- ts(sapply(1:nlevels_,
                                 function(idx)
                                   conf_ints[[idx]][, 2]),
                          start = start_preds)

    res_obj$lower <- ts(ans_mean + residuals_lower,
                        start = start_preds)
    res_obj$upper <- ts(ans_mean + residuals_upper,
                        start = start_preds)
  }

  if (ci == "garch")
  {
    resid_fcast <- after::garch11f(ans_residuals, h = h)
    res_obj$mean <- ans_mean + resid_fcast$mean
    res_obj$lower <- ans_mean + resid_fcast$lower
    res_obj$upper <- ans_mean + resid_fcast$upper
  }

  res_obj$weights <- weights

  res_obj$level <- level

  class(res_obj) <- "forecast"

  return (res_obj)
}


check_list_obj <- function(fitted_objects)
{
  x <- fitted_objects[[1]]$x
  for (i in 2:length(fitted_objects))
  {
    if (!all.equal(fitted_objects[[i]]$x, x))
      return (FALSE)
  }
  return(TRUE)
}
