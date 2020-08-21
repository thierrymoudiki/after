# 0 - Import functions ---------------------------------------------------

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir("afterfuncs")

# 1 - Individual models ---------------------------------------------------

#' Title
#'
#' @param y
#' @param h
#' @param ... additional parameters to be passed to `stats::arima0`
#'
#' @return
#' @export
#'
#' @examples
#'
#' try(arimaf(AirPassengers, h=7), silent=TRUE)
#'
#' try(arimaf(Nile, h=6), silent=TRUE)
#'
#' try(arimaf(WWWusage, h=8), silent=TRUE)
#'
#' try(arimaf(USAccDeaths, h=10), silent=TRUE)
#'
arimaf <- function(y, h=5, level=c(80, 95),
                   ci = c("gaussian", "E", "A", "T",
                          "garch"),
                   ...)
{
  n <- length(y)
  ci <- match.arg(ci)
  arima_combns_grid <- expand.grid(c(0, 1, 2),
                 c(0, 1, 2),
                 c(0, 1, 2))

  objective <- function(xx, ...)
  {
    k <- sum(xx)
    res <- try(stats::arima0(x,
                             order = c(xx[1], xx[2], xx[3]),
                             ...)$aic + (2*k*(k+1))/(n-k-1),
               silent = TRUE)
    if (class(res) == "try-error")
      return(.Machine$integer.max)
    return(res)
  }

  df <- cbind.data.frame(arima_combns_grid,
                         sapply(1:nrow(arima_combns_grid),
                                function (i) objective(as.numeric(arima_combns_grid[i, ]))))

  order_ <- as.numeric(df[which.min(df[, 4]), ][, c(1, 2, 3)])

  fit <- stats::arima0(y, order = order_, ...)
  out <- list()
  out$x <- y
  out$residuals <- residuals(fit)
  out$model <- fit
  out$method <- paste0("arima(", order_[1], ",",
                       order_[2], ",", order_[3], ")")
  print(fit)
  preds <- predict(fit, n.ahead = h)
  fcast  <- preds$pred
  out$sigma <- preds$se

  # Compute prediction intervals

  if (ci == "E")
  {
    resid_fcast <- forecast::forecast(forecast::ets(out$residuals),
                                      h = h)
    out$mean <- fcast + resid_fcast$mean
    out$lower <- fcast + resid_fcast$lower
    out$upper <- fcast + resid_fcast$upper
  }

  if (ci == "T")
  {
    resid_fcast <- forecast::thetaf(out$residuals,
                                    h = h)
    out$mean <- fcast + resid_fcast$mean
    out$lower <- fcast + resid_fcast$lower
    out$upper <- fcast + resid_fcast$upper
  }

  if (ci == "garch")
  {
    resid_fcast <- garch11f(out$residuals,
                                   h = h)
    out$mean <- fcast + resid_fcast$mean
    out$lower <- fcast + resid_fcast$lower
    out$upper <- fcast + resid_fcast$upper
  }

  if (ci == "gaussian")
  {
    qts <- sapply(level, function (x) qnorm(1-(1-(100-x)/200)))
    nlevels_ <- length(qts)
    out$mean <- fcast
    out$lower <- ts(matrix(nrow = h, ncol = nlevels_),
                    start = start(out$mean),
                    frequency = frequency(out$mean))
    colnames(out$lower) <- paste("Lo", level)
    out$upper <- ts(matrix(nrow = h, ncol = nlevels_),
                    start = start(out$mean),
                    frequency = frequency(out$mean))
    colnames(out$upper) <- paste("Hi", level)
    for (i in 1:nlevels_)
    {
      out$lower[, i] <-  fcast - qts[i]*out$sigma
      out$upper[, i] <-  fcast + qts[i]*out$sigma
    }
  }

  out$level <- level

  return(structure(out, class = "forecast"))

}


#' Title
#'
#' @param y
#' @param h
#' @param level
#' @param fit_func
#' @param predict_func
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' par(mfrow=c(3, 2))
#' plot(dynrmf(USAccDeaths, level=c(80, 90, 95, 99)))
#' plot(dynrmf(AirPassengers, level=c(80, 90, 95, 99)))
#' plot(dynrmf(lynx, level=c(80, 90, 95, 99)))
#' plot(dynrmf(WWWusage, level=c(80, 90, 95, 99)))
#' plot(dynrmf(Nile, level=c(80, 90, 95, 99)))
#' plot(dynrmf(fdeaths, level=c(80, 90, 95, 99)))
#'
#'
dynrmf <- function(y, h = 5,
                   level = c(80, 95),
                   fit_func = fit_ridge,
                   predict_func = predict_ridge,
                   ...)
{
  obj <- dynrm_fit(y,
                          fit_func = fit_func,
                          predict_func = predict_func,
                          ...)

  return(dynrm_predict(obj, h=h, level=level, ...))
}


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
#'print(eatf(WWWusage, method = "EAT",
#'weights = c(0.5, 0, 0.5)))
#'
#'print(eatf(WWWusage, method = "EAT"))
#'
#'print(eatf(WWWusage, method = "T"))
#'
#'
#'obj <- eatf(WWWusage, method = "EAT",
#'weights = c(0, 0.5, 0.5),
#'ci = "T")
#'plot(obj)
#'
#'
#'obj <- eatf(WWWusage, method = "EAT",
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
  resids_matrix <- matrix(0, nrow = n_objects,
                          ncol = length(fitted_objects[[1]]$residuals))

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
  res_obj$method <- paste0(c(paste0(sapply(1:n_objects,
                                           function(i) fitted_objects[[i]]$method),
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
    qts <- sapply(level, function (x) qnorm(1-(1-(100-x)/200)))
    nlevels_ <- length(qts)
    res_obj$mean <- ans_mean
    res_obj$lower <- ts(matrix(nrow = h, ncol = nlevels_),
                    start = start(res_obj$mean),
                    frequency = frequency(res_obj$mean))
    colnames(res_obj$lower) <- paste("Lo", level)
    res_obj$upper <- ts(matrix(nrow = h, ncol = nlevels_),
                    start = start(res_obj$mean),
                    frequency = frequency(res_obj$mean))
    colnames(res_obj$upper) <- paste("Hi", level)
    for (i in 1:nlevels_)
    {
      res_obj$lower[, i] <-  fcast - qts[i]*res_obj$sigma
      res_obj$upper[, i] <-  fcast + qts[i]*res_obj$sigma
    }
  }



  if (ci == "garch")
  {
    resid_fcast <- garch11f(ans_residuals, h = h)
    res_obj$mean <- ans_mean + resid_fcast$mean
    res_obj$lower <- ans_mean + resid_fcast$lower
    res_obj$upper <- ans_mean + resid_fcast$upper
  }

  res_obj$weights <- weights

  res_obj$level <- level

  class(res_obj) <- "forecast"

  return (res_obj)
}


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
#' # (obj <- garch11f(x))
#'
garch11f <- function(x, h=5, level=c(80, 95))
{

  # container for the results
  fit <- fGarch::garchFit(~ garch(1,1), data = x, trace = FALSE)
  preds <- fGarch::predict(fit, h)
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


#' Title
#'
#' @param y
#' @param h
#' @param level
#' @param method
#' @param fan
#' @param lambda
#' @param biasadj
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' par(mfrow=c(2, 2))
#' plot(mmf(Nile))
#' plot(mmf(Nile, method = "median"))
#' plot(mmf(AirPassengers, method = "median"))
#' plot(mmf(AirPassengers))
#'
mmf <- function(y, h = 5,
                level = c(80, 95),
                method=c("mean", "median"),
                fan = FALSE,
                lambda = NULL,
                biasadj = FALSE,
                ...)
{
  method <- match.arg(method)

  x <- diff(y)

  # adapted from forecast::meanf
  n <- length(x)
  if (!is.null(lambda)) {
    origx <- x
    x <- BoxCox(x, lambda)
    lambda <- attr(x, "lambda")
  }

  if (method == "mean")
  {
    meanx <- mean(x, na.rm = TRUE)
  } else {
    meanx <- median(x, na.rm = TRUE)
  }

  fits <- cumsum(c(y[1], rep(meanx, length(x))))
  res <- y - fits
  f <- as.numeric(tail(y, 1)) + cumsum(rep(meanx, h))

  if (fan) {
    level <- seq(51, 99, by = 3)
  }
  else {
    if (min(level) > 0 && max(level) < 1) {
      level <- 100 * level
    }
    else if (min(level) < 0 || max(level) > 99.99) {
      stop("Confidence limit out of range")
    }
  }
  nconf <- length(level)
  s <- sd(x, na.rm = TRUE)

  lower <- upper <- matrix(NA, nrow = h, ncol = nconf)
  for (i in 1:nconf) {
    if (n > 1) {
      tfrac <- qt(0.5 - level[i]/200, n - 1)
    }
    else {
      tfrac <- -Inf
    }
    w <- -tfrac * s * sqrt(1 + 1/n)
    lower[, i] <- f - w
    upper[, i] <- f + w
  }

  colnames(lower) <- colnames(upper) <- paste(level, "%", sep = "")

  if (is.ts(x)) {
    fits <- ts(fits)
    res <- ts(res)
    tsp(fits) <- tsp(res) <- tsp(y)
    freq <- frequency(y)
    f <- ts(f, start = tsp(y)[2] + 1/freq, frequency = freq)
    lower <- ts(lower, start = tsp(y)[2] + 1/freq, frequency = freq)
    upper <- ts(upper, start = tsp(y)[2] + 1/freq, frequency = freq)
  }

  if (!is.null(lambda)) {
    fits <- InvBoxCox(fits, lambda)
    x <- origx
    f <- InvBoxCox(f, lambda, biasadj,
                   list(level = level,
                        upper = upper, lower = lower))
    lower <- InvBoxCox(lower, lambda)
    upper <- InvBoxCox(upper, lambda)
  }

  out <- list(method = "mmf",
              level = level,
              x = y,
              series = deparse(substitute(y)),
              mean = f,
              lower = lower, upper = upper,
              model = structure(list(mu = f[1],
                                     mu.se = s/sqrt(length(x)),
                                     sd = s),
                                class = "meanf"),
              lambda = lambda,
              fitted = fits,
              residuals = res)

  out$model$call <- match.call()
  return(structure(out, class = "forecast"))

}


#' Title
#'
#' @param x
#' @param h
#' @param degree
#' @param lambda
#' @param level
#' @param fan
#'
#' @return
#' @export
#'
#' @examples
polythetaf <- function (y, h = ifelse(frequency(x) > 1, 2 * frequency(x), 10),
                        lambda=10^seq(from=-5, to=4,
                                      length.out = 100),
                        level = c(80, 95), fan = FALSE)
{
  if (fan) {
    level <- seq(51, 99, by = 3)
  }
  else {
    if (min(level) > 0 && max(level) < 1) {
      level <- 100 * level
    }
    else if (min(level) < 0 || max(level) > 99.99) {
      stop("Confidence limit out of range")
    }
  }

  x <- ts(diff(y), end = end(y), frequency = frequency(y))
  n <- length(x)
  x <- as.ts(x)
  m <- frequency(x)
  if (m > 1 && !is.constant(x) && n > 2 * m) {
    r <- as.numeric(acf(x, lag.max = m, plot = FALSE)$acf)[-1]
    stat <- sqrt((1 + 2 * sum(r[-m]^2))/n)
    seasonal <- (abs(r[m])/stat > qnorm(0.95))
  }
  else {
    seasonal <- FALSE
  }

  origx <- x
  if (seasonal) {
    decomp <- decompose(x, type = "multiplicative")
    if (any(abs(seasonal(decomp)) < 1e-10)) {
      warning("Seasonal indexes equal to zero. Using non-seasonal Theta method")
    }
    else {
      x <- seasadj(decomp)
    }
  }

  fcast <- ses(x, h = h)
  training_trend <- 0:(n - 1)
  test_trend <- 0:(h - 1)

  training_data <- cbind.data.frame(y = x, trend = training_trend)
  test_data <- as.matrix(test_trend)

  regr <- polytheta::fit_ridge(x=training_data[,-1], y=training_data$y,
                               lambda=lambda)

  index_lambda <- which.min(regr$GCV)

  tmp2 <- regr$coef[which.min(regr$GCV)]/2

  alpha <- pmax(1e-10, fcast$model$par["alpha"])

  fcast$mean <- fcast$mean +  predict_ridge(regr, test_data)[,index_lambda]  + tmp2 * (1 - (1 -
                                                                                              alpha)^n)/alpha

  if (seasonal) {
    fcast$mean <- fcast$mean * rep(tail(decomp$seasonal,
                                        m), trunc(1 + h/m))[1:h]
    fcast$fitted <- fcast$fitted * decomp$seasonal
  }

  fcast$residuals <- origx - fcast$fitted


  fcast.se <- sqrt(fcast$model$sigma2) * sqrt(test_trend *
                                                alpha^2 + 1)
  nconf <- length(level)

  fcast$mean <- ts(cumsum(c(tail(y, 1), fcast$mean))[1:h],
                   start = start(fcast$mean), frequency = frequency(y))

  fcast$lower <- fcast$upper <- ts(matrix(NA, nrow = h, ncol = nconf))
  tsp(fcast$lower) <- tsp(fcast$upper) <- tsp(fcast$mean)
  for (i in 1:nconf) {
    zt <- -qnorm(0.5 - level[i]/200)
    fcast$lower[, i] <- fcast$mean - zt * fcast.se
    fcast$upper[, i] <- fcast$mean + zt * fcast.se
  }
  fcast$x <- y
  fcast$fitted <- ts(cumsum(c(y[1], fcast$fitted)),
                     end = end(y), frequency = frequency(y))
  fcast$residuals <- y - fcast$fitted
  fcast$level <- level
  fcast$method <- "Theta"
  fcast$model <- list(alpha = alpha, drift = tmp2, sigma = fcast$model$sigma2)
  fcast$model$call <- match.call()
  return(fcast)
}



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
#' (res <- smoothf(y, h=10))
#' plot(res)
#'
#' y <- AirPassengers
#' (res <- smoothf(y, h = 10))
#' plot(res)
#'
#' y <- USAccDeaths
#' (res <- smoothf(y, h = 10))
#' plot(res)
#'
#' y <- WWWusage
#' (res <- smoothf(y, h = 10))
#' plot(res)
#'
#' y <- WWWusage
#' (res <- smoothf(y, h = 10, lags = 2L))
#' plot(res)
#'
#'
#' y <- WWWusage
#' (res <- smoothf(y, h = 10, lags = 3L))
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


  Y <- embedc(x = y, lags = lags)
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


  y <- c(y, tail(fitted, 1))
  for (i in 1:(h-1))
  {
    Y <- embedc(x = y, lags = lags)
    fitted <- switch(method,
                     "mean" = rowMeans(tcrossprod(Y, W)),
                     "median" = apply(tcrossprod(Y, W), 1,
                                      median))
    y <- c(y, tail(fitted, 1))
  }


  ans_mean <- ts(tail(y, h), start = start_preds,
                 frequency = freq_y)


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
    resid_fcast <- garch11f(res_obj$residuals,
                                   level=level, h = h)
    res_obj$mean <- ans_mean + resid_fcast$mean
    res_obj$lower <- ans_mean + resid_fcast$lower
    res_obj$upper <- ans_mean + resid_fcast$upper
  }

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
