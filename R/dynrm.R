



# Defaults:
# For non-seasonal data, p chosen using AIC from linear AR(p) model
# For seasonal data, p chosen using AIC from linear AR(p) model after
#    seasonally adjusting with STL decomposition, and P=1
# size set to average of number of inputs and number of outputs: (p+P+1)/2
# if xreg is included then size = (p+P+ncol(xreg)+1)/2



#' Neural Network Time Series Forecasts
#'
#' Feed-forward neural networks with a single hidden layer and lagged inputs
#' for forecasting univariate time series.
#'
#' A feed-forward neural network is fitted with lagged values of \code{y} as
#' inputs and a single hidden layer with \code{size} nodes. The inputs are for
#' lags 1 to \code{p}, and lags \code{m} to \code{mP} where
#' \code{m=frequency(y)}. If \code{xreg} is provided, its columns are also
#' used as inputs. If there are missing values in \code{y} or
#' \code{xreg}, the corresponding rows (and any others which depend on them as
#' lags) are omitted from the fit. A total of \code{repeats} networks are
#' fitted, each with random starting weights. These are then averaged when
#' computing forecasts. The network is trained for one-step forecasting.
#' Multi-step forecasts are computed recursively.
#'
#' For non-seasonal data, the fitted model is denoted as an NNAR(p,k) model,
#' where k is the number of hidden nodes. This is analogous to an AR(p) model
#' but with nonlinear functions. For seasonal data, the fitted model is called
#' an NNAR(p,P,k)[m] model, which is analogous to an ARIMA(p,0,0)(P,0,0)[m]
#' model but with nonlinear functions.
#'
#' @aliases print.dynrm print.dynrmmodels
#'
#' @param y A numeric vector or time series of class \code{ts}.
#' @param p Embedding dimension for non-seasonal time series. Number of
#' non-seasonal lags used as inputs. For non-seasonal time series, the default
#' is the optimal number of lags (according to the AIC) for a linear AR(p)
#' model. For seasonal time series, the same method is used but applied to
#' seasonally adjusted data (from an stl decomposition).
#' @param P Number of seasonal lags used as inputs.
#' @param size Number of nodes in the hidden layer. Default is half of the
#' number of input nodes (including external regressors, if given) plus 1.
#' @param repeats Number of networks to fit with different random starting
#' weights. These are then averaged when producing forecasts.
#' @param xreg Optionally, a vector or matrix of external regressors, which
#' must have the same number of rows as \code{y}. Must be numeric.
#' @param model Output from a previous call to \code{dynrm}. If model is
#' passed, this same model is fitted to \code{y} without re-estimating any
#' parameters.
#' @param subset Optional vector specifying a subset of observations to be used
#' in the fit. Can be an integer index vector or a logical vector the same
#' length as \code{y}. All observations are used by default.
#' @param scale_inputs If TRUE, inputs are scaled by subtracting the column
#' means and dividing by their respective standard deviations. If \code{lambda}
#' is not \code{NULL}, scaling is applied after Box-Cox transformation.
#' @param x Deprecated. Included for backwards compatibility.
#' @param \dots Other arguments passed to \code{\link[nnet]{nnet}} for
#' \code{dynrm}.
#' @inheritParams forecast
#'
#' @return Returns an object of class "\code{dynrm}".
#'
#' The function \code{summary} is used to obtain and print a summary of the
#' results.
#'
#' The generic accessor functions \code{fitted.values} and \code{residuals}
#' extract useful features of the value returned by \code{dynrm}.
#'
#' \item{model}{A list containing information about the fitted model}
#' \item{method}{The name of the forecasting method as a character string}
#' \item{x}{The original time series.}
#' \item{xreg}{The external regressors used in fitting (if given).}
#' \item{residuals}{Residuals from the fitted model. That is x minus fitted values.}
#' \item{fitted}{Fitted values (one-step forecasts)}
#' \item{...}{Other arguments}
#'
#' @author Rob J Hyndman and Gabriel Caceres
#' @keywords ts
#' @examples
#' #fit <- dynrm(lynx)
#' #fcast <- forecast(fit)
#' #plot(fcast)
#'
#' ## Arguments can be passed to nnet()
#' #fit <- dynrm(lynx, decay=0.5, maxit=150)
#' #plot(forecast(fit))
#' #lines(lynx)
#'
#' ## Fit model to first 100 years of lynx data
#' #fit <- dynrm(window(lynx,end=1920), decay=0.5, maxit=150)
#' #plot(forecast(fit, h=14))
#' #lines(lynx)
#'
#' ## Apply fitted model to later data, including all optional arguments
#' #fit2 <- dynrm(window(lynx,start=1921), model=fit)
#'
#' @export
dynrm_fit <- function(y,
                   p,
                   P = 1,
                   xreg = NULL,
                   fit_func = stats::lm.fit,
                   predict_func = predict.lm,
                   fit_params = NULL,
                   lambda = NULL,
                   scale_inputs = TRUE,
                   x = y,
                   seed = 123,
                   ...) {

  yname <- deparse(substitute(y))

  # transform data
  if (!is.null(lambda)) {
    xx <- forecast::BoxCox(x, lambda)
    lambda <- attr(xx, "lambda")
  } else {
    xx <- x
  }

  # scale series
  scalex <- NULL
  if (scale_inputs) {
    tmpx <- scale(xx, center = TRUE, scale = TRUE)
    scalex <- list(
      center = attr(tmpx, "scaled:center"),
      scale = attr(tmpx, "scaled:scale")
    )
    xx <-
      after::scale(xx, center = scalex$center, scale = scalex$scale)
    xx <- xx[, 1]
  }


  # check xreg class & dim
  xxreg <- NULL
  scalexreg <- NULL
  if (!is.null(xreg)) {
    xxreg <- xreg <- as.matrix(xreg)

    if (length(x) != NROW(xreg)) {
      stop("Number of rows in xreg does not match series length")
    }

    # scale xreg
    if (scale_inputs) {
      tmpx <- scale(xxreg, center = TRUE, scale = TRUE)
      scalexreg <- list(
        center = attr(tmpx, "scaled:center"),
        scale = attr(tmpx, "scaled:scale")
      )

      xxreg <- scale(xxreg,
                     center = scalexreg$center,
                     scale = scalexreg$scale)
    }

  }


  # Set up lagged matrix
  n <- length(xx)
  xx <- as.ts(xx)
  m <- max(round(frequency(xx)), 1L)

  if (m == 1) {
    if (missing(p)) {
      p <- max(length(stats::ar(xx)$ar), 1)
    }

    if (p >= n) {
      warning("Reducing number of lagged inputs due to short series")
      p <- n - 1
    }

    lags <- 1:p

    if (P > 1) {
      warning("Non-seasonal data, ignoring seasonal lags")
    }

    P <- 0

  } else {
    if (missing(p)) {
      if (n > 2 * m) {
        x.sa <- forecast::seasadj(forecast::mstl(forecast::na.interp(xx)))
      } else {
        x.sa <- forecast::na.interp(xx)
      }

      p <- max(length(stats::ar(x.sa)$ar), 1)

    }

    if (p >= n) {
      warning("Reducing number of lagged inputs due to short series")
      p <- n - 1
    }

    if (P > 0 && n >= m * P + 2) {
      lags <- sort(unique(c(1:p, m * (1:P))))
    } else {
      lags <- 1:p
      if (P > 0) {
        warning("Series too short for seasonal lags")
        P <- 0
      }
    }

  }

  maxlag <- max(lags)
  nlag <- length(lags)
  y <- xx[-(1:maxlag)]
  lags.X <- matrix(NA_real_, ncol = nlag, nrow = n - maxlag)

  for (i in 1:nlag)
    lags.X[, i] <- xx[(maxlag - lags[i] + 1):(n - lags[i])]

  cat("lags.X: ", "\n")
  print(head(lags.X))
  print(tail(lags.X))
  cat("\n")

  # Add xreg into lagged matrix
  lags.X <- cbind(lags.X, xxreg[-(1:maxlag), ])


  # Remove missing values if present
  j <- complete.cases(lags.X, y)


  ## Stop if there's no data to fit (e.g. due to NAs or NaNs)
  if (NROW(lags.X[j, , drop = FALSE]) == 0) {
    stop("No data to fit (possibly due to NA or NaN)")
  }

  cat("y: ", "\n")
  print(head(y))
  print(tail(y))
  cat("\n")

  # fit
  set.seed(seed) # in case the algo is randomized
  fit <- do.call(what = fit_func,
                 args = c(list(x = lags.X[j, , drop = FALSE],
                               y = y[j]), fit_params))


  # Return results
  out <- list()
  out$x <- as.ts(x)
  out$m <- m
  out$p <- p
  out$P <- P
  out$scalex <- scalex
  out$scalexreg <- scalexreg
  out$xreg <- xreg
  out$lambda <- lambda
  out$model <- fit
  out$dynrmgs <- list(...)


  # Fitted values
  fits <- try(predict_func(fit,
                           newdata = lags.X[j, , drop = FALSE]),
              silent = TRUE)

  if (class(fits) == "try-error" || is.null(fits))
    fits <- try(predict_func(fit,
                             newx = lags.X[j, , drop = FALSE]),
                silent = TRUE)

  if (scale_inputs) {
    fits <- fits * scalex$scale + scalex$center
  }

  fits <- ts(fits, end = end(out$x),
             frequency = frequency(out$x))
  if (!is.null(lambda)) {
    fits <- forecast::InvBoxCox(fits, lambda)
  }

  out$fitted <- fits
  out$residuals <- out$x - out$fitted
  out$lags <- lags
  out$predict_func <- predict_func
  out$series <- yname
  out$method <- paste("DynRM ", p, sep = "")

  if (P > 0) {
    out$method <- paste(out$method, ",", P, sep = "")
  }

  if (P > 0) {
    out$method <- paste(out$method, "[", m, "]", sep = "")
  }

  out$call <- match.call()

  return(structure(out, class = c("dynrm")))
}


#' Forecasting using neural network models
#'
#' Returns forecasts and other information for univariate neural network
#' models.
#'
#' Prediction intervals are calculated through simulations and can be slow.
#' Note that if the network is too complex and overfits the data, the residuals
#' can be arbitrarily small; if used for prediction interval calculations, they
#' could lead to misleadingly small values. It is possible to use out-of-sample
#' residuals to ameliorate this, see examples.
#'
#' @param object An object of class "\code{dynrm}" resulting from a call to
#' \code{\link{dynrm}}.
#' @param h Number of periods for forecasting. If \code{xreg} is used, \code{h}
#' is ignored and the number of forecast periods is set to the number of rows
#' of \code{xreg}.
#' @param PI If TRUE, prediction intervals are produced, otherwise only point
#' forecasts are calculated. If \code{PI} is FALSE, then \code{level},
#' \code{fan}, \code{bootstrap} and \code{npaths} are all ignored.
#' @param level Confidence level for prediction intervals.
#' @param fan If \code{TRUE}, level is set to \code{seq(51,99,by=3)}. This is
#' suitable for fan plots.
#' @param xreg Future values of external regressor variables.
#' @param ... Additional arguments passed to \code{\link{simulate.dynrm}}
#' @inheritParams forecast
#'
#' @return An object of class "\code{forecast}".
#'
#' The function \code{summary} is used to obtain and print a summary of the
#' results, while the function \code{plot} produces a plot of the forecasts and
#' prediction intervals.
#'
#' The generic accessor functions \code{fitted.values} and \code{residuals}
#' extract useful features of the value returned by \code{after.dynrm}.
#'
#' An object of class "\code{forecast}" is a list containing at least the
#' following elements:
#'   \item{model}{A list containing information about the fitted model}
#'   \item{method}{The name of the forecasting method as a character string}
#'   \item{mean}{Point forecasts as a time series}
#'   \item{lower}{Lower limits for prediction intervals}
#'   \item{upper}{Upper limits for prediction intervals}
#'   \item{level}{The confidence values associated with the prediction intervals}
#'   \item{x}{The original time series (either \code{object} itself or the time series
#'            used to create the model stored as \code{object}).}
#'   \item{xreg}{The external regressors used in fitting (if given).}
#'   \item{residuals}{Residuals from the fitted model. That is x minus fitted values.}
#'   \item{fitted}{Fitted values (one-step forecasts)}
#'   \item{...}{Other arguments}
#'
#' @author Rob J Hyndman and Gabriel Caceres
#' @seealso \code{\link{dynrm}}.
#' @keywords ts
#' @examples
#' ## Fit & forecast model
#' #fit <- dynrm(USAccDeaths, size=2)
#' #fcast <- forecast(fit, h=20)
#' #plot(fcast)
#'
#' \dontrun{
#' ## Include prediction intervals in forecast
#' #fcast2 <- forecast(fit, h=20, PI=TRUE, npaths=100)
#' #plot(fcast2)
#'
#' ## Set up out-of-sample innovations using cross-validation
#' #fit_cv <- CVar(USAccDeaths,  size=2)
#' #res_sd <- sd(fit_cv$residuals, na.rm=TRUE)
#' #myinnovs <- rnorm(20*100, mean=0, sd=res_sd)
#' ## Forecast using new innovations
#' #fcast3 <- forecast(fit, h=20, PI=TRUE, npaths=100, innov=myinnovs)
#' #plot(fcast3)
#' #}
#'
#' @export
dynrm_predict <-
  function(out,
           h = ifelse(out$m > 1, 2 * out$m, 10),
           level = c(80, 95),
           fan = FALSE,
           xreg = NULL,
           lambda = out$lambda,
           ci = c("E", "A", "T",
                  "garch", "gaussian"),
           ...) {

    tspx <- tsp(out$x)

    ci <- match.arg(ci)

    if (fan) {
      level <- seq(51, 99, by = 3)
    } else {
      if (min(level) > 0 && max(level) < 1) {
        level <- 100 * level
      } else if (min(level) < 0 || max(level) > 99.99) {
        stop("Confidence limit out of range")
      }
    }


    # Check if xreg was used in fitted model
    if (is.null(out$xreg)) {
      if (!is.null(xreg)) {
        warning("External regressors were not used in fitted model, xreg will be ignored")
      }
      xreg <- NULL
    }
    else {
      if (is.null(xreg)) {
        stop("No external regressors provided")
      }
      xreg <- as.matrix(xreg)
      if (NCOL(xreg) != NCOL(out$xreg)) {
        stop("Number of external regressors does not match fitted model")
      }
      if (!identical(colnames(xreg), colnames(out$xreg))) {
        warning(
          "xreg contains different column names from the xreg used in training. Please check that the regressors are in the same order."
        )
      }
      h <- NROW(xreg)
    }


    fcast <- numeric(h)
    xx <- out$x
    xxreg <- xreg
    if (!is.null(lambda)) {
      xx <- forecast::BoxCox(xx, lambda)
      lambda <- attr(xx, "lambda")
    }


    # Check and apply scaling of fitted model
    if (!is.null(out$scalex)) {
      xx <- after::scale(xx,
                         center = out$scalex$center,
                         scale = out$scalex$scale)
      if (!is.null(xreg)) {
        xxreg <- stats::scale(xreg,
                              center = out$scalexreg$center,
                              scale = out$scalexreg$scale)
      }
    }


    # Get lags used in fitted model
    lags <- out$lags
    maxlag <- max(lags)
    flag <- rev(tail(xx, n = maxlag))

    cat("flag", "\n")
    print(flag)
    cat("\n")

    # Iterative 1-step forecast
    for (i in 1:h)
    {
      newdata <- c(flag[lags], xxreg[i, ])

      cat("newdata", "\n")
      print(newdata)
      cat("\n")

      newdata_ <- as.matrix(rbind(newdata,
                                  rnorm(length(newdata))))

      preds <- try(as.numeric(out$predict_func(out$model,
                                               newdata = newdata_)[1]), # do not change (ever)
                   silent = TRUE)
      if (class(preds) == "try-error")
      {
        preds <- try(as.numeric(out$predict_func(out$model,
                                                 newx = newdata_)[1]), # do not change (ever)
                     silent = TRUE)
        if (class(preds) == "try-error")
        {
          preds <- NA
        }
      }

      cat("preds", "\n")
      print(preds)
      cat("\n")

      fcast[i] <- preds

      flag <- c(fcast[i], flag[-maxlag])
    }


    # Re-scale point forecasts
    if (!is.null(out$scalex)) {
      fcast <- fcast * out$scalex$scale + out$scalex$center
    }


    # Add ts properties
    fcast <- ts(fcast, start = tspx[2] + 1 / tspx[3],
                frequency = tspx[3])


    # Back-transform point forecasts
    if (!is.null(lambda)) {
      fcast <- forecast::InvBoxCox(fcast, lambda)
    }


    # Compute prediction intervals

    if (ci == "E")
    {
      resid_fcast <- forecast::forecast(forecast::ets(out$residuals),
                                        h = h)
      out$mean <- fcast + resid_fcast$mean
      out$lower <- fcast + resid_fcast$lower
      out$upper <- fcast + resid_fcast$upper
    }

    if (ci == "A")
    {
      resid_fcast <- forecast::forecast(forecast::auto.arima(out$residuals),
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
      resid_fcast <- after::garch11f(out$residuals,
                                     h = h)
      out$mean <- fcast + resid_fcast$mean
      out$lower <- fcast + resid_fcast$lower
      out$upper <- fcast + resid_fcast$upper
    }

    out$level <- level

    return(structure(out, class = "forecast"))
  }
