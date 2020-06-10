
# Defaults:
# For non-seasonal data, p chosen using AIC from linear AR(p) model
# For seasonal data, p chosen using AIC from linear AR(p) model after
#    seasonally adjusting with STL decomposition, and P=1
# size set to average of number of inputs and number of outputs: (p+P+1)/2
# if xreg is included then size = (p+P+ncol(xreg)+1)/2



#' Fit a dynrm
#'
#' @param y
#' @param p
#' @param P
#' @param xreg
#' @param fit_func
#' @param predict_func
#' @param fit_params
#' @param lambda
#' @param scale_inputs
#' @param x
#' @param seed
#' @param ...
#' @example
#'
#'
#'
#' dat <- list(lynx, USAccDeaths, Nile, WWWusage, fdeaths, AirPassengers)
#'
#' h <- 10
#'
#'
#' # Example 1: glmnet
#'
#' par(mfrow=c(3, 2))
#'
#' for (i in 1:length(dat))
#' {
#'  x <- dat[[i]]
#'  cat("Dataset #", i, ":", length(x), "points \n")
#'
#'  obj <- after::dynrm_fit(x, fit_func = function(x, y) glmnet::glmnet(x, y, lambda = 0.01),
#'  p=1, predict_func = glmnet::predict.glmnet, xreg = as.matrix(1:length(x)))
#'
#'  plot(after::dynrm_predict(obj, xreg = as.matrix((length(x)+1):(length(x)+h-1)),
#'  ci="A", h=h))
#' }
#'
#'
#'
#' # Example 2: ridge
#'
#' par(mfrow=c(3, 2))
#'
#' for (i in 1:length(dat))
#' {
#'  x <- dat[[i]]
#'  cat("Dataset #", i, ":", length(x), "points \n")
#'
#'  obj <- after::dynrm_fit(x, fit_func = after::fit_ridge,
#'  p=1, predict_func = after::predict_ridge, xreg = as.matrix(1:length(x)))
#'
#'  plot(after::dynrm_predict(obj, xreg = as.matrix((length(x)+1):(length(x)+h-1)),
#'  ci="A", h=h))
#' }
#'
dynrm_fit <- function(y,
                   p,
                   P = 1,
                   xreg = NULL,
                   fit_func = stats::lm.fit,
                   predict_func = predict.lm,
                   fit_params = NULL,
                   lambda = NULL, alpha = NULL,
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
      tmpx <- base::scale(xxreg, center = TRUE, scale = TRUE)
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

   # cat("lags.X: ", "\n")
   # print(head(lags.X))
   # print(tail(lags.X))
   # cat("\n")

  if (is.null(alpha) == FALSE)
  {
    n_alphas <- nrow(lags.X)
    alphas <- alpha*((1-alpha)^((n_alphas + p + P - 2):0))
     # cat("alphas: ", "\n")
     # print(dim(exp(embedc(alphas, p + P -1))))
     # cat("\n")
    # Add xreg into lagged matrix
    lags.X <- cbind(lags.X*exp(embedc(alphas, ncol(lags.X)-1)),
                    xxreg[-(1:maxlag), ])
  } else {
    # Add xreg into lagged matrix
    lags.X <- cbind(lags.X, xxreg[-(1:maxlag), ])
  }

  # Remove missing values if present
  j <- complete.cases(lags.X, y)


  ## Stop if there's no data to fit (e.g. due to NAs or NaNs)
  if (NROW(lags.X[j, , drop = FALSE]) == 0) {
    stop("No data to fit (possibly due to NA or NaN)")
  }

   # cat("y: ", "\n")
   # print(head(y))
   # print(tail(y))
   # cat("\n")

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
  out$sigma <- sd(out$residuals)

  if (P > 0) {
    out$method <- paste(out$method, ",", P, sep = "")
  }

  if (P > 0) {
    out$method <- paste(out$method, "[", m, "]", sep = "")
  }

  out$call <- match.call()

  # return
  invisible(structure(out, class = c("dynrm")))
}



#' Predict from dynrm
#'
#' @param out
#' @param h
#' @param level
#' @param fan
#' @param xreg
#' @param lambda
#' @param ci
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' dat <- list(lynx, USAccDeaths, Nile, WWWusage, fdeaths, AirPassengers)
#'
#' h <- 10
#'
#' # Example 2: ridge
#'
#' par(mfrow=c(3, 2))
#'
#' for (i in 1:length(dat))
#' {
#'  x <- dat[[i]]
#'  cat("Dataset #", i, ":", length(x), "points \n")
#'
#'  obj <- after::dynrm_fit(x, fit_func = after::fit_ridge,
#'  predict_func = after::predict_ridge, xreg = as.matrix(1:length(x)))
#'
#'   plot(after::dynrm_predict(obj, xreg = as.matrix((length(x)+1):(length(x)+h-1)),
#'   ci="A", h=h))
#' }
#'
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
        xxreg <- base::scale(xreg,
                              center = out$scalexreg$center,
                              scale = out$scalexreg$scale)
      }
    }


    # Get lags used in fitted model
    lags <- out$lags
    maxlag <- max(lags)
    flag <- rev(tail(xx, n = maxlag))

    # cat("flag", "\n")
    # print(flag)
    # cat("\n")

    # Iterative 1-step forecast
    for (i in 1:h)
    {
      newdata <- c(flag[lags], xxreg[i, ])

      # cat("newdata", "\n")
      # print(newdata)
      # cat("\n")

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

      # cat("preds", "\n")
      # print(preds)
      # cat("\n")

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
