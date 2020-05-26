


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
