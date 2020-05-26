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
