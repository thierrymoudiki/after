
#' Outliers treatment
#'
#' @param x a time series
#' @param plot a boolean; plot time series or not
#' @param replace a boolean; replace outliers or not? (if TRUE, creates a new time series)
#'
#' @return
#' @export
#'
#' @examples
#'
#' x <- rnorm(100)
#' x[5] <- 10
#' x[10] <- 10
#' x[20] <- 10
#' x[70] <- 10
#' after::treat_outliers(x, plot = TRUE)
#'
treat_outliers <- function(x, plot=FALSE, replace=TRUE)
{
  # adapted from: https://stats.stackexchange.com/questions/1142/simple-algorithm-for-online-outlier-detection-of-a-generic-time-series

  x <- as.ts(x)
  median_x <- median(x)
  if(frequency(x)>1)
    resid <- stl(x,s.window="periodic",robust=TRUE)$time.series[,3]
  else
  {
    tt <- 1:length(x)
    resid <- residuals(loess(x ~ tt))
  }
  resid.q <- quantile(resid,prob=c(0.25,0.75))
  iqr <- diff(resid.q)
  limits <- resid.q + 1.5*iqr*c(-1,1)
  score <- abs(pmin((resid-limits[1])/iqr,0) + pmax((resid - limits[2])/iqr,0))
  criterion1 <- score > 0
  diff_x <- diff(x)
  high_changes <- boxplot.stats(diff_x)$out
  criterion2 <- pmatch(high_changes, diff_x)
  criterion2 <- criterion2[!is.na(criterion2)]

  if(plot)
  {
    if (replace)
      par(mfrow=c(2, 1))

    plot(x)
    x2 <- ts(rep(NA,length(x)))
    x2[criterion1] <- x[criterion1]
    x2[criterion2] <- x[criterion2]
    tsp(x2) <- tsp(x)
    points(x2, pch=19, col="red")

    if (!replace)
    {
      return(invisible(list(pct=(sum(criterion1)+length(criterion2))/length(x))))
    } else {
      x_ <- x
      x_[criterion1] <- median_x
      x_[criterion2] <- median_x
      plot(x_, ylim=c(min(x), max(x)),
           lwd=2, col="gray60")
      return(invisible(list(x_= x_,
                            pct=(sum(criterion1)+length(criterion2))/length(x))))
    }

  } else {

    if (!replace)
    {
      return(list(pct=(sum(criterion1)+length(criterion2))/length(x)))
    } else {
      x_ <- x
      x_[criterion1] <- median_x
      x_[criterion2] <- median_x
      return(list(x_ = x_,
                  pct=(sum(criterion1)+length(criterion2))/length(x)))
    }

  }

}


#' Scale a univariate time series
#'
#' @param x
#' @param center
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
scale.ts <- function(x, center = TRUE, scale = TRUE) {
  tspx <- tsp(x)
  x <- as.ts(scale.default(x, center = center, scale = scale))
  tsp(x) <- tspx
  return(x)
}

