
#' Title
#'
#' @param x
#' @param h
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' after::arimaf(AirPassengers, h=7)
#'
#' after::arimaf(Nile, h=6)
#'
#' after::arimaf(WWWusage, h=8)
#'
#' after::arimaf(USAccDeaths, h=10)
#'
arimaf <- function(x, h=5, ...)
{
  n <- length(x)

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
  list(order=order_,
       preds=predict(stats::arima0(x,
                        order = order_,
                        ...), n.ahead = h, ...))
}


