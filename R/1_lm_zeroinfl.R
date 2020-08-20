

#' Title
#'
#' @param x
#' @param y
#' @param alpha
#' @param lambda
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fit_zeroinfl <- function(x, y, ...)
{
  dfyx <- cbind.data.frame(y, x)
  print(head(dfyx))
  print(tail(dfyx))
  return(pscl::zeroinfl(y ~ .,
                        data = dfyx,
                        ...))
}


#' Title
#'
#' @param obj
#' @param newx
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict_zeroinfl <- function(obj, newx, ...)
{

  return(as.numeric(predict(obj,
                            newx = as.matrix(newx),
                            type = "response",
                            ...)))
}
