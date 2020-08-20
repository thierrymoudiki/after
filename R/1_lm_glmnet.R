

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
fit_glmnet <- function(x, y,
                       alpha=0.5,
                       lambda=10^seq(-10, 10, length.out = 100),
                       ...)
{
  fit_obj <- glmnet::glmnet(x=as.matrix(x),
                            y=y,
                            alpha=alpha,
                            lambda=lambda,
                            ...)

  residuals_train <- glmnet::predict.glmnet(fit_obj,
                                            newx=as.matrix(x),
                                            ...) - y
  n_train <- nrow(x)
  (bic <- n_train*log(colMeans(residuals_train^2)) + (fit_obj$df + 2)*log(n_train))

  fit_obj$lambda <- lambda

  fit_obj$bic <- bic

  return(fit_obj)
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
predict_glmnet <- function(obj, newx, ...)
{

  return(as.numeric(glmnet::predict.glmnet(obj,
                                    newx = as.matrix(newx),
                                    s = obj$lambda[which.min(obj$bic)],
                                    ...)))
}
