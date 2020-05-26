#' Title
#'
#' @param x
#' @param y
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
fit_ridge <- function(x, y, lambda)
{
  x <- as.matrix(x)
  y <- as.vector(y)
  nlambda <- length(lambda)

  ym <- mean(y)
  centered_y <- y - ym

  x_scaled <- base::scale(x)
  attrs <- attributes(x_scaled)
  X <- as.matrix(x_scaled[,])

  Xs <- La.svd(X)
  rhs <- crossprod(Xs$u, centered_y)
  d <- Xs$d
  nb_di <- length(d)
  div <- d ^ 2 + rep(lambda, rep(nb_di, nlambda))
  a <- drop(d * rhs) / div
  dim(a) <- c(nb_di, nlambda)
  n <- nrow(X)

  coef <- crossprod(Xs$vt, a)
  colnames(coef) <- lambda
  centered_y_hat <- X %*% coef

  fitted_values <- drop(ym +  centered_y_hat)
  colnames(fitted_values) <- lambda
  residuals <- centered_y - centered_y_hat
  GCV <- colSums(residuals^2)/(n - colSums(matrix(d^2/div, nb_di)))^2
  BIC <- n*log(colMeans(residuals^2)) + (ncol(X) + 2)*log(n)

  return(
    list(
      coef = drop(coef),
      ym = ym,
      xm = attrs$`scaled:center`,
      xsd = attrs$`scaled:scale`,
      lambda = lambda,
      fitted_values = fitted_values,
      residuals = centered_y - centered_y_hat,
      GCV = GCV,
      BIC = BIC,
      x = x,
      y = y
    )
  )
}


#' Title
#'
#' @param obj
#' @param newx
#'
#' @return
#' @export
#'
#' @examples
predict_ridge <- function(obj, newx)
{
  return(base::scale(newx, center=obj$xm,
                     scale=obj$xsd)%*%obj$coef + obj$ym)
}


