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
#'
#' n <- 1000 ; p <- 100
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#'
#' (fit_obj <- after::fit_ridge(X, y))
#'
#' matplot(fit_obj$lambda, t(fit_obj$coef), type = 'l')
#'
fit_ridge <- function(x, y, lambda=10^seq(-5, 4,
                                          length.out = 100))
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
      best_lam = lambda[which.min(GCV)],
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
#'
#'
#' n <- 10000 ; p <- 100
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#'
#' fit_obj <- after::fit_ridge(X, y)
#'
#' n_test <- 10
#'
#' predict_ridge(fit_obj, newx=matrix(rnorm(n_test * p), n_test, p),
#' cv=TRUE)
#'
#' predict_ridge(fit_obj, newx=matrix(rnorm(n_test * p), n_test, p),
#' cv=FALSE)
#'
#'
predict_ridge <- function(obj, newx, cv=TRUE)
{
  if (cv){
    return(drop(base::scale(newx, center=obj$xm,
                       scale=obj$xsd)%*%obj$coef[,which.min(obj$GCV)] + obj$ym))
  }

  return(base::scale(newx, center=obj$xm,
                     scale=obj$xsd)%*%obj$coef + obj$ym)
}


