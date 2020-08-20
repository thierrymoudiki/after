
# Data transformation -----------------------------------------------------

# calculate std's of columns
my_sd <- function(x)
{
  n <- dim(x)[1]
  return(drop(rep(1/(n-1), n) %*% (x - tcrossprod(rep.int(1, n), colMeans(x)))^2)^0.5)
}
my_sd <- compiler::cmpfun(my_sd)

# scaling matrixes
my_scale <- function(x, xm = NULL, xsd = NULL)
{
  rep_1_n <- rep.int(1, dim(x)[1])

  # centering and scaling, returning the means and sd's
  if(is.null(xm) && is.null(xsd))
  {
    xm <- colMeans(x)
    xsd <- my_sd(x)
    return(list(res = (x - tcrossprod(rep_1_n, xm))/tcrossprod(rep_1_n, xsd),
                xm = xm,
                xsd = xsd))
  }

  # centering and scaling
  if(is.numeric(xm) && is.numeric(xsd))
  {
    return((x - tcrossprod(rep_1_n, xm))/tcrossprod(rep_1_n, xsd))
  }

  # centering only
  if(is.numeric(xm) && is.null(xsd))
  {
    return(x - tcrossprod(rep_1_n, xm))
  }

  # scaling only
  if(is.null(xm) && is.numeric(xsd))
  {
    return(x/tcrossprod(rep_1_n, xsd))
  }
}
my_scale <- compiler::cmpfun(my_scale)

# re-scale a matrix
re_scale <- function(x, xm, xsd)
{
  n <- dim(x)[1]
  p <- dim(x)[2]
  stopifnot(length(xm) == p && length(xsd) == p)

  mat_xm <- tcrossprod(rep(1, n), xm)
  mat_xsd <- tcrossprod(rep(1, n), xsd)

  return (x * mat_xsd + mat_xm)
}

remove_zero_cols <- function(x, with_index = FALSE)
{
  if (with_index == FALSE)
  {
    return(x[, colSums(x == 0) != nrow(x)])
  } else {
    index <- colSums(x == 0) != nrow(x)
    return(list(mat = x[, index],
                index = index))
  }
}

is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

sort_df <- function(df, by_col)
{
  df[with(df, order(by_col)), ]
}

# From MASS::ginv
my_ginv <- function(X, tol = sqrt(.Machine$double.eps))
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
  {
    stop("'X' must be a numeric or complex matrix")
  }

  Xsvd <- La.svd(X)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive))
  {
    return(crossprod(Xsvd$vt, (1/Xsvd$d * t(Xsvd$u))))
  }
  else if(!any(Positive))
  {
    return(array(0, dim(X)[2L:1L]))
  }
  else {
    return(crossprod(Xsvd$vt[, Positive, drop = FALSE], ((1/Xsvd$d[Positive]) *
                                                   t(Xsvd$u[, Positive, drop = FALSE]))))
  }
}
my_ginv <- compiler::cmpfun(my_ginv)

log_returns <- function(x)
{
  time_x <- time(x)

  if (!is.null(dim(x)))
  {
    n <- nrow(x)
    res <- ts(log(x[-1,]/x[1:(n-1),]),
              start = time_x[2], end = time_x[n],
              frequency = frequency(x))
    return(res)
  } else {
    n <- length(x)
    res <- ts(log(x[-1]/x[1:(n-1)]),
              start = time_x[2], end = time_x[n],
              frequency = frequency(x))
    return(res)
  }

}

# Error on the validation set -----------------------------------------------------

validation_error <- function(preds, obs)
{
  difference <- preds - obs
  return(list(difference = difference,
         rmse = sqrt(colMeans(difference^2)),
         mape = 100*colMeans(abs(difference/obs))))
}

# Bagging confidence intervals -----------------------------------------------------

bag_conf_int <- function(u, level = 95)
{
  nblevels <- length(level)
  if (nblevels > 1)
  {
    lapply(1:nblevels, function(i) {
      apply(u, 1,
            function(x) {res <- t.test(x, conf.level = level[i]/100);
            return(c(res$conf.int[1], res$est, res$conf.int[2]))})
    })
  } else {
    apply(u, 1,
          function(x) {res <- t.test(x, conf.level = level/100);
          return(c(res$conf.int[1], res$est, res$conf.int[2]))})
  }
}

# AICc -----------------------------------------------------

calculate_aicc <- function(residuals, nb_predictors)
{
  n <- nrow(residuals)
  sse <- colSums(residuals^2)
  return(n*log(sse/n) + 2*(nb_predictors + 2)*(1 + (nb_predictors + 3)/(n - nb_predictors - 3)))
}
