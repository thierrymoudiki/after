#' Simulate a correlated variable with a given correlation
#'
#' @param x response variable
#' @param rho target correlation
#' @param sim_fun simulation function
#'
#' @return
#' @export
#'
#' @examples
#'
#' B <- 250
#' res <- rep(0, B)
#'
#' for (i in 1:250)
#' {
#'   y <- rnorm(1000)
#'   res[i] <- get_correlated_var(y, rho = 0.1)$cors
#' }
#'
#' hist(res)
#'
get_correlated_var <- function(x, rho, sim_fun=rnorm, seed=123) {

  set.seed(seed)

  n <- length(x)
  C <- matrix(rep(rho, 4), nrow = 2, ncol = 2)
  diag(C) <- 1
  C <- chol(C)

  X2 <- sim_fun(n)
  X <- cbind(x, X2)

  # induce correlation (does not change X1)
  X_ <- X %*% C

  return(list(vars=X_ %*% C,
              cors = cor(X_)[1, 2],
              check=all.equal(x, X_[,1])))
}


#' Target-based encoder
#'
#' @param x matrix; explanatory variables
#' @param y vector; response
#' @param rho float; desired correlation
#'
#' @return a list, encoded regressors and codes
#' @export
#'
#' @examples
#'
#' n <- 100
#' X <- cbind.data.frame(as.factor(sample(x = c(1, 2, 3),
#' size = n, replace = TRUE)), as.factor(sample(x = c(0, 1),
#' size = n, replace = TRUE)))
#' X$X3 <- rnorm(n)
#' colnames(X) <- c("X1", "X2", "X3")
#'
#' y <- rt(nrow(X), df=2)
#'
#' z <- after::target_based_encoder(X, y)$newx
#'
#' head(X)
#' head(z)
#' tail(X)
#' tail(z)
#'
#' cor(y, z[,2])
#' cor(y, z[,1])
#'
#' sd(z[,2])
#' sd(z[,1])
#'
target_based_encoder <- function(x, y, rho=0, seed=123)
{
  n <- nrow(x)
  p <- ncol(x)
  newx <- matrix(0, nrow = nrow(x), ncol = ncol(x))

  y_ <- get_correlated_var(y, rho, seed = seed)$vars

  pb <- txtProgressBar(min = 0, max = p, style = 3)

  codes <- vector("list", 0)

  col_names <- colnames(x)

  for (j in 1:p)
  {

    x_j <- x[, j]

    if (is.factor(x_j))
    {

      levels_ <- unique(x_j)
      codes[[col_names[j]]] <- vector("list", 0)
      for (l in levels_)
      {
        filter_ <- (x_j == l)
        replace_val <- sum(y_[filter_])
        codes[[col_names[j]]][[l]] <- replace_val
        newx[filter_, j] <- replace_val
      }

    } else {

      newx[, j] <- x[, j]

    }
    utils::setTxtProgressBar(pb, j)
  }

  colnames(newx) <- colnames(x)

  return(list(newx=newx,
              codes=codes))
}
