#' (used in ensemblef)
#'
#' @param fitted_objects
#'
#' @return
#' @export
#'
#' @examples
check_list_obj <- function(fitted_objects)
{
  x <- fitted_objects[[1]]$x
  for (i in 2:length(fitted_objects))
  {
    if (!all.equal(fitted_objects[[i]]$x, x))
      return (FALSE)
  }
  return(TRUE)
}

#' Embed time series
#'
#' @param x
#' @param y
#' @param lags
#' @param model_matrix
#'
#' @return
#' @export
#'
#' @examples
#'
#' x <- 1:10
#' set.seed(2020)
#' df <- data.frame(x1 = x,
#'                  x2 = c("a", "b", "c", "b", "a",
#'                         "a", "a", "c", "b", "c"),
#'                  x3=sort(rnorm(10)))
#' df$x2 <- as.factor(df$x2)
#'
#' lags <- 3L
#'
#' cat("1 - data -----")
#' cat("\n")
#' print(df)
#' cat("\n")
#'
#' cat("2 - training (model_matrix==FALSE) -----")
#' cat("\n")
#' cat("\n")
#' print(after::embed_reg(y=df$x1, x=df[, c("x2", "x3")],
#'                        lags = lags))
#' cat("\n")
#'
#' cat("2 - training (model_matrix==TRUE) -----")
#' cat("\n")
#' cat("\n")
#' print(after::embed_reg(y=df$x1, x=df[, c("x2", "x3")],
#'                        lags = lags, encoding = "model_matrix"))
#' cat("\n")
#'
#' cat("3 - training target-based -----")
#' cat("\n")
#' cat("\n")
#' print(after::embed_reg(y=df$x1, x=df[, c("x2", "x3")],
#'                        lags = lags, encoding = "target_based"))
#' cat("\n")
#'
embed_reg <- function(x, y, lags = 1,
                      encoding=c('None',
                                 'model_matrix',
                                 'target_based'),
                      testing=FALSE)
{
  stopifnot(is.data.frame(x))

  encoding <- match.arg(encoding)

  if (lags < 1) # but useless
  {
    return(return(list(
      y = y,
      x = x,
      df = cbind.data.frame(y, x)
    )))
  }

  lags_plus <- lags + 1


  temp_y <- stats::embed(as.vector(y),
                         lags_plus)

  if (!testing)
  {
    # At training time
    response_y <- temp_y[, 1]
    covariates_y <- temp_y[,-1]
    if (lags > 1)
    {
      colnames(covariates_y) <- paste0("y_l", 1:lags)
    } else {
      names(covariates_y) <- paste0("y_l", 1)
    }

  } else {
    # At testing time
    covariates_y <- temp_y[, 1:lags]
    if (lags > 1)
    {
      colnames(covariates_y) <- paste0("y_l", 1:lags)
    } else {
      names(covariates_y) <- paste0("y_l", 1)
    }
  }


  `%op%` <- foreach::`%do%`
  j <- NULL


  if(encoding=="None")
  {
    x_ <- as.matrix(x)
  }

  if(encoding=="model_matrix")
  {
    df_temp <- cbind.data.frame(y, x)
    x_ <- model.matrix(y ~ ., data = df_temp)[,-1]
  }

  if (encoding=="target_based")
  {
    encoder <- after::target_based_encoder(x = x, y = y)
    x_ <- encoder$newx
  }


  if (!testing)
  {
    covariates_x <- foreach::foreach(j = 1:ncol(x_),
                                     .combine = cbind) %op% {
                                       stats::embed(as.vector(x_[, j]), lags_plus)[, -1]
                                     }
    names_x_ <- colnames(x_)
    colnames(covariates_x) <- unlist(sapply(1:ncol(x_), function(col) paste0(names_x_[col],
                                                                             "_l", 1:lags)))
  } else {
    covariates_x <- foreach::foreach(j = 1:ncol(x_),
                                     .combine = cbind) %op% {
                                       temp <- stats::embed(as.vector(x_[, j]),
                                                            lags_plus)
                                       temp[, -ncol(temp)]
                                     }
    names_x_ <- colnames(x_)
    colnames(covariates_x) <- unlist(sapply(1:ncol(x_), function(col) paste0(names_x_[col],
                                                                             "_l", 1:lags)))
  }

  if (!testing)
  {
    covariates <- cbind.data.frame(covariates_y,
                                   covariates_x)
    #colnames(covariates) <-

    if(encoding == "target_based")
    {
      return(list(
        y = response_y,
        x = covariates,
        df = cbind.data.frame(y=response_y, covariates),
        codes = encoder$codes
      ))
    }

    # At training time
    return(list(
      y = response_y,
      x = covariates,
      df = cbind.data.frame(y=response_y, covariates)
    ))

  } else {
    # At testing time
    covariates <- cbind.data.frame(covariates_y,
                                   covariates_x)
    #colnames(covariates) <-
    return(list(
      y = rep(NA, nrow(covariates)),
      x = covariates))
  }

}


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

