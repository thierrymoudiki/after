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


  # 1 - response ------------------------------------------------------------
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


  # 2 - covariates ------------------------------------------------------------

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

  # 3 - colnames ------------------------------------------------------------
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
