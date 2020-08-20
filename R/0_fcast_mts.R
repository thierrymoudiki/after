#' @useDynLib regularizedmts
#' @importFrom Rcpp sourceCpp
fcast_obj_mts <- function(fit_obj,
                          h = 5,
                          type_ci = "none",
                          type_forecast = c("recursive", "direct"),
                          level = 95)
{
  # Fit method
  fit_method <- fit_obj$class_res
  # Method for forecasting the residuals
  stopifnot(!is.na(match(
    type_ci, c("none", "arima",
               "ets", "theta",
               "VAR")
  )))
  type_forecast <- match.arg(type_forecast)
  # index of chosen predictors (TRUE if all, NumericVector if col_sample < 1)
  col_index <- fit_obj$col_index
  # non-zero columns index in hidden layer
  hidden_layer_index <- fit_obj$hidden_layer_index

  # 1 - recursive forecasts -------------------------------------------------

  # recursive forecasts
  if (type_forecast == "recursive")
  {
    if (fit_method == "ridge2")
    {
      # if col_sample < 1
      if(is.numeric(col_index) == FALSE)
      {
        # observed values (minus lagged)
        y <- fit_obj$y
        # observed values (minus selected columns)
        lags <- fit_obj$lags
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd
        w <- fit_obj$w
        g <- fit_obj$activ


        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))

          } else {

            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }

          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
        }
      } else {
        # observed values (minus lagged)
        y <- fit_obj$y
        # observed values (minus selected columns)
        y_reduced <- as.matrix(fit_obj$y[ , col_index])
        lags <- fit_obj$lags
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd

        w <- fit_obj$w
        g <- fit_obj$activ

        for (i in 1:h)
        {
          newx <- reformat_cpp(y_reduced, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))
          } else {
            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }

          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
          y_reduced <- as.matrix(y[ , col_index])
        }
      }
    }

    if (fit_method == "ridge")
    {
      y <- fit_obj$y
      lags <- fit_obj$lags

      if (fit_obj$nb_hidden > 0)
      {
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd
        w <- fit_obj$w
        g <- fit_obj$activ

        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))

          } else {

            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }

          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
        }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)
          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
        }
      }
    }

    if (fit_method == "lm" || fit_method == "nnls")
    {
      y <- fit_obj$y
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ

      if (fit_obj$nb_hidden > 0)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))

          } else {

            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }

          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
        }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)
          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
        }
      }
    }

    if (fit_method == "mgaussian")
    {
      y <- fit_obj$y
      lags <- fit_obj$lags

      if (fit_obj$nb_hidden > 0)
      {
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd
        w <- fit_obj$w
        g <- fit_obj$activ
        fit_glmnet <- fit_obj$fit_obj

        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w))
          } else {
            newx <- cbind(newx, g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)
            ) %*% w))
          }

          preds <- predict(fit_glmnet,
                           type = "response",
                           s = fit_obj$s,
                           newx = newx)[, , 1]
          y <- rbind_vecmat_cpp(preds, y)
        }
      } else {
        newx <- reformat_cpp(y, lags)
        fit_glmnet <- fit_obj$fit_obj
        preds <- predict(fit_glmnet,
                         type = "response",
                         s = fit_obj$s,
                         newx = newx)[, , 1]
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    if (fit_method == "glmboost" || fit_method == "xgboost")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ
      fit_glmboost <- fit_obj$fit_obj
      nb_series <- fit_obj$nb_series
      ym <- fit_obj$ym
      xm <- fit_obj$xm
      scales <- fit_obj$scales
      direct_link <- fit_obj$direct_link

      if (direct_link == TRUE)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

            if(fit_obj$hidden_layer_bias == FALSE)
            {
              newx <- cbind(newx, matrix(g(my_scale(
                newx, xm = nn_xm,
                xsd = nn_xsd
              ) %*% w)[, hidden_layer_index], nrow = 1))
            } else {
              newx <- cbind(newx, matrix(g(my_scale(
                cbind(1, newx), xm = c(1, nn_xm),
                xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
            }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          y <- rbind_vecmat_cpp(preds, y)
        }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

            if(fit_obj$hidden_layer_bias == FALSE)
            {
              newx <- matrix(g(my_scale(
                newx, xm = nn_xm,
                xsd = nn_xsd
              ) %*% w)[, hidden_layer_index], nrow = 1)
            } else {
              newx <- matrix(g(my_scale(
                cbind(1, newx), xm = c(1, nn_xm),
                xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1)
            }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          y <- rbind_vecmat_cpp(preds, y)
        }
      }
    }

    if (fit_method == "gaussian" ||
        fit_method == "matern32" || fit_method == "matern52")
    {
      y <- fit_obj$y
      xm <- fit_obj$xm
      ym <- fit_obj$ym
      xsd <- fit_obj$xsd
      lags <- fit_obj$lags
      xreg <- fit_obj$xreg
      sigma <- fit_obj$sigma
      l <- fit_obj$l
      mat_coefs <- fit_obj$coef
      class_res <- fit_obj$class_res
      series_names <- fit_obj$series_names
      Kxy <- switch(
        fit_method,
        "gaussian" = gaussian_kxy_cpp,
        "matern32" = matern32_kxy_cpp,
        "matern52" = matern52_kxy_cpp
      )

      for (i in 1:h)
      {
        newx <- reformat_cpp(y, lags)
        K_star <- Kxy(
          x = xreg,
          y = my_scale(x = newx, xm = xm,
                       xsd = xsd),
          sigma = sigma,
          l = l
        )
        preds <- drop(crossprod(K_star, mat_coefs)) + ym
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    if (fit_method == "polynomial")
    {
      y <- fit_obj$y
      xm <- fit_obj$xm
      ym <- fit_obj$ym
      xsd <- fit_obj$xsd
      lags <- fit_obj$lags
      xreg <- fit_obj$xreg
      sigma <- fit_obj$sigma
      l <- fit_obj$l
      d <- fit_obj$d
      mat_coefs <- fit_obj$mat_coefs
      class_res <- fit_obj$class_res
      series_names <- fit_obj$series_names

      for (i in 1:h)
      {
        # enlever l'intercept
        newx <- reformat_cpp(y, lags)
        K_star <-
          poly_kxy_cpp(
            x = xreg,
            y = my_scale(x = newx, xm = xm,
                         xsd = xsd),
            sigma = sigma,
            d = d,
            l = l
          )
        preds <- drop(crossprod(K_star, mat_coefs)) + ym
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    if (fit_method == "VAR")
    {
      y <- fit_obj$y
      p <- length(fit_obj$series_names)
      if (type_ci == "none")
      {
        preds <- predict(fit_obj$fit_obj, n.ahead = h)
        res <- sapply(1:p,
                      function(i)
                        preds$fcst[[i]][, 1])
        colnames(res) <- fit_obj$series_names
      } else {
        type_ci <- "VAR"
        res <- predict(fit_obj$fit_obj, n.ahead = h,
                       ci = level / 100)
      }
    }

    if (fit_method == "lassoVAR")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      ym <- fit_obj$ym
      xm <- fit_obj$xm
      xsd <- fit_obj$scales
      p <- length(fit_obj$series_names)

      for (i in 1:h)
      {
        newinfo <- my_scale(y[1:(lags + 1), ],
                            xm = xm, xsd = xsd)
        newx <- reformat_cpp(newinfo, lags)
        # consider using parSapply when p is high
        preds <- sapply(1:p, function(i) newx%*%fit_obj$fit_lasso[[i]]$beta) + ym
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    if (fit_method == "scn")
    {
      y <- fit_obj$y
      lags <- fit_obj$lags

      for (i in 1:h)
      {
        newx <- reformat_cpp(y, lags)
        preds <- predict_SCN(fit_obj, newx = newx)
        y <- rbind_vecmat_cpp(preds, y)
      }
    }

    if (fit_method == "pls")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ
      ym <- fit_obj$ym

      if (fit_obj$direct_link == TRUE)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

             if(fit_obj$hidden_layer_bias == FALSE)
             {
                 newx <- cbind(newx, matrix(g(my_scale(
                   newx, xm = nn_xm,
                   xsd = nn_xsd
                 ) %*% w)[, hidden_layer_index], nrow = 1))
            } else {
                newx <- cbind(newx, matrix(g(my_scale(
                  cbind(1, newx), xm = c(1, nn_xm),
                  xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
            }

          preds <- predict_pls(fit_obj, newx = newx) + ym
          y <- rbind_vecmat_cpp(preds, y)
          }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1)
          } else {
            newx <- matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1)
          }

          preds <- predict_pls(fit_obj, newx = newx) + ym
          y <- rbind_vecmat_cpp(preds, y)
        }
      }
    }

    if (fit_method == "pcr")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      ym <- fit_obj$ym
      ncomp <- fit_obj$ncomp
      p <- length(fit_obj$series_names)

        for (i in 1:h)
        {
          newx <- my_scale(reformat_cpp(y, lags),
                           xm = fit_obj$xm, xsd = fit_obj$scales)
          preds <- sapply(1:p, function(i) predict(fit_obj$fit_obj[[i]],
                        newdata = newx)[, , ncomp]) + ym
          y <- rbind_vecmat_cpp(preds, y)
        }
    }
  }

  # 2 - direct forecasts -------------------------------------------------

  # direct forecasts
  if (type_forecast == "direct")
  {
    if (fit_method == "ridge2")
    {
      # if col_sample < 1
      if(is.numeric(col_index) == FALSE)
      {
        # observed values (minus lagged)
        y <- fit_obj$y
        # observed values (minus selected columns)
        lags <- fit_obj$lags
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd
        w <- fit_obj$w
        g <- fit_obj$activ


        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))

          } else {
            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }

          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_ridge_mts(x = newtrainingx, lags = fit_obj$lags,
                                   nb_hidden = fit_obj$nb_hidden, fit_method = "ridge2",
                                   nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                   hidden_layer_bias = fit_obj$hidden_layer_bias,
                                   col_sample = fit_obj$col_sample, row_sample = fit_obj$row_sample,
                                   a = fit_obj$a, lambda_1 = fit_obj$lambda_1, lambda_2 = fit_obj$lambda_2,
                                   seed = fit_obj$seed)
        }
      } else {
        # observed values (minus lagged)
        y <- fit_obj$y
        # observed values (minus selected columns)
        y_reduced <- as.matrix(fit_obj$y[ , col_index])
        lags <- fit_obj$lags
        nn_xm <- fit_obj$nn_xm
        nn_xsd <- fit_obj$nn_xsd

        w <- fit_obj$w
        g <- fit_obj$activ

        for (i in 1:h)
        {
          newx <- reformat_cpp(y_reduced, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))
          } else {
            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }

          preds <- predict_myridge(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_ridge_mts(x = newtrainingx, lags = fit_obj$lags,
                                   nb_hidden = fit_obj$nb_hidden, fit_method = "ridge2",
                                   nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                   hidden_layer_bias = fit_obj$hidden_layer_bias,
                                   col_sample = fit_obj$col_sample, row_sample = fit_obj$row_sample,
                                   a = fit_obj$a, lambda_1 = fit_obj$lambda_1, lambda_2 = fit_obj$lambda_2,
                                   seed = fit_obj$seed)
          y_reduced <- as.matrix(y[ , col_index])
        }
      }
    }

    if (fit_method == "lassoVAR")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      ym <- fit_obj$ym
      xm <- fit_obj$xm
      xsd <- fit_obj$scales
      lambda <- fit_obj$lambda
      p <- length(fit_obj$series_names)

      for (i in 1:h)
      {
        newinfo <- my_scale(y[1:(lags + 1), ],
                            xm = xm, xsd = xsd)
        newx <- reformat_cpp(newinfo, lags)
        preds <- sapply(1:p, function(i) newx%*%fit_obj$fit_lasso[[i]]$beta) + ym
        y <- rbind_vecmat_cpp(preds, y)
        newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
        fit_obj <- fit_var_mts(x = newtrainingx, penalization = "l1",
                               lags = lags, lambda = lambda)
      }
    }

    if (fit_method == "scn")
    {
      y <- fit_obj$y
      lags <- fit_obj$lags

      if (fit_obj$method == "greedy")
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)
          preds <- predict_SCN(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_scn_mts(x = newtrainingx, method = "greedy",
                                 lags = fit_obj$lags, activ = fit_obj$activ,
                                 hidden_layer_bias = fit_obj$hidden_layer_bias,
                                 B = ncol(fit_obj$betas_opt),
                                 nu = fit_obj$nu, lam = fit_obj$lam, r = fit_obj$r,
                                 tol = fit_obj$tol, col_sample = fit_obj$col_sample,
                                 verbose = FALSE,
                                 type_optim = fit_obj$type_optim)
        }
      }

      if (fit_obj$method == "direct")
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)
          preds <- predict_SCN(fit_obj, newx = newx)
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_scn_mts(x = newtrainingx, method = "direct",
                                 lags = fit_obj$lags, activ = fit_obj$activ,
                                 hidden_layer_bias = fit_obj$hidden_layer_bias,
                                 B = ncol(fit_obj$betas_opt),
                                 nu = fit_obj$nu, lam = fit_obj$lam, r = fit_obj$r,
                                 tol = fit_obj$tol, col_sample = fit_obj$col_sample,
                                  verbose = FALSE,
                                 type_optim = fit_obj$type_optim)
        }
      }
    }

    if (fit_method == "pls")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ
      ym <- fit_obj$ym

      if (fit_obj$direct_link == TRUE)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))
          } else {
            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
          }
          preds <- predict_pls(fit_obj, newx = newx) + ym
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_pls_mts(x = newtrainingx,
                                 lags = fit_obj$lags,
                                 activ = fit_obj$activ_name,
                                 hidden_layer_bias = fit_obj$hidden_layer_bias,
                                 nb_hidden = fit_obj$nb_hidden,
                                 nodes_sim = fit_obj$method,
                                 direct_link = TRUE,
                                 B = fit_obj$ncomp)
        }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1)
          } else {
            newx <- matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1)
          }

          preds <- predict_pls(fit_obj, newx = newx) + ym
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_pls_mts(x = newtrainingx,
                                 lags = fit_obj$lags,
                                 activ = fit_obj$activ_name,
                                 hidden_layer_bias = fit_obj$hidden_layer_bias,
                                 nb_hidden = fit_obj$nb_hidden,
                                 nodes_sim = fit_obj$method,
                                 direct_link = FALSE,
                                 B = fit_obj$ncomp)
        }
      }
    }

    if (fit_method == "pcr")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      ym <- fit_obj$ym
      ncomp <- fit_obj$ncomp
      p <- length(fit_obj$series_names)

      for (i in 1:h)
      {
        newx <- my_scale(reformat_cpp(y, lags),
                         xm = fit_obj$xm, xsd = fit_obj$scales)
        preds <- sapply(1:p, function(i) predict(fit_obj$fit_obj[[i]],
                                                 newdata = newx)[, , ncomp]) + ym
        y <- rbind_vecmat_cpp(preds, y)
        newtrainingx <- rbind(fit_obj$x, preds)[-1, ]
        fit_obj <- fit_pcr_mts(x = newtrainingx,
                               lags = lags, ncomp = ncomp)
      }
    }

    if (fit_method == "glmboost")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ
      fit_glmboost <- fit_obj$fit_obj
      nb_series <- fit_obj$nb_series
      ym <- fit_obj$ym
      xm <- fit_obj$xm
      scales <- fit_obj$scales
      direct_link <- fit_obj$direct_link

      if (direct_link == TRUE)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

            if(fit_obj$hidden_layer_bias == FALSE)
            {
              newx <- cbind(newx, matrix(g(my_scale(
                newx, xm = nn_xm,
                xsd = nn_xsd
              ) %*% w)[, hidden_layer_index], nrow = 1))
            } else {
              newx <- cbind(newx, matrix(g(my_scale(
                cbind(1, newx), xm = c(1, nn_xm),
                xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))
            }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          #cat("i = ", i, "\n")
          #cat("preds = ", preds, "\n")
          y <- rbind_vecmat_cpp(preds, y)
          #cat("y = ", y, "\n")
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_glmboost_mts(x = newtrainingx, B = fit_obj$B, eta = fit_obj$eta,
                                      lags = fit_obj$lags, nb_hidden = fit_obj$nb_hidden,
                                      nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                      hidden_layer_bias = fit_obj$hidden_layer_bias,
                                      direct_link = TRUE, a = fit_obj$a, seed = fit_obj$seed)
        }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

            if(fit_obj$hidden_layer_bias == FALSE)
            {
              newx <- matrix(g(my_scale(
                newx, xm = nn_xm,
                xsd = nn_xsd
              ) %*% w)[, hidden_layer_index], nrow = 1)
            } else {
              newx <- matrix(g(my_scale(
                cbind(1, newx), xm = c(1, nn_xm),
                xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1)
            }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_glmboost_mts(x = newtrainingx, B = fit_obj$B, eta = fit_obj$eta,
                                      lags = fit_obj$lags, nb_hidden = fit_obj$nb_hidden,
                                      nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                      hidden_layer_bias = fit_obj$hidden_layer_bias,
                                      direct_link = FALSE, a = fit_obj$a, seed = fit_obj$seed)
        }
      }
    }

    if (fit_method == "xgboost")
    {
      # observed values (minus lagged)
      y <- fit_obj$y
      # observed values (minus selected columns)
      lags <- fit_obj$lags
      nn_xm <- fit_obj$nn_xm
      nn_xsd <- fit_obj$nn_xsd
      w <- fit_obj$w
      g <- fit_obj$activ
      fit_xgboost <- fit_obj$fit_obj
      nb_series <- fit_obj$nb_series
      ym <- fit_obj$ym
      xm <- fit_obj$xm
      scales <- fit_obj$scales
      direct_link <- fit_obj$direct_link

      if (direct_link == TRUE)
      {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- cbind(newx, matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1))
          } else {
            newx <- cbind(newx, matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1))

          }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x

          fit_obj <- fit_xgboost_mts(x = newtrainingx, B = fit_obj$B, eta = fit_obj$eta,
                                     lambda = fit_obj$lambda, alpha = fit_obj$alpha,
                                     hidden_layer_bias = fit_obj$hidden_layer_bias,
                                      lags = fit_obj$lags, nb_hidden = fit_obj$nb_hidden,
                                      nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                      direct_link = TRUE, a = fit_obj$a, seed = fit_obj$seed)
        }
      } else {
        for (i in 1:h)
        {
          newx <- reformat_cpp(y, lags)

          if(fit_obj$hidden_layer_bias == FALSE)
          {
            newx <- matrix(g(my_scale(
              newx, xm = nn_xm,
              xsd = nn_xsd
            ) %*% w)[, hidden_layer_index], nrow = 1)
          } else {
            newx <- matrix(g(my_scale(
              cbind(1, newx), xm = c(1, nn_xm),
              xsd = c(1, nn_xsd)) %*% w)[, hidden_layer_index], nrow = 1)
          }

          newx <- my_scale(x = newx, xm = xm, xsd = scales)
          preds <- ym + sapply(1:nb_series,
                               function (i) predict(fit_obj$fit_obj[[i]],
                                                    newdata = newx))
          y <- rbind_vecmat_cpp(preds, y)
          newtrainingx <- rbind(fit_obj$x, preds)[-1, ] # same window length as x
          fit_obj <- fit_xgboost_mts(x = newtrainingx, B = fit_obj$B, eta = fit_obj$eta,
                                     hidden_layer_bias = fit_obj$hidden_layer_bias,
                                      lags = fit_obj$lags, nb_hidden = fit_obj$nb_hidden,
                                      nodes_sim = fit_obj$method, activ = fit_obj$activ_name,
                                      direct_link = FALSE, a = fit_obj$a, seed = fit_obj$seed)
        }
      }
    }
  }

  # Forecast method
  # if no confidence interval is required
  if (type_ci == "none")
  {
    if (fit_method == "VAR") {
      res2 <- rbind(as.matrix(y), res)
      n <- nrow(res2)
      res <- res2[(n - h + 1):n,]
      colnames(res) <- fit_obj$series_names
      return(res)
    } else {
      res2 <- rev_matrix_cpp(y)
      n <- nrow(res2)
      res <- res2[(n - h + 1):n,]
      colnames(res) <- fit_obj$series_names
      return(res)
    }
  } else {
    # if a confidence interval is required
    if (fit_method == "VAR") {
      res2 <- rbind(as.matrix(y),
                    sapply(1:p,
                           function(i)
                             res$fcst[[i]][, 1]))

      ans <- lapply(1:p,
                    function(i) {
                      resid_fcast_matrix <- res$fcst[[i]][, c(1, 2, 3)]
                      colnames(resid_fcast_matrix) <-
                        c("Point Forecast",
                          paste0("Lo ", level),
                          paste0("Hi ", level))
                      resid_fcast_matrix
                    })
      names(ans)  <- fit_obj$series_names

      return(list(obs = res2, fcst = ans))
    } else {
      res2 <- rev_matrix_cpp(y)
      n <- nrow(res2)
      res <- res2[(n - h + 1):n,]
      colnames(res) <- fit_obj$series_names
      p <- length(fit_obj$series_names)
      ans <- vector("list", p)
      names(ans) <- fit_obj$series_names
      resid_fcast <- switch(
        type_ci,
        "arima" = lapply(1:p,
                         function(x)
                           forecast::forecast(
                             forecast::auto.arima(fit_obj$resid[, x]),
                             h = h,
                             level = level
                           )),
        "ets" = lapply(1:p,
                       function(x)
                         forecast::forecast(
                           forecast::ets(fit_obj$resid[, x]),
                           h = h,
                           level = level
                         )),
        "theta" = lapply(1:p,
                         function(x)
                           forecast::thetaf(fit_obj$resid[, x],
                                            h = h, level = level))
      )

      nb <- 1 + 2 * length(level)

      for (i in 1:p)
      {
        resid_fcast_matrix <- cbind(0, resid_fcast[[i]]$lower,
                                    resid_fcast[[i]]$upper)
        colnames(resid_fcast_matrix) <-
          c("Point Forecast",
            paste0("Lo ", level),
            paste0("Hi ", level))
        ans[[i]] <- resid_fcast_matrix + matrix(rep(res[, i], nb),
                                                ncol = nb, byrow = FALSE)
      }

      colnames(res2) <- fit_obj$series_names
      return(list(obs = res2, fcst = ans))
    }
  }

}
