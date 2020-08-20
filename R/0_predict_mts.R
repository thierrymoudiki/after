# ridge regression
predict_myridge <- function(fit_obj, newx)
{
  my_scale(x = newx, xm = fit_obj$xm,
           xsd = fit_obj$scales)%*%fit_obj$coef + fit_obj$ym
}

# scn
predict_SCN <- function(fit_obj, newx)
{
  # if a bias is used in the hidden layers
  hidden_layer_bias <- fit_obj$hidden_layer_bias
  # maxL
  maxL <- max(1, ncol(fit_obj$betas_opt))

  # columns' shifting when bias term is (not) included
  col_shift <- 0

    # fit_obj$ contains ym, matrix_betasL_opt, matrix_w_opt, matrix_b_opt, nu, activation
    if(is.vector(newx))
    {
      newx_scaled <- my_scale(x = t(newx), xm = fit_obj$xm,
                                              xsd = fit_obj$xsd)
      # initial fit
      fitted_xL <- fit_obj$ym
    } else {
      newx_scaled <- my_scale(x = newx, xm = fit_obj$xm,
                                              xsd = fit_obj$xsd)
      # initial fit
      fitted_xL <- tcrossprod(rep(1, nrow(newx)),
                              fit_obj$ym)
    }

  if (fit_obj$col_sample < 1) { # if columns' subsampling is used

    if (is.vector(newx_scaled))
    {
      newx_scaled <- t(newx_scaled)
    }

        # not all the boosting iterations, but the ones before early stopping
        for (L in 1:maxL)
        {
          if (hidden_layer_bias == FALSE)
          {
            xreg_scaled <- newx_scaled[, fit_obj$col_sample_indices[, L]]
          } else {
            if(dim(newx_scaled)[1] == 1)
            {
              xreg_scaled <- c(1, newx_scaled[, fit_obj$col_sample_indices[, L]])
            } else {
              xreg_scaled <- cbind(1, newx_scaled[, fit_obj$col_sample_indices[, L]])
            }
          }

          if (is.vector(xreg_scaled))
          {
            xreg_scaled <- t(xreg_scaled)
          } else {
            xreg_scaled <- matrix(xreg_scaled,
                                  nrow = nrow(fitted_xL))
          }

          fitted_xL <- fitted_xL + calculate_fittedeL(betasL = fit_obj$betas_opt[, L],
                                                hL = calculate_hL(x = xreg_scaled,
                                                                  w = as.vector(fit_obj$ws_opt[, L]),
                                                                  activation = fit_obj$activ),
                                                nu = fit_obj$nu)
        }

  } else { # if columns' subsampling is not used

    if(hidden_layer_bias == TRUE)#here
    {
      newx_scaled <- cbind(1, newx_scaled)
    }

      # not all the boosting iterations, but the ones before early stopping
      for (L in 1:max(1, ncol(fit_obj$betas_opt)))
      {
        fitted_xL <- fitted_xL + calculate_fittedeL(betasL = fit_obj$betas_opt[, L],
                                                            hL = calculate_hL(x = newx_scaled,
                                                                              w = as.vector(fit_obj$ws_opt[, L]),
                                                                              activation = fit_obj$activ),
                                                            nu = fit_obj$nu)
      }
  }

  return(fitted_xL)
}

# partial least squares
predict_pls <- function(fit_obj, newx)
{
  ncomp <- fit_obj$ncomp

   #inspired from caret::predict.PLS
   if(is.vector(newx)) newx <- t(newx)

  # from coef.mvr in pls package
  B <- fit_obj$coefficients[, , 1:ncomp, drop = FALSE]
  dB <- dim(B)
  dB[1] <- dB[1] + 1

  BInt <- array(dim = dB)
  BInt[-1, , ] <- B
  for (i in seq(along = 1:ncomp)) BInt[1, , i] <- fit_obj$Ymeans - fit_obj$Xmeans %*% B[, , i]
  B <- BInt
  # stop

  mat_coeffs <- B[-1, , ncomp]

  # from predict.mvr in pls package
  #return(sweep(newx[, 1:nrow(mat_coeffs)] %*% mat_coeffs, 2, B[1, , ncomp], "+"))
  return(sweep(newx %*% mat_coeffs, 2, B[1, , ncomp], "+"))
}

# pcr
predict_pcr <- function(fit_obj, newx)
{
  pls::predict(my_scale(x = newx, xm = fit_obj$xm,
                           xsd = fit_obj$scales) + fit_obj$ym)
}
