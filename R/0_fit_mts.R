# 0 - individual models' fitting (not time series) -------------------------------------------

# scn -----
fit_SCN <- function(x,
                    y,
                    B = 10,
                    nu = 0.1,
                    col_sample = 1,
                    seed = 123,
                    lam = 0.1,
                    r = 0.3,
                    tol = 1e-10,
                    type_optim = c("nlminb", "nmkb"),
                    activation = c("sigmoid", "tanh"),
                    method = c("greedy", "direct"),
                    hidden_layer_bias = FALSE,
                    verbose = FALSE)
{
  if (is.ts(x))
  {
    freq_x <- frequency(x)
    start_preds <- tsp(x)[2] + 1 / freq_x
  } else {
    freq_x <- NULL
    start_preds <- NULL
  }

  stopifnot(nu > 0 && nu <= 1)
  stopifnot(r > 0 && r < 1)
  stopifnot(B > 1)
  stopifnot(col_sample >= 0.5 || col_sample <= 1) #must be &&
  d <- ncol(x)
  # d_reduced <- 0 # for col_sample < 1
  # dd <- 0 # for hidden_layer_bias = TRUE
  # dd_reduced <- 0 # for col_sample < 1 && for hidden_layer_bias = TRUE

  if (hidden_layer_bias == TRUE)
  {
    dd <- d + 1
  }
  N <- nrow(x)
  m <- ncol(y)
  stopifnot(nrow(y) == nrow(x))
  ym <- colMeans(y)
  xscales <- my_scale(x)
  xm <- xscales$xm
  xsd <- xscales$xsd
  x_scaled <- xscales$res
  centered_y <- my_scale(x = y, xm = ym)
  type_optim <- match.arg(type_optim)
  L <- NULL
  col_sample_indices <- NULL
  d_reduced <- NULL
  # choice between Algo SC-I = "greedy" and Algo SC-III = "direct"
  # from Wang et Li (2017)
  method <- match.arg(method)
  # choice of activation function
  activation <- match.arg(activation)

  # columns' subsampling
  if (col_sample < 1)
  {
    d_reduced <- max(1, floor(col_sample * d))
    if (hidden_layer_bias == TRUE)
    {
      dd_reduced <- d_reduced + 1
    }

    if (d_reduced > 1)
    {
      set.seed(seed)
      col_sample_indices <-
        sapply(1:B, function (i)
          sort(sample.int(n = d,
                          size = d_reduced)))
    } else {
      set.seed(seed)
      col_sample_indices <-
        t(sapply(1:B, function (i)
          sort(sample.int(
            n = d,
            size = 1
          ))))
    }
  }

  # columns' names
  names_L <- paste0("L", 1:B)
  names_m <- paste0("m", 1:m)
  names_d <- paste0("d", 1:d)
  names_N <- paste0("N", 1:N)
  if (col_sample < 1)
  {
    names_d_reduced <- paste0("d", 1:d_reduced)
    if (hidden_layer_bias == TRUE)
      names_dd_reduced <- paste0("d", 1:dd_reduced)
  }

  matrix_betas_opt <- matrix(0, nrow = m, ncol = B)
  colnames(matrix_betas_opt) <- names_L
  rownames(matrix_betas_opt) <- names_m

  if (method == "direct")
  {
    matrix_hL_opt <- matrix(0, nrow = N, ncol = B)
    colnames(matrix_hL_opt) <- names_L
    rownames(matrix_hL_opt) <- names_N
  }

  if (col_sample < 1)
  {
    if (hidden_layer_bias == FALSE)
    {
      matrix_ws_opt <- matrix(0, nrow = d_reduced, ncol = B)
      colnames(matrix_ws_opt) <- names_L
      rownames(matrix_ws_opt) <- names_d_reduced
    } else {
      matrix_ws_opt <- matrix(0, nrow = dd_reduced, ncol = B)
      #colnames(matrix_ws_opt) <- names_L
      #rownames(matrix_ws_opt) <- names_d_reduced
    }
  } else {
    if (hidden_layer_bias == FALSE)
    {
      matrix_ws_opt <- matrix(0, nrow = d, ncol = B)
      colnames(matrix_ws_opt) <- names_L
      rownames(matrix_ws_opt) <- names_d
    } else {
      matrix_ws_opt <- matrix(0, nrow = dd, ncol = B)
      #colnames(matrix_ws_opt) <- names_L
      #rownames(matrix_ws_opt) <- names_d
    }
  }

  # beginning of the algorithm
  current_error <- centered_y
  colnames(current_error) <- names_m
  rownames(current_error) <- names_N
  current_error_norm <- norm(current_error, type = "F")

  # Main boosting loop
  L <- 1
  if (col_sample < 1)
    # columns' subsampling
  {
    if (hidden_layer_bias == FALSE)
    {
      # 0 - 1 - col_sample < 1 && hidden_layer_bias == FALSE

      # inequality objective function
      InequalityOF <- function(w) {
        # calculate hL
        # calculate xsi = (xsi_1, ..., xsi_m)
        xsi_vec <- calculate_xsiL(
          eL = current_error,
          hL = calculate_hL(
            x = as.matrix(x_scaled[, col_sample_indices[, L]]),
            w = w,
            activation = activation
          ),
          nu = nu,
          r = r,
          L = L
        )
        # calculate xsiL = sum(xsi)
        # return -xsiL*(min(xsi) > 0)
        return(-sum(xsi_vec) * (min(xsi_vec) > 0))
      }
      InequalityOF <- compiler::cmpfun(InequalityOF)

    } else {
      # hidden_layer_bias == TRUE

      # inequality objective function
      InequalityOF <- function(w) {
        # calculate hL
        # calculate xsi = (xsi_1, ..., xsi_m)
        xsi_vec <- calculate_xsiL(
          eL = current_error,
          hL = calculate_hL(
            x = as.matrix(cbind(1, x_scaled[, col_sample_indices[, L]])),
            w = w,
            activation = activation
          ),
          nu = nu,
          r = r,
          L = L
        )
        # calculate xsiL = sum(xsi)
        # return -xsiL*(min(xsi) > 0)
        return(-sum(xsi_vec) * (min(xsi_vec) > 0))
      }
      InequalityOF <- compiler::cmpfun(InequalityOF)

    }

    while (L <= B && current_error_norm > tol) {
      if (verbose)
      {
        cat("L = ", L, "\n")
        cat("\n")
      }

      if (hidden_layer_bias == FALSE)
      {
        lower <- rep(-lam, d_reduced)
        upper <- rep(lam, d_reduced)
      } else {
        lower <- rep(-lam, dd_reduced)
        upper <- rep(lam, dd_reduced)
      }

      if (type_optim == "nlminb")
      {
        set.seed(L)
        out_opt <-
          stats::nlminb(
            start = lower + (upper - lower) * runif(length(lower)),
            objective = InequalityOF,
            lower = lower,
            upper = upper
          )
      }

      if (type_optim == "nmkb")
      {
        if (length(lower) > 1)
        {
          set.seed(L)
          out_opt <-
            dfoptim::nmkb(
              par = lower + (upper - lower) * runif(length(lower)),
              fn = InequalityOF,
              lower = lower,
              upper = upper
            )
        } else {
          set.seed(L)
          out_opt <-
            suppressWarnings(
              stats::optim(
                par = lower + (upper - lower) * runif(length(lower)),
                fn = InequalityOF,
                method = "Nelder-Mead"
              )
            )
        }
      }

      w_opt <- out_opt$par
      matrix_ws_opt[, L] <- w_opt

      if (verbose)
      {
        if (hidden_layer_bias == FALSE)
        {
          names(w_opt) <- paste0("w", 1:d_reduced)
        } else {
          names(w_opt) <- paste0("w", 1:dd_reduced)
        }
        cat("w_opt", "\n")
        print(w_opt)
        cat("\n")
      }

      # calculate hL_opt
      if (hidden_layer_bias == FALSE)
      {
        hL_opt <-
          calculate_hL(x = as.matrix(x_scaled[, col_sample_indices[, L]]),
                       w = w_opt,
                       activation = activation)
      } else {
        hL_opt <-
          calculate_hL(
            x = cbind(1, as.matrix(x_scaled[, col_sample_indices[, L]])),
            w = w_opt,
            activation = activation
          )
      }

      # calculate betaL_opt with Algo SC-I
      if (method == "greedy")
      {
        betaL_opt <- calculate_betasL(current_error, hL_opt)
        matrix_betas_opt[, L] <- betaL_opt
        # update the error
        current_error <-
          current_error - calculate_fittedeL(betasL = betaL_opt,
                                             hL = hL_opt,
                                             nu = nu)
      }

      # calculate betaL_opt with Algo SC-III
      if (method == "direct")
      {
        # cat("----- L = ", "\n")
        # print(L)
        # cat("\n")
        matrix_hL_opt[, L] <- hL_opt
        # cat("matrix_hL_opt", "\n")
        # print(matrix_hL_opt)
        # cat("\n")
        betaL_opt <- .lm.fit(x = as.matrix(matrix_hL_opt[, 1:L]),
                             y = centered_y)$coef
        # cat("beta", "\n")
        # print(betaL_opt)
        # cat("\n")
        matrix_betas_opt[, L] <- betaL_opt[L,]
        # update the error
        current_error <-
          current_error - calculate_fittedeL(
            betasL = as.vector(matrix_betas_opt[, L]),
            hL = hL_opt,
            nu = nu
          )
      }

      # update the norm of the error matrix
      current_error_norm <- norm(current_error, type = "F")

      L <- L + 1
    } # end while(L <= B && current_error_norm > tol) for col_sample < 1

  } else {
    # if col_sample == 1 (no subsampling of the columns)

    # 0-2 - col_sample == 1 && hidden_layer_bias == FALSE

    # inequality objective function
    if (hidden_layer_bias == FALSE)
    {
      InequalityOF <- function(w) {
        # calculate hL
        # calculate xsi = (xsi_1, ..., xsi_m)
        xsi_vec <- calculate_xsiL(
          eL = current_error,
          hL = calculate_hL(
            x = x_scaled,
            w = w,
            activation = activation
          ),
          nu = nu,
          r = r,
          L = L
        )
        # calculate xsiL = sum(xsi)
        # return -xsiL*(min(xsi) > 0)
        return(-sum(xsi_vec) * (min(xsi_vec) > 0))
      }
      InequalityOF <- compiler::cmpfun(InequalityOF)

    } else {
      # hidden_layer_bias == TRUE

      InequalityOF <- function(w) {
        # calculate hL
        # calculate xsi = (xsi_1, ..., xsi_m)
        xsi_vec <- calculate_xsiL(
          eL = current_error,
          hL = calculate_hL(
            x = cbind(1, x_scaled),
            w = w,
            activation = activation
          ),
          nu = nu,
          r = r,
          L = L
        )
        # calculate xsiL = sum(xsi)
        # return -xsiL*(min(xsi) > 0)
        return(-sum(xsi_vec) * (min(xsi_vec) > 0))
      }
      InequalityOF <- compiler::cmpfun(InequalityOF)

    }

    while (L <= B && current_error_norm > tol) {
      if (verbose)
      {
        cat("L = ", L, "\n")
        cat("\n")
      }

      if (hidden_layer_bias == FALSE)
      {
        lower <- rep(-lam, d)
        upper <- rep(lam, d)

        if (type_optim == "nlminb")
        {
          set.seed(L)
          out_opt <-
            stats::nlminb(
              start = lower + (upper - lower) * runif(length(lower)),
              objective = InequalityOF,
              lower = lower,
              upper = upper
            )
        }

        if (type_optim == "nmkb")
        {
          if (length(lower) > 1)
          {
            set.seed(L)
            out_opt <-
              dfoptim::nmkb(
                par = lower + (upper - lower) * runif(length(lower)),
                fn = InequalityOF,
                lower = lower,
                upper = upper
              )
          } else {
            set.seed(L)
            out_opt <-
              suppressWarnings(
                stats::optim(
                  par = lower + (upper - lower) * runif(length(lower)),
                  fn = InequalityOF,
                  method = "Nelder-Mead"
                )
              )
          }
        }

      } else {
        # hidden_layer_bias == TRUE

        lower <- rep(-lam, dd)
        upper <- rep(lam, dd)

        if (type_optim == "nlminb")
        {
          set.seed(L)
          out_opt <-
            stats::nlminb(
              start = lower + (upper - lower) * runif(length(lower)),
              # dd <- d + 1
              objective = InequalityOF,
              lower = lower,
              upper = upper
            )
        }

        if (type_optim == "nmkb")
        {
          if (length(lower) > 1)
          {
            set.seed(L)
            out_opt <-
              dfoptim::nmkb(
                par = lower + (upper - lower) * runif(length(lower)),
                # dd <- d + 1
                fn = InequalityOF,
                lower = lower,
                upper = upper
              )
          } else {
            set.seed(L)
            out_opt <-
              suppressWarnings(
                stats::optim(
                  par = lower + (upper - lower) * runif(length(lower)),
                  fn = InequalityOF,
                  method = "Nelder-Mead"
                )
              )
          }
        }

      }

      w_opt <- out_opt$par
      matrix_ws_opt[, L] <- w_opt
      if (verbose)
      {
        names(w_opt) <- paste0("w", 1:max(d, dd))
        cat("w_opt", "\n")
        print(w_opt)
        cat("\n")
      }

      # calculate hL_opt
      if (hidden_layer_bias == FALSE)
      {
        hL_opt <- calculate_hL(x = x_scaled,
                               w = w_opt,
                               activation = activation)
      } else {
        hL_opt <- calculate_hL(x = cbind(1, x_scaled),
                               w = w_opt,
                               activation = activation)
      }

      # calculate betaL_opt
      if (method == "greedy")
      {
        betaL_opt <- calculate_betasL(current_error, hL_opt)
        matrix_betas_opt[, L] <- betaL_opt
        # update the error
        current_error <-
          current_error - calculate_fittedeL(betasL = betaL_opt,
                                             hL = hL_opt,
                                             nu = nu)
      }

      if (method == "direct")
      {
        # cat("----- L = ", "\n")
        # print(L)
        # cat("\n")
        matrix_hL_opt[, L] <- hL_opt
        # cat("matrix_hL_opt", "\n")
        # print(matrix_hL_opt)
        # cat("\n")
        betaL_opt <- .lm.fit(x = as.matrix(matrix_hL_opt[, 1:L]),
                             y = centered_y)$coef
        # cat("beta", "\n")
        # print(betaL_opt)
        # cat("\n")
        matrix_betas_opt[, L] <- betaL_opt[L,]

        # update the error
        current_error <-
          current_error - calculate_fittedeL(
            betasL = as.vector(matrix_betas_opt[, L]),
            hL = hL_opt,
            nu = nu
          )
      }

      # update the norm of the error matrix
      current_error_norm <- norm(current_error, type = "F")

      L <- L + 1
    } # end while(L <= B && current_error_norm > tol) for col_sample == 1

  } # end main boosting loop

  bool_non_zero_betas <- colSums(matrix_betas_opt) != 0
  names(ym) <- names_m
  names(xm) <- names_d
  names(xsd) <- names_d

  if (!is.null(d_reduced) && d_reduced == 1)
  {
    return(
      list(
        y = y,
        x = x,
        ym = ym,
        xm = xm,
        xsd = xsd,
        col_sample = col_sample,
        betas_opt = matrix_betas_opt[, bool_non_zero_betas],
        ws_opt = t(matrix_ws_opt[, bool_non_zero_betas]),
        type_optim = type_optim,
        col_sample_indices = col_sample_indices,
        activ = activation,
        hidden_layer_bias = hidden_layer_bias,
        nu = nu,
        current_error = current_error,
        current_error_norm = current_error_norm
      )
    )
  } else {
    return(
      list(
        y = y,
        x = x,
        ym = ym,
        xm = xm,
        xsd = xsd,
        col_sample = col_sample,
        betas_opt = as.matrix(matrix_betas_opt[, bool_non_zero_betas]),
        ws_opt = as.matrix(matrix_ws_opt[, bool_non_zero_betas]),
        type_optim = type_optim,
        col_sample_indices = col_sample_indices,
        activ = activation,
        hidden_layer_bias = hidden_layer_bias,
        nu = nu,
        current_error = current_error,
        current_error_norm = current_error_norm
      )
    )
  }

}

# partial least squares -----
fit_pls <- function(x, y, ncomp = 1)
{
  fit_obj_pls <- pls::simpls.fit(
    X = x,
    Y = y,
    ncomp = ncomp,
    stripped = TRUE
  ) # to avoid lengthy calculations
  fit_obj_pls$ncomp <- ncomp



  return(fit_obj_pls)
}

# principal components regression -----
fit_pcr <- function(x, y, ncomp = 1)
{
  nb_series <- ncol(y)

  df_list <- lapply(1:nb_series,
                    function (i)
                      cbind.data.frame(y = y[, i], x))

  fit_obj_pcr <- lapply(1:nb_series,
                        function (i)
                          pls::pcr(y ~ ., data = df_list[[i]],
                                   scale = TRUE))

  fit_obj_pcr$ncomp <- ncomp

  return(fit_obj_pcr)
}


# 1 - individual fitting models -------------------------------------------
# regularized regression models
fit_ridge_mts <- function(x,
                          lags = 1,
                          nb_hidden = 5,
                          fit_method = c("ridge2", "ridge", "mgaussian"),
                          nodes_sim = c("sobol", "halton", "unif"),
                          activ = c("relu", "sigmoid", "tanh",
                                    "leakyrelu", "elu", "linear"),
                          hidden_layer_bias = FALSE,
                          col_sample = 1,
                          row_sample = 1,
                          a = 0.01,
                          lambda = 0.1,
                          alpha = 0.5,
                          lambda_1 = 0.1,
                          lambda_2 = 0.1,
                          seed = 1)
{
  stopifnot(col_sample >= 0 && col_sample <= 1)
  stopifnot(floor(col_sample * ncol(x)) >= 1)
  stopifnot(row_sample >= 0.5 && row_sample <= 1)
  stopifnot(is.wholenumber(nb_hidden))
  fit_method <- match.arg(fit_method)
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  blocks_index <- TRUE
  col_index <- TRUE

  series_names <- colnames(x)

  if (row_sample < 1)
  {
    nb_rows_reduced <- max(4, floor(row_sample * nrow(x)))
    # !!! because the ts object is in reverse order
    set.seed(seed)
    x <- rev_matrix_cpp(x)[1:sample(4:nb_rows_reduced, size = 1),]
  } else {
    # !!! because the ts object is in reverse order
    x <- rev_matrix_cpp(x)
  }

  if (col_sample <= 1)
  {
    nb_cols_reduced <- max(1, floor(col_sample * ncol(x)))
    set.seed(seed)
    col_index <- sort(sample(1:ncol(x), size = nb_cols_reduced))
    blocks_index <- lapply(col_index,
                           function(i)
                             seq(
                               from = (i - 1) * lags + 1,
                               to = i * lags,
                               by = 1
                             ))
    names(blocks_index) <- series_names[col_index]
    blocks_index <- unlist(blocks_index)
  }

  y_x <- create_train_inputs_cpp(x, lags)

  if (nb_hidden > 0)
  {
    list_xreg <- create_new_predictors(
      as.matrix(y_x$xreg[, blocks_index]),
      nb_hidden = nb_hidden,
      hidden_layer_bias = hidden_layer_bias,
      method = nodes_sim,
      activ = activ,
      a = a,
      seed = seed
    )
    xreg <- list_xreg$predictors
  } else {
    xreg <-  y_x$xreg[, TRUE]
  }

  # observed values, minus the lags (!)(beware)(!)
  observed_values <- y <- y_x$y

  if (fit_method == "ridge2")
  {
    stopifnot(nb_hidden > 0)
    stopifnot(lambda_1 > 0 && lambda_2 > 0)

    ym <- colMeans(y)
    centered_y <- my_scale(x = y, xm = ym)
    x_scaled <- my_scale(xreg)
    xreg <- x_scaled$res

    xm <- x_scaled$xm
    xsd <- x_scaled$xsd

    k_p <-
      ifelse(col_sample == 1, lags * ncol(x), lags * nb_cols_reduced)
    index <- 1:k_p

    # original predictors (scaled)
    X <- as.matrix(xreg[, index])

    # transformed predictors (scaled)
    Phi_X <- xreg[,-index]

    B <- crossprod(X) + lambda_1 * diag(k_p)
    C <- crossprod(Phi_X, X)
    D <- crossprod(Phi_X) + lambda_2 * diag(ncol(Phi_X))
    B_inv <- my_ginv(B)
    W <- C %*% B_inv
    S_mat <- D - tcrossprod(W, C)
    S_inv <- my_ginv(S_mat)
    Y <- S_inv %*% W
    inv <- rbind(cbind(B_inv + crossprod(W, Y),-t(Y)),
                 cbind(-Y, S_inv))

    lscoef <- inv %*% crossprod(xreg, centered_y)
    colnames(lscoef) <- series_names
    rownames(lscoef) <- colnames(xreg)

    lsfit <- xreg %*% lscoef
    fitted_values <-
      rev_matrix_cpp(lsfit + matrix(rep(ym, each = nrow(lsfit)),
                                    ncol = ncol(lsfit)))
    resid <- rev_matrix_cpp(observed_values) - fitted_values
    colnames(resid) <- series_names

    return(
      list(
        y = y,
        x = x,
        lags = lags,
        col_index = col_index,
        series_names = series_names,
        class_res = fit_method,
        nb_hidden = nb_hidden,
        method = nodes_sim,
        w = list_xreg$w,
        activ = list_xreg$activ,
        activ_name = activ,
        hidden_layer_bias = hidden_layer_bias,
        col_sample = col_sample,
        row_sample = row_sample,
        a = a,
        lambda_1 = lambda_1,
        lambda_2 = lambda_2,
        hidden_layer_index = list_xreg$hidden_layer_index,
        seed = seed,
        nn_xm = list_xreg$xm,
        nn_xsd = list_xreg$xsd,
        ym = ym,
        xm = xm,
        scales = xsd,
        coef = lscoef,
        resid = resid
      )
    )
  }

  if (fit_method == "ridge")
  {
    stopifnot(lambda > 0)
    ym <- colMeans(y)
    centered_y <- my_scale(x = y, xm = ym)
    x_scaled <- my_scale(xreg)
    xreg <- x_scaled$res

    xm <- x_scaled$xm
    xsd <- x_scaled$xsd

    inv <- my_ginv(crossprod(xreg) + lambda * diag(ncol(xreg)))
    lscoef <- inv %*% crossprod(xreg, centered_y)
    colnames(lscoef) <- series_names
    rownames(lscoef) <- colnames(xreg)

    lsfit <- xreg %*% lscoef
    fitted_values <-
      rev_matrix_cpp(lsfit + matrix(rep(ym, each = nrow(lsfit)),
                                    ncol = ncol(lsfit)))
    resid <- rev_matrix_cpp(observed_values) - fitted_values

    if (nb_hidden > 0)
    {
      return(
        list(
          y = y,
          lags = lags,
          series_names = series_names,
          class_res = fit_method,
          nb_hidden = nb_hidden,
          method = nodes_sim,
          w = list_xreg$w,
          activ = list_xreg$activ,
          hidden_layer_bias = hidden_layer_bias,
          hidden_layer_index = list_xreg$hidden_layer_index,
          seed = seed,
          nn_xm = list_xreg$xm,
          nn_xsd = list_xreg$xsd,
          ym = ym,
          xm = xm,
          scales = xsd,
          coef = lscoef,
          resid = resid
        )
      )
    } else {
      return(
        list(
          y = y,
          lags = lags,
          series_names = series_names,
          class_res = fit_method,
          nb_hidden = 0,
          ym = ym,
          xm = xm,
          scales = xsd,
          coef = lscoef,
          resid = resid
        )
      )
    }
  }

  if (fit_method == "mgaussian")
  {
    stopifnot(lambda > 0)

    fit_glmnet <-
      glmnet::glmnet(
        x = xreg,
        y = y,
        standardize = FALSE,
        intercept = FALSE,
        alpha = alpha,
        family = "mgaussian"
      )

    fitted_values <-
      predict(fit_glmnet, s = lambda, newx = xreg)[, , 1]

    resid <-
      rev_matrix_cpp(observed_values - fitted_values) # looks false (!)

    if (nb_hidden > 0)
    {
      return(
        list(
          y = y,
          lags = lags,
          series_names = series_names,
          class_res = fit_method,
          nb_hidden = nb_hidden,
          method = nodes_sim,
          fit_obj = fit_glmnet,
          s = lambda,
          w = list_xreg$w,
          activ = list_xreg$activ,
          hidden_layer_bias = hidden_layer_bias,
          hidden_layer_index = list_xreg$hidden_layer_index,
          seed = seed,
          nn_xm = list_xreg$xm,
          nn_xsd = list_xreg$xsd,
          resid = resid
        )
      )
    } else {
      # nb_hidden <= 0
      return(
        list(
          y = y,
          lags = lags,
          series_names = series_names,
          class_res = fit_method,
          nb_hidden = 0,
          fit_obj = fit_glmnet,
          s = lambda,
          seed = seed,
          resid = resid
        )
      )
    }
  }

}


# least squares regression  model to multiple time series
fit_lm_mts <- function(x,
                       lags = 1,
                       nb_hidden = 5,
                       nodes_sim = c("sobol", "halton", "unif"),
                       activ = c("relu", "sigmoid", "tanh",
                                 "leakyrelu", "elu", "linear"),
                       hidden_layer_bias = FALSE,
                       col_sample = 1,
                       a = 0.01,
                       seed = 1)
{
  stopifnot(col_sample > 0 || col_sample <= 1)
  stopifnot(is.wholenumber(nb_hidden))

  series_names <- colnames(x)
  # !!! because the ts object is in reverse order
  x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(x, lags)
  observed_values <- y <- y_x$y
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)

  list_xreg <- create_new_predictors(
    y_x$xreg,
    nb_hidden = nb_hidden,
    hidden_layer_bias = hidden_layer_bias,
    method = nodes_sim,
    activ = activ,
    a = a,
    seed = seed
  )
  xreg <- list_xreg$predictors

  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)
  x_scaled <- my_scale(xreg)
  xreg <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  k_p <- lags * ncol(x)
  index <- 1:k_p

  fit_lm <- .lm.fit(x = xreg, y = centered_y)

  lscoef <- fit_lm$coefficients
  colnames(lscoef) <- series_names
  rownames(lscoef) <- colnames(xreg)

  lsfit <- xreg %*% lscoef
  fitted_values <-
    rev_matrix_cpp(lsfit + matrix(rep(ym, each = nrow(lsfit)),
                                  ncol = ncol(lsfit)))
  resid <- rev_matrix_cpp(observed_values) - fitted_values

  if (nb_hidden > 0)
  {
    return(
      list(
        y = y,
        lags = lags,
        series_names = series_names,
        class_res = "lm",
        nb_hidden = nb_hidden,
        method = nodes_sim,
        w = list_xreg$w,
        activ = list_xreg$activ,
        hidden_layer_bias = hidden_layer_bias,
        hidden_layer_index = list_xreg$hidden_layer_index,
        seed = seed,
        nn_xm = list_xreg$xm,
        nn_xsd = list_xreg$xsd,
        ym = ym,
        xm = xm,
        scales = xsd,
        coef = lscoef,
        resid = resid
      )
    )
  } else {
    return(
      list(
        y = y,
        lags = lags,
        series_names = series_names,
        class_res = "lm",
        nb_hidden = nb_hidden,
        method = nodes_sim,
        seed = seed,
        nn_xm = list_xreg$xm,
        nn_xsd = list_xreg$xsd,
        ym = ym,
        xm = xm,
        scales = xsd,
        coef = lscoef,
        resid = resid
      )
    )
  }
}


# Fitting a constrained regression model
# to multiple time series (with > 0 coeffs)
fit_nnls_mts <- function(x,
                         lags = 1,
                         nb_hidden = 5,
                         nodes_sim = c("sobol", "halton", "unif"),
                         activ = c("relu", "sigmoid", "tanh",
                                   "leakyrelu", "elu", "linear"),
                         hidden_layer_bias = FALSE,
                         col_sample = 1,
                         a = 0.01,
                         seed = 1)
{
  stopifnot(col_sample > 0 || col_sample <= 1)
  stopifnot(is.wholenumber(nb_hidden))

  series_names <- colnames(x)
  # !!! because the ts object is in reverse order
  x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(x, lags)
  observed_values <- y <- y_x$y
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)

  list_xreg <- create_new_predictors(
    y_x$xreg,
    nb_hidden = nb_hidden,
    method = nodes_sim,
    activ = activ,
    hidden_layer_bias = hidden_layer_bias,
    a = a,
    seed = seed
  )
  xreg <- list_xreg$predictors

  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)
  x_scaled <- my_scale(xreg)
  xreg <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  k_p <- lags * ncol(x)
  index <- 1:k_p

  lscoef <- sapply(1:ncol(centered_y),
                   function (i)
                     coefficients(nnls::nnls(xreg, centered_y[, i])))
  colnames(lscoef) <- series_names
  rownames(lscoef) <- colnames(xreg)

  lsfit <- xreg %*% lscoef
  fitted_values <-
    rev_matrix_cpp(lsfit + matrix(rep(ym, each = nrow(lsfit)),
                                  ncol = ncol(lsfit)))
  resid <- rev_matrix_cpp(observed_values) - fitted_values

  if (nb_hidden > 0)
  {
    return(
      list(
        y = y,
        lags = lags,
        series_names = series_names,
        class_res = "nnls",
        nb_hidden = nb_hidden,
        method = nodes_sim,
        w = list_xreg$w,
        activ = list_xreg$activ,
        hidden_layer_bias = hidden_layer_bias,
        hidden_layer_index = list_xreg$hidden_layer_index,
        seed = seed,
        nn_xm = list_xreg$xm,
        nn_xsd = list_xreg$xsd,
        ym = ym,
        xm = xm,
        scales = xsd,
        coef = lscoef,
        resid = resid
      )
    )
  } else {
    return(
      list(
        y = y,
        lags = lags,
        series_names = series_names,
        class_res = "nnls",
        nb_hidden = nb_hidden,
        method = nodes_sim,
        seed = seed,
        nn_xm = list_xreg$xm,
        nn_xsd = list_xreg$xsd,
        ym = ym,
        xm = xm,
        scales = xsd,
        coef = lscoef,
        resid = resid
      )
    )
  }
}

# Fitting kernel Ridge regression  model to multiple time series
fit_krls_mts <- function(x,
                         lags = 1,
                         lambda_krls = 0.1,
                         l = 0.1,
                         sigma = 2,
                         d = 1,
                         kernel_type = c("gaussian", "polynomial",
                                         "matern32", "matern52"),
                         inv_method = c("chol", "ginv"))
{
  kernel_type <- match.arg(kernel_type)
  inv_method <- match.arg(inv_method)

  series_names <- colnames(x)

  x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(x, lags)
  xreg <- y_x$xreg
  x_scaled <- my_scale(xreg)
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd
  y <-  observed_values <- y_x$y
  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)

  n_k <- nrow(xreg)
  p <- ncol(y)

  # Fit with kernel
  K <- switch(
    kernel_type,
    "gaussian" = gaussian_kxx_cpp(x = xreg,
                                  sigma = sigma,
                                  l = l),
    "polynomial" = poly_kxx_cpp(
      x = xreg,
      sigma = sigma,
      d = d,
      l = l
    ),
    "matern32" = matern32_kxx_cpp(x = xreg,
                                  sigma = sigma,
                                  l = l),
    "matern52" = matern52_kxx_cpp(x = xreg,
                                  sigma = sigma,
                                  l = l)
  )

  mat_coefs <- switch(
    inv_method,
    "chol" = chol2inv(chol(K + lambda_krls * diag(n_k))) %*%
      centered_y,
    "ginv" = my_ginv(K + lambda_krls * diag(n_k)) %*% centered_y
  )

  lsfit <- drop(crossprod(K, mat_coefs))
  fitted_values <-
    rev_matrix_cpp(lsfit  + matrix(rep(ym, each = nrow(lsfit)),
                                   ncol = ncol(lsfit)))

  resid <- rev_matrix_cpp(observed_values) - fitted_values

  class_res <- kernel_type

  if (class_res == "gaussian" ||
      class_res == "matern32" || class_res == "matern52")
  {
    return(
      list(
        y = y,
        series_names = series_names,
        lags = lags,
        xreg = xreg,
        sigma = sigma,
        l = l,
        lambda_krls = lambda_krls,
        mat_coefs = mat_coefs,
        class_res = class_res,
        ym = ym,
        xm = xm,
        xsd = xsd,
        resid = resid
      )
    )
  }

  if (class_res == "polynomial")
  {
    return(
      list(
        y = y,
        series_names = series_names,
        lags = lags,
        xreg = xreg,
        sigma = sigma,
        d = d,
        l = l,
        lambda_krls = lambda_krls,
        mat_coefs = mat_coefs,
        class_res = class_res,
        ym = ym,
        xm = xm,
        xsd = xsd
      )
    )
  }
}

# Fitting a VAR model to multiple time series
fit_var_mts <- function(x,
                        lags = 1,
                        penalization = c("none", "l1"),
                        lambda = 0.1,
                        # for penalization == "l1" only
                        type_VAR = c("const", "trend",
                                     "both", "none"))
  # for penalization == "none" only
{
  series_names <- colnames(x)
  lag_names <- as.vector(outer(paste0("lag", 1:lags, "_"),
                               series_names, FUN = "paste0"))
  penalization <- match.arg(penalization)

  # unrestricted VAR algo from package 'vars'
  if (penalization == "none")
  {
    type_VAR <- match.arg(type_VAR)
    fit_obj <- vars::VAR(y = x, p = lags, type = type_VAR)

    resids <- resid(fit_obj$fit_obj)

    return(
      list(
        y = x[-(1:lags), ],
        fit_obj = fit_obj,
        resid = resids,
        series_names = series_names,
        class_res = "VAR"
      )
    )
  }

  # Fu (1998) algo for coordinate descent
  if (penalization == "l1")
  {
    nb_series <- ncol(x)
    # !!! because the ts object is in reverse order
    rev_x <- rev_matrix_cpp(x)

    y_x <- create_train_inputs_cpp(rev_x, 1)
    xreg <- y_x$xreg
    x_scaled <- my_scale(xreg)
    xm <- x_scaled$xm
    xsd <- x_scaled$xsd
    observed_values <- y <- y_x$y
    ym <- colMeans(y)
    centered_y <- my_scale(x = y, xm = ym)

    #X <- reformat_cpp(x_scaled$res,
    #                  n_k = lags)
    X <- x_scaled$res

    XX <- crossprod(X)
    XX2 <- 2 * XX

    # Xy2 <- lapply(1:nb_series,
    #               function (i)
    #                 2 * X * centered_y[1, i])
    Xy2 <- lapply(1:nb_series,
                  function (i)
                    2 * X * centered_y[, i])

    # initial beta parameter for the lasso algo
    # beta0 <- lapply(1:nb_series, function(i) {
    #   coeff <-
    #     as.numeric(solve(XX + lambda * diag(ncol(XX))) %*% (t(X) * centered_y[1, i]))
    #
    #   names(coeff) <- lag_names
    #
    #   return(coeff)
    # })
    beta0 <- lapply(1:nb_series, function(i) {
      coeff <-
        as.numeric(solve(XX + lambda * diag(ncol(XX))) %*% (t(X) * centered_y[, i]))

      names(coeff) <- lag_names

      return(coeff)
    })

    # naming the columns
    names(beta0) <- paste0("y_", series_names)

    # apply lasso algo (shooting)
    # envisage to do this in parallel if there are a lot of series
    fit_lasso <-
      lapply(1:nb_series, function (i)
        lasso_shoot_cpp(
          beta = beta0[[i]],
          XX2 = XX2,
          Xy2 = Xy2[[i]],
          lambda = lambda,
          tol = 1e-05,
          max_iter = 10000
        ))
    names(fit_lasso) <- series_names
    # naming the coefficients
    for (i in 1:nb_series)
      names(fit_lasso[[i]]$beta) <- lag_names

    return(
      list(
        y = y,
        x = as.matrix(x),
        lags = lags,
        lambda = lambda,
        series_names = series_names,
        class_res = "lassoVAR",
        ym = ym,
        xm = xm,
        scales = xsd,
        #resid = NULL,
        fit_lasso = fit_lasso
      )
    )
  }
}

# fitting a pls
fit_pls_mts <- function(x,
                        lags = 1,
                        nb_hidden = 5,
                        nodes_sim = c("sobol", "halton", "unif"),
                        activ = c("relu", "sigmoid", "tanh",
                                  "leakyrelu", "elu", "linear"),
                        hidden_layer_bias = FALSE,
                        direct_link = FALSE,
                        a = 0.01,
                        B = 5,
                        seed = 1)
{
  stopifnot(is.wholenumber(nb_hidden))
  stopifnot(nb_hidden > 0)
  p <- ncol(x)

  series_names <- colnames(x)
  nb_series <- ncol(x)
  # !!! because the ts object is in reverse order
  rev_x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(rev_x, lags)
  observed_values <- y <- y_x$y
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)

  if (direct_link == FALSE)
  {
    list_xreg <- create_new_predictors(
      as.matrix(y_x$xreg),
      nb_hidden = nb_hidden,
      method = nodes_sim,
      activ = activ,
      hidden_layer_bias = hidden_layer_bias,
      a = a,
      seed = seed
    )
    xreg <- as.matrix(list_xreg$predictors[,-(1:(lags * nb_series))])
  } else {
    list_xreg <- create_new_predictors(
      as.matrix(y_x$xreg),
      nb_hidden = nb_hidden,
      method = nodes_sim,
      activ = activ,
      hidden_layer_bias = hidden_layer_bias,
      a = a,
      seed = seed
    )
    xreg <- list_xreg$predictors
  }

  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)
  x_scaled <- my_scale(xreg)
  xreg <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  ym_mat <- tcrossprod(rep(1, nrow(xreg)), ym)

  fit_obj_pls <- fit_pls(x = xreg, y = centered_y,
                         ncomp = B)
  fitted_values <-
    rev_matrix_cpp(ym_mat + predict_pls(fit_obj = fit_obj_pls,
                                        newx = xreg))
  resids <- fitted_values - rev_matrix_cpp(y)

  return(
    list(
      y = y,
      x = x,
      lags = lags,
      series_names = series_names,
      class_res = "pls",
      nb_hidden = nb_hidden,
      method = nodes_sim,
      w = list_xreg$w,
      activ = list_xreg$activ,
      hidden_layer_bias = hidden_layer_bias,
      hidden_layer_index = list_xreg$hidden_layer_index,
      activ_name = activ,
      direct_link = direct_link,
      seed = seed,
      nn_xm = list_xreg$xm,
      nn_xsd = list_xreg$xsd,
      ym = ym,
      xm = xm,
      scales = xsd,
      coefficients = fit_obj_pls$coefficients,
      Xmeans = fit_obj_pls$Xmeans,
      Ymeans = fit_obj_pls$Ymeans,
      ncomp = B,
      fitted_values = fitted_values,
      resid = resids
    )
  )
}

# fitting a pcr
fit_pcr_mts <- function(x, lags = 1,
                        ncomp = 5)
{
  stopifnot(ncomp <= lags * ncol(x))

  p <- ncol(x)

  series_names <- colnames(x)
  # !!! because the ts object is in reverse order
  rev_x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(rev_x, lags)
  observed_values <- y <- y_x$y

  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)
  x_scaled <- my_scale(y_x$xreg)
  xreg <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  ym_mat <- tcrossprod(rep(1, nrow(xreg)), ym)

  fit_obj_pcr <- fit_pcr(x = xreg, y = centered_y,
                         ncomp = ncomp)

  fitted_values <- rev_matrix_cpp(ym_mat + sapply(1:p,
                                                  function (i)
                                                    predict(fit_obj_pcr[[i]],
                                                            newdata = xreg)[, , ncomp]))

  resids <- fitted_values - rev_matrix_cpp(y)

  return(
    list(
      y = y,
      x = x,
      lags = lags,
      series_names = series_names,
      class_res = "pcr",
      ym = ym,
      xm = xm,
      scales = xsd,
      ncomp = ncomp,
      fit_obj = fit_obj_pcr,
      fitted_values = fitted_values,
      resid = resids
    )
  )
}



# 2 - ensemble fitting models -------------------------------------------
# Fitting a(n unconstrained) regression  model to multiple time series
fit_scn_mts <- function(x,
                        lags = 1,
                        activ = c("sigmoid", "tanh"),
                        hidden_layer_bias = FALSE,
                        B = 100,
                        nu = 0.1,
                        lam = 100,
                        r = 0.3,
                        tol = 1e-10,
                        col_sample = 1,
                        method = c("greedy", "direct"),
                        type_optim = c("nlminb", "nmkb"),
                        verbose = FALSE)
{
  stopifnot(col_sample >= 0.5 || col_sample <= 1)
  series_names <- colnames(x)
  # !!! because the ts object is in reverse order
  rev_x <- rev_matrix_cpp(as.matrix(x))
  y_x <- create_train_inputs_cpp(rev_x, lags)
  observed_values <- y <- y_x$y
  activ <- match.arg(activ)
  method <- match.arg(method)
  type_optim <- match.arg(type_optim)

  fit_obj <- fit_SCN(
    x = y_x$xreg,
    y =  y_x$y,
    B = B,
    nu = nu,
    col_sample = col_sample,
    lam = lam,
    r = r,
    tol = tol,
    activation = activ,
    hidden_layer_bias = hidden_layer_bias,
    type_optim = type_optim,
    method = method,
    verbose = verbose
  )

  resids <- rev_matrix_cpp(fit_obj$current_error)

  return(
    list(
      y = y,
      x = x,
      lags = lags,
      betas_opt = fit_obj$betas_opt,
      ws_opt = fit_obj$ws_opt,
      activ = activ,
      hidden_layer_bias = hidden_layer_bias,
      nu = fit_obj$nu,
      lam = lam,
      r = r,
      tol = tol,
      col_sample = col_sample,
      col_sample_indices = fit_obj$col_sample_indices,
      method = method,
      type_optim = fit_obj$type_optim,
      series_names = series_names,
      class_res = "scn",
      ym = fit_obj$ym,
      xm = fit_obj$xm,
      xsd = fit_obj$xsd,
      fitted_values = rev_matrix_cpp(y_x$y) - resids,
      resid = resids
    )
  )
}

# glmboosting
fit_glmboost_mts <- function(x,
                             B = 10,
                             eta = 0.1,
                             lags = 1,
                             nb_hidden = 1,
                             nodes_sim = c("sobol", "halton", "unif"),
                             activ = c("relu", "sigmoid", "tanh",
                                       "leakyrelu", "elu", "linear"),
                             hidden_layer_bias = FALSE,
                             direct_link = FALSE,
                             a = 0.01,
                             seed = 1)
{
  stopifnot(is.wholenumber(nb_hidden))
  nb_series <- ncol(x)
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  blocks_index <- TRUE
  col_index <- TRUE

  series_names <- colnames(x)

  rev_x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(rev_x, lags)

  if (nb_hidden > 0)
  {
    if (direct_link == FALSE)
    {
      list_xreg <- create_new_predictors(
        as.matrix(y_x$xreg),
        nb_hidden = nb_hidden,
        method = nodes_sim,
        activ = activ,
        hidden_layer_bias = hidden_layer_bias,
        a = a,
        seed = seed
      )
      xreg <-
        as.matrix(list_xreg$predictors[,-(1:(lags * nb_series))])
    } else {
      list_xreg <- create_new_predictors(
        as.matrix(y_x$xreg),
        nb_hidden = nb_hidden,
        method = nodes_sim,
        activ = activ,
        hidden_layer_bias = hidden_layer_bias,
        a = a,
        seed = seed
      )
      xreg <- list_xreg$predictors
    }
  } else {
    xreg <-  y_x$xreg[, TRUE]
  }

  # observed values, minus the lags (!)(beware)(!)
  observed_values <- y <- y_x$y

  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)
  x_scaled <- my_scale(xreg)
  xreg <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  fit_control <- mboost::boost_control(
    mstop = B,
    nu = eta,
    risk = "inbag",
    stopintern = FALSE,
    center = TRUE,
    trace = FALSE
  )

  fit_obj <- suppressWarnings(lapply(1:nb_series,
                                     function (i)
                                       mboost::glmboost(
                                         x = xreg,
                                         y = centered_y[, i],
                                         center = TRUE,
                                         control = fit_control
                                       )))
  names(fit_obj) <- series_names

  resid <- rev_matrix_cpp(sapply(1:nb_series,
                  function (i)
                    residuals(fit_obj[[i]])))

  return(
    list(
      y = y,
      x = x,
      lags = lags,
      col_index = col_index,
      fit_obj = fit_obj,
      series_names = series_names,
      nb_series = nb_series,
      class_res = "glmboost",
      nb_hidden = nb_hidden,
      B = B,
      eta = eta,
      method = nodes_sim,
      seed = seed,
      w = list_xreg$w,
      activation_name = activ,
      activ = list_xreg$activ,
      hidden_layer_bias = hidden_layer_bias,
      hidden_layer_index = list_xreg$hidden_layer_index,
      direct_link = direct_link,
      seed = seed,
      nn_xm = list_xreg$xm,
      nn_xsd = list_xreg$xsd,
      ym = ym,
      xm = xm,
      scales = xsd,
      resid = resid)
    )


}

# xgboosting
fit_xgboost_mts <- function(x,
                            B = 10,
                            eta = 0.5,
                            lambda = 0.1,
                            alpha = 0.5,
                            lags = 1,
                            nb_hidden = 1,
                            nodes_sim = c("sobol", "halton", "unif"),
                            activ = c("relu", "sigmoid", "tanh",
                                      "leakyrelu", "elu", "linear"),
                            hidden_layer_bias = FALSE,
                            direct_link = FALSE,
                            a = 0.01,
                            seed = 1)
{
  stopifnot(is.wholenumber(nb_hidden))
  stopifnot(nb_hidden > 0)
  nb_series <- ncol(x)
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  blocks_index <- TRUE
  col_index <- TRUE

  series_names <- colnames(x)

  rev_x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(rev_x, lags)

  if (nb_hidden > 0)
  {
    if (direct_link == FALSE)
    {
      list_xreg <- create_new_predictors(
        as.matrix(y_x$xreg),
        nb_hidden = nb_hidden,
        method = nodes_sim,
        activ = activ,
        hidden_layer_bias = hidden_layer_bias,
        a = a,
        seed = seed
      )
      xreg <-
        as.matrix(list_xreg$predictors[,-(1:(lags * nb_series))])
    } else {
      list_xreg <- create_new_predictors(
        as.matrix(y_x$xreg),
        nb_hidden = nb_hidden,
        method = nodes_sim,
        activ = activ,
        hidden_layer_bias = hidden_layer_bias,
        a = a,
        seed = seed
      )
      xreg <- list_xreg$predictors
    }
  } else {
    xreg <-  y_x$xreg[, TRUE]
  }

  # observed values, minus the lags (!)(beware)(!)
  observed_values <- y <- y_x$y

  ym <- colMeans(y)
  centered_y <- my_scale(x = y, xm = ym)
  x_scaled <- my_scale(xreg)
  xreg <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  fit_obj <- lapply(1:nb_series,
                    function (i) {
                      data.train <- cbind(y = centered_y[, i], xreg)

                      dtrain <-
                        xgboost::xgb.DMatrix(data.train, label = data.train[, 1])

                      xgboost::xgboost(
                        data = dtrain,
                        nrounds = B,
                        params = list(
                          eta = eta,
                          lambda = lambda,
                          alpha = alpha,
                          objective = 'reg:squarederror'
                        ),
                        verbose = 0
                      )
                    })
  names(fit_obj) <- series_names

  resid <- rev_matrix_cpp(sapply(1:nb_series,
                  function(i){
                    observed_values[, i] - predict(fit_obj[[i]], xreg)
                  }))

  return(
    list(
      y = y,
      x = x,
      lags = lags,
      col_index = col_index,
      fit_obj = fit_obj,
      series_names = series_names,
      nb_series = nb_series,
      class_res = "xgboost",
      nb_hidden = nb_hidden,
      B = B,
      eta = eta,
      lambda = lambda,
      alpha = alpha,
      method = nodes_sim,
      seed = seed,
      w = list_xreg$w,
      activation_name = activ,
      activ = list_xreg$activ,
      hidden_layer_bias = hidden_layer_bias,
      hidden_layer_index = list_xreg$hidden_layer_index,
      direct_link = direct_link,
      seed = seed,
      nn_xm = list_xreg$xm,
      nn_xsd = list_xreg$xsd,
      ym = ym,
      xm = xm,
      scales = xsd,
      resid = resid
    )
  )

}

# 3 - generic fitting models -------------------------------------------

# generic function
fit_generic_mts <- function(x,
                            fit_func,
                            predict_func,
                            lags = 1,
                            nb_hidden = 5,
                            nodes_sim = c("sobol", "halton", "unif"),
                            activ = c("relu", "sigmoid", "tanh",
                                      "leakyrelu", "elu", "linear"),
                            hidden_layer_bias = FALSE,

                            col_sample = 1,
                            a = 0.01,
                            seed = 1)
{
  stopifnot(col_sample > 0 || col_sample <= 1)
  stopifnot(is.wholenumber(nb_hidden))

  series_names <- colnames(x)
  # !!! because the ts object is in reverse order
  x <- rev_matrix_cpp(x)
  y_x <- create_train_inputs_cpp(x, lags)
  observed_values <- y <- y_x$y
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)

  list_xreg <- create_new_predictors(
    y_x$xreg,
    nb_hidden = nb_hidden,
    hidden_layer_bias = hidden_layer_bias,
    method = nodes_sim,
    activ = activ,
    a = a,
    seed = seed
  )
  xreg <- list_xreg$predictors

  ym <- colMeans(y)
  centered_y <- after::my_scale(x = y, xm = ym)
  x_scaled <- after::my_scale(xreg)
  xreg <- x_scaled$res
  xm <- x_scaled$xm
  xsd <- x_scaled$xsd

  k_p <- lags * ncol(x)
  index <- 1:k_p

  # model fitting
  fit_gen <- lapply(1:ncol(centered_y),
                    function(i)
                      fit_func(x = xreg, y = centered_y[, i]))


  genfit <- try(sapply(1:ncol(centered_y), function (i)
      predict_func(fit_gen[[i]],
                  newx = xreg)), silent = TRUE)

  if (class(genfit) == "try-error")
  {
    genfit <- try(sapply(1:ncol(centered_y), function (i)
                  predict_func(fit_gen[[i]],
                  newx = xreg)), silent = TRUE)
  }

  fitted_values <-
    after::rev_matrix_cpp(genfit + matrix(rep(ym, each = nrow(genfit)),
                                   ncol = ncol(genfit)))
  resid <- after::rev_matrix_cpp(observed_values) - fitted_values

  if (nb_hidden > 0)
  {
    return(
      list(
        y = y,
        lags = lags,
        series_names = series_names,
        class_res = "gen",
        nb_hidden = nb_hidden,
        method = nodes_sim,
        w = list_xreg$w,
        activ = list_xreg$activ,
        hidden_layer_bias = hidden_layer_bias,
        hidden_layer_index = list_xreg$hidden_layer_index,
        seed = seed,
        nn_xm = list_xreg$xm,
        nn_xsd = list_xreg$xsd,
        ym = ym,
        xm = xm,
        scales = xsd,
        resid = resid,
        fit_obj = fit_gen,
        predict_func = predict_func
      )
    )
  } else { #nb_hidden <= 0
    return(
      list(
        y = y,
        lags = lags,
        series_names = series_names,
        class_res = "gen",
        nb_hidden = nb_hidden,
        method = nodes_sim,
        seed = seed,
        nn_xm = list_xreg$xm,
        nn_xsd = list_xreg$xsd,
        ym = ym,
        xm = xm,
        scales = xsd,
        coef = lscoef,
        resid = resid,
        fit_obj = fit_gen,
        predict_func = predict_func
      )
    )
  }

}
