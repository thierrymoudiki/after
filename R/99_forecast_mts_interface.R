# 0 - Import functions ---------------------------------------------------

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir("regmtsfuncs")

# 1 - Individual models ---------------------------------------------------

# Forecasts from ridge2
ridge2f <- function(x,
                    lags = 1,
                    nb_hidden = 5,
                    nodes_sim = c("sobol", "halton", "unif"),
                    activ = c("relu", "sigmoid", "tanh",
                              "leakyrelu", "elu", "linear"),
                    hidden_layer_bias = FALSE,
                    col_sample = 1,
                    row_sample = 1,
                    a = 0.01,
                    lambda_1 = 0.1,
                    lambda_2 = 0.1,
                    seed = 1,
                    h = 5,
                    type_ci = "none",
                    type_forecast = c("recursive", "direct"),
                    level = 95)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  stopifnot(col_sample >= 0.5 && col_sample <= 1)
  stopifnot(row_sample >= 0.5 && row_sample <= 1)
  type_forecast <- match.arg(type_forecast)

  # Fitting a regularized regression  model to multiple time series
  fit_obj <- fit_ridge_mts(
    x,
    lags = lags,
    nb_hidden = nb_hidden,
    fit_method = "ridge2",
    nodes_sim = nodes_sim,
    activ = activ,
    hidden_layer_bias = hidden_layer_bias,
    col_sample = col_sample,
    row_sample = row_sample,
    a = a,
    lambda_1 = lambda_1,
    lambda_2 = lambda_2,
    seed = seed
  )

  preds <- ts(data = fcast_obj_mts(
                    fit_obj,
                    h = h,
                    type_ci = type_ci,
                    type_forecast = type_forecast,
                    level = level),
              start = start_preds, frequency = freq_x)

  resids <- ts(data = fit_obj$resid,
               start = start_preds, frequency = freq_x)

  # Forecast from fit_obj
  return(list(preds = preds,
              resid = resids))
}

# Forecasts from ridge
ridgef <- function(x,
                   lags = 1,
                   nb_hidden = 5,
                   fit_method = c("ridge", "mgaussian"),
                   nodes_sim = c("sobol", "halton", "unif"),
                   activ = c("relu", "sigmoid", "tanh",
                             "leakyrelu", "elu", "linear"),
                   hidden_layer_bias = FALSE,
                   a = 0.01,
                   lambda = 0.1,
                   alpha = 0.5,
                   seed = 1,
                   h = 5,
                   type_ci = "none",
                   type_forecast = c("recursive", "direct"),
                   level = 95)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  fit_method <- match.arg(fit_method)
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  type_forecast <- match.arg(type_forecast)

  # Fitting a regularized regression  model to multiple time series
  fit_obj <- fit_ridge_mts(
    x,
    lags = lags,
    nb_hidden = nb_hidden,
    fit_method = fit_method,
    nodes_sim = nodes_sim,
    activ = activ,
    hidden_layer_bias = hidden_layer_bias,
    a = a,
    lambda = lambda,
    alpha = alpha,
    seed = seed
  )

  # Forecast from fit_obj

  preds <- ts(data = fcast_obj_mts(
    fit_obj,
    h = h,
    type_ci = type_ci,
    type_forecast = type_forecast,
    level = level),
    start = start_preds, frequency = freq_x)

  resids <- ts(data = fit_obj$resid,
               start = start_preds, frequency = freq_x)

  # Forecast from fit_obj
  return(list(preds = preds,
              resid = resids))
}

# Forecasts from unsconstrained linear model
lmf <- function(x,
                lags = 1,
                nb_hidden = 5,
                nodes_sim = c("sobol", "halton", "unif"),
                activ = c("relu", "sigmoid", "tanh",
                          "leakyrelu", "elu", "linear"),
                hidden_layer_bias = FALSE,
                a = 0.01,
                seed = 1,
                h = 5,
                type_ci = "none",
                type_forecast = c("recursive", "direct"),
                level = 95)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  type_forecast <- match.arg(type_forecast)

  # Fitting a(n unconstrained) regression  model to multiple time series
  fit_obj <- fit_lm_mts(
    x,
    lags = lags,
    nb_hidden = nb_hidden,
    nodes_sim = nodes_sim,
    activ = activ,
    hidden_layer_bias = hidden_layer_bias,
    a = a,
    seed = seed
  )

  # Forecast from fit_obj
  preds <- ts(data = fcast_obj_mts(
    fit_obj,
    h = h,
    type_ci = type_ci,
    type_forecast = type_forecast,
    level = level),
    start = start_preds, frequency = freq_x)

  resids <- ts(data = fit_obj$resid,
               start = start_preds, frequency = freq_x)

  # Forecast from fit_obj
  return(list(preds = preds,
              resid = resids))
}

# Forecasts from constrained linear model (positive coefficients)
nnlsf <- function(x,
                  lags = 1,
                  nb_hidden = 5,
                  nodes_sim = c("sobol", "halton", "unif"),
                  activ = c("relu", "sigmoid", "tanh",
                            "leakyrelu", "elu", "linear"),
                  hidden_layer_bias = FALSE,
                  a = 0.01,
                  seed = 1,
                  h = 5,
                  type_ci = "none",
                  type_forecast = c("recursive", "direct"),
                  level = 95)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  type_forecast <- match.arg(type_forecast)

  # Fitting a(n unconstrained) regression  model to multiple time series
  fit_obj <- fit_nnls_mts(
    x,
    lags = lags,
    nb_hidden = nb_hidden,
    nodes_sim = nodes_sim,
    activ = activ,
    hidden_layer_bias = hidden_layer_bias,
    a = a,
    seed = seed
  )

  # Forecast from fit_obj
  preds <- ts(data = fcast_obj_mts(
    fit_obj,
    h = h,
    type_ci = type_ci,
    type_forecast = type_forecast,
    level = level),
    start = start_preds, frequency = freq_x)

  resids <- ts(data = fit_obj$resid,
               start = start_preds, frequency = freq_x)

  # Forecast from fit_obj
  return(list(preds = preds,
              resid = resids))
}

# RENAME # RENAME # RENAME
# Forecasts from krls models
krlsf <- function(x,
                lags = 1,
                lambda_krls = 0.1,
                l = 0.5,
                sigma = 10,
                d = 1,
                kernel_type = c("gaussian", "polynomial",
                                "matern32", "matern52"),
                inv_method = c("chol", "ginv"),
                h = 5,
                type_ci = "none",
                type_forecast = c("recursive", "direct"),
                level = 95)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  kernel_type <- match.arg(kernel_type)
  inv_method  <- match.arg(inv_method)
  type_forecast <- match.arg(type_forecast)

  # Fitting a gaussian process regression  model to multiple time series
  fit_obj <- fit_krls_mts(
    x,
    lags = lags,
    lambda_krls = lambda_krls,
    l = l,
    sigma = sigma,
    d = d,
    kernel_type = kernel_type,
    inv_method = inv_method
  )

  # Forecast from fit_obj
  preds <- ts(data = fcast_obj_mts(
    fit_obj,
    h = h,
    type_ci = type_ci,
    type_forecast = type_forecast,
    level = level),
    start = start_preds, frequency = freq_x)

  resids <- ts(data = fit_obj$resid,
               start = start_preds, frequency = freq_x)

  # Forecast from fit_obj
  return(list(preds = preds,
              resid = resids))
}

# Forecasts from unsconstrained VAR models
varf <- function(x,
                 lags = 1,
                 penalization = c("none", "l1"),
                 lambda = 0.1, # only for penalization == "l1"
                 type_VAR = c("const", "trend",
                              "both", "none"), # currently: only for penalization == "none"
                 type_forecast = c("recursive", "direct"),
                 h = 5,
                 type_ci = "none",
                 level = 95)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  type_VAR <- match.arg(type_VAR)
  penalization <- match.arg(penalization)
  type_forecast <- match.arg(type_forecast)

  if (penalization == "none")
  {
    # Fitting a VAR model to multiple time series
    fit_obj <- fit_var_mts(x,
                           lags = lags,
                           penalization = penalization,
                           type_VAR = type_VAR)

    # Forecast from fit_obj

    preds <- ts(data = fcast_obj_mts(
      fit_obj,
      h = h,
      type_ci = type_ci,
      type_forecast = type_forecast,
      level = level),
      start = start_preds, frequency = freq_x)

    resids <- ts(data = fit_obj$resid,
                 start = start_preds, frequency = freq_x)

    # Forecast from fit_obj
    return(list(preds = preds,
                resid = resids))
  }

  if (penalization == "l1")
  {
    # Fitting a VAR model to multiple time series
    fit_obj <- fit_var_mts(x,
                           lags = lags,
                           penalization = penalization,
                           lambda = lambda)

    # Forecast from fit_obj

    preds <- ts(data = fcast_obj_mts(
      fit_obj,
      h = h,
      # type_ci = type_ci,
      type_forecast = type_forecast),#,
      #level = level),
      start = start_preds, frequency = freq_x)

    resids <- ts(data = fit_obj$resid,
                 start = start_preds, frequency = freq_x)

    # Forecast from fit_obj
    return(list(preds = preds,
                resid = resids))
  }

}

# principal components regression
pcrf <- function(x, ncomp = ncol(x), lags = 1,
                 h = 5,
                 type_ci = "none",
                 type_forecast = c("recursive", "direct"),
                 level = 95)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  type_forecast <- match.arg(type_forecast)

  fit_obj <- fit_pcr_mts(x = x, lags = lags,
                                         ncomp = ncomp)

  preds <- ts(data = fcast_obj_mts(
    fit_obj,
    h = h,
    type_ci = type_ci,
    type_forecast = type_forecast,
    level = level),
    start = start_preds, frequency = freq_x)

  resids <- ts(data = fit_obj$resid,
               start = start_preds, frequency = freq_x)

  # Forecast from fit_obj
  return(list(preds = preds,
              resid = resids))
}


# 2 - Ensemble models ---------------------------------------------------

# Bootstrap aggregating
bag_ridge2f <- function(x,
                        B = 10,
                        agg = c("mean", "median"),
                        lags = 1,
                        nb_hidden = 5,
                        nodes_sim = c("sobol", "halton", "unif"),
                        activ = c("relu", "sigmoid", "tanh",
                                  "leakyrelu", "elu", "linear"),
                        hidden_layer_bias = FALSE,
                        verbose = TRUE,
                        col_sample = 0.5,
                        row_sample = 1,
                        a = 0.01,
                        lambda_1 = 0.1,
                        lambda_2 = 0.1,
                        h = 5)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  nb_series <- ncol(x)
  series_names <- colnames(x)
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  stopifnot(col_sample >= 0.5 && col_sample <= 1)
  stopifnot(row_sample >= 0.5 && row_sample <= 1)
  agg <- match.arg(agg)

  if ((nodes_sim == "sobol" ||
       nodes_sim == "halton") && (col_sample == 1 &&
                                  row_sample == 1))
  {
    stop(
      "For nodes_sim == 'sobol' || nodes_sim == 'halton', we must have \n
      col_sample < 1 && row_sample < 1"
    )
  }

  # if (verbose == TRUE)
  # {
  #   pb <- txtProgressBar(min = 0, max = B, style = 3)
  # }

  # CONSIDER USING A LOOP or lapply
  i <- NULL
  `%op%` <-  foreach::`%do%`
  list_preds <- foreach::foreach(i = 1:B,
                                 .errorhandling = "remove") %op% {

                                   # if (verbose == TRUE)
                                   # {
                                   #    setTxtProgressBar(pb, i)
                                   # }

                                   fit_obj_ridge2f <- ridge2f(
                                     x = x,
                                     lags = lags,
                                     nb_hidden = nb_hidden,
                                     nodes_sim = nodes_sim,
                                     activ = activ,
                                     hidden_layer_bias = hidden_layer_bias,
                                     col_sample = col_sample,
                                     row_sample = row_sample,
                                     a = a,
                                     lambda_1 = lambda_1,
                                     lambda_2 = lambda_2,
                                     seed = i,
                                     h = h,
                                     type_ci = "none"
                                   )$preds

                                 }
  # if (verbose == TRUE)
  # {
  #   close(pb)
  # }

  base_preds <- lapply(1:nb_series, function(j)
    sapply(1:length(list_preds), function(i)
      list_preds[[i]][, j]))
  names(base_preds) <- series_names

  preds <- switch(agg,
                  "mean" = sapply(1:nb_series,
                                  function(i)
                                    rowMeans(base_preds[[i]])),
                  "median" = sapply(1:nb_series,
                                    function(i)
                                      apply(base_preds[[i]], 1, median)))
  #colnames(mean_preds) <- series_names
  #colnames(median_preds) <- series_names

  return(list(
    preds = preds,
    all_preds = base_preds
  ))
  }

# Stacked generalization
stack_data_ridge2f <- function(x,
                               B = 10,
                               lags = 1,
                               nb_hidden = 5,
                               nodes_sim = c("sobol", "halton", "unif"),
                               activ = c("relu", "sigmoid", "tanh",
                                         "leakyrelu", "elu", "linear"),
                               hidden_layer_bias = FALSE,
                               col_sample = 0.9,
                               row_sample = 0.8,
                               a = 0.01,
                               lambda_1 = 0.1,
                               lambda_2 = 0.1)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  nb_series <- ncol(x)
  series_names <- colnames(x)
  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  nb_points <- nrow(x)
  index_train <- 1:floor(nb_points / 2)
  x_train <- x[index_train, ]
  x_test <- x[-index_train, ]
  h <- nrow(x_test)

  if ((nodes_sim == "sobol" ||
       nodes_sim == "halton") && (col_sample == 1 &&
                                  row_sample == 1))
  {
    stop(
      "For nodes_sim == 'sobol' || nodes_sim == 'halton', we must have \n
      col_sample < 1 && row_sample < 1"
    )
  }

  i <- NULL
  `%op%` <-  foreach::`%do%`
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  list_preds <- foreach::foreach(i = 1:B,
                                 .errorhandling = "remove") %op% {
                                   utils::setTxtProgressBar(pb, i)
                                   ridge2f(
                                     x = x_train,
                                     lags = lags,
                                     nb_hidden = nb_hidden,
                                     nodes_sim = nodes_sim,
                                     activ = activ,
                                     hidden_layer_bias = hidden_layer_bias,
                                     col_sample = col_sample,
                                     row_sample = row_sample,
                                     a = a,
                                     lambda_1 = lambda_1,
                                     lambda_2 = lambda_2,
                                     seed = i,
                                     h = h,
                                     type_ci = "none"
                                   )$preds
                                 }
  close(pb)

  base_preds <- lapply(1:nb_series, function(j)
    sapply(1:length(list_preds), function(i)
      list_preds[[i]][, j]))
  names(base_preds) <- series_names

  new_preds <- foreach::foreach (i = 1:nb_series, .combine = cbind)%op%
  {
    colnames(base_preds[[i]]) <- paste0(series_names[i], 1:B)
    base_preds[[i]]
  }

  newdata <- cbind(x_test, new_preds)
  index <- duplicated(t(newdata))

  return(list(
    base_preds = base_preds,
    new_preds = new_preds,
    newdata = newdata[, !index]
  ))
  }

# glmboosting
glmboostf <- function(x, B = 10, eta = 0.1, lags = 1, nb_hidden = 1,
                      nodes_sim = c("sobol", "halton", "unif"),
                      activ = c("relu", "sigmoid", "tanh", "leakyrelu", "elu", "linear"),
                      hidden_layer_bias = FALSE,
                      a = 0.01, direct_link = FALSE, seed = 1,
                      h = 5,
                      type_ci = "none",
                      type_forecast = c("recursive", "direct"),
                      level = 95)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x
  type_forecast <- match.arg(type_forecast)

  # fitting the model
  fit_obj <- fit_glmboost_mts(x = x, B = B, eta = eta,
                                              lags = lags, nb_hidden = nb_hidden,
                                              nodes_sim = match.arg(nodes_sim),
                                              activ = match.arg(activ),
                                              hidden_layer_bias = hidden_layer_bias,
                                              direct_link = direct_link,
                                              a = a, seed = seed)

  type_forecast <- match.arg(type_forecast)

  preds <- ts(data = fcast_obj_mts(
    fit_obj,
    h = h,
    type_ci = type_ci,
    type_forecast = type_forecast,
    level = level),
    start = start_preds, frequency = freq_x)

  resids <- ts(data = fit_obj$resid,
               start = start_preds, frequency = freq_x)

  # Forecast from fit_obj
  return(list(preds = preds,
              resid = resids))
}

# xgboosting
xgboostf <- function(x, B = 10, eta = 0.1,
                     lambda = 0.1, alpha = 0.5,
                     lags = 1, nb_hidden = 5,
                      nodes_sim = c("sobol", "halton", "unif"),
                      activ = c("relu", "sigmoid", "tanh", "leakyrelu", "elu", "linear"),
                      hidden_layer_bias = FALSE,
                      a = 0.01, direct_link = FALSE, seed = 1,
                      h = 5,
                      type_ci = "none",
                      type_forecast = c("recursive", "direct"),
                      level = 95)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  type_forecast <- match.arg(type_forecast)

  # fitting the model
  fit_obj <- fit_xgboost_mts(x = x, B = B, eta = eta,
                                             lambda = lambda, alpha = alpha,
                                              lags = lags, nb_hidden = nb_hidden,
                                              nodes_sim = match.arg(nodes_sim),
                                              activ = match.arg(activ),
                                              hidden_layer_bias = hidden_layer_bias,
                                              direct_link = direct_link,
                                              a = a, seed = seed)

  preds <- ts(data = fcast_obj_mts(
    fit_obj,
    h = h,
    type_ci = type_ci,
    type_forecast = type_forecast,
    level = level),
    start = start_preds, frequency = freq_x)

  resids <- ts(data = fit_obj$resid,
               start = start_preds, frequency = freq_x)

  # Forecast from fit_obj
  return(list(preds = preds,
              resid = resids))
}

# scn
scnf <- function(x, B = 100, lags = 1, activ = c("tanh", "sigmoid"),
                 hidden_layer_bias = FALSE,
                 eta = 0.1, lam = 100, r = 0.3, tol = 1e-10, col_sample = 1,
                 method = c("greedy", "direct"), type_optim = c("nlminb", "nmkb"),
                 verbose = FALSE,
                 h = 5,
                 type_ci = "none",
                 type_forecast = c("recursive", "direct"),
                 level = 95
                 )
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  activ <- match.arg(activ)
  method <- match.arg(method)
  type_forecast <- match.arg(type_forecast)
  type_optim <- match.arg(type_optim)

    fit_obj <- fit_scn_mts(x = x, lags = lags, activ = activ,
                                           hidden_layer_bias = hidden_layer_bias,
                                           B = B, nu = eta, lam = lam, r = r, tol = tol,
                                           col_sample = col_sample,
                                           method = method, type_optim = type_optim,
                                           verbose = verbose)

  # Forecast from fit_obj
  preds <- ts(data = fcast_obj_mts(
    fit_obj,
    h = h,
    type_ci = type_ci,
    type_forecast = type_forecast,
    level = level),
    start = start_preds, frequency = freq_x)

  resids <- ts(data = fit_obj$resid,
               start = start_preds, frequency = freq_x)

  # Forecast from fit_obj
  return(list(preds = preds,
              resid = resids))
}

# partial least squares
plsf <- function(x, B = 5, lags = 1, nb_hidden = 5,
                 nodes_sim = c("sobol","halton", "unif"),
                 activ = c("relu", "sigmoid", "tanh",
                           "leakyrelu", "elu", "linear"),
                 hidden_layer_bias = FALSE,
                 direct_link = FALSE, a = 0.01, seed = 1,
                 h = 5,
                 type_ci = "none",
                 type_forecast = c("recursive", "direct"),
                 level = 95)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  nodes_sim <- match.arg(nodes_sim)
  activ <- match.arg(activ)
  type_forecast <- match.arg(type_forecast)

  fit_obj <- fit_pls_mts(x = x, lags = lags, nb_hidden = nb_hidden,
                                         nodes_sim = nodes_sim, activ = activ,
                                         hidden_layer_bias = hidden_layer_bias,
                                         direct_link = direct_link,
                                         a = a, B = B, seed = seed)

  preds <- ts(data = fcast_obj_mts(
    fit_obj,
    h = h,
    type_ci = type_ci,
    type_forecast = type_forecast,
    level = level),
    start = start_preds, frequency = freq_x)

  resids <- ts(data = fit_obj$resid,
               start = start_preds, frequency = freq_x)

  # Forecast from fit_obj
  return(list(preds = preds,
              resid = resids))
}

# 3 - Simple benchmarks (mean, median, random walk) ---------------------------------------------------

#' Title
#'
#' @param x
#' @param h
#'
#' @return
#' @export
#'
#' @examples
mtsmeanf <- function(x, h = 5)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  series_names <- colnames(x)
  xm <- colMeans(x)
  preds <- tcrossprod(rep(1, h), xm)
  colnames(preds) <- series_names

  return(list(preds = ts(data=preds,
                         start = start_preds,
                         frequency = freq_x)))
}

#' Title
#'
#' @param x
#' @param h
#'
#' @return
#' @export
#'
#' @examples
mtsmedianf <- function(x, h = 5)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  series_names <- colnames(x)
  xm <- apply(x, 2, median)
  preds <- tcrossprod(rep(1, h), xm)
  colnames(preds) <- series_names

  return(list(preds = ts(data=preds,
                         start = start_preds,
                         frequency = freq_x)))
}

#' Title
#'
#' @param x
#' @param h
#'
#' @return
#' @export
#'
#' @examples
mtsrwf <- function(x, h = 5)
{
  if (!is.ts(x))
  {
    x <- ts(x)
  }
  freq_x <- frequency(x)
  start_fits <- start(x)
  start_preds <- tsp(x)[2] + 1/freq_x

  series_names <- colnames(x)
  last_obs <- x[nrow(x), ]
  preds <- tcrossprod(rep(1, h), last_obs)
  colnames(preds) <- series_names

  return(list(preds = ts(data=preds,
                         start = start_preds,
                         frequency = freq_x)))
}
