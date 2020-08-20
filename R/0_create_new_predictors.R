# create new predictors
create_new_predictors <- function(x,
                                  nb_hidden = 5,
                                  hidden_layer_bias = FALSE,
                                  method = c("sobol", "halton", "unif"),
                                  activ = c("relu", "sigmoid", "tanh",
                                            "leakyrelu", "elu", "linear"),
                                  a = 0.01,
                                  seed = 123)
{
  n <- nrow(x)

  if (nb_hidden > 0)
  {
    p <- ncol(x)
    method <- match.arg(method)

    # Activation function
    g <- switch(
      match.arg(activ),
      "relu" = function(x)
        x * (x > 0),
      "sigmoid" = function(x)
        1 / (1 + exp(-x)),
      "tanh" = function(x)
        tanh(x),
      "leakyrelu" = function(x)
        x * (x > 0) + a * x * (x <= 0),
      "elu" = function(x)
        x * (x >= 0) + a * (exp(x) - 1) * (x < 0),
      "linear" = function(x)
        x
    )

      if (hidden_layer_bias == FALSE)
      {
          # used for columns sample and for 'method == unif'
          set.seed(seed + 1)
          w <- remove_zero_cols(switch(
            method,
            "sobol" = 2 * t(randtoolbox::sobol(nb_hidden + 1, p)) - 1,
            "halton" = 2 * t(randtoolbox::halton(nb_hidden, p)) - 1,
            "unif" = matrix(
              runif(nb_hidden * p, min = -1, max = 1),
              nrow = p,
              ncol = nb_hidden
            )
          ))
          scaled_x <- my_scale(x)
          hidden_layer_obj <- remove_zero_cols(g(scaled_x$res %*% w),
                                                               with_index = TRUE)
          hidden_layer <- hidden_layer_obj$mat

      } else { # hidden_layer_bias == TRUE
          pp <- p + 1
          # used for columns sample and for 'method == unif'
          set.seed(seed + 1)
          w <- remove_zero_cols(switch(
            method,
            "sobol" = 2 * t(sobol(nb_hidden + 1, pp)) - 1,
            "halton" = 2 * t(halton(nb_hidden, pp)) - 1,
            "unif" = matrix(
              runif(nb_hidden * pp, min = -1, max = 1),
              nrow = pp,
              ncol = nb_hidden
            )
          ))

          scaled_x <- my_scale(x)
          hidden_layer_obj <- remove_zero_cols(g(cbind(1, scaled_x$res) %*% w),
                                                               with_index = TRUE)
          hidden_layer <- hidden_layer_obj$mat
      }

    res <- cbind(x, hidden_layer)
    nb_nodes <- ncol(hidden_layer)
    if (!is.null(nb_nodes))
      colnames(res) <- c(paste0("x", 1:p), # maybe use the real names
                         paste0("h", 1:nb_nodes))


    # if nb_hidden > 0 && (nb_predictors >= 2 && col_sample < 1)
        return(
          list(
            activ = g,
            xm = scaled_x$xm,
            xsd = scaled_x$xsd,
            w = w,
            predictors = res,
            hidden_layer_index = hidden_layer_obj$index
          )
        )
  } else {# if nb_hidden <= 0
    scaled_x <- my_scale(x)
      return(
        list(
          xm = scaled_x$xm,
          xsd = scaled_x$xsd,
          predictors = x,
          hidden_layer_index = hidden_layer_obj$index
        )
      )
  }
}
