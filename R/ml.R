#' Title
#'
#' @param x
#' @param y
#' @param lags
#' @param h
#' @param encoding
#' @param fit_func
#' @param predict_func
#' @param level
#' @param ci
#'
#' @return
#' @export
#'
#' @examples
#'
#'x_train <- data.frame(matrix(rnorm(20), ncol=2))
#'x_train$v <- as.factor(sample(c("a", "b", "c"), size=10, replace=TRUE))
#'x_test <- data.frame(matrix(rnorm(20), ncol=2))
#'x_test$v <- as.factor(sample(c("a", "b", "c"), size=10, replace=TRUE))
#'y_train <- rnorm(10)
#'#mlf(x_train, y_train, x_test)
#'
#'
#'(obj1 <- after::mlf(x_train = x_train, y_train = y_train,
#'x_test = x_test, lags = 2, fit_func = randomForest::randomForest,
#'predict_func = predict, ntree=2, mtry=2, h=5))
#'
#'(obj2 <- after::mlf(x_train = x_train, y_train = y_train, x_test = x_test,
#'fit_func = glm.fit, predict_func = predict.glm, lags = 1, h=5))
#'
#'(obj3 <- after::mlf(x_train = x_train, y_train = y_train, x_test = x_test,
#'fit_func = glmnet::glmnet, predict_func = glmnet::predict.glmnet, lags = 2,
#'h=5, alpha=0.5, lambda=0.1))
#'
mlf <- function(x_train, y_train, x_test,
                lags=1, h=5,
                fit_func=after::fit_lm,
                predict_func=after::predict_lm,
                level = c(80, 95),
                ci = c("E", "A", "T",
                       "garch", "gaussian"),
                seed=123, encoding=c("target", "model_matrix"),
                ...)

{
  stopifnot(length(y_train) == nrow(x_train))
  stopifnot(is.data.frame(x_train))

  n_x_train <- nrow(x_train)
  p_x_train <- ncol(x_train)
  n_x_test <- nrow(x_test)
  p_x_test <- ncol(x_test)
  encoding <- match.arg(encoding)


  # embed ----------------------------------------------------


  # embed x_train ----------------------------------------------------

  xy_train <- after::embed_reg(x=x_train, y=y_train,
                               testing=FALSE, lags=lags,
                               encoding="target")
  y_train_ <- xy_train$y
  x_train_ <- as.matrix(xy_train$x)
  names(y_train_) <- rownames(x_train_) <- paste0("train", 1:nrow(x_train_))
  codes <- xy_train$codes
  training_set_factors <- names(codes)


  # embed x_test ----------------------------------------------------

  # encode
  set.seed(seed)
  xy_test <- after::embed_reg(x=x_test,
                              #y=rnorm(nrow(x_test))+median(y_train),
                              y=rep(median(y_train), nrow(x_test)),
                              testing=FALSE, lags=lags,
                              encoding="target")
  x_test_ <- x_test <- xy_test$x
  for (f in training_set_factors) # loop on training set factor columns
  {
    f_ <- paste0(f, "_")
    # x_test_[[f_]] <- rep(-1, length(x_test[[f]]))
    idx_factors_test <- which(f == colnames(x_test))
    for (level in levels(x_test[[f]])) # loop on levels
    {
        x_test_[[f_]][x_test[[f]] == level] <- codes[[f]][[level]]
    }
  }
  x_test_[[f]] <- x_test_[[f_]]
  x_test_[[f_]] <- NULL
  rownames(x_test_) <- paste0("test", 1:nrow(x_test_))
  x_test_ <- rbind(tail(x_train_, 1), x_test_)
  x_train_ <- x_train_[-nrow(x_train_),]
  y_train_ <- y_train_[-length(y_train_)]
  rm(x_test); gc()


  # fit ----------------------------------------------------

  # get y_new, residual
  `%op%` <- foreach::`%do%`
  ans <- foreach::foreach(i = 1:h, .combine = c, .errorhandling = "pass")%op%
  {
    cat("\n")
    cat("i: ", i, "---------- \n")

    cat("\n")
    cat("x_train_", "\n")
    print(x_train_)
    cat("\n")

    cat("\n")
    cat("y_train_", "\n")
    print(y_train_)
    cat("\n")

    cat("x_test_", "\n")
    print(x_test_)
    cat("\n")

    set.seed(seed) # in case the algo is randomized
    fit_obj <-
      do.call(what = fit_func,
              args = c(list(x = as.matrix(x_train_),
                            y = y_train_,
                            ...)))

     # cat("fit_obj", "\n")
     # print(fit_obj)
     # cat("\n")


    # predict ----------------------------------------------------

    y_new <- try(as.numeric(predict_func(fit_obj,
                              newdata = as.matrix(rbind(x_test_[1, ],
                                              rnorm(ncol(x_test_)))))[1]), # do not change (ever)
                 silent = TRUE)
    if (class(y_new) == "try-error")
    {
      y_new <- try(as.numeric(predict_func(fit_obj,
                                newx = as.matrix(rbind(x_test_[1, ],
                                             rnorm(ncol(x_test_)))))[1]), # do not change (ever)
                   silent = TRUE)
      if (class(y_new) == "try-error")
      {
        y_new <- NA
      }
    }

    cat("y_new", "\n")
    print(y_new)
    cat("\n")


  # new xregs ----------------------------------------------------

    y_train_ <- c(y_train_, y_new)[-1]
    x_train_ <- rbind(x_train_, x_test_[1, ])[-1,]
    x_test_ <- x_test_[-1, ]
    x_test_[1, 1:lags] <- rev(tail(y_train_, lags))

    y_new
  }


  # residuals ----------------------------------------------------

  # NOPE

  # residual <- try(residuals(fit_obj), silent = TRUE)
  # # print("here0")
  # if (class(residual) == "try-error" || is.null(residual))
  # {
  #   print("here1")
  #   residual <- try(predict_func(fit_obj,
  #                                newdata = x_train) - y_train,
  #                   silent = TRUE)
  #
  #   if (class(residual) == "try-error")
  #   {
  #     print("here2")
  #     residual <- try(predict_func(fit_obj,
  #                                  newx = x_train) - y_train,
  #                     silent = TRUE)
  #
  #     if (class(training_residuals) == "try-error")
  #     {
  #       print("here3")
  #       residual <- NA
  #     }
  #   }
  # }
  #
  # return(list(f=ans, resid=residual[1:n_x_train]))

  return(list(f=ans))

  # confidence intervals ----------------------------------------------------

  # res_obj$ci <- ci
  # training_residuals <- try(residuals(fit_obj), silent = TRUE)
  # # print("here0")
  # if (class(training_residuals) == "try-error" || is.null(training_residuals))
  # {
  #   print("here1")
  #   training_residuals <- try(predict_func(fit_obj,
  #                                          newdata = x_train)
  #                             silent = TRUE)
  #   training_residuals <- training_residuals - y_train
  #   if (class(training_residuals) == "try-error")
  #   {
  #     print("here2")
  #     training_residuals <- try(predict_func(fit_obj,
  #                                            newx = x_train)
  #                               silent = TRUE)
  #     training_residuals <- training_residuals - y_train
  #     if (class(training_residuals) == "try-error")
  #     {
  #       print("here3")
  #       training_residuals <- rep(NA, length(train_index))
  #     }
  #   }
  # }
  # training_residuals <- ts(training_residuals, start = start_y,
  #                          frequency = freq_y)
  #
  # cat("\n")
  # cat("training_residuals", "\n")
  # print(training_residuals)
  # cat("\n")
  #
  # if (ci == "E")
  # {
  #   resid_fcast <- forecast::forecast(forecast::ets(training_residuals),
  #                        h = h)
  #   res_obj$mean <- ans_mean + resid_fcast$mean
  #   res_obj$lower <- ans_mean + resid_fcast$lower
  #   res_obj$upper <- ans_mean + resid_fcast$upper
  # }
  #
  # if (ci == "A")
  # {
  #   resid_fcast <- forecast::forecast(forecast::auto.arima(training_residuals),
  #                                   h = h)
  #   res_obj$mean <- ans_mean + resid_fcast$mean
  #   res_obj$lower <- ans_mean + resid_fcast$lower
  #   res_obj$upper <- ans_mean + resid_fcast$upper
  # }
  #
  # if (ci == "T")
  # {
  #   resid_fcast <- forecast::thetaf(training_residuals,
  #                                   h = h)
  #   res_obj$mean <- ans_mean + resid_fcast$mean
  #   res_obj$lower <- ans_mean + resid_fcast$lower
  #   res_obj$upper <- ans_mean + resid_fcast$upper
  # }
  #
  # if (ci == "garch")
  # {
  #   resid_fcast <- after::garch11f(training_residuals, h = h)
  #   res_obj$mean <- ans_mean + resid_fcast$mean
  #   res_obj$lower <- ans_mean + resid_fcast$lower
  #   res_obj$upper <- ans_mean + resid_fcast$upper
  #
  # }
  #
  # if (ci == "gaussian")
  # {
  #   rep_1_h <- rep(1, h)
  #   conf_ints <- lapply(level,
  #                       function (x)
  #                         tcrossprod(rep_1_h, t.test(training_residuals,
  #                                                    conf.level = x /
  #                                                      100)$conf))
  #   res_obj$mean <- ts(ans_mean + mean(training_residuals),
  #                      start = start_preds)
  #
  #   residuals_lower <- ts(sapply(1:nlevels_,
  #                                function(idx)
  #                                  conf_ints[[idx]][, 1]),
  #                         start = start_preds)
  #   residuals_upper <- ts(sapply(1:nlevels_,
  #                                function(idx)
  #                                  conf_ints[[idx]][, 2]),
  #                         start = start_preds)
  #
  #   res_obj$lower <- ts(ans_mean + residuals_lower,
  #                       start = start_preds)
  #   res_obj$upper <- ts(ans_mean + residuals_upper,
  #                       start = start_preds)
  # }
  #
  # res_obj$method <- try(class(fit_obj))
  #
  # res_obj$level <- level
  #
  # class(res_obj) <- "forecast"
  #
  # return (res_obj) # ets obj, change $mean, $lower, $upper, $residuals, $method
}
