#
# # dataset
# response_index <- 4
# x <- data.matrix(quakes[, -response_index])
# y <- quakes[, response_index]
# (n <- nrow(x))
# (p <- ncol(x))
#
# set.seed(1225)
# (ix <- sample(n, size=floor(0.8*n)))
# x_train <- x[-ix, ]
# y_train <- y[-ix]
# x_test <- x[ix, ]
# y_test <- y[ix]
#
# (n_train <- nrow(x_train))
# lambdas <- 10^seq(from=-10, to=10,
#                   length.out = 100)
#
#
# # fit
# fit1 <- glmnet::glmnet(x_train, y_train,
#                        alpha=0.5,
#                        lambda=lambdas)
# residuals_train <- predict(fit1, newx=x_train) - y_train
# (bic <- n_train*log(colMeans(residuals_train^2)) + (fit1$df + 2)*log(n_train))
# plot(bic, type='l')
#
#
# # preds
# (s_opt <- lambdas[which.min(bic)])
# (preds <- as.numeric(glmnet::predict.glmnet(fit1,
#                                             newx = x_test,
#                                             s = s_opt)))
# sqrt(mean((preds - y_test)^2))
# hist(quakes$mag)
# summary(quakes$mag)
