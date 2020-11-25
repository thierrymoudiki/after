
# 0 - imports ----------------------------------------------------------

library(glmnet)





# 1 - univariate time series test ----------------------------------------------------------

x <- AirPassengers

test_that("test univariate time series methods", {
  #expect_equal(round(after::arimaf(x)$mean[3]), 498)
  expect_equal(after::dynrmf(x)$mean[3], 459.68)
  expect_equal(round(after::dynrmf(x, fit_func=glmnet::glmnet,
                                   predict_func=predict)$mean[3]), 294)
  expect_equal(round(after::eatf(x)$mean[3]), 497)
  expect_equal(round(after::mmf(x)$mean[3]), 439)
  expect_equal(round(after::polythetaf(x)$mean[3]), 432)
  expect_equal(round(after::smoothf(x)$mean[3]), 433)
})

# 2 - univariate time series test ----------------------------------------------------------

x <- fpp::insurance

test_that("test multivariate time series methods", {
  # varf, xgboost remain to be tested
  expect_equal(check_closeness(as.numeric(after::ridge2f(x)$preds[3, 1]), 14.56234), TRUE)
  expect_equal(check_closeness(as.numeric(after::ridgef(x)$preds[3, 1]), 14.56234), TRUE)
  expect_equal(check_closeness(as.numeric(after::krlsf(x)$preds[3, 1]), 13.8922), TRUE)
  expect_equal(check_closeness(as.numeric(after::lmf(x)$preds[3, 1]), 14.94347), TRUE)
  expect_equal(check_closeness(as.numeric(after::mtsmeanf(x)$preds[3, 1]), 13.60435), TRUE)
  expect_equal(check_closeness(as.numeric(after::mtsmedianf(x)$preds[3, 1]), 13.28744), TRUE)
  expect_equal(check_closeness(as.numeric(after::mtsmeanf(x)$preds[3, 1]), 13.60435), TRUE)
  expect_equal(check_closeness(as.numeric(after::mtsrwf(x)$preds[3, 1]), 14.49168), TRUE)
  expect_equal(check_closeness(as.numeric(after::nnlsf(x)$preds[3, 1]), 12.62374), TRUE)
  expect_equal(check_closeness(as.numeric(after::pcrf(x)$preds[3, 1]), 13.77753), TRUE)
  expect_equal(check_closeness(as.numeric(after::plsf(x)$preds[3, 1]), 13.71616 ), TRUE)
  #expect_equal(check_closeness(as.numeric(after::glmboostf(x)$preds[3, 1]), 13.71616), TRUE)
  expect_equal(check_closeness(as.numeric(after::scnf(x)$preds[3, 1]), 15.63736), TRUE)

  })
