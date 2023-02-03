test_that("linear model's analytical gradient is close to numerical one", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = 'linear', seed = 1918)
  loglik_func <- function (coef) { 
    calc_linear_loglik(coef, data$design, data$outcome)
  }
  set.seed(615)
  n_test <- 10
  grads_are_close <- TRUE
  for (i in 1:n_test) {
    if (!grads_are_close) break
    regcoef <- rnorm(n_pred)
    analytical_grad <- calc_linear_grad(regcoef, data$design, data$outcome)
    numerical_grad <- approx_grad_via_finite_diff(loglik_func, regcoef)
    grads_are_close <- are_all_close(
      analytical_grad, numerical_grad, abs_tol = Inf, rel_tol = 1e-3
    )
  }
  expect_true(grads_are_close)
})

test_that("logit model's analytical gradient is close to numerical one", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = 'logit', seed = 1918)
  loglik_func <- function (coef) { 
    calc_logit_loglik(coef, data$design, data$outcome)
  }
  set.seed(615)
  n_test <- 10
  grads_are_close <- TRUE
  for (i in 1:n_test) {
    if (!grads_are_close) break
    regcoef <- rnorm(n_pred)
    analytical_grad <- calc_logit_grad(regcoef, data$design, data$outcome)
    numerical_grad <- approx_grad_via_finite_diff(loglik_func, regcoef)
    grads_are_close <- are_all_close(
      analytical_grad, numerical_grad, abs_tol = Inf, rel_tol = 1e-3
    )
  }
  expect_true(grads_are_close)
})