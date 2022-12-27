test_that("linalg and optim least-sq coincide", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = 'linear', seed = 1918)
  design <- data$design; outcome <- data$outcome
  via_linalg_out <- hiper_glm(design, outcome, model = 'linear')
  via_bfgs_out <- hiper_glm(
    design, outcome, model = 'linear', option = list(mle_solver = 'BFGS')
  )
  expect_true(are_all_close(
    coef(via_linalg_out), coef(via_bfgs_out), abs_tol = 1e-2, rel_tol = 1e-2
  ))
})

test_that("newton and bfgs outputs coincide on logit model", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = 'logit', seed = 1918)
  design <- data$design; outcome <- data$outcome
  via_newton_out <- hiper_glm(design, outcome, model = 'logit')
  via_bfgs_out <- hiper_glm(
    design, outcome, model = 'logit', option = list(mle_solver = 'BFGS')
  )
  expect_true(are_all_close(
    coef(via_newton_out), coef(via_bfgs_out), abs_tol = 1e-2, rel_tol = 1e-2
  ))
})

test_that("vanilla/weighted least-sq Newton updates coincide", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model = 'logit', seed = 1918)
  design <- data$design; outcome <- data$outcome
  set.seed(615)
  init_coef <- rnorm(n_pred)
  wls_updated_coef <- 
    take_one_newton_step(init_coef, design, outcome, solver = "weighted-leqst-sq")
  ne_updated_coef <- 
    take_one_newton_step(init_coef, design, outcome, solver = "normal-eq")
  expect_true(are_all_close(wls_updated_coef, ne_updated_coef))
})
