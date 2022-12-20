test_that("`hglm` linear model est coincides with `lm` one", {
  data <- simulate_data(32, 4, intercept = 1, seed = 1918)
  design <- data$design; outcome <- data$outcome
  lm_out <- stats::lm(outcome ~ design + 0)
  hglm_out <- hiper_glm(design, outcome)
  expect_true(are_all_close(coef(hglm_out), coef(lm_out)))
  expect_true(are_all_close(  
    as.vector(vcov(hglm_out)), as.vector(vcov(lm_out))
  ))
})

test_that("`hglm` logit model est coincides with `glm` one", {
  data <- simulate_data(32, 4, model = 'logit', seed = 1918)
  design <- data$design; outcome <- data$outcome
  glm_out <- stats::glm(outcome ~ design + 0, family = binomial('logit'))
  hglm_out <- hiper_glm(
    design, outcome, model = 'logit', option = list(mle_solver = 'BFGS')
  )
  expect_true(are_all_close(
    coef(hglm_out), coef(glm_out), abs_tol = Inf, rel_tol = 1e-3
  ))
  # TODO: test covariance estimate
})