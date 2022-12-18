test_that("`hglm` coef/vcov estimate coincides with `lm` one", {
  data <- simulate_data(32, 4, intercept = 1, seed = 1918)
  design <- data$design; outcome <- data$outcome
  lm_out <- stats::lm(outcome ~ design + 0)
  hglm_out <- hiper_glm(design, outcome)
  expect_true(are_all_close(coef(hglm_out), coef(lm_out)))
  expect_true(are_all_close(  
    as.vector(vcov(hglm_out)), as.vector(vcov(lm_out))
  ))
})