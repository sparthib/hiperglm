test_that("`hglm` coef estimate coincides with `lm` one", {
  data <- simulate_data(32, 4, intercept = 1, seed = 1918)
  lm_out <- stats::lm(data$outcome ~ data$design + 0)
  hglm_out <- hiper_glm(data$design, data$outcome)
  expect_true(are_all_close(coef(hglm_out), coef(lm_out)))
})