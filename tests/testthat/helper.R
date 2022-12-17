are_all_close <- function(v, w, abs_tol = 1e-6, rel_tol = 1e-6) {
  abs_diff <- abs(v - w)
  are_all_within_atol <- all(abs_diff < abs_tol)
  are_all_within_rtol <- all(abs_diff < rel_tol * pmax(abs(v), abs(w)))
  return(are_all_within_atol && are_all_within_rtol)
}

approx_grad_via_finite_diff <- function(func, x, dx = 1e-6) {
  numerical_grad <- rep(0, length(x))
  for (i in 1:length(x)) {
    x_plus <- x 
    x_plus[i] <- x[i] + dx
    x_minus <- x
    x_minus[i] <- x[i] - dx
    numerical_grad[i] <- (func(x_plus) - func(x_minus)) / (2 * dx)
  }
  return(numerical_grad)
}

simulate_data <- function(
    n_obs, n_pred, model = "linear", intercept = NULL, 
    coef_true = NULL, design = NULL, seed = NULL, signal_to_noise = NULL
  ) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if ((model != "linear")  && !is.null(signal_to_noise)) {
    warning(paste(
      "The `signal_to_noise` arg is currently unsupported for non-linear models",
      "and will be ignored."
    ))
    signal_to_noise <- NULL
  }
  if (is.null(signal_to_noise)) {
    signal_to_noise <- 0.1
  } 
  if (is.null(coef_true)) {
    coef_true <- rnorm(n_pred, sd = 1 / sqrt(n_pred))
  } 
  if (is.null(design)) {
    design <- matrix(rnorm(n_obs * n_pred), nrow = n_obs, ncol = n_pred)
  }
  if (!is.null(intercept)) {
    if (!is.numeric(intercept)) {
      stop("The intercept argument must be numeric.")
    }
    coef_true <- c(intercept, coef_true)
    design <- cbind(rep(1, n_obs), design)
  }
  expected_mean <- as.vector(design %*% coef_true)
  if (model == 'linear') {
    noise_magnitude <- sqrt(var(expected_mean) / signal_to_noise^2)
    noise <- noise_magnitude * rnorm(n_obs)
    outcome <- expected_mean + noise 
  } else {
    n_trial <- 1
    prob <- 1 / (1 + exp(-expected_mean))
    outcome <- rbinom(n_obs, n_trial, prob)
  }
  return(list(design = design, outcome = outcome, coef_true = coef_true))
}