calc_linear_loglik <- function(reg_coef, design, outcome, noise_var = 1) {
  predicted_val <- design %*% reg_coef
  loglik <- - 0.5 * sum((outcome - predicted_val)^2) / noise_var
  return(loglik)
}

calc_linear_grad <- function(reg_coef, design, outcome, noise_var = 1) {
  predicted_val <- design %*% reg_coef
  grad <- t(design) %*% (outcome - predicted_val)
  grad <- as.vector(grad)
  return(grad)
}

calc_logit_loglik <- function(
    reg_coef, design, n_success, n_trial = NULL
) {
  if (is.null(n_trial)) {
    # Assume binary outcome unless otherwise specified.
    n_trial <- rep(1, length(n_success))
  }
  logit_prob <- design %*% reg_coef
  loglik <- sum(n_success * logit_prob - n_trial * log(1 + exp(logit_prob)))
    # TODO: improve numerical stability for logit_prob >> 1
  return(loglik)
}

calc_logit_grad <- function(
    reg_coef, design, n_success, n_trial = NULL
) {
  if (is.null(n_trial)) {
    n_trial <- rep(1, length(n_success))
  }
  logit_prob <- design %*% reg_coef
  predicted_prob <- 1 / (1 + exp(-logit_prob))
  grad <- t(design) %*% (n_success - n_trial * predicted_prob)
  grad <- as.vector(grad)
  return(grad)
}

calc_logit_hessian <- function(
    reg_coef, design, n_success, n_trial = NULL
) {
  if (is.null(n_trial)) {
    n_trial <- rep(1, length(n_success))
  }
  logit_prob <- as.vector(design %*% reg_coef)
  predicted_prob <- 1 / (1 + exp(-logit_prob))
  weight <- n_trial * predicted_prob * (1 - predicted_prob)
  hess <- - t(design) %*% (outer(weight, rep(1, ncol(design))) * design)
  return(hess)
}