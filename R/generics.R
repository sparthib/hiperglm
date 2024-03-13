calc_loglik <- function(x, outcome, model_name, reg_coef, ...) {
  UseMethod("calc_loglik")
}

calc_loglik.matrix <- function(x, outcome, model_name, reg_coef, noise_var = 1) { 
  design <- x 
  
  if (model_name == "linear") {
    predicted_val <- design %*% reg_coef
    loglik <- -0.5 * sum((outcome - predicted_val)^2) / noise_var
  } else if (model_name == "logit") {
    if (is.list(outcome)) {
      n_success <- outcome$n_success
      n_trial <- outcome$n_trial
    } else {
      n_success <- outcome
      n_trial <- rep(1, length(n_success)) # Assume binary outcome
    }
    logit_prob <- design %*% reg_coef
    loglik <- sum(n_success * logit_prob - n_trial * log(1 + exp(logit_prob)))
  }
  return(loglik)
}

calc_loglink_deriv <- function(x, outcome, model_name, reg_coef, ...) {
  UseMethod("calc_loglink_deriv")
}

calc_loglink_deriv.matrix <- function(x, outcome, model_name, reg_coef, order = 1, ...) { 
  design <- x 
  if(model_name == "logit"){ 
    if (is.list(outcome)) {
      n_success <- outcome$n_success
      n_trial <- outcome$n_trial
    } else {
      n_success <- outcome
      n_trial <- rep(1, length(n_success)) # Assume binary outcome
    }
    logit_prob <- as.vector(design %*% reg_coef)
    predicted_prob <- 1 / (1 + exp(-logit_prob))
    if (order == 1) {
      deriv <- n_success - n_trial * predicted_prob
    } else if (order == 2) {
      deriv <- n_trial * predicted_prob * (1 - predicted_prob)
    } else {
      stop("3rd+ order derivative calculations are not supported")
    }
    deriv <- as.vector(deriv)
  }
  
  return(deriv)
}

calc_grad <- function(design, outcome, method, reg_coef, noise_var = 1, ...) { 
  if (method == "linear"){ 
    predicted_val <- design %*% reg_coef
    grad <- t(design) %*% (outcome - predicted_val) / noise_var
    grad <- as.vector(grad)
  }
  else if(method == "logit"){ 
    loglink_grad <- calc_logit_loglink_deriv(reg_coef, design, outcome, order = 1)
    grad <- t(design) %*% loglink_grad
    grad <- as.vector(grad)
  }
  return(grad)
}
