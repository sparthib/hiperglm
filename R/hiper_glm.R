#' @export
hiper_glm <- function(design, outcome, model = "linear", option = list()) {
  supported_model <- c("linear", "logit")
  if (!(model %in% supported_model)) {
    stop(sprintf("The model %s is not supported.", model))
  }
  warning("`hiper_glm` is yet to be implemented.")
  mle <- find_mle(design, outcome, model, option)
  hglm_out <- list(coef = mle$coef, info_mat = mle$info_mat)
  class(hglm_out) <- "hglm"
  return(hglm_out)
}

find_mle <- function(design, outcome, model, option) {
  if (model == 'linear') {
    if (is.null(option$mle_solver)) {
      result <- solve_via_least_sq(design, outcome)
    } else {
      result <- solve_via_optim(design, outcome, option$mle_solver)
    }
  } else {
    # TODO: implement iteratively reweighted least-sq
    stop("Not yet implemented.")
  }
  return(result)
}

solve_via_least_sq <- function(design, outcome) {
  XTX <- t(design) %*% design 
  mle_coef <- solve(XTX, t(design) %*% outcome)
  mle_coef <- as.vector(mle_coef)
  noise_var <- mean((outcome - design %*% mle_coef)^2)
  n_obs <- nrow(design); n_pred <- ncol(design)
  noise_var <- noise_var / (1 - n_pred / n_obs) 
  # Use the same nearly-unbiased estimator as in `stats::lm`
  info_mat <- XTX / noise_var
  return(list(coef = mle_coef, info_mat = info_mat))
}

solve_via_optim <- function(design, outcome, method) {
  init_coef <- rep(0, ncol(design))
  obj_fn <- function (coef) {
    calc_linear_loglik(coef, design, outcome) 
  }
  obj_grad <- function (coef) {
    calc_linear_grad(coef, design, outcome)
  }
  optim_result <- stats::optim(
    init_coef, obj_fn, obj_grad, method = method,
    control = list(fnscale = -1) # Maximize the function
  )
  optim_converged <- (optim_result$convergence == 0L)
  if (!optim_converged) {
    warning("Optimization did not converge. The estimates may be meaningless.")
  }
  return(list(coef = optim_result$par))
}
