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
  if (is.null(option$mle_solver)) {
    if (model == 'linear') {
      result <- solve_via_least_sq(design, outcome)
    } else {
      result <- solve_via_newton(
        design, outcome, option$n_max_iter, option$rel_tol, option$abs_tol
      )
    }
  } else {
    result <- solve_via_optim(design, outcome, model, option$mle_solver)
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

solve_via_newton <- function(design, outcome, n_max_iter, rel_tol, abs_tol) {
  if (is.null(n_max_iter)) { n_max_iter <- 25L }
  if (is.null(rel_tol)) { rel_tol <- 1e-6 }
  if (is.null(abs_tol)) { abs_tol <- 1e-6 }
  coef_est <- rep(0, ncol(design))
  n_iter <- 0L
  max_iter_reached <- FALSE
  converged <- FALSE
  curr_loglik <- calc_logit_loglik(coef_est, design, outcome)
  while (!(converged || max_iter_reached)) {
    prev_loglik <- curr_loglik
    grad <- calc_logit_grad(coef_est, design, outcome)
    hess <- calc_logit_hessian(coef_est, design, outcome)
    coef_est <- coef_est - solve(hess, grad)
    curr_loglik <- calc_logit_loglik(coef_est, design, outcome)
    converged <- (
      2 * abs(curr_loglik - prev_loglik) < (abs_tol + rel_tol * abs(curr_loglik))
    )
    n_iter <- n_iter + 1L
    max_iter_reached <- (n_iter == n_max_iter)
  }
  if (max_iter_reached && !converged) {
    warning("Newton's method did not converge. The estimates may be meaningless.")
  }
  return(list(coef = coef_est, info_mat = - hess, converged = converged))
}

solve_via_optim <- function(design, outcome, model, method) {
  init_coef <- rep(0, ncol(design))
  if (model == 'linear') {
    obj_fn <- function (coef) {
      calc_linear_loglik(coef, design, outcome) 
    }
    obj_grad <- function (coef) {
      calc_linear_grad(coef, design, outcome)
    }
  } else {
    obj_fn <- function (coef) {
      calc_logit_loglik(coef, design, outcome) 
    }
    obj_grad <- function (coef) {
      calc_logit_grad(coef, design, outcome)
    }
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
