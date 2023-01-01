# Solve for the L2-norm minimizer of sqrt(weight) * (X %*% beta - y)
solve_least_sq_via_qr <- function(X, y, weight = NULL) {
  if (!is.null(weight)) {
    X <- outer(sqrt(weight), rep(1, ncol(X))) * X
    y <- sqrt(weight) * y
  }
  out <- qr_wrapper(X)
  qr_decomp <- out$qr_decomp; R <- out$R
  solution <- solve(qr_decomp, y)
  solution <- as.vector(solution)
  return(list(solution = solution, R = R))
}

# Find an inverse of the gram matrix t(X) %*% X from the QR factor of X
invert_gram_mat_from_qr <- function(R) {
  return(chol2inv(R))
}

# Return un-pivoted (Q, R) factors so that X == Q %*% R
qr_wrapper <- function(X, require_Q = FALSE) {
  qr_decomp <- qr(X)
  R <- qr.R(qr_decomp)
  R <- R[ , order(qr_decomp$pivot)]
  result <- list(qr_decomp = qr_decomp, R = R)
  if (require_Q) {
    result$Q <- qr.Q(qr_decomp)
  }
  return(result)
}