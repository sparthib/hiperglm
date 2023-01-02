# Solve for the L2-norm minimizer of sqrt(weight) * (X %*% beta - y)
solve_least_sq_via_qr <- function(X, y, weight = NULL, use_eigen = TRUE) {
  if (!is.null(weight)) {
    X <- outer(sqrt(weight), rep(1, ncol(X))) * X
    y <- sqrt(weight) * y
  }
  y <- as.numeric(y) # Ensure `double`
  if (use_eigen) {
    ls_result <- solve_least_sq_via_qr_cpp_eig(X, y)
  } else {
    out <- qr_wrapper(X)
    qr_decomp <- out$qr_decomp; R <- out$R
    solution <- as.vector(solve(qr_decomp, y))
    ls_result <- list(solution = solution, R = R)
  }
  return(ls_result)
}

# Find an inverse of the gram matrix t(X) %*% X from the QR factor of X
invert_gram_mat_from_qr <- function(R) {
  return(chol2inv(R))
}

# Return the QR decomposition along with the unpivoted "R" factor so that X == Q %*% R
qr_wrapper <- function(X) {
  qr_decomp <- qr(X)
  R <- qr.R(qr_decomp)
  R <- R[ , order(qr_decomp$pivot)]
  return(list(qr_decomp = qr_decomp, R = R))
}