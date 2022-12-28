# Solve for the L2-norm minimizer of sqrt(weight) * (X %*% beta - y)
solve_least_sq_via_qr <- function(X, y, weight = NULL) {
  if (!is.null(weight)) {
    X <- outer(sqrt(weight), rep(1, ncol(X))) * X
    y <- sqrt(weight) * y
  }
  qr_decomp <- qr_wrapper(X) 
  Q <- qr_decomp$Q; R <- qr_decomp$R
  solution <- backsolve(R, t(Q) %*% y)
  solution <- as.vector(solution)
  return(list(solution = solution, Q = Q, R = R))
}

# Find an inverse of the gram matrix t(X) %*% X from the QR factor of X
invert_gram_mat_from_qr <- function(R) {
  return(chol2inv(R))
}

# Returns un-pivoted (Q, R) factors so that X == Q %*% R
qr_wrapper <- function(X) {
  qr_result <- qr(X)
  Q <- qr.Q(qr_result)
  R <- qr.R(qr_result)
  R <- R[ , order(qr_result$pivot)]
  return(list(Q = Q, R = R))
}