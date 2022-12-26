# Solve for the L2-norm minimizer of X %*% beta - y
solve_least_sq_via_qr <- function(X, y) {
  qr_decomp <- qr_wrapper(X) 
  Q <- qr_decomp$Q; R <- qr_decomp$R
  solution <- backsolve(R, t(Q) %*% y)
  solution <- as.vector(solution)
  return(list(solution = solution, Q = Q, R = R))
}

# Returns un-pivoted (Q, R) factors so that X == Q %*% R
qr_wrapper <- function(X) {
  qr_result <- qr(X)
  Q <- qr.Q(qr_result)
  R <- qr.R(qr_result)
  R <- R[ , order(qr_result$pivot)]
  return(list(Q = Q, R = R))
}