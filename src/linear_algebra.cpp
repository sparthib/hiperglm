#include "hiperglm_types.h"

// The `depends` attribute tells Rcpp to create hooks for RcppEigen
// [[Rcpp::depends(RcppEigen)]]

// Solve for the l2-norm minimizer of X %*% beta - y
// [[Rcpp::export]]
List solve_leqst_sq_via_qr_cpp_eig(
    const Map<MatrixXd> X, const Map<VectorXd> y
) {
  Eigen::HouseholderQR<MatrixXd> qr(X);
  VectorXd solution = qr.solve(y);
  int n_col = X.cols();
  const Eigen::TriangularView<const Eigen::Block<const MatrixXd>, Eigen::Upper> R \
    = qr.matrixQR().topRows(n_col).triangularView<Eigen::Upper>();
  MatrixXd R_inv = R.solve(MatrixXd::Identity(n_col, n_col));
  MatrixXd inverse_gram_mat = R_inv * R_inv.transpose();
  return List::create(
    _["solution"] = solution,
    _["inverse_gram"] = inverse_gram_mat
  );
}