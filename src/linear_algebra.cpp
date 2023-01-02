#include "hiperglm_types.h"

// The `depends` attribute tells Rcpp to create hooks for RcppEigen
// [[Rcpp::depends(RcppEigen)]]

// Solve for the l2-norm minimizer of X %*% beta - y
// [[Rcpp::export]]
List solve_least_sq_via_qr_cpp_eig(
    const Map<MatrixXd> X, const Map<VectorXd> y, 
    bool require_inverse_gram = false
) {
  Eigen::HouseholderQR<MatrixXd> qr(X);
  VectorXd solution = qr.solve(y);
  MatrixXd inverse_gram_mat;
  if (require_inverse_gram) {
    int n_col = X.cols();
    const Eigen::TriangularView<const Eigen::Block<const MatrixXd>, Eigen::Upper> R \
      = qr.matrixQR().topRows(n_col).triangularView<Eigen::Upper>();
    MatrixXd R_inv = R.solve(MatrixXd::Identity(n_col, n_col));
    inverse_gram_mat = R_inv * R_inv.transpose();
  } else {
    inverse_gram_mat = MatrixXd::Zero(0, 0);
  }
  return List::create(
    _["solution"] = solution,
    _["inverse_gram"] = inverse_gram_mat
  );
}