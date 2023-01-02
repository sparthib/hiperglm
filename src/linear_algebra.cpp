#include "hiperglm_types.h"

// The `depends` attribute tells Rcpp to create hooks for RcppEigen
// [[Rcpp::depends(RcppEigen)]]

// Solve for the l2-norm minimizer of X %*% beta - y
// [[Rcpp::export]]
List solve_least_sq_via_qr_cpp_eig(
    const Map<MatrixXd> X, const Map<VectorXd> y
) {
  Eigen::HouseholderQR<MatrixXd> qr(X);
  VectorXd solution = qr.solve(y);
  int n_col = X.cols();
  MatrixXd R = qr.matrixQR().topRows(n_col).triangularView<Eigen::Upper>();
  return List::create(
    _["solution"] = solution,
    _["R"] = R
  );
}