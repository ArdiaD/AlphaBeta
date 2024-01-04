
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Eigen::Lower;

// [[Rcpp::export]]
VectorXd f_eigenvalues_eigen(Map<MatrixXd> M) {
  SelfAdjointEigenSolver<MatrixXd> es(M);
  return es.eigenvalues();
}

// [[Rcpp::export]]
MatrixXd f_crossprod_eigen(Map<MatrixXd> X) {

const int n(X.cols());
MatrixXd XtX(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(X.adjoint()));

return XtX;
}

// [[Rcpp::export]]
MatrixXd f_tcrossprod_eigen(Map<MatrixXd> X) {
  
  const int m(X.rows());
  MatrixXd XXt(MatrixXd(m, m).setZero().selfadjointView<Lower>().rankUpdate(X));
  
  return XXt;
}

