// Files to look to get MM version: ReciprocityModel for basic functions,
// BinaryReciprocityModel.cpp for version written already,
// MMBinaryReciprocityModel.cpp for more eleborate version
// #define ARMA_64BIT_WORD 1;
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat eigenvectors_sparse(
    const arma::sp_mat& X,
    int n_vec) {
  arma::vec eigval;
  arma::mat eigvec;

  arma::eigs_sym(eigval, eigvec, X, n_vec);
  return(eigvec);
}
