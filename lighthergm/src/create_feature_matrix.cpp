// Files to look to get MM version: ReciprocityModel for basic functions,
// BinaryReciprocityModel.cpp for version written already,
// MMBinaryReciprocityModel.cpp for more eleborate version
// #define ARMA_64BIT_WORD 1;
#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_max_threads() 0
#endif
#include <RcppArmadillo.h>
#include "helper.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


// [[Rcpp::export]]
arma::sp_mat get_sparse_feature_adjmat(const arma::vec& x) {
  int numOfVertices = x.size();
  arma::sp_mat S(numOfVertices, numOfVertices);
  // When x[i] == x[j] and i != j, S[i,j] = 1.
  #pragma omp parallel for
  for (int j = 0; j < numOfVertices; j++) {
    for (int i = 0; i < numOfVertices; i++) {
      if (i != j) {
        if (x[i] == x[j]) {
          S(i, j) = 1;
        }
      }
    }
  }
  return S;
}

// [[Rcpp::export]]
arma::sp_mat get_sparse_feature_adjmat_from_string(const Rcpp::StringVector& x) {
  int numOfVertices = x.size();
  arma::sp_mat S(numOfVertices, numOfVertices);
  // When x[i] == x[j] and i != j, S[i,j] = 1.
  #pragma omp parallel for
  for (int j = 0; j < numOfVertices; j++) {
    for (int i = 0; i < numOfVertices; i++) {
      if (i != j) {
        if (x[i] == x[j]) {
          S(i, j) = 1;
        }
      }
    }
  }
  return S;
}

// Compute something like X := - (S + T + U) + (S % T + T % U + U % S) - S % T % U.
// [[Rcpp::export]]
arma::sp_mat get_matrix_for_denominator(int numOfVertices, const Rcpp::List& list_feature_adjmat)
{
  int n_feature = list_feature_adjmat.length();
  int n_item = pow(2, n_feature);
  arma::sp_mat output(numOfVertices, numOfVertices);

  // It is difficult to explain this part...
  for (int s = 1; s < n_item; s++) {
    arma::vec index = decimal_to_binary_vector(s, n_feature);
    int k = sum(index);
    arma::sp_mat X(numOfVertices, numOfVertices);
    // Set a counter
    int counter = 0;
    for (int t = 0; t < n_feature; t++) {
      if (index[t] == 1) {
        arma::sp_mat S = list_feature_adjmat[t];
        counter += 1;
        if (counter == 1) {
          X = S;
        } else {
          X = X % S;
        }
      }
    }
    // X = arma::trimatu(X);
    output += pow(-1, k) * X;
  }
  return output;
}


// [[Rcpp::export]]
Rcpp::List get_elementwise_multiplied_matrices(const arma::sp_mat& adjmat,
                                               const Rcpp::List& list_feature_adjmat) {
  // Number of nodes
  int n_node = adjmat.n_rows;
  // Number of feature matrices
  int n_feature = list_feature_adjmat.length();

  // Append all the matrices in a single list
  Rcpp::List list_mat(n_feature+1);
  list_mat[0] = adjmat;
  for (int i = 0; i < n_feature; i++) {
    list_mat[i+1] = list_feature_adjmat[i];
  }

  // Create a list to store multiplied matrices
  int n_matrix = list_mat.length();
  int length_output = pow(2, n_matrix);
  Rcpp::List output(length_output);

  // The first element of the output list should contain the matrix for the denominator of pi_d0x0.
  arma::sp_mat denom = get_matrix_for_denominator(n_node, list_feature_adjmat);
  output[0] = denom;

  // Element-wise matrix multiplication without breaking sparsity
  for (int s = 1; s < length_output; s++) {
    // Convert an integer to a binary numeric vector
    arma::vec index = decimal_to_binary_vector(s, n_matrix);
    // Initialize a sparse matrix
    arma::sp_mat X(n_node, n_node);
    // Set a counter
    int counter = 0;

    // Start element-wise matrix multiplication
    for (int t = 0; t < n_matrix; t++) {
      // Prepare a matrix to be multiplied
      arma::sp_mat S = list_mat[t];
      // First, multiply matrices that don't need subtraction like (one - X).
      if (index[t] == 1) {
        counter += 1;
        if (counter == 1) {
          X = S;
        } else {
          X = X % S;
        }
      }
    }

    // Then multiply the rest of the matrices
    for (int t = 0; t < n_matrix; t++) {
      // Prepare a matrix to be multiplied
      arma::sp_mat S = list_mat[t];
      if (index[t] == 0) {
        X = X - S % X;
      }
    }
    // Lastly, store the multiplied matrix in the output list
    output[s] = X;
  }
  // Return the output
  return output;
}

