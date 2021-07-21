#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec decimal_to_binary_vector(int decimal,
                                   int vec_length) {
  int n = decimal;
  arma::vec output(vec_length);

  for (int i = 0; i < vec_length; i++) {
    output[i] = n % 2;
    n = n/2;
  }
  return output;
}

// Summing up a matrix by row
arma::rowvec sumDoubleMatrixByRow(const arma::mat& matrix) {
  arma::rowvec vector = arma::sum(matrix, 0);
  return vector;
}


// Normalizing tau
void normalizeTau(arma::mat& tau,
                  double minValue) {
  int numOfVertices = tau.n_rows;
  int numOfClasses = tau.n_cols;
  // normalize
  for (int i = 0; i < numOfVertices; i++) {
    double denominator = 0;
    for (int k = 0; k < numOfClasses; k++) {
      denominator += tau(i, k);
    }
    bool again = false;
    for (int k = 0; k < numOfClasses; k++) {
      tau(i, k) /= denominator;
      if (tau(i, k) < minValue) {
        tau(i, k) = minValue;
        again = true;
      }
    }
    if (again) {
      denominator = 0;
      for (int k = 0; k < numOfClasses; k++) {
        denominator += tau(i, k);
      }
      for (int k = 0; k < numOfClasses; k++)
        tau(i, k) /= denominator;
    }
  }
}


// [[Rcpp::export]]
arma::mat compute_sumTaus(int numOfVertices,
                          int numOfClasses,
                          const arma::mat& tau,
                          int verbose = 0) {

  if (verbose >= 5) {
    Rcpp::Rcout  << "find_sumTaus: sum by row";
  }
  arma::rowvec tauL = sumDoubleMatrixByRow(tau);

  if (verbose >= 5) {
    Rcpp::Rcout  << "find_sumTaus: calculating sumTaus";
  }
  arma::mat sumTaus = tau.t() * -(tau.each_row() - tauL); // Check computation speed here

  if (verbose >= 5) {
    Rcpp::Rcout  << "find_sumTaus: returning";
  }

  return sumTaus;
}

// A naive implementation of quadratic coefficient computation
// [[Rcpp::export]]
arma::mat compute_quadratic_term_naive(int numOfVertices,
                                       int numOfClasses,
                                       const arma::mat& pi,
                                       const arma::mat& tau,
                                       const arma::sp_mat& network) {
  arma::mat pi1 = pi;
  arma::mat pi0 = 1 - pi;
  arma::mat logPi0 = arma::log(pi0);
  arma::mat logPi1 = arma::log(pi1);

  arma::mat A(numOfVertices, numOfClasses);
  A.zeros();
  for (int i = 0; i < numOfVertices; i++) {
    for (int k = 0; k < numOfClasses; k++) {
      for (int j = 0; j < numOfVertices; j++) {
        if (i != j) {
          for (int l = 0; l < numOfClasses; l++) {
            if (network(i, j) == 0 ) {
              A(i, k) += tau(j, l) * logPi0(k, l);
            } else {
              A(i, k) += tau(j, l) * logPi1(k, l);
            }
          }
        }
      }
    }
  }

  // Finalize by subtracting half of from 1 dividing tau_{ik}
  for (int i = 0; i < numOfVertices; i++) {
    for (int k = 0; k < numOfClasses; k++) {
      // In theory, A(i, k) must be negative or 0.
      if (A(i, k) > 0) { // In reality, A(i, k) can be greater than 0 because of numerical precision.
        A(i, k) = 0;     // Therefore, we cut it off to 0 in this case
      }
      A(i, k) = 1 - A(i, k) / 2;
      A(i, k) /= tau(i, k);
    }
  }
  return A;
}

