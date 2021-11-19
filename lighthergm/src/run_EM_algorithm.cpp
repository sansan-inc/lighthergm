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

const double minTau = 1e-6;
double minPi =  1e-4;


// Quadratic programming problem solver
void solveQP
  (const arma::mat& m,
   const arma::mat& s,
   arma::mat& tau,
   double precision,
   int verbose = 0) {

  if (m.n_rows <= 0)
    return;

  int n = m.n_cols;


  bool * J_k = new bool[n]();
  bool * J_lambda_a = new bool[n]();

  bool * J_lambda_k = new bool[n]();
  bool * J_lambda_k_a = new bool[n]();
  bool * J_lambda_k_b = new bool[n]();

  double lambda_k;

  int nrows = m.n_rows;

  if (verbose >= 5) {
    Rcpp::Rcout << "solveQP: about to enter QP iteration" << "\n";
  }

  for (int rowIndex = 0; rowIndex < nrows; rowIndex++) {
    // Step 1:
    for (int j = 0; j < n; j++) {
      J_k[j] = true;
      J_lambda_a[j] = false;
    }

    bool done = false;
    // we never loop over n times.
    int count = 0;
    while (!done) {
      count++;

      // Step 2:
      double value1 = 0, value2 = 0;
      for (int j = 0; j < n; j++) {
        if (J_k[j]) {
          value1 += 1 / m(rowIndex, j);
          value2 += s(rowIndex, j) / m(rowIndex, j);
        }
      }
      lambda_k = (value2 - 2) / value1;

        // Step3:
        for (int j = 0; j < n; j++) {
          if (J_k[j]) {
            if (lambda_k >= s(rowIndex, j)) {
              J_lambda_k[j] = false;
              J_lambda_k_a[j] = true;
              J_lambda_k_b[j] = false;
            } else if (lambda_k
                         > (-2 * m(rowIndex, j) + s(rowIndex, j))) {
              J_lambda_k[j] = true;
              J_lambda_k_a[j] = false;
              J_lambda_k_b[j] = false;
            } else {
              J_lambda_k[j] = false;
              J_lambda_k_a[j] = false;
              J_lambda_k_b[j] = true;
            }
          } else {
            J_lambda_k[j] = false;
            J_lambda_k_a[j] = false;
            J_lambda_k_b[j] = false;
          }
        }

        // Step 4:
        double delta = 0;
        double value3 = 0, value4 = 0;
        bool J_lambda_k_empty = true;
        for (int j = 0; j < n; j++) {
          delta += J_lambda_k_b[j];
          if (J_lambda_k[j]) {
            value3 += s(rowIndex, j) / m(rowIndex, j);
            value4 += 1 / m(rowIndex, j);
            J_lambda_k_empty = false;
          }
        }
        delta += (value3 / 2 - lambda_k * value4 / 2 - 1);

        if (fabs(delta) < precision || J_lambda_k_empty || count >= n) {
          for (int j = 0; j < n; j++) {
            if (J_lambda_a[j] || J_lambda_k_a[j])
              tau(rowIndex, j) = 0;
            else if (J_lambda_k_b[j])
              tau(rowIndex, j) = 1;
            else
              tau(rowIndex, j) = (s(rowIndex, j) - lambda_k) / (2 * m(rowIndex, j));
          }
          done = true;
        } else if (delta > 0) {
          for (int j = 0; j < n; j++) {
            if (J_lambda_k_a[j]) {
              J_lambda_a[j] = true;
              J_k[j] = false;
            }
          }
        } else {
          for (int j = 0; j < n; j++) {
            if (J_lambda_k_b[j])
              tau(rowIndex, j) = 1;
            else
              tau(rowIndex, j) = 0;
          }
          done = true;
        }
    }
  }
  // Don't forget to delete the bool arrays
  delete [] J_k;
  delete [] J_lambda_a;

  delete [] J_lambda_k;
  delete [] J_lambda_k_a;
  delete [] J_lambda_k_b;

}

// Compute coefficients on the linear terms of the surrogate function
// [[Rcpp::export]]
arma::mat compute_linear_term(
    int numOfVertices,
    int numOfClasses,
    const arma::vec& alpha,
    const arma::mat& tau,
    double& LB) {
  arma::mat s(numOfVertices, numOfClasses);
  // Calculate the linear coefficients
  arma::vec logAlpha = arma::log(alpha);
  #pragma omp parallel for
  for (int k = 0; k < numOfClasses; k++) {
    for (int i = 0; i < numOfVertices; i++) {
      s(i, k) = 1 + logAlpha(k) - log(tau(i, k));
    }
  }
  // Compute the lower bound
  #pragma omp parallel for reduction(+:LB)
  for (int k = 0; k < numOfClasses; k++) {
    for (int i = 0; i < numOfVertices; i++) {
      LB += tau(i, k) * (logAlpha(k) - log(tau(i, k)));
    }
  }
  return s;
}


// Compute \pi matrix without features
// [[Rcpp::export]]
arma::mat compute_pi(
    int numOfVertices,
    int numOfClasses,
    const arma::sp_mat& stat,
    const arma::mat& tau)
{
  // Compute pi
  arma::mat sumTaus = compute_sumTaus(numOfVertices, numOfClasses, tau);
  arma::mat pi = tau.t() * (stat * tau);
  pi /= sumTaus;
  // Remove extremely small elements in the denominator
  for (auto& val : pi) {
    if (val < minPi) {
      val = minPi;
    }
  }
  return pi;
}

// Compute coefficients on the quadratic terms of the surrogate function without features
// [[Rcpp::export]]
arma::mat compute_quadratic_term(
    int numOfVertices,
    int numOfClasses,
    const arma::vec& alpha,
    const arma::mat& tau,
    const arma::sp_mat& network,
    double& LB,
    int verbose = 0) {

  //Calculate pi's for reciprocity model
  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: calculating pi11";
  }

  arma::mat sumTaus = compute_sumTaus(numOfVertices, numOfClasses, tau);

  arma::mat pi11 = (tau.t() * network * tau) / sumTaus;
  // Remove extremely small elements in the denominator
  for (auto& val : pi11) {
    if (val < minPi) {
      val = minPi;
    }
  }

  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: calculating pi00";
  }
  arma::mat pi00 = 1 - pi11;
  // Remove extremely small elements in the denominator
  for (auto& val : pi00) {
    if (val < minPi) {
      val = minPi;
    }
  }

  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: calculating logPi00";
  }
  arma::mat logPi00 = arma::log(pi00);

  // Calculate the quadratic coefficients
  // Compute the norm term, i.e. \pi_kl^0
  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: calculating tauL";
  }
  arma::rowvec tauL = sumDoubleMatrixByRow(tau);

  // When D_ij = 0
  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: calculating A";
  }
  arma::mat tau_nought = -(tau.each_row() - tauL);
  arma::mat A = (logPi00 * (tau_nought.t())).t();
  LB += accu((tau.t() * tau_nought) % logPi00);

  // When D_ij = 1
  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: calculating logPi11";
  }
  arma::mat logPi11 = arma::log(pi11/pi00);

  // If the matrix is symmetric, only need to iterate over the positive values,
  // because the condition if (network(i,j) == 1 && network(j,i) == 1) will always be met.
  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: updating A";
  }
  arma::mat tauG = tau.t() * network;
  A += (logPi11 * tauG).t();
  LB += accu((tauG * tau) % logPi11);

  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: subtract from A";
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

// A wrapper function to implement the EM algorithm without features
// [[Rcpp::export]]
Rcpp::List run_EM_without_features(
    int numOfVertices,
    int numOfClasses,
    const arma::vec& alpha,
    arma::mat& tau,
    const arma::sp_mat& network,
    int verbose = 0) {

  // Lower bound
  double LB = 0;
  // Compute quadratic term
  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: compute_quadratic_term";
  }
  arma::mat A = compute_quadratic_term(numOfVertices, numOfClasses, alpha, tau, network, LB, verbose);

  // Compute linear term
  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: compute_linear_term";
  }
  arma::mat s = compute_linear_term(numOfVertices, numOfClasses, alpha, tau, LB);

  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: solveQP";
  }
  solveQP(A, s, tau, minTau, verbose);

  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: normalizeTau";
  }
  normalizeTau(tau, minTau);

  if (verbose >= 5) {
    Rcpp::Rcout  << "runFixedPointEstimationEStepMM_sparse: returning";
  }
  Rcpp::List output(2);
  output[0] = tau;
  output[1] = LB;
  return output;
}


// A function that computes the denominator of \pi_{d=1, x=0}.
// [[Rcpp::export]]
arma::mat compute_denominator_for_pi_d1x0(
    int numOfVertices,
    double numOfClasses,
    const arma::sp_mat& matrix_for_denominator,
    const arma::mat& tau,
    int verbose
    ) {
  // We would like to compute something like t(alpha) * (one - S) % (one - T) % (one - U) * alpha without breaking matrix sparsity.
  // Note that (one - S) % (one - T) % (one - U) = one - (S + T + U) + (S % T + T % U + U % S) - S % T % U
  // and that the diagonals of the matrix must be zeros.
  // Utilizing this, we compute as follows.
  // 1a. Get a matrix whose elements in jth column is sum of alpha's ith column. Call it sum_alpha.
  // 1b. Compute A1 := t(sum_alpha) * alpha. This corresponds to t(alpha) * one * alpha
  // 2a. Compute X := - (S + T + U) + (S % T + T % U + U % S) - S % T % U.
  // 2a. Compute A2 := t(alpha) * X * alpha
  //  3. Compute A1 + A2, which equals t(alpha) * (one - S) % (one - T) % (one - U) * alpha.
  // Below we write a code that implements the above.

  // Step 1: Compute A1 = t(sum_alpha) * alpha.
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_denominator_for_pi_d1x0: computing A1." << "\n";
  }
  arma::rowvec tauL = sumDoubleMatrixByRow(tau);
  arma::mat A1 = tau.t() * -(tau.each_row() - tauL);

  // Step 2: Compute A2 = t(tau) * (- (S + T + U) + (S % T + T % U + U % S) - S % T % U) * alpha.
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_denominator_for_pi_d1x0: computing A2." << "\n";
  }
  /*
  arma::mat tau_copy(tau);
  for (auto& val : tau_copy) {
    if (val <= 1 / (998 + numOfClasses)) {
      val = 0;
    }
  }
  arma::sp_mat tau_sp(tau_copy);
  */
  arma::mat A2 = tau.t() * matrix_for_denominator * tau;
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_denominator_for_pi_d1x0: computing A1 + A2." << "\n";
  }
  arma::mat output = A1 + A2;
  // Remove extremely small elements in the denominator
  for (auto& val : output) {
    if (val < 1e-10) {
      val = 1;
    }
  }
  // Return the output
  return output;
}


// A function that computes \pi_{d=1, x=0}.
// [[Rcpp::export]]
arma::mat compute_pi_d1x0(
    int numOfVertices,
    int numOfClasses,
    const Rcpp::List& list_multiplied_feature_adjmat,
    const arma::mat& tau,
    int verbose) {

  arma::mat tau_copy(tau);
  for (auto& val : tau_copy) {
    if (val <= 1 / (998 + numOfClasses)) {
      val = 0;
    }
  }
  arma::sp_mat tau_sp(tau_copy);

  // Compute the denominator
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_pi_d1x0: started computing the denominator of pi_d1x0." << "\n";
  }
  arma::mat sumTaus_feature = compute_denominator_for_pi_d1x0(numOfVertices, numOfClasses, list_multiplied_feature_adjmat[0], tau, verbose);

  // Compute Pi_{d=1;X=0}
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_pi_d1x0: started computing pi_d1x0." << "\n";
  }

  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: multiplying matrices."<< "\n";
  }
  arma::sp_mat S = list_multiplied_feature_adjmat[1];
  arma::mat pi_d1x0 = (tau.t() * S * tau) / sumTaus_feature;

  // Remove extremely small elements in pi_d0x0
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: removing extremely small elements in pi_d0x0."<< "\n";
  }
  for (auto& val : pi_d1x0) {
    if (val < minPi) {
      val = minPi;
      }
    }
  // Return the output
  return pi_d1x0;
}

// Compute coefficients on the quadratic terms of the surrogate function
// [[Rcpp::export]]
arma::mat compute_quadratic_term_with_features(
    int numOfVertices,
    int numOfClasses,
    const Rcpp::List& list_multiplied_feature_adjmat,
    const arma::mat& tau,
    double& LB,
    int verbose = 0)
{
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: started computing the quadratic coeffcients." << "\n";
  }

  // Calculate pi_d1x0
  arma::mat pi_d1x0 = compute_pi_d1x0(numOfVertices, numOfClasses, list_multiplied_feature_adjmat, tau, verbose);

  // Compute pi_d0x0
  arma::mat pi_d0x0 = 1 - pi_d1x0;
  // Remove extremely small elements in pi_d0x0
  for (auto& val : pi_d0x0) {
    if (val < minPi) {
      val = minPi;
    }
  }

  // When D_ij = 0 and X_ij = 0
  arma::mat logPi_d0x0 = arma::log(pi_d0x0);
  arma::rowvec tauL = sumDoubleMatrixByRow(tau);
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: computing the baseline quadratic coeffcient." << "\n";
  }

  arma::mat tau_nought = -(tau.each_row() - tauL);
  arma::mat A = (logPi_d0x0 * (tau_nought.t())).t();

  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: computing the lower bound for D_ij = 0 and X_ij = 0." << "\n";
  }
  LB += accu((tau.t() * tau_nought) % logPi_d0x0);

  // When D_ij = 1 and X_ij = 0
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: computing the quadratic coeffcient for D_ij = 1 and X_ij = 0." << "\n";
  }
  arma::mat logPi_d1x0 = arma::log(pi_d1x0/pi_d0x0);
  arma::sp_mat S0 = list_multiplied_feature_adjmat[1];
  arma::mat tau_one = tau.t() * S0;
  A += (logPi_d1x0 * tau_one).t();
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: computing the lower bound for D_ij = 1 and X_ij = 0." << "\n";
  }
  LB += accu((tau_one * tau) % logPi_d1x0);

  // Update A for the other combinations of features
  double length_list = list_multiplied_feature_adjmat.length();
  int numOfSteps = length_list/2;

  for (int s = 1; s < numOfSteps; s++) {
    // First, compute the conditional probability for D_ij = 0 and D_ij = 1 given X_ij = x
    // Create indices to extract appropriate multiplied matrices
    int index = 2 * s;
    // Get the multiplied matrix for D = 0
    arma::sp_mat D0X = list_multiplied_feature_adjmat[index];
    // Get the multiplied matrix for D = 1
    arma::sp_mat D1X = list_multiplied_feature_adjmat[index+1];
    // Compute the denominator for pi
    if (verbose >= 5) {
      auto time = std::chrono::system_clock::now();
      std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
      Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: computing D0X + D1X in loop " << s << " of " << numOfSteps<< "\n";
    }
    arma::sp_mat X = D0X + D1X;
    if (verbose >= 5) {
      auto time = std::chrono::system_clock::now();
      std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
      Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: computing pi_denom in loop " << s << " of " << numOfSteps << "\n";
    }
    arma::mat pi_denom = tau.t() * X * tau;
    // Compute the numerator for D_ij = 1
    if (verbose >= 5) {
      auto time = std::chrono::system_clock::now();
      std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
      Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: computing the numerator for D_ij = 1 in loop " << s << " of " << numOfSteps << "\n";
    }
    arma::mat pi_num = tau.t() * D1X * tau;
    // Compute pi for D_ij = 1
    if (verbose >= 5) {
      auto time = std::chrono::system_clock::now();
      std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
      Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: computing pi for D_ij = 1 in loop " << s << " of " << numOfSteps << "\n";
    }
    arma::mat pi_d1 = pi_num/pi_denom;
    // Remove extremely small elements in pi_d1
    // Remove extremely small elements in the denominator
    for (auto& val : pi_d1) {
      if (val < minPi) {
        val = minPi;
      }
    }
    // Compute pi for D_ij = 0
    if (verbose >= 5) {
      auto time = std::chrono::system_clock::now();
      std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
      Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: computing pi for D_ij = 0 in loop " << s << " of " << numOfSteps << "\n";
    }
    arma::mat pi_d0 = 1 - pi_d1;
    // Remove extremely small elements in pi_d0
    // Remove extremely small elements in the denominator
    for (auto& val : pi_d0) {
      if (val < minPi) {
        val = minPi;
      }
    }

    // Second, using the computed pis above, update the quadratic term
    // When D_ij = 0
    if (verbose >= 5) {
      auto time = std::chrono::system_clock::now();
      std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
      Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: updating A for D_ij = 0 in loop " << s << " of " << numOfSteps << "\n";
    }
    arma::mat logPi_d0 = arma::log(pi_d0/pi_d0x0);
    A += (logPi_d0 * (tau.t() * D0X)).t();
    if (verbose >= 5) {
      auto time = std::chrono::system_clock::now();
      std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
      Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: computing the lower bound D_ij = 0 in loop " << s << " of " << numOfSteps << "\n";
    }
    arma::mat tauX0tau = pi_denom - pi_num;
    LB += accu(tauX0tau % logPi_d0);

    // When D_ij = 1
    if (verbose >= 5) {
      auto time = std::chrono::system_clock::now();
      std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
      Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: updating A for D_ij = 1 in loop " << s << " of " << numOfSteps << "\n";
    }
    arma::mat logPi_d1 = arma::log(pi_d1/pi_d0x0);
    A += (logPi_d1 * (tau.t() * D1X)).t();
    if (verbose >= 5) {
      auto time = std::chrono::system_clock::now();
      std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
      Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: computing the lower bound D_ij = 1 in loop " << s << " of " << numOfSteps << "\n";
    }
    LB += accu(pi_num % logPi_d1);
  }

  // Finalize by subtracting half of from 1 dividing tau_{ik}
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: Finalizing A."<< "\n";
  }
  for (auto& val : A) {
    // In reality, A(i, k) can be greater than 0 because of numerical precision.
    // Therefore, we cut it off to 0 in this case
    if (val > 0) {
      val = 0;
    }
    val = 1 - val/2;
  }
  A /= tau;

  // Return the quadratic term
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "compute_quadratic_term_with_features: Finished."<< "\n";
  }
  return A;
}


// Compute \pi matrix using features
// [[Rcpp::export]]
Rcpp::List compute_pi_with_features(
    int numOfVertices,
    int numOfClasses,
    const Rcpp::List& list_multiplied_feature_adjmat,
    const arma::mat& tau) {

  // This function computes pi for each combination of features. It yields pi only for D = 1

  // Prepare an Rcpp::List to which pis will be appended.
  int length_list_multiplied_feature_adjmat = list_multiplied_feature_adjmat.length();
  int length_list = length_list_multiplied_feature_adjmat / 2;
  Rcpp::List list_pi(length_list);

  // Calculate pi_d1x0
  arma::mat pi_d1x0 = compute_pi_d1x0(numOfVertices, numOfClasses, list_multiplied_feature_adjmat, tau, 0);
  list_pi[0] = pi_d1x0;

  // Compute pi for the other combinations of features
  int n_comb = list_multiplied_feature_adjmat.length() / 2;

  for (int s = 1; s < n_comb; s++) {
    // First, compute the conditional probability for D_ij = 0 and D_ij = 1 given X_ij = x
    // Create indices to extract appropriate multiplied matrices
    int index = 2 * s;
    // Get the multiplied matrix for D = 0
    arma::sp_mat D0X = list_multiplied_feature_adjmat[index];
    // Get the multiplied matrix for D = 1
    arma::sp_mat D1X = list_multiplied_feature_adjmat[index+1];
    // Compute the denominator for pi
    arma::sp_mat X = D0X + D1X;
    arma::mat pi_denom = tau.t() * X * tau;
    // Remove extremely small elements in the denominator
    #pragma omp parallel for
    for (int k = 0; k < numOfClasses; k++) {
      for (int l = 0; l < numOfClasses; l++) {
        if (pi_denom(k, l) < 1e-10) {
          pi_denom(k, l) = 1;
        }
      }
    }
    // Compute the numerator for D_ij = 1
    arma::mat pi_num = tau.t() * D1X * tau;
    // Compute pi for D_ij = 1
    arma::mat pi_d1 = pi_num/pi_denom;
    // Remove extremely small elements in pi_d1
    #pragma omp parallel for
    for (int k = 0; k < numOfClasses; k++) {
      for (int l = 0; l < numOfClasses; l++) {
        if (pi_d1(k, l) < minPi) {
          pi_d1(k, l) = minPi;
        }
      }
    }
    // Store pi
    list_pi[s] = pi_d1;
  }
    // Return the output
    return list_pi;
}

// A wrapper function to implement the EM algorithm using features
// [[Rcpp::export]]
Rcpp::List run_EM_with_features
  (int numOfVertices,
   int numOfClasses,
   const arma::vec& alpha,
   const Rcpp::List& list_multiplied_feature_adjmat,
   arma::mat& tau,
   int verbose = 0) {
  // Initialize lower bound
  double LB = 0;

  // Compute quadratic term with multiple features
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "runFixedPointEstimationEStepMM_sparse: compute quadratic term with features."<< "\n";
  }
  arma::mat A = compute_quadratic_term_with_features
    (numOfVertices,
     numOfClasses,
     list_multiplied_feature_adjmat,
     tau,
     LB,
     verbose);

  // Compute linear term
  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "runFixedPointEstimationEStepMM_sparse: compute linear."<< "\n";
  }
  arma::mat s = compute_linear_term(numOfVertices, numOfClasses, alpha, tau, LB);

  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "runFixedPointEstimationEStepMM_sparse: solveQP."<< "\n";
  }
  solveQP(A, s, tau, minTau, verbose);

  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "runFixedPointEstimationEStepMM_sparse: normalizeTau."<< "\n";
  }
  normalizeTau(tau, minTau);

  if (verbose >= 5) {
    auto time = std::chrono::system_clock::now();
    std::time_t timestamp = std::chrono::system_clock::to_time_t(time);
    Rcpp::Rcout << std::ctime(&timestamp) << "runFixedPointEstimationEStepMM_sparse: returning."<< "\n";
  }

  // Return \alpha and the lower bound as a list
  Rcpp::List output(2);
  output[0] = tau;
  output[1] = LB;
  return output;
}
