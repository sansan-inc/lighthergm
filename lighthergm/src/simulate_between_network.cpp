#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_max_threads() 0
#endif
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// Function that simulates a between-block network.
// The first element of `coef_between` must be the edges parameter.
// [[Rcpp::export]]
arma::sp_mat simulate_between_network
  (int numOfVertices,
   const Rcpp::List& list_feature_adjmat,
   const arma::vec& coef_between,
   const arma::vec& block_membership,
   bool directed
   )
{
  // Number of covariates
  int numOfCovariates = list_feature_adjmat.length();
  // Initialize a sparse adjacency matrix for the between-block network
  arma::sp_mat between_adjmat = arma::sp_mat(numOfVertices, numOfVertices);
  // Prepare a sparse adjacency cube
  arma::field<arma::sp_mat> feature_cube(numOfCovariates);
  for (int p = 0; p < numOfCovariates; p++) {
    feature_cube(p) = Rcpp::as<arma::sp_mat>(list_feature_adjmat[p]);
  }
  // Necessary for R random number generator
  GetRNGstate();

  #pragma omp parallel
  {
    // Simulate between-block links
    #pragma omp for
    for (int j = 0; j < numOfVertices; j++) {
      for (int i = 0; i < numOfVertices; i++) {
        // Skip as many unnecessary calculations as possible in this nested loop, which makes the computation faster.
        if (block_membership[i] != block_membership[j] && ((directed && i != j) || (!directed && i < j))) {
          double x = unif_rand();
          double u = coef_between[0];
          for (int p = 0; p < numOfCovariates; p++) {
            double elem = feature_cube(p)(i, j);
            double elem_coef = coef_between[p+1];
            u += elem_coef * elem;
          }
          //std::printf("x: %f, Thread: %d, Loop: (%d, %d), u: %f\n", x, omp_get_thread_num(), i, j, u);
          if (u > log(x/(1-x))) {
            between_adjmat(i, j) = 1;
            }
          }
        }
      }
    }
  // This must be called after GetRNGstate before returning to R.
  PutRNGstate();
  return between_adjmat;
}
