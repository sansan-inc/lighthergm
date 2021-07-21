#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute Yule's Φ-coefficient
//' @param z_star a true block membership
//' @param z an estimated block membership
//' @export
// [[Rcpp::export]]
double compute_yule_coef(
    const arma::vec& z_star,
    const arma::vec& z) {

  double n_00 = 0;
  double n_01 = 0;
  double n_10 = 0;
  double n_11 = 0;

  int numOfVertices = z.n_elem;

  // Remove missing values from z_star and z
  arma::vec z_star_new = z_star;
  arma::vec z_new = z;
  for (int i = 0; i < numOfVertices; i++) {
    if (std::isnan(z_star_new(i))) {
      z_star_new[i] = -100;
    }
    if (std::isnan(z_new(i))) {
      z_new[i] = -1;
    }
  }

  for (int i = 0; i < numOfVertices; i++) {
    for (int j = i+1; j < numOfVertices; j++) {
      if (z_star_new(i) == z_star_new(j)) {
        if (z_new(i) == z_new(j)) {
          n_11 += 1;
        } else {
          n_10 += 1;
        }
      } else {
        if (z_new(i) == z_new(j)) {
          n_01 += 1;
        } else {
          n_00 += 1;
        }
      }
    }
  }

  // Compute Yule's φ-coefficnent
  double num = (n_00 * n_11) - (n_01 * n_10);
  double denom = sqrt((n_00 + n_01) * (n_10 + n_11) * (n_00 + n_10) * (n_01 + n_11));
  double phi = num / denom;

  // Return the output
  return phi;
}
