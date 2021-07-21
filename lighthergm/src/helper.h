#ifndef __HELPER__
#define __HELPER__

arma::vec decimal_to_binary_vector(
    int decimal,
    int vec_length);

arma::mat compute_sumTaus(
    int numOfVertices,
    int numOfClasses,
    const arma::mat& tau,
    int verbose = 0);

void normalizeTau(
    arma::mat& tau,
    double minValue);

arma::rowvec sumDoubleMatrixByRow(
    const arma::mat& matrix);

#endif // __HELPER__
