#ifndef GAUSSIANDIST_HPP_INCLUDED
#define GAUSSIANDIST_HPP_INCLUDED
#define MPI 3.14159265358979323846
#define ME 2.718281828459045235360

/*
 * Computes gaussian entropy from given covariance matrix.
 *  covM: Covariance matrix
 *  n: Size of correlation matrix.
 *
 *  returns: Gaussian entropy.
 */
double gaussianEntropy(double** covM, int n);

/*
 * Computes gaussian entropy from given covariance matrix.
 *  covM: Covariance matrix
 *  n: Size of correlation matrix.
 *
 *  returns: Gaussian entropy.
 */
double gaussianMutualInformation(double** corrM, int n);

#endif // GAUSSIANDIST_HPP_INCLUDED
