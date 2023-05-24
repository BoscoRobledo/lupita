#ifndef GAUSSIANDISTRIBUTION_HPP_INCLUDED
#define GAUSSIANDISTRIBUTION_HPP_INCLUDED

/** \brief Computes gaussian entropy from given covariance matrix.
 *
 * \param covM double** Covariance matrix
 * \param n int Size of matrix.
 * \returns Gaussian Entropy.
 *
 */
double gaussianEntropy(double** covM, int n);

/** \brief Computes gaussian mutual information from given correlation matrix.
 *
 * \param covM double** Correlation matrix
 * \param n int Size of matrix.
 * \returns Mutual information.
 *
 */
double gaussianMutualInformation(double** corrM, int n);

#endif // GAUSSIANDISTRIBUTION_HPP_INCLUDED
