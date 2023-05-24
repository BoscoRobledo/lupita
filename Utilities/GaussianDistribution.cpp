#define _USE_MATH_DEFINES

#include <cmath>
#include "Determinant.h"
#include "GaussianDistribution.h"

/** \brief Computes gaussian entropy from given covariance matrix.
 *
 * \param covM double** Covariance matrix
 * \param n int Size of matrix.
 * \returns Gaussian Entropy.
 *
 */
double gaussianEntropy(double** covM, int n)
{
    if(n==1)
        return 0.5*log(2.0*M_PI*M_E*covM[0][0]);
    else
        return (double)n/2.0*log(2.0*M_PI*M_E)+0.5*log(determinant(covM,n));
}


/** \brief Computes gaussian mutual information from given correlation matrix.
 *
 * \param covM double** Correlation matrix
 * \param n int Size of matrix.
 * \returns Mutual information.
 *
 */
double gaussianMutualInformation(double** corrM, int n)
{
    if(n==1)
        return 0.0;
    return -0.5*log(determinant(corrM,n));
}
