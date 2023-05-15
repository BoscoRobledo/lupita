#include <cmath>
#include "Determinant.h"
#include "GaussianDistribution.h"

/*
 * Computes gaussian entropy from given covariance matrix.
 *  covM: Covariance matrix
 *  n: Size of correlation matrix.
 *
 *  returns: Gaussian entropy.
 */
double entropy(double** covM, int n)
{
    if(n==1)
        return 0.5*log(2.0*MPI*ME*covM[0][0]);
    else
        return (double)n/2.0*log(2.0*MPI*ME)+0.5*log(determinant(covM,n));
}

/*
 * Computes gaussian mutual information from given correlation matrix.
 *  corrM: Correlation matrix
 *  n: Size of correlation matrix.
 *
 *  returns: Mutual information.
 */
double mutualInformation(double** corrM, int n)
{
    if(n==1)
        return 0.0;
    return -0.5*log(determinant(corrM,n));
}
