#include <cmath>
#include "Determinant.hpp"
#include "GaussianDist.hpp"
double entropy(double** covM, int n)
{
    if(n==1)
        return 0.5*log(2.0*MPI*ME*covM[0][0]);
    else
        return (double)n/2.0*log(2.0*MPI*ME)+0.5*log(determinant(covM,n));
}

double mutualInfo(double** corrM, int n)
{
    if(n==1)
        return 0.0;
    return -0.5*log(determinant(corrM,n));
}
