#include <iomanip>
#include "Matrix.h"

/** \brief Computes Vector- Square Matrix product
 *
 * \param x double* Vector (factor)
 * \param A double** Matrix (factor)
 * \param n int Size of matrix/vector.
 * \param b double* (reference) vector.
 * \return void
 * \todo Check if it can be optimised by using block products.
 *
 */
void vectorSqMatrixProduct(double* x,double** A,int n, double* b)
{
    for(int c=0;c<n;c++)
    {
        b[c]=0;
        for(int d=0;d<n;d++)
            b[c]+=x[d]*A[d][c];
    }
}

/** \brief Computes Vector dot product
 *
 * \param x double* Vector (factor)
 * \param y double* Vector (factor)
 * \param n int Size of vectors.
 * \return double Dot Product.
 *
 */
double vectorDotProduct(double* x, double* y, int n)
{
    double res=0.0;
    for(int c=0;c<n;c++)
        res+=x[c]*y[c];
    return res;
}
