#include <iomanip>
#include "Matrix.h"

/*
 * Computes Vector- Square Matrix product
 * x Vector (factor)
 * A: Matrix (factor)
 * n: Size of matrix/vector.
 * b: Output (reference) vector.
 * TODO: Check if it can be optimised by using block products.
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

/*
 * Computes Vector dot product
 * x Vector (factor)
 * y: Vector (factor)
 * n: Size of vectors.
 *
 * returns: Dot Product.
 */
double vectorDotProduct(double* x, double* y, int n)
{
    double res=0.0;
    for(int c=0;c<n;c++)
        res+=x[c]*y[c];
    return res;
}
