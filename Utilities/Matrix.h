#ifndef MATRIX_HPP_INCLUDED
#define MATRIX_HPP_INCLUDED

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
void vectorSqMatrixProduct(double* x,double** A,int n, double* b);

/** \brief Computes Vector dot product
 *
 * \param x double* Vector (factor)
 * \param y double* Vector (factor)
 * \param n int Size of vectors.
 * \return double Dot Product.
 *
 */
double vectorDotProduct(double* x, double* y, int n);

#endif // MATRIX_HPP_INCLUDED
