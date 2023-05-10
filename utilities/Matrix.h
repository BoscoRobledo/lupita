#ifndef MATRIX_HPP_INCLUDED
#define MATRIX_HPP_INCLUDED

/*
 * Computes Vector- Square Matrix product
 * x Vector (factor)
 * A: Matrix (factor)
 * n: Size of matrix/vector.
 * b: Output (reference) vector.
 *
 */
void vectorSqMatrixProduct(double* x, double** A, int n, double* b);

/*
 * Computes Vector dot product
 * x Vector (factor)
 * y: Vector (factor)
 * n: Size of vectors.
 *
 * returns: Dot Product.
 */
double vectorDotProduct(double* x, double* y, int n);

#endif // MATRIX_HPP_INCLUDED
