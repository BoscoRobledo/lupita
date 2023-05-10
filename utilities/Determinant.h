#ifndef DETERMINANT_HPP_INCLUDED
#define DETERMINANT_HPP_INCLUDED

/*
 * Computes determinant from given matrix using the best strategy for smaller matrices
 *  D: matrix
 *  n: Size of matrix.
 *
 *  returns: determinant.
 */
double determinant(double** D, int order);
#endif // DETERMINANT_HPP_INCLUDED
