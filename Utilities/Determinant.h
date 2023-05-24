#ifndef DETERMINANT_HPP_INCLUDED
#define DETERMINANT_HPP_INCLUDED

/** \brief Computes determinant from given matrix using the best strategy for smaller matrices
 *
 * \param D double** matrix
 * \param order int Size of matrix.
 * \returns determinant.
 * \todo Check if log(det) works as well, in order to get smaller data
 *
 */
double determinant(double** D, int order);
#endif // DETERMINANT_HPP_INCLUDED
