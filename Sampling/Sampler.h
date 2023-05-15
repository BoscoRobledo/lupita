#ifndef SIMULATE_HPP_INCLUDED
#define SIMULATE_HPP_INCLUDED
#include <iostream>
#include <random>
#include <chrono>
#include <Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

namespace Eigen {
namespace internal {
template<typename Scalar>
struct scalar_normal_dist_op
{
  static boost::mt19937 rng;    // The uniform pseudo-random algorithm
  mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

  EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

  template<typename Index>
  inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
};

template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;

template<typename Scalar>
struct functor_traits<scalar_normal_dist_op<Scalar> >
{ enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };
} // end namespace internal
} // end namespace Eigen

void getSamplesMarginal(double* mu, double* var, double ** sample, int D, int n);
void getSamplesChowLiu(double* mu, double** covM, double** corrM, double ** sample, int D, int n);
/** \brief Get samples from a Gaussian Network built from a t-cherry junction tree.
 *
 * \param mu double* Unconditional means vector
 * \param covM double** Covariance Matrix to build the Cherry tree from
 * \param corrM double** Correlation Matrix to build the Cherry tree from
 * \param sample double** Memory to fill with samples
 * \param D int Sample dimensionality
 * \param n int Sample size
 * \param k int Desired order of the Cherry tree
 * \return double Weight of the Cherry Tree modeled from population
 *
 */
double getSamplesTChJT(double* mu, double** covM, double** corrM, double ** sample, int D, int n, int k, int evals, int func);
void getSamplesMVN(double* mu, double** covM,  double ** sample, int D, int n);

#endif // SIMULATE_HPP_INCLUDED
