
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_lyap.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_LYAP_H_
#define IN_FRA_LYAP_H_


#include "fra_type.h"



/* This file contains function declarations for routines related */
/* to computing Lyapunov exponents.                              */

#ifdef __cplusplus
extern "C" {
#endif


/** Lyapunov scaling curves.
 * This function returns the Lyapunov scaling function over a
 * specified range for a given set of embedding dimensions and neighborhood
 * distances (epsilons). The maximal Lyapunov exponent may then be estimated
 * by the average slope of the curves.
 *
 * For each given dimension and epsilon pair (d,eps), this function creates
 * a d-dimensional delay embedding of the given time series data and attempts
 * to find {\bf nref_points} reference points among the resulting delay
 * vectors from which to compute the Lyapunov scaling function.
 *
 * The output of the function is a matrix set which contains two matrices.
 * The first matrix has dimensions ({\bf nfuture}+1)$\times DE$, where $D$
 * is the length of input {\bf dims} and $E$ is the length of {\bf epsilons}.
 * Each column holds a Lyapunov scaling function curve for
 * a specific (d,eps), with eps running faster than d. The (d,eps) pairs are
 * taken in order from the {\bf dims} and {\bf epsilons} input vectors. The
 * second matrix of the output matrix set contains a column vector holding
 * the values 0,..nfuture. This may be used as the absicca in subsequent
 * plotting of the Lyapunov scaling function curves.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_lyap.h
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = frauniv_lyapunov_maximal(&ts,dly,&dims,&eps,FRA_DISTANCE_L1,nref,mneig,olag,nfut,iptr,&scl);#
 * @return Standard mutils error/OK code.
 * @param time_series   Pointer to a universal matrix of type MUTIL_DOUBLE
 *                      containing the scalar time series to be embedded. May
 *                      be a row or column vector.
 * @param delay         Positive integer specifying the delay to use when
 *                      embedding the given time series.
 * @param dims          Pointer to a universal matrix of type MUTIL_SINT32
 *                      containing the desired dimensions of the time
 *                      series embeddings. May be a row or column vector with
 *                      all elements greater than 1.
 * @param epsilons      Pointer to a universal matrix of type MUTIL_DOUBLE
 *                      containing the maximum distances from a reference
 *                      point, with respect to the given metric, that
 *                      neighbors may be found.
 * @param metric        Specifies the distance measure used when searching
 *                      neighbors. See \Ref{fra_distance_metric} for possible
 *                      values.
 * @param nref_points   Positive integer specifiying the number of reference
 *                      points to find for each delay embedding. This is a
 *                      request and gives an upper bound on the number of
 *                      reference points to find.
 * @param min_nneig     Positive integer specifying the minimim number of
 *                      neighbors that each reference point must have.
 * @param orbital_lag   Positive integer specifying the minimum temporal
 *                      distance between reference points. No two reference
 *                      points will be less than this number of lags apart.
 * @param nfuture       Number of future time steps in which to compute the
 *                      scaling function.
 * @param intrp_ptr     Pointer for implementing user interrupts.
 * @param result        A matrix set containing two universal matrices of
 *                      type MUTIL_DOUBLE. See descriptions of the universal
 *                      matrices above.
 *
 * @see frauniv_lyapunov_scaling_function
 * @see _fra_distance_metric
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_lyapunov_maximal(
  const univ_mat             *time_series,
  const sint32                delay,
  const univ_mat             *dims,
  const univ_mat             *epsilons,
  const fra_distance_metric   metric,
  const sint32                nref_points,
  const sint32                min_nneig,
  const sint32                orbital_lag,
  const sint32                nfuture,
  void                       *intrp_ptr,
  mat_set                    *result );




/** Lyapunov scaling function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_lyap.h
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = frauniv_lyapunov_scaling_function(&ts,dly,dim,eps,FRA_DISTANCE_L2,nref,mneig,olag,nfut,iptr,&ridx,&rnneig,&scl);#
 * @return Standard mutils error/OK code.
 * @param time_series   Pointer to a universal matrix of type MUTIL_DOUBLE
 *                      containing the scalar time series to be embedded. May
 *                      be a row or column vector.
 * @param delay         Positive integer specifying the delay to use when
 *                      embedding the given time series.
 * @param dim           A positive integer specifying the desired dimension
 *                      of the time series embedding. Must be greater than 1.
 * @param epsilon       A positive number specifying the radius, with respect
 *                      to the given metric, of the neighborhood in which to
 *                      find reference point neighbors.
 * @param metric        Specifies the distance measure used when searching
 *                      neighbors. See \Ref{fra_distance_metric} for possible
 *                      values.
 * @param nref_points   Positive integer specifiying the desired number of
 *                      reference points to be found. This is a request and
 *                      provieds an upper bound.
 * @param min_nneig     Positive integer specifying the minimim number of
 *                      neighbors that each reference point must have.
 * @param orbital_lag   Positive integer specifying the minimum temporal
 *                      distance between reference points. No two reference
 *                      points will be less than this number of lags apart.
 * @param nfuture       Number of future time steps in which to compute the
 *                      scaling function.
 * @param intrp_ptr     Pointer used for user interrupts.
 * @param nref_found    The number of reference points actually found by the
 *                      function. Will be less than or equal to
 *                      {\bf nref_points}.
 * @param ref_indices   A universal matrix that holds the indices to the
 *                      found reference points. These are zero-based
 *                      indices with respect to the embedding matrix created
 *                      by the function.
 * @param ref_nneig     A universal matrix the same length as {\bf ref_indices}
 *                      which contains the number of neighbors found for
 *                      each reference point.
 * @param scaling       A universal matrix of type MUTIL_DOUBLE and length
 *                      {\bf nfuture}+1 that contains the Lyapunov scaling
 *                      function.
 *
 * @see _fra_distance_metric
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_lyapunov_scaling_function(
  const univ_mat             *time_series,
  const sint32                delay,
  const sint32                dim,
  const double                epsilon,
  const fra_distance_metric   metric,
  const sint32                nref_points,
  const sint32                min_nneig,
  const sint32                orbital_lag,
  const sint32                nfuture,
  void                       *intrp_ptr,
  sint32                     *nref_found,
  univ_mat                   *ref_indices,
  univ_mat                   *ref_nneig,
  univ_mat                   *scaling );


/** Local Lyapunov spectrum estimation.
 * Estimates the local Lyapunov exponents over a range
 * of user supplied scales and dimensions. The local Lyapunov
 * spectrum is calculated as follows: (1) A delayed embedding
 * of the input time series is formed by embedding the series
 * in dimension embedding\_dimension using a time delay of
 * time\_lag. (2) For each global reference point (specified
 * by an intger index in the reference matrix) a local Lyapunov
 * spectrum is calculated, one exponent for each dimension
 * from 1 to local\_dimension and for each (integer) scale
 * specified in the scale matrix. As the scales grow larger,
 * the Lyapunov exponent estimates tend toward asymptotic
 * values corresponding to the global Lyapunov exponents. The
 * details of how each local spectrum is estimated is given below.
 * (3) The local spectra are then averaged over each global reference
 * point to stabilize the results.
 *
 * Each local spectrum is obtained by estimating the eigenvalues
 * of the so-called Oseledec matrix, which is formed through
 * a matrix product of successive local Jacobians with the transpose
 * of the Jacobians. The number of Jacobians in the product is equivalent
 * to the scale. Each Jacobian is formed by fitting a local neighborhood
 * of points (relative to a some reference point) with a multidimensional
 * polynomial of order polynomial\_order. The number of neighbors found
 * for each reference point in the embedding is chosen to be twice the
 * polynomial order for numerical stability. To further stabilize the results,
 * a local Lyapunov spectrum is formed for each local reference point of which
 * there are n\_reference\_local total local reference points.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_lyap.c
 * @library fractal
 * @usage #  err = frauniv_local_lyapunov_spectrum(
time_series,embedding_dimension,time_lag,orbital_lag,sampling_interval,local_dimensionl,polynomial_order,n_reference_local,scale,intrp_ptr,&result );#
 * @param time_series Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series.
 * @param embedding_dimension An integer representing the maximal embedding dimension.
 *   This dimension is used to embed the data and for finding nearest neighbors.
 * @param time_lag  An integer representing the delay between coordinates.
 * @param orbital_lag The number of points along the trajectory of the
 *                    current point that must be exceeded in order for
 *                    another point in the phase space to be considered
 *                    a neighbor candidate. This argument is used
 *                    to help attenuate temporal correlation in the
 *                    the embedding which can lead to spuriously low
 *                    information dimension estimates. The orbital lag
 *                    must be positive or zero.
 * @param sampling_interval The sampling interval of the time series.
 * @param local_dimension An integer representing the dimension (number of)
 *   local Lyapunov exponents to estimate. This value msut be less than
 *   or equal to the embedding dimension.
 * @param polynomial_order The order of the polynomial to use in fitting data
 *                    around reference points in the phase space. This poloynomial
 *                    fit will be used to form the Jacobians which are
 *                    in  turn used to calcualte the Lypaunov exponents.
 * @param global_reference   A pointer to a single-column or single row universal
 *  matrix of type MUTIL\_SINT32 containing a set of global reference point indices
 *  to use in estimating the local Lyapunov spectrum. A local spectrum is estimated
 *  around each global reference point, and all the local spectra are then
 *  averaged to stabilize the results. These global reference points should be
 *  chosen such that they are far apart in time.
 * @param n_reference_local The number of neighbors to use in
 *                    in developing the kd-tree (used as a quick means
 *                    of finding nearest neighbors in the phase space).
 *                    These neighbors are colelcted relative to the reference
 *                    points.
 * @param metric      The metric used to define the distance between
 *                    points in the embedding. The argument is of
 *                    type \Ref{_fra_distance_metric}.
 * @param scale       A pointer to a single-column or single-row
 *   universal matrix of type MUTIL\_SINT32. Each element of this vector
 *   represents the number of consecutive points along a trajectory
 *   over which the local Lyapunov exponents are averaged. As this scale
 *   increases, one expects the local Lyapunov exponent estimates to converge
 *   towards the global estimates. All scales must be great than one.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a matrix set which (upon return) will
 *                    contain local\_dimension single-column univeral
 *                    matrices of type MUTIL\_DOUBLE. The length of each
 *                    vector corresponds to the number of scales analyzed
 *                    (number of elements in the scale matrix). Each vector
 *                    thus contains the local Lyapunov exponent estimates
 *                    across the supplied scales for a given dimension.
 *                    The memory for the
 *                    matrix set is automatically allocated within the
 *                    function and will containing two universal
 *                    matrices of type MUTIL\_DOUBLE.
 * @see frauniv_embed
 * @see frauniv_dimension_correlation_summation
 * @see frauniv_dimension_information
 */

MUTIL_LIBEXPORT mutil_errcode frauniv_local_lyapunov_spectrum(
  const univ_mat *time_series,
  const sint32    embedding_dimension,
  const sint32    time_lag,
  const sint32    orbital_lag,
  const double    sampling_interval,
  const sint32    local_dimension,
  const sint32    polynomial_order,
  const univ_mat *global_reference,
  const sint32    n_reference_local,
  const fra_distance_metric metric,
  const univ_mat *scale,
  void           *intrp_ptr,
  mat_set        *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_LYAP_H_ */
