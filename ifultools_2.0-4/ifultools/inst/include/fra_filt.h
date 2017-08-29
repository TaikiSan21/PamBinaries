
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_filt.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_FILT_H_
#define IN_FRA_FILT_H_

#include "ut_plat.h"
#include "mat_type.h"
#include "fra_type.h"

/* This file contains function declarations for
   nonlinear denoising.
   the MUTILS fractal library. */

#ifdef __cplusplus
extern "C" {
#endif

/** Median filtering of a time series.
 *
 * Given a filter order P, the output of applying a median filter to a
 * time series X[t] for t = 0, ..., N - 1 is
 *
 * P odd: Y[ k ] = median( X[ k - ( P - 1 ) / 2, ...,   k + ( P - 1 ) / 2 ] )
 *
 * P even: Y[ k ] = median( X[ k - P / 2, ...,   k + P / 2 - 1 ] )
 *
 * for k = 0, ..., N - 1.
 * Thus, median filtering replaces the k-th
 * value of the time series with the median of the time series
 * over an P-point window centered about point k.
 * In the case where a portion of the window exceeds the boundaries
 * of the time series, the values outside the boundaries are ignored in
 * the median value calculation.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_filt.h
 * @source fra\_filt.c
 * @library fractal
 * @usage #err = frauniv_filter_median( &time_series, order, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of any type except MUTIL\_DCOMPLEX
 *                    containing the time series.
 * @param order       An integer representing the filter order.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a universal matrix of type MUTIL\_DOUBLE
 *                    which (upon return) will contain the median filtered
 *                    time series.The memory for this matrix is automatically
 *                    allocated within the function. The size of the returned
 *                    matrix is equivalent to that of the input time series.
 * @see frauniv_filter_local_projection
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_filter_median(
  const univ_mat *time_series,
  const sint32    order,
  void           *intrp_ptr,
  univ_mat       *result );



/** Nonlinear denoising via local projection.
 *
 * Performs one interation of nonlinear denoising of a given time series. The
 * technique is intended for denoising deterministic (chaotic) systems
 * contaminated with additive measurement noise as opposed to system/dynamic
 * noise. The algorithm works by creating a delay embedding of the give time
 * series and attenuating the noise in every point of the embedding by (1)
 * finding the nearest neighbors, (2) linearizing the flow in the
 * neighborhood, (3) estimating the noise component of the linearized flow by
 * projecting it onto a subspace defined by the weakest eigenvalues and
 * eigenvectors of the covariance matrix for the linearized flow neighborhood,
 * and (4) subtracting the noise component from the original point. In
 * general, this technique is effective if the embedding dimension is much
 * higher than the manifold containing the dynamics. The dimension of the
 * noise space must be estimated a priori and used as an input argument to
 * this function.
 *
 * References:
 *
 * 1. Peter Grassberger, Rainer Hegger, Holger Kantz, and Carsten Schaffrath,
 * ``On noise reduction methods for chaotic data'', Chaos 3(2), pp. 127-141,
 * 1993.
 *
 * 2. Holger Kantz and Thomas Schreiber, ``Nonlinear Time Series Analysis'',
 * Cambridge University Press, 1997.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_filt.h
 * @source fra\_filt.c
 * @library fractal
 * @usage #err = frauniv_filter_nonlinear_local_projection(&ts,edim,nneig,eps,metric,ndim,TRUE,intrp_ptr,&result );#
 * @return Standard mutils error/OK code.
 * @param time_series     Pointer to a pre-allocated row or column universal
 *                        matrix of type MUTIL\_DOUBLE containing the time
 *                        series to be denoised.
 * @param dim             The desired delay embedding dimension.
 * @param delay           The delay to used when embedding the input time
 *                        series.
 * @param min_nneig       The minimum number of neighbors to include
 *                        in the set of neighbor candidates for a given
 *                        point in the embedding. This input must be greater
 *                        than dim.
 * @param max_radius      Used to obtain the subset of neighbor candidates
 *                        to form the so-called admissible neighbors.
 *                        By definition, an admissible neighbor is one
 *                        which is separated from the original point
 *                        by a distance less than max\_radius.
 *                        The metric used to measure the separation is defined
 *                        by the distance\_metric input argument.
 * @param distance_metric The metric used to define the distance between
 *                        points in the embedding. The argument is of
 *                        type \Ref{_fra_distance_metric}.
 * @param noise_dim       The estimated dimension of the noise. This value
 *                        obviously must be less than input dim.
 * @param correct_curvature A logical flag. If TRUE, the center of mass
 *                        estimates, used in linearizing the local flow, are
 *                        corrected for local curvature effects. This results
 *                        in a higher computational burden with the benefit of
 *                        providing a more accurate portryal of the local
 *                        dynamics and thus a more accurate noise reduction.
 * @param intrp_ptr       Pointer for implementation of interrupt checking.
 * @param result          Pointer to a universal matrix of type MUTIL\_DOUBLE
 *                        which, upon return, will contain the denoised time
 *                        series. The memory for this matrix is automatically
 *                        allocated within the function. The size of the
 *                        returned matrix is equivalent to that of input
 *                        time\_series.
 *
 * @see _fra_distance_metric
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_filter_nonlinear_local_projection(
  const univ_mat            *time_series,
  const sint32               dim,
  const sint32               delay,
  const sint32               min_nneig,
  const double               max_radius,
  const fra_distance_metric  distance_metric,
  const sint32               noise_dim,
  const boolean              correct_curvature,
  void                      *intrp_ptr,
  univ_mat                  *result );


#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_FILT_H_ */
