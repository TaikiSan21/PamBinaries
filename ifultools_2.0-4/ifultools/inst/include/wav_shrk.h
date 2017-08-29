
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_shrk.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_WAV_SHRK_H_
#define IN_WAV_SHRK_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "wav_type.h"

/* This file contains function declarations the wavelet mutils library.
*/

#ifdef __cplusplus
extern "C" {
#endif

/** Nonlinear noise reduction in a time series via wavelet shrinkage.
 * Performs nonlinear noise reduction by defining a threshold
 * based on the estimated noise variance, shrinking (reducing
 * in amplitude toward zero) the discrete wavelet coefficients,
 * followed by an inverse discrete wavelet transform.
 *
 * References:
 *
 * 1. D. B. Percival and A. T. Walden, {\it Wavelet Methods for
 * Time Series Analysis}, Cambridge University Press, 2000.
 *
 * 2. Donoho, D. and Johnstone, I. (1992a). Ideal Spatial Adaptation
 * by Wavelet Shrinkage. Technical report, Department of Statistics,
 * Stanford University.
 *
 * 3. Donoho, D. and Johnstone, I. (1992b). Adapting to Unknown Smoothness
 * via Wavelet Shrinkage. Technical report, Department of Statistics,
 * Stanford University.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_shrk.h
 * @source wav\_shrk.c
 * @library wavelets
 * @usage #err = wavuniv_shrink( &time_series, &filters, &threshold, n_level, shrink_fun, threshold_type, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param  time_series Pointer to a pre-allocated universal matrix
 *    containing the time series to analyze.
 *    This input must a single-row or single-column
 *    matrix and can be of any type with the exception of MUTIL\_DCOMPLEX.
 * @param  filters Pointer to a matrix set containing two pre-allocated
 *    universal matrices of type MUTIL\_DOUBLE. The matrices
 *    must contain (in order) valid wavelet and scaling
 *    filters. The filter matrices must be a single-column
 *    or single-row and each must contain the same number of elements.
 * @param  threshold Pointer to a pre-allocated universal matrix of type 
 *    MUTIL\_DOUBLE containing the wavelet threshold(s).
 *    This matrix must be a single-column or single-row
 *    and can contain either a single element or the number
 *    of elements must be equal to the number of decomposition
 *    levels (n\_level). If the thresholds are to be estimated
 *    by a model instead (see arguments below), then set this
 *    argument to (univ_mat *) NULL (a NULL pointer).
 * @param threshold_function An enum of type \Ref{_wav_shrink_threshold}
 *    representing the threshold estimator. Choices
 *    are universal, minimax, and adaptive. Adaptive
 *    uses Stein's Unbiased Risk Estimator (SURE)
 *    for the threshold at each decomposition level
 *    and forces soft thresholding. This argument is only used
 *    if the threshold matrix is set to a NULL pointer.
 * @param threshold_scale A positive double value. The model-based
 *    threshold estimates are multiplied by this value to 
 *    either smooth or a roughen the resulting shrinkage result.
 *    This argument is only used if the threshold matrix is set 
 *    to a NULL pointer.
 * @param noise_variance A double value representing the estimated
 *    variance of the additive Gaussian noise. If this argument
 *    is zero or negative, then the noise variance is estimated 
 *    by the median absolute deviation of the scale 1 wavelet
 *    coefficients (normalized by 0.6745 for Gaussian distributions).
 *    This argument is only used if the threshold matrix is set 
 *    to a NULL pointer.
 * @param  shrink_function An enum of type \Ref{_wav_shrink_function} 
 *    representing the shrinkage function. Choices are hard, soft, and
 *    mid (hard/soft) thresholding.
 * @param  n_level The number of decomposition levels. This input should
 *    not exceed $2^{\lfloor \log_2{N}\rfloor}$.
 * @param  decimated  A boolean (logical) value. If TRUE, a decimated
 *    wavelet transform is performed. Otherwise, the
 *    maximal overlap DWT is used and the threshold
 *    values are automatically adjusted to accommodate
 *    larger sample sizes.
 * @param  intrp_ptr Pointer for implementation of interrupt checking.
 * @param  result Pointer to a universal matrix of type MUTIL\_DOUBLE
 *    which upon return will contain the result. The
 *    dimensions of this matrix are the same as that of the
 *    original time series. Memory for this matrix is
 *    allocated within the function.
 * @see _wav_shrink_function
 * @see _wav_shrink_threshold
 * @see wavuniv_shrink_threshold
 * @see wavuniv_transform_discrete_wavelet_convolution
 * @see wavuniv_transform_discrete_wavelet_convolution_inverse
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_transform_maximum_overlap_inverse
 * @see wavuniv_filters_daubechies
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_shrink(
  const univ_mat             *time_series,
  const mat_set              *filters,
  const univ_mat             *threshold,
  const wav_shrink_threshold  threshold_function,
  const double                threshold_scale,
  const double                noise_variance,
  const wav_shrink_function   shrink_function,
  const sint32                n_level,
  const boolean               decimated,
  void                       *intrp_ptr,
  univ_mat                   *result );

/** Waveshrink thresholds.
 * Calculates a waveshrink threshold for each specified
 * wavelet transform decomposition level.
 *
 * References:
 *
 * 1. D. B. Percival and A. T. Walden, {\it Wavelet Methods for
 * Time Series Analysis}, Cambridge University Press, 2000.
 *
 * 2. Donoho, D. and Johnstone, I. (1992a). Ideal Spatial Adaptation
 * by Wavelet Shrinkage. Technical report, Department of Statistics,
 * Stanford University.
 *
 * 3. Donoho, D. and Johnstone, I. (1992b). Adapting to Unknown Smoothness
 * via Wavelet Shrinkage. Technical report, Department of Statistics,
 * Stanford University.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_shrk.h
 * @source wav\_shrk.c
 * @library wavelets
 * @usage #err = wavuniv_shrink_threshold( &wavelet_transform, decimated, threshold_type, shrink_fun, noise_variance, intrp_ptr, &result);#
 * @return Standard mutils error/OK code.
 * @param  wavelet_transform Pointer to a pre-allocated matrix set
 *   containing either a decimated wavelet transform (DWT) or 
 *   a maximal overlap DWT (MODWT).
 * @param  decimated   A boolean (logical) value. If TRUE (FALSE), the 
 *   wavelet\_transform matrix set is assumed to contain a DWT (MODWT). 
 * @param  threshold_type An enum of type \Ref{_wav_shrink_threshold}
 *                     representing the threshold estimator. Choices
 *                     are universal, minimax, and adaptive. Adaptive
 *                     uses Stein's Unbiased Risk Estimator (SURE)
 *                     for the threshold at each decomposition level
 *                     and forces soft thresholding.
 * @param  shrink_fun  An enum of type \Ref{_wav_shrink_function} representing
 *                     the shrinkage function. Choices are hard, soft, and
 *                     mid (hard/soft) thresholding.
 * @param  noise_variance  A double value representing (an estimate of)
 *   the additive Gaussian white noise variance. If unknown, setting
 *   this value to 0.0 (or less) will prompt the function to automatically
 *   estimate the noise variance based on the median absolute deviation (MAD)
 *   of the scale one wavelet coefficients.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  result      Pointer to a single-row universal matrix of type 
 *   MUTIL\_DOUBLE containing the wavelet threshold(s). Memory for this
 *   matrix is automatically allocated within the function.
 *                     
 * @see _wav_shrink_function
 * @see _wav_shrink_threshold
 * @see wavuniv_transform_discrete_wavelet_convolution
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_shrink
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_shrink_threshold(
  const mat_set              *wavelet_transform,
  const boolean               decimated,
  const wav_shrink_threshold  threshold_type,
  const wav_shrink_function   shrink_fun,
  const double                noise_variance,
  void                       *intrp_ptr,
  univ_mat                   *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_SHRK_H_*/

