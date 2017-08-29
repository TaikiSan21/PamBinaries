
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_sdf.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_SDF_H_
#define IN_FRA_SDF_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "fra_type.h"

/* This file contains function declarations the mutils fractal library. */

#ifdef __cplusplus
extern "C" {
#endif

/** Nonparametric cross-spectral density function estimation.
 * Estimates the cross-spectral density function for a real-valued
 * uniformly sampled time series via a direct spectral
 * estimation scheme. Using a rectangular window with unit energy will yield
 * a periodogram.
 *
 * References:
 *
 * 1. D.B. Percival and A.T. Walden,
 * ``Spectral Analysis for Physical Applications: Multitaper and
 * Conventional Univariate Techniques'', Cambridge University Press,
 * 1993..
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_sdf.h
 * @source fra\_sdf.c
 * @library fractal
 * @usage #err = frauniv_spectral_density_function_direct( &time_series, &taper, recenter, single_sided, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series data, one time series per column,
 * @param taper       Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the taper. The number of elements must
 *                    be the same as that of the time series, but can
 *                    have transposed dimensions.
 * @param center      A logical flag. If TRUE, the mean is subtracted from the
 *                    time series prior to estimating the SDF.
 * @param recenter    A logical flag. If TRUE, the mean of tapered
 *                    time series is removed prior to estimating the SDF.
 * @param single_sided A logical flag. If TRUE, the single-sided SDF
 *   estimate is returned, i.e. over normalized frequencies
 *   f(k) = k / (npad) for k = 0, ..., floor(npad/2)/npad.
 *   Otherwise, the double-sided (symmetric) SDF is returned.
 * @param npad        The total length of each time series to analyze 
 *                    after padding with zeros (npad >= N).
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a single-column universal matrix of type
 *                    MUTIL\_DOUBLE with twice the number of elements as the
 *                    input time series and containing the (double-sided)
 *                    direct spectral density function estimate. The memory for this
 *                    matrix is automatically allocated by the function.
 * @see sig_taper
 * @see _fra_sdf_taper
 * @see fra_spectral_density_function_lag_window
 * @see fra_spectral_density_function_wosa
 * @see fra_spectral_density_function_multitaper
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_spectral_density_function_direct(
  const univ_mat *time_series,
  const univ_mat *taper,
  const boolean   center,
  const boolean   recenter,
  const boolean   single_sided,
  const sint32    npad,
  void           *intrp_ptr,
  univ_mat       *result );

/** Nonparametric cross-spectral density function estimation.
 * Estimates the cross-spectral density function for a real-valued
 * uniformly sampled time series via a lag window spectral
 * estimation scheme.
 *
 * References:
 *
 * 1. D.B. Percival and A.T. Walden,
 * ``Spectral Analysis for Physical Applications: Multitaper and
 * Conventional Univariate Techniques'', Cambridge University Press,
 * 1993..
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_sdf.h
 * @source fra\_sdf.c
 * @library fractal
 * @usage #err = frauniv_spectral_density_function_lag_window( &time_series, &lag_window, &taper, recenter, single_sided, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series data, one time series per column,
 * @param lag_window  Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the taper. The number of elements must
 *                    be the same as that of the time series, but can
 *                    have transposed dimensions. The window is expected
 *                    to contain values correpsonding to lags 0,...,N-1
 *                    where N is the number of samples in the original time
 *                    series. Recommended lag windows are the Parzen,
 *                    Paopoulis, and Daniell windows.
 * @param taper       Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the taper. The number of elements must
 *                    be the same as that of the time series, but can
 *                    have transposed dimensions.
 * @param recenter    A logical flag. If TRUE, the mean is subtracted from the
 *                    time series prior to estimating the SDF.
 * @param single_sided A logical flag. If TRUE, the single-sided SDF
 *   estimate is returned, i.e. over normalized frequencies
 *   f(k) = k / (npad) for k = 0, ..., floor(npad/2)/npad.
 *   Otherwise, the double-sided (symmetric) SDF is returned.
 * @param npad        The total length of each time series to analyze 
 *                    after padding with zeros (npad >= N).
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a single-column universal matrix of type
 *                    MUTIL\_DOUBLE with twice the number of elements as the
 *                    input time series and containing the (double-sided)
 *                    direct spectral density function estimate. The memory for this
 *                    matrix is automatically allocated by the function.
 * @see sig_taper
 * @see _fra_sdf_taper
 * @see fra_spectral_density_function_direct
 * @see fra_spectral_density_function_wosa
 * @see fra_spectral_density_function_multitaper
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_spectral_density_function_lag_window(
  const univ_mat *time_series,
  const univ_mat *lag_window,
  const univ_mat *taper,
  const boolean   center,
  const boolean   recenter,
  const boolean   single_sided,
  const sint32    npad,
  void           *intrp_ptr,
  univ_mat       *result );

/** Nonparametric cross-spectral density function estimation.
 * Estimates the cross-spectral density function for a real-valued
 * uniformly sampled time series via Welch's Overlapped
 * Segment Averaging (WOSA) spectral estimation scheme.
 * The segmentation of the time series is done in such a way
 * as to include all of the coefficients in the original time
 * series. Inasmuch, the number of overlapping points between
 * segments is adjusted accordingly. Use the
 *
 * References:
 *
 * 1. D.B. Percival and A.T. Walden,
 * ``Spectral Analysis for Physical Applications: Multitaper and
 * Conventional Univariate Techniques'', Cambridge University Press,
 * 1993..
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_sdf.h
 * @source fra\_sdf.c
 * @library fractal
 * @usage #err = frauniv_spectral_density_function_wosa( &time_series, &taper, overlap, recenter, single_sided, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series data, one time series per column,
 * @param taper       Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the taper. The number of elements in
 *                    this taper is used to define the width of the segmentation
 *                    window, Ns < N where N is the number of samples
 *                    in the original time series.
 * @param overlap     The suggested proportion of overlap (0 <= overlap < 1).
 *                    The actual number of overlap points is
 *                    No = Ns - floor( ( N - Ns ) / ( Nb - 1 ) ) + 1
 *                    where Nb = 1 + floor( ( N - Ns ) / ( Ns * ( 1 - overlap )))
 *                    is the number of blocks used to segment the time series.
 *                    This scheme esnures that all of the points in the
 *                    time series are used in the SDF estimation.
 * @param recenter    A logical flag. If TRUE, the mean is subtracted from the
 *                    time series prior to estimating the SDF.
 * @param single_sided A logical flag. If TRUE, the single-sided SDF
 *   estimate is returned, i.e. over normalized frequencies
 *   f(k) = k / (npad) for k = 0, ..., floor(npad/2)/npad.
 *   Otherwise, the double-sided (symmetric) SDF is returned.
 * @param npad        The total length of each time series to analyze 
 *                    after padding with zeros (npad >= N).
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a single-column universal matrix of type
 *                    MUTIL\_DOUBLE containing the (double-sided)
 *                    direct spectral density function estimate. The number
 *                    of elements in this matrix is 2*Ns. The memory for this
 *                    matrix is automatically allocated by the function.
 * @see sig_taper
 * @see _fra_sdf_taper
 * @see fra_spectral_density_function_direct
 * @see fra_spectral_density_function_lag_window
 * @see fra_spectral_density_function_multitaper
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_spectral_density_function_wosa(
  const univ_mat *time_series,
  const univ_mat *taper,
  const double    overlap,
  const boolean   center,
  const boolean   recenter,
  const boolean   single_sided,
  const sint32    npad,
  void           *intrp_ptr,
  univ_mat       *result );


/** Nonparametric cross-spectral density function estimation.
 * Estimates the cross-spectral density function for a real-valued
 * uniformly sampled time series via a multitaper
 * spectral estimation scheme.
 *
 * References:
 *
 * 1. D.B. Percival and A.T. Walden,
 * ``Spectral Analysis for Physical Applications: Multitaper and
 * Conventional Univariate Techniques'', Cambridge University Press,
 * 1993..
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_sdf.h
 * @source fra\_sdf.c
 * @library fractal
 * @usage #err = frauniv_spectral_density_function_multitaper( &time_series, &tapers, recenter, single_sided, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series data, one time series per column,
 * @param tapers      Pointer to a pre-allocated K x N
 *                    universal matrix of type MUTIL\_DOUBLE.
 *                    Each row of this matrix contains an N-point taper
 *                    where N is the number of samples in the original
 *                    time series. The number of tapers K, is limited
 *                    to K <= N. It is legal for K to be unity. The
 *                    tapers should be orthogonal to one another.
 *                    Seggested tapers are the sinusoidal tapers and
 *                    discrete prolate spheroidal sequences.
 * @param recenter    A logical flag. If TRUE, the mean is subtracted from the
 *                    time series prior to estimating the SDF.
 * @param single_sided A logical flag. If TRUE, the single-sided SDF
 *   estimate is returned, i.e. over normalized frequencies
 *   f(k) = k / (npad) for k = 0, ..., floor(npad/2)/npad.
 *   Otherwise, the double-sided (symmetric) SDF is returned.
 * @param npad        The total length of each time series to analyze 
 *                    after padding with zeros (npad >= N).
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a single-column universal matrix of type
 *                    MUTIL\_DOUBLE containing the (double-sided)
 *                    direct spectral density function estimate. The number
 *                    of elements in this matrix is 2*N. The memory for this
 *                    matrix is automatically allocated by the function.
 * @see sig_taper
 * @see _fra_sdf_taper
 * @see fra_spectral_density_function_direct
 * @see fra_spectral_density_function_lag_window
 * @see fra_spectral_density_function_multitaper
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_spectral_density_function_multitaper(
  const univ_mat *time_series,
  const univ_mat *taper,
  const boolean   center,
  const boolean   recenter,
  const boolean   single_sided,
  const sint32    npad,
  void           *intrp_ptr,
  univ_mat       *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_SDF_H_ */

