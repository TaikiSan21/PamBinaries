
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_surr.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_SURR_H_
#define IN_FRA_SURR_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "fra_type.h"

/* This file contains function declarations the mutils fractal library. */

#ifdef __cplusplus
extern "C" {
#endif

/** Generate surrogate time series via Theiler's methods.
 * Create surrogates of a time series using Fourier based techniques.
 * Choices are (i) phase randomization, (ii) Theiler's amplitude
 * adjusted Fourier transform (AAFT).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_surr.h
 * @source fra\_surr.c
 * @library fractal
 * @usage #err = mutil_errcode frauniv_bootstrap_theiler( &time_series, method, seed, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series.
 * @param method      An enum of type \Ref{_fra_surrogate} denoting the
 *                    method to use in creating the surrogate data.
 * @param seed Unsigned long integer representing the initial
 *             random seed value. If zero, the seed is reset
 *             based on the current time.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a single-column universal matrix of type
 *                    MUTIL\_DOUBLE and of the same length of the original
 *                    time series. Upon return, this matrix will contain
 *                    the surrogate realization. The memory for this
 *                    matrix is automatically allocated within the function.
 *
 * @see frauniv_bootstrap_davison_hinkley
 * @see frauniv_bootstrap_circulant_embedding
 * @see _fra_surrogate
 * @see frauniv_determinism_delta_epsilon
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_bootstrap_theiler(
  const univ_mat      *time_series,
  const fra_surrogate  method,
  const uint32         seed,
  void                *intrp_ptr,
  univ_mat            *result );

/** Generate surrogate time series via Davison-Hinkley method.
 * Create surrogates of a time series using Davison-Hinkley
 * method which randomizes both the phases and the amplitudes
 * of original series.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_surr.h
 * @source fra\_surr.c
 * @library fractal
 * @usage #err = mutil_errcode frauniv_bootstrap_davison_hinkley( &time_series, seed, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series.
 * @param seed Unsigned long integer representing the initial
 *             random seed value. If zero, the seed is reset
 *             based on the current time.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a single-column universal matrix of type
 *                    MUTIL\_DOUBLE and of the same length of the original
 *                    time series. Upon return, this matrix will contain
 *                    the surrogate realization. The memory for this
 *                    matrix is automatically allocated within the function.
 *
 * @see frauniv_bootstrap_theiler
 * @see frauniv_bootstrap_circulant_embedding
 * @see _fra_surrogate
 * @see frauniv_determinism_delta_epsilon
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_bootstrap_davison_hinkley(
  const univ_mat      *time_series,
  const uint32         seed,
  void                *intrp_ptr,
  univ_mat            *result );

/** Generate surrogate time series via circulant embedding.
 * Create surrogates of a time series using a single-sided
 * spectral density function (estimate) in a circulant embedding
 * scheme.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_surr.h
 * @source fra\_surr.c
 * @library fractal
 * @usage #err = mutil_errcode frauniv_bootstrap_circulant_embedding( &sdf, seed, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param sdf         Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the single-sided spectral density function
 *                    estimate for normalized frequencies f(k) = k/(2N).
 *                    for k = 0, ..., N where N is the number of samples in the
 *                    original time series. Thus, this is an N + 1 element
 *                    vector.
 * @param seed Unsigned long integer representing the initial
 *             random seed value. If zero, the seed is reset
 *             based on the current time.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a single-column universal matrix of type
 *                    MUTIL\_DOUBLE and of the same length of the original
 *                    time series. Upon return, this matrix will contain
 *                    the surrogate realization. The memory for this
 *                    matrix is automatically allocated within the function.
 *
 * @see frauniv_bootstrap_theiler
 * @see frauniv_bootstrap_davison_hinkley
 * @see _fra_surrogate
 * @see frauniv_spectral_density_function_direct
 * @see frauniv_spectral_density_function_lag_window
 * @see frauniv_spectral_density_function_wosa
 * @see frauniv_spectral_density_function_multitaper
 * @see frauniv_determinism_delta_epsilon
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_bootstrap_circulant_embedding(
  const univ_mat      *sdf,
  const uint32         seed,
  void                *intrp_ptr,
  univ_mat            *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_SURR_H_ */




