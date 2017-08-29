
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_cwt.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_WAV_CWT_H_
#define IN_WAV_CWT_H_

#include "mat_type.h"
#include "ut_err.h"
#include "ut_plat.h"
#include "ut_type.h"
#include "wav_type.h"

/* This file contains function declarations the wavelet mutils library. */

#ifdef __cplusplus
extern "C" {
#endif

/** The continuous wavelet transform.
 * The continuous wavelet transform (CWT) is a highly redundant
 * transformation of a real-valued or complex-valued function $f(x)$, mapping it
 * from the time domain to the so-called time-scale domain. Loosely,
 * speaking the CWT coefficients are proportional to the variability of a
 * function at a given time and scale.
 *
 * The CWT is defined by a
 * complex correlation of a scaled and time-shifted mother wavelet with a
 * function $f(x)$.
 * Let $\psi(x)$ be a real- or complex-valued function representing a
 * \textit{mother} wavelet, i.e. a function which meets the standard
 * mathemtical criteria for a wavelet (see \Ref{wavuniv_filters_continuous})
 * and one that can be used to generate
 * all other wavelets within the same family.
 * Let $\psi^*(\cdot)$ be the complex conjugate of $\psi(\cdot)$. The CWT of
 * $f(x)$ is defined as
 *
 * \[
 * W_f(a,b) \equiv \frac{1}{\sqrt{a}} \int_{-\infty}^\infty f(x) \psi^* \Bigl( \frac{x-b}{a}
 * \Bigr) \; dx \qquad \mbox{for } (a,b) \in \mathbf{R} \mbox{ and } a > 0
 * \]
 *
 * where $a$ is the scale of the wavelet and $b$ is the shift of the
 * wavelet in time. It can be shown the that the above complex correlation
 * maintains a duality with the Fourier transform defined by the relation
 *
 * \[ W_f(a,b) \equiv \frac{1}{\sqrt{a}} \int_{-\infty}^\infty f(x) \psi^* \Bigl( \frac{x-b}{a}
 * \Bigr) \; dx \longleftrightarrow \sqrt{a} \, F(\omega) \,\Psi^*(a\omega) \]
 *
 * where $F(\cdot)$ is the Fourier transform of $f(x)$ and $\omega$ is the frequency in radians.
 * This function calculates the CWT in the Fourier domain followed by an inverse Fourier transform.
 *
 * References:
 * D. B. Percival and  A. T. Walden, {\it Wavelet Methods for
 * Time Series Analysis}, Cambridge University Press, 2000.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_cwt.h
 * @source wav\_cwt.c
 * @library wavelets
 * @usage #err = wavuniv_transform_continuous_wavelet( &time_series, sampling_interval, filter_type, filter_arg, n_octave, n_voice, scale_min, intrp_ptr, &scale, &cwt );#
 * @return Standard mutils error/OK code.
 * @param  time_series Pointer to a pre-allocated universal matrix
 *                     containing the time series to analyze.
 *                     This input must a single-row or single-column
 *                     matrix and can be of any type with the
 *                     exception of MUTIL\_DCOMPLEX.
 * @param  sampling_interval The sampling interval of the time series.
 * @param  filter_type Filter type. Only the wav\_filter\_type
 *                     WAV\_FILTER\_HAAR, WAV\_FILTER\_GAUSSIAN\_I,
 *                     WAV\_FILTER\_GAUSSIAN\_II, and WAV\_FILTER\_MORLET
 *                     are supported.
 * @param  filter_arg  A double value representing a secondary argument to be
 *                     passed to the filter function. If the filter
 *                     is of type WAV\_FILTER\_GAUSSIAN\_I or
 *                     WAV\_FILTER\_GAUSSIAN\_II then this parameter
 *                     represents the standard deviation $\sigma$ of a Guassian PDF.
 *                     If the filter is of type WAV\_FILTER\_MORLET, then
 *                     this argument represents the frequency shift.
 * @param  scale       Pointer to a pre-allocated universal matrix of type MUTIL\_DOUBLE
 *                     containing the scales over which to calculate the CWT.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  cwt         Pointer to a universal matrix of type MUTIL\_DCOMPLEX
 *                     and of size [N x Ns] where N is the length of the input
 *                     time series and Ns is the number of scales.
 *                     This matrix contains the CWT coefficients.
 *                     The memory for this matrix is allocated by the function.
 * @see _wav_filter_type
 * @see wavuniv_filters_continuous
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_continuous_wavelet(
  const univ_mat        *time_series,
  const double           sampling_interval,
  const wav_filter_type  filter_type,
  const double           filter_arg,
  const univ_mat        *scale,
  void                  *intrp_ptr,
  univ_mat              *cwt );

#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_CWT_H_*/
