
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_boot.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_WAV_BOOT_H_
#define IN_WAV_BOOT_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "wav_type.h"

/* This file contains function declarations the wavelet mutils library. */

#ifdef __cplusplus
extern "C" {
#endif

/** Find the whitest set of DWPT crystals.
 * Performs a top down (level zero on up) analyis of a discrete wavelet packet
 * transform (using natural crystal indexing) and
 * finds the whitest transform of the many transsforms
 * available in a DWPT.
 *
 * References:
 *
 * 1. D. B. Percival, S. Sardy and A. C. Davison, {\it Wavestrapping Time Series:
 * Adaptive Wavelet-Based Bootstrapping}, In W. J. Fitzgerald, R. L. Smith,
 * A. T. Walden and P. C. Young (Eds.), Nonlinear and Nonstationary Signal
 * Processing, Cambridge, England: Cambridge University Press.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_boot.h
 * @source wav\_boot.c
 * @library wavelets
 * @usage #err = mutil_errcode wavuniv_transform_packet_whitest( &dwpt, significance, white_noise_test, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param  dwpt Pointer to a pre-allocated matrix set of type MUTIL\_DOUBLE
 *              containing a discrete wavelet packet transform
 *              such as that retuned by \Ref{wavuniv_transform_packet}.
 * @param  significance  The significance level to use in calculating comparative chi-square distribution
 *                     p x 100 percentage points where p = 1 - significance (the chi-square degrees of
 *                     freedom are estimated automatically within the specified white noise test).
 * @param  white_noise_test An enum of type \Ref{_wav_white_test} denoting the white noise
 *                     test to use in testing wavelet packet crystals.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  result      Pointer to a two-row universal matrix of
 *                     type MUTIL\_SINT32 which upon return will contain
 *                     the sequency ordered coordinates of the whitest set of
 *                     crystals in the DWPT bssed on the supplied
 *                     white\_noise\_test.test. The rows of this output
 *                     contain the decomposition levels and the oscillation
 *                     indices for the whitest set. Each column thus defines
 *                     the location of one such crystal in the DWPT in sequency
 *                     order. The memory for this matrix is allocated within the
 *                     function.
 * @see _wav_white_test
 * @see wavuniv_bootstrap
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_packet_inverse
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_packet_whitest(
  const mat_set  *dwpt,
  const double    significance,
  const wav_white_test  white_noise_test,
  void           *intrp_ptr,
  univ_mat       *result );


/** Adaptive wavelet-based bootstrapping.
 * Generates boostrapped version(s) of a time series
 * by (1) finding the whitest transform of the many
 * transforms in a discrete wavelet packet transform,
 * (2) randomly shuffling (with replacement) the coefficients within
 * each crystal of the whitest transform, and (3)
 * performing an inverse transform of the result.

 * References:
 *
 * 1. D. B. Percival, S. Sardy and A. C. Davison, {\it Wavestrapping Time Series:
 * Adaptive Wavelet-Based Bootstrapping}, In W. J. Fitzgerald, R. L. Smith,
 * A. T. Walden and P. C. Young (Eds.), Nonlinear and Nonstationary Signal
 * Processing, Cambridge, England: Cambridge University Press.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_boot.h
 * @source wav\_boot.c
 * @library wavelets
 * @usage #err = mutil_errcode wavuniv_bootstrap( &dwpt, &filters, &white_indices, n_realization, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param  dwpt Pointer to a pre-allocated matrix set of type MUTIL\_DOUBLE
 *              containing a discrete wavelet packet transform
 *              such as that retuned by \Ref{wavuniv_transform_packet}.
 * @param white_indices Pointer to a pre-allocated universal matrix of type
 *         MUTIL\_SINT32 containing the indices that define the
 *         whitest transform of the DWPT such as that returned
 *         by the \Ref{wavuniv_transform_packet_whitest} function.
 *         The wavelet packet crystals are assumed to be arranged in the so-called
 *         natural order ala [W\_{0,0} , W\_{1,0}, W\_{1,1}, W\_{2,0}, W\_{2,1},
 *         W\_{2,2}, W\_{2,3}, ... , W\_{J,0}, ...  W\_{J, N\_J}]
 *         where J  is the number of decomposition levels and $N\_J = 2^J$.
 *         By definition, W\_{0,0} is the original time series. A given crystal
 *         is identified in the W\_{j,n} form above or by a flattened index.
 *         For example, the DWPT crystal in level 2 at oscillation index 1 is
 *         identified either by {j,n} = {2,1} or by its flattened index 4
 *         (with zero based indexing, 4 represents the fifth crystal of the
 *         wavelet packet transform in natural order). If the crystals
 *         are identified using flattened indices, then the white\_indices
 *         matrix must be a single-column or single-row matrix. Conversely, if the
 *         crystals are identified ala the {j,n} style, then the white\_indices
 *         matrix must be either a two-column matrix with columns representing
 *         the j and n indices, respectively, or a two-row matrix with rows
 *         representing j and n indices, respectively.
 * @param n_realization The numer of realizations (bootstrapped series) to create.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result  Pointer to a matrix set containing n\_realization universal matrices
 *                 of type MUTIL\_DOUBLE and dimension 1 x N where N
 *                 is the number of samples in the original time series.
 *                 each universal matrix represents one bootstrapped
 *                 version of the original time series.
 *                 The memory for this matrix set is allocated within the function.
 * @see wavuniv_transform_packet_whitest
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_packet_inverse
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_bootstrap(
  const mat_set    *dwpt,
  const mat_set    *filters,
  const univ_mat   *white_indices,
  const sint32      n_realization,
  void             *intrp_ptr,
  mat_set          *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_BOOT_H_*/


