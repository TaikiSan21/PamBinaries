
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mth_var.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_MATH_VAR_H
#define IN_MATH_VAR_H

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
   This file contains proptypes for functions used to estimate
   distributional properties of a time series.
*/


/** One-sided autocovariance.
 * A fast implementation of the autocovariance sequence
 * using the forward and inverse forms of the discrete Fourier transform.
 * Given the sequence $\{X_t\}_{t=0}^{N-1}$, an estimate of the autocovariance sequence is given by
 *
 * \[\hat{s}_k \equiv  { 1 \over \mbox{norm} }\sum_{t = 0}^{N-|k| - 1} (X_t - \bar{X})(X_{t+|k|}-\bar{X})\]
 *
 * for $k = 0, \pm 1, \ldots, \pm (N-1)$. For the biased case norm = N while for the unbiased case norm = N - |k|.
 * A fast Fourier domain technique is used to calculate the acvs of the input sequence for lags k = 0, ..., N-1.
 *
 * References:
 *
 * 1. D.B. Percival and  A.T. Walden, ``Spectral Analysis for Physical Applications'', Cambridge University Press,
 * 1993, pp. 196-197.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_var.h
 * @source mth\_var.c
 * @library math
 * @usage #err = mthuniv_acvs( &time_series, biased, recenter, intrp_ptr, &result  );#
 * @return Standard mutils error/OK code.
 * @param time_series  A pointer to a single-column or single-row
 *                     universal matrix of type MUTIL\_DOUBLE
 *                     containing an evenly sampled time series.
 * @param biased       A boolean flag. If TRUE, the result is the biased
 *                     estimator (normalized by N, the number of samples
 *                     in the time series). If FALSE, the result is
 *                     the unbiased estimator (the t'thacvs value is normalized by
 *                     N - |t| for the unbiased case where t = 0, ..., N - 1.
 * @param recenter     A boolean flag. If TRUE, the data is recentered
 *                     by subtracting the sample mean. If FALSE, the
 *                     mean is set to zero.
 * @param intrp_ptr    Pointer for implementation of interrupt checking.
 * @param result       A pointer to a universal matrix the same and type
 *                     as the input time series which upon return will contain
 *                     the ACVS. The memory for this matrix is automatically
 *                     allocated within the function.
 * @see fra_spectral_density_function_direct
 */
MUTIL_LIBEXPORT mutil_errcode mthuniv_acvs(
  const univ_mat   *time_series,
  const boolean     biased,
  const boolean     recenter,
  void             *intrp_ptr,
  univ_mat         *result );


#ifdef __cplusplus
}
#endif

#endif /* #ifdef IN_MATH_VAR_H */
