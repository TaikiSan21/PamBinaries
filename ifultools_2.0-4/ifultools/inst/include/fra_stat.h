
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_stat.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_STAT_H_
#define IN_FRA_STAT_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "fra_type.h"

/* This file contains function declarations the mutils fractal library. */

#ifdef __cplusplus
extern "C" {
#endif



/** Stationarity tests for a time series.
 *
 *
 * References:
 *
 * 1.Priestley, Subba Rao
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_stat.h
 * @source fra\_stat.c
 * @library fractal
 * @usage #err = frauniv_stationarity_priestley_subba_rao( &time_series, sampling_interval, n_taper, n_block, significance, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series.
 * @param sampling_interval  The sampling interval of the time series.
 * @param n_taper     The number of tapers to use in forming the eigenspectra.
 * @param n_block     The number of nonoverlapping uniform blocks over which the
 *                    time series is divided.
 * @param significance The significance is the number of times you expect the
 *                    underlying hypothesis of stationarity to fail even
 *                    though stationarity remains true. Essentially, you are
 *                    allowing for error in the result. A significance of
 *                    0.05 means that you are allowing 5% error.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a matrix set which upon return will contain
 *                    three universal matrices of type MUTIL\_DOUBLE.
 *                    The first matrix is a 2 x 1 vector containing the
 *                    'between time' and 'interaction + residual' sum of
 *                    squares statistics. The second is a ( n\_block + 1 ) x
 *                    ( n\_freq + 1 ) matrix containing the Priestley-Subba Rao
 *                    ANOVA results (n\_freq is automatically determined within
 *                    the program based on the bandwidth estimations.
 *                    The n\_freq + 1 column contains the
 *                    row means, the n\_block + 1 row contains the column means,
 *                    and the last element in the matrix contains the grand mean.
 *                    The final matrix is a n\_freq x 1 vector containing the
 *                    Fourier frequencies corresponding to the first n\_freq
 *                    columns of the ANOVA matrix. The memory for the
 *                    matrix set is automatically allocated within the function.
 *
 * @see frauniv_determinism_delta_epsilon
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_stationarity_priestley_subba_rao(
  const univ_mat      *time_series,
  const double         sampling_interval,
  const sint32         n_taper,
  const sint32         n_block,
  const double         significance,
  const boolean        center,
  const boolean        recenter,
  void                *intrp_ptr,
  mat_set             *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_STAT_H_ */




