
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/src/wav_math.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_WAV_MATH_H
#define IN_WAV_MATH_H

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
   This file contains mathematical operations for wavelet functions.
*/

/** Table lookup using linear interpolation for 1D matrices.
 * This function uses linear interpolation to find the yi, which are
 * the values of the underlying function y at the points indicated by
 * the vector xi.  The vector x specifies the points at which the data
 * y is given. The vector x must be monotonic and evenly spaced.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_intp.h
 * @source wav\_intp.c
 * @library wavelets
 * @usage #err = wavuniv_statistic_interpolation_linear(&y, &x, &xi, intrp_ptr, &yi);#
 * @return Standard mutils error/OK code.
 * @param  y           Pointer to a single column or row matrix of type
 *                     MUTIL\_DOUBLE
 *                     containing the dependent variable (table) data.
 * @param  x           Pointer to a single column or row matrix of type
 *                     MUTIL\_DOUBLE
 *                     containing the independent variable (table) data.
 *                     This matrix must be the
 *                     same size and type as the y matrix.
 * @param  xi          Pointer to a single column or row matrix of type
 *                     MUTIL\_DOUBLE
 *                     containing the independent variable (lookup) data.
 *                     This matrix must be the
 *                     same size and type as the yi matrix.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  yi          Pointer to a single column or row matrix of type
 *                     MUTIL\_DOUBLE containing the
 *                     interpolated dependent variable (lookup) data.
 * @private
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_statistic_interpolation_linear(
  const univ_mat *y,
  const univ_mat *x,
  const univ_mat *xi,
  void           *intrp_ptr,
  univ_mat       *yi );

#ifdef __cplusplus
}
#endif

#endif /* #ifdef IN_WAV_MATH_H */








