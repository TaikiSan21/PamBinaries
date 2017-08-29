
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mth_stat.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_MATH_STAT_H
#define IN_MATH_STAT_H

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
   This file contains proptypes for (scalar) mathematical functions.
*/

/** The gamma function.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_stat.h
 * @source mth\_stat.c
 * @library math
 * @usage #y = mth_gamma( x );#
 * @return Double value.
 * @param  x Numeric value of type double.
 * @see MUTIL_GAMMA
 * @see MUTIL_DIGAMMA
 * @see MUTIL_TRIGAMMA
 * @see mth_digamma
 * @see mth_trigamma
 */
MUTIL_LIBEXPORT double mth_gamma( double x );

/** The digamma function.
 * The digamma function is defined as the first derivative
 * of the natural log of the gamma function.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_stat.h
 * @source mth\_stat.c
 * @library math
 * @usage #y = mth_digamma( x );#
 * @return Double value.
 * @param  x Numeric value of type double.
 * @see MUTIL_GAMMA
 * @see MUTIL_DIGAMMA
 * @see MUTIL_TRIGAMMA
 * @see mth_gamma
 * @see mth_trigamma
 */
MUTIL_LIBEXPORT double mth_digamma( double x );

/** The trigamma function.
 * The trigamma function is defined as the first derivative
 * of the natural log of the gamma function.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_stat.h
 * @source mth\_stat.c
 * @library math
 * @usage #y = mth_trigamma( x );#
 * @return Double value.
 * @param  x Numeric value of type double.
 * @see MUTIL_GAMMA
 * @see MUTIL_DIGAMMA
 * @see MUTIL_TRIGAMMA
 * @see mth_gamma
 * @see mth_digamma
 */
MUTIL_LIBEXPORT double mth_trigamma( double x );

/** The error function.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_stat.h
 * @source mth\_stat.c
 * @library math
 * @usage #y = mth_erf( x );#
 * @return Double value.
 * @param  x Numeric value of type double over which the function is evaluated.
 * @see MUTIL_ERF
 * @see MUTIL_ERFC
 * @see mth_erfc
 */
MUTIL_LIBEXPORT double mth_erf( double x );

/** The complimentary error function.
 * Defined as mth\_erfc( x ) = 1 - mth\_erf( x ).
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_stat.h
 * @source mth\_stat.c
 * @library math
 * @usage #y = mth_erfc( x );#
 * @return Double value.
 * @param  x Numeric value of type double over which the function is evaluated.
 * @see MUTIL_ERF
 * @see MUTIL_ERFC
 * @see mth_erf
 */
MUTIL_LIBEXPORT double mth_erfc( double x );

#ifdef __cplusplus
}
#endif

#endif /* #ifdef IN_MATH_STAT_H */
