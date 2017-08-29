
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mth_mac.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_MATH_MAC_H
#define IN_MATH_MAC_H

#include "mth_stat.h"
#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
   This file contains macro wrappers for (scalar)
   mathematical functions.
*/


/** Macro wrapper for the gamma function.
 * This macro casts the input value to double in the call to the
 * \Ref{mth_gamma} function.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_mac.h
 * @source mth\_mac.h
 * @library math
 * @usage #y = MUTIL_GAMMA( x );#
 * @return Double-precision value.
 * @param  x Numeric value over which the function is evaluated. This
 *           input is cast to type double in the function call.
 * @see mth_gamma
 */
#define MUTIL_GAMMA( x ) mth_gamma( (double) ( x ) )


/** Macro wrapper for the digamma function.
 * This macro casts the input value to double in the call to the
 * \Ref{mth_digamma} function.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_mac.h
 * @source mth\_mac.h
 * @library math
 * @usage #y = MUTIL_DIGAMMA( x );#
 * @return Double-precision value.
 * @param  x Numeric value over which the function is evaluated. This
 *           input is cast to type double in the function call.
 * @see mth_digamma
 */
#define MUTIL_DIGAMMA( x ) mth_digamma( (double) ( x ) )


/** Macro wrapper for the trigamma function.
 * This macro casts the input value to double in the call to the
 * \Ref{mth_trigamma} function.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_mac.h
 * @source mth\_mac.h
 * @library math
 * @usage #y = MUTIL_TRIGAMMA( x );#
 * @return Double-precision value.
 * @param  x Numeric value over which the function is evaluated. This
 *           input is cast to type double in the function call.
 * @see mth_trigamma
 */
#define MUTIL_TRIGAMMA( x ) mth_trigamma( (double) ( x ) )


/** Macro wrapper for the cumulative distribution function for a random normal variable.
 * This macro casts the input value to double in the call to the
 * \Ref{mth_pnorm} function.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_mac.h
 * @source mth\_mac.h
 * @library math
 * @usage #y = MUTIL_PNORM( x );#
 * @return Double-precision value.
 * @param  x Numeric value over which the function is evaluated. This
 *           input is cast to type double in the function call.
 * @return Double-precision value.
 * @see mth_pnorm
 */
#define MUTIL_PNORM( x ) mth_pnorm( (double) ( x ) )


/** Macro wrapper for the error function.
 * This macro casts the input value to double in the call to the
 * \Ref{mth_erf} function.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_mac.h
 * @source mth\_mac.h
 * @library math
 * @usage #y = MUTIL_ERF( x );#
 * @param  x Numeric value over which the function is evaluated. This
 *           input is cast to type double in the function call.
 * @return Double-precision value.
 * @see mth_erf
 */
#define MUTIL_ERF( x ) mth_erf( (double) ( x ) )


/** Macro wrapper for the complimentary error function.
 * This macro casts the input value to double in the call to the
 * \Ref{mth_erfc} function.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_mac.h
 * @source mth\_mac.h
 * @library math
 * @usage #y = MUTIL_ERFC( x );#
 * @param  x Numeric value over which the function is evaluated. This
 *           input is cast to type double in the function call.
 * @return Double-precision value.
 * @see mth_erfc
 */
#define MUTIL_ERFC( x ) mth_erfc( (double) ( x ) )


#ifdef __cplusplus
}
#endif

#endif /* #ifdef IN_MATH_MAC_H */








