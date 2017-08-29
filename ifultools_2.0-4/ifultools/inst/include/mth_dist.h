
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mth_dist.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_MATH_DIST_H
#define IN_MATH_DIST_H

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
   This file contains proptypes for characterizing random distributions.
*/

/** The cumulative distribution function for a random normal variable.
 * Returns the cumulative probability for a random normal variable.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_dist.h
 * @source mth\_dist.c
 * @library math
 * @usage #y = mth_pnorm( x );#
 * @return Double value.
 * @param  x Numeric value of type double representing the quantile over which
 *           the function is evaluated.
 * @see MUTIL_PNORM
 * @see MUTIL_ERF
 * @see MUTIL_ERFC
 * @see mth_erf
 * @see mth_erfc
 * @see mth_qnorm
 */
MUTIL_LIBEXPORT double mth_pnorm( double x );

/** Quantiles for a Gaussian distribution.
 * The quantile Q( p ) is the p x 100 percentage
 * point for a Gaussian random variable.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_dist.h
 * @source mth\_dist.c
 * @library math
 * @usage #result = mth_qnorm( probability );#
 * @return A double value denoting the quantile..
 * @param  probability The probability.
 * @see mth_pnorm
 * @see mth_qchisq
 */
MUTIL_LIBEXPORT double mth_qnorm( double probability );

/** Quantiles for a chi-square distribution.
 * The quantile Q( p, v ) is the p x 100 percentage
 * point for a chi-square random variable with v degrees of freedom.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mth\_dist.h
 * @source mth\_dist.c
 * @library math
 * @usage #mth_qchisq( probability, dof );#
 * @return The chi-square quantile..
 * @param  probability The probability.
 * @param  dof The degrees of freedom.
 * @see mth_qnorm
 */
MUTIL_LIBEXPORT double mth_qchisq( double probability, double dof );

#ifdef __cplusplus
}
#endif

#endif /* #ifdef IN_MATH_DIST_H */
