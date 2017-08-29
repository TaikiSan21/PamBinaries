
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_math.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_MATH_H
#define IN_UT_MATH_H

#include <math.h>

/*
   This file contains generally useful mathematical operations.
*/

#ifdef __cplusplus
extern "C" {
#endif


/** Value of pi.
 * Constant to define pi, to 20 decimal places.  The value is taken
 * from the Solaris version of the math.h c header file.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 */
#define MUTIL_PI 3.14159265358979323846


/** Round to the nearest whole number.
 * Macro that calls floor(x + 0.5) to round a double-precision
 * number to the nearest whole number.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #y = MUTIL_ROUND(x);#
 * @return Double-precision rounded number.
 * @param   x  Number to round.
 */
#define MUTIL_ROUND(x) floor(( x ) + 0.5 )


/** Find maximum of two numbers.
 * Macro that uses the conditional expression operator (?:) to
 * take the maximum of two numbers.  If the arguments are
 * expressions, they may be evaluated twice.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #z = MUTIL_MAX(x, y);#
 * @return The maximum of the two numbers.
 * @param    x  First number.
 * @param    y  Second number.
 * @see MUTIL_MIN
 */
#define MUTIL_MAX(x, y) (( (x) > (y) ) ? (x) : (y) )


/** Find minimum of two numbers.
 * Macro that uses the conditional expression operator (?:) to
 * take the minimum of two numbers.  If the arguments are
 * expressions, they may be evaluated twice.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #z = MUTIL_MIN(x, y);#
 * @return The minimum of the two numbers.
 * @param x  First number.
 * @param y  Second number.
 * @see MUTIL_MAX
 */
#define MUTIL_MIN(x, y) (( (x) > (y) ) ? (y) : (x) )


/** Find minimum of three numbers.
 * Macro that uses the conditional expression operator (?:) to
 * take the minimum of three numbers.  If the arguments are
 * expressions, they may be evaluated twice.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #z = MUTIL_MIN(x, y, z);#
 * @return The minimum of the three numbers.
 * @param x  First number.
 * @param y  Second number.
 * @param z  Third number.
 * @see MUTIL_MIN
 * @see MUTIL_MAX
 */
#define MUTIL_MIN3(x, y, z) MUTIL_MIN( MUTIL_MIN( x, y ), z )


/** Check approximate equality of floating point numbers.
 * Checks the relative equality of two floating point numbers by checking
 * if their absolute difference is less than or equal to the sum of their
 * absolute values multiplied by 1e-N, where N is an integer provided by the
 * user.
 *
 * The formula used is |a - b| <= (|a| + |b|) x 1e-N
 *
 * This macro should NOT be used to decide if a number is equal or very close
 * to zero.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage  To check whether x and y are equal to one part per million:
 *     #equal = MUTIL_EQUAL_TO(x, y, 6);#
 * @return TRUE or FALSE.
 * @param  x  First number.
 * @param  y  Second number.
 * @param  N  The exponent in 1e-N. It should be an integer number, not a
 *    variable.
 */
#define MUTIL_EQUAL_TO(x, y, N) \
  ( fabs( (x) - (y) ) <= (( fabs(x) + fabs(y) ) * (1e-##N) ))


/** Square of a number.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #y = MUTIL_SQR(x);#
 * @return The square of the number.
 * @param   x  The input number.
 */
#define MUTIL_SQR(x) ( (x) * (x) )


/** Absolute value of a number.
 * Absolute value of a number, defined so that it works
 * for all data types.  May produce a compiler warning,
 * which can be ignored, if used on unsigned integer data.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #y = MUTIL_ABS(x);#
 * @return The absolute value of the number.
 * @param    x  The input number.
 */
#define MUTIL_ABS(x) (( (x) >= 0 ) ? (x) : -(x) )


/** Sign of a number.
 * Sign of a number, defined so that it works
 * for all data types. Returns a +1 if the number is greater than or equal to
 * zero and -1 otherwise.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #y = MUTIL_SGN(x);#
 * @return The sign of the number (+1 or -1).
 * @param    x  The input number.
 */
#define MUTIL_SGN(x) (( (x) >= 0 ) ? 1 : -1 )


/** Assign a sign to a number x given the sign of number y.
 * Assign of a sign to a number x given the sign of number y,
 * defined so that it works for all data types.
 * This function mimics the fortran function ``sign''.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #y = MUTIL_SIGN(x, y);#
 * @return If y is positive, MUTIL\_ABS(x) is returned.
 *     If y is negative, -MUTIL\_ABS(x) is returned.
 * @param    x  The input number.
 * @param    y  The input signed number.
 */
#define MUTIL_SIGN(x,y) ((y) >= 0.0 ? MUTIL_ABS(x) : -MUTIL_ABS(x))


/** Base-2 logarithm of a number.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #y = MUTIL_LOG2(x);#
 * @return The base-2 logarithm of the number.
 * @param    x  The input number.
 */
#define MUTIL_LOG2(x) ( log10( (double) (x) ) / log10( 2.0 ) )


/** Test if a number is a power of two.
 * This macro uses the \Ref{MUTIL_LOG2} and
 * \Ref{MUTIL_EQUAL_TO} macros to determine if a number is a power of two,
 * to fifteen decimal places.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #istwo = MUTIL_TEST_POW2(x);#
 * @return   TRUE or FALSE.
 * @param  x  The input number.
 * @see MUTIL_EQUAL_TO
 * @see MUTIL_LOG2
 */
#define MUTIL_TEST_POW2(x) \
    ( x > 0 && ( ceil( MUTIL_LOG2( x ) ) == floor( MUTIL_LOG2( x ) ) ) )


/** Power function $x^{y}$.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #z = MUTIL_POW(x, y);#
 * @return x raised to the power y.
 * @param    x  The input number.
 * @param    y  The input power.
 */
#define MUTIL_POW(x, y) ( pow( (double) (x), (double) (y) ) )


/** Exponential function $\exp{x}$.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #y = MUTIL_EXP(x);#
 * @return e raised to the power x.
 * @param    x  The input number.
 */
#define MUTIL_EXP(x) ( exp( (double) (x) ) )


/** Add two complex numbers.
 * Takes two complex numbers, adds them together, and puts the result in
 * a third complex number.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #MUTIL_CPX_ADD(num1, num2, result);#
 *
 * @return none
 * @param a The first complex number to be added.
 * @param b The second complex number to be added.
 * @param c The result complex number.
 * @see MUTIL_CPX_SUB
 * @see MUTIL_CPX_MULT
 * @see MUTIL_CPX_ABS
 * @see _dcomplex
 */
#define MUTIL_CPX_ADD(a, b, c) \
      (c).re = (a).re + (b).re; \
      (c).im = (a).im + (b).im



/** Subtract two complex numbers.
 * Take tow complex numbers, and subtracts the second from the first, putting
 * the result in a third complex number.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #MUTIL_CPX_SUB(num1, num2, result);#
 *
 * @return none
 * @param a The first complex number.
 * @param b The complex number to be subtracted from the first.
 * @param c The result complex number.
 * @see MUTIL_CPX_ADD
 * @see MUTIL_CPX_ABS
 * @see MUTIL_CPX_MULT
 * @see _dcomplex
 */
#define MUTIL_CPX_SUB(a, b, c) \
      (c).re = (a).re - (b).re; \
      (c).im = (a).im - (b).im



/** Multiply two complex numbers.
 * Takes two complex numbers and multiplies them together, putting the
 * result in a third complex number.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #MUTIL_CPX_MULT(num1, num2, result);#
 *
 * @return none
 * @param a The first complex number.
 * @param b The second complex number.
 * @param c The result complex number.
 * @see MUTIL_CPX_ADD
 * @see MUTIL_CPX_SUB
 * @see MUTIL_CPX_ABS
 * @see _dcomplex
 */
#define MUTIL_CPX_MULT(a, b, c) \
      (c).re = (a).re * (b).re - (a).im * (b).im; \
      (c).im = (a).re * (b).im + (a).im * (b).re



/** Absolute value of a complex number (modulus).
 * Takes a complex number and returns its absolute value (modulus),
 * defined as the square root of the sum of the squares of the real and
 * imaginary parts of the number. This macro calls the standard C function
 * sqrt.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_math.h
 * @source ut\_math.h
 * @library
 * @usage #abs_value = MUTIL_CPX_ABS(num);#
 *
 * @return A double precision number containing the absolute value.
 * @param a The complex number whose absolute value must be calculated.
 * @see MUTIL_CPX_ADD
 * @see MUTIL_CPX_SUB
 * @see MUTIL_CPX_MULT
 * @see _dcomplex
 */
#define MUTIL_CPX_ABS(a) \
      sqrt( (a).re * (a).re + (a).im * (a).im )


#ifdef __cplusplus
}
#endif

#endif /* #ifdef IN_UT_MATH_H */
