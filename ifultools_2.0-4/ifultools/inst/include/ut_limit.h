
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_limit.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_LIMIT_H_
#define IN_UT_LIMIT_H_

/* This file contains ranges for the basic data types used
   in the mutils library.
*/

/** Maximum allowed value for uint8.
 * The maximum value (255) that the mutils library will store in a
 * uint8 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see uint8
 */
#define MUTIL_UINT8_MAX 255


/** Minimum allowed value for sint8.
 * The minimum value (-127) that the mutils library will store in an
 * sint8 value.  This is the minimum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see sint8
 * @see MUTIL_SINT8_MAX
 */
#define MUTIL_SINT8_MIN -127


/** Maximum allowed value for sint8.
 * The maximum value (127) that the mutils library will store in an
 * sint8 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see sint8
 * @see MUTIL_SINT8_MIN
 */
#define MUTIL_SINT8_MAX 127


/** Maximum allowed value for uint16.
 * The maximum value (65535) that the mutils library will store in a
 * uint16 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see uint16
 */
#define MUTIL_UINT16_MAX 65535


/** Minimum allowed value for sint16.
 * The minimum value (-32767) that the mutils library will store in an
 * sint16 value.  This is the minimum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see sint16
 * @see MUTIL_SINT16_MAX
 */
#define MUTIL_SINT16_MIN -32767


/** Maximum allowed value for sint16.
 * The maximum value (32767) that the mutils library will store in an
 * sint16 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see sint16
 * @see MUTIL_SINT16_MIN
 */
#define MUTIL_SINT16_MAX 32767


/** Maximum allowed value for uint32.
 * The maximum value (4294967295) that the mutils library will store in a
 * uint32 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see uint32
 */
#define MUTIL_UINT32_MAX 4294967295UL


/** Minimum allowed value for sint32.
 * The minimum value (-2147483647) that the mutils library will store in an
 * sint32 value.  This is the minimum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see sint32
 * @see MUTIL_SINT32_MAX
 */
#define MUTIL_SINT32_MIN -2147483647L


/** Maximum allowed value for sint32.
 * The maximum value (2147483647) that the mutils library will store in an
 * sint32 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see sint32
 * @see MUTIL_SINT32_MIN
 */
#define MUTIL_SINT32_MAX 2147483647L


/** Maximum allowed value for floating-point number.
 * The maximum absolute value (1E+37) that the mutils library will store in
 * a float value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 */
#define MUTIL_FLOAT_MAX 1E+37


/** Maximum allowed value for double floating-point number.
 * The maximum value (1E+37) that the mutils library will store in an
 * double value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 */
#define MUTIL_DOUBLE_MAX 1E+37


/** Smallest number such that, for double numbers, 1 + x does not equal 1.
 * The smallest value (1E-9) that, for any double number x greater than
 * or equal to IE-9, 1.0 + x does not equal 1.0.  This is the value that
 * is guaranteed to work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 */
#define MUTIL_DOUBLE_EPSILON 1E-9


/** Smallest number such that, for float numbers, 1 + x does not equal 1.
 * The smallest value (1E-9) that, for any float number x greater than
 * or equal to IE-9, 1.0 + x does not equal 1.0.  This is the value that
 * is guaranteed to work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 */
#define MUTIL_FLOAT_EPSILON 1E-9


/** Number of bits in the uint8 type.
 * The number of bits (8) that the mutils library will store in a
 * uint8 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see uint8
 */
#define MUTIL_UINT8_NBIT 8


/** Number of unsigned bits in the sint8 type.
 * The number of unsigned bits (7) that the mutils library will store in a
 * sint8 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see sint8
 */
#define MUTIL_SINT8_NBIT 7


/** Number of bits in the uint16 type.
 * The number of bits (16) that the mutils library will store in a
 * uint16 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see uint16
 */
#define MUTIL_UINT16_NBIT 16


/** Number of unsigned bits in the sint16 type.
 * The number of unsigned bits (15) that the mutils library will store in a
 * sint16 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see sint16
 */
#define MUTIL_SINT16_NBIT 15


/** Number of bits in the uint32 type.
 * The number of bits (32) that the mutils library will store in a
 * uint32 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see uint32
 */
#define MUTIL_UINT32_NBIT 32


/** Number of unsigned bits in the sint32 type.
 * The number of unsigned bits (31) that the mutils library will store in a
 * sint32 value.  This is the maximum value that is guaranteed to
 * work in an ANSI-compatible compiler.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 * @see sint32
 */
#define MUTIL_SINT32_NBIT 31


/** Maximum allowed value for size\_t.
 * The maximum value (4294967295) that the mutils library will store in a
 * size\_t value.  This assumes that size\_t is at least a 32 bit integer.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 */
#define MUTIL_SIZE_T_MAX 4294967295UL


/** Maximum allowed value for double angle, in radians.
 * The maximum value that the mutils library will accept for an
 * angle input, in radians.
 * This is the value above which the smallest allowable changes
 * (up to machine precision) in angle are of order 2 pi, so that
 * trigonometric functions no longer make sense.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_limit.h
 * @source ut\_limit.h
 * @library
 */
#define MUTIL_DOUBLE_ANGLE_MAX ( 1.0 / MUTIL_DOUBLE_EPSILON )

#endif /*ifndef IN_UT_LIMIT_H*/
