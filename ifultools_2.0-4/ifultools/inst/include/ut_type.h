
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_type.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_TYPE_H_
#define IN_UT_TYPE_H_

/* This file contains typedefs and structs for the basic data types used
   in the mutils library. */

#ifdef __cplusplus
extern "C" {
#endif

/* USING S.h to differentiate R and SPLUS */

#ifdef USE_RINTERNALS
  #include "R.h"
  #undef boolean
#endif

/** Enum of basic data types for the mutils library.
 * This lists the scalar-like types that can be used for
 * universal matrices and scalars in the mutils library.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @same
 *   #typedef enum _mutil_data_type mutil_data_type;#
 * @see _univ_scalar
 * @see Matrix Data Types
 */
enum _mutil_data_type {

    /** 8-bit unsigned integer */
    MUTIL_UINT8,

    /** 8-bit signed integer */
    MUTIL_SINT8,

    /** 16-bit unsigned integer */
    MUTIL_UINT16,

    /** 16-bit signed integer */
    MUTIL_SINT16,

    /** 32-bit unsigned integer */
    MUTIL_UINT32,

    /** 32-bit signed integer */
    MUTIL_SINT32,

    /** single-precision floating-point number */
    MUTIL_FLOAT,

    /** double-precision floating-point number */
    MUTIL_DOUBLE,

    /** double-precision complex number */
    MUTIL_DCOMPLEX
};


/* See the documentation _mutil_data_type above for a
 * description of the enum. */
typedef enum _mutil_data_type mutil_data_type;


/** 8-bit unsigned integer.
 * This type is defined as the ANSI
 * unsigned char data type, which is guaranteed to be one byte in
 * ANSI C.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @see sint8
 * @see _uint8_mat
 * @see Data Type Characterization
 */
typedef unsigned char uint8;


/** 8-bit signed integer.
 * This type is defined as the ANSI
 * signed char data type, which is guaranteed to be one byte in
 * ANSI C.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @see uint8
 * @see _sint8_mat
 * @see Data Type Characterization
 */
typedef signed char sint8;


/** 16-bit unsigned integer.
 * This type is defined to be the ANSI
 * unsigned short data type, which is guaranteed to be at least 2 bytes.
 * It could actually be larger on some platforms.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @see sint16
 * @see _uint16_mat
 * @see Data Type Characterization
 */
typedef unsigned short uint16;


/** 16-bit signed integer.
 * This type is defined to be the ANSI
 * signed short data type, which is guaranteed to be at least 2 bytes.
 * It could actually be larger on some platforms.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @see uint16
 * @see _sint16_mat
 * @see Data Type Characterization
 */
typedef signed short sint16;


/** 32-bit unsigned integer.
 * This type is defined to be the ANSI
 * unsigned Sint data type, which is guaranteed to be at least 4 bytes.
 * It could actually be larger on some platforms.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @see sint32
 * @see _uint32_mat
 * @see Data Type Characterization
 */
#ifdef USING_R
	typedef unsigned int uint32;
#else 
	typedef unsigned long uint32;
#endif

/** 32-bit signed integer.
 * This type is defined to be the ANSI
 * signed Sint data type, which is guaranteed to be at least 4 bytes.
 * It could actually be larger on some platforms.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @see uint32
 * @see _sint32_mat
 * @see Data Type Characterization
 */
#ifdef USING_R
	typedef signed int sint32;
#else
  typedef signed long sint32;
#endif

/** Struct for complex number.
 * The complex number is a struct with double real and imaginary parts.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @same
 *   #typedef struct _dcomplex dcomplex;#
 * @see _dcomplex_mat
 */
struct _dcomplex {
  /** real part */
  double re;

  /** imaginary part */
  double im;
};


/* See the documentation for _dcomplex above for a
 * description of the struct. */
typedef struct _dcomplex dcomplex;


/** Boolean (true/false) type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @see TRUE
 * @see FALSE
 */
typedef unsigned char boolean;


/** Scalar union for the universal scalar.
 * Union of basic scalar types used in the universal scalar.  This
 * union is not meant to be used by itself, but is part of the
 * universal scalar struct.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @see _mutil_data_type
 * @see _univ_scalar
 * @see Basic Matrix Functions
 */
union _mutil_univscal_union
{
  /** 8-bit unsigned integer */
  uint8  u8;

  /** 8-bit signed integer */
  sint8  s8;

  /** 16-bit unsigned integer */
  uint16 u16;

  /** 16-bit signed integer */
  sint16 s16;

  /** 32-bit unsigned integer */
  uint32 u32;

  /** 32-bit signed integer */
  sint32 s32;

  /** single-precision floating-point number */
  float  flt;

  /** double-precision floating-point number */
  double dbl;

  /** double-precision complex number */
  dcomplex cpx;
};


/** Struct for the universal scalar.
 * The universal scalar is a wrapper structure that can contain any of
 * the scalar data types defined in the mutils library.  It is used
 * to pass scalar data to and from universal matrix functions that
 * need different types of scalar data for different input matrix types.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @same
 *   #typedef struct _univ_scalar univ_scalar;#
 * @see _mutil_univscal_union
 * @see _mutil_data_type
 * @see Basic Matrix Functions
 */
struct _univ_scalar {
    /** the wrapped scalar */
    union _mutil_univscal_union num;

    /** the data type of the wrapped scalar */
    mutil_data_type  type;
};


/* See the documentation for _univ_scalar above for a description
 * of the struct */
typedef struct _univ_scalar univ_scalar;


#ifdef DEF_TF


/** Boolean true value (1).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @see FALSE
 * @see boolean
 * @see Compile-Time Options
 */
#define TRUE 1


/** Boolean false value (0).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 * @see TRUE
 * @see boolean
 * @see Compile-Time Options
 */
#define FALSE 0


#endif /* ifdef DEF_TF */


/** Code for invalid length (-1).
 * This number designates an invalid number of rows, columns, or elements.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_type.h
 * @source ut\_type.h
 * @library
 */
#define MUTIL_INVALID_LENGTH -1


#ifdef __cplusplus
}
#endif

#endif /*ifndef IN_UT_TYPE_H*/
