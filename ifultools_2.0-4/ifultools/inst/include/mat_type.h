
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_type.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_TYPE_H_
#define IN_MAT_TYPE_H_

#include "ut_type.h"
#include "ut_err.h"

/* This file contains typedefs and structs for the basic matrix types
   for the mutils library.
   All of the matrices are stored in row-major order.
*/

#ifdef __cplusplus
extern "C" {
#endif


/** Struct for matrix of 8-bit unsigned integers.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same
 *   #typedef struct _uint8_mat uint8_mat;#
 * @see uint8
 */
struct _uint8_mat {
    /** number of rows */
    sint32 nrow;

    /** number of columns */
    sint32 ncol;

    /** total number of elements */
    sint32 nelem;

    /** pointer to flat data array in row-major order */
    uint8 *data;
};


/* See documentation for _uint8_mat (above) for description */
typedef struct _uint8_mat uint8_mat;


/** Struct for matrix of 8-bit signed integers.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same
 *   #typedef struct _sint8_mat sint8_mat;#
 * @see sint8
 */
struct _sint8_mat {
    /** number of rows */
    sint32 nrow;

    /** number of columns */
    sint32 ncol;

    /** total number of elements */
    sint32 nelem;

    /** pointer to flat data array in row-major order */
    sint8 *data;
};


/* See documentation for _sint8_mat (above) for description */
typedef struct _sint8_mat sint8_mat;


/** Struct for matrix of 16-bit unsigned integers.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same
 *   #typedef struct _uint16_mat uint16_mat;#
 * @see uint16
 */
struct _uint16_mat {
    /** number of rows */
    sint32 nrow;

    /** number of columns */
    sint32 ncol;

    /** total number of elements */
    sint32 nelem;

    /** pointer to flat data array in row-major order */
    uint16 *data;
};


/* See documentation for _uint16_mat (above) for description */
typedef struct _uint16_mat uint16_mat;


/** Struct for matrix of 16-bit signed integers.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same
 *   #typedef struct _sint16_mat sint16_mat;#
 * @see sint16
 */
struct _sint16_mat {
    /** number of rows */
    sint32 nrow;

    /** number of columns */
    sint32 ncol;

    /** total number of elements */
    sint32 nelem;

    /** pointer to flat data array in row-major order */
    sint16 *data;
};


/* See documentation for _sint16_mat (above) for description */
typedef struct _sint16_mat sint16_mat;


/** Struct for matrix of 32-bit unsigned integers.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same
 *   #typedef struct _uint32_mat uint32_mat;#
 * @see uint32
 */
struct _uint32_mat {
    /** number of rows */
    sint32 nrow;

    /** number of columns */
    sint32 ncol;

    /** total number of elements */
    sint32 nelem;

    /** pointer to flat data array in row-major order */
    uint32 *data;
};


/* See documentation for _uint32_mat (above) for description */
typedef struct _uint32_mat uint32_mat;


/** Struct for matrix of 32-bit signed integers.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same
 *   #typedef struct _sint32_mat sint32_mat;#
 * @see sint32
 */
struct _sint32_mat {
    /** number of rows */
    sint32 nrow;

    /** number of columns */
    sint32 ncol;

    /** total number of elements */
    sint32 nelem;

    /** pointer to flat data array in row-major order */
    sint32 *data;
};


/* See documentation for _sint32_mat (above) for description */
typedef struct _sint32_mat sint32_mat;


/** Struct for matrix of single-precision floating-point numbers.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same
 *   #typedef struct _float_mat float_mat;#
 */
struct _float_mat {
    /** number of rows */
    sint32 nrow;

    /** number of columns */
    sint32 ncol;

    /** total number of elements */
    sint32 nelem;

    /** pointer to flat data array in row-major order */
    float *data;
};


/* See documentation for _float_mat (above) for description */
typedef struct _float_mat float_mat;


/** Struct for matrix of double-precision floating-point numbers.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same
 *   #typedef struct _double_mat double_mat;#
 */
struct _double_mat {
    /** number of rows */
    sint32 nrow;

    /** number of columns */
    sint32 ncol;

    /** total number of elements */
    sint32 nelem;

    /** pointer to flat data array in row-major order */
    double *data;
};


/* See documentation for _double_mat (above) for description */
typedef struct _double_mat double_mat;


/** Struct for matrix of complex numbers.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @see _dcomplex
 * @same
 *   #typedef struct _dcomplex_mat dcomplex_mat;#
 */
struct _dcomplex_mat {
    /** number of rows */
    sint32 nrow;

    /** number of columns */
    sint32 ncol;

    /** total number of elements */
    sint32 nelem;

    /** pointer to flat data array in row-major order */
    dcomplex *data;
};


/* See documentation for _dcomplex_mat above for information
 * about this struct. */
typedef struct _dcomplex_mat dcomplex_mat;


/** Matrix union for the universal matrix.
 * Union of basic matrix types used in the universal matrix.  This
 * union is not meant to be used by itself.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @see _uint8_mat
 * @see _sint8_mat
 * @see _uint16_mat
 * @see _sint16_mat
 * @see _uint32_mat
 * @see _sint32_mat
 * @see _float_mat
 * @see _double_mat
 * @see _dcomplex_mat
 * @see _mutil_data_type
 * @see _univ_mat
 * @see Basic Matrix Functions
 */
union _mutil_univ_union
{
  /** matrix of 8-bit unsigned integers */
  uint8_mat  u8mat;

  /** matrix of 8-bit signed integers */
  sint8_mat  s8mat;

  /** matrix of 16-bit unsigned integers */
  uint16_mat u16mat;

  /** matrix of 16-bit signed integers */
  sint16_mat s16mat;

  /** matrix of 32-bit unsigned integers */
  uint32_mat u32mat;

  /** matrix of 32-bit signed integers */
  sint32_mat s32mat;

  /** matrix of single-precision floating-point numbers */
  float_mat  fltmat;

  /** matrix of double-precision floating-point numbers */
  double_mat dblmat;

  /** matrix of double-precision complex numbers */
  dcomplex_mat cpxmat;
};


/** Struct for the universal matrix.
 * The universal matrix is a wrapper structure that can contain any of
 * the matrix data types defined in the mutils library.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same
 *   #typedef struct _univ_mat univ_mat;#
 * @see _mutil_univ_union
 * @see _mutil_data_type
 * @see Basic Matrix Functions
 */
struct _univ_mat {
    /** the wrapped matrix */
    union _mutil_univ_union mat;

    /** the data type of the wrapped matrix */
    mutil_data_type  type;
};


/* See documentation for _univ_mat (above) for description */
typedef struct _univ_mat univ_mat;


/** Struct for a matrix set.
 * The matrix set is a wrapper structure that contains a single-
 * or multi-dimensional set of universal matrices, which could
 * represent frequency bands or colors for images, time-indexed images,
 * depth-indexed images, or other structures.  The matrix set
 * is stored in a flat array, and if it is multi-dimensional, the
 * order is such that the last dimension varies the fastest (as in
 * C multi-dimensional arrays).
 *
 * For example, three-hundred grayscale slices of a 200 row by 100
 * column image would be stored in a matrix set with ndim equal to 1,
 * dims set to be the array [300],  and nelem set to 300.
 * The underlying universal matrices
 * would all have 200 rows and 100 columns.
 *
 * For example, five RGB color slices of a 128 row by 128
 * column image would be stored in a matrix set with ndim equal to 2,
 * dims set to be the array [5, 3], and nelem set to 15.
 * The underlying universal matrices
 * would all have 128 rows and 128 columns.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same
 *   #typedef struct _mat_set mat_set;#
 * @see _univ_mat
 * @see Matrix Set Functions
 */
struct _mat_set {
  /** number of dimensions in array */
  sint32 ndim;

  /** dimensions of the array */
  sint32 *dims;

  /** pointer to flat array of matrices */
  univ_mat *mats;

  /** total number of matrices in flat array, which is equivalent to the
   product of dims */
  sint32 nelem;

  /** if TRUE, the matrices are allocated in a contiguous block of memory */
  boolean contiguous;
};


/* See documentation for _mat_set (above) for description */
typedef struct _mat_set mat_set;

/** Enum of boundary conditions for matrices.
 * These represent various choices of
 * what to do at the boundaries of matrices.
 * They are useful for image and signal processing
 * operations, such as convolution, correlation, and wavelet transforms.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same #typedef enum _mutil_boundary_type mutil_boundary_type;#
 */
enum _mutil_boundary_type
{
  /** zero boundary condition: pixels beyond matrix boundary are
    treated as zeros */
  MUTIL_BOUNDARY_ZERO,

  /** periodic boundary condition: pixels beyond matrix boundary are treated
    as a periodic (tiled) continuation of the matrix */
  MUTIL_BOUNDARY_PERIODIC,

  /** reflective boundary condition: pixels beyond the matrix boundary are
    as if reflected by mirrors along the edges of the matrix */
  MUTIL_BOUNDARY_REFLECT,

  /** continuation boundary condition: pixels beyond the matrix boundary are
    taken to be the last value inside the matrix */
  MUTIL_BOUNDARY_CONTINUE
};

/* See above for documentation on this enum. */
typedef enum _mutil_boundary_type mutil_boundary_type;


/** Enum for interpolation algorithms.
 * These represent various choices of
 * interpolation algorithms.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same #typedef enum _mutil_interpolation_type mutil_interpolation_type;#
 */
enum _mutil_interpolation_type
{
  /** no interpolation */
  MUTIL_INTERPOLATION_NONE,

  /** linear interpolation */
  MUTIL_INTERPOLATION_LINEAR
};

/* See above for documentation on this enum. */
typedef enum _mutil_interpolation_type mutil_interpolation_type;

/** Enum of relation comparison operators.
 * These represent the supported relational operators
 * used by matrix comparison functions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_type.h
 * @source mat\_type.h
 * @library matrix
 * @same #typedef enum _mutil_relation mutil_relation;#
 */
enum _mutil_relation{
  /** less than operator */
  MUTIL_RELATION_LESS_THAN,

  /** less than or equal operator */
  MUTIL_RELATION_LESS_THAN_OR_EQUAL,

  /** equivalence operator */
  MUTIL_RELATION_EQUAL,

  /** non-equivalence operator */
  MUTIL_RELATION_NOT_EQUAL,

  /** greater than operator */
  MUTIL_RELATION_GREATER_THAN,

  /** greater than or equal operator */
  MUTIL_RELATION_GREATER_THAN_OR_EQUAL
};

/* See above for documentation on this enum. */
typedef enum _mutil_relation mutil_relation;


#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_TYPE_H */
