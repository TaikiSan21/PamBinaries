
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_err.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_ERR_H_
#define IN_UT_ERR_H_

#include "ut_plat.h"

/* This file contains the standard error codes for the mutils library */

#ifdef __cplusplus
extern "C" {
#endif


/** Enum of standard error codes for mutils library.
 * These are the standard error codes for the mutils library.
 * They are used as the return values for all mutils functions.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_err.h
 * @source ut\_err.h
 * @library matrix
 * @same
 *   #typedef enum _mutil_errcode mutil_errcode;#
 * @see mutil_err_string
 */
enum _mutil_errcode
{
  /** No error: function succeeded. */
  MUTIL_ERR_OK = 0,

  /** Memory allocation error: out of dynamic memory or too much requested. */
  MUTIL_ERR_MEM_ALLOC,

  /** Passed in a null pointer where not allowed. */
  MUTIL_ERR_NULL_POINTER,

  /** Illegal memory address for one or more parameter arguments. */
  MUTIL_ERR_ILLEGAL_ADDRESS,

  /** Illegal or inconsistent data size. */
  MUTIL_ERR_ILLEGAL_SIZE,

  /** Illegal data type for this operation. */
  MUTIL_ERR_ILLEGAL_TYPE,

  /** Illegal data value. */
  MUTIL_ERR_ILLEGAL_VALUE,

  /** Attempted access past data bounds. */
  MUTIL_ERR_OUT_OF_BOUNDS,

  /** Attempted to divide by zero. */
  MUTIL_ERR_DIVIDE_BY_ZERO,

  /** Input/output error. */
  MUTIL_ERR_INPUT_OUTPUT,

  /** Data overflow: result cannot be stored in given data type. */
  MUTIL_ERR_OVERFLOW,

  /** User interrupted calculation. */
  MUTIL_ERR_INTERRUPT,

  /** Singular or nearly-singular matrix. */
  MUTIL_ERR_SINGULAR_MATRIX,

  /** Iterative algorithm is not converging. */
  MUTIL_ERR_NOT_CONVERGING,

  /** Cumulative roundoff errors. */
  MUTIL_ERR_CUM_ROUNDOFF,

  /** Singularity encountered in data. */
  MUTIL_ERR_SINGULARITY,

  /** Tree-structured algorithm has error in structure. */
  MUTIL_ERR_TREE_STRUCTURE,

  /** No neighbors found in a nearest neighbor search. */
  MUTIL_ERR_ZERO_NEIGHBORS_FOUND,

  /** Feature has not been implemented. */
  MUTIL_ERR_FEATURE_NOT_IMPLEMENTED
};


/* See the documentation for _mutil_errcode (above) for a
 * description of the enum. */
typedef enum _mutil_errcode mutil_errcode;


/** Character strings for standard errors.
 * Look up and return a pointer to a character string corresponding to
 * one of the standard error codes.  The character string cannot
 * be freed or modified by the calling function.
 *
 *
 * @return   Standard mutils error/OK code.
 * @param  msg_code  Error code to look up.
 * @param  message   Character string for error code.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_err.h
 * @source ut\_err.c
 * @library matrix
 * @see _mutil_errcode
 */
MUTIL_LIBEXPORT mutil_errcode mutil_err_string(mutil_errcode msg_code,
  const char **message);


#ifdef __cplusplus
}
#endif

#endif /*ifndef IN_UT_ERR_H_*/
