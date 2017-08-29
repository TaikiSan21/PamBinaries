
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_nnls.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_NNLS_H_
#define IN_IMAT_NNLS_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

/* This file contains function declarations for matrix nonnegative
  least squares.  The functions are defined in mat_nnls.h.  */

#ifdef __cplusplus
extern "C" {
#endif

/** Nonnegative least squares.
 * Nonnegative least squares (NNLS) solver for ax=b
 * where a in an m by n matrix, b is a vector of length m,
 * and x is a vector of length n, subject to x >= 0.
 *
 * This code is based on Fortran code developed by
 * Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
 * June 15, 1973.  See
 * {\em Solving Least Squares Problems}, Prentice-Hall, 1974.
 * Revised February 1995 to accompany reprinting of the book by SIAM.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_nnls.h
 * @source mat\_nnls.c
 * @library matrix
 * @usage  #err_code = matdbl_non_negative_least_squares(&a, &b, &intrp_ptr, &x);#
 * @return           Standard mutils error/OK code.
 * @param  a         Pointer to m by n matrix a.
 * @param  b         Pointer to vector of length m.
 * @param  intrp_ptr Pointer for implementation of interrupt checking.
 * @param  x         Pointer to the solution of length n.
 * @see Matrix Data Types
 * @see Interrupt Handling
*/
MUTIL_LIBEXPORT mutil_errcode matdbl_non_negative_least_squares(
  const double_mat *a,
  const double_mat *b,
  void             *intrp_ptr,
  double_mat       *x );


#ifdef __cplusplus
}
#endif

#endif /* IN_IMG_SMA_H_*/
