
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_arit.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_ARIT_H
#define IN_MAT_ARIT_H

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "mat_any.h"

/* This file contains functions and macros for universal matrix operations
   normally associated with the arithmetic component of an ALU.
 */

#ifdef __cplusplus
extern "C" {
#endif


/** Add two matrices.
 * Calculate the sum of two matrices, which must have equal
 * dimensions, and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for any or all of the three matrix pointers to point
 * to the same matrix (e.g., to assign the result to one of the operands),
 * and the addition will still work correctly.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data type
 * is capable of holding the arithmetic result.  If you are unsure about
 * the range of your input data, it is safest to use double data types
 * to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_arit.h
 * @source mat\_arit.c
 * @library matrix
 * @usage  #err_code = matuniv_add(&mat1, &mat2, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param   mat1       Pointer to the first matrix operand.
 * @param   mat2       Pointer to the second matrix operand.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   result     Pointer to resulting matrix sum.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_add(const double_mat *mat1, const double_mat *mat2, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_add(const float_mat *mat1, const float_mat *mat2, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_add(const uint8_mat *mat1, const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_add(const uint16_mat *mat1, const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_add(const uint32_mat *mat1, const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_add(const sint16_mat *mat1, const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_add(const sint32_mat *mat1, const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_add_scalar
 * @see matuniv_subtract
 * @see matuniv_multiply_elem
 * @see matuniv_multiply
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_add( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_add in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_add( const double_mat *mat1,
  const double_mat *mat2, void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_add in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matflt_add( const float_mat *mat1,
  const float_mat *mat2, void *intrp_ptr, float_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats32_add( const sint32_mat *mat1,
  const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats16_add( const sint16_mat *mat1,
  const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu32_add( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu16_add( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_add in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu8_add( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result );


/** Add a scalar and a matrix.
 * Perform an addition between a matrix and a scalar,
 * and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for both of the input and output matrix pointers to point
 * to the same matrix (i.e., to assign the result to the matrix operand),
 * and the addition will still work correctly.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data type
 * is capable of holding the arithmetic result.  If you are unsure about
 * the range of your input data, it is safest to use double data types
 * to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_arit.h
 * @source mat\_arit.c
 * @library matrix
 * @usage  #err_code = matuniv_add_scalar( &mat, mynum, intrp_ptr, &result );#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix operand.
 * @param scalar     The scalar operand.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to resulting product matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_add_scalar(const double_mat *mat, double scalar, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_add_scalar(const float_mat *mat, float scalar, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_add_scalar(const uint8_mat *mat, uint8 scalar, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_add_scalar(const uint16_mat *mat, uint16 scalar, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_add_scalar(const uint32_mat *mat, uint32 scalar, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_add_scalar(const sint16_mat *mat, sint16 scalar, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_add_scalar(const sint32_mat *mat, sint32 scalar, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_add
 * @see matuniv_multiply_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_add_scalar( const univ_mat *mat,
  univ_scalar scalar, void *intrp_ptr, univ_mat *result);


/* This function is documented under matuniv_add_scalar in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_add_scalar( const double_mat *mat,
  double scalar, void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_add_scalar in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matflt_add_scalar( const float_mat *mat,
  float scalar, void *intrp_ptr, float_mat *result );

/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats32_add_scalar( const sint32_mat *mat,
  sint32 scalar, void *intrp_ptr, sint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats16_add_scalar( const sint16_mat *mat,
  sint16 scalar, void *intrp_ptr, sint16_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu32_add_scalar( const uint32_mat *mat,
  uint32 scalar, void *intrp_ptr, uint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu16_add_scalar( const uint16_mat *mat,
  uint16 scalar, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_add_scalar in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu8_add_scalar( const uint8_mat *mat,
  uint8 scalar, void *intrp_ptr, uint8_mat *result );


/** Subtract two matrices.
 * Calculate the difference between two matrices, which must have equal
 * dimensions, and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for any or all of the three matrix pointers to point
 * to the same matrix (e.g., to assign the result to one of the operands),
 * and the subtraction will still work correctly.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data type
 * is capable of holding the arithmetic result.  If you are unsure about
 * the range of your input data, it is safest to use double data types
 * to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_arit.h
 * @source mat\_arit.c
 * @library matrix
 * @usage  #err_code = matuniv_subtract(&mat1, &mat2, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param   mat1       Pointer to the first matrix operand.
 * @param   mat2       Pointer to the second matrix operand.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   result     Pointer to resulting matrix difference.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_subtract(const double_mat *mat1, const double_mat *mat2, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_subtract(const float_mat *mat1, const float_mat *mat2, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_subtract(const uint8_mat *mat1, const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_subtract(const uint16_mat *mat1, const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_subtract(const uint32_mat *mat1, const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_subtract(const sint16_mat *mat1, const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_subtract(const sint32_mat *mat1, const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_add
 * @see matuniv_multiply_elem
 * @see matuniv_multiply
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_subtract( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_subtract in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_subtract( const double_mat *mat1,
  const double_mat *mat2, void *intrp_ptr, double_mat *result);


/* This function is documented under matuniv_subtract in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matflt_subtract( const float_mat *mat1,
  const float_mat *mat2, void *intrp_ptr, float_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats32_subtract( const sint32_mat *mat1,
  const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats16_subtract( const sint16_mat *mat1,
  const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu32_subtract( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu16_subtract( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_subtract in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu8_subtract( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result );


/** Matrix multiplication.
 * Perform matrix multiplication on two matrices, which
 * must have compatible dimensions, and put the
 * result into a previously allocated matrix with the
 * correct dimensions for the matrix product.
 * The output matrix must not have the same data location as
 * either input matrix.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data type
 * is capable of holding the arithmetic result.  If you are unsure about
 * the range of your input data, it is safest to use double data types
 * to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_arit.h
 * @source mat\_arit.c
 * @library matrix
 * @usage  #err_code = matuniv_multiply(&mat1, &mat2, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param mat1       Pointer to the first matrix operand.
 * @param mat2       Pointer to the second matrix operand.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to resulting matrix product.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_multiply(const double_mat *mat1, const double_mat *mat2, void *intrp_ptr, double_mat *result);#
 * \item #MUTIL_LIBEXPORT mutil_errcode matflt_multiply(const float_mat *mat1, const float_mat *mat2, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_multiply(const uint8_mat *mat1, const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_multiply(const uint16_mat *mat1, const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_multiply(const uint32_mat *mat1, const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_multiply(const sint16_mat *mat1, const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_multiply(const sint32_mat *mat1, const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_add
 * @see matuniv_subtract
 * @see matuniv_multiply_elem
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_multiply( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_multiply in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_multiply(const double_mat *mat1,
  const double_mat *mat2, void *intrp_ptr, double_mat *result);


/* This function is documented under matuniv_multiply in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matflt_multiply( const float_mat *mat1,
  const float_mat *mat2, void *intrp_ptr, float_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats32_multiply( const sint32_mat *mat1,
  const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats16_multiply( const sint16_mat *mat1,
  const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu32_multiply( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu16_multiply( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_multiply in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu8_multiply( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result );


/** Element-wise multiplication of two matrices.
 * Perform element-by-element multiplication of two
 * matrices,  which must have equal
 * and positive dimensions, and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for any or all of the three matrix pointers to point
 * to the same matrix (e.g., to assign the result to one of the operands),
 * and the multiplication will still work correctly.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data type
 * is capable of holding the arithmetic result.  If you are unsure about
 * the range of your input data, it is safest to use double data types
 * to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_arit.h
 * @source mat\_arit.c
 * @library matrix
 * @usage  #err_code = matuniv_multiply_elem(&mat1, &mat2, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param mat1       Pointer to the first matrix operand.
 * @param mat2       Pointer to the second matrix operand.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to resulting product matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_multiply_elem(const double_mat *mat1, const double_mat *mat2, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_multiply_elem(const float_mat *mat1, const float_mat *mat2, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_multiply_elem(const uint8_mat *mat1, const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_multiply_elem(const uint16_mat *mat1, const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_multiply_elem(const uint32_mat *mat1, const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_multiply_elem(const sint16_mat *mat1, const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_multiply_elem(const sint32_mat *mat1, const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_multiply
 * @see matuniv_multiply_scalar
 * @see matuniv_add
 * @see matuniv_subtract
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_multiply_elem( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_multiply_elem in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_multiply_elem( const double_mat *mat1,
  const double_mat *mat2, void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_multiply_elem in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matflt_multiply_elem( const float_mat *mat1,
  const float_mat *mat2, void *intrp_ptr, float_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats32_multiply_elem( const sint32_mat *mat1,
  const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats16_multiply_elem( const sint16_mat *mat1,
  const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu32_multiply_elem( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu16_multiply_elem( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_multiply_elem in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu8_multiply_elem( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result );


/** Multiply a scalar and a matrix.
 * Perform multiplication between a matrix and a scalar,
 * and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for both of the input and output matrix pointers to point
 * to the same matrix (i.e., to assign the result to the matrix operand),
 * and the multiplication will still work correctly.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data type
 * is capable of holding the arithmetic result.  If you are unsure about
 * the range of your input data, it is safest to use double data types
 * to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_arit.h
 * @source mat\_arit.c
 * @library matrix
 * @usage  #err_code = matuniv_multiply_scalar( &mat, mynum, intrp_ptr, &result );#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix operand.
 * @param scalar     The scalar operand.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to resulting product matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_multiply_scalar(const double_mat *mat, double scalar, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_multiply_scalar(const float_mat *mat, float scalar, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_multiply_scalar(const uint8_mat *mat, uint8 scalar, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_multiply_scalar(const uint16_mat *mat, uint16 scalar, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_multiply_scalar(const uint32_mat *mat, uint32 scalar, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_multiply_scalar(const sint16_mat *mat, sint16 scalar, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_multiply_scalar(const sint32_mat *mat, sint32 scalar, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_multiply
 * @see matuniv_add_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_multiply_scalar( const univ_mat *mat,
  univ_scalar scalar, void *intrp_ptr, univ_mat *result);


/* This function is documented under matuniv_multiply_scalar in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_multiply_scalar( const double_mat *mat,
  double scalar, void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_multiply_scalar in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matflt_multiply_scalar( const float_mat *mat,
  float scalar, void *intrp_ptr, float_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats32_multiply_scalar( const sint32_mat *mat,
  sint32 scalar, void *intrp_ptr, sint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats16_multiply_scalar( const sint16_mat *mat,
  sint16 scalar, void *intrp_ptr, sint16_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu32_multiply_scalar( const uint32_mat *mat,
  uint32 scalar, void *intrp_ptr, uint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu16_multiply_scalar( const uint16_mat *mat,
  uint16 scalar, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_multiply_scalar in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu8_multiply_scalar( const uint8_mat *mat,
  uint8 scalar, void *intrp_ptr, uint8_mat *result );


/** Divide a scalar by a matrix or a matrix by a scalar.
 * If all elements of the denominator are nonzero,
 * then each output element is the requested scalar divided by the corresponding
 * input element or input element divided by the scalar.
 * The result is put into a previously allocated matrix of the same dimensions.
 * It is legal for both of the input and output matrix pointers to point
 * to the same matrix (i.e., to assign the result to the matrix operand),
 * and the division will still work correctly.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data type
 * is capable of holding the arithmetic result.  If you are unsure about
 * the range of your input data, it is safest to use double data types
 * to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_arit.h
 * @source mat\_arit.c
 * @library matrix
 * @usage  #err_code = matuniv_divide_scalar( &mat, scalar, TRUE, intrp_ptr, &result );#
 * @return        Standard mutils error/OK code.
 * @param mat           Pointer to the matrix operand.
 * @param scalar        The scalar operand.
 * @param mat_numerator If TRUE, the input matrix elements are divided by the scalar.
 *     If FALSE the scalar is divided by each input matrix element.
 * @param intrp_ptr     Pointer for implementation of interrupt checking.
 * @param result        Pointer to resulting product matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_divide_scalar(const double_mat *mat, double scalar, boolean mat_numerator, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_divide_scalar(const float_mat *mat, float scalar, boolean mat_numerator, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_divide_scalar(const uint8_mat *mat, uint8 scalar, boolean mat_numerator, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_divide_scalar(const uint16_mat *mat, uint16 scalar, boolean mat_numerator, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_divide_scalar(const uint32_mat *mat, uint32 scalar, boolean mat_numerator, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_divide_scalar(const sint16_mat *mat, sint16 scalar, boolean mat_numerator, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_divide_scalar(const sint32_mat *mat, sint32 scalar, boolean mat_numerator, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_divide_elem
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_divide_scalar( const univ_mat *mat,
  univ_scalar scalar, boolean mat_numerator, void *intrp_ptr,
  univ_mat *result);


/* This function is documented under matuniv_divide_scalar in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_divide_scalar( const double_mat *mat,
  double scalar, boolean mat_numerator, void *intrp_ptr,
  double_mat *result );


/* This function is documented under matuniv_divide_scalar in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matflt_divide_scalar( const float_mat *mat,
  float scalar, boolean mat_numerator, void *intrp_ptr,
  float_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats32_divide_scalar( const sint32_mat *mat,
  sint32 scalar, boolean mat_numerator, void *intrp_ptr,
  sint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats16_divide_scalar( const sint16_mat *mat,
  sint16 scalar, boolean mat_numerator, void *intrp_ptr,
  sint16_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu32_divide_scalar( const uint32_mat *mat,
  uint32 scalar, boolean mat_numerator, void *intrp_ptr,
  uint32_mat *result );


/* Documented with the corresponding universal matrix function in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu16_divide_scalar( const uint16_mat *mat,
  uint16 scalar, boolean mat_numerator, void *intrp_ptr,
  uint16_mat *result );


/* This function is documented under matuniv_divide_scalar in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu8_divide_scalar( const uint8_mat *mat,
  uint8 scalar, boolean mat_numerator, void *intrp_ptr,
  uint8_mat *result );


/** Divide element by element a matrix by another matrix.
 * If all elements of the denominator are nonzero,
 * then each output element is the corresponding
 * element of the first input matrix divided by the corresponding element
 * of the second input matrix.
 * The result is put into a previously allocated matrix of the same dimensions.
 * It is legal for both of the input and output matrix pointers to point
 * to the same matrix (i.e., to assign the result to the matrix operand),
 * and the division will still work correctly.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data type
 * is capable of holding the arithmetic result.  If you are unsure about
 * the range of your input data, it is safest to use double data types
 * to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_arit.h
 * @source mat\_arit.c
 * @library matrix
 * @usage  #err_code = matuniv_divide_scalar( &mat, scalar, intrp_ptr, &result );#
 * @return        Standard mutils error/OK code.
 * @param mat1          Pointer to the first matrix operand.
 * @param mat2          Pointer to the second matrix operand.
 * @param intrp_ptr     Pointer for implementation of interrupt checking.
 * @param result        Pointer to resulting product matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_divide_elem(const double_mat *mat1, const double_mat *mat2, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_divide_elem(const float_mat *mat1, const float_mat *mat2, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_divide_elem(const uint8_mat *mat1, const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_divide_elem(const uint16_mat *mat1, const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_divide_elem(const uint32_mat *mat1, const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_divide_elem(const sint16_mat *mat1, const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_divide_elem(const sint32_mat *mat1, const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_divide_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_divide_elem( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, univ_mat *result);


/* This function is documented under matuniv_divide_elem in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_divide_elem( const double_mat *mat1,
  const double_mat *mat2, void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_divide_elem in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matflt_divide_elem( const float_mat *mat1,
  const float_mat *mat2, void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_divide_elem in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats32_divide_elem( const sint32_mat *mat1,
  const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_divide_elem in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats16_divide_elem( const sint16_mat *mat1,
  const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_divide_elem in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu32_divide_elem( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result );


/* This function is documented under matuniv_divide_elem in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu16_divide_elem( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_divide_elem in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu8_divide_elem( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result );


/** Absolute value.
 * Computes the absolute value of a matrix,
 * considering the matrix data as a single flat vector.
 * The output matrix may have the same data location as
 * the input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_arit.h
 * @source mat\_arit.c
 * @library matrix
 * @usage  #err_code = matuniv_abs(&mat_in, intrp_ptr, &result_mat);#
 * @return          Standard mutils error/OK code.
 * @param mat_in   Pointer to the input matrix.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param mat_out  Pointer to the resulting matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_abs(const double_mat *mat_in, void *intrp_ptr, double_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_abs(const float_mat *mat_in, void *intrp_ptr, float_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_abs(const sint16_mat *mat_in, void *intrp_ptr, sint16_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_abs(const sint32_mat *mat_in, void *intrp_ptr, sint32_mat *mat_out);#
 *  \end{itemize}
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_abs(
    const univ_mat *mat_in,
    void           *intrp_ptr,
    univ_mat       *mat_out );


/* This function is documented under matuniv_abs in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_abs( const double_mat *mat_in,
  void *intrp_ptr, double_mat *mat_out );


/* This function is documented under matuniv_abs in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matflt_abs( const float_mat *mat_in,
  void *intrp_ptr, float_mat *mat_out );


/* This function is documented under matuniv_abs in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats32_abs( const sint32_mat *mat_in,
  void *intrp_ptr, sint32_mat *mat_out );


/* This function is documented under matuniv_abs in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats16_abs( const sint16_mat *mat_in,
  void *intrp_ptr, sint16_mat *mat_out );


/** Rescale a matrix.
 * Perform a rescaling of the data in a matrix,
 * and put the result into a previously
 * allocated matrix of the same dimensions and type.
 * The elements of the input matrix with the minimum value
 * are rescaled to min\_val in the resulting matrix;
 * the elements of the input matrix with the maximum value
 * are rescaled to max\_val in the resulting matrix;
 * all other values are proportionally rescaled
 * to values between min\_val and max\_val.
 * It is legal for both of the input and output matrix pointers to point
 * to the same matrix (i.e., to assign the result to the matrix operand).
 *
 * If all values in the input matrix are the same they are rescaled
 * to min\_val.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_arit.h
 * @source mat\_arit.c
 * @library matrix
 * @usage  #err_code = matuniv_rescale(&mat, min_val, max_val, intrp_ptr, &mat_out);#
 * @return           Standard mutils error/OK code.
 * @param mat_in     Pointer to the matrix operand.
 * @param min_val    The minimum value in the resulting matrix,
 *    of the same data type as mat\_in.
 * @param max_val    The maximum value in the resulting matrix,
 *    of the same data type as mat\_in.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param mat_out     Pointer to resulting matrix,
 *    of the same data type as mat\_in.
 *
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_rescale(const double_mat *mat_in, const double min_val, const double max_val, void *intrp_ptr, double_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_rescale(const float_mat *mat_in, const float min_val, const float max_val, void *intrp_ptr, float_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_rescale(const sint16_mat *mat_in, const sint16 min_val, const sint16 max_val, void *intrp_ptr, sint16_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_rescale(const sint32_mat *mat_in, const sint32 min_val, const sint32 max_val, void *intrp_ptr, sint32_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_rescale(const uint8_mat *mat_in, const uint8 min_val, const uint8 max_val, void *intrp_ptr, uint8_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_rescale(const uint16_mat *mat_in, const uint16 min_val, const uint16 max_val, void *intrp_ptr, uint16_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_rescale(const uint32_mat *mat_in, const uint32 min_val, const uint32 max_val, void *intrp_ptr, uint32_mat *mat_out);#
 *  \end{itemize}
 *
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_rescale(
  const univ_mat *mat_in, const univ_scalar min_val,
  const univ_scalar max_val, void *intrp_ptr,
  univ_mat *mat_out );

/* This function is documented under matuniv_rescale in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_rescale( const double_mat *mat_in,
  const double min_val, const double max_val,
  void *intrp_ptr, double_mat *mat_out );


/* This function is documented under matuniv_rescale in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matflt_rescale( const float_mat *mat_in,
  const float min_val, const float max_val,
  void *intrp_ptr, float_mat *mat_out );


/* This function is documented under matuniv_rescale in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu8_rescale( const uint8_mat *mat_in,
  const uint8 min_val, const uint8 max_val,
  void *intrp_ptr, uint8_mat *mat_out );


/* This function is documented under matuniv_rescale in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu16_rescale( const uint16_mat *mat_in,
  const uint16 min_val, const uint16 max_val,
  void *intrp_ptr, uint16_mat *mat_out );


/* This function is documented under matuniv_rescale in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode matu32_rescale( const uint32_mat *mat_in,
  const uint32 min_val, const uint32 max_val,
  void *intrp_ptr, uint32_mat *mat_out );


/* This function is documented under matuniv_rescale in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats16_rescale( const sint16_mat *mat_in,
  const sint16 min_val, const sint16 max_val,
  void *intrp_ptr, sint16_mat *mat_out );


/* This function is documented under matuniv_rescale in mat_arit.h */
MUTIL_LIBEXPORT mutil_errcode mats32_rescale( const sint32_mat *mat_in,
  const sint32 min_val, const sint32 max_val,
  void *intrp_ptr, sint32_mat *mat_out );


#ifdef __cplusplus
}
#endif

#endif /* ifdef IN_MAT_ARIT_H */
