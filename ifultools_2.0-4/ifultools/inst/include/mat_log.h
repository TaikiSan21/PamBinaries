
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_log.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_LOG_H
#define IN_MAT_LOG_H

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "mat_any.h"

/* This file contains functions and macros for universal matrix operations.
   normally associated with the logic component of an ALU.
 */

#ifdef __cplusplus
extern "C" {
#endif


/*
********************************
Logical Functions
********************************
*/

/** Bitwise AND two matrices.
 * Calculate the bitwise AND of two matrices, which must have same
 * data types and equal dimensions, and put the result into a previously
 * allocated matrix of the same data type and dimensions.
 * It is legal for any or all of the three matrix pointers to point
 * to the same matrix (e.g., to assign the result to one of the operands),
 * and the operation will still work correctly.
 *
 * @limits These functions are limited to unsigned integers data types.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_log.h
 * @source mat\_log.c
 * @library matrix
 * @usage  #err_code = matuniv_bit_and(&mat1, &mat2, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param   mat1       Pointer to the first matrix operand.
 * @param   mat2       Pointer to the second matrix operand.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   result     Pointer to resulting matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_bit_and(const uint8_mat *mat1, const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_bit_and(const uint16_mat *mat1, const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_bit_and(const uint32_mat *mat1, const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_bit_or
 * @see matuniv_bit_not
 * @see matuniv_bit_xor
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_and( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_bit_and in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu8_bit_and( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result );


/* This function is documented under matuniv_bit_and in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu16_bit_and( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result );


/*  This function is documented under matuniv_bit_and in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu32_bit_and( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result );


/** Bitwise OR two matrices.
 * Calculate the bitwise OR of two matrices, which must have same
 * data types and equal dimensions, and put the result into a previously
 * allocated matrix of the same data type and dimensions.
 * It is legal for any or all of the three matrix pointers to point
 * to the same matrix (e.g., to assign the result to one of the operands),
 * and the operation will still work correctly.
 *
 * @limits These functions are limited to unsigned integers data types.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_log.h
 * @source mat\_log.c
 * @library matrix
 * @usage  #err_code = matuniv_bit_or(&mat1, &mat2, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param   mat1       Pointer to the first matrix operand.
 * @param   mat2       Pointer to the second matrix operand.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   result     Pointer to resulting matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_bit_or(const uint8_mat *mat1, const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_bit_or(const uint16_mat *mat1, const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_bit_or(const uint32_mat *mat1, const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_bit_and
 * @see matuniv_bit_not
 * @see matuniv_bit_xor
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_or( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_bit_or in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu8_bit_or( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result );


/* This function is documented under matuniv_bit_or in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu16_bit_or( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result );


/*  This function is documented under matuniv_bit_or in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu32_bit_or( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result );


/** Bitwise XOR two matrices.
 * Calculate the bitwise XOR of two matrices, which must have same
 * data types and equal dimensions, and put the result into a previously
 * allocated matrix of the same data type and dimensions.
 * It is legal for any or all of the three matrix pointers to point
 * to the same matrix (e.g., to assign the result to one of the operands),
 * and the operation will still work correctly.
 *
 * @limits These functions are limited to unsigned integers data types.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_log.h
 * @source mat\_log.c
 * @library matrix
 * @usage  #err_code = matuniv_bit_xor(&mat1, &mat2, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param   mat1       Pointer to the first matrix operand.
 * @param   mat2       Pointer to the second matrix operand.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   result     Pointer to resulting matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_bit_xor(const uint8_mat *mat1, const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_bit_xor(const uint16_mat *mat1, const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_bit_xor(const uint32_mat *mat1, const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_bit_and
 * @see matuniv_bit_not
 * @see matuniv_bit_or
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_xor( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_bit_xor in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu8_bit_xor( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result );


/* This function is documented under matuniv_bit_xor in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu16_bit_xor( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result );


/*  This function is documented under matuniv_bit_xor in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu32_bit_xor( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result );


/** Bitwise NOT a matrix.
 * Calculate the bitwise NOT of a matrix, and put the result into a previously
 * allocated matrix of the same data type and dimensions.
 * It is legal for two matrix pointers to point
 * to the same matrix (e.g., to assign the result to the operand),
 * and the operation will still work correctly.
 *
 * @limits These functions are limited to unsigned integers data types.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_log.h
 * @source mat\_log.c
 * @library matrix
 * @usage  #err_code = matuniv_bit_not(&mat, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param   mat       Pointer to the first matrix operand.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   result     Pointer to resulting matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_bit_not(const uint8_mat *mat, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_bit_not(const uint16_mat *mat, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_bit_not(const uint32_mat *mat, void *intrp_ptr, uint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_bit_and
 * @see matuniv_bit_or
 * @see matuniv_bit_xor
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_not( const univ_mat *mat,
  void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_bit_not in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu8_bit_not( const uint8_mat *mat,
  void *intrp_ptr, uint8_mat *umat_result );


/* This function is documented under matuniv_bit_not in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu16_bit_not( const uint16_mat *mat,
  void *intrp_ptr, uint16_mat *result );


/*  This function is documented under matuniv_bit_not in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu32_bit_not( const uint32_mat *mat,
  void *intrp_ptr, uint32_mat *result );


/** Bitwise AND a scalar and a matrix.
 * Perform bitwise AND between a matrix and a scalar,
 * and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for both of the input and output matrix pointers to point
 * to the same matrix (i.e., to assign the result to the matrix operand),
 * and the operation will still work correctly.
 *
 * @limits These functions are limited to unsigned integers data types.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_log.h
 * @source mat\_log.c
 * @library matrix
 * @usage  #err_code = matuniv_bit_and_scalar( &mat, scalar, intrp_ptr, &result );#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix operand.
 * @param scalar     The scalar operand.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to resulting product matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_bit_and_scalar(const uint8_mat *mat, uint8 scalar, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_bit_and_scalar(const uint16_mat *mat, uint16 scalar, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_bit_and_scalar(const uint32_mat *mat, uint32 scalar, void *intrp_ptr, uint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_bit_or_scalar
 * @see matuniv_bit_xor_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_and_scalar( const univ_mat *mat,
  univ_scalar scalar, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_bit_and_scalar in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu8_bit_and_scalar( const uint8_mat *mat,
  uint8 scalar, void *intrp_ptr, uint8_mat *result );


/* This function is documented under matuniv_bit_and_scalar in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu16_bit_and_scalar( const uint16_mat *mat,
  uint16 scalar, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_bit_and_scalar in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu32_bit_and_scalar( const uint32_mat *mat,
  uint32 scalar, void *intrp_ptr, uint32_mat *result );


/** Bitwise OR a scalar and a matrix.
 * Perform bitwise OR between a matrix and a scalar,
 * and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for both of the input and output matrix pointers to point
 * to the same matrix (i.e., to assign the result to the matrix operand),
 * and the operation will still work correctly.
 *
 * @limits These functions are limited to unsigned integers data types.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_log.h
 * @source mat\_log.c
 * @library matrix
 * @usage  #err_code = matuniv_bit_or_scalar( &mat, scalar, intrp_ptr, &result );#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix operand.
 * @param scalar     The scalar operand.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to resulting product matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_bit_or_scalar(const uint8_mat *mat, uint8 scalar, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_bit_or_scalar(const uint16_mat *mat, uint16 scalar, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_bit_or_scalar(const uint32_mat *mat, uint32 scalar, void *intrp_ptr, uint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_bit_and_scalar
 * @see matuniv_bit_xor_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_or_scalar( const univ_mat *mat,
  univ_scalar scalar, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_bit_or_scalar in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu8_bit_or_scalar( const uint8_mat *mat,
  uint8 scalar, void *intrp_ptr, uint8_mat *result );


/* This function is documented under matuniv_bit_or_scalar in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu16_bit_or_scalar( const uint16_mat *mat,
  uint16 scalar, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_bit_or_scalar in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu32_bit_or_scalar( const uint32_mat *mat,
  uint32 scalar, void *intrp_ptr, uint32_mat *result );


/** Bitwise XOR a scalar and a matrix.
 * Perform bitwise XOR between a matrix and a scalar,
 * and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for both of the input and output matrix pointers to point
 * to the same matrix (i.e., to assign the result to the matrix operand),
 * and the operation will still work correctly.
 *
 * @limits These functions are limited to unsigned integers data types.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_log.h
 * @source mat\_log.c
 * @library matrix
 * @usage  #err_code = matuniv_bit_xor_scalar( &mat, scalar, intrp_ptr, &result );#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix operand.
 * @param scalar     The scalar operand.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to resulting product matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_bit_xor_scalar(const uint8_mat *mat, uint8 scalar, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_bit_xor_scalar(const uint16_mat *mat, uint16 scalar, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_bit_xor_scalar(const uint32_mat *mat, uint32 scalar, void *intrp_ptr, uint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_bit_and_scalar
 * @see matuniv_bit_or_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_xor_scalar( const univ_mat *mat,
  univ_scalar scalar, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_bit_xor_scalar in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu8_bit_xor_scalar( const uint8_mat *mat,
  uint8 scalar, void *intrp_ptr, uint8_mat *result );


/* This function is documented under matuniv_bit_xor_scalar in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu16_bit_xor_scalar( const uint16_mat *mat,
  uint16 scalar, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_bit_xor_scalar in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu32_bit_xor_scalar( const uint32_mat *mat,
  uint32 scalar, void *intrp_ptr, uint32_mat *result );


/** Bitwise shift matrix elements to the left.
 * Perform left shift on a matrix,
 * and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for both of the input and output matrix pointers to point
 * to the same matrix (i.e., to assign the result to the matrix operand),
 * and the operation will still work correctly.
 *
 * @limits These functions are limited to unsigned integers data types.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_log.h
 * @source mat\_log.c
 * @library matrix
 * @usage  #err = matuniv_bit_shift_left( &mat, 3, intrp_ptr, &result );#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix operand.
 * @param scalar     The number of places to shift.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to resulting bit-shifted matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_bit_shift_left(const uint8_mat *mat, uint8 scalar, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_bit_shift_left(const uint16_mat *mat, uint16 scalar, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_bit_shift_left(const uint32_mat *mat, uint32 scalar, void *intrp_ptr, uint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_bit_shift_right
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_shift_left( const univ_mat *mat,
  uint8 scalar, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_bit_shift_left in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu8_bit_shift_left( const uint8_mat *mat,
  uint8 scalar, void *intrp_ptr, uint8_mat *result);


/* This function is documented under matuniv_bit_shift_left in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu16_bit_shift_left( const uint16_mat *mat,
  uint16 scalar, void *intrp_ptr, uint16_mat *result);


/* This function is documented under matuniv_bit_shift_left in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu32_bit_shift_left( const uint32_mat *mat,
  uint32 scalar, void *intrp_ptr, uint32_mat *result);


/** Bitwise shift matrix elements to the right.
 * Perform right shift on a matrix,
 * and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for both of the input and output matrix pointers to point
 * to the same matrix (i.e., to assign the result to the matrix operand),
 * and the operation will still work correctly.
 *
 * @limits These functions are limited to unsigned integers data types.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_log.h
 * @source mat\_log.c
 * @library matrix
 * @usage  #err = matuniv_bit_shift_right( &mat, 3, intrp_ptr, &result );#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix operand.
 * @param scalar     The number of places to shift.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to resulting bit-shifted matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_bit_shift_right(const uint8_mat *mat, uint8 scalar, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_bit_shift_right(const uint16_mat *mat, uint16 scalar, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_bit_shift_right(const uint32_mat *mat, uint32 scalar, void *intrp_ptr, uint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_bit_shift_left
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_shift_right( const univ_mat *mat,
  uint8 scalar, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_bit_shift_right in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu8_bit_shift_right( const uint8_mat *mat,
  uint8 scalar, void *intrp_ptr, uint8_mat *result);


/* This function is documented under matuniv_bit_shift_right in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu16_bit_shift_right( const uint16_mat *mat,
  uint16 scalar, void *intrp_ptr, uint16_mat *result);


/* This function is documented under matuniv_bit_shift_right in mat_log.h */
MUTIL_LIBEXPORT mutil_errcode matu32_bit_shift_right( const uint32_mat *mat,
  uint32 scalar, void *intrp_ptr, uint32_mat *result);


#ifdef __cplusplus
}
#endif

#endif /* ifdef IN_MAT_LOG_H */
