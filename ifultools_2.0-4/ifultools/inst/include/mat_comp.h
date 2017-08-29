
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_comp.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_COMP_H
#define IN_MAT_COMP_H

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "mat_any.h"

/* This file contains functions and macros for universal matrix operations
   dealing with comparisons such as number_equal, ==, <, >.
 */

#ifdef __cplusplus
extern "C" {
#endif


/*
********************************
Macros for universal matrices
********************************
*/


/** Find the element by element minimum of two matrices.
 * Find the element by element minimum of two matrices, which must have equal
 * dimensions, and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for any or all of the three matrix pointers to point
 * to the same matrix (e.g., to assign the result to one of the operands),
 * and the minimum will still work correctly.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data
 * type is capable of holding the result. For example, if
 * mat1 is of type sint8, and mat2 and result are of type uint16, then
 * negative values in the sint8 matrix which must be less than the
 * uint16 matrix will not be preserved.  The function casts both
 * inputs to the result type (if the result type is equal to or
 * ``larger'' based on MUTIL\_DATA\_TYPE).  Clipping is not permitted
 * during the cast, thus detecting the problem just described. If you
 * are unsure about the range of your input data, it is safest to use
 * double data types to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_comp.h
 * @source mat\_comp.c
 * @library matrix
 * @usage  #err_code = matuniv_min(&mat1, &mat2, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param   mat1       Pointer to the first matrix operand.
 * @param   mat2       Pointer to the second matrix operand.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   result     Pointer to resulting matrix difference.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_min(const double_mat *mat1, const double_mat *mat2, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_min(const float_mat *mat1, const float_mat *mat2, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_min(const uint8_mat *mat1, const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_min(const uint16_mat *mat1, const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_min(const uint32_mat *mat1, const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_min(const sint16_mat *mat1, const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_min(const sint32_mat *mat1, const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_max
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_min( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_min in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_min( const double_mat *mat1,
  const double_mat *mat2, void *intrp_ptr, double_mat *result);


/* This function is documented under matuniv_min in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matflt_min( const float_mat *mat1,
  const float_mat *mat2, void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_min in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats32_min( const sint32_mat *mat1,
  const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_min in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats16_min( const sint16_mat *mat1,
  const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_min in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu32_min( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);


/* This function is documented under matuniv_min in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu16_min( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_min in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu8_min( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result );


/** Find the element by element maximum of two matrices.
 * Find the element by element maximum of two matrices, which must have equal
 * dimensions, and put the result into a previously
 * allocated matrix of the same dimensions.
 * It is legal for any or all of the three matrix pointers to point
 * to the same matrix (e.g., to assign the result to one of the operands),
 * and the maximum will still work correctly.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data
 * type is capable of holding the result. For example, if
 * mat1 is of type uint8, and mat2 and result are of type sint16, then
 * very large values in the uint8 matrix which must be less than the
 * sint8 matrix will not be preserved.  The function casts both inputs
 * to the ``larger'' type based on MUTIL\_DATA\_TYPE.  Clipping is not
 * permitted during the cast, thus detecting the problem just
 * described. If you are unsure about the range of your input data, it
 * is safest to use double data types to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_comp.h
 * @source mat\_comp.c
 * @library matrix
 * @usage  #err_code = matuniv_max(&mat1, &mat2, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param   mat1       Pointer to the first matrix operand.
 * @param   mat2       Pointer to the second matrix operand.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   result     Pointer to resulting matrix difference.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_max(const double_mat *mat1, const double_mat *mat2, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_max(const float_mat *mat1, const float_mat *mat2, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_max(const uint8_mat *mat1, const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_max(const uint16_mat *mat1, const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_max(const uint32_mat *mat1, const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_max(const sint16_mat *mat1, const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_max(const sint32_mat *mat1, const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_min
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_max( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_max( const double_mat *mat1,
  const double_mat *mat2, void *intrp_ptr, double_mat *result);


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matflt_max( const float_mat *mat1,
  const float_mat *mat2, void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats32_max( const sint32_mat *mat1,
  const sint32_mat *mat2, void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats16_max( const sint16_mat *mat1,
  const sint16_mat *mat2, void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu32_max( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, uint32_mat *result);


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu16_max( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu8_max( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, uint8_mat *result );


/** Find the number of equal elements in two matrices.
 * Find the number of equal elements in two matrices, which must have equal
 * dimensions and the same type. It is legal for the two matrix pointers
 * to point to the same matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_comp.h
 * @source mat\_comp.c
 * @library matrix
 * @usage  #err_code = matuniv_number_equal(&mat1, &mat2, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param   mat1       Pointer to the first matrix operand.
 * @param   mat2       Pointer to the second matrix operand, of the same
 *                     type as mat1.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   result     Pointer for number of identical elements.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_number_equal(const double_mat *mat1, const double_mat *mat2, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_number_equal(const float_mat *mat1, const float_mat *mat2, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_number_equal(const uint8_mat *mat1, const uint8_mat *mat2, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_number_equal(const uint16_mat *mat1, const uint16_mat *mat2, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_number_equal(const uint32_mat *mat1, const uint32_mat *mat2, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_number_equal(const sint16_mat *mat1, const sint16_mat *mat2, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_number_equal(const sint32_mat *mat1, const sint32_mat *mat2, void *intrp_ptr, sint32 *result);#
 *  \end{itemize}
 * @see matuniv_number_equal_scalar
 * @see matuniv_number_less_than_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_number_equal( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_number_equal( const double_mat *mat1,
  const double_mat *mat2, void *intrp_ptr, sint32 *result);


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matflt_number_equal( const float_mat *mat1,
  const float_mat *mat2, void *intrp_ptr, sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats32_number_equal( const sint32_mat *mat1,
  const sint32_mat *mat2, void *intrp_ptr, sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats16_number_equal( const sint16_mat *mat1,
  const sint16_mat *mat2, void *intrp_ptr, sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu32_number_equal( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, sint32 *result);


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu16_number_equal( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu8_number_equal( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, sint32 *result );


/** Find the number of elements in a matrix that are equal to the given scalar.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_comp.h
 * @source mat\_comp.c
 * @library matrix
 * @usage  #err_code = matuniv_number_equal_scalar(&mat, scalar, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param   mat        Pointer to the matrix operand.
 * @param   scalar     Pointer to the scalar, the type must by compatible with
 *    that of mat.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   result     Pointer for number of identical elements.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_number_equal_scalar(const double_mat *mat, const double scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_number_equal_scalar(const float_mat *mat, const float scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_number_equal_scalar(const uint8_mat *mat, const uint8 scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_number_equal_scalar(const uint16_mat *mat, const uint16 scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_number_equal_scalar(const uint32_mat *mat, const uint32 scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_number_equal_scalar(const sint16_mat *mat, const sint16 scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_number_equal_scalar(const sint32_mat *mat, const sint32 scalar, void *intrp_ptr, sint32 *result);#
 *  \end{itemize}
 * @see matuniv_number_equal
 * @see matuniv_number_less_than_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_number_equal_scalar( const univ_mat *mat,
  const univ_scalar scalar, void *intrp_ptr, sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_number_equal_scalar( const double_mat *mat,
  const double scalar, void *intrp_ptr, sint32 *result);


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matflt_number_equal_scalar( const float_mat *mat,
  const float scalar, void *intrp_ptr, sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats32_number_equal_scalar( const sint32_mat *mat,
  const sint32 scalar, void *intrp_ptr, sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats16_number_equal_scalar( const sint16_mat *mat,
  const sint16 scalar, void *intrp_ptr, sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu32_number_equal_scalar( const uint32_mat *mat,
  const uint32 scalar, void *intrp_ptr, sint32 *result);


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu16_number_equal_scalar( const uint16_mat *mat,
  const uint16 scalar, void *intrp_ptr, sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu8_number_equal_scalar( const uint8_mat *mat,
  const uint8 scalar, void *intrp_ptr, sint32 *result );


/** Find the number of elements in a matrix that are less than
 * the given scalar.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_comp.h
 * @source mat\_comp.c
 * @library matrix
 * @usage  #err_code = matuniv_number_less_than_scalar(&mat, scalar, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param   mat        Pointer to the matrix operand.
 * @param   scalar     Pointer to the scalar, the type must by compatible with
 *    that of mat.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   result     Pointer for number of identical elements.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_number_less_than_scalar(const double_mat *mat, const double scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_number_less_than_scalar(const float_mat *mat, const float scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_number_less_than_scalar(const uint8_mat *mat, const uint8 scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_number_less_than_scalar(const uint16_mat *mat, const uint16 scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_number_less_than_scalar(const uint32_mat *mat, const uint32 scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_number_less_than_scalar(const sint16_mat *mat, const sint16 scalar, void *intrp_ptr, sint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_number_less_than_scalar(const sint32_mat *mat, const sint32 scalar, void *intrp_ptr, sint32 *result);#
 *  \end{itemize}
 * @see matuniv_number_equal
 * @see matuniv_number_equal_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_number_less_than_scalar(
  const univ_mat *mat, const univ_scalar scalar, void *intrp_ptr,
  sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_number_less_than_scalar(
  const double_mat *mat, const double scalar, void *intrp_ptr,
  sint32 *result);


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matflt_number_less_than_scalar(
  const float_mat *mat, const float scalar, void *intrp_ptr,
  sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats32_number_less_than_scalar(
  const sint32_mat *mat, const sint32 scalar, void *intrp_ptr,
  sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats16_number_less_than_scalar(
  const sint16_mat *mat, const sint16 scalar, void *intrp_ptr,
  sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu32_number_less_than_scalar(
  const uint32_mat *mat, const uint32 scalar, void *intrp_ptr,
  sint32 *result);


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu16_number_less_than_scalar(
  const uint16_mat *mat, const uint16 scalar, void *intrp_ptr,
  sint32 *result );


/* This function is documented under matuniv_max in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu8_number_less_than_scalar(
  const uint8_mat *mat, const uint8 scalar, void *intrp_ptr,
  sint32 *result );


/** Return the values and/or flattened index locations of a matrix which
 * satisfy the given scalar relation.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_comp.h
 * @source mat\_comp.c
 * @library matrix
 * @usage  #err_code = matuniv_compare_scalar(&mat, MUTIL_RELATION_LESS_THAN_OR_EQUAL, scalar, intrp_ptr, &match_index, &match_value);#
 * @return        Standard mutils error/OK code.
 * @param   mat        Pointer to the matrix operand.
 * @param   relation   An enumerated value of type \Ref{_mutil_relation} denoting
 *   the relation operator.
 * @param   scalar     Comparative scalar, the type must by compatible with
 *    that of mat.
 * @param   intrp_ptr  Pointer for implementation of interrupt checking.
 * @param   match_index Pointer to a sint32 matrix  which, upon
 *                      return, will contain the flattened indices of the input
 *                      matrix which satsify the given scalar relation. This
 *                      pointer may be set to NULL if the storage of the
 *                      matching indices is not desired. In this case, however,
 *                      the pointer to match\_value must not be NULL. The
 *                      memory for this matrix is allocated by this function
 *                      as a row vector. If no matches occur, no memory is
 *                      allocated and the matrix header elements (nelem, nrow,
 *                      and ncol) are set to zero.
 * @param   match_value Pointer to a universal matrix of the same type as
 *                      the input matrix which, upon
 *                      return, will contain the values of the input
 *                      matrix which satsify the given scalar relation. This
 *                      pointer may be set to NULL if the storage of the
 *                      matching values is not desired. In this case, however,
 *                      the pointer to match\_index must not be NULL. The
 *                      memory for this matrix is allocated by this function
 *                      as a row vector. If no matches occur, no memory is
 *                      allocated and the matrix header elements (nelem, nrow,
 *                      and ncol) are set to zero.
 * @limits Complex input matrices are not supported.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_compare_scalar(const double_mat *mat, const mutil_relation relation, const double scalar, void *intrp_ptr, sint32_mat *match_index, double_mat *match_value);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_compare_scalar(const float_mat *mat, const mutil_relation relation, const float scalar, void *intrp_ptr, sint32_mat *match_index, float_mat *match_value);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_compare_scalar(const uint8_mat *mat, const mutil_relation relation, const uint8 scalar, void *intrp_ptr, sint32_mat *match_index, uint8_mat *match_value);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_compare_scalar(const uint16_mat *mat, const mutil_relation relation, const uint16 scalar, void *intrp_ptr, sint32_mat *match_index, uint16_mat *match_value);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_compare_scalar(const uint32_mat *mat, const mutil_relation relation, const uint32 scalar, void *intrp_ptr, sint32_mat *match_index, unit32_mat *match_value);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_compare_scalar(const sint16_mat *mat, const mutil_relation relation, const sint16 scalar, void *intrp_ptr, sint32_mat *match_index, sint16_mat *match_value);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_compare_scalar(const sint32_mat *mat, const mutil_relation relation, const sint32 scalar, void *intrp_ptr, sint32_mat *match_index, sint32_mat *match_value);#
 *  \end{itemize}
 * @see matuniv_number_equal
 * @see matuniv_number_equal_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_compare_scalar(
  const univ_mat       *mat,
  const mutil_relation  relation,
  const univ_scalar     scalar,
  void                 *intrp_ptr,
  sint32_mat           *match_index,
  univ_mat             *match_value );

/* This function is documented under matuniv_compare_scalar in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_compare_scalar(
  const double_mat *mat, const mutil_relation relation, const double scalar,
  void *intrp_ptr, sint32_mat *match_index, double_mat *match_value);


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matflt_compare_scalar(
  const float_mat *mat, const mutil_relation  relation, const float scalar,
  void *intrp_ptr, sint32_mat *match_index, float_mat *match_value );


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats32_compare_scalar(
  const sint32_mat *mat, const mutil_relation  relation, const sint32 scalar,
  void *intrp_ptr, sint32_mat *match_index, sint32_mat *match_value );


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode mats16_compare_scalar(
  const sint16_mat *mat, const mutil_relation  relation, const sint16 scalar,
  void *intrp_ptr, sint32_mat *match_index, sint16_mat *match_value );


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu32_compare_scalar(
  const uint32_mat *mat, const mutil_relation  relation, const uint32 scalar,
  void *intrp_ptr, sint32_mat *match_index, uint32_mat *match_value);


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu16_compare_scalar(
  const uint16_mat *mat, const mutil_relation  relation, const uint16 scalar,
  void *intrp_ptr, sint32_mat *match_index, uint16_mat *match_value );


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
MUTIL_LIBEXPORT mutil_errcode matu8_compare_scalar(
  const uint8_mat *mat, const mutil_relation  relation, const uint8 scalar,
  void *intrp_ptr, sint32_mat *match_index, uint8_mat *match_value );


#ifdef __cplusplus
}
#endif

#endif /* ifdef IN_MAT_COMP_H */
