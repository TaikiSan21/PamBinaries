
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_assn.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_ASSN_H
#define IN_MAT_ASSN_H

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "mat_any.h"

/* This file contains functions and macros for universal matrix operations.
 */

#ifdef __cplusplus
extern "C" {
#endif


/** Assign data from one matrix to another.
 * Copy the data values in one matrix to another
 * matrix previously allocated to the same size and type.
 * The operation is performed using the ANSI C memmove function,
 * so it will work even if the arrays overlap.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_assn.h
 * @source mat\_assn.c
 * @library matrix
 * @usage  #err_code = matuniv_assign(&mat, intrp_ptr, &result);#
 * @return       Standard mutils error/OK code.
 * @param mat       Pointer to the matrix to be copied.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param result    Pointer to the copy of the matrix of the same size
 *     and type.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_assign(const double_mat *mat, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_assign(const float_mat *mat, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_assign(const uint8_mat *mat, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_assign(const uint16_mat *mat, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_assign(const uint32_mat *mat, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_assign(const sint16_mat *mat, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_assign(const sint32_mat *mat, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_cast
 * @see matuniv_assign_zeropad
 * @see matuniv_assign_submat
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_assign( const univ_mat *mat,
  void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_assign in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_assign( const double_mat *mat,
  void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_assign in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matflt_assign( const float_mat *mat,
  void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_assign in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats32_assign( const sint32_mat *mat,
  void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_assign in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats16_assign( const sint16_mat *mat,
  void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_assign in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu32_assign( const uint32_mat *mat,
  void *intrp_ptr, uint32_mat *result );


/* This function is documented under matuniv_assign in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu16_assign( const uint16_mat *mat,
  void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_assign in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu8_assign( const uint8_mat *mat,
  void *intrp_ptr, uint8_mat *result );


/** Assign a scalar to every element of a matrix.
 * Copy a single scalar value to every element of a previously
 * allocated matrix of the same type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_assn.h
 * @source mat\_assn.c
 * @library matrix
 * @usage  #err_code = matuniv_assign_scalar(uscal, intrp_ptr, &mymat);#
 * @return        Standard mutils error/OK code.
 * @param scalar     Scalar value to copy.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param mat        Pointer to matrix to copy scalar into.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_assign_scalar(double scalar, void *intrp_ptr, double_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_assign_scalar(float scalar, void *intrp_ptr, float_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_assign_scalar(uint8 scalar, void *intrp_ptr, uint8_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_assign_scalar(uint16 scalar, void *intrp_ptr, uint16_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_assign_scalar(uint32 scalar, void *intrp_ptr, uint32_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_assign_scalar(sint16 scalar, void *intrp_ptr, sint16_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_assign_scalar(sint32 scalar, void *intrp_ptr, sint32_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matcpx_assign_scalar(dcomplex scalar, void *intrp_ptr, dcomplex_mat *mat);#
 *  \end{itemize}
 * @see matuniv_assign
 * @see _univ_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_assign_scalar( univ_scalar scalar,
  void *intrp_ptr, univ_mat *mat );


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_assign_scalar( double scalar,
  void *intrp_ptr, double_mat *mat );


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matflt_assign_scalar( float scalar,
  void *intrp_ptr, float_mat *mat );


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats32_assign_scalar( sint32 scalar,
  void *intrp_ptr, sint32_mat *mat );


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats16_assign_scalar( sint16 scalar,
  void *intrp_ptr, sint16_mat *mat );


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu32_assign_scalar( uint32 scalar,
  void *intrp_ptr, uint32_mat *mat );


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu16_assign_scalar( uint16 scalar,
  void *intrp_ptr, uint16_mat *mat );


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu8_assign_scalar( uint8 scalar,
  void *intrp_ptr, uint8_mat *mat );


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_assign_scalar( dcomplex scalar,
  void *intrp_ptr, dcomplex_mat *mat );


/** Assign sub-matrix values into a matrix.
 * Take the values in a smaller matrix, and assign them into
 * a larger matrix of the same type, putting the upper left corner
 * of the sub-matrix values at given row and column indices.
 * The output matrix must not have the same data location as
 * the input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_assn.h
 * @source mat\_assn.c
 * @library matrix
 * @usage  #err_code = matuniv_assign_submat(&smallmat, row, col, intrp_ptr, &largemat);#
 * @return          Standard mutils error/OK code.
 * @param smallmat  Pointer to the input matrix.
 * @param row       Row index of upper left corner of assignment location.
 * @param col       Column index of upper left corner of assignment location.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param result    Pointer to the resulting larger matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_assign_submat(const double_mat *smallmat, sint32 row, sint32 col, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_assign_submat(const float_mat *smallmat, sint32 row, sint32 col, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_assign_submat(const uint8_mat *smallmat, sint32 row, sint32 col, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_assign_submat(const uint16_mat *smallmat, sint32 row, sint32 col, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_assign_submat(const uint32_mat *smallmat, sint32 row, sint32 col, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_assign_submat(const sint16_mat *smallmat, sint32 row, sint32 col, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_assign_submat(const sint32_mat *smallmat, sint32 row, sint32 col, void *intrp_ptr, sint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matcpx_assign_submat(const dcomplex_mat *smallmat, sint32 row, sint32 col, void *intrp_ptr, dcomplex_mat *result);#
 *  \end{itemize}
 * @see matuniv_assign
 * @see matuniv_assign_zeropad
 * @see matuniv_extract
 * @see matuniv_translate
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_assign_submat( const univ_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_assign_submat in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_assign_submat( const double_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_assign_submat in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matflt_assign_submat( const float_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_assign_submat in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats32_assign_submat( const sint32_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_assign_submat in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats16_assign_submat( const sint16_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_assign_submat in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu32_assign_submat( const uint32_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, uint32_mat *result );


/* This function is documented under matuniv_assign_submat in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu16_assign_submat( const uint16_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_assign_submat in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu8_assign_submat( const uint8_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, uint8_mat *result );


/* This function is documented under matuniv_assign_submat in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_assign_submat(
  const dcomplex_mat *smallmat, sint32 row, sint32 col, void *intrp_ptr,
  dcomplex_mat *result );


/** Assign matrix values and zeros into a matrix.
 * Take the values in a smaller matrix, and assign them into
 * a larger matrix of the same type at the same row and column indices.  Fill
 * the rest of the larger matrix with zeros.
 * The output matrix must not have the same data location as
 * the input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_assn.h
 * @source mat\_assn.c
 * @library matrix
 * @usage  #err_code = matuniv_assign_zeropad(&smallmat, intrp_ptr, &largemat);#
 * @return          Standard mutils error/OK code.
 * @param smallmat  Pointer to the input matrix.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param result    Pointer to the resulting larger matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_assign_zeropad(const double_mat *smallmat, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_assign_zeropad(const float_mat *smallmat, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_assign_zeropad(const uint8_mat *smallmat, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_assign_zeropad(const uint16_mat *smallmat, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_assign_zeropad(const uint32_mat *smallmat, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_assign_zeropad(const sint16_mat *smallmat, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_assign_zeropad(const sint32_mat *smallmat, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_assign
 * @see matuniv_assign_submat
 * @see matuniv_extract
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_assign_zeropad( const univ_mat *smallmat,
  void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_assign_zeropad( const double_mat *smallmat,
  void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matflt_assign_zeropad( const float_mat *smallmat,
  void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats32_assign_zeropad( const sint32_mat *smallmat,
  void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats16_assign_zeropad( const sint16_mat *smallmat,
  void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu32_assign_zeropad( const uint32_mat *smallmat,
  void *intrp_ptr, uint32_mat *result );


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu16_assign_zeropad( const uint16_mat *smallmat,
  void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu8_assign_zeropad( const uint8_mat *smallmat,
  void *intrp_ptr, uint8_mat *result );


/** Extract a submatrix from a matrix.
 * Extract a continuous, rectangular submatrix from a matrix, starting
 * from the given row and column number, and
 * put the result into a previously allocated matrix of the same type.
 * The size of the extracted part of the matrix is taken from
 * the dimensions of the passed-in storage for the result.
 * The output matrix must not have the same data location as
 * the input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_assn.h
 * @source mat\_assn.c
 * @library matrix
 * @usage  #err_code = matuniv_extract(&mymat, start_row, start_col, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param mat       Pointer to the input matrix.
 * @param start_row The starting row number.
 * @param start_col The starting column number.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param result    Pointer to the resulting submatrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_extract(const double_mat *mat, sint32 start_row, sint32 start_col, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_extract(const float_mat *mat, sint32 start_row, sint32 start_col, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_extract(const uint8_mat *mat, sint32 start_row, sint32 start_col, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_extract(const uint16_mat *mat, sint32 start_row, sint32 start_col, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_extract(const uint32_mat *mat, sint32 start_row, sint32 start_col, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_extract(const sint16_mat *mat, sint32 start_row, sint32 start_col, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_extract(const sint32_mat *mat, sint32 start_row, sint32 start_col, void *intrp_ptr, sint32_mat *result);#
 * \item #MUTIL_LIBEXPORT mutil_errcode matcpx_extract(const dcomplex_mat *mat, sint32 start_row, sint32 start_col, void *intrp_ptr, dcomplex_mat *result);#
 *  \end{itemize}
 * @see matuniv_assign
 * @see matuniv_assign_submat
 * @see matuniv_assign_zeropad
 * @see matuniv_translate
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_extract( const univ_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, univ_mat *result);


/* This function is documented under matuniv_extract in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_extract( const double_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_extract in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matflt_extract( const float_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_extract in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats32_extract( const sint32_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_extract in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats16_extract( const sint16_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_extract in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu32_extract( const uint32_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, uint32_mat *result );


/* This function is documented under matuniv_extract in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu16_extract( const uint16_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_extract in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu8_extract( const uint8_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, uint8_mat *result );

/* This function is documented under matuniv_extract in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_extract( const dcomplex_mat *mat,
    sint32 start_row, sint32 start_col, void *intrp_ptr, dcomplex_mat *result );

/** Find the transpose of a matrix.
 * Find the transpose of a matrix, and put the result
 * into a previously allocated matrix of the same type with transposed
 * dimensions.  Except for single-row and single-column matrices,
 * the output matrix must not have the same data location as the
 * input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_assn.h
 * @source mat\_assn.c
 * @library matrix
 * @usage  #err_code = matuniv_transpose(&mat, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to input matrix.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to the transposed matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_transpose(const double_mat *mat, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_transpose(const float_mat *mat, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_transpose(const uint8_mat *mat, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_transpose(const uint16_mat *mat, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_transpose(const uint32_mat *mat, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_transpose(const sint16_mat *mat, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_transpose(const sint32_mat *mat, void *intrp_ptr, sint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matcpx_transpose(const dcomplex_mat *cmat, void *intrp_ptr, dcomplex_mat *result);#
 *  \end{itemize}
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_transpose( const univ_mat *mat,
  void *intrp_ptr, univ_mat *result);


/* This function is documented under matuniv_transpose in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_transpose( const double_mat *mat,
  void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_transpose in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matflt_transpose( const float_mat *mat,
  void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_transpose in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats32_transpose( const sint32_mat *mat,
  void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_transpose in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats16_transpose( const sint16_mat *mat,
  void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_transpose in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu32_transpose( const uint32_mat *mat,
  void *intrp_ptr, uint32_mat *result );


/* This function is documented under matuniv_transpose in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu16_transpose( const uint16_mat *mat,
  void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_transpose in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu8_transpose( const uint8_mat *mat,
  void *intrp_ptr, uint8_mat *result );


/* This function is documented under matuniv_transpose in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_transpose( const dcomplex_mat *cmat,
  void *intrp_ptr, dcomplex_mat *result );


/** Flip a matrix left/right.
 * Flip the columns of matrix from left to right, and put the result
 * into a previously allocated matrix of the same size and type.
 * The output matrix may have the same data location as the
 * input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_assn.h
 * @source mat\_assn.c
 * @library matrix
 * @usage  #err_code = matuniv_flip_left_right(&mat, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to input matrix.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to the column flipped matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_flip_left_right(const double_mat *mat, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_flip_left_right(const float_mat *mat, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_flip_left_right(const uint8_mat *mat, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_flip_left_right(const uint16_mat *mat, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_flip_left_right(const uint32_mat *mat, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_flip_left_right(const sint16_mat *mat, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_flip_left_right(const sint32_mat *mat, void *intrp_ptr, sint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matcpx_flip_left_right(const dcomplex_mat *cmat, void *intrp_ptr, dcomplex_mat *result);#
 *  \end{itemize}
 * @see matuniv_flip_up_down
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_flip_left_right( const univ_mat *mat,
  void *intrp_ptr, univ_mat *result);


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_flip_left_right( const double_mat *mat,
  void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matflt_flip_left_right( const float_mat *mat,
  void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats32_flip_left_right( const sint32_mat *mat,
  void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats16_flip_left_right( const sint16_mat *mat,
  void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu32_flip_left_right( const uint32_mat *mat,
  void *intrp_ptr, uint32_mat *result );


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu16_flip_left_right( const uint16_mat *mat,
  void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu8_flip_left_right( const uint8_mat *mat,
  void *intrp_ptr, uint8_mat *result );


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_flip_left_right( const dcomplex_mat *cmat,
  void *intrp_ptr, dcomplex_mat *result );


/** Flip a matrix up/down.
 * Flip the rows of matrix from top to bottom, and put the result
 * into a previously allocated matrix of the same size and type.
 * The output matrix may have the same data location as the
 * input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_assn.h
 * @source mat\_assn.c
 * @library matrix
 * @usage  #err_code = matuniv_flip_up_down(&mat, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to input matrix.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to the column flipped matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_flip_up_down(const double_mat *mat, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_flip_up_down(const float_mat *mat, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_flip_up_down(const uint8_mat *mat, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_flip_up_down(const uint16_mat *mat, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_flip_up_down(const uint32_mat *mat, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_flip_up_down(const sint16_mat *mat, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_flip_up_down(const sint32_mat *mat, void *intrp_ptr, sint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matcpx_flip_up_down(const dcomplex_mat *cmat, void *intrp_ptr, dcomplex_mat *result);#
 *  \end{itemize}
 * @see matuniv_flip_left_right
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_flip_up_down( const univ_mat *mat,
  void *intrp_ptr, univ_mat *result);


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_flip_up_down( const double_mat *mat,
  void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matflt_flip_up_down( const float_mat *mat,
  void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats32_flip_up_down( const sint32_mat *mat,
  void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats16_flip_up_down( const sint16_mat *mat,
  void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu32_flip_up_down( const uint32_mat *mat,
  void *intrp_ptr, uint32_mat *result );


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu16_flip_up_down( const uint16_mat *mat,
  void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu8_flip_up_down( const uint8_mat *mat,
  void *intrp_ptr, uint8_mat *result );


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_flip_up_down( const dcomplex_mat *cmat,
  void *intrp_ptr, dcomplex_mat *result );


/** Translate a matrix.
 * Translate a matrix by a provided number of rows and
 * columns, and set the element positions previously occupied
 * by the translated elements to a value provided by the user.
 * Store the result in a previously allocated matrix of the same
 * dimensions and type as the input matrix (some matrix values will
 * be lost).
 * The output matrix must not have the same data location as
 * the input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_assn.h
 * @source mat\_assn.c
 * @library matrix
 * @usage  #err_code = matuniv_translate(&mat, row_shift, col_shift, pad_value, intrp_ptr, &result);#
 * @return          Standard mutils error/OK code.
 * @param mat       Pointer to the input matrix.
 * @param row_shift The number of rows to shift.
 * @param col_shift The number of columns to shift.
 * @param pad_value Value to assign to the elements previously occupied by the
 *    translated elements.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param result    Pointer to the resulting matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_translate(const double_mat *mat, sint32 row_shift, sint32 col_shift, double pad_value, void *intrp_ptr, double_mat *result);#
 * \item #MUTIL_LIBEXPORT mutil_errcode matflt_translate( const float_mat *mat, sint32 row_shift, sint32 col_shift, float pad_value, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_translate(const uint8_mat *mat, sint32 row_shift, sint32 col_shift, uint8 pad_value, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_translate(const uint16_mat *mat, sint32 row_shift, sint32 col_shift, uint16 pad_value, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_translate(const uint32_mat *mat, sint32 row_shift, sint32 col_shift, uint32 pad_value, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_translate(const sint16_mat *mat, sint32 row_shift, sint32 col_shift, sint16 pad_value, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_translate(const sint32_mat *mat, sint32 row_shift, sint32 col_shift, sint32 pad_value, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_assign_zeropad
 * @see matuniv_extract
 * @see _univ_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_translate( const univ_mat *mat,
  sint32 row_shift, sint32 col_shift, univ_scalar pad_value,
  void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_translate in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_translate( const double_mat *mat,
  sint32 row_shift, sint32 col_shift, double pad_value,
  void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_translate in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matflt_translate( const float_mat *mat,
  sint32 row_shift, sint32 col_shift, float pad_value,
  void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_translate in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats32_translate( const sint32_mat *mat,
  sint32 row_shift, sint32 col_shift, sint32 pad_value,
  void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_translate in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode mats16_translate( const sint16_mat *mat,
  sint32 row_shift, sint32 col_shift, sint16 pad_value,
  void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_translate in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu32_translate( const uint32_mat *mat,
  sint32 row_shift, sint32 col_shift, uint32 pad_value,
  void *intrp_ptr, uint32_mat *result );


/* This function is documented under matuniv_translate in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu16_translate( const uint16_mat *mat,
  sint32 row_shift, sint32 col_shift, uint16 pad_value,
  void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_translate in mat_assn.h */
MUTIL_LIBEXPORT mutil_errcode matu8_translate( const uint8_mat *mat,
  sint32 row_shift, sint32 col_shift, uint8 pad_value,
  void *intrp_ptr, uint8_mat *result );


#ifdef __cplusplus
}
#endif

#endif /* ifdef IN_MAT_ASSN_H */
