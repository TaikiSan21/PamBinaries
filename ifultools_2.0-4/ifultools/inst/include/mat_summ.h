
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_summ.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_SUMM_H
#define IN_MAT_SUMM_H

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "mat_any.h"

/* This file contains functions and macros for universal matrix
   summary operations.
 */

#ifdef __cplusplus
extern "C" {
#endif


/** Sum the elements of a matrix.
 * Take the sum of all the elements in a matrix, and output
 * the results into a \Ref{_univ_scalar} argument (or a
 * regular scalar for specific-type functions).
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data type
 * is capable of holding the result.  If you are unsure about
 * the range of your input data, it is safest to use double data types
 * to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_summ.h
 * @source mat\_summ.c
 * @library matrix
 * @usage  #err_code = matuniv_sum(&mat, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix to sum.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to returning the result.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_sum(const double_mat *mat, void *intrp_ptr, double *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_sum(const float_mat *mat, void *intrp_ptr, float *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_sum(const uint8_mat *mat, void *intrp_ptr, uint8 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_sum(const uint16_mat *mat, void *intrp_ptr, uint16 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_sum(const uint32_mat *mat, void *intrp_ptr, uint32 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_sum(const sint16_mat *mat, void *intrp_ptr, sint16 *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_sum(const sint32_mat *mat, void *intrp_ptr, sint32 *result);#
 *  \end{itemize}
 * @see matuniv_sum_rows
 * @see matuniv_sum_cols
 * @see _univ_scalar
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_sum( const univ_mat *mat,
  void *intrp_ptr, univ_scalar *result );


/* This function is documented under matuniv_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_sum( const double_mat *mat, void *intrp_ptr,
  double *result );


/* This function is documented under matuniv_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matflt_sum( const float_mat *mat, void *intrp_ptr,
  float *result );


/* This function is documented under matuniv_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats32_sum( const sint32_mat *mat, void *intrp_ptr,
  sint32 *result );


/* This function is documented under matuniv_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats16_sum( const sint16_mat *mat, void *intrp_ptr,
  sint16 *result );


/* This function is documented under matuniv_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu32_sum( const uint32_mat *mat, void *intrp_ptr,
  uint32 *result );


/* This function is documented under matuniv_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu16_sum( const uint16_mat *mat, void *intrp_ptr,
  uint16 *result );


/* This function is documented under matuniv_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu8_sum( const uint8_mat *mat, void *intrp_ptr,
  uint8 *result );


/** Cumulative sum.
 * Computes the cumulative sum of a matrix
 * considering the matrix data as a single flat vector.
 * Stores the result in a previously allocated matrix of the same
 * size and type as the input matrix (overflow/underflow is
 * not checked).
 * The output matrix may have the same data location as
 * the input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_summ.h
 * @source mat\_summ.c
 * @library matrix
 * @usage  #err_code = matuniv_cumulative_sum(&in_mat, intrp_ptr, &result_mat);#
 * @return          Standard mutils error/OK code.
 * @param mat_in   Pointer to the input matrix.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param mat_out  Pointer to the resulting matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_cumulative_sum(const double_mat *mat_in, void *intrp_ptr, double_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_cumulative_sum(const float_mat *mat_in, void *intrp_ptr, float_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_cumulative_sum(const uint8_mat *mat_in, void *intrp_ptr, uint8_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_cumulative_sum(const uint16_mat *mat_in, void *intrp_ptr, uint16_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_cumulative_sum(const uint32_mat *mat_in, void *intrp_ptr, uint32_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_cumulative_sum(const sint16_mat *mat_in, void *intrp_ptr, sint16_mat *mat_out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_cumulative_sum(const sint32_mat *mat_in, void *intrp_ptr, sint32_mat *mat_out);#
 *  \end{itemize}
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_cumulative_sum(
    const univ_mat *mat_in,
    void           *intrp_ptr,
    univ_mat       *mat_out );


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_cumulative_sum( const double_mat *mat_in,
  void *intrp_ptr, double_mat *mat_out );


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matflt_cumulative_sum( const float_mat *mat_in,
  void *intrp_ptr, float_mat *mat_out );


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats32_cumulative_sum( const sint32_mat *mat_in,
  void *intrp_ptr, sint32_mat *mat_out );


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats16_cumulative_sum( const sint16_mat *mat_in,
  void *intrp_ptr, sint16_mat *mat_out );


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu32_cumulative_sum( const uint32_mat *mat_in,
  void *intrp_ptr, uint32_mat *mat_out );


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu16_cumulative_sum( const uint16_mat *mat_in,
  void *intrp_ptr, uint16_mat *mat_out );


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu8_cumulative_sum( const uint8_mat *mat_in,
  void *intrp_ptr, uint8_mat *mat_out );


/** Sum the rows of a matrix.
 * Take the sum of each row of a matrix, and put
 * the results into a previously allocated matrix
 * with the same number of rows and a single column.
 * If the input matrix has only one column, it may be the same
 * matrix as the output matrix.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data type
 * is capable of holding the result.  If you are unsure about
 * the range of your input data, it is safest to use double data types
 * to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_summ.h
 * @source mat\_summ.c
 * @library matrix
 * @usage  #err_code = matuniv_sum_rows(&mat, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix to sum.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to the row sum matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_sum_rows(const double_mat *mat, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_sum_rows(const float_mat *mat, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_sum_rows(const uint8_mat *mat, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_sum_rows(const uint16_mat *mat, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_sum_rows(const uint32_mat *mat, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_sum_rows(const sint16_mat *mat, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_sum_rows(const sint32_mat *mat, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_sum_cols
 * @see matuniv_sum
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_sum_rows( const univ_mat *mat,
  void *intrp_ptr, univ_mat *result);


/* This function is documented under matuniv_sum_rows in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_sum_rows( const double_mat *mat,
  void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_sum_rows in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matflt_sum_rows( const float_mat *mat,
  void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_sum_rows in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats32_sum_rows( const sint32_mat *mat,
  void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_sum_rows in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats16_sum_rows( const sint16_mat *mat,
  void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_sum_rows in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu32_sum_rows( const uint32_mat *mat,
  void *intrp_ptr, uint32_mat *result );


/* This function is documented under matuniv_sum_rows in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu16_sum_rows( const uint16_mat *mat,
  void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_sum_rows in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu8_sum_rows( const uint8_mat *mat,
  void *intrp_ptr, uint8_mat *result );


/** Sum the columns of a matrix.
 * Take the sum of each column of a matrix, and put
 * the results into a previously allocated matrix
 * with the same number of columns and a single row.
 * If the input matrix has only one row, it may be the same
 * matrix as the output matrix.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * This function assumes, but does not verify, that the output data type
 * is capable of holding the result.  If you are unsure about
 * the range of your input data, it is safest to use double data types
 * to avoid overflow.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_summ.h
 * @source mat\_summ.c
 * @library matrix
 * @usage  #err_code = matuniv_sum_cols(&mat, intrp_ptr, &result);#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix to sum.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to the column sum matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_sum_cols(const double_mat *mat, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_sum_cols(const float_mat *mat, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_sum_cols(const uint8_mat *mat, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_sum_cols(const uint16_mat *mat, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_sum_cols(const uint32_mat *mat, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_sum_cols(const sint16_mat *mat, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_sum_cols(const sint32_mat *mat, void *intrp_ptr, sint32_mat *result);#
 *  \end{itemize}
 * @see matuniv_sum_rows
 * @see matuniv_sum
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_sum_cols( const univ_mat *mat,
  void *intrp_ptr, univ_mat *result );


/* This function is documented under matuniv_sum_cols in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_sum_cols( const double_mat *mat,
  void *intrp_ptr, double_mat *result );


/* This function is documented under matuniv_sum_cols in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matflt_sum_cols( const float_mat *mat,
  void *intrp_ptr, float_mat *result );


/* This function is documented under matuniv_sum_cols in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats32_sum_cols( const sint32_mat *mat,
  void *intrp_ptr, sint32_mat *result );


/* This function is documented under matuniv_sum_cols in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats16_sum_cols( const sint16_mat *mat,
  void *intrp_ptr, sint16_mat *result );


/* This function is documented under matuniv_sum_cols in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu32_sum_cols( const uint32_mat *mat,
  void *intrp_ptr, uint32_mat *result );


/* This function is documented under matuniv_sum_cols in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu16_sum_cols( const uint16_mat *mat,
  void *intrp_ptr, uint16_mat *result );


/* This function is documented under matuniv_sum_cols in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu8_sum_cols( const uint8_mat *mat,
  void *intrp_ptr, uint8_mat *result );


/** Find the range of values in a matrix.
 * Find the minimum and maximum values in a matrix.
 * Complex matrices are not supported, and the value is always
 * returned as the same type as the input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_summ.h
 * @source mat\_summ.c
 * @library matrix
 * @usage  #err_code = matuniv_range(&mat, intrp_ptr, &min, &max);#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix to find range of.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param minval     Pointer to return the minimum value.
 * @param maxval     Pointer to return the maximum value.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_range(const double_mat *mat, void *intrp_ptr, double *minval, double *maxval);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_range(const float_mat *mat, void *intrp_ptr, float *minval, float *maxval);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_range(const uint8_mat *mat, void *intrp_ptr, uint8 *minval, uint8 *maxval);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_range(const uint16_mat *mat, void *intrp_ptr, uint16 *minval, uint16 *maxval);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_range(const uint32_mat *mat, void *intrp_ptr, uint32 *minval, uint32 *maxval);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_range(const sint16_mat *mat, void *intrp_ptr, sint16 *minval, sint16 *maxval);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_range(const sint32_mat *mat, void *intrp_ptr, sint32 *minval, sint32 *maxval);#
 *  \end{itemize}
 * @see matuniv_range_robust
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_range( const univ_mat *mat,
  void *intrp_ptr, univ_scalar *minval, univ_scalar *maxval );


/* This function is documented under matuniv_range in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_range( const double_mat *mat,
  void *intrp_ptr, double *minval, double *maxval );


/* This function is documented under matuniv_range in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matflt_range( const float_mat *mat,
  void *intrp_ptr, float *minval, float *maxval );


/* This function is documented under matuniv_range in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats32_range( const sint32_mat *mat,
  void *intrp_ptr, sint32 *minval, sint32 *maxval );


/* This function is documented under matuniv_range in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats16_range( const sint16_mat *mat,
  void *intrp_ptr, sint16 *minval, sint16 *maxval );


/* This function is documented under matuniv_range in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu32_range( const uint32_mat *mat,
  void *intrp_ptr, uint32 *minval, uint32 *maxval );


/* This function is documented under matuniv_range in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu16_range( const uint16_mat *mat,
  void *intrp_ptr, uint16 *minval, uint16 *maxval );


/* This function is documented under matuniv_range in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu8_range( const uint8_mat *mat,
  void *intrp_ptr, uint8 *minval, uint8 *maxval );


/** Find a robust range of values in a matrix.
 * Find the minimum and maximum values in a matrix,
 * excluding values that lie outside a given range of values.
 * Typically, this would be used to exclude the top and bottom
 * few percent of the range and see whether the matrix still
 * spanned the entire range.  If there are no values within the
 * allowed range, an error value is returned.
 *
 * Complex matrices are not supported, and the value is always
 * returned as the same type as the input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_summ.h
 * @source mat\_summ.c
 * @library matrix
 * @usage  #err_code = matuniv_range_robust(&mat, excl1, excl2, intrp_ptr, &min, &max);#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix to find range of.
 * @param lowexc     Exclude values below this value.
 * @param highexc    Exclude values above this value.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param minval     Pointer to return the minimum value.
 * @param maxval     Pointer to return the maximum value.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_range_robust(const double_mat *mat, double lowexc, double highexc, void *intrp_ptr, double *minval, double *maxval);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_range_robust(const float_mat *mat, float lowexc, float highexc, void *intrp_ptr, float *minval, float *maxval);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_range_robust(const uint8_mat *mat, uint8 lowexc, uint8 highexc, void *intrp_ptr, uint8 *minval, uint8 *maxval);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_range_robust(const uint16_mat *mat, uint16 lowexc, uint16 highexc, void *intrp_ptr, uint16 *minval, uint16 *maxval);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_range_robust(const uint32_mat *mat, uint32 lowexc, uint32 highexc, void *intrp_ptr, uint32 *minval, uint32 *maxval);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_range_robust(const sint16_mat *mat, sint16 lowexc, sint16 highexc, void *intrp_ptr, sint16 *minval, sint16 *maxval);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_range_robust(const sint32_mat *mat, sint32 lowexc, sint32 highexc, void *intrp_ptr, sint32 *minval, sint32 *maxval);#
 *  \end{itemize}
 * @see matuniv_range
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_range_robust( const univ_mat *mat,
  double lowexc, double highexc, void *intrp_ptr,
  univ_scalar *minval, univ_scalar *maxval );


/* This function is documented under matuniv_range_robust in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_range_robust( const double_mat *mat,
  double lowexc, double highexc, void *intrp_ptr,
  double *minval, double *maxval );


/* This function is documented under matuniv_range_robust in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matflt_range_robust( const float_mat *mat,
  float lowexc, float highexc, void *intrp_ptr,
  float *minval, float *maxval );


/* This function is documented under matuniv_range_robust in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats32_range_robust( const sint32_mat *mat,
  sint32 lowexc, sint32 highexc, void *intrp_ptr,
  sint32 *minval, sint32 *maxval );


/* This function is documented under matuniv_range_robust in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode mats16_range_robust( const sint16_mat *mat,
  sint16 lowexc, sint16 highexc, void *intrp_ptr,
  sint16 *minval, sint16 *maxval );


/* This function is documented under matuniv_range_robust in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu32_range_robust( const uint32_mat *mat,
  uint32 lowexc, uint32 highexc, void *intrp_ptr,
  uint32 *minval, uint32 *maxval );


/* This function is documented under matuniv_range_robust in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu16_range_robust( const uint16_mat *mat,
  uint16 lowexc, uint16 highexc, void *intrp_ptr,
  uint16 *minval, uint16 *maxval );


/* This function is documented under matuniv_range_robust in mat_summ.h */
MUTIL_LIBEXPORT mutil_errcode matu8_range_robust( const uint8_mat *mat,
  uint8 lowexc, uint8 highexc, void *intrp_ptr,
  uint8 *minval, uint8 *maxval );


#ifdef __cplusplus
}
#endif

#endif /* ifdef IN_MAT_SUMM_H */
