
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_stat.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_STAT_H_
#define IN_MAT_STAT_H_

#include "mat_type.h"

#include "ut_plat.h"
#include "ut_type.h"
#include "ut_err.h"
#include "ut_math.h"

/* This file contains function declarations for
   statistical functions on matrices for the mutils library.
*/

#ifdef __cplusplus
extern "C" {
#endif


/*
***************************
 Universal matrix functions
***************************
 */


/** Median of a matrix.
 * Find the median of the values in a matrix, using a partial
 * quicksort to find the middle value(s), and then averaging them
 * if the number of elements in the matrix is even.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_stat.h
 * @source mat\_stat.c
 * @library matrix
 * @usage  #err_code = matuniv_median(&mymat, intrp_ptr, &out);#
 * @return        Standard mutils error/OK code.
 * @param mat       Pointer to the matrix to find median of.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param median    Pointer to returned median value.
 * @see matuniv_sort_index_partial
 * @see matuniv_quantiles
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_median( const univ_mat *mat,
  void *intrp_ptr, double *median );


/** Mean and variance of a matrix.
 * Compute the mean and variance of the values of a matrix containing
 * non-complex data. The parameter {\tt unbiased} controls the method of
 * estimation for the variance. If TRUE, the variance is computed as
 * $var = \frac{\sum ( x - \mbox{mean}(x) )^{2}}{(N-1)}$, which is unbiased
 * if the matrix
 * elements (x) are obtained by simple random sampling. Otherwise, the
 * variance is computed as $var = \frac{\sum ( x - \mbox{mean}(x) )^{2}}{N}$.
 * The mean and variance are returned in separate variables of type double.
 *
 * This function uses the standard mutils casting conventions
 * (see \Ref{Library Basics}).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_stat.h
 * @source mat\_stat.c
 * @library matrix
 * @usage  #err_code = matuniv_mean_variance(&mymat, unbiased, intrp_ptr, &mean, &variance);#
 * @return        Standard mutils error/OK code.
 * @param mat       Pointer to the matrix for which the mean and variance
 *    are computed.
 * @param unbiased  TRUE for unbiased variance estimate, FALSE otherwise.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param mean      Pointer to return value for the mean.
 *   (NULL to omit this output).
 * @param variance  Pointer to return value for the variance.
 *   (NULL to omit this output).
 * @see matuniv_median
 * @see matuniv_quantiles
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_mean_variance( const univ_mat *mat,
  boolean unbiased, void *intrp_ptr, double *mean, double *variance );


/** Quantiles of a matrix.
 * Find quantile values of a matrix with given probabilities
 * (which must be between 0 and 1), using a partial
 * quicksort to find the values in the matrix with the appropriate
 * sorted ranks.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_stat.h
 * @source mat\_stat.c
 * @library matrix
 * @usage  #err_code = matuniv_quantiles(&mymat, &probs, intrp_ptr, &out);#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix to find quantiles of.
 * @param probs      Pointer to matrix of probabilities to calculate
 *    quantiles for (treated as flat array, of type double).
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param quantiles  Pointer for returned quantiles (same dimensions as
 *    probs, of type double).
 * @see matuniv_sort_index_partial
 * @see matuniv_median
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_quantiles( const univ_mat *mat,
  const univ_mat *probs, void *intrp_ptr, univ_mat *quantiles );


/** Histogram of values in a matrix set.
 * Count the number of values in a matrix set that lie in
 * each of the equally-spaced histogram bins between given
 * starting and ending values.  Values that lie on a bin boundary
 * will be counted as belonging to the lower bin, and the starting
 * value itself will not be counted unless the include\_end parameter
 * is set to TRUE.
 *
 * Matrices of the matrix set must all be of the same type and size.
 * The output histogram must be a 1-column universal matrix
 * of type MUTIL\_UINT32.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_stat.h
 * @source mat\_stat.c
 * @library matrix
 * @usage  #err_code = matset_histogram(&matset, minval, maxval, TRUE, intrp_ptr, &result);#
 * @return    Standard mutils error/OK code.
 * @param matset      Pointer to the input matrix set.
 * @param start_val   Low end of histogram range.
 * @param end_val     High end of histogram range.
 * @param include_end If TRUE, count start\_val values in lowermost bin.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param histogram   Pointer to previously-allocated 1-column vector for
 *    result, of type MUTIL\_UINT32, as long as the number of desired bins.
 * @see matuniv_histogram
 * @see matuniv_range
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matset_histogram(
  const mat_set *matset, const double start_val, const double end_val,
  const boolean include_end, void *intrp_ptr, univ_mat *histogram );


/** Histogram of values in a matrix.
 * Count the number of values in a matrix that lie in
 * each of the equally-spaced histogram bins between given
 * starting and ending values.  Values that lie on a bin boundary
 * will be counted as belonging to the lower bin, and the starting
 * value itself will not be counted unless the include\_end parameter
 * is set to TRUE.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_stat.h
 * @source mat\_stat.c
 * @library matrix
 * @usage  #err_code = matuniv_histogram(&mat, minval, maxval, TRUE, intrp_ptr, &result);#
 * @return    Standard mutils error/OK code.
 * @param mat         Pointer to the input matrix.
 * @param start_val   Low end of histogram range.
 * @param end_val     High end of histogram range.
 * @param include_end If TRUE, count start\_val values in lowermost bin.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param histogram   Pointer to previously-allocated 1-column vector for
 *    result, of type MUTIL\_UINT32, as long as the number of desired bins.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_histogram(const double_mat *mat, double start_val, double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_histogram(const float_mat *mat, double start_val, double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_histogram(const uint8_mat *mat, double start_val, double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_histogram(const uint16_mat *mat, double start_val, double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_histogram(const uint32_mat *mat, double start_val, double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_histogram(const sint16_mat *mat, double start_val, double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_histogram(const sint32_mat *mat, double start_val, double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram);#
 * \end{itemize}
 * @see matset_histogram
 * @see matuniv_range
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_histogram( const univ_mat *mat,
  double start_val, double end_val, boolean include_end,
  void *intrp_ptr, univ_mat *histogram );


/* This function is documented under matuniv_histogram in mat_stat.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_histogram( const double_mat *mat,
  double start_val, double end_val, boolean include_end,
  void *intrp_ptr, uint32_mat *histogram );


/* This function is documented under matuniv_histogram in mat_stat.h */
MUTIL_LIBEXPORT mutil_errcode matflt_histogram( const float_mat *mat,
  double start_val, double end_val, boolean include_end,
  void *intrp_ptr, uint32_mat *histogram );


/* This function is documented under matuniv_histogram in mat_stat.h */
MUTIL_LIBEXPORT mutil_errcode matu8_histogram( const uint8_mat *mat,
  double start_val, double end_val, boolean include_end,
  void *intrp_ptr, uint32_mat *histogram );


/* This function is documented under matuniv_histogram in mat_stat.h */
MUTIL_LIBEXPORT mutil_errcode matu16_histogram( const uint16_mat *mat,
  double start_val, double end_val, boolean include_end,
  void *intrp_ptr, uint32_mat *histogram );


/* This function is documented under matuniv_histogram in mat_stat.h */
MUTIL_LIBEXPORT mutil_errcode matu32_histogram( const uint32_mat *mat,
  double start_val, double end_val, boolean include_end,
  void *intrp_ptr, uint32_mat *histogram );


/* This function is documented under matuniv_histogram in mat_stat.h */
MUTIL_LIBEXPORT mutil_errcode mats16_histogram( const sint16_mat *mat,
  double start_val, double end_val, boolean include_end,
  void *intrp_ptr, uint32_mat *histogram );


/* This function is documented under matuniv_histogram in mat_stat.h */
MUTIL_LIBEXPORT mutil_errcode mats32_histogram( const sint32_mat *mat,
  double start_val, double end_val, boolean include_end,
  void *intrp_ptr, uint32_mat *histogram );


#ifdef __cplusplus
}
#endif

#endif /* #ifndef IN_MAT_STAT_H_ */
