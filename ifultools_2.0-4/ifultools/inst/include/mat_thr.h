
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_thr.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_THR_H_
#define IN_MAT_THR_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

/* This file contains function declarations  for the
   basic matrix thresholding operation for the mutils library.
 */

#ifdef __cplusplus
extern "C" {
#endif


/** Thresholding.
 * Takes as input a matrix, two threshold values, and two
 * output values.  Returns a matrix of the same size, whose
 * elements assume either of two values, depending on whether or not the
 * original matrix element values lie between the two
 * thresholds (values equal to either threshold are considered to be
 * between the thresholds).  Depending on the inputs, implements classical
 * thresholding (binarize) or may be used to
 * emphasize/deemphasize values within a certain range. It is also
 * possible to assign the same value to the entire matrix.
 * The input and output data may point to the same location in memory.
 * If the data types of in\_data and results are the same, the specific
 * subtype function (e.g. matdbl\_threshold) will be called,
 * resulting in a faster implementation.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_thr.h
 * @source mat\_thr.c
 * @library matrix
 * @usage  #err_code = matuniv_threshold(&in_data, low_t, high_t, out_val, in_val, intrp_ptr, &result);#
 * @return     Standard mutils error/OK code.
 * @param in_data       Pointer to the input matrix.
 * @param low_t     Lower threshold.
 * @param high_t    Upper threshold.
 * @param out_val   New value for all values not between the two thresholds.
 * @param in_val    New value for all values between the two thresholds.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param result    Pointer to a previously allocated matrix of the same
 *    size for output.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_threshold(const double_mat *in_data, double low_t, double high_t, double out_val, double in_val, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_threshold(const float_mat *in_data, double low_t, double high_t, float out_val, float in_val, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_threshold(const uint8_mat *in_data, double low_t, double high_t, uint8 out_val, uint8 in_val, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_threshold(const uint16_mat *in_data, double low_t, double high_t, uint16 out_val, uint16 in_val, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_threshold(const uint32_mat *in_data, double low_t, double high_t, uint32 out_val, uint32 in_val, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_threshold(const sint16_mat *in_data, double low_t, double high_t, sint16 out_val, sint16 in_val, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_threshold(const sint32_mat *in_data, double low_t, double high_t, sint32 out_val, sint32 in_val, void *intrp_ptr, sint32_mat *result);#
 * \end{itemize}
 * @see clsuniv_scalar_quantize_linear
 * @see clsuniv_scalar_quantize_angle
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_threshold( const univ_mat *in_data,
  double low_t, double high_t, univ_scalar out_val, univ_scalar	in_val,
  void *intrp_ptr, univ_mat *result );

/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matdbl_threshold( const double_mat *in_data,
  double low_t, double high_t, double out_val, double in_val,
  void  *intrp_ptr, double_mat *result );

/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matflt_threshold( const float_mat *in_data,
  double low_t, double high_t, float out_val, float in_val,
  void  *intrp_ptr, float_mat *result );

/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu8_threshold( const uint8_mat *in_data,
  double low_t, double high_t, uint8 out_val, uint8 in_val,
  void  *intrp_ptr, uint8_mat *result );

/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu16_threshold( const uint16_mat *in_data,
  double low_t, double high_t, uint16 out_val, uint16 in_val,
  void  *intrp_ptr, uint16_mat *result );

/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu32_threshold( const uint32_mat *in_data,
  double low_t, double high_t, uint32 out_val, uint32 in_val,
  void  *intrp_ptr, uint32_mat *result );

/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats16_threshold( const sint16_mat *in_data,
  double low_t, double high_t, sint16 out_val, sint16 in_val,
  void  *intrp_ptr, sint16_mat *result );

/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats32_threshold( const sint32_mat *in_data,
  double low_t, double high_t, sint32 out_val, sint32 in_val,
  void  *intrp_ptr, sint32_mat *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_THR_H_*/
