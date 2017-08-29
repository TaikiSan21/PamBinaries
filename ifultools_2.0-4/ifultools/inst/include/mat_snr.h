
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_snr.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_SNR_H_
#define IN_MAT_SNR_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

/* This file contains function declarations for matrix characterization
   related to signal to noise ratio and entropy.
*/

#ifdef __cplusplus
extern "C" {
#endif

/*
***************************
 Universal matrix functions
***************************
 */


/** Mean squared error.
 * This function calculates the mean squared sum of the
 * pixel-wise difference between two matrices.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_snr.h
 * @source mat\_snr.c
 * @library matrix
 * @usage #errcode = matuniv_mean_squared_error(&matrix1, &matrix2, &intrp_ptr, &result);#
 * @return   Standard mutils error/OK code.
 * @param  matrix1    Pointer to first input matrix.
 * @param  matrix2    Pointer to second input matrix.
 * @param  intrp_ptr  Pointer for implementation of interrupt checking.
 * @param  result     Pointer for returning result.
 * @see matuniv_signal_to_noise_ratio
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_mean_squared_error(
    const univ_mat *matrix1,
    const univ_mat *matrix2,
    void           *intrp_ptr,
    double         *result );


/** Signal-to-noise ratio.
 * This function takes two matrices, one representing a noiseless
 * matrix, the other representing the same matrix corrupted by noise, and
 * returns the signal-to-noise ratio of the corrupted matrix.
 * It also returns a flag to indicate infinite signal-to-noise ratio --
 * if the noise power is zero, the result
 * is set to MUTIL\_DOUBLE\_MAX, and the flag is set to TRUE.  The flag,
 * and not the result value, should be used to test whether the
 * signal-to-noise ratio was infinite.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_snr.h
 * @source mat\_snr.c
 * @library matrix
 * @usage #errcode = matuniv_signal_to_noise_ratio(&orig_mat, &noisy_mat, &intrp_ptr, &result, &inf_flag);#
 * @return    Standard mutils error/OK code.
 * @param  orig_mat   Pointer to the noiseless matrix.
 * @param  noisy_mat  Pointer to the noisy matrix.
 * @param  intrp_ptr  Pointer for implementation of interrupt checking.
 * @param  result     Pointer for returning result.
 * @param  inf_flag   Pointer for returning flag indicating infinite SNR.
 * @see matuniv_mean_squared_error
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_signal_to_noise_ratio(
    const univ_mat *orig_mat,
    const univ_mat *noisy_mat,
    void           *intrp_ptr,
    double         *result,
    boolean        *inf_flag );


/** Entropy calculation.
 * The entropy, or average value of self information, is defined as
 * $-\sum_{i=0}^{N-1} P[i] \log_{2} P[i]$.
 * P[i] is the probability distribution of source symbols from
 * a discrete memoryless source. The entropy represents the
 * minimum number of binary digits (bits) that can be used, on average,
 * to represent the source symbols.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_snr.h
 * @source mat\_snr.c
 * @library matrix
 * @usage  #err_code = matuniv_entropy(&frequency, intrp_ptr, &entropy);#
 * @return   Standard mutils error/OK code.
 * @param frequency   Pointer to the histogram or frequency.  If the
 *     data is of integer type, then frequency is treated as a histogram
 *     where P[i] is computed as histogram[i]/sum(histogram).
 *     If the data is of floating point type, then frequency is
 *     is treated as a set of probabilities P[i] that must sum to 1.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param entropy     Pointer to type double result.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_entropy(const double_mat *frequency, void *intrp_ptr, double *entropy);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_entropy(const double_mat *frequency, void *intrp_ptr, double *entropy);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_entropy(const uint8_mat *frequency, void *intrp_ptr, double *entropy);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_entropy(const uint16_mat *frequency, void *intrp_ptr, double *entropy);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_entropy(const uint32_mat *frequency, void *intrp_ptr, double *entropy);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_entropy(const sint16_mat *frequency, void *intrp_ptr, double_mat *entropy);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_entropy(const sint32_mat *frequency, void *intrp_ptr, double *entropy);#
 *  \end{itemize}
 * @see Matrix Data Types
 * @see Interrupt Handling
*/
MUTIL_LIBEXPORT mutil_errcode matuniv_entropy(
    const univ_mat *frequency,
    void           *intrp_ptr,
    double         *entropy );


/* This function is documented under matuniv_entropy in mat_snr.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_entropy(
    const double_mat *frequency,
    void             *intrp_ptr,
    double           *entropy );


/* This function is documented under matuniv_entropy in mat_snr.h */
MUTIL_LIBEXPORT mutil_errcode matflt_entropy(
    const float_mat *frequency,
    void            *intrp_ptr,
    double          *entropy );


/* This function is documented under matuniv_entropy in mat_snr.h */
MUTIL_LIBEXPORT mutil_errcode mats32_entropy(
    const sint32_mat *frequency,
    void             *intrp_ptr,
    double           *entropy );


/* This function is documented under matuniv_entropy in mat_snr.h */
MUTIL_LIBEXPORT mutil_errcode mats16_entropy(
    const sint16_mat *frequency,
    void             *intrp_ptr,
    double           *entropy );


/* This function is documented under matuniv_entropy in mat_snr.h */
MUTIL_LIBEXPORT mutil_errcode matu32_entropy(
    const uint32_mat *frequency,
    void             *intrp_ptr,
    double           *entropy );


/* This function is documented under matuniv_entropy in mat_snr.h */
MUTIL_LIBEXPORT mutil_errcode matu16_entropy(
    const uint16_mat *frequency,
    void             *intrp_ptr,
    double           *entropy );


/* This function is documented under matuniv_entropy in mat_snr.h */
MUTIL_LIBEXPORT mutil_errcode matu8_entropy(
    const uint8_mat *frequency,
    void            *intrp_ptr,
    double          *entropy );


#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_SNR_H_*/
