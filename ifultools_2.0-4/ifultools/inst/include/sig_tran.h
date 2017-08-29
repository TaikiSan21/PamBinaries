
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/sig_tran.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_SIG_TRAN_H
#define IN_SIG_TRAN_H

#include "ut_plat.h"
#include "mat_type.h"
#include "ut_err.h"
#include "ut_type.h"

/* This file contains function declarations for signal processing
   involving transformations such as discrete cosine
   and Fourier transforms.
*/

#ifdef __cplusplus
extern "C" {
#endif

/*
**************************
Universal matrix Functions
**************************
*/

/** Discrete cosine transform of type II.
 * Computes the discrete cosine transform of type II for
 * a signal of any length. For 1D signals, the data should be
 * arranged in a column vector. If the input matrix is 2D, it
 * computes the transform of each column. This is useful for
 * multi-channel signals, with one channel per column.
 *
 * Uses a discrete Fourier transform algorithm to compute the transform,
 * exploiting the close relationship between the DFT and the DCT-II.
 * The algorithm calculates
 * the DFT of a sequence which is twice as long as the original one,
 * and zero padded at the end. Then to go from the DFT to the DCT, the
 * following formula is used for the first N elements of the DFT, where
 * N is the original signal length:
 *
 * $G(n) = \sqrt{\frac{2}{N}} (\cos(\alpha) \Re(G_{dft}(n)) +
 *          \sin(\alpha) \Im(G_{dft}(n)))$
 *
 * where $\alpha$ is a sequence of angles with values
 * $\frac{\pi}{2N}$
 *
 * Reference: K.R.Rao and P.Yip, ``Discrete Cosine Transform,''
 *            Academic Press, 1990.
 *
 * @limits Not all data types are supported.
 * @usage #err = siguniv_transform_discrete_cosine_II(&signal, &intrp_ptr, &result);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_tran.h
 * @source sig\_tran.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param sig Pointer to matrix containing data to be transformed.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param result Pointer to matrix that will contain the transform.
 *               Must be pre-allocated and of the same size and type as sig.
 * @same \begin{itemize}
 *    \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_transform_discrete_cosine_II(const double_mat *sig, void *intrp_ptr, double_mat *result);#
 * \end{itemize}
 * @see siguniv_transform_discrete_cosine_II_inverse
 * @see siguniv_transform_discrete_cosine_II_2d
 * @see siguniv_transform_discrete_cosine_II_pow2
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_transform_discrete_cosine_II(
  const univ_mat *sig,
  void *intrp_ptr,
  univ_mat *result );


/* Function documented with universal matrix function, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_transform_discrete_cosine_II(
  const double_mat *sig,
  void *intrp_ptr,
  double_mat *result );


/** Inverse discrete cosine transform of type II.
 * Computes the inverse discrete cosine transform of type II for
 * a signal of any length. For 1D signals, the data should be
 * arranged in a column vector. If the input matrix is 2D, it
 * computes the transform of each column. This is useful for
 * multi-channel signals, with one channel per column.
 *
 * Uses a discrete Fourier transform algorithm to compute the transform,
 * exploiting the close relationship between the DFT and the DCT-II.
 * The algorithm calculates
 * the DFT of a sequence which is four times as long as the original one,
 * and zero padded at the end. Then, to go from the DFT to the original
 * signal, the following formula is used for the first N elements of
 * the DFT, where N is the original signal length:
 *
 * $g(n) = \sqrt{\frac{2}{N}} \Re(G_{dft}(2n+1))$
 *
 * Reference: K.R.Rao and P.Yip, ``Discrete Cosine Transform,''
 *            Academic Press, 1990.
 * @limits Not all data types are supported.
 * @usage #err = siguniv_transform_discrete_cosine_II_inverse(&signal, &intrp_ptr, &result);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_tran.h
 * @source sig\_tran.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param sig Pointer to matrix containing data to be transformed.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param result Pointer to matrix that will contain the transform.
 *               Must be pre-allocated and of the same size and type as sig.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_transform_discrete_cosine_II_inverse(const double_mat *sig, void *intrp_ptr, double_mat *result);#
 * \end{itemize}
 * @see siguniv_transform_discrete_cosine_II
 * @see siguniv_transform_discrete_cosine_II_2d
 * @see siguniv_transform_discrete_cosine_II_inverse_pow2
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_transform_discrete_cosine_II_inverse(
  const univ_mat *sig,
  void *intrp_ptr,
  univ_mat *result );


/* Function documented with universal matrix function, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_transform_discrete_cosine_II_inverse(
  const double_mat *sig,
  void *intrp_ptr,
  double_mat *result);


/** Two-dimensional discrete cosine transform of type II.
 * Computes the 2D discrete cosine transform of type II. The inverse
 * transform is computed if the parameter inverse\_flag is TRUE.
 *
 * The 2D transform is obtained by applying the 1D version of
 * the algorithm to the rows and columns independently.
 * For the direct transform, the 1D algorithm is
 * first applied to the rows of the input matrix, and then to the columns
 * of the row-transformed matrix. For the inverse transform, the 1D
 * algorithm is applied first to the columns of the input
 * matrix and then to the rows of the column-transformed matrix.
 * See \Ref{siguniv_transform_discrete_cosine_II} and \Ref{siguniv_transform_discrete_cosine_II_inverse} for more details.
 *
 * Reference: K.R.Rao and P.Yip, ``Discrete Cosine Transform,''
 *            Academic Press, 1990.
 *
 * @limits Not all data types are supported.
 *
 * @usage #err = siguniv_transform_discrete_cosine_II_2d(&signal, inverse_flag, &intrp_ptr, &result);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_tran.h
 * @source sig\_tran.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param sig Pointer to matrix containing the original data.
 * @param inverse_flag Boolean flag set to TRUE to compute inverse
 *                     transform, otherwise FALSE.
 * @param intrp_ptr Pointer to implement user interrupt.
 * @param result Pointer to matrix for result. Must be
 *               pre-allocated and of the same size and type as sig
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_transform_discrete_cosine_II_2d( const double_mat *sig, boolean inverse_flag, void *intrp_ptr, double_mat *result);#
 * \end{itemize}
 * @see siguniv_transform_discrete_cosine_II
 * @see siguniv_transform_discrete_cosine_II_inverse
 * @see siguniv_transform_discrete_cosine_II_2d_pow2
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_transform_discrete_cosine_II_2d(
  const univ_mat *sig, boolean inverse_flag, void *intrp_ptr,
  univ_mat *result);


/* Function documented with universal matrix function, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_transform_discrete_cosine_II_2d(
  const double_mat *image, boolean inverse_flag, void *intrp_ptr,
  double_mat *result );


/** Discrete Fourier transform and its inverse.
 * Calculates the DFT (and inverse DFT) of a signal.
 * If the input matrix has multiple columns, calculates the DFT on each
 * column separately. This is useful for multi-channel signals, with one
 * channel per column. The inverse transform is calculated if the boolean
 * parameter inverse\_flag is set to TRUE.
 *
 * The forward transform does not scale the data; the inverse transform scaled
 * the data by the inverse of the number of elements (or number of rows in the
 * case of a multichannel signal)
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_tran.h
 * @source sig\_tran.c
 * @library signal
 * @usage #err = siguniv_transform_discrete_fourier(&signal, inverse_flag, &intrp_ptr, &result);#
 * @limits For a single-channel, 1D signal, the data must be arranged in
 * a single-column complex matrix. For multi-channel signals, each channel
 * must be contained in a column of a matrix.
 * @return Standard mutils error/OK code.
 * @param sig Pointer to matrix containing the original data.
 * @param inverse_flag Boolean flag. Set to TRUE to calculate inverse DFT,
 *                     otherwise set to FALSE.
 * @param intrp_ptr Pointer to implement user interrupt.
 * @param result Pointer to matrix for the result. Must be
 *               pre-allocated, of the same size as sig, and of type dcomplex.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_transform_discrete_fourier(const double_mat *sig, boolean inverse_flag, void *intrp_ptr, dcomplex_mat *result);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigcpx_transform_discrete_fourier(const dcomplex_mat *sig, boolean inverse_flag, void *intrp_ptr, dcomplex_mat *result);#
 * \end{itemize}
 * @see siguniv_transform_discrete_cosine_II
 * @see siguniv_transform_discrete_cosine_II_inverse
 * @see siguniv_transform_discrete_cosine_II_2d
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_transform_discrete_fourier(
  const univ_mat *sig, boolean inverse_flag, void *intrp_ptr,
  univ_mat *result);


/* Function documented with universal matrix function, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_transform_discrete_fourier(
  const double_mat *sig, boolean inverse_flag, void *intrp_ptr,
  dcomplex_mat *result );


/* Function documented with universal matrix function, above */
MUTIL_LIBEXPORT mutil_errcode sigcpx_transform_discrete_fourier( const dcomplex_mat *sig,
  boolean inverse_flag, void *intrp_ptr, dcomplex_mat *result);


#ifdef __cplusplus
}
#endif

#endif /* #ifndef IN_SIG_TRAN_H_ */
