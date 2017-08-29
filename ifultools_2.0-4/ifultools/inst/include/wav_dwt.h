
/* $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_dwt.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file */

#ifndef IN_WAV_DWT_H
#define IN_WAV_DWT_H

#include "mat_type.h"


/* This file contains function declarations for the discrete
 * wavelet transform and its inverse, for 1D and 2D.
 * The functions are defined in wav_dwt.c
 */


#ifdef __cplusplus
extern "C" {
#endif


/** Discrete wavelet transform.
 * Compute the discrete wavelet transform (DWT) with the pyramid algorithm.
 * Given the one dimensional column vectors for the input and the coefficients
 * of the FIR lowpass and highpass analysis filters, computes the DWT by
 * convolving the input sequence with the filters using the function
 * \Ref{siguniv_convolve}. Downsampling is carried out within the convolution
 * function.
 * The result is a matrix set containing num\_levels + 1 matrices.
 * Each matrix is a column vector of details coefficients at each decomposition
 * level, ordered from finest to coarsest scale. The last matrix of the
 * matrix set contains the averages coefficients
 * at the last decomposition level (coarsest scale).
 *
 * The boundary extension rules supported by this function are
 * MUTIL\_BOUNDARY\_ZERO, MUTIL\_BOUNDARY\_CONTINUE, and
 * MUTIL\_BOUNDARY\_PERIODIC. The boundary rules are applied during
 * convolution at each decomposition level.
 * The nature of the input sequence and the type of wavelet used dictate which
 * boundary condition is most suitable for a specific case. In particular, the
 * periodic extension rule will lead to a perfectly reconstructible
 * decomposition only if the sequence length is a multiple of $2^J$, where J is
 * the decomposition level. The other boundary extension rules are not
 * guaranteed to lead
 * to perfect reconstruction, independently of the signal length.
 * The function does not check whether perfect reconstruction will
 * be possible or not.
 *
 * References:
 *
 * 1. A. Bruce, H. Y. Gao, ``Applied Wavelet Analysis with S-PLUS'', Springer,
 * 1996.
 *
 * 2. {\em S+Wavelets User's Manual}, MathSoft, Inc. , 1994.
 *
 * 3. A. Bruce, H. Y. Gao, and D. Ragozin, ``S+Wavelets: Algorithms and
 * Technical Details'', Technical Report, MathSoft Inc., January 16, 1995.
 *
 * @limits Only matrices of type MUTIL\_DOUBLE are supported.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwt.h
 * @source wav\_dwt.c
 * @library wavelets
 * @usage #err = wavuniv_transform_discrete_wavelet(&input, &lowpass, &highpass, number_of_levels, MUTIL_BOUNDARY_PERIODIC, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param invec      Pointer to sequence to be transformed. Must be a column
 *    vector.
 * @param lowpass    Pointer to column vector containing the FIR lowpass filter
 *   coefficients associated with the wavelet to be used.
 * @param highpass   Pointer to column vector containing the FIR highpass
 *   filter coefficients associated with the wavelet to be used.
 * @param num_levels Number of decomposition levels to be computed. Must be
 *   at least 1 and cannot be greater than (floor(log2(sequence length)).
 * @param intrp_ptr  Pointer to implement interrupt checking.
 * @param result     Pointer to matrix set for the result. The memory for the
 *    set and the matrices in it is allocated by this function automatically.
 *    Each matrix is a column vector of details coefficients at each
 *    decomposition level, ordered from finest to coarsest scale. The last
 *    matrix of the matrix set contains the averages coefficients
 *    at the last decomposition level (coarsest scale).
 * @same \begin{itemize}
 * \item #MUTIL_LIBEXPORT mutil_errcode wavdbl_transform_discrete_wavelet( const double_mat* invec, const double_mat* lowpass, const double_mat* highpass, sint32 num_levels, mutil_boundary_type boundary, void* intrp_ptr, mat_set* result );#
 * \end{itemize}
 * @see wavuniv_transform_discrete_wavelet_inverse
 * @see wavuniv_transform_discrete_wavelet_convolution
 * @see siguniv_convolve
 * @see _mutil_boundary_type
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_discrete_wavelet(
  const univ_mat      *invec,
  const univ_mat      *lowpass,
  const univ_mat      *highpass,
  sint32               num_levels,
  mutil_boundary_type  boundary,
  void                *intrp_ptr,
  mat_set             *result );


/* Function documented with universal matrix version above */
MUTIL_LIBEXPORT mutil_errcode wavdbl_transform_discrete_wavelet(
  const double_mat    *invec,
  const double_mat    *lowpass,
  const double_mat    *highpass,
  sint32               num_levels,
  mutil_boundary_type  boundary,
  void                *intrp_ptr,
  mat_set             *result );


/** Two-dimensional discrete wavelet transform.
 * Compute the two-dimensional discrete wavelet transform (DWT) with
 * the pyramid algorithm.
 * Given a two-dimensional matrix and the lowpass and highpass filters to
 * be applied to its rows and columns, the function computes the DWT by
 * convolving each row of the input matrix with the row lowpass and highpass
 * filters, and then convolving the results with the column lowpass and
 * highpass filters. The function \Ref{siguniv_convolve} computes convolution
 * with dyadic downsampling.
 *
 * The result is a matrix set containing 3 * num\_levels + 1 matrices. The
 * matrices in the matrix set are ordered in increasing decomposition level,
 * and within each level the matrices are ordered as follows:
 * highpass along the rows followed by lowpass along the columns,
 * lowpass along the rows followed by highpass along the columns,
 * highpass along the rows followed by highpass along the columns.
 * This arrangement is repeated for each wavelet decomposition level, except
 * for the very last one, which has an additional matrix with the results of
 * lowpass along the rows followed by lowpass along the columns.
 *
 * The boundary extension rules supported by this function are
 * MUTIL\_BOUNDARY\_ZERO, MUTIL\_BOUNDARY\_CONTINUE, and
 * MUTIL\_BOUNDARY\_PERIODIC. The boundary rules are applied during
 * convolution at each decomposition level.
 * The nature of the input sequence and the type of wavelet used dictate which
 * boundary condition is most suitable for a specific case. In particular, the
 * periodic extension rule will lead to a perfectly reconstructible
 * decomposition only if the sequence length is a multiple of $2^J$, where J is
 * the decomposition level. The other boundary extension rules are not
 * guaranteed to lead
 * to perfect reconstruction, independently of the signal length.
 * The function does not check whether perfect reconstruction will
 * be possible or not.
 *
 * References:
 *
 * 1. A. Bruce, H. Y. Gao, ``Applied Wavelet Analysis with S-PLUS'', Springer,
 * 1996.
 *
 * 2. {\em S+Wavelets User's Manual}, MathSoft, Inc. , 1994.
 *
 * 3. A. Bruce, H. Y. Gao, and D. Ragozin, ``S+Wavelets: Algorithms and
 * Technical Details'', Technical Report, MathSoft Inc., January 16, 1995.
 *
 * @limits Only matrices of type MUTIL\_DOUBLE are supported.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwt.h
 * @source wav\_dwt.c
 * @library wavelets
 * @usage #err = wavuniv_transform_discrete_wavelet_2d(&input, &lowpass_rows, &lowpass_cols, &highpass_rows, &highpass_cols, number_of_levels, MUTIL_BOUNDARY_PERIODIC, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param inmat         Pointer to matrix to be transformed.
 * @param lowpass_rows  Pointer to row vector containing the FIR lowpass filter
 *   coefficients associated with the wavelet to be used.
 * @param lowpass_cols  Pointer to column vector containing the FIR lowpass
 *   filter coefficients associated with the wavelet to be used.
 * @param highpass_rows Pointer to row vector containing the FIR highpass
 *   filter coefficients associated with the wavelet to be used.
 * @param highpass_cols Pointer to column vector containing the FIR highpass
 *   filter coefficients associated with the wavelet to be used.
 * @param num_levels    Number of decomposition levels to be computed. Must be
 *   at least 1 and cannot be greater than (floor(log2(min(input dimensions))).
 * @param intrp_ptr     Pointer to implement interrupt checking.
 * @param result        Pointer to matrix set for the result. The memory for
 *   the set and the matrices in it are allocated by this function
 *   automatically. The matrix set will contain 3 * num\_levels + 1 matrices.
 * @same \begin{itemize}
 * \item #MUTIL_LIBEXPORT mutil_errcode wavdbl_transform_discrete_wavelet_2d( const double_mat* inmat, const double_mat* lowpass_rows, const double_mat* lowpass_cols, const double_mat* highpass_rows, const double_mat* highpass_cols, sint32 num_levels, mutil_boundary_type boundary, void* intrp_ptr, mat_set* result );#
 * \end{itemize}
 * @see wavuniv_transform_discrete_wavelet
 * @see siguniv_convolve
 * @see _mutil_boundary_type
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_discrete_wavelet_2d(
  const univ_mat      *inmat,
  const univ_mat      *lowpass_rows,
  const univ_mat      *lowpass_cols,
  const univ_mat      *highpass_rows,
  const univ_mat      *highpass_cols,
  sint32               num_levels,
  mutil_boundary_type  boundary,
  void                *intrp_ptr,
  mat_set             *result );


/* Function documented with universal matrix version above */
MUTIL_LIBEXPORT mutil_errcode wavdbl_transform_discrete_wavelet_2d(
  const double_mat    *inmat,
  const double_mat    *lowpass_rows,
  const double_mat    *lowpass_cols,
  const double_mat    *highpass_rows,
  const double_mat    *highpass_cols,
  sint32                num_levels,
  mutil_boundary_type   boundary,
  void                *intrp_ptr,
  mat_set             *result );


/** Inverse discrete wavelet transform.
 * Compute the inverse discrete wavelet transform with the pyramid algorithm.
 * Given a matrix set containing the wavelet coefficients of a one dimensional
 * column vector and the FIR lowpass and highpass synthesis filters,
 * reconstructs the
 * original sequence up to a desired level, specified by the function argument
 * {\em level}. The result is a one dimensional column vector containing the
 * reconstructed sequence at the user-specified decomposition level.
 * For example, to reconstruct the original sequence,
 * set level = 0. To reconstruct the sequence in the first subspace, set
 * level = 1, etc. Memory for the result
 * is allocated by this function and the upsampling operation required
 * by the pyramid algorithm is performed by \Ref{siguniv_sample}.
 *
 * The matrices in the matrix set are assumed to be ordered in the same way
 * as returned by \Ref{wavuniv_transform_discrete_wavelet}:
 * Each matrix in the matrix set contains the details coefficients at a
 * decomposition level, and the matrices are ordered from finest to coarsest
 * scale. The last matrix in the matrix set contains the averages coefficients
 * at the coarsest scale.
 *
 * The boundary extension rules supported by this function are
 * MUTIL\_BOUNDARY\_ZERO, MUTIL\_BOUNDARY\_PERIODIC,
 * MUTIL\_BOUNDARY\_CONTINUE. However, perfect reconstruction is
 * achievable only for the periodic extension and for sequences whose
 * original length is a multiple of $2^J$, where J the maximum decomposition
 * level computed by the forward transform.
 * This function does not check whether perfect reconstruction is possible or
 * not.
 *
 * References:
 *
 * 1. A. Bruce, H. Y. Gao, ``Applied Wavelet Analysis with S-PLUS'', Springer,
 * 1996.
 *
 * 2. {\em S+Wavelets User's Manual}, MathSoft, Inc. , 1994.
 *
 * 3. A. Bruce, H. Y. Gao, and D. Ragozin, ``S+Wavelets: Algorithms and
 * Technical Details'', Technical Report, MathSoft, Inc, January 16, 1995.
 *
 * @limits Only matrices of type MUTIL\_DOUBLE are supported.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwt.h
 * @source wav\_dwt.c
 * @library wavelets
 * @usage #err = wavuniv_transform_discrete_wavelet_inverse(&coeffs, &lowpass, &highpass, level, MUTIL_BOUNDARY_PERIODIC, intrp_ptr, &reconstructed );#
 * @return Standard mutils error/OK code.
 * @param coeffs    Pointer to matrix set containing the wavelet decomposition
 *    coefficients. See the documentation above for information on how the
 *    matrices in the set should be ordered.
 * @param lowpass   Pointer to column vector of FIR lowpass synthesis filter
 *    coefficients.
 * @param highpass  Pointer to column vector of FIR highpass synthesis filter
 *    coefficients.
 * @param level     The decomposition level at which reconstruction should
 *    stop. To recover a sequence at its original scale, set level = 0.
 * @param boundary  The type of boundary extension rule to use.
 * @param intrp_ptr Pointer to implement user interrupt checking.
 * @param result    The reconstructed sequence. Memory is allocated by this
 *    function.
 * @same \begin{itemize}
 * \item #MUTIL_LIBEXPORT mutil_errcode wavdbl_transform_discrete_wavelet_inverse( const mat_set* coeffs, const double_mat* lowpass, const double_mat* highpass, sint32 level, mutil_boundary_type boundary, void* intrp_ptr, double_mat* result );#
 * \end{itemize}
 * @see wavuniv_transform_discrete_wavelet
 * @see siguniv_convolve
 * @see siguniv_sample
 * @see _mutil_boundary_type
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_discrete_wavelet_inverse(
  const mat_set       *coeffs,
  const univ_mat      *lowpass,
  const univ_mat      *highpass,
  sint32               level,
  mutil_boundary_type  boundary,
  void                *intrp_ptr,
  univ_mat            *result );


/* Function documented with universal matrix version above */
MUTIL_LIBEXPORT mutil_errcode wavdbl_transform_discrete_wavelet_inverse(
  const mat_set       *coeffs,
  const double_mat    *lowpass,
  const double_mat    *highpass,
  sint32                level,
  mutil_boundary_type   boundary,
  void                *intrp_ptr,
  double_mat          *result );


/** Two-dimensional inverse discrete wavelet transform.
 * Compute the two-dimensional inverse discrete wavelet transform with
 * the pyramid algorithm.
 * Given a matrix set containing the wavelet coefficients of a two dimensional
 * matrix and the FIR lowpass and highpass synthesis filters to be applied to
 * the rows and columns,
 * reconstructs the
 * original matrix up to a desired level, specified by the function argument
 * {\em level}. The result is a matrix containing the
 * reconstructed sequence at the user-specified decomposition level.
 * For example, to reconstruct the original sequence,
 * set level = 0. To reconstruct the sequence in the first subspace, set
 * level = 1, etc. Memory for the result
 * is allocated by this function and the upsampling operation required
 * by the pyramid algorithm is performed by \Ref{siguniv_sample}.
 *
 * The matrices in the matrix set are assumed to be ordered in the same way
 * as returned by \Ref{wavuniv_transform_discrete_wavelet_2d}:
 * Each matrix in the matrix set contains the details coefficients at a
 * decomposition level, and the matrices are ordered from finest to coarsest
 * scale. The last matrix in the matrix set contains the averages coefficients
 * at the coarsest scale.
 *
 * The boundary extension rules supported by this function are
 * MUTIL\_BOUNDARY\_ZERO, MUTIL\_BOUNDARY\_PERIODIC,
 * MUTIL\_BOUNDARY\_CONTINUE. However, perfect reconstruction is
 * achievable only for the periodic extension and for sequences whose
 * original length is a multiple of $2^J$, where J is the maximum decomposition
 * level computed by the forward transform.
 * This function does not check whether perfect reconstruction is possible or
 * not.
 *
 * References:
 *
 * 1. A. Bruce, H. Y. Gao, ``Applied Wavelet Analysis with S-PLUS'', Springer,
 * 1996.
 *
 * 2. {\em S+Wavelets User's Manual}, MathSoft, Inc. , 1994.
 *
 * 3. A. Bruce, H. Y. Gao, and D. Ragozin, ``S+Wavelets: Algorithms and
 * Technical Details'', Technical Report, MathSoft, Inc, January 16, 1995.
 *
 * @limits Only matrices of type MUTIL\_DOUBLE are supported.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwt.h
 * @source wav\_dwt.c
 * @library wavelets
 * @usage #err = wavuniv_transform_discrete_wavelet_inverse_2d(&coeffs, &lowpass_rows, &lowpass_cols, &highpass_rows, &highpass_cols, level, MUTIL_BOUNDARY_PERIODIC, intrp_ptr, &reconstructed );#
 * @return Standard mutils error/OK code.
 * @param coeffs    Pointer to matrix set containing the wavelet decomposition
 *    coefficients. See the documentation above for information on how the
 *    matrices in the set should be ordered.
 * @param lowpass_rows  Pointer to row vector containing the FIR lowpass filter
 *   coefficients associated with the wavelet to be used.
 * @param lowpass_cols  Pointer to column vector containing the FIR lowpass
 *   filter coefficients associated with the wavelet to be used.
 * @param highpass_rows Pointer to row vector containing the FIR highpass
 *   filter coefficients associated with the wavelet to be used.
 * @param highpass_cols Pointer to column vector containing the FIR highpass
 *   filter coefficients associated with the wavelet to be used.
 * @param level         Decomposition level to be reconstructed. Set level = 0
 *   to reconstruct the original matrix.
 * @param intrp_ptr     Pointer to implement interrupt checking.
 * @param result        Pointer to reconstructed matrix. Memory for the result
 *   is allocated by this function.
 * @same \begin{itemize}
 * \item #MUTIL_LIBEXPORT mutil_errcode wavdbl_transform_discrete_wavelet_inverse_2d( const double_mat* inmat, const double_mat* lowpass_rows, const double_mat* lowpass_cols, const double_mat* highpass_rows, const double_mat* highpass_cols, sint32 level, mutil_boundary_type boundary, void* intrp_ptr, mat_set* result );#
 * \end{itemize}
 * @see wavuniv_transform_discrete_wavelet_2d
 * @see wavuniv_transform_discrete_wavelet
 * @see wavuniv_transform_discrete_wavelet_inverse
 * @see siguniv_convolve
 * @see _mutil_boundary_type
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_discrete_wavelet_inverse_2d(
  const mat_set*       coeffs,
  const univ_mat*      lowpass_rows,
  const univ_mat*      lowpass_cols,
  const univ_mat*      highpass_rows,
  const univ_mat*      highpass_cols,
  sint32               level,
  mutil_boundary_type  boundary,
  void*                intrp_ptr,
  univ_mat*            result );


/* Function documented with universal matrix version above */
MUTIL_LIBEXPORT mutil_errcode wavdbl_transform_discrete_wavelet_inverse_2d(
  const mat_set*      coeffs,
  const double_mat*   lowpass_rows,
  const double_mat*   lowpass_cols,
  const double_mat*   highpass_rows,
  const double_mat*   highpass_cols,
  sint32              level,
  mutil_boundary_type boundary,
  void*               intrp_ptr,
  double_mat*         result );


#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_DWT_H*/
