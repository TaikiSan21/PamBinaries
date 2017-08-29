
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/sig_conv.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_SIG_CONV_H_
#define IN_SIG_CONV_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "ut_err.h"
#include "mat_type.h"
#include "mat_assn.h"

/* This file contains function declarations for the convolution and correlation
   functions.
*/

#ifdef __cplusplus
extern "C" {
#endif


/** One- or two-dimensional signal convolution with arbitrary step and phase.
 * This function convolves a kernel with a signal,
 * and writes the result to a previously-allocated matrix.
 * Note that convolution involves flipping the kernel about its midpoint
 * for the one-dimensional case, or a 180 degree rotation of the kernel for
 * the two-dimensional case.
 * The step parameters control the size of the
 * lag in the convolution, and thus allow for convolution with
 * downsampling. For example, to convolve and downsample by 2, use
 * row\_step = 2 and col\_step = 2. For regular convolution, set both step
 * parameters to 1.
 *
 * The overlap parameters control the initial overlap of the kernel with the
 * input matrix. For a zero phase convolution, they should be set to the
 * number of rows and columns of the kernel minus their halves. For
 * example, row\_overlap = num\_rows - (num\_rows / 2), keeping in mind that
 * these are integer operations. For the traditional convolution
 * the overlap parameters should be set to 1, and the phase shift will be
 * num\_rows - (num\_rows / 2) - row\_overlap (and similarly for the columns).
 *
 * The matrix boundaries can be handled by zero padding, continuation,
 * periodic padding, or by reflection; see \Ref{_mutil_boundary_type}
 * for more information.
 *
 * For one-dimensional convolution, the input and the kernel matrices should
 * both be either row or column vectors, otherwise the function will perform
 * a two-dimensional the convolution, and extend the input according to the
 * boundary type provided.
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
 * @include sig\_conv.h
 * @source sig\_conv.c
 * @library signal
 * @usage  #err_code = siguniv_convolve(&in_sig, &kernel, row_step, col_step, row_overlap, col_overlap, MUTIL_BOUNDARY_ZERO, intrp_ptr, &out_sig);#
 * @return     Standard mutils error/OK code.
 * @param in_sig      Pointer to the input signal matrix.
 * @param kernel      Pointer to matrix representing the kernel.
 * @param row_step    Step size for downsampling the number of rows.
 * @param col_step    Step size for downsampling the number of columns.
 * @param row_overlap Number of kernel rows initially overlapping the
 *   input matrix. Must be between 1 and the number of rows of the kernel.
 * @param col_overlap Number of kernel columns initially overlapping the
 *   input matrix. Must be between 1 and the number of columns of the kernel.
 * @param boundary    Boundary condition to use.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param out_sig     Pointer to the output signal matrix, which may not be the
 *   same as the input matrix; its dimensions cannot be bigger than M + L - 1,
 *   where M is the numbers of rows (or columns) of the input, and L is the
 *   number of rows (or columns) of the kernel.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_convolve(const double_mat *in_sig, const double_mat *kernel, sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, double_mat *out_sig);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigflt_convolve(const float_mat *in_sig, const float_mat *kernel, sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, float_mat *out_sig);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigu8_convolve(const uint8_mat *in_sig, const uint8_mat *kernel, sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, uint8_mat *out_sig);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigu16_convolve(const uint16_mat *in_sig, const uint16_mat *kernel, row_step, col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, uint16_mat *out_sig);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigu32_convolve(const uint32_mat *in_sig, const uint32_mat *kernel,  sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, uint32_mat *out_sig);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigs32_convolve(const sint32_mat *in_sig, const sint32_mat *kernel,  sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, sint32_mat *out_sig);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigcpx_convolve( const dcomplex_mat *in_sig, const dcomplex_mat *kernel, sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, dcomplex_mat *out_sig );#
 * \end{itemize}
 * @see siguniv_correlate
 * @see _mutil_boundary_type
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_convolve( const univ_mat *in_sig,
    const univ_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, univ_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_convolve( const double_mat *in_sig,
    const double_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, double_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigflt_convolve( const float_mat *in_sig,
    const float_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, float_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigu8_convolve( const uint8_mat *in_sig,
    const uint8_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, uint8_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigu16_convolve( const uint16_mat *in_sig,
    const uint16_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, uint16_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigu32_convolve( const uint32_mat *in_sig,
    const uint32_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, uint32_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigs32_convolve( const sint32_mat *in_sig,
    const sint32_mat *kernel,  sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, sint32_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigcpx_convolve( const dcomplex_mat *in_sig,
    const dcomplex_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, dcomplex_mat *out_sig );

/** One- or two-dimensional signal correlation, with arbitrary step and phase.
 * This function correlates the input kernel directly
 * with the input matrix, and writes the result to a previously-allocated
 * output matrix.  For kernels that are both horizontally and vertically
 * symmetric, the convolution (\Ref{siguniv_convolve})
 * and correlation functions should give the same result.
 *
 * The step parameters control the size of the
 * lag in the correlation, and thus allow for correlation with
 * downsampling. For example, to correlate and downsample by 2, use
 * row\_step = 2 and col\_step = 2. For regular correlation, set both step
 * parameters to 1.
 *
 * The overlap parameters control the initial overlap of the kernel with the
 * input matrix. For a zero phase correlation, they should be set to the
 * number of rows and columns of the kernel minus their halves. For
 * example, row\_overlap = num\_rows - (num\_rows / 2), keeping in mind that
 * these are integer operations. For the traditional correlation,
 * the overlap parameters should be set to 1, and the phase shift will be
 * num\_rows - (num\_rows / 2) - row\_overlap (and similarly for the columns).
 *
 * The matrix boundaries can be handled by zero padding, continuation,
 * periodic padding, or by reflection; see \Ref{_mutil_boundary_type}
 * for more information.
 *
 * For one-dimensional correlation, the input and the kernel matrices should
 * both be either row or column vectors, otherwise the function will perform
 * a two-dimensional the correlation, and extend the input according to the
 * boundary type provided.
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
 * @include sig\_conv.h
 * @source sig\_conv.c
 * @library signal
 * @usage  #err_code = siguniv_correlate(&in_sig, &kernel, row_step, col_step, row_overlap, col_overlap, MUTIL_BOUNDARY_ZERO, intrp_ptr, &out_sig);#
 * @return           Standard mutils error/OK code.
 * @param in_sig      Pointer to the input signal matrix.
 * @param kernel      Pointer to matrix representing the kernel.
 * @param row_step    Step size for downsampling the number of rows.
 * @param col_step    Step size for downsampling the number of columns.
 * @param row_overlap Number of kernel rows initially overlapping the
 *   input matrix. Must be between 1 and the number of rows of the kernel.
 * @param col_overlap Number of kernel columns initially overlapping the
 *   input matrix. Must be between 1 and the number of columns of the kernel.
 * @param boundary   Boundary condition to use.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param out_sig    Pointer to the output signal matrix, which may not be the
 *   same as the input matrix; its dimensions cannot be bigger than M + L - 1,
 *   where M is the numbers of rows (or columns) of the input, and L is the
 *   number of rows (or columns) of the kernel.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_correlate(const double_mat *in_sig, const double_mat *kernel, sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, double_mat *out_sig);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigcpx_correlate(const dcomplex_mat *in_sig, const dcomplex_mat *kernel, sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, dcomplex_mat *out_sig);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode sigu8_correlate(const uint8_mat *in_sig, const uint8_mat *kernel, sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, uint8_mat *out_sig);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode sigu16_correlate(const uint16_mat *in_sig, const uint16_mat *kernel, sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, uint16_mat *out_sig);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigu32_correlate(const uint32_mat *in_sig, const uint32_mat *kernel, sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, uint32_mat *out_sig);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigs32_correlate(const sint32_mat *in_sig, const sint32_mat *kernel, sint32 row_step, sint32 col_step, sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary, void *intrp_ptr, sint32_mat *out_sig);#
 *  \end{itemize}
 * @see siguniv_convolve
 * @see _mutil_boundary_type
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_correlate( const univ_mat *in_sig,
    const univ_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, univ_mat *out_sig);


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigcpx_correlate( const dcomplex_mat *in_sig,
    const dcomplex_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, dcomplex_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_correlate( const double_mat *in_sig,
    const double_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, double_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigflt_correlate( const float_mat *in_sig,
    const float_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, float_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigu8_correlate( const uint8_mat *in_sig,
    const uint8_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap,  mutil_boundary_type boundary,
    void *intrp_ptr, uint8_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigu16_correlate( const uint16_mat *in_sig,
    const uint16_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, uint16_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigu32_correlate( const uint32_mat *in_sig,
    const uint32_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, uint32_mat *out_sig );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigs32_correlate( const sint32_mat *in_sig,
    const sint32_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, sint32_mat *out_sig );


/** Compute correlation of matrices in matrix sets.
 * This functions takes two matrix sets containing the same number of
 * matrices and computes the correlation of each matrix in one set
 * with each matrix in the other set. The result is another matrix set,
 * with the same number of matrices as the other two sets, and with each
 * matrix containing the result of each correlation operation.
 *
 * The output matrix set and its matrices must be pre-allocated by the
 * calling function, and all sets must have the same number of matrices.
 * However, it is not required that the matrix set dimensions match.
 *
 * Restrictions on the dimensions and data types of the matrices are
 * the same as for \Ref{siguniv_correlate}.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_conv.h
 * @source sig\_conv.c
 * @library signal
 * @usage #err_code = sigset_correlate(&in_set, &kernel_set, &row_step, &col_step, &row_overlap, &col_overlap, MUTIL_BOUNDARY_ZERO, intrp_ptr, &out_set);#
 * @return Standard mutils error/OK code
 * @param in_set Pointer to input matrix set.
 * @param kernel_set Pointer to matrix set containing kernels.
 * @param row_step   Pointer to matrix of sint32 containing the downsampling
 *   step along the rows for each correlation. Its number of
 *   elements must be equal to the number of matrices in the sets.
 * @param col_step   Pointer to matrix of sint32 containing the downsampling
 *   step along the columns for each correlation. Its number
 *   of elements must be equal to the number of matrices in the sets.
 * @param row_overlap Pointer to matrix of sint32 containing the row overlap
 *   for each correlation. Its number of elements must be equal to the number
 *   of matrices in the sets.
 * @param col_overlap Pointer to matrix of sint32 containing the columns
 *   overlap for each correlation. Its number of elements must be equal to
 *   the number of matrices in the sets.
 * @param boundary   Boundary condition to use.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param out_set    Pointer to matrix set containing the results of the
 *   correlations.
 * @see sigset_convolve
 * @see siguniv_correlate
 * @see _mutil_boundary_type
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode sigset_correlate( const mat_set *in_set,
    const mat_set *kernel_set, const sint32_mat *row_step,
    const sint32_mat *col_step, const sint32_mat *row_overlap,
    const sint32_mat *col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, mat_set *out_set );


/** Compute convolution of matrices in matrix sets.
 * This functions takes two matrix sets containing the same number of
 * matrices and computes the convolution of each matrix in one set
 * with each matrix in the other set. The result is another matrix set,
 * with the same number of matrices as the other two sets, and with each
 * matrix containing the result of each convolution operation.
 *
 * The output matrix set and its matrices must be pre-allocated by the
 * calling function, and all sets must have the same number of matrices.
 * However, it is not required that the matrix set dimensions match.
 *
 * Restrictions on the dimensions and data types of the matrices are
 * the same as for \Ref{siguniv_convolve}.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_conv.h
 * @source sig\_conv.c
 * @library signal
 * @usage #err_code = sigset_convolve(&in_set, &kernel_set, &row_step, &col_step, &row_overlap, &col_overlap, MUTIL_BOUNDARY_ZERO, intrp_ptr, &out_set);#
 * @return            Standard mutils error/OK code
 * @param in_set      Pointer to input matrix set.
 * @param kernel_set  Pointer to matrix set containing kernels.
 * @param row_step    Pointer to matrix of sint32 containing the downsampling
 *   step along the rows for each convolution. Its number of
 *   elements must be equal to the number of matrices in the sets.
 * @param col_step    Pointer to matrix of sint32 containing the downsampling
 *   step along the columns for each convolution. Its number
 *   of elements must be equal to the number of matrices in the sets.
 * @param row_overlap Pointer to matrix of sint32 containing the row overlap
 *   for each convolution. Its number of elements must be equal to the number
 *   of matrices in the sets.
 * @param col_overlap Pointer to matrix of sint32 containing the columns
 *   overlap for each convolution. Its number of elements must be equal to
 *   the number of matrices in the sets.
 * @param boundary    Boundary condition to use.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param out_set     Pointer to matrix set containing the results of the
 *   convolutions.
 * @see sigset_correlate
 * @see siguniv_convolve
 * @see _mutil_boundary_type
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode sigset_convolve( const mat_set *in_set,
    const mat_set *kernel_set, const sint32_mat *row_step,
    const sint32_mat *col_step, const sint32_mat *row_overlap,
    const sint32_mat *col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, mat_set *out_set );


#ifdef __cplusplus
}
#endif

#endif /* IN_SIG_CONV_H */
