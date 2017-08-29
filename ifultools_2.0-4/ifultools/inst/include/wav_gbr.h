
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_gbr.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file */

#ifndef IN_WAV_GBR_H
#define IN_WAV_GBR_H

#include "mat_type.h"
#include "wav_type.h"

/* This file contains function declarations for the Gabor
   wavelet and transform. The functions are defined in wav_gbr.c
 */

#ifdef __cplusplus
extern "C" {
#endif


/*
************************************
Functions for creating Gabor kernels
************************************
*/


/** Create a two-dimensional Gabor kernel.
 * Fills a pre-allocated universal matrix with samples of a
 * two-dimensional Gabor kernel, centered about the middle of the
 * matrix.  The two-dimensional Gabor kernel consists of a
 * two-dimensional Gaussian envelope, rotated by an angle with respect
 * to the row (r) and column (c) directions, and modulated by a
 * complex exponential along the column direction:
 *
 * \[ G(r,c)=(\frac{1}{\sqrt{\pi \lambda \sigma^{2}}})\exp\left( - \frac{(c/\lambda)^{2}+r^{2}}{2\sigma^{2}}\right) \exp\left( -j \omega_{o} c \right) \]
 *
 * Lambda is the aspect ratio (columns/rows) of the elliptical
 * Gaussian envelope, sigma the standard deviation, and omega the
 * frequency, measured in radians per matrix, in the complex
 * exponential. The coordinates (r,c) are obtained by rotating the
 * ($\hat{r},\hat{c}$) axis by an angle, so that the resulting modulated
 * Gaussian is rotated as well by the same angle:
 *
 * \[r = -\hat{c} \sin(\theta) + \hat{r} \cos(\theta)\]
 * \[c =  \hat{c} \cos(\theta) + \hat{r} \sin(\theta)\]
 *
 * Note: This function shifts the values of the real part of the Gabor
 * kernel so that its sum is zero. This is required for the kernel to
 * be an admissible wavelet. The imaginary part is always
 * antisymmetric around the center of the Gaussian, so it always sums
 * to zero.
 *
 * References:
 *
 * 1. Tai Sing Lee, ``Image Representation Using 2D Gabor Wavelets,''
 * {\em IEEE Trans. Pattern Analysis and Machine Intelligence}, vol. 18,
 * no. 10, Oct. 1996.
 *
 * 2. A. C. Bovik, M. Clark, and W. Geisler. ``Multichannel Texture
 * Analysis Using Localized Spatial Filters,'' {\em IEEE Trans. Pattern
 * Analysis and Machine Intelligence}, vol. 12, no. 1, Jan. 1990.
 *
 * @limits Only matrices of type MUTIL\_DCOMPLEX are supported because the
 *    Gabor wavelet is complex.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_gbr.h
 * @source wav\_gbr.c
 * @library wavelets
 * @usage #err_code = wavuniv_kernel_gabor_2d(sigma, aspect, freq, angle, step, norm_factor, intrp_ptr, &kernel);#
 * @return            Standard mutils error/OK code.
 * @param sigma       Standard deviation in the two-dimensional Gaussian,
 *    must be greater than 0.
 * @param aspect      Aspect ratio of the elliptical two-dimensional envelope,
 *    (columns/rows), must be greater than 0
 * @param freq        Frequency in complex exponential, must be positive.
 * @param angle       Angle of rotation of the kernel, in radians, must be
 *    in (-pi/2, pi/2].
 * @param step        The step size between successive pixels in the
 *    matrix, must be greater than 0.
 * @param norm_factor Normalization factor in front of kernel. All elements of
 *    the kernel matrix G(r,c) are multiplied by this number. Useful for
 *    creating kernels for fractionally dilated (sub-octave) Gabor wavelets,
 *    which have a different scaling (see reference 1).
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param kernel      Pointer to pre-allocated universal matrix containing
 *    result. It must contain a matrix of type MUTIL\_DCOMPLEX.
 *    For suggestions on how to allocate a big enough matrix, see
 *    \Ref{wavset_kernel_gabor_2d}.
 * @same \begin{itemize}
 * \item #MUTIL_LIBEXPORT mutil_errcode wavcpx_kernel_gabor_2d(double sigma, double aspect, double freq, double angle, double step, double norm_factor, void *intrp_ptr, dcomplex_mat *kernel);#
 * \end{itemize}
 * @see wavuniv_transform_gabor_2d
 * @see wavset_kernel_gabor_2d
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_kernel_gabor_2d( double sigma,
    double aspect, double freq, double angle, double step, double norm_factor,
    void *intrp_ptr, univ_mat *kernel );


/* Function documented with universal matrix version above */
MUTIL_LIBEXPORT mutil_errcode wavcpx_kernel_gabor_2d( double sigma,
    double aspect, double freq, double angle, double step, double norm_factor,
    void *intrp_ptr, dcomplex_mat *kernel );


/** Fill a matrix set with two-dimensional Gabor kernel matrices.
 * This function allocates memory for a matrix set and its matrices
 * and fills the matrices with samples of two-dimensional Gabor
 * kernels at various wavelet decomposition levels, voices per level,
 * and rotation angles. The matrix set dimensions will be: dims[0] =
 * number of levels, dims[1]= number of voices, dims[2] = number of
 * angles. See \Ref{_mat_set} for details on matrix sets, and
 * \Ref{wavuniv_kernel_gabor_2d} for details on the Gabor wavelet.
 *
 * Convolving each matrix in the matrix set with another matrix
 * will yield the Gabor wavelet transform of that matrix, as will
 * \Ref{wavuniv_transform_gabor_2d}, provided the same parameters are
 * used. This function offers a way to inspect the kernels used in the
 * wavelet transform.
 *
 * Memory for the matrices is allocated with
 * \Ref{matset_malloc_matrices_arbitrary_size}.  Matrix sizes depend
 * on the values of the decomposition level, the voice, and the
 * rotation angle. All matrices will be of type MUTIL\_DCOMPLEX,
 * because the Gabor kernel is complex. Each kernel is computed by a
 * call to the wavcpx\_kernel\_gabor\_2d function.
 *
 * References:
 *
 * 1. Tai Sing Lee, ``Image Representation Using 2D Gabor Wavelets,''
 * {\em IEEE Trans. Pattern Analysis and Machine Intelligence}, vol. 18,
 * no. 10, Oct. 1996.
 *
 * 2. A. C. Bovik, M. Clark, and W. Geisler. ``Multichannel Texture
 * Analysis Using Localized Spatial Filters,'' {\em IEEE Trans. Pattern
 * Analysis and Machine Intelligence}, vol. 12, no. 1, Jan. 1990.
 *
 * @algorithm
 * {\bf kernel construction:} This function tries to minimize the size
 * of the kernels while guaranteeing that the kernels are big enough
 * to contain a modulated Gaussian without truncating the tails
 * significantly. First, the number of rows is assigned a value of six
 * times the standard deviation and the number of columns is assigned
 * the number of rows times the aspect ratio (columns/rows):

 * \[ \mbox{NROW} = \frac{6 \sigma}{\mbox{step}} \]
 * \[ \mbox{NCOL} = \mbox{NROW} \times \mbox{aspect} \]

 * Then the dimensions are adjusted using the cosine and sine of the angle of
 * rotation. Odd dimensions are ensured.
 *
 * Example: Consider a kernel whose initial dimensions, as determined
 * from the standard deviation, aspect ratio, and step size are 11
 * rows and 5 columns, with the axes oriented horizontally and
 * vertically.  When rotated, the dimensions will change, for example
 * to 7 columns and 7 rows.
 * \begin{verbatim}
 * -----       -------
 * | | |       |     |
 * | | |       |-   ||
 * | | |       | - | |
 * | | |       |  |  |
 * | | |       | | - |
 * | | |       ||   -|
 * |-|-|       |     |
 * | | |       -------
 * | | |
 * | | |
 * | | |
 * -----
 * \end{verbatim}
 *
 * {\bf voices:} It is possible to create Gabor wavelets that are
 * related by a non-dyadic scale by specifying a number of voices
 * greater than 1 in the {\em params} structure. The wavelets thus
 * created are also called {\em fractionally dilated wavelets},
 * because the scales are not related by a power of 2, but rather by a
 * {\em fraction} of a power of 2. For example, a Gabor kernel G(r,c)
 * at a fixed orientation angle and at scale m can be written as:
 * \[ G(r,c) = 2^{-m} G(2^{-m}r, 2^{-m}c). \]
 * For a fractionally dilated wavelet, with N voices per scale, the
 * expression for the Gabor kernel is:
 * \[ G(r,c) = 2^{-2n/N} G(2^{-n/N}r, 2^{-n/N}c) \]
 * where n = 0, 1, ... N - 1. The consequence of using fractionally
 * dilated wavelets is that the Gabor wavelet transform coefficients
 * are defined on a finer grid than the dyadic ones, and this tends to
 * make the Gabor frame tighter. Also note that the scale factor in
 * front of the wavelet is different for the fractionally dilated
 * wavelets. For more information, see the references and the
 * documentation for the functions listed at the end.
 *
 * {\bf angles of rotation:} The desired angles of rotation of the
 * Gabor wavelets are specified within the params structure by the
 * fields params.angle\_init, num\_angles, and params.angle\_step. The
 * angle\_init field specifies the angle of rotation of the first
 * wavelet in radians, num\_angles the total number of rotation
 * angles, and angle\_step the separation between two successive
 * rotation angles. These parameters allow finer or coarser coverage
 * of the entire (-pi/2, pi/2] range, as well as the ability to limit
 * the range to a sub-interval. All angles must be in the interval
 * (-pi/2, pi/2].
 *
 * Examples: To cover the entire (-pi/2, pi/2] range with 8 angles,
 * specify angle\_init = pi/2, num\_angles = 8, and angle\_step =
 * pi/8. To cover only the (-pi/4, pi/4] range with 4 angles,
 * specify angle\_init = pi/4, num\_angles = 4, and angle\_step =
 * pi/8.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_gbr.h
 * @source wav\_gbr.c
 * @library wavelets
 * @usage #err_code = wavset_kernel_gabor_2d(&params, intrp_ptr, &kernel_set);#
 * @return            Standard mutils error/OK code.
 * @param params      Pointer to structure containing parameter values needed
 *     to compute the kernel.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param kernel_set  Pointer to matrix set for the result that has not
 *     been previously allocated.
 * @see wavuniv_kernel_gabor_2d
 * @see wavuniv_transform_gabor_2d
 * @see _wav_params_gabor_2d
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode wavset_kernel_gabor_2d( const wav_params_gabor_2d
    *params, void *intrp_ptr, mat_set *kernel_set );


/*
***************************************
Functions for calculating Gabor wavelet
***************************************
*/


/** Compute the two-dimensional Gabor wavelet transform.
 * Given a two-dimensional matrix, this function calculates the Gabor
 * wavelet transform by computing the convolution of the Gabor wavelet
 * kernel with the input matrix at each dyadic level, for each voice,
 * and for each rotation angle.  This function calls
 * \Ref{siguniv_convolve} to compute the convolution of each Gabor
 * wavelet with the input matrix.  The result is returned in a
 * \Ref{_mat_set} structure. The matrices in the result matrix set
 * will contain the results of the convolution operations, and the
 * dimensions of the matrix set will be as follows: dims[0] =
 * num\_levels, dims[1] = num\_voices, dims[2] = num\_angles.  The
 * downsample parameter controls whether the convolution is computed
 * with or without downsampling.
 *
 * References:
 *
 * 1. Tai Sing Lee, ``Image Representation Using 2D Gabor Wavelets,''
 * {\em IEEE Trans. Pattern Analysis and Machine Intelligence}, vol. 18,
 * no. 10, Oct. 1996.
 *
 * 2. A. C. Bovik, M. Clark, and W. Geisler. ``Multichannel Texture
 * Analysis Using Localized Spatial Filters,'' {\em IEEE Trans. Pattern
 * Analysis and Machine Intelligence}, vol. 12, no. 1, Jan. 1990.
 *
 * 3. I. Daubechies. {\em Ten Lectures in Wavelets}, SIAM, 1992.
 *
 * @algorithm Refer to the documentation for \Ref{wavset_kernel_gabor_2d} and
 * \Ref{wavuniv_kernel_gabor_2d} before reading further.
 *
 * {\bf Relationship between standard deviation, frequency, and
 * bandwidth:} Although the standard deviation (sigma) and frequency
 * (omega) are independent of each other as defined in the
 * \Ref{_wav_params_gabor_2d} structure, they are related by:
 * \[ \sigma = \frac{K}{\omega \lambda} \]
 * where K is a constant related to the bandwidth of the Fourier
 * transform of the kernel. If the half-amplitude bandwidth is B, then
 * K is
 * \[ K = \sqrt{2\ln{2}}(\frac{2^{B} + 1}{2^{B} - 1}). \]
 * The calling function may use these
 * equations to pre-compute the value of the frequency, based on the value
 * of sigma and the bandwidth, or it may use any positive value. However,
 * in the latter case some properties of the wavelet kernels may not be
 * satisfied. For details see references 1 and 2.
 *
 * {\bf Tightness of wavelet frames:} The Gabor wavelets are not
 * orthonormal, but may form a frame if the decomposition parameters
 * are chosen appropriately. For the Gabor wavelets to form a frame,
 * the step size and the frequency must be such that step * frequency
 * <= 2pi. However, for the frame to be tight, it is usually necessary
 * to have step * frequency << 2pi. For a fixed frequency, this can be
 * achieved by lowering the step parameter. In turn, this will lead to
 * smoother kernels, at the expense of a smaller kernel size (with all
 * other parameters kept constant). Another way to make the frame
 * tighter, and thus approximate an orthogonal base, is to oversample,
 * that is, to use a larger number of rotations, voices, and levels
 * than strictly necessary.  It is possible to reconstruct the
 * original matrix from a Gabor wavelet decomposition that forms a
 * tight frame by simply adding the wavelets themselves, scaled by the
 * coefficients of the decomposition and by a factor proportional to
 * the frame bounds. For details, see references 3 and 1. Reference 1
 * reports on combinations of parameter values that lead to a (quasi)
 * tight frame. A very tight frame is formed by the Gabor wavelets generated
 * by a step < 0.8, 3 voices per level, 20 angles of rotation, for an
 * bandwidth of 1.5. An ``acceptably tight'' frame is formed even with only
 * 1 voice per level, 8 angles of rotation, and step = 0.8. The frequency and
 * sigma values were computed from these parameters using the formulas above.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_gbr.h
 * @source wav\_gbr.c
 * @library wavelets
 * @usage #err_code = wavuniv_kernel_gabor_2d(&in_mat, &params, MUTIL_BOUNDARY_ZERO, downsample, intrp_ptr, &result);#
 * @return           Standard mutils error/OK code.
 * @param in_mat     Pointer to input universal matrix to be transformed.
 * @param params     Pointer to data structure containing parameters
 *     necessary for calculating Gabor kernels and the Gabor transform.
 * @param boundary   Type of boundary condition to use in convolution.
 * @param downsample TRUE or FALSE for downsampling during convolution.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to matrix set structure for the result. Memory
 *     is allocated by this function.
 * @see wavuniv_kernel_gabor_2d
 * @see wavset_kernel_gabor_2d
 * @see _wav_params_gabor_2d
 * @see _mutil_boundary_type
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_gabor_2d( const univ_mat
    *in_mat, const wav_params_gabor_2d *params, mutil_boundary_type boundary,
    boolean downsample, void *intrp_ptr, mat_set *result);

#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_GBR_H */
