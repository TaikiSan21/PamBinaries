
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/sig_win.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file */

#ifndef IN_SIG_WIN_H_
#define IN_SIG_WIN_H_

#include "ut_plat.h"
#include "sig_type.h"
#include "mat_type.h"
#include "ut_err.h"

/* This file contains declarations of functions for creating
   window tapers frequently used for windowing signals
   in signal processing.
*/

#ifdef __cplusplus
extern "C" {
#endif

/*
**************************
Universal Matrix Functions
**************************
*/


/** Create rectangular window.
 * Takes a pre-allocated one-dimensional matrix and
 * fills it with values corresponding to a rectangular window
 * (all ones).
 *
 * @usage #err = siguniv_window_rectangle(&intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_rectangle(void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_triangle
 * @see siguniv_window_raised_cosine
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_rectangle( void *intrp_ptr,
    univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_rectangle( void *intrp_ptr,
  double_mat *win );


/** Create triangular window.
 * Takes a pre-allocated one-dimensional matrix and
 * fills it with values corresponding to a triangular window.
 *
 * The triangular window is symmetric, so if the length of the
 * window is odd the mid-point actually reaches the value of 1,
 * otherwise it does not. The end points do not reach zero.
 *
 * @usage #err = siguniv_window_triangle(&intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_triangle(void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_rectangle
 * @see siguniv_window_raised_cosine
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_triangle( void *intrp_ptr,
    univ_mat *win);


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_triangle( void *intrp_ptr,
    double_mat *win);


/** Create raised cosine window.
 * Takes a pre-allocated one-dimensional matrix and
 * fills it with values corresponding to a raised cosine window.
 *
 * The raised cosine window is symmetric, with values equal
 * or close to 1 in the center, and values closer to 0 at the
 * tails. It is effectively a rectangular window with
 * sidewalls that decay smoothly instead of abruptly.
 *
 * A value of 1 for the percent parameter is equivalent to the
 * Hanning window. A value of 0 is equivalent to a rectangular window.
 *
 * @usage #err = siguniv_window_raised_cosine(&intrp_ptr, percent, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param percent Percent of the window to be raised cosine. Must be
 *                between 0 and 1 inclusive.
 * @param win Pointer to matrix to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_raised_cosine(double percent, void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_triangle
 * @see siguniv_window_hanning
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_raised_cosine( double percent,
    void *intrp_ptr,
    univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_raised_cosine( double percent,
    void *intrp_ptr,
    double_mat *win );


/** Create Hanning window.
 * Takes a pre-allocated one-dimensional matrix and
 * fills it with values corresponding to a Hanning window.
 *
 * The Hanning window is symmetric. For odd-length, the peak value of
 * window is 1. For even-length, the peak does not reach 1.
 * The end points do not reach zero.
 *
 * @usage #err = siguniv_window_hanning(&intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_hanning(void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_triangle
 * @see siguniv_window_raised_cosine
 * @see siguniv_window_hamming
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_hanning( void *intrp_ptr,
    univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_hanning( void *intrp_ptr,
    double_mat *win );


/** Create Hamming window.
 * Takes a pre-allocated one-dimensional matrix and
 * fills it with values corresponding to a Hamming window.
 *
 * The Hamming window is symmetric. For odd-length, the peak value of
 * window is 1. For even-length, the peak does not reach 1.
 * The end points do not reach zero.
 *
 * @usage #err = siguniv_window_hamming(&intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_hamming(void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_raised_cosine
 * @see siguniv_window_hanning
 * @see siguniv_window_blackman
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_hamming( void *intrp_ptr,
    univ_mat *win );


/* Function documented with universal matrix version, below */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_hamming( void *intrp_ptr,
    double_mat *win );


/** Create Blackman window.
 * Takes a pre-allocated one-dimensional matrix and
 * fills it with values corresponding to a Blackman window.
 *
 * The Blackman window is symmetric. For odd-length, the peak value of
 * window is 1. For even-length, the peak does not reach 1.
 * The end points do not reach zero.
 *
 * @usage #err = siguniv_window_blackman(&intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_blackman(void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_hanning
 * @see siguniv_window_hamming
 * @see siguniv_window_nuttall
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_blackman( void *intrp_ptr,
    univ_mat *win );


/* Function documented with universal matrix version, below */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_blackman( void *intrp_ptr,
    double_mat *win );


/** Create Nuttall window.
 * Takes a pre-allocated one-dimensional matrix and
 * fills it with values corresponding to a Nuttall window.
 *
 * The Nuttall window is symmetric. For odd-length, the peak value of
 * window is 1. For even-length, the peak does not reach 1.
 * The end points reach zero.
 *
 * @usage #err = siguniv_window_nuttall(&intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_nuttall(void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_triangle
 * @see siguniv_window_rectangle
 * @see siguniv_window_kaiser
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_nuttall( void *intrp_ptr,
    univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_nuttall( void *intrp_ptr,
    double_mat *win );


/** Create Gaussian window.
 * Takes a pre-allocated one-dimensional matrix and
 * fills it with values corresponding to a Gaussian window.
 *
 * The argument alpha is the reciprocal of the standard deviation. So, for
 * alpha = 4, the window will be a Gaussian taper between
 * two standard deviations above and below the mean.
 *
 * The Gaussian window is symmetric. For odd-length, the peak value of the
 * window is 1. For even-length, the peak does not reach 1.
 * The end points do not reach zero.
 *
 * @usage #err = siguniv_window_gaussian(alpha, &intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param alpha Reciprocal of standard deviation.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_gaussian(double alpha, void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_blackman
 * @see siguniv_window_nuttall
 * @see siguniv_window_kaiser
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_gaussian( double alpha,
    void *intrp_ptr,
    univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_gaussian( double alpha,
    void *intrp_ptr,
    double_mat *win );


/** Create Kaiser window.
 * Takes a pre-allocated one-dimensional matrix and
 * fills it with values corresponding to a Kaiser window.
 *
 * This function uses a closed form expression for
 * the Kaiser window in terms of a Bessel function.
 * \[ w(x,\theta)=\frac{J_{0}(\theta*\sqrt{(1-x^{2})})}{J_0(\theta)} \]
 *
 * The Kaiser window is symmetric. For odd-length, the peak value of
 * window is 1. For even-length, the peak does not reach 1.
 * The end points do not reach zero.
 *
 * @usage #err = siguniv_window_kaiser(res, &intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param res Resolution parameter for the taper.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix to be filled.
 *  @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_kaiser(double res, void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_triangle
 * @see siguniv_window_rectangle
 * @see siguniv_window_chebyshev
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_kaiser( double res,
    void *intrp_ptr,
    univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_kaiser(
  double beta,
  void *intrp_ptr,
  double_mat *win );


/** Create a Chebyshev Window.
 * Create a Chebyshev window.
 *
 * @usage #err = siguniv_window_chebyshev(sidelobe, &intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param sidelobe Level of the sidelove in dB.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_chebyshev(double sidelobe, void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_triangle
 * @see siguniv_window_rectangle
 * @see siguniv_window_born_jordan
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_chebyshev( double sidelobe,
    void *intrp_ptr,
    univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_chebyshev(
  double sidelobe,
  void *intrp_ptr,
  double_mat *win );

/** Create a Born-Jordan Window.
 * Create a Born-Jordan window taper. Takes a pre-allocated, one-dimensional
 * matrix and fills it with the appropriate values.
 *
 * @usage #err = siguniv_window_born_jordan(&intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_born_jordan(void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_triangle
 * @see siguniv_window_rectangle
 * @see siguniv_window_kaiser
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_born_jordan( void *intrp_ptr,
  univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_born_jordan( void *intrp_ptr,
  double_mat *win);

/** Create multiple sinusoidal tapers.
 * Takes a pre-allocated, K x N
 * matrix and fills the kth row with an N point kth ordered sinusoidal taper.
 * These tapers are typically used in a multitaper spectral density
 * function estimation scheme.
 *
 * @usage #err = siguniv_window_sinusoidal(&intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix of type MUTIL\_DOUBLE to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_sinusoidal(void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_rectangle
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_sinusoidal( void *intrp_ptr,
  univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_sinusoidal( void *intrp_ptr,
  double_mat *win);


/** Create a Parzen window.
 * Takes a pre-allocated, one-dimensional
 * matrix and fills it with the appropriate values.
 * This window is typically used in a lag window spectral density
 * function estimation scheme.
 *
 * @usage #err = siguniv_window_parzen( cutoff, &intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param cutoff An integer denoting a truncation point. The window
 *               tapers down smoothly towards the cutoff
 *               and beyond it the window values are zero.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix of type MUTIL\_DOUBLE to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_parzen( sint32 cutoff, void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_papoulis
 * @see siguniv_window_daniell
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_parzen( sint32 cutoff, void *intrp_ptr,
  univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_parzen( sint32 cutoff, void *intrp_ptr,
  double_mat *win);

/** Create a Papoulis window.
 * Takes a pre-allocated, one-dimensional
 * matrix and fills it with the appropriate values.
 * This window is typically used in a lag window spectral density
 * function estimation scheme.
 *
 * @usage #err = siguniv_window_papoulis(&intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param cutoff An integer denoting a truncation point. The window
 *               tapers down smoothly towards the cutoff
 *               and beyond it the window values are zero.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix of type MUTIL\_DOUBLE to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_papoulis(void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_parzen
 * @see siguniv_window_daniell
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_papoulis( sint32 cutoff, void *intrp_ptr,
  univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_papoulis( sint32 cutoff, void *intrp_ptr,
  double_mat *win);


/** Create a Daniell window.
 * Takes a pre-allocated, one-dimensional
 * matrix and fills it with the appropriate values.
 * This window is typically used in a lag window spectral density
 * function estimation scheme.
 *
 * @usage #err = siguniv_window_daniell(&intrp_ptr, &win);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param roughness Controls the degree of averaging that is performed
 *               on the preliminary direct spectral density function estimate.
 *               The smaller the roughness, the greater the amount of smoothing.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param win Pointer to matrix of type MUTIL\_DOUBLE to be filled.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_window_daniell(void *intrp_ptr, double_mat *win);#
 * \end{itemize}
 * @see siguniv_window_parzen
 * @see siguniv_window_papoulis
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_window_daniell( double roughness, void *intrp_ptr,
  univ_mat *win );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_window_daniell( double roughness, void *intrp_ptr,
  double_mat *win);

/** Oracle function for obtaining a particular taper/window.
 * Used to obtain any of the supported tapers.
 * For the single taper case, set either nrow or ncol to unity. For the multitaper case,
 * the number of rows denotes the number of tapers and the number of columns
 * denotes the number of samples per taper. Use the param argument to pass in
 * any additional information needed to create the taper based on the following chart:
 *
 * Taper                          param description
 * ---------------------------    -----------------------------
 * SIG\_TAPER\_RAISED\_COSINE        Fractional span of flat portion ( [0.0,1.0] ).
 * SIG\_TAPER\_GAUSSIAN             Reciprocal of standard deviation ( > 0.0 ).
 * SIG\_TAPER\_KAISER               Resolution. Shape factor used to trade sidelobe
 *                                amplitude for mainlobe width. As the resolution
 *                                increases, so does the mainlobe width ( >= 0.0 ).
 * SIG\_TAPER\_CHEBYSHEV            Sidelobe level in decibels ( > 0.0 ).
 * SIG\_TAPER\_RAISED\_COSINE    p x 100 percentage cosine taper ( [0.0,1.0] ).
 * SIG\_TAPER\_PARZEN               Cutoff index, m. Beyond m, the taper is zero ( > 1 ).
 * SIG\_TAPER\_PAPOULIS             Cutoff index, m. Beyond m, the taper is zero ( > 0 ).
 * SIG\_TAPER\_DANIELL              Roughness, m. As m increases, so does the window smoothness ( > 0.0 ).
 *
 * The param argument is ignored otherwise.
 * See the help documentation for each taper function for more information.
 *
 * @usage #err = siguniv_taper( SIG_TAPER_PAPOULIS, 1, ncol, ncol / 2, TRUE, &intrp_ptr, &result );#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_win.h
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param taper An enum of type \Ref{_sig_taper} used to denoting the taper type.
 * @param nrow The number of rows in the resulting matrix.
 * @param ncol The number of columns in the resulting matrix.
 * @param param Addition parameter needed to create the taper. See the chart
 *              above for more information.
 * @param normalize A logical flag. If TRUE, the taper is normalized to have unit energy.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @param result Pointer to universal matrix of type MUTIL\_DOUBLE containing the taper.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode sigdbl_taper( const sig_taper taper, const sint32 nrow, const sint32 ncol, const double param, const boolean normalize, void *intrp_ptr, double_mat *result );#
 * \end{itemize}
 * @see siguniv_window_rectangle
 * @see siguniv_window_triangle
 * @see siguniv_window_raised_cosine
 * @see siguniv_window_hanning
 * @see siguniv_window_hamming
 * @see siguniv_window_blackman
 * @see siguniv_window_nuttall
 * @see siguniv_window_gaussian
 * @see siguniv_window_kaiser
 * @see siguniv_window_chebyshev
 * @see siguniv_window_born_jordan
 * @see siguniv_window_sinusoidal
 * @see siguniv_window_parzen
 * @see siguniv_window_papoulis
 * @see siguniv_window_daniell
 */
MUTIL_LIBEXPORT mutil_errcode siguniv_taper(
  const sig_taper_type  taper,
  const sint32          nrow,
  const sint32          ncol,
  const double          param,
  const boolean         normalize,
  void                 *intrp_ptr,
  univ_mat             *result );

/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode sigdbl_taper(
  const sig_taper_type  taper,
  const sint32          nrow,
  const sint32          ncol,
  const double          param,
  const boolean         normalize,
  void                 *intrp_ptr,
  double_mat      *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_SIG_WIN_H_ */
