
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_type.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_WAV_TYPE_H
#define IN_WAV_TYPE_H

#include "mat_type.h"

/* This file contains data structures, typedefs and
   declarations for the wavelet data types
 */

#ifdef __cplusplus
extern "C" {
#endif

/** Structure containing wavelet coefficient energy and associated properties.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_fdp.h
 * @source wav\_fdp.h
 * @library wavelets
 * @same
 *   #typedef struct _wave_energy wave_energy;#
 *
 * @see wavuniv_fdp_estimator_instantaneous
 * @see wavuniv_fdp_estimator_block
 */
struct _wave_energy{

  /** the energy of wavelet coefficient(s)
      in each decomposition level */

  double_mat *energy;

  /** lower limit of the user specified
      range of the FD parameters (delta) */

  double delta_min;

  /** the decomposition levels */

  sint32_mat *level;

  /** number of wavelet coefficients used to
      calculate the energy in each decomposition
      level */

  sint32_mat *n_coeff;

  /** boolean specifying whether or not
      biased estimates were used in forming
      the energy vector */

  boolean biased;

  /** number of decomposition levels */

  sint32 n_level;

  /** number of samples in time series */

  sint32 n_sample;

  /** number of levels usable for
      innovation variance and reduced
      log-likelihood estimations */

  sint32 n_level_usable;

  /** boolean specifying whether or not
      a stationary model for block MLE
      of delta should be used */

  boolean stationary;

  /** The mid-octave SDF of an FD process
      for the wavelet octave associated with
      the scaling coefficients */

  double Cprime_scaling;
};

/* See above documentation for _wave_energy */
typedef struct _wave_energy wave_energy;

/** Struct for two-dimensional Gabor wavelet transform parameters.
 * This structure contains the parameters that are necessary to compute
 * the forward two-dimensional Gabor wavelet transform.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_type.h
 * @source wav\_type.h
 * @library wavelets
 * @same \begin{itemize}
 * \item #typedef struct _wav_params_gabor_2d wav_params_gabor_2d;#
 * \end{itemize}
 */
struct _wav_params_gabor_2d {
    /** standard deviation of Gaussian at decomposition level 0 along the
     * vertical (rows) direction
     */
    double sigma;

    /** aspect ratio of two-dimensional Gaussian (columns/rows) */
    double aspect;

    /** frequency in complex exponential at decomposition level 0 */
    double freq;

    /** sampling period for kernel */
    double step;

    /** number of decomposition levels */
    sint32 num_levels;

    /** number of voices per decomposition level */
    sint32 num_voices;

    /** number of rotations of the kernel per voice */
    sint32 num_angles;

    /** starting rotation angle */
    double angle_init;

    /** angle step size separating each rotation angle */
    double angle_step;

};

/* see above documentation for _wav_params_gabor_2d for description */
typedef struct _wav_params_gabor_2d wav_params_gabor_2d;


/** Enum of two-dimensional discrete wavelet transform subband types.
 * These types identify the subbands resulting from applying the
 * two-dimensional discrete wavelet transform to a matrix, for example
 * by using the function \Ref{wavuniv_transform_discrete_wavelet_2d}.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_type.h
 * @source wav\_type.h
 * @library wavelets
 * @same #typedef enum _wav_subband_dwt_2d wav_subband_dwt_2d#
 */
enum _wav_subband_dwt_2d
{

  /** Detail coefficients along the rows: highpass along the rows followed by
   * lowpass along the columns. */
  MUTIL_SUBBAND_ROWS,

  /** Detail coefficients along the columns: lowpass along the rows followed
   * by highpass along the columns. */
  MUTIL_SUBBAND_COLS,

  /** Diagonal detail coefficients: highpass along the rows followed by
   * highpass along the columns. */
  MUTIL_SUBBAND_DIAGONAL,

  /** Smooth coefficients: lowpass along the rows followed by lowpass along
   * the columns. */
  MUTIL_SUBBAND_SMOOTH
};


/* See above documentation for _wav_dwt_2d_subband_type for explanation. */
typedef enum _wav_subband_dwt_2d wav_subband_dwt_2d;

/** Enum of Daubechies filter types.
 * These types identify the class of Daubechies filters.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_type.h
 * @source wav\_type.h
 * @library wavelets
 * @same #typedef enum _wav_filter_type#
 */
enum _wav_filter_type
{
  /** Daubechies extremal (minimum) phase filters.
   * The build-up of its partial energy sequence (cumulative summation
   * over the squares of the filter's impulse response) is the most rapid
   * of any Daubechies filter with width L/2 (where L is the filter length)
   * and with the same squared gain function. */
  WAV_FILTER_EXTREMAL_PHASE,

  /** Daubechies least asymmetric (symmlet) filters.
   * Closest to linear phase for filters having the same squared gain
   * function (as do the extremal phase and best localized families). */
  WAV_FILTER_LEAST_ASYMMETRIC,

  /** Daubechies best localized filters.
   * Represents a best localized factorization of the squared gain
   * function for Daubechies least asymmetric (LA) scaling filter. The LA
   * filters are designed to minimize the maximum deviation from the
   * phase function over all frequencies. However, since the scaling
   * filter is low-pass, the phases at high frequencies are not as important
   * as low frequencies.  The best localized family is designed to
   * penalize departures from linear phase in the lower frequencies more
   * than those of higher frequencies. The result is that the best localized
   * family has better linear phase approximation than the LA family for
   * certain filter lengths (e.g. 14 tap), but is worse for others
   * (e.g. 18 and 20 tap). */
  WAV_FILTER_BEST_LOCALIZED,

  /** Daubechies Coiflet filters.
   * The extremal phase, least asymmetric, and best localized families
   * have the same squared gain function, built by specifying certain
   * ``vanishing moment'' conditions on a wavelet function that is
   * entirely determined by the associated scaling filter. Coiflets
   * are designed by specifying vanishing moment conditions on both
   * the wavelet and scaling functions. The benefit of Coiflets is
   * that they have an excellent approximate linear phase property
   * (typically better than the LA family) but at the cost of imposing
   * somewhat jagged artifacts in resulting transform coefficients.
   * Smoother sequences are best transformed using the (smoother) LA family. */
  WAV_FILTER_COIFLET,

  /** Continuous Gaussian filter (first derivative).
   * Real-valued first derivative of a Gaussian probability
   * density function scaled to meet the square integrability condition. */
  WAV_FILTER_GAUSSIAN_I,

  /** Continuous Gaussian filter (second derivative).
   * Real-valued second derivative of a Gaussian probability
   * density function scaled to meet the square integrability condition. This
   * filter is also known as the Mexican Hat or Sombrero wavelet. */
  WAV_FILTER_GAUSSIAN_II,

  /** Continuous Morlet filters.
   * Complex-valued Morlet filters. The mother wavelet
   * is a complex exponential modulated by a Gaussian window. */
  WAV_FILTER_MORLET,

  /** Piecewise continuous Haar filter.
   * Real-valued piecewise continuous Haar filter. The mother wavelet
   * is a negated step function followed directly by a positive step function. */
  WAV_FILTER_HAAR
};

/* See above documentation for _wav_filter_type for explanation. */
typedef enum _wav_filter_type wav_filter_type;

/** Enum of discrete wavelet transform types.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_type.h
 * @source wav\_type.h
 * @library wavelets
 * @same #typedef enum _wav_transform#
 */
enum _wav_transform
{
  /** Maximum overlap discrete wavelet transform. */
  WAV_TRANSFORM_MODWT,

  /** Maximum overlap discrete wavelet packet transform. */
  WAV_TRANSFORM_MODWPT,

  /** Discrete wavelet transform. */
  WAV_TRANSFORM_DWT,

  /** Discrete wavelet packet transform. */
  WAV_TRANSFORM_DWPT,

  /** Two-dimensional discrete wavelet transform. */
  WAV_TRANSFORM_DWT2D
};

/* See above documentation for _wav_filter_type for explanation. */
typedef enum _wav_transform wav_transform;

/** Enum of wavelet transform local extrema type
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_type.h
 * @source wav\_type.h
 * @library wavelets
 * @same #typedef enum _wav_transform_peak#
 */
enum _wav_transform_peak
{
  /** Local maxima and minima. */
  WAV_TRANSFORM_PEAK_EXTREMA,

  /** Local maxima. */
  WAV_TRANSFORM_PEAK_MAXIMA,

  /** Local minima. */
  WAV_TRANSFORM_PEAK_MINIMA
};

/* See above documentation for _wav_transform_peak for explanation. */
typedef enum _wav_transform_peak wav_transform_peak;

/** Struct for parameters in a hidden Markov tree (HMT) model.
 * After declaring a \Ref{_wav_hmt_model} structure variable it
 * should initialized using the \Ref{WAVHMT_ZERO_MODEL} macro.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_type.h
 * @source wav\_type.h
 * @library wavelet
 * @same
 *   #typedef struct _wav_hmt_model wav_hmt_model;#
 */
struct _wav_hmt_model {
  /** Number of signals or images used to train the model */
  sint32                num_train;

  /** Number of hidden states for each wavelet coefficient */
  sint32                num_states;

  /** Number of scale levels in the wavelet trees */
  sint32                num_scales;

  /** Boundary condition used in computing the discrete wavelet transform */
  mutil_boundary_type   boundary;

  /** Hidden state probability mass function values at each scale and
      subband */
  double_mat            state_pmfs;

  /** Hidden state means at each scale and subband */
  double_mat            means;

  /** Hidden state variances at each scale and subband */
  double_mat            vars;

  /** Transition probabilities for each child state, parent state,
      scale and subband */
  double_mat            trans_probs;

  /** Total log-likelihood of training images */
  double                final_like;

  /** Total number of EM iterations used to train the model */
  sint32                iters;

  /** Wavelet filter type used for the disrete wavelet tranform */
  wav_filter_type       filt_type;

  /** Wavelet filter length */
  sint32                filt_taps;

  /** Row dimensions of matrices in training images(s) DWT matrix set */
  sint32_mat            dwt_row_dims;

  /** Column dimensions of matrices in training images(s) DWT matrix set */
  sint32_mat            dwt_col_dims;

  /** Sample mean of original image(s) DWT scaling coefficients */
  double                scl_mean;

  /** Sample variance of original image(s) DWT scaling coefficients */
  double                scl_var;
};

/* See above documentation for _wav_hmt_model */
typedef struct _wav_hmt_model wav_hmt_model;

/** Enum of fractionally differenced (FD) process parameter estimator.
 * These types identify the method used to estimate FD model parameters.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_type.h
 * @source wav\_type.h
 * @library wavelets
 * @same #typedef enum _wav_fdp_estimator#
 */
enum _wav_fdp_estimator{

  /** Maximum likelihood estimator. */
  WAV_FDP_MAXIMUM_LIKELIHOOD,

  /** Least squares estimator. */
  WAV_FDP_LEAST_SQUARES
};

typedef enum _wav_fdp_estimator wav_fdp_estimator;

/** Enum of wavelet shrink functions.
 * These types identify the method use for waveshrinking.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_type.h
 * @source wav\_type.h
 * @library wavelets
 * @same #typedef enum _wav_shrink_function#
 */
enum _wav_shrink_function{

  /** Hard thresholding. */
  WAV_SHRINK_FUNCTION_HARD,

  /** Soft thresholding. */
  WAV_SHRINK_FUNCTION_SOFT,

  /** Mid (hard/soft) thresholding. */
  WAV_SHRINK_FUNCTION_MID
};

typedef enum _wav_shrink_function wav_shrink_function;

/** Enum of wavelet shrink threshold estimators.
 * These types identify the method use to estimate waveshrink threshold(s).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_type.h
 * @source wav\_type.h
 * @library wavelets
 * @same #typedef enum _wav_shrink_threshold#
 */
enum _wav_shrink_threshold{

  /** Universal threshold. */
  WAV_SHRINK_THRESHOLD_UNIVERSAL,

  /** Minimax thresholding. */
  WAV_SHRINK_THRESHOLD_MINIMAX,

  /** Adaptive thresholding. */
  WAV_SHRINK_THRESHOLD_ADAPTIVE
};

typedef enum _wav_shrink_threshold wav_shrink_threshold;


/** Enum of whiteness tests for discrete wavelet packet transform crystals.
 * These types identify the method use to test the whitness (in
 * a spectral sense) of a discrete wavelet packet crystal or (in general)
 * a uniformly sampled time series.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_type.h
 * @source wav\_type.h
 * @library wavelets
 * @same #typedef enum _wav_white_test#
 */
enum _wav_white_test{

  /** Portmanteau I: Box-Pierce. */
  WAV_WHITE_TEST_PORTMANTEAU_I,

  /** Portmanteau II: Lung-Box-Pierce. */
  WAV_WHITE_TEST_PORTMANTEAU_II,

  /** Portmanteau III: Lung Box-Pierce on sample
      autocorrelations for the squares of the discrete
      wavelet packet coefficients. */
  WAV_WHITE_TEST_PORTMANTEAU_III,

  /** Cumulative periodogram test. */
  WAV_WHITE_TEST_CUMULATIVE_PERIODOGRAM
};

typedef enum _wav_white_test wav_white_test;


/** Structure for storing and tracking extra DWPT coefficients.
 * An 'extra' DWPT coefficient is defined as the last atom of a
 * given DWPT crystal if and only if the length of that crystal is odd.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_type.h
 * @source wav\_type.h
 * @library wavelet
 * @same
 *   #typedef struct _wav_dwpt_extra wav_dwpt_extra;#
 */
struct _wav_dwpt_extra {
  /** A vector to store a concatenation of all extra atoms */
  double_mat atoms;

  /** Number of extra atoms that are stored */
  sint32 nelem;

  /** A vector whose jth element denotes the index of the atoms
      vector that contains the first extra coefficient for the
      jth DWPT decomposition level, i.e., the storage location in
      atoms where the last element of crystal W(j,0) is stored. */
  sint32_mat levelmap;
};

/* See above documentation for _wav_dwpt_extra */
typedef struct _wav_dwpt_extra wav_dwpt_extra;


#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_TYPE_H */










