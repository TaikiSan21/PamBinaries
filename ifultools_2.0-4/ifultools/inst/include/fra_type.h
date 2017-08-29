
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_type.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_TYPE_H
#define IN_FRA_TYPE_H

#include "ut_type.h"
#include "ut_mem.h"

/* This file contains data structures, typedefs and
   declarations for fractal data types
 */

#ifdef __cplusplus
extern "C" {
#endif

/** Structure containing machine dependent arithmetic constants.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_type.h
 * @source fra\_type.h
 * @library fractal
 * @same
 *   #typedef struct _machine machine;#
 *
 * @see fraflt_machine_constants
 */

struct _machine{

  /** radix for the floating-point representation */

  sint32 radix;

  /** number of base radix digits in the floating-point significand */

  sint32 ndigits;

  /** The rounding action:
      (0) if floating-point addition chops,
      (1) if floating-point addition rounds, but not in the IEEE style,
      (2) if floating-point addition rounds in the IEEE style,
      (3) if floating-point addition chops, and there is partial underflow,
      (4) if floating-point addition rounds, but not in the
      IEEE style, and there is partial underflow,
      (5) if floating-point addition rounds in the IEEE style,
      and there is partial underflow  */

  sint32 round;

  /** number of guard digits for multiplication with truncating arithmetic.
      (0) if floating-point arithmetic rounds, or if it
      truncates and only digit base radix digits
      participate in the post-normalization shift of the
      floating-point significand in multiplication,
      (1) if floating-point arithmetic truncates and more
      than digits base  radix  digits participate in the
      post-normalization shift of the floating-point
      significand in multiplication. */

  sint32 nguard;

  /** the largest negative integer such that
      1.0 + FLOAT(radix)**ulp_digits .NE. 1.0, except that
      ulp_digits is bounded below by  -(ndigits+3) */

  sint32 ulp_digits;

  /** the largest negative integer such that
      1.0-FLOAT(radix)**neg_ulp_digitss .NE. 1.0, except that
      neg_ulp_digits is bounded below by  -(ndigits+3) */

  sint32 neg_ulp_digits;

  /** the number of bits (decimal places if radix = 10)
      reserved for the representation of the exponent
      (including the bias or sign) of a floating-point
      number */

  sint32 nexponent;

  /** the largest in magnitude negative integer such that
      FLOAT(radix)**minexp is positive and normalized */

  sint32 minexp;

  /** the smallest positive power of BETA that overflows */

  sint32 maxexp;

  /** the smallest positive floating-point number such
      that  1.0+eps .NE. 1.0. In particular, if either
      radix = 2  or  ROUND = 0, eps = FLOAT(radix)**ulp_digits.
      Otherwise,  eps = (FLOAT(radix)**ulp_digits)/2 */

  float eps;

  /** A small positive floating-point number such that
      1.0-epsneg .NE. 1.0. In particular, if radix = 2
      or  ROUND = 0, epsneg = FLOAT(radix)**neg_ulp_digitss.
      Otherwise,  epsneg = (radix**neg_ulp_digitss)/2.  Because
      neg_ulp_digitss is bounded below by -(ndigits+3), epsneg may not
      be the smallest number that can alter 1.0 by
      subtraction. */

  float epsneg;

  /** the smallest non-vanishing normalized floating-point
      power of the radix, i.e.,  xmin = FLOAT(radix)**minexp */

  float xmin;

  /** the largest finite floating-point number.  In
      particular  xmax = (1.0-epsneg)*FLOAT(radix)**maxexp
      Note - on some machines  xmax  will be only the
      second, or perhaps third, largest number, being
      too small by 1 or 2 units in the last digit of
      the significand. */

  float xmax;
};

/* See above documentation for _machine */

typedef struct _machine machine;

/** Enum of nearest neighbor search metric.
 * These types specifiy the distance metric to use
 * in finding nearest neighbors for multidimensional
 * phase space embeddings.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_type.h
 * @source fra\_type.h
 * @library fractal
 * @same #typedef enum _fra_distance_metric#
 */
enum _fra_distance_metric
{
  /** The L1 norm.
   * This norm is defined by: sum | x[i] - y[i] |, where x and y are vectors
   * and i = 1,...,E (where E is the embedding dimension) are their components. */
  FRA_DISTANCE_L1,

  /** The L2 norm.
   * This norm is defined by: sqrt(sum ( x[i] - y[i] )^2), where x and y are
   * vectors and i = 1,...,E (where E is the embedding dimension) are their
   * components. */
  FRA_DISTANCE_L2,

  /** The L-infinity norm.
   * This norm is defined by: max | x[i] - y[i] |, where x and y are
   * vectors and i = 1,...,E (where E is the embedding dimension) are their
   * components. */
  FRA_DISTANCE_LINFINITY,

  /** The Mahalanobis norm.
   * This norm is defined by: sqrt((x-y)' * M * (x-y)), where x and y are
   * embedding vectors, ' denotes transpose, and M is a postivie definite
   * matrix. */
  FRA_DISTANCE_MAHALANOBIS
};

/* See above documentation for _fra_distance_metric for explanation. */
typedef enum _fra_distance_metric fra_distance_metric;


/** Enum of surrogate data method.
 * These types specifiy the method to use
 * in creating surrogate data for uniformly sampled
 * time series.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_type.h
 * @source fra\_type.h
 * @library fractal
 * @same #typedef enum _fra_surrogate#
 */
enum _fra_surrogate
{
  /** The phase randomization method.
   * A conjugate-symmetric complex-valued sequence is formed based on
   * the product of weighted spectral density function coefficients
   * followed by a forward DFT, resulting in a purely real sequence.
   * The weights are created via exp(i * phi), where phi are random
   * variables uniformly distributed on [0,2*PI]. */
  FRA_SURROGATE_RANDOM_PHASE,

  /** The Amplitude Adjusted Fourier Transform method.
   * This method works by (1) creating a zero mean unit variance
   * Gaussian white noise sequence and sorting it based upon the rank of the
   * indices of the original time series, (2) phase-randomizing the rank sorted
   * Gaussian sequence resulting in a sequency Y, and (3) creating a surrogate
   * series by sorting the original time series based on the rank of Y. */
  FRA_SURROGATE_AAFT,

  /** Davison-Hinkley phase and amplitude randomization method.
   * The DFT of the original time series is calcualted and modulated
   * with uniformaly distributed random phase son [0,2*PI]. Then
   * an averaging operation is used to also modify the amplitude
   * information. The result is a conjugate-symmetric complex-valued
   * sequence which is iverted via an inverse DFT. */
  FRA_SURROGATE_DAVISON_HINKLEY,

  /** The circulant embedding.
   * Davies and Harte's circulant embedding technique.
   * A conjugate-symmetric complex-valued sequence is formed based on
   * the product of weighted spectral density function coefficients
   * followed by a forward DFT, resulting in a purely real sequence.
   * The weights are created in the form X = i*Y, where X and Y are
   * zero mean, unit variance Gaussian random variables. */
  FRA_SURROGATE_CIRCULANT_EMBEDDING,

  /** Adaptive wavelet-based bootstrapping.
   * Coefficients within a subset transform of a DWPT are
   * randomly shuffled with replacement followed by an inversion
   * process to create surrogate data. The subset transform
   * should be optimized for whiteness, i.e., flatness of the SDF within
   * each subband. */
  FRA_SURROGATE_WAVELET_PACKET_BOOTSTRAP,

  /** The wavelet transform shuffle method.
   * Discrete wavelet transform coefficients are randomly shuffled
   * within each scale followed by an inverse DWT operation to
   * produce the surrogate. */
  FRA_SURROGATE_WAVELET_SHUFFLE
};

/* See above documentation for _fra_surrogate for explanation. */
typedef enum _fra_surrogate fra_surrogate;


/** Enum for types of extrema.
 * These specifiy a type of extrema to find in a given time series.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_type.h
 * @source fra\_type.h
 * @library fractal
 * @same #typedef enum _fra_extrema_type#
 */
enum _fra_extrema_type
{
  /** Time series minima. */
  FRA_EXTREMA_MINIMA = 0,

  /** Time series maxima. */
  FRA_EXTREMA_MAXIMA,

  /** All time series extrema, i.e., minima and maxima */
  FRA_EXTREMA_ALL
};

/* See above documentation for _fra_extrema_type for explanation. */
typedef enum _fra_extrema_type fra_extrema_type;



#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_TYPE_H */
