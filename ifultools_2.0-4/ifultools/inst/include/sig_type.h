
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/sig_type.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_SIG_TYPE_H
#define IN_SIG_TYPE_H

#include "ut_type.h"
#include "ut_mem.h"

/* This file contains data structures, typedefs and
   declarations for signal data types
 */

#ifdef __cplusplus
extern "C" {
#endif

/** Enum for types for various tapers/windows. Some of them
 * are mainly used in nonparametric spectral density function
 * estimation (SDF) schemes.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include sig\_type.h
 * @source sig\_type.h
 * @library signal
 * @same #typedef enum _sig_taper_type#
 */
enum _sig_taper_type
{
  /** Rectangular. In combination with direct SDF estimator, this yields the periodogram */
  SIG_TAPER_RECTANGULAR,

  /** Triangular window. */
  SIG_TAPER_TRIANGLE,

  /** Raised cosine window. */
  SIG_TAPER_RAISED_COSINE,

  /** Hanning window. */
  SIG_TAPER_HANNING,

  /** Hamming window. */
  SIG_TAPER_HAMMING,

  /** Blackman window. */
  SIG_TAPER_BLACKMAN,

  /** Nuttall window. */
  SIG_TAPER_NUTTALL,

  /** Gaussian window. */
  SIG_TAPER_GAUSSIAN,

  /** Kaiser window. */
  SIG_TAPER_KAISER,

  /** Chebyshev window. */
  SIG_TAPER_CHEBYSHEV,

  /** Born Jordan window. */
  SIG_TAPER_BORN_JORDAN,

  /** SINUSOIDAL. Typically used in multitaper SDF estimation schemes. */
  SIG_TAPER_SINUSOIDAL,

  /** Parzen window. Typically used in lag window SDF estimation schemes. */
  SIG_TAPER_PARZEN,

  /** Papoulis window. Typically used in lag window SDF estimation schemes. */
  SIG_TAPER_PAPOULIS,

  /** Daniell window. Typically used in lag window SDF estimation schemes. */
  SIG_TAPER_DANIELL
};

/* See above documentation for _sig_taper_type for explanation. */
typedef enum _sig_taper_type sig_taper_type;

  /** Split cosine bell. */
/*   SIG_TAPER_SPLIT_COSINE_BELL, */



#ifdef __cplusplus
}
#endif

#endif /* IN_SIG_TYPE_H */
