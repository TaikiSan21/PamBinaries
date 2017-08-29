
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_filt.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_WAV_FILT_H_
#define IN_WAV_FILT_H_

#include "mat_type.h"
#include "ut_err.h"
#include "ut_plat.h"
#include "ut_type.h"
#include "wav_type.h"

/* This file contains function declarations the wavelet mutils library. */

#ifdef __cplusplus
extern "C" {
#endif

/** Daubechies wavelet and scaling filters.
 * Ingrid Daubechies, a noted pioneer in wavelet theory, has
 * established a number of wavelet filter types, each with different
 * mathematical properties.  The wavuniv\_filters\_daubechies()
 * function calculates the wavelet coefficients and scaling
 * coefficients for a given filter type.  The wavelet coefficients
 * $h_k$ for $k=0,\ldots,L-1$ where $L$ is the filter length are
 * related to the scaling coefficients through the quadrature mirror
 * filter (QMF) relation
 * \[h_k = (-1)^{k-L} g_{L-1-k} \]
 *
 * References:
 * D. B. Percival and A. T. Walden, {\it Wavelet Methods for
 * Time Series Analysis}, Cambridge University Press, 2000.
 * I. Daubechies, ``Orthonormal Bases of Compactly Supported Wavelets'',
 * {\it Communications on Pure and, Applied Mathematics}, 41, 909-96.
 *
 * @limits Only relevant for Daubechies filter types. Inconsistent ordering
 * of the coefficients in Daubechies' book was recognized and corrected by
 * Percival (above). The``correct'' order is given here.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_filt.h
 * @source wav\_filt.c
 * @library wavelets
 * @usage #err = wavuniv_filters_daubechies(filter_length, filter_type, normalize, intrp_ptr, &filters);#
 * @return Standard mutils error/OK code.
 * @param  filter_length     The length of the wavelet filter.
 * @param  filter_type       The type of Daubechies filter.
 * @param  normalize         If TRUE, the filters are normalized by dividing
 *                           each filter coefficient by the sqrt(2)
 *                           (useful for maximum overlap wavelet transforms).
 *                           If FALSE, no normalization is used.
 * @param  intrp_ptr         Pointer for implementation of interrupt checking.
 * @param  result            Pointer to a matrix set which (upon return)
 *                           will contain two universal matrices of size
 *                           [1 x filter\_length] and of type MUTIL\_DOUBLE.
 *                           The first and second vectors of the matrix set
 *                           contain the wavelet and scaling filter,
 *                           respectively. The memory for the matrix set
 *                           is automatically allocated by the function.
 *
 * @see _wav_filter_type
 * @see wavuniv_filters_daubechies_verify
 * @see wavuniv_transform_discrete_wavelet_convolution
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_filters_daubechies_gain
 * @see wavuniv_filters_zero_phase
 * @see wavuniv_transform_coefficient_boundaries
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_filters_daubechies(
  sint32           filter_length,
  wav_filter_type  filter_type,
  boolean          normalize,
  void            *intrp_ptr,
  mat_set         *result );

/** Daubechies wavelet and scaling filter verification.
 * Verifies the necessary mathematical properties for
 * Daubechies wavelet and scaling filters. Checks are made
 * for (1) the sum of the filter coefficients, (2) the
 * energy of the filter, and (3) even shift orthogonality
 * between and within the scaling and wavelet filters.
 *
 * References:
 * D. B. Percival and  A. T. Walden, {\it Wavelet Methods for
 * Time Series Analysis}, Cambridge University Press, 2002.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_filt.h
 * @source wav\_filt.c
 * @library wavelets
 * @usage #err = wavuniv_filters_daubechies_verify( &wavelet_filter, &scaling_filter, normalize );#
 * @return Standard mutils error/OK code.
 * @param  wavelet_filter Pointer to wavelet filter.
 * @param  scaling_filter Pointer to scaling filter.
 * @param  normalize      A boolean representing the normalization used
 *                        in developing the Daubechies filters.
 *                        If TRUE, the filters are assumed to be normalized
 *                        such that the sum over the scaling filter
 *                        coefficients is equal to one, and the energy of
 *                        the wavelet or scaling filter is one half.
 *                        If FALSE, the filters are assumed to be normalized
 *                        such that the sum over the scaling filter
 *                        coefficients is equal to the square root of 2,
 *                        and the energy of the wavelet or scaling filter
 *                        is equal to one. In either case the
 *                        wavelet filter is expected to have zero mean.
 *
 * @see _wav_filter_type
 * @see wavuniv_filters_daubechies
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_filters_daubechies_verify(
  const univ_mat  *wavelet_filter,
  const univ_mat  *scaling_filter,
  boolean          normalize);


/** The gain functions for Daubechies wavelet and scaling filters.
 * Given $\{g\}$ and $\{h\}$ are the impulse responses for the scaling
 * and wavelet filters, respectively, and G1(f) and H1(f) are their
 * corresponding gain functions, then the gain functions for decomposition
 * level j > 1 is calculated in a recursive algorithm based
 * on the equations
 *
 * \begin{center}
 * \TEX{
 * \begin{eqnarray}
 *  G_j(f) &=& H(2^{j-1} f)G_{j-1}(f)\nonumber\\
 *  H_j(f) &=& H(2^{j-1} f)G_{j-1}(f)\nonumber\\\nonumber
 * \end{eqnarray}
 * }
 * \end{center}
 *
 * References:
 *
 * 1. D.B. Percival and  A.T. Walden, ``Wavelet Methods for Time Series
 * Analysis'', Cambridge University Press, 2000.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_filt.h
 * @source wav\_filt.c
 * @library wavelets
 * @usage #err = wavuniv_filters_daubechies_gain( WAV_FILTER_TYPE, filter_length, num_levels, num_fft,  normalize, intrp_ptr, &gain_frequency, &gain_wavelet, &gain_scaling );#
 * @return Standard mutils error/OK code.
 * @param  filter_type       The type of Daubechies filter
 *                           (WAV\_FILTER\_TYPE enum).
 * @param  filter_length     The length of the wavelet filter.
 * @param  num_levels        The number of decomposition levels.
 * @param  num_fft           The number of Fourier coefficients to use in
 *                           approximating the gain functions.
 * @param  normalize         Boolean indicating normalization mode. If TRUE,
 *                           the filters are normalized by sqrt(2),
 *                           otherwise no normalization is used.
 * @param  intrp_ptr         Pointer for implementation of interrupt checking.
 * @param  gain_frequency    Pointer to a universal matrix of
 *                           type MUTIL\_DOUBLE to store the normalized
 *                           frequency vector corresponding to the
 *                           wavelet and scaling filter gain functions.
 *                           The memory for this matrix
 *                           is automatically allocated by the function.
 * @param  gain_wavelet      Pointer to a universal matrix of
 *                           type MUTIL\_DCOMPLEX to store the gain function
 *                           for the wavelet filter for levels
 *                           1, ..., num\_levels.
 *                           The size of this matrix is
 *                           [num\_levels x N] where
 *                           N = max(filter\_length, num\_fft).
 *                           The memory for this matrix
 *                           is automatically allocated by the function.
 * @param  gain_scaling      Pointer to a universal matrix of
 *                           type MUTIL\_DCOMPLEX to store the gain function
 *                           for the scaling filter for levels
 *                           1, ..., num\_levels.
 *                           The size of the matrix is
 *                           [num\_levels x N] where
 *                           N = max(filter\_length, num\_fft).
 *                           The memory for this matrix
 *                           is automatically allocated by the function.
 *
 * @see _wav_filter_type
 * @see wavuniv_filters_daubechies
 * @see wavuniv_filters_zero_phase
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_filters_daubechies_gain(
  wav_filter_type   filter_type,
  sint32            filter_length,
  sint32            num_levels,
  sint32            num_fft,
  boolean           normalize,
  void             *intrp_ptr,
  univ_mat         *gain_frequency,
  univ_mat         *gain_wavelet,
  univ_mat         *gain_scaling);


/** Boundary and interior wavelet coefficient identification for the
 * DWT and MODWT.
 * The boundary wavelet and scaling coefficients are those subject to
 * circular filtering operations. The wavuniv\_coefficient\_boundaries()
 * function calculates the range of indices which span the interior
 * (or nonboundary) wavelet and scaling coefficients.
 *
 * References:
 *
 * 1. D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press, 2000.
 *
 * @limits Only relevant for DWT and MODWT definitions as given in
 * the above reference.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_filt.h
 * @source wav\_filt.c
 * @library wavelets
 * @usage #err = wavuniv_transform_coefficient_boundaries(num_levels, filter_length, num_points, transform_type, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param  num_levels       The number of decomposition levels.
 * @param  filter_length    The length of the wavelet filter.
 * @param  num_points       The number of points in the time series.
 * @param  transform_type   Specifies the type of wavelet transform.
 * @param  intrp_ptr        Pointer for implementation of interrupt checking.
 * @param  result           Pointer to a universal matrix set which (upon
 *                          return) will contain five matrices, each of
 *                          type MUTIL\_SINT32 and of size [1 x num\_levels].
 *                          The contents of the matrix set are in order:
 *                          (0) indices denoting the start of the interior
 *                          wavelet coefficients (one per level).
 *                          (1) indices denoting the end of the interior
 *                          wavelet coefficients (one per level),
 *                          (2) the total number of interior
 *                          wavelet coefficients (one per level),
 *                          (3) the total number of boundary
 *                          wavelet coefficients (one per level),
 *                          (4) the total number of all
 *                          wavelet coefficients (one per level).
 *                          The memory for this matrix set
 *                          is automatically allocated by the function.
 *
 * @see wavuniv_filters_zero_phase
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_filters_daubechies
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_coefficient_boundaries(
  sint32         n_level,
  sint32         filter_length,
  sint32         n_sample,
  wav_transform  transform_type,
  void          *intrp_ptr,
  mat_set       *result );

/** Zero phase shift factors for Daubechies symmlet and Coiflet filters.
 * Daubechies Coiflet and symmlet filters are approximate linear phase
 * filters. Consequently, the wavelet and scaling coefficients of the DWT
 * (using convolution style filtering) and MODWT can be circularly
 * shifted for approximate zero phase alignment with the original time
 * series. This function calculates the
 * circular shift factors needed to bring the wavelet and scaling
 * coefficients to approximate zero phase.
 *
 * References:
 *
 * 1. D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press, 2000.
 *
 * @limits Only relevant for DWT and MODWT definitions as given in the
 *         above reference and is valid only for Daubechies symmlet and
 *         Coiflet filters.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_filt.h
 * @source wav\_filt.c
 * @library wavelets
 * @usage #err = wavuniv_filters_zero_phase(filter_type, filter_length, n_level, intrp_ptr, &result);#
 * @return Standard mutils error/OK code.
 * @param  filter_type   Specifies the type of Daubechies filter.
 *                       Allowable type are WAV\_FILTER\_COIFLET
 *                       and WAV\_FILTER\_LEAST\_ASYMMETRIC, WAV\_FILTER\_EXTREMAL\_PHASE,
 *                       and WAV\_FILTER\_HAAR.
 * @param  filter_length The length of the wavelet (or scaling) filter.
 * @param  n_level       The number of decomposition levels.
 * @param  intrp_ptr     Pointer for implementation of interrupt checking.
 * @param  result        Pointer to a universal matrix set which (upon
 *                       return) will contain two matrices, each of
 *                       type MUTIL\_SINT32 and of size
 *                       [1 x ( num\_levels * 2 )].
 *                       The contents of the matrix set are (in order)
 *                       the zero phase shifts for the DWT and the MODWT.
 *                       MUTIL\_SINT32. The columns in each matrix correspond to
 *                       [W1 | ... | WJ | V1 | ... | VJ] where W and V are
 *                       the wavelet and scaling coefficient vectors,
 *                       respectively, and J = num\_levels is the number
 *                       of decomposition levels. A negative shift factor
 *                       implies an advance (circular shift to the left)
 *                       of the wavelet transform coefficient vectors W or V.
 *                       The memory for this matrix set
 *                       is automatically allocated by the function.
 *
 * @see _wav_filter_type
 * @see wavuniv_transform_coefficient_boundaries
 * @see wavuniv_filters_daubechies
 * @see wavuniv_transform_maximum_overlap
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_filters_zero_phase(
  wav_filter_type  filter_type,
  sint32           filter_length,
  sint32           n_level,
  void            *intrp_ptr,
  mat_set         *result );


/** Creates frequency domain filters for the continuous wavelet transform.
 *
 * Returns the frequency response of a continuous
 * wavelet filter. The choices for filters are limited to
 * Haar, Gaussian, and Morlet families. The following table
 * provides a summary of supported wavelet filters.
 * \TEX{
 * \hspace{-3cm}
 *  \begin{tabular}{|l||c|c|}
 *  \hline
 *  \textbf{Family} & \textbf{\large Time Domain} & \textbf{\large
 *    Frequency Domain} \\
 *  \hline \hline
 *  \textbf{Haar} &
 *  \begin{minipage}[t]{2.5in}
 *  \vspace{0.04in}
 *  $\psi_H(x) \equiv \left\{
 *  \begin{array}{ll}
 *  -1/\sqrt{2}, & x \in (-1,0]; \\
 *  \hphantom{-}1/\sqrt{2}, & x \in (0, 1]; \\
 *  \hphantom{-}0, & \mbox{otherwise}
 *  \end{array} \right.$
 *  \vspace{0.12in}
 *  \end{minipage} &
 *  \begin{minipage}[t]{2.5in}
 *  \vspace{0.12in}
 *  $\begin{array}{lll}
 *  \Psi_H(\omega) & = & \frac{\sqrt{2}\;i}{\omega} \, \sin^2(\omega/2) \\
 *  \Psi_H(0) & = & 0
 *  \end{array}$
 *  \end{minipage} \\
 *  \hline
 *  \textbf{Gaussian} &
 *  $\psi^{(1)}_G(x) \equiv \frac{ \sqrt{2}\, x \, e^{-x^2/2\sigma^2} }{ \sigma^{3/2}\,\pi^{1/4}}$
 *  &
 *  $\Psi^{(1)}_G(\omega) = -i\,2\,\sigma^{3/2}\,\pi^{1/4}\,\omega\,e^{-\omega^2\sigma^2/2}$
 *  \\
 *  &
 *  $\psi^{(2)}_G(x) \equiv
 *  \frac{ 2 \bigl( 1 - \frac{x^2}{\sigma^2} \bigr) \,e^{-x^2/2\sigma^2
 *    } } {\pi^{1/4}\sqrt{3\sigma}}$
 *  & $\Psi^{(2)}_G(\omega) = \sqrt{8/3} \, \sigma^{5/2} \,\pi^{1/4} \,\omega^2
 *  \, e^{-\omega^2\sigma^2/2}$ \\
 *  \hline
 *  \textbf{Morlet}
 *  &
 *  $\psi_M(x) \equiv C_{\omega_0} \, e^{-i \omega_0 x} \Bigl( e^{-x^2/2} - \sqrt{2}
 *  \, e^{-\omega_0^2/4} \, e^{-x^2} \Bigr)$
 *  &
 *  $\Psi_M(\omega) = C_{\omega_0} \, \sqrt{2\pi} \, e^{-(\omega + \omega_0)^2/2} \Bigl(\; 1 -
 *   e^{(\omega^2 + 2 \omega \omega_0 ) / 4}\;  \Bigr)$ \\
 *  & & $C_{\omega_0} = \frac{\pi^{-1/4}}{\Bigl( 1 - \frac{4}{\sqrt{3}} \,
 *   e^{-\omega_0^2/4} + \sqrt{2} \, e^{-\omega_0^2/2}\Bigr)^{1/2} }$ \\
 *   \hline
 *  \end{tabular}
 * }
 *
 * References:
 *
 * 1. D.B. Percival and  A.T. Walden, ``Wavelet Methods for Time Series
 * Analysis'', Cambridge University Press, 2000.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_filt.h
 * @source wav\_filt.c
 * @library wavelets
 * @usage #err = wavuniv_filters_continuous( filter_type, filter_arg, &frequency, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param  filter_type        Filter type. Only the wav\_filter\_type
 *                            WAV\_FILTER\_HAAR, WAV\_FILTER\_GAUSSIAN\_I,
 *                            WAV\_FILTER\_GAUSSIAN\_II, and WAV\_FILTER\_MORLET
 *                            are supported.
 * @param  filter_arg         A double value representing a secondary argument to be
 *                            passed to the filter function. If the filter
 *                            is of type WAV\_FILTER\_GAUSSIAN\_I or
 *                            WAV\_FILTER\_GAUSSIAN\_II then this parameter
 *                            represents the standard deviation $\sigma$ of a Gaussian PDF.
 *                            If the filter is of type WAV\_FILTER\_MORLET, then
 *                            this value denotes the frequency shift variable w0.
 * @param  frequency          A pointer to a pre-allocated universal matrix of type
 *                            MUTIL\_DOUBLE with either a single column or row
 *                            containing the frequencies over which the frequency
 *                            response is evaluated.
 * @param  intrp_ptr          Pointer for implementation of interrupt checking.
 * @param  result             Pointer to a pre-allocated single-row or single-column
 *                            universal matrix of type MUTIL\_DCOMPLEX containing
 *                            the same number of elements as the frequency vector.
 *                            Upon return this matrix will contain the frequency
 *                            response of the specified wavelet filter.
 *
 * @see _wav_filter_type
 * @see wavuniv_transform_continuous_wavelet
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_filters_continuous(
  const wav_filter_type  filter_type,
  const double           filter_arg,
  const univ_mat        *frequency,
  void                  *intrp_ptr,
  univ_mat              *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_FILT_H_*/

/* adding this to offset a missing \right\} */
