
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_var.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_WAV_VARIANCE_H_
#define IN_WAV_VARIANCE_H_

#include "mat_set.h"
#include "mat_type.h"
#include "ut_err.h"
#include "ut_plat.h"
#include "ut_type.h"
#include "wav_type.h"

/* This file contains function declarations the wavelet mutils library.
*/

#ifdef __cplusplus
extern "C" {
#endif

/**  Discrete wavelet variance estimation.
 *
 * The discrete wavelet variance is a useful alternative to the spectral
 * density function (SDF) and is seen as an octave-band regularization of
 * the SDF. The wavelet variance decomposes the variance of certain stochastic
 * processes on a scale-by-scale basis, and thus, is very appealing
 * to the analyst studying physical phenomena which fluctuate both within
 * and across a wide range of scale.
 *
 * @algorithm
 *
 * {\bf Discrete wavelet variance estimates using the MODWT:}
 *
 * Let N be the the number of samples in a time series $\{X_t\}$,
 * L be the length of the wavelet filter,
 * $L_j \equiv (2^j - 1)(L - 1) + 1$ be the equivalent filter width
 * at level j in a MODWT, and $\tau_j \equiv 2^{j-1}$ be the scale
 * of the data at level j for j = 1,...,J. Then the unbiased wavelet
 * variance is defined as
 *
 * \[ \hat\nu_X^2(\tau_j ) \equiv {1 \over M_j} \sum_{t = L_j - 1}^{N-1}
 * \widetilde{W}_{j,t}^2 \]
 *
 * where $\widetilde{W}_{j,t}$ are the MODWT coefficients at level j and
 * time t, and $M_j \equiv N - L_j + 1$. The unbiased wavelet variance
 * estimator avoids so-called boundary coefficients which are those
 * coefficients subject to circular filter operations in a discrete
 * wavelet transform. The biased estimator is defined as
 *
 * \[ \widetilde\nu_X^2(\tau_j ) \equiv {1 \over N} \sum_{t = 0}^{N-1}
 * \widetilde{W}_{j,t}^2, \]
 *
 * and includes the effects of the boundary coefficients.
 *
 * {\bf Discrete wavelet variance estimates using the DWT:}
 *
 * The DWT can also be used to calculate wavelet variance estimates,
 * but is not preferred over the MODWT due to its poor statistical
 * properties. Let $N_j \equiv \lfloor N / 2^j \rfloor$ be the number
 * of DWT wavelet coefficients at level j, and $L'_j \equiv \lceil
 * (L-2)(1 - 2^{-j}) \rceil$ be the number of DWT boundary coefficients
 * at level j (assuming $N_j > L'j$). Then the DWT version of the
 * unbiased wavelet variance is defined as
 *
 * \[ \hat{\hat\nu}_X^2(\tau_j ) \equiv {1 \over N - 2^jL'_j}
 * \sum_{t = L'_j - 1}^{N_j-1} W_{j,t}^2 \]
 *
 * where $W_{j,t}$ are the DWT coefficients at level j and time t.
 * Similarly, the DWT version of the biased wavelet variance is defined as
 *
 * \[ \widetilde{\widetilde\nu}_X^2(\tau_j ) \equiv {1 \over N}
 * \sum_{t = 0}^{N_j-1} W_{j,t}^2. \]
 *
 * References:
 *
 * 1. D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press, 2000.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_var.h
 * @source wav\_var.c
 * @library wavelets
 * @usage #err = wavuniv_variance( &time_series, transform_type, filter_type, filter_length, n_level, &sdf, intrp_ptr, &variance_time, &variance_block, &confidence, &edof );#
 * @return Standard mutils error/OK code.
 * @param time_series     Pointer to a pre-allocated single-row or
 *                        single-column universal matrix of type
 *                        MUTIL\_DOUBLE containing the time series.
 * @param transform_type  Specifies the type of wavelet transform
 *                        (\Ref{_wav_transform}).
 * @param filter_type     Daubechies filter type \Ref{_wav_filter_type}.
 * @param filter_length   The length of the wavelet filter.
 * @param n_level         The number of decomposition levels.
 * @param sdf          	  Pointer to pre-allocated single-row or
 *                        single-column universal matrix of type
 *                     	  MUTIL\_DOUBLE containing a discretized approximation
 *                     	  of the process spectral density function. The
 *                     	  coefficients of this argument should correspond
 *                     	  exactly with the normalized Fourier frequencies
 *                     	  f = [0, 1/P , 2/P, 3/P, ..., (M-1)/P] where
 *                     	  P = 2*(M-1) and M is the number of points in the sdf
 *                     	  vector. For example, if the sdf vector contains 5
 *                     	  elements, the corresponding frequencies will be
 *                     	  f = [0, 1/8, 1/4, 3/8, 1/2]. This argument is
 *                        is used only for calculating mode 2 EDOF.
 *                        If the EDOF mode 2 estimates are not desired,
 *                        send in a NULL pointer ( ala (univ\_mat *) NULL )
 *                        for this argument and the EDOF mode 2 and
 *                        correpsonding confidence intervals will not
 *                        be calculated.
 * @param intrp_ptr       Pointer for implementation of interrupt checking.
 * @param variance_time   Pointer to a matrix set which (upon return) will
 *                        contain the time-dependent wavelet variance
 *                        estimates. The memory for the matrix set
 *                        is automatically allocated within the function
 *                        and will containing two universal matrices of
 *                        type MUTIL\_DOUBLE and size [J x N]. The two
 *                        matrices contain (in order) the biased and
 *                        unbiased time-dependent estimates.
 * @param variance_block  Pointer to a matrix set which (upon return) will
 *                        contain the block-dependent wavelet variance
 *                        estimates. The memory for the matrix set
 *                        is automatically allocated within the function
 *                        and will contain two universal matrices of
 *                        type MUTIL\_DOUBLE and size [1 x J]. The two
 *                        matrices contain (in order) the biased and
 *                        unbiased block-dependent estimates.
 * @param confidence      Pointer to matrix set which (upon return) will
 *                        contain the confidence intervals for the
 *                        block-dependent wavelet variance
 *                        estimates. The memory for the matrix set
 *                        is automatically allocated within the function
 *                        and will contain three universal matrices of
 *                        type MUTIL\_DOUBLE and size [2 x J]. The
 *                        first and second row of each matrix hold
 *                        the low and high confidence limits, respectively.
 *                        The three matrices contain (in order) the
 *                        confidence intervals for EDOF mode 1, EDOF mode 2,
 *                        and EDOF mode 3, respectively. If confidence
 *                        intervals are not desired, send in a
 *                        NULL pointer ( ala (mat\_set *) NULL ) for
 *                        this argument and they will not be calculated.
 * @param edof            Pointer to a matrix set which (upon return) will
 *                        contain the chi-squared equivalent degrees of
 *                        freedom (EDOF) used to estimate the confidence
 *                        intervals for the block-dependent wavelet variance
 *                        estimates. The memory for the matrix set
 *                        is automatically allocated within the function
 *                        and will contain three universal matrices of
 *                        type MUTIL\_DOUBLE and size [1 x J].
 *                        The three matrices contain (in order) the
 *                        EDOF mode 1, EDOF mode 2, and EDOF mode 3,
 *                        respectively. If the EDOF are not desired,
 *                        send in a NULL pointer ( ala (mat\_set *) NULL )
 *                        for this argument and neither they nor the
 *                        confidence intervals will not be calculated.
 *
 * @see wavuniv_variance_confidence
 * @see wavuniv_variance_edof
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_variance(
 const univ_mat  *time_series,
 wav_transform    transform_type,
 wav_filter_type  filter_type,
 sint32           filter_length,
 sint32           n_level,
 const univ_mat  *sdf,
 void            *intrp_ptr,
 mat_set         *variance_time,
 mat_set         *variance_block,
 mat_set         *confidence,
 mat_set         *edof );

/** Confidence intervals for the unbiased and blocked averaged discrete wavelet
 * variance estimates.
 *
 * Given $\hat\nu_X^2(\tau_j)$ are the time independent unbiased wavelet
 * variance estimates at scale $\tau_j\equiv 2^{j-1}$ where j is the
 * decomposition level, the approximate 100(1-2p)\% confidence interval
 * is given by
 *
 * \[ \biggl[ { n\hat\nu_X^2(\tau_j) \over Q_n(1-p) } ,
 * { n\hat\nu_X^2(\tau_j) \over Q_n(p) } \biggr] \]
 *
 * where $Q_n(p)$ is the p x 100\% percentage point for the $\chi_n^2$
 * distribution.
 *
 * References:
 *
 * 1. D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press, 2000.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_var.h
 * @source wav\_var.c
 * @library wavelets
 * @usage #err = wavuniv_variance_confidence( &variance, &edof, probability, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param  variance    Pointer to pre-allocated universal matrix of type
 *                     MUTIL\_DOUBLE containing the blocked unbiased
 *                     discrete wavelet variance estimates (one for each
 *                     decomposition level). This matrix must be a
 *                     row or column vector with J elements where J is
 *                     the number of decomposition levels.
 * @param  edof        Pointer to pre-allocated universal matrix of type
 *                     MUTIL\_DOUBLE containing the equivalent
 *                     degrees of freedom estimates (one for each
 *                     decomposition level). This matrix must be a
 *                     a row or column vector with J elements where J is
 *                     the number of decomposition levels.
 * @param  probability The probability desired for the confidence intervals.
 *                     Allowable values are
 *                     [0.005, .025, .05, .95, .975, .995].
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  result      Pointer to matrix set containing two
 *                     universal matrices of type MUTIL\_DOUBLE and of size
 *                     [1 x J]. The two vectors contain (in order) the.
 *                     low and high confidence interval limits for levels
 *                     1,..., J. The memory for the matrix set is
 *                     automatically allocated within the function.
 *
 * @see wavuniv_variance
 * @see wavuniv_variance_edof
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_variance_confidence(
  const univ_mat *variance,
  const univ_mat *edof,
  double          probability,
  void           *intrp_ptr,
  mat_set        *result);

/** Equivalent degrees of freedom (EDOF) estimates for a chi-squared
 * distribution assumption on the interior wavelet coefficients. This program
 * produces three estimates of the EDOF for each level of a discrete wavelet
 * transform. The three modes are described as follows for the MODWT of an
 * an input sequence $\{X_t\}_{t=0}^{N-1}$:
 *
 * EDOF $\eta_1$ (large sample approximation that requires an SDF estimation
 * via wavelet coefficients):
 *
 * \[\eta_1 = { M_j (\hat{s}_{j,0})^2 \over \hat{A}_j} \]
 *
 * where $\hat{s}_{j,\tau}$ is the autocovariance sequence defined by
 *
 * \[ \hat{s}_{j,\tau} \equiv {1 \over M_j} \sum_{t=0}^{M_j - 1}
 * \widetilde{W}_{j,t}^{(int)} \widetilde{W}_{j,t + |\tau|}^{(int)}
 * \;\; 0\le|\tau|\le M_j-1 \]
 *
 * and $\widetilde{W}_{j,t}^{(int)}$ are the $M_j$ jth level interior
 * MODWT wavelet coefficients.
 *
 * EDOF $\eta_2$ (large sample approximation where the SDF is known a priori):
 *
 * \[\eta_2 = { 2 {\biggl( \sum_{k=1}^{\lfloor (M_j - 1)/ 2 \rfloor} C_j(f_k)
 *  \biggr)}^2  \over \sum_{k=1}^{\lfloor (M_j - 1) / 2 \rfloor} C_j^2(f_k)} \]
 *
 * where $f_k \equiv k /M_j$ and
 *
 * \[ C_j \equiv \widetilde{\mathcal{H}}_j^{(D)}(f) S_X(f)\]
 *
 * is the product of Daubechies wavelet filter squared gain function and the
 * spectral density function of $X_t$.
 *
 * EDOF $\eta_3$ (large sample approximation using a band-pass approximation
 * for the SDF):
 *
 * \[\eta_3 = \max\{M_j/2^j, 1\}. \]
 *
 * References:
 *
 * 1. D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press, 2000.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_var.h
 * @source wav\_var.c
 * @library wavelets
 * @usage #err = wavuniv_variance_edof( &interior, &num_coefs, &var_block_unbiased, &level, &sdf, filter_length, filter_type, &intrp_ptr, &result)#
 * @return Standard mutils error/OK code.
 * @param interior       Pointer to pre-allocated single-row or single-column
 *                       universal matrix of type MUTIL\_DOUBLE containing a
 *                       concatenated collection of (interior) wavelet
 *                       coefficients in order of decomposition level.
 * @param num_coefs      Pointer to pre-allocated single-row or single-column
 *                       universal matrix of type MUTIL\_SINT32 containing
 *                       the number of (interior) wavelet coefficents in each level.
 * @param var_block_unbiased Pointer to pre-allocated single-row or single-column
 *                       universal matrix of type MUTIL\_DOUBLE containing the
 *                       (unbiased) block wavelet variance estimates.
 * @param level          Pointer to pre-allocated 1D universal matrix of type
 *                       MUTIL\_SINT32 containing the decomposition levels.
 *                       This argument coordinates with the interior,
 *                       num\_coefs, and var\_block\_unbiased arguments.
 *                       It therefore does not need to contain a monotonic
 *                       continuous progression of levels.
 * @param sdf          	 Pointer to pre-allocated single-row or
 *                       single-column universal matrix of type
 *                     	 MUTIL\_DOUBLE containing a discretized approximation
 *                     	 of the process spectral density function. The
 *                     	 coefficients of this argument should correspond
 *                     	 exactly with the normalized Fourier frequencies
 *                     	 f = [0, 1/P , 2/P, 3/P, ..., (M-1)/P] where
 *                     	 P = 2*(M-1) and M is the number of points in the sdf
 *                     	 vector. For example, if the sdf vector contains 5
 *                     	 elements, the corresponding frequencies will be
 *                     	 f = [0, 1/8, 1/4, 3/8, 1/2]. This argument is
 *                       is used only for calculating mode 2 EDOF.
 *                       If the EDOF mode 2 estimates are not desired,
 *                       send in a NULL pointer ( ala (univ\_mat *) NULL )
 *                       for this argument and the EDOF mode 2 and
 *                       correpsonding confidence intervals will not
 *                       be calculated.
 * @param filter_type    Daubechies filter type \Ref{_wav_filter_type}.
 * @param filter_length  The length of the wavelet filter.
 * @param intrp_ptr      Pointer for implementation of interrupt checking.
 * @param result         Pointer to a matrix set which (upon return) will
 *                       contain the chi-squared equivalent degrees of
 *                       freedom (EDOF). The memory for the matrix set
 *                       is automatically allocated within the function
 *                       and will contain three universal matrices of
 *                       type MUTIL\_DOUBLE and size [1 x J].
 *                       The three matrices contain (in order) the
 *                       EDOF mode 1, EDOF mode 2, and EDOF mode 3,
 *                       respectively. If the EDOF are not desired,
 *                       send in a NULL pointer ( ala (mat\_set *) NULL )
 *                       for this argument and they will not be calculated.
 *
 * @see wavuniv_variance
 * @see wavuniv_variance_confidence
 * @see wavuniv_filters_daubechies
 */

MUTIL_LIBEXPORT mutil_errcode wavuniv_variance_edof(
  const univ_mat  *interior,
  const univ_mat  *num_coefs,
  const univ_mat  *variance,
  const univ_mat  *level,
  const univ_mat  *sdf,
  wav_filter_type  filter_type,
  sint32           filter_length,
  void            *intrp_ptr,
  mat_set         *result);


#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_VARIANCE_H_*/
