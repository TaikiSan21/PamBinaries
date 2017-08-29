
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_fdp.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_WAV_FDP_H_
#define IN_WAV_FDP_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "wav_type.h"

/* This file contains function declarations the wavelet mutils library. */

#ifdef __cplusplus
extern "C" {
#endif

/** Block-independent (instantaneous) estimation of fractionally
 * differenced (FD) model parameters.
 * The MODWT is used to calculate estimates of the FD parameter delta,
 * the variance of the FD parameter estimates and the innovation variance.
 * The user can select between maximum likelihood and least squares
 * estimators. Localized estimates may also be formed by using multiple
 * chi-squared degrees of freedom in estimating the FD model parameters.
 *
 * References:
 * D. B. Percival and  A. T. Walden, {\it Wavelet Methods for
 * Time Series Analysis}, Cambridge University Press, 2000, 340-392.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_fdp.h
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = wavuniv_fdp_estimator_instantaneous(&time_series, &levels, filter_type, filter_length, estimator, biased, dof_order, &delta_range, intrp_ptr, &delta, &variance_delta, &innovation_variance);#
 * @return Standard mutils error/OK code.
 * @param  time_series    Pointer to a pre-allocated single-row or
 *                        single-column universal matrix of type
 *                        MUTIL\_DOUBLE containing the time series.
 * @param  levels         Pointer to a pre-allocated single-row or
 *                        single-column universal matrix of type
 *                        MUTIL\_SINT32 containing the decomposition levels.
 *                        The levels can be given in any order,
 *                        but must be positive.
 * @param  filter_type    Daubechies filter type \Ref{_wav_filter_type}.
 * @param  filter_length  The length of the wavelet filter.
 * @param  estimator      The method to estimate the FD model parameters.
 *                        This argument is an enumerated type
 *                        \Ref{_wav_fdp_estimator}.
 * @param  biased         Boolean flag denoting biased or unbiased
 *                        estimates. Biased estimates are those which
 *                        use all available levels in calculating
 *                        the FD model parameters. Unbiased estimates
 *                        are calculated with only those wavelet
 *                        coefficients not subject to circular filter
 *                        operations, i.e. only the interior wavelet
 *                        coefficients are used in calculating unbiased
 *                        estimates.
 * @param  dof_order      The degree of freedom (DOF) order. The number of
 *                        chi-squared DOF used in estimating the FD
 *                        parameters is equal to 2K + 1, where K is
 *                        the dof\_order such that (K > 0). As the order
 *                        increases, the estimates will become smoother
 *                        but less localized in time.
 * @param  delta_range    Pointer to a pre-allocated two element universal
 *                        matrix of type MUTIL\_DOUBLE containing
 *                        the minimum and maximum search range
 *                        to use in estimating the FD parameter.
 * @param  intrp_ptr      Pointer for implementation of interrupt checking.
 * @param  delta          Pointer to a universal matrix to contain the
 *                        estimates of the FD parameter delta. The memory
 *                        for this output is automatically generated within
 *                        the function and (upon return) will contain a
 *                        single-column universal matrix of type MUTIL\_DOUBLE
 *                        with the same length as the original time series.
 * @param  variance_delta Pointer to a universal matrix to contain the variance
 *                        of the FD parameter estimates. The memory
 *                        for this output is automatically generated within
 *                        the function and (upon return) will contain a
 *                        single-column universal matrix of type MUTIL\_DOUBLE
 *                        with the same length as the original time series.
 * @param  innovation_variance Pointer a to universal matrix to contain the
 *                        the FD innovation variances. The memory
 *                        for this output is automatically generated within
 *                        the function and (upon return) will contain a
 *                        single-column universal matrix of type MUTIL\_DOUBLE
 *                        with the same length as the original time series.
 *
 * @see _wav_fdp_estimator
 * @see _wav_filter_type
 * @see wavuniv_fdp_estimator_block
 * @see wavuniv_fdp_bandpass_variance
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_filters_daubechies
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_fdp_estimator_instantaneous(
  const univ_mat    *time_series,
  const univ_mat    *levels,
  wav_filter_type    filter_type,
  sint32             filter_length,
  wav_fdp_estimator  estimator,
  boolean            biased,
  sint32             dof_order,
  const univ_mat    *delta_range,
  void              *intrp_ptr,
  univ_mat          *delta,
  univ_mat          *variance_delta,
  univ_mat          *innovation_variance);

/** Block-dependent estimation of fractionally differenced (FD) model
 * parameters.
 * A discrete wavelet transform is used to estimate the
 * FD parameter, the variance of the FD parameter and the
 * innovation variance for a given time series. Both a
 * maximum likelihood estimation (MLE) and weighted least
 * squares estimation (WLSE) scheme are available.
 * If an MLE scheme is chosen, then the DWT
 * \Ref{wavuniv_transform_discrete_wavelet_convolution}
 * is used for its ability to decorrelate long-memory
 * processes. If a WLSE scheme is chosen, then the MODWT
 * \Ref{wavuniv_transform_maximum_overlap}
 * is used for its known statistical wavelet variance properties.
 *
 * References:
 * D. B. Percival and  A. T. Walden, {\it Wavelet Methods for
 * Time Series Analysis}, Cambridge University Press, 2000, 340-392.
 *
 * @limits For the MLE scheme, the scaling coefficients may be
 *         used (in addition to the wavelet coefficients) by
 *         setting the boundary\_mode = TRUE (no detrending)
 *         and by setting the levels vector = {1,2,...,J} where J
 *         is the maximum number of levels in a FULL DWT. Use of
 *         the scaling coefficients in FD parameter estimation
 *         is done under an implicit assumption that the FD
 *         process is stationary ( delta < 1/2 ) and should
 *         not be used otherwise. The WLSE scheme should be
 *         limited to unbiased estimates (boundary\_mode =
 *         FALSE) since the confidence intervals for the
 *         biased estimator have not been sufficiently studied.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_fdp.h
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = wavuniv_fdp_estimator_block(&time_series, &levels, filter_type, filter_length, estimator, boundary_mode, edof_mode, &sdf, intrp_ptr, &delta, &variance_delta, &innovation_variance );#
 * @return Standard mutils error/OK code.
 * @param  time_series    Pointer to a pre-allocated single-row or
 *                        single-column universal matrix of type
 *                        MUTIL\_DOUBLE containing the time series.
 * @param  levels         Pointer to a pre-allocated single-row or
 *                        single-column universal matrix of type
 *                        MUTIL\_SINT32 containing the decomposition levels.
 *                        The levels can be given in any order,
 *                        but must be positive.
 * @param  filter_type    Daubechies filter type \Ref{_wav_filter_type}.
 * @param  filter_length  The length of the wavelet filter.
 * @param  estimator      The method to estimate the FD model parameters.
 *                        This argument is an enumerated type
 *                        \Ref{_wav_fdp_estimator}.
 * @param  boundary_mode  Boolean flag representing the different methods by
 *                        which boundary wavelet (and scaling) coefficients
 *                        are handled in calculating the FD model parameters.
 *                        If the estimator is WAV\_FDP\_MAXIMUM\_LIKELIHOOD,
 *                        this argument is used to specify a detrending
 *                        scheme (TRUE = no detrending, FALSE = detrending).
 *                        Detrending a time series means to avoid the scaling
 *                        coefficients and boundary wavelet coefficients
 *                        since they are subject to (polynomial) trend
 *                        contamination. If the estimator is
 *                        WAV\_FDP\_LEAST\_SQUARES, then this argument
 *                        is used to specify whether biased or unbiased
 *                        estimates are to be used (TRUE = biased, FALSE =
 *                        unbiased). Unbiased estimates are those which
 *                        exclude the wavelet boundary coefficients (the
 *                        scaling coefficients are always rejected in
 *                        the WLSE of FD parameters).
 * @param edof_mode       The mode by which the equivalent degrees of
 *                        freedom are calculated. This argument is
 *                        limited to 1,2, or 3 and is used only for the
 *                        WLSE scheme. See \Ref{wavuniv_variance_edof}
 *                        for details.
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
 *                        is used only for the WLSE scheme.
 * @param  delta_range    Pointer to a pre-allocated two element universal
 *                        matrix of type MUTIL\_DOUBLE containing
 *                        the minimum and maximum search range
 *                        to use in estimating the FD parameter.
 * @param  intrp_ptr      Pointer for implementation of interrupt checking.
 * @param  delta          Pointer to a double value which (upon return)
 *                        will contain a block-dependent estimate of the
 *                        FD parameter.
 * @param  variance_delta Pointer to a double value which (upon return)
 *                        will contain the variance of a block-dependent
 *                        estimate of the FD parameter.
 * @param  innovation_variance Pointer to a double value which (upon return)
 *                        will contain a block-dependent estimate of the
 *                        FD innovation variance.
 *
 * @see _wav_fdp_estimator
 * @see _wav_filter_type
 * @see wavuniv_fdp_estimator_instantaneous
 * @see wavuniv_variance_edof
 * @see wavuniv_fdp_bandpass_variance
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_filters_daubechies
 * @see wavuniv_variance
 */

MUTIL_LIBEXPORT mutil_errcode wavuniv_fdp_estimator_block(
  const univ_mat    *time_series,
  const univ_mat    *levels,
  wav_filter_type    filter_type,
  sint32             filter_length,
  wav_fdp_estimator  estimator,
  boolean            boundary_mode,
  sint32             edof_mode,
  const univ_mat    *sdf,
  const univ_mat    *delta_range,
  void              *intrp_ptr,
  double            *delta,
  double            *variance_delta,
  double            *innovation_variance );

/** Mid-octave spectral density function (SDF) estimation.
 * The wavelet and scaling filters used for wavelet decompositions
 * are nominally associated with approximate bandpass filters.
 * Specifically, at decomposition level j, the wavelet transform
 * coefficients correspond approximately to the normalized
 * frequency range of $[ 1/2^(j+1), 1 /2^j ]$. The square of the
 * wavelet coefficients are used to form the so-called wavelet
 * variance (or wavelet spectrum) which is seen as a regularization
 * of the SDF. Under an assumed fractionally differenced process (FDP)
 * model (used to model stochastic fractal processes which
 * have a power law SDF), the wavuniv\_fdp\_bandpass\_variance
 * function estimates the mid-octave SDF values. The estimates are
 * calculated assuming that the wavelet transform filters form
 * perfect (rectangular) passbands. Decomposition levels 1 and 2
 * are calculated using a second order Taylor series expansion
 * about the mid-octave frequencies while, for levels greater
 * than 2, a small angle approximation ($sin(\pi f) \approx \pi f$ )
 * is used to develop a closed form solution which is a function of FDP
 * model parameters as well as the mid-octave frequencies.
 *
 * For more information about the algorithm, see
 * D. B. Percival and  A. T. Walden, {\it Wavelet Methods for
 * Time Series Analysis}, Cambridge University Press, 2000, 343-54.
 *
 * @limits Estimates are made for the scaling filter band based upon
 *         an implicit assumption that the FD process is stationary,
 *         i.e. that delta < 0.5.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_fdp.h
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = wavuniv_fdp_bandpass_variance( &levels, delta, n_sample, intrp_ptr, &result);#
 * @return Standard mutils error/OK code.
 * @param  levels    Pointer to a pre-allocated single-row or
 *                   single-column universal matrix of type
 *                   MUTIL\_SINT32 containing the decomposition levels.
 *                   The levels can be given in any order,
 *                   but must be positive. However, for n\_sample > 0,
 *                   the levels vector must contain the values
 *                   1, 2, 3, ..., J where J is the maximum
 *                   wavelet transform decomposition level.
 * @param  delta     The fractional difference parameter. If the
 *                   scaling band estimates are desired (prompted by
 *                   setting n\_sample > 0), then delta must be less
 *                   than 0.5 since the formulae for calculating the
 *                   scaling band estimates implicitly assume
 *                   stationarity (i.e., a stationary FD process
 *                   has an FD parameter delta < 0.5).
 * @param  n_sample  The number of samples in the time series.
 *                   Although no time series is actually passed to
 *                   the wavuniv\_fdp\_bandpass\_variance function,
 *                   the n\_sample argument is used in estimating the
 *                   mid-octave SDF value over the band of frequencies
 *                   which are nominally associated with the scaling
 *                   filter in a wavelet transform.  If n\_sample > 0,
 *                   this function will append the estimate of the
 *                   average SDF value over the scaling band to the
 *                   wavelet octave estimates. If n\_sample <= 0,
 *                   only the wavelet octave estimates are returned.
 * @param  intrp_ptr Pointer for implementation of interrupt checking.
 * @param  result    Pointer to a universal matrix which (upon return)
 *                   will contain the mid-octave SDF values. The
 *                   memory for the matrix is automatically allocated
 *                   within the function. The universal matrix will be
 *                   a single-column of type MUTIL\_DOUBLE. In the case
 *                   where n\_sample is negative or equal to zero, the
 *                   result matrix will have the same number of
 *                   elements as does the levels matrix.  For the
 *                   case where n\_sample is positive, the length of
 *                   the result vector will be equal to J+1 where
 *                   J is the wavelet transform decomposition level.
 *
 * @see wavuniv_fdp_estimator_instantaneous
 * @see wavuniv_variance_confidence
 * @see wavuniv_variance_edof
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_fdp_bandpass_variance(
  const univ_mat *levels,
  double          delta,
  sint32          n_sample,
  void           *intrp_ptr,
  univ_mat       *result );


#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_FDP_H_*/

