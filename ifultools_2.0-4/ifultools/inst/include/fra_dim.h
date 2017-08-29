
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_dim.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_DIMENSION_H_
#define IN_FRA_DIMENSION_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "fra_type.h"

/* This file contains function declarations the mutils fractal library. */

#ifdef __cplusplus
extern "C" {
#endif

/** Delay embedding of a univariate time series.
 * Given a uniformly sampled discrete time series x[t]
 * for t = 0, ..., N - 1, an integer delay T,
 * and an embedding dimension E, a delay embedding
 * is created by forming a matrix whose columns are
 * defined by:
 *
 * x[t], x[t + T], x[t + 2*T]..., x[t + (E - 1)*T].
 *
 * The number of columns in the resulting matrix
 * is E and the number of rows is N - (E - 1)*T.
 * A delay embedding is used to emulate the the phase
 * space of a multivariate dynamical system when
 * only a single variable of that system is measured.
 *
 * References:
 *
 * 1.Takens, F.
 * ``On the numerical determination of the dimension of an attractor'',
 * Dynamical Systems and Vibrations, oecture Notes in Math. Vol. 1125,
 * Springer, Heidelberg, NY, 1985.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_dim.h
 * @source fra\_dim.c
 * @library fractal
 * @usage #err = frauniv_embed( &time_series, embedding_dimension, time_lag, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series.
 * @param embedding_dimension An integer representing the embedding dimension.
 * @param time_lag  An integer representing the delay between coordinates.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a universal matrix of type MUTIL\_DOUBLE
 *                    which (upon return) will contain the delay embedding.
 *                    The memory for this matrix is allocated within the
 *                    function.
 * @see frauniv_dimension_correlation_summation
 * @see frauniv_dimension_information
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_embed(
  const univ_mat *time_series,
  sint32          embedding_dimension,
  const sint32    time_lag,
  void           *intrp_ptr,
  univ_mat       *result );

/** Information dimension statistics.
 * The information dimension (D1) is one of an infinite
 * number of fractal dimensions of a chaotic system.
 * For generalalized fractal dimension estimations,
 * correlation integral
 * moments are determined as an average of the contents of
 * neighbohoods in the phase space of equal radius eps.
 * Using this approach. the information dimension
 * for a given embedding dimension E is estimated
 * via D1( E ) = <ln( p )> / ln( eps ) in the limit as eps approaches zero,
 * where eps is the radius of an E-dimensional hypersphere,
 * p is the density (also known as the mass fraction),
 * and <ln( p )> is the average Shannon information needed
 * to specify an arbitrary point in the phase space with
 * accuracy eps.
 *
 * Alternatively, the neighborhoods can be constructed with
 * variable radii but with constant density. The scaling behavior
 * of the average radii of these neighborhoods as a function of density
 * is then used to estimate the fractal dimensions.
 * In this function, we use this constant density approach
 * to calculate the statistics for estimating the information
 * dimension.
 *
 * For single variable time series, the
 * phase space is approximated with a delay embedding and
 * D1(E) is thus estimated over statistics gathered for
 * dimensions 1,...,E. For chaotic
 * systems, these statistics will `saturate' at a finite
 * embedding dimension, revealing both the (estimated)
 * information dimension and an appropriate
 * embedding dimension for the system. A linear regression
 * scheme should be to estimate the D1(E) using the
 * statistics returned by this function.
 *
 * References:
 *
 * 1. Peter Grassberger and Itamar Procaccia,
 * ``Measuring the strangeness of strange attractors'',
 * Physica 9D (1983) 189.
 *
 * 2. Holger Kantz and Thomas Schreiber,
 * ``Nonlinear Time Series Analysis'',
 * Cambridge University Press, 1997.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_dim.h
 * @source fra\_dim.c
 * @library fractal
 * @usage #err = frauniv_dimension_information( &time_series, embedding_dimension, time_lag, orbital_lag, n_density, distance_metric, max_neighbors, n_average, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series.
 * @param embedding_dimension An integer representing the maximal embedding dimension.
 * @param time_lag  An integer representing the delay between coordinates.
 * @param orbital_lag The number of points along the trajectory of the
 *                    current point that must be exceeded in order for
 *                    another point in the phase space to be considered
 *                    a neighbor candidate. This argument is used
 *                    to help attenuate temporal correlation in the
 *                    the embedding which can lead to spuriously low
 *                    information dimension estimates. The orbital lag
 *                    must be positive or zero.
 * @param n_density   The number of densities to investigate.
 *                    The densities will be
 *                    logarithmically distributed
 *                    over the range [ 1 / N, 1 ], where N is the number
 *                    of points in the embedding.
 * @param distance_metric The metric used to define the distance between
 *                    points in the embedding. The argument is of
 *                    type \Ref{_fra_distance_metric}.
 * @param max_neighbors The maximum number of neighbors to use in
 *                    in developing the kd-tree (used as a quick means
 *                    of finding nearest neighbors in the phase space).
 *                    Given a density p, the number of neighbors
 *                    to find is k ~ p * N where N is the number of
 *                    points to use to develop the kd-tree. If
 *                    k <= max\_neighbors, N is the maximum number
 *                    of points in the embedding. If, however,
 *                    k > max\_neighbors, then k is reset to
 *                    k = max\_neighbors
 *                    and N is reduced accordingly to achieve the
 *                    specified density. While this process reduces
 *                    the computational and memory burden, it increases
 *                    the variability of the result since we are only
 *                    using a subset of the all available points to develop
 *                    the kd-tree. Set max\_neighbors = 100 for typical
 *                    applications, higher if time and memory are not an issue.
 *                    Note: in all cases, all the points in the embedding
 *                    are considered as neighbor candidates
 *                    and it is only the kd-tree size (and number of
 *                    neighbors to find) that are reduced.
 * @param n_reference The number of points in the embedding over which the
 *                    average neighborhood radius should be estimated
 *                    given a specified density.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a matrix set which (upon return) will
 *                    contain the information dimension statistics.
 *                    The memory for the
 *                    matrix set is automatically allocated within the
 *                    function and will containing two universal
 *                    matrices of type MUTIL\_DOUBLE.
 *
 *                    The first matrix has dimension P x E, where
 *                    P is the number of densities and E
 *                    is the maximal embedding dimension, and
 *                    contains the log2 of the average neighborhood radii.
 *
 *                    The second matrix is a single-column vector of double
 *                    values containing the log2 of densities
 *                    used to calculate the average radii. For
 *                    convenience, these densities will be
 *                    logarithmically distributed
 *                    over the range [ 1 / N, 1 ], where N is the number
 *                    of points in the embedding.
 *
 * @see frauniv_embed
 * @see frauniv_dimension_correlation_summation
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_dimension_information(
  const univ_mat *time_series,
  const sint32    embedding_dimension,
  const sint32    time_lag,
  const sint32    orbital_lag,
  const sint32    n_density,
  const fra_distance_metric distance_metric,
  const sint32    max_neighbors,
  const sint32    n_reference,
  void           *intrp_ptr,
  mat_set        *result );

/** Discrete correlation integral estimation using a delay embedding
 * of a single variable time series.
 *
 * This function calculates the correlation summations over
 * all available scales and for a range of finite
 * embedding dimensions. Following Taken's Theorem, the
 * embedding is created using a single variable time series.
 * Given the time series (X[t]), the embedding dimension (emb\_dim),
 * and the time delay (time\_lag), the embedding coordinates are
 * defined as X[t], X[t + time\_lag], ... , X[t + (emb\_dim - 1)*time\_lag].
 * For a given scale (eps) and embedding dimension (E), the
 * correlation summation C2(eps,E) is formed by calculating the
 * fraction of all points in the phase space that are separated
 * by a distance less than or equal to eps. The algorithm is
 * made computationally efficient by
 *
 * (1) using the L-infinity norm to calculate
 * the distance between neighbors in the phase space as
 * opposed to (say) the L2 norm which involves taking
 * computationally intense square root and power of two
 * operations. The L-infinity norm of the distance between
 * two points in the phase space is the absolute value of
 * the maximal difference between any of the points' respective
 * coordinates.
 *
 * (2) using bitwise masking and shift operations to reveal
 * the radix-2 exponent of the L-infinity norm. This direct
 * means of obtaining the exponent immediately yields the
 * asociated scale of the distance between neighbors in the
 * phase space while avoiding costly log operations. The bitwise
 * mask and shift factors are based on the IEEE standard 754
 * for binary floating-point arithmetic. Initial tests are
 * performed in the code to verify that the current machine
 * follows this standard.
 *
 * (3) using a computationally efficient routine to calculate
 * the resulting value of a float raised to a positive integer power.
 * The L-infinity norm is raised to an integer power (p) to
 * effectively increase the usable scales (spatial resolution)
 * by a factor of p.
 *
 * References:
 *
 * 1. Peter Grassberger and Itamar Procaccia,
 * ``Measuring the strangeness of strange attractors'',
 * Physica 9D (1983) 189.
 *
 * 2. Holger Kantz and Thomas Schreiber,
 * ``Nonlinear Time Series Analysis'',
 * Cambridge University Press, 1997.
 *
 * 3. Peter Grassberger and Itamar Procaccia,
 * ``Characterization of strange attractors'',
 * Physical Rewiew Letters, Vol.50(5), 1983.
 *
 * @limits Currently, no efficient nearest neighbor algorithm is used to
 * speed up the algorithm for small scale computations. Only delay coordinate
 * embeddings are supported for this routine.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_dim.h
 * @source fra\_dim.c
 * @library fractal
 * @usage #err = frauniv_dimesion_correlation_summation( &time_series, emb_dim, time_lag, orbital_lag, resolution, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param data            Pointer to a pre-allocated (a) single-row or
 *                        single-column universal matrix of type
 *                        MUTIL\_DOUBLE containing the time series or
 *                        (b) an N x E universal matrix of type
 *                        MUTIL\_DOUBLE containing the embedding where
 *                        N is the number of points and
 *                        E is the dimension of the embedding.
 *                        If the latter case, the emb\_dim and
 *                        time\_lag inputs are ignored.
 *
 * @param emb_dim         The maximal embedding dimension.
 * @param time_lag        The time lag used to create the delay coordinates
 *                        using the original time series. The time lag
 *                        must be positive.
 * @param orbital_lag  The number of points along the trajectory of the
 *                        current point that must be exceeded in order for
 *                        another point in the phase space to be considered
 *                        a neighbor candidate. This argument is used
 *                        to help attenuate temporal correlation in the
 *                        the embedding which can lead to spuriously low
 *                        correlation dimension estimates. The orbital lag
 *                        must be positive or zero.
 * @param resolution      An integer representing the spatial resolution factor.
 *                        A value of P increases the number of effective
 *                        scales by a factor of P at a cost of raising
 *                        the L-infinity norm to the Pth power.
 *                        For example, setting the resolution to 2
 *                        will double the number of scales while imposing
 *                        and additional multiplication operation. The
 *                        resolution must exceed unity.
 * @param intrp_ptr       Pointer for implementation of interrupt checking.
 * @param result          Pointer to a matrix set which (upon return) will
 *                        contain the correlation results. The memory for the
 *                        matrix set is automatically allocated within the
 *                        function and will containing two universal
 *                        matrices of type MUTIL\_DOUBLE.
 *
 *                        The first matrix has dimension [N x emb\_dim], where
 *                        N is the number of effective scales and emb\_dim
 *                        is the maximal embedding dimension, and
 *                        contains the log base-2 of the correlation summations.
 *                        The number of effective scales N is both a function
 *                        of the data and the resolution. value greater than
 *                        one indicates that no data exists at the corresponding
 *                        scale and embedding dimension.
 *                        The second matrix is a [N x 1] vector containing the
 *                        log base-2 of the effective scales.
 * @see frauniv_embed
 * @see frauniv_dimension_information
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_dimension_correlation_summation(
 const univ_mat *time_series,
 sint32          emb_dim,
 sint32          time_lag,
 sint32          orbital_lag,
 sint32          n_scale_octave,
 void           *intrp_ptr,
 mat_set        *result );


/** Estimation of the proper embedding dimension for a single-variable time series.
 * The method of False Nearest Neighbors (FNN) is a direct means of
 * estimating the proper embedding dimension using a delay embedding scheme
 * for a single-variable time series. The method works by finding the nearest
 * admissible neighbor for each point of the data embedded in a given dimension.
 * The data is then embedded in the next dimension and if the increase in
 * distance between the original neighbors grows sufficiently then those
 * neighbors are "false neighbors". As the base
 * embedding dimension increases, it is expected that the percentage of FNN
 * will approach zero for deterministic series. The base dimension in which this
 * occurs is seen as the proper embedding dimension for the time series.
 *
 * The statistics used for determining a FNN are based on two Euclidean
 * tolerances supplied by the user, namely, rtol and atol.
 *
 * Test 1:
 *
 * Let R(d) be the
 * Euclidean distance between two nearest neighbors in embedding dimension d,
 * and dR(d+1) be the Euclidean distance gained by embedding the neighbors in
 * dimension d + 1. These neighbors are declared false if dR(d+1) / R(d) > rtol.
 *
 * Test 2:
 *
 * Let A be the estimated size the attractor in one dimension. Two neighbors
 * are in embedding dimension d are declared false if R(d + 1) / A > atol.
 *
 * Highly correlated time series can induce spurious results on the FNN method
 * by driving the percentage of FNN toward zero with an increase in embedding
 * dimension, and may do so even if the series is a realization of a
 * highly correlated random process (such as a random-walk). To avoid this
 * pitfall, we define an "admissible neighbor candidate" as one in which
 * is close spatially but not temporally. The minimum separation time for admissible
 * neighbors is supplied by the user through the orbital\_lag input argument.
 *
 * References:
 *
 * 1. M. B. Kennel, R. Brown, and H. D. I. Abarbanel (1992),
 * "Determining embedding dimension for phase-space reconstruction
 * using a geometrical construction",
 * Physical Review A, 45(6), 3403-3411.
 *
 * 2. Fredkin, D. R., and Rice, J. A. (1995),
 * "Method of false nearest neighbors: A cautionary note",
 * Physical Review E, 51(4), 2950-2954.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_dim.h
 * @source fra\_dim.c
 * @library fractal
 * @usage #err = frauniv_dimension_false_nearest_neighbors( &time_series, embedding_dimension, time_lag, orbital_lag, rtol, atol, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series.
 * @param embedding_dimension An integer representing the maximal embedding dimension.
 * @param time_lag  An integer representing the delay between coordinates.
 * @param orbital_lag The number of points along the trajectory of the
 *                    current point that must be exceeded in order for
 *                    another point in the phase space to be considered
 *                    a neighbor candidate. This argument is used
 *                    to help attenuate temporal correlation in the
 *                    the embedding which can lead to spuriously low
 *                    FNN estimates. The orbital lag
 *                    must be positive or zero.
 * @param rtol        The rtol tolerance (see above for details).
 * @param atol        The atol tolerance (see above for details).
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a universal matrix of type MUTIL\_DOUBLE
 *                    and of size [3 x embedding\_dimension]. The first
 *                    two rows contain the percentage of FNN based on the atol
 *                    and rtol criterion, respectively. The last row is a
 *                    the percentage of FNN based on the joint (atol and rtol)
 *                    tests. The joint test states that two nearest neighbors are
 *                    declared as false neighbors if either the atol or rtol tests fail.
 *
 * @see frauniv_embed
 * @see frauniv_dimension_information
 * @see frauniv_dimension_correlation_summation
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_dimension_false_nearest_neighbors(
  const univ_mat *time_series,
  const sint32    embedding_dimension,
  const sint32    time_lag,
  const sint32    orbital_lag,
  const double    rtol,
  const double    atol,
  void           *intrp_ptr,
  univ_mat       *result );


/** Estimation of the proper embedding dimension for a single-variable time series.
 * The method of False Nearest Strands (FNS) is a direct means of
 * estimating the proper embedding dimension using a time delayed embedding scheme
 * for a single-variable time series. The method works by finding the nearest
 * admissible neighbor for each point of the data embedded in a given dimension.
 * [original, neighbor] index pairs are formed and agglomerated into strands.
 * Each is labeled either as a false or true strand based on a 
 * Euclidean distance metric.
 * If the percentage of false strands remains sufficiently small at and
 * beyond embedding dimension E, then E is considered to be a sufficient
 * embedding dimension for the data. This method is an improvement to
 * the False Nearest Neighbors method, which can produce spuriously
 * low dimension estimates due to high temporal correlation in the data.
 * Linking the nearest neighbor pairs into strands attenuates this
 * adverse effect.
 *
 * The statistic used for determining a false nearest strand (FNS)
 * is based on a Euclidean
 * tolerance supplied by the user (atol).
 * Let S(d) be the mean Euclidean distance in the projected (d+1)th coordinate 
 * between strand pairs found to be nearest neighbors in embedding dimension d. 
 * If S(d) / Ra > atol, where Ra is the estimated attractor size, then
 * the strand is considered to be a false strand. Ra is typically calculated
 * to be the sample standard deviation of the original time series and atol
 * is typically chosen to be 1 or 2. The S(d) metric is a measure of the average
 * additional Euclidean distance we gain by embedding the strand in the next
 * dimension, and the FNS statistic is used to assess when this extra distance
 * has grown too large, indicating a false strand. 
 *
 * References:
 *
 * 1. M. B. Kennel and Henry D.I. Abarbanel (2002),
 * "False neighbors and false strands: A reliable minimum embedding dimension
 * algorithm", Physical Review E, (66), 026209, 1--19.
 *
 * 2. M. B. Kennel, R. Brown, and H. D. I. Abarbanel (1992),
 * "Determining embedding dimension for phase-space reconstruction
 * using a geometrical construction",
 * Physical Review A, 45(6), 3403-3411.
 *
 * 3. Fredkin, D. R., and Rice, J. A. (1995),
 * "Method of false nearest neighbors: A cautionary note",
 * Physical Review E, 51(4), 2950-2954.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_dim.h
 * @source fra\_dim.c
 * @library fractal
 * @usage #err = frauniv_dimension_false_nearest_strands( &time_series, embedding_dimension, time_lag, orbital_lag, iterate_tolerance, atol, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series.
 * @param embedding_dimension An integer representing the maximal embedding dimension.
 * @param time_lag  An integer representing the delay between coordinates.
 * @param orbital_lag The number of points along the trajectory of the
 *                    current point that must be exceeded in order for
 *                    another point in the phase space to be considered
 *                    a neighbor candidate. This argument is used
 *                    to help attenuate temporal correlation in the
 *                    the embedding. The orbital lag
 *                    must be positive or zero.
 * @param iterate_tolerance  An integer defining the so-called iterate tolerance. Nearest
 *    neighbor pairs (i,J(i)) are separated in time by a point index span dindex = |i-J(i)|,
 *    where J(i) represents the index of the nearest neighbor to point i. If a point near i,
 *    say k points away also has a nearest neighbor such that |k - J(k)| = dindex +/- M,
 *    where M is the iterate tolerance, then the pair (k, J(k)) is added to the current
 *    strand. Typically, M=0 or M=1. If M=0, then the difference in index must be exactly
 *    the same for each pair included in the strand. If M=1, the index difference is allowed
 *    to be 1 point off from the reference pair.  
 * @param atol        The atol tolerance (see above for details).
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a universal matrix of type MUTIL\_DOUBLE
 *                    and of size [embedding\_dimension x 1] containing the
 *                    FNS percentages as a function of embedding dimension. The first
 *                    element corresponds to the FNS percentage for dimension 1.
 *
 * @see frauniv_embed
 * @see frauniv_false_nearest_neighbors
 * @see frauniv_dimension_information
 * @see frauniv_dimension_correlation_summation
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_dimension_false_nearest_strands(
  const univ_mat *time_series,
  const sint32    embedding_dimension,
  const sint32    time_lag,
  const sint32    orbital_lag,
  const sint32    iterate_tolerance,
  const double    atol,
  void           *intrp_ptr,
  univ_mat       *result );


/** Poincare map.
 * This function creates a map from a given scalar time series, s(t),
 * assumed to have underlying non-linear dynamical structure. It is
 * known that the coordinates s(t), s'(t), s''(t), ... may be used to
 * reconstruct the underlying phase space dynamics, i.e., the flow. A Poincare
 * surface of sections may then be formed by the surface given by s'(t) = 0.
 * The collection of intersections of the flow trajectories in a chosen
 * direction with this surface of section form an N-1 dimensional map, where
 * N is the dimension of the phase space. This so-called Poincare map may
 * subsequently be used in a delay embedding to visualize the underlying
 * dynamics or to estimate invariant measures of the underlying dynamics,
 * e.g., correlation dimension and Lyapunov exponents.
 *
 * @algorithm
 *
 * If the denoise parameter is TRUE, waveshrink is applied to the original
 * series in order to attenuate Gaussian white noise contamination.
 * The first and second derivatives of the resulting series
 * are approximated via the continuous
 * wavelet transform using the first derivative of a Gaussian as a mother
 * wavelet filter (see references for details). The locations of the
 * local extrema are then estimated using the standard first and second
 * derivative tests on the CWT coefficients
 * at a single and appropriate scale (an appropriate scale is one that is
 * large enough to smooth out noisy components but not so large as to
 * the oversmooth the data). The extrema locations are then fit with a quadratic
 * interpolation scheme to estimate the amplitude of the extrema using the
 * (possibly waveshrunk) original time series.
 *
 * References:
 *
 * 1. Nie L, Wu S, Lin X, Zheng L, Rui L.,
 * ``Approximate derivative calculated by using continuous wavelet transform'',
 *  J. Chem. Inf. Comput. Sci. (2002) Mar-Apr;42(2):274-83.
 *
 * 2. Holger Kantz and Thomas Schreiber,
 * ``Nonlinear Time Series Analysis'',
 * Cambridge University Press, 1997.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_dim.h
 * @source fra\_dim.c
 * @library fractal
 * @usage #err = frauniv_poincare_map( &time_series,FRA_EXTREMA_MIN, FALSE, intrp_ptr, &locs, &vals );#
 * @return Standard mutils error/OK code.
 * @param  time_series   Pointer to a pre-allocated universal matrix containing
 *                       the time series to analyze. This input must be a row
 *                       or column vector of type MUTIL\_DOUBLE.
 * @param  ex_type       An enumerated value of type \Ref{_fra_extrema_type}
 *                       which determines the type of extrema used to form the
 *                       resulting Poincare map. FRA\_EXTREMA\_MINIMA uses the
 *                       time series minima only, FRA\_EXTREMA\_MAXIMA uses only
 *                       the time series maxima and FRA\_EXTREMA\_ALL uses both
 *                       minima and maxima.
 * @param  denoise       A boolean indicating whether the input time series
 *                       will be denoised prior to finding its extrema. If
 *                       TRUE then the WaveShrink algorithm will be performed
 *                       on the input time series.
 * @param  intrp_ptr     Pointer for implementation of interrupt checking.
 * @param  ex_locs       Pointer to a universal matrix of type MUTIL\_DOUBLE
 *                       which upon return will contain the abscissa locations
 *                       of the found extrema. This assumes the abscissa of
 *                       the input time series is 0, 1, 2, ... The orientation
 *                       of the matrix (row or column vector) will match that
 *                       of the time series input. Memory for this matrix is
 *                       allocated within the function.
 * @param  ex_vals       Pointer to a universal matrix of type MUTIL\_DOUBLE
 *                       which upon return will contain the resulting Poincare
 *                       map, i.e., the interpolated time series values at the
 *                       abscissa locations given in {\bf ex\_locs}. The
 *                       orientation of the matrix (row or column vector) will
 *                       match that of the time series input. Memory for this
 *                       matrix is allocated within the function.
 * @see _fra_extrema_type
 * @see wavuniv_shrink
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_poincare_map(
  const univ_mat             *time_series,
  const fra_extrema_type      extrema_type,
  const boolean               denoise,
  void                       *intrp_ptr,
  univ_mat                   *extrema_location,
  univ_mat                   *extrema_amplitude );

/** Space time separation plot.
 * This function takes a time series input and computes contour lines as
 * viewed in a space time separation (STS) plot. Each contour, C(p,.), in
 * an STS plot corresponds to a particular probability, p, and gives spatial
 * distance as a function of temporal separation between pairs of points in
 * the phase space. In particular, any pair of points seperated in time by
 * delta_t are separated by C(p,delta_t) in spatial distance with
 * probabilty p.
 *
 * The inputs to the funtion are the time series to be embedded, along with
 * the dembedding dimension and the desired delay. In addition, a matrix of
 * orbital lags (integers representing time separations in phase space) and
 * a fractional probability (frac_prob) at which to compute the first
 * countour. An M x N matrix is returned where N is the number of orbital
 * lags given. M is the largest integer such that M*frac_prob <= 1. Hence,
 * smaller frac_prob produces more contour lines. Each row of the output
 * matrix may be plotted versus the given time separations to produce the
 * STS plot.
 *
 * NOTE: All distances are computed using the infinity norm.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_visu.h
 * @source fra\_visu.c
 * @library fractal
 * @usage #err = frauniv_space_time_separation_plot(&ts,edim,dly,delta_t,0.01,int_ptr,&plotdata);#
 * @return Standard mutils error/OK code.
 * @param time_series  Pointer to a universal matrix of type MUTIL_DOUBLE,
 *                     which contains a row or column vector representing a
 *                     time series.
 * @param dim          An integer representing the embedding dimension.
 * @param delay        An integer representing the delay between coordinates.
 * @param orbital_lags A pointer to a row or column universal matrix of type
 *                     MUTIL_SINT32 specifying the time separations at which
 *                     to compute the countours;
 * @param frac_prob    A positive number between 0 and 1 which gives the
 *                     probability represented by the first countour. The
 *                     smaller this number the more countour lines the will
 *                     be generated.
 * @param intrp_ptr    Parameter used for user interrupts.
 * @param max_eps      Pointer to a double. For the user's convenience the
 *                     maximum possible distance between any pair of embedding
 *                     vectors is return through this parameter.
 * @param result       Pointer to an unallocated universal matrix. Upon
 *                     return, an N x M  matrix of type MUTIL\_DOUBLE which
 *                     holds the computed contours. Here, M is the number of
 *                     elements in input orbital_lags and N is the number of
 *                     generated countours (see details above).
 * @see frauniv_embed
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_space_time_separation_plot(
  const univ_mat   *time_series,
  const sint32      dim,
  const sint32      delay,
  const univ_mat   *orbital_lags,
  const double      frac_prob,
  void             *intrp_ptr,
  double           *max_eps,
  univ_mat         *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_DIMENSION_H_ */
