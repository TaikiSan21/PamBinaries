
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_modl.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_MODL_H_
#define IN_FRA_MODL_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "fra_type.h"

/* This file contains function declarations the mutils fractal library. */

#ifdef __cplusplus
extern "C" {
#endif



/** Detecting determinism in a time series. The delta-epsilon
 * test for detecting
 * deterministic structure in a time series works by exploiting
 * the continuity of the data in a delay coordinate embedding.
 * This structure is non-existent in randomized versions of
 * the original time series (surrogate data) and differences
 * of an appropriate discriminating statistic between an
 * ensemble of surrogate data and the original time series
 * indicates the presence of determinstic structure in the
 * the original time series.
 * The discriminating statistic is calculated as follows:
 * \TEX{
 * \begin{eqnarray}
 * \delta _{j,k} = |z_{j} - z_{k}| \\
 * \epsilon _{j,k} = |z_{j+ \kappa} - z_{k+\kappa}| \\
 * e(r) \equiv \overline{\epsilon _{j,k}} \qquad\hbox{for $j,k$ s.t. } r
 * \leq \delta_{j,k} < r + \Delta r
 * \end{eqnarray}
 * }
 *
 * where $\Delta r$ is the Euclidean bin size.  A cumulative average of $e(r)$ is formed by
 * \TEX{
 * \[ E(r) \equiv \sum \overline{e(r)}. \]
 * }
 *
 * If there exists a distinct separation of $E(r)$ from surrogate data, it implies that the signal
 * is deterministic.
 * The image separation $\kappa$ should be chosen large enough to sufficiently
 * decorrelate the points along an orbit, i.e. $\kappa$ should be considered as
 * an orbital lag as no points $|j - k| < \kappa$ are considered in the calculations.
 *
 * References:
 *
 * 1.Kaplan, D.,
 * ``Exceptional Events as Evidence for Determinism'',
 * Physica D 73 (1994) 38-48.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_modl.h
 * @source fra\_modl.c
 * @library fractal
 * @usage #err = frauniv_determinism_delta_epsilon( &time_series, embedding_dimension, time_lag, orbital_lag, scale_min, scale_max, resolution, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a pre-allocated single-row or
 *                    single-column universal matrix of type MUTIL\_DOUBLE
 *                    containing the time series.
 * @param embedding_dimension An integer representing the maximal embedding
 *                    dimension.
 * @param time_lag    An integer representing the delay between coordinates.
 * @param orbital_lag The number of points along the trajectory of the
 *                    current point that must be exceeded in order for
 *                    another point in the phase space to be considered
 *                    a neighbor candidate. This argument is used
 *                    to help attenuate temporal correlation in the
 *                    the embedding which can lead to spuriously low
 *                    FNN estimates. The orbital lag
 *                    must be positive or zero.
 * @param scale_min   The miminum scale (delta).
 * @param scale_max   The maximum scale (delta).
 * @param resolution  The separation between adjacent scales.
 * @param minimize    A boolean flag. If TRUE, only the nonzero portions of
 *                    the epsilon and cumsum epsilon matrices are returned,
 *                    essentially trimming the data of superfluous scales.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a matrix set containing three universal
 *                    matrices of type MUTIL\_DOUBLE.
 *                    The first is a single-column matrix of length N
 *                    containing the scales definining the bin boundaries.
 *                    The second and third matrices are of size
 *                    N x embedding\_dimension whose columns contain
 *                    epsilon and the cumulative sum of the epislon statistic,
 *                    respectively. The rows of these last two matrices
 *                    correspond to the scales vector. If minmize is FALSE,
 *                    N is the number of scales. Otherwise, N is data dependent,
 *                    and is on the range of [1, n_scale]. In this case,
 *                    The lowest scale returned is the smallest scale
 *                    in which there was at least one delta recorded in
 *                    dimension 1. The highest scale returned is the
 *                    maximum scale at which there was at least one delta
 *                    recorded for dimension embedding\_dimension.
 *                    The memory for the matrix
 *                    set is automatically allocated within the function.
 *
 * @see frauniv_surrogate
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_determinism_delta_epsilon(
 const univ_mat *time_series,
 const sint32    embedding_dimension,
 const sint32    time_lag,
 const sint32    orbital_lag,
 const sint32    image_lag,
 const double    scale_min,
 const double    scale_max,
 const double    resolution,
 const boolean   minimize,
 void           *intrp_ptr,
 mat_set        *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_MODL_H_ */

