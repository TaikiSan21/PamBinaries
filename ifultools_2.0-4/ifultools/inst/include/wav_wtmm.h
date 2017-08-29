
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_wtmm.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_WAV_WTMM_H_
#define IN_WAV_WTMM_H_

#include "mat_type.h"
#include "ut_err.h"
#include "ut_plat.h"
#include "ut_type.h"
#include "wav_type.h"

/* This file contains function declarations the wavelet mutils library. */

#ifdef __cplusplus
extern "C" {
#endif


/** The modulus maxima of a continuous wavelet transform.
 * Returns the index locations in time and scale of the continuous wavelet
 * transform modulus maxima. A point in the CWT X(t,j) is defined
 * as a maximum if |X(t-1,j)| + tol < |X(t,j)| and |X(t+1,j)| + tol < |X(t,j)|
 * where tol is a (scale-dependent) tolerance specified by the user.
 * The search algorithm is also adpated to identify plateaus in the data,
 * and will select the the middle of the plateau as a maximum location
 * when encountered. The data |X(t,j)| is first scaled so that its
 * maximum value is 1.0, so the tolerances should be adjusted accordingly.
 * Since the CWT coefficients are in effect a result band-pass filtering operations,
 * the large scale coefficients form a smoother curve than do the small
 * scale coefficients. Thus, the tolerance vector allows the user to specify
 * scale-dependent tolerances, helping to weed out undesirable local maxima.
 * It is recommended that the tolerance be set proportional to the scale,
 * e.g., tolerance = C / sqrt( scale) where C is a constant 0 < C < 1.
 *
 * References:
 * J.F. Muzy, E. Bacry, and A. Arneodo., "The multifractal formalism revisited with wavelets.",
 * International Journal of Bifurcation and Chaos 4, 245-302 (1994).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_wtmm.h
 * @source wav\_wtmm.c
 * @library wavelets
 * @usage #err = wavuniv_transform_continuous_wavelet_modulus_maxima( &cwt, &tolerance, peaktype intrp_ptr, &itime, &iscale );#
 * @return Standard mutils error/OK code.
 * @param  cwt         Pointer to a pre-allocated universal matrix containing
 *                     a continuous wavelet transform. This matrix must be
 *                     of type MUTIL\_DCOMPLEX and of size [N x Ns] where N
 *                     is the length of the original time series and Ns is
 *                     the number of scales.
 * @param  tolerance   A pointer to a universal matrix of type MUTIL\_DOUBLE
 *                     containing Ns tolerances to use in finding local
 *                     CWT modulus maxima. The jth element defines the
 *                     tolerance to use in finding modulus maxima at
 *                     the jth scale of the CWT.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  itime       Pointer to a sint32 matrix which (upon return)
 *                     will contain a concatenated list of the index locations
 *                     in time where the CWT modulus maxima were found.
 *                     The memory for this matrix is allocated by the function.
 * @param  iscale      Pointer to a sint32 matrix which (upon return)
 *                     will contain a concatenated list of the index locations
 *                     in scale where the CWT modulus maxima were found.
 *                     The memory for this matrix is allocated by the function.
 * @param peaktype     Enum defining the peak: extrema, minima, extrema.
 *
 * @see wavuniv_transform_continuous_wavelet
 * @see wavuniv_transform_continuous_wavelet_modulus_maxima_tree
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_continuous_wavelet_modulus_maxima(
  const univ_mat        *cwt,
  const univ_mat        *tolerance,
  const wav_transform_peak peak_type,
  void                  *intrp_ptr,
  sint32_mat            *itime,
  sint32_mat            *iscale );

/** The modulus maxima tree of a continuous wavelet transform.
 * Converts a list of continuous wavelet transform modulus maxima
 * (WTMM) index locations in time and in scale into a list of
 * branches. Each branch represents a collection of WTMM that
 * correspond to the same ridge in the WTMM time-scale plane.
 * A coarse-to-fine scale strategy is used to identify the members
 * of each branch as follows: (i) a single WTMM at the coarsest
 * scale is selected as the start of a given branch,
 * (ii) the closest neighboring WTMM in time at the next finest
 * scale is then added to the branch, (iii) step ii is repeated
 * until the smallest scale is reached or an apparent break
 * occurs in the branch across scale, and (iv) steps i-iii are
 * repeated until all WTMM have been accounted. Branches are
 * allowed to wrap off one end of the time-scale plane in time
 * and back onto the other. Gaps encountered in the search for
 * additional branch members may be bridged by setting the
 * bridge\_gap boolean to TRUE. A branch is not grown unless
 * the nearest neighbor candidate at the next finest scale is
 * close in time to the last recorded branch member, where "close"
 * is defined as being less than the current scale of the neighbor
 * candidate. This means that the window in time for admissible
 * neighbor WTMM candidates (at the next finest scale) shrinks
 * proportionally with scale.
 *
 * References:
 * J.F. Muzy, E. Bacry, and A. Arneodo., "The multifractal formalism revisited with wavelets.",
 * International Journal of Bifurcation and Chaos 4, 245-302 (1994).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_wtmm.h
 * @source wav\_wtmm.c
 * @library wavelets
 * @usage #err = wavuniv_transform_continuous_wavelet_modulus_maxima_tree( &wtmm_time_index, &wtmm_scale_index, &wtmm_scale, n_scale, bridge_gaps, n_octave_min, wtmm_strength_min, intrp_ptr, &tree );#
 * @return Standard mutils error/OK code.
 * @param  wtmm_time_index   Pointer to a pre-allocated sint32 matrix containing
 *                           a concatenated list of the index locations
 *                           in time where the CWT modulus maxima were found.
 *                           The matrix must be a single column or single row
 *                           and must coordinate (element by element) with the
 *                           wtmm\_scale\_index input matrix.
 * @param  wtmm_scale_index  Pointer to a pre-allocated sint32 matrix containing
 *                           a concatenated list of the index locations
 *                           in scale where the CWT modulus maxima were found.
 *                           The matrix must be a single column or single row
 *                           and must coordinate (element by element) with the
 *                           wtmm\_time\_index input matrix.
 * @param  wtmm_scale        Pointer to a pre-allocated double matrix containing
 *                           the scales used in developing the original
 *                           continuous wavelet transform.
 * @param  n_sample          The number of samples in the original time series.
 * @param  bridge_gaps       A logical flag. If TRUE, an attempt is made to bridge
 *                           any gaps (across scale) that are encountered in linking
 *                           WTMM into a single branch. Often, these gaps are produced
 *                           as an artifact from setting too low a tolerance in the
 *                           \Ref{wavuniv_transform_continuous_wavelet_modulus_maxima}
 *                           function, which is used to develop the wtmm\_time\_index
 *                           wtmm\_scale\_index vectors.
 * @param  n_octave_min      A pruning factor for excluding non-persistent branches. If
 *                           a WTMM branch does not span this number of octaves, it is
 *                           excluded from the tree.
 * @param  wtmm_strength_min A pruning factor for excluding weak branches. A given WTMM
 *                           is excluded from the current branch if it less than or equal
 *                           to the product of the maximum WTMM at the current scale and
 *                           wtmm\_strength\_min. This parameter must be in the range [0,1].
 *                           If set to zero, all WTMM are allowed. If set to unity,
 *                           only the maximum WTMM at the current scale is allowed.
 * @param  intrp_ptr         Pointer for implementation of interrupt checking.
 * @param  tree              Pointer to a matrix set which (upon return) will
 *                           contain the WTMM branches. Each branch is returned as
 *                           a [5 x Nj] double matrix where Nj is the length of
 *                           branch j. The rows of each matrix correspond
 *                           to the (1) time indices, (2) scale indices, (3) times,
 *                           (4) scales, and (5) WTMM for a single branch.The memory
 *                           for this matrix set and its constitutive matrices is
 *                           allocated by the function.
 *
 * @see wavuniv_transform_continuous_wavelet
 * @see wavuniv_transform_continuous_wavelet_modulus_maxima
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_continuous_wavelet_modulus_maxima_tree(
  const sint32_mat   *wtmm_time_index,
  const sint32_mat   *wtmm_scale_index,
  const dcomplex_mat *cwt,
  const double_mat   *cwt_time,
  const double_mat   *cwt_scale,
  const boolean       bridge_gaps,
  const double        n_octave_min,
  const double        wtmm_strength_min,
  void               *intrp_ptr,
  mat_set            *tree );

#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_WTMM_H_*/
