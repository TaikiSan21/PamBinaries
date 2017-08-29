
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_dwtc.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_WAV_DWTC_H_
#define IN_WAV_DWTC_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "wav_type.h"

/* This file contains function declarations the wavelet mutils library.
*/

#ifdef __cplusplus
extern "C" {
#endif


/** The discrete wavelet transform using convolution style filtering.
 * This is a restricted version of \Ref{wavuniv_transform_discrete_wavelet}.
 * Given j, t are the decomposition level, and time index,
 * respectively, and V(0,t) = X(t) where X(t) is an N point time
 * series, and $N_j = \lfloor N / 2^j \rfloor$, the convolution style
 * DWT is given by
 *
 * \begin{center}
 * \TEX{
 * \begin{eqnarray}
 *  W_{j,t} & \equiv & \sum_{l=0}^{L-1} h_l V_{j-1, 2t + 1 - l\;
 *     \mbox{ mod }N_{j-1}}\nonumber\\
 *  V_{j,t} & \equiv & \sum_{l=0}^{L-1} g_l V_{j-1, 2t + 1 - l\;
 *     \mbox{ mod }N_{j-1}}\nonumber\\\nonumber
 * \end{eqnarray}
 * }
 * \end{center}
 *
 * for t = 0, ..., Nj - 1. The variable L is the length of both the
 * scaling filter (g) and wavelet filter (h).  W(j,t) are the wavelet
 * coefficients at decomposition level j and time index t while V(j,t)
 * are the scaling coefficients.
 *
 * A note on odd length time series:
 *
 * Odd length time series are allowed, but require a special (ad hoc)
 * method in order to preserve certain statistical properties of the
 * transform.  The method used here is energy preserving (and is
 * therefore usable for wavelet variance estimation). The method works
 * by preserving (at most) one scaling coefficient per decomposition
 * level. Specifically, the last scaling coefficient at a given level
 * is preserved if the number of scaling coefficients (in that level)
 * is odd. For simplicity, the collection of preserved scaling
 * coefficients are appended to the (returned) transform
 * coefficients. The length of the original series is also preserved
 * so that for an N point time series, N transform coefficients are
 * returned. For dyadic length sequences, there are no preserved
 * scaling coefficients because at every level the number of scaling
 * coefficients is even.
 *
 * For more information about the algorithm, see D.B. Percival and
 * A.T. Walden, ``Wavelet Methods for Time Series Analysis,''
 * Cambridge University Press, 2000.
 *
 * @limits Only relevant for the DWT as defined in the above
 * reference. Periodic extension is used to handle boundary
 * conditions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwtc.h
 * @source wav\_dwtc.c
 * @library wavelets
 * @usage #err = wavuniv_transform_discrete_wavelet_convolution( &time_series, &filters, n_level, intrp_ptr, &result);#
 * @return Standard mutils error/OK code.
 * @param  time_series      Pointer to universal matrix
 *                          containing the time series to analyze.
 *                          This input must have only one column and be
 *                          of type MUTIL\_DOUBLE.
 * @param  filters          Pointer to a matrix set containing
 *                          two pre-allocated universal matrices
 *                          of type MUTIL\_DOUBLE and containing
 *                          (in order) the wavelet and scaling
 *                          filter coefficients. The filter matrices
 *                          must be a single-column or single-row
 *                          and each must contain the same number of elements.
 * @param  n_level          The number of decomposition
 *                          levels. This input should not
 *                          exceed $2^{\lfloor \log_2{N}\rfloor}$ where N
 *                          is the length of the time series.
 * @param  intrp_ptr        Pointer for implementation of interrupt checking.
 * @param  result           Pointer to matrix set which (upon return)
 *                          will contain universal matrices of
 *                          type MUTIL\_DOUBLE to store transform coefficients.
 *                          The number of matrices returned in the matrix set
 *                          is dependent upon whether or not there are any
 *                          extra scaling coefficients to return. In the case
 *                          where there are no extra scaling coefficients,
 *                          the number of matrices is equal to n\_levels + 1.
 *                          In the case where there are extra scaling
 *                          coefficients to be returned, the number of
 *                          matrices is equal to n\_levels + 2.
 *                          The order of matrices in the matrix set is
 *                          given by W\_1,  W\_2, ... , W\_J , V\_J,
 *                          V\_{preserved}],
 *                          where W\_j are the N\_j wavelet coefficients at
 *                          level j,  V\_J are the N\_J scaling coefficients
 *                          at level J = num\_levels, and V\_{preserved} are
 *                          the collection of preserved scaling coefficients
 *                          stored in order of increasing scale (the number
 *                          of preserved coefficients can range anywhere
 *                          from 0 to J+1).
 *
 * @see _wav_filter_type
 * @see wavuniv_transform_discrete_wavelet
 * @see wavuniv_transform_discrete_wavelet_inverse
 * @see wavuniv_transform_discrete_wavelet_convolution_inverse
 * @see wavuniv_coefficient_zero_phase
 * @see wavuniv_coefficient_boundaries
 * @see wavuniv_transform_discrete_wavelet
 * @see wavuniv_transform_maximum_overlap
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_discrete_wavelet_convolution(
  const univ_mat *time_series,
  const mat_set  *filters,
  sint32          n_level,
  void           *intrp_ptr,
  mat_set        *result );

/** The discrete wavelet packet transform using convolution style filtering.
 * The discrete wavelet packet transform (DWPT) is regarded as any
 * one of a collection of orthonormal transforms, each of which can be readily
 * computed using a simple modificaiton to the DWT scheme. Each DWPT
 * is associtaed with a level j, and the jth level DWPT decomposes the
 * normalized frequency interval [0,1/2] into $2^j$ equal and
 * individual intervals. Each interval is referred to as a "crystal"
 * (or "node") denoted by the variable "W".
 * The location of a particular crystal is given by the so-called "natural order" of a
 * wavelet packet tree (j,n), where where j is the decomposition level and n is
 * the local node index. By definition, the original time series is located at
 * crystal (0,0) in the wavelet packet tree. The first level MODWPT contains
 * crystals (1,0) and (1,1), the second level contains crystals (2,0),(2,1),(2,2)
 * and (2,3), and so on. Each level of a wavelet packet transform
 * contains exactly $2^j$ crystals with the local node index range given by
 * $0 \le n \le 2^j - 1$. The wavelet transform coefficients at crystal (j,n)
 * are nominally associated with the normalized frequency band
 * $f \in [n/2^{j+1}, (n+1)/2^{j+1}]$. For example, crystal (2,3)
 * nominally corresponds to frequencies [3/8, 1/2]. Thus, an increase in the
 * local node index (n) corresponds to an increase in frequency content.
 *
 * Let g and h be the scaling fiter and wavelet filter respectively.
 * Given j, n, t are the decomposition level, local  node index, and time index,
 * respectively, and W(0,0,t) = X(t) where X(t) is an N point time
 * series, and $N_j = \lfloor N / 2^j \rfloor$, the convolution style
 * DWPT is given by
 *
 * \[ W_{j,n,t} \equiv \sum_{l=0}^{L-1} u_{n,l}
 * W_{j-1, \lfloor n/2 \rfloor, 2t + 1 - l \bmod N_{j-1}} \]
 *
 * where
 *
 * \begin{eqnarray}
 *   u_{n,l} \equiv &g_l& \mbox{\hspace{1cm}if $n \bmod 4 = 0$ or $3$}\nonumber \\
 * &h_l&  \mbox{\hspace{1cm} if $n \bmod 4 = 1$ or $2$}. \nonumber
 * \end{eqnarray}
 *
 * A note on odd length time series:
 *
 * Odd length time series are allowed, but require a special (ad hoc)
 * method in order to preserve certain statistical properties of the
 * transform.  The method used here is energy preserving (and is
 * therefore usable for wavelet variance estimation). The method works
 * by preserving (at most) one scaling coefficient per decomposition
 * level. Specifically, the last scaling coefficient at a given level
 * is preserved if the number of scaling coefficients (in that level)
 * is odd. For simplicity, the collection of preserved scaling
 * coefficients are appended to the (returned) transform
 * coefficients. The length of the original series is also preserved
 * so that for an N point time series, N transform coefficients are
 * returned. For dyadic length sequences, there are no preserved
 * scaling coefficients because at every level the number of scaling
 * coefficients is even.
 *
 * For more information about the algorithm, see D.B. Percival and
 * A.T. Walden, ``Wavelet Methods for Time Series Analysis,''
 * Cambridge University Press, 2000.
 *
 * @limits Periodic extension is used to handle boundary
 * conditions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwtc.h
 * @source wav\_dwtc.c
 * @library wavelets
 * @usage #err = wavuniv_transform_packet( &time_series, &filters, n_level, intrp_ptr, &result);#
 * @return Standard mutils error/OK code.
 * @param  time_series      Pointer to universal matrix
 *                          containing the time series to analyze.
 *                          This input must have only one column and be
 *                          of type MUTIL\_DOUBLE.
 * @param  filters          Pointer to a matrix set containing
 *                          two pre-allocated universal matrices
 *                          of type MUTIL\_DOUBLE and containing
 *                          (in order) the wavelet and scaling
 *                          filter coefficients. The filter matrices
 *                          must be a single-column or single-row
 *                          and each must contain the same number of elements.
 * @param  n_level          The number of decomposition
 *                          levels. This input should not
 *                          exceed $2^{\lfloor \log_2{N}\rfloor}$ where N
 *                          is the length of the time series.
 * @param  intrp_ptr        Pointer for implementation of interrupt checking.
 * @param  result      Pointer to a matrix set to store the result.
 *                     The matrix set header and matrices are automatically
 *                     allocated by this function and (upon return) will
 *                     contain $2^{n_level + 1 } - 1$ single-row universal
 *                     matrices of type MUTIL\_DOUBLE. The order of the
 *                     matrices is given by
 *                     [W\_{0,0} , W\_{1,0}, W\_{1,1}, W\_{2,0}, W\_{2,1},
 *                     W\_{2,2}, W\_{2,3}, ... , W\_{J,0}, ...  W\_{J, N\_J}]
 *                     where J  = n\_level and $N\_J = 2^J$. By definition,
 *                     W\_{0,0} is the original time series.
 *
 * @see _wav_filter_type
 * @see wavuniv_transform_discrete_wavelet
 * @see wavuniv_transform_maximum_overlap_packet
 * @see wavuniv_transform_packet_detail
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_packet(
  const univ_mat *time_series,
  const mat_set  *filters,
  sint32          n_level,
  void           *intrp_ptr,
  mat_set        *result );

/** The inverse discrete wavelet transform using convolution style
 * filtering.  This is a restricted version of
 * \Ref{wavuniv_transform_discrete_wavelet_inverse}.  Performs a
 * synthesis of the convolution style DWT. In the synthesis, the
 * so-called extra scaling coefficients (which are preserved in levels
 * where the number of scaling coefficients is odd) are added at each
 * stage of reconstruction.
 *
 * For more information about the algorithm, see D.B. Percival and
 * A.T. Walden, ``Wavelet Methods for Time Series Analysis,''
 * Cambridge University Press, 2000.
 *
 * @limits Only relevant for the DWT as defined in the above
 * reference. Periodic extension is used to handle boundary
 * conditions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwtc.h
 * @source wav\_dwtc.c
 * @library wavelets
 * @usage #err = wavuniv_transform_discrete_wavelet_convolution_inverse( &dwt, &filters, n_level, intrp_ptr, &result);#
 * @return Standard mutils error/OK code.
 * @param  dwt        Pointer to a pre-allocated universal matrix set
 *                    containing the convolution style DWT coefficients
 *                    stored in the [1 x N\_j] universal matrices:
 *                    [W\_1 | W\_2 | ... | W\_J | V\_J | V\_{preserved},
 *                    where W\_j are the N\_j wavelet coefficients at
 *                    level j, V\_J are the N\_J scaling coefficients
 *                    at level J = num\_levels, and V\_{preserved} are
 *                    the colleciton of preserved scaling coefficients
 *                    stored in order of increasing scale (the number
 *                    of preserved coefficients can range anywhere
 *                    from 0 to J+1). For certain cases (such as the
 *                    the DWT of a dyadic length sequence) there
 *                    are no extra scaling coefficients to return and
 *                    hence the V\_preseved matrix is not present
 *                    in the dwt matrix set.
 * @param  filters    Pointer to a matrix set containing
 *                    two pre-allocated universal matrices
 *                    of type MUTIL\_DOUBLE and containing
 *                    (in order) the wavelet and scaling
 *                    filter coefficients. The filter matrices
 *                    must be a single-column or single-row
 *                    and each must contain the same number of elements.
 * @param  intrp_ptr  Pointer for implementation of interrupt checking.
 * @param  result     Pointer to a universal matrix which (upon return)
 *                    will contain the synthesis. The return matrix will
 *                    be of of type MUTIL\_DOUBLE and size [1 x N] where
 *                    N is the number of points in the original time
 *                    series. The memory for the return matrix is
 *                    automatically allocated in the function.
 *
 * @see _wav_filter_type
 * @see wavuniv_transform_discrete_wavelet
 * @see wavuniv_transform_discrete_wavelet_inverse
 * @see wavuniv_transform_discrete_wavelet_convolution
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_discrete_wavelet_convolution_inverse(
  const mat_set  *dwt,
  const mat_set  *filters,
  void           *intrp_ptr,
  univ_mat       *result );

/** Conversion and validation of wavelet packet indices.
 * The (maximal overlap) discrete wavelet packet transform
 * contains many different transforms, each a result of a
 * projection of the original time series onto a unique
 * basis. The crystals which define an individiual transform
 * comprise a subset of all such crystals in a wavelet packet transform.
 * This function examines a 'list' of such crystals, checking
 * for (1) redundant entries, (2) whether the oscillation index
 * exceeds its limits at the current decomposition level, and
 * (3) if the normalized bandwidth of the union of all crystals
 * is 1/2 (otherwise, the subset represents only a partial projection
 * of the original time series onto the corresponding basis).
 * If all tests are successful, the transform indices are returned
 * in all the supported forms.
 *
 * For more information about the algorithm, see D.B. Percival and
 * A.T. Walden, ``Wavelet Methods for Time Series Analysis,''
 * Cambridge University Press, 2000.
 *
 * @limits Only relevant for the DWPT and MODWPT.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwtc.h
 * @source wav\_dwtc.c
 * @library wavelets
 * @usage #err = wavuniv_transform_packet_convert_indices( &transform_indices, intrp_ptr, &flat, &level, &osc );#
 * @return Standard mutils error/OK code.
 * @param  transform_indices  Pointer to a pre-allocated universal matrix
 *         of type MUTIL\_SINT32 containing the representative indices of the
 *         crystals which define the subset of a wavelet packet transform to analyze.
 *         The wavelet packet crystals are assumed to be arranged in the so-called
 *         natural order ala [W\_{0,0} , W\_{1,0}, W\_{1,1}, W\_{2,0}, W\_{2,1},
 *         W\_{2,2}, W\_{2,3}, ... , W\_{J,0}, ...  W\_{J, N\_J}]
 *         where J  is the number of decomposition levels and $N\_J = 2^J$.
 *         By definition, W\_{0,0} is the original time series. A given crystal
 *         is identified in the W\_{j,n} form above or by a flattened index.
 *         For example, the DWPT crystal in level 2 at oscillation index 1 is
 *         identified either by {j,n} = {2,1} or by its flattened index 4
 *         (with zero based indexing, 4 represents the fifth crystal of the
 *         wavelet packet transform in natural order). If the crystals
 *         are identified using flattened indices, then the transform\_indices
 *         matrix must be a single-column or single-row matrix. Conversely, if the
 *         crystals are identified ala the {j,n} style, then the transform\_indices
 *         matrix must be either a two-column matrix with columns representing
 *         the j and n indices, respectively, or a two-row matrix with rows
 *         representing j and n indices, respectively.
 * @param  intrp_ptr  Pointer for implementation of interrupt checking.
 * @param  flat  Pointer to a universal matrix of type MUTIL\_SINT32.
 *               The memory for this single-dimensional matrix is
 *               automatically allocated within the function and
 *               (upon return) will contain the flattened transform indices.
 * @param  level Pointer to a universal matrix of type MUTIL\_SINT32.
 *               The memory for this single-dimensional matrix is
 *               automatically allocated within the function and
 *               (upon return) will contain the decomposition levels
 *               coresponding to specified transform indices.
 * @param  osc   Pointer to a universal matrix of type MUTIL\_SINT32.
 *               The memory for this single-dimensional matrix is
 *               automatically allocated within the function and
 *               (upon return) will contain the oscillation indices
 *               coresponding to specified transform indices.
 * @same \begin{itemize}
 * \item #MUTIL_LIBEXPORT mutil_errcode wavs32_transform_packet_convert_indices(const sint32_mat *transform_indices, void *intrp_ptr);#
 * \end{itemize}
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_maximum_overlap_packet
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_packet_convert_indices(
  const univ_mat *transform_indices,
  void           *intrp_ptr,
  univ_mat       *flat,
  univ_mat       *level,
  univ_mat       *osc );

/* Function documented with universal matrix version above */
MUTIL_LIBEXPORT mutil_errcode wavs32_transform_packet_convert_indices(
  const sint32_mat *transform_indices,
  void             *intrp_ptr,
  sint32_mat       *flat,
  sint32_mat       *level,
  sint32_mat       *osc );

/** Extracts a discrete wavelet packet transform subset.
 * Extracts a subset of a wavelet packet transform matrix set
 * produced either by \Ref{wavuniv_transform_packet}
 * or \Ref{wavuniv_transform_maximum_overlap_packet}.
 * The subset is itself a matrix set and represents a projection
 * of the original time series onto one of many possible bases
 * in the wavelet packet transform.
 *
 * An extra DWPT atoms storage system is also returned and is explained as follows:
 *
 * If, for any given parent crystal in a DWPT, the number of
 * coefficients is odd, then then some sort of boundary treatment
 * must be applied to extend the crystal to an even length.
 * As this is somehwat arbitrary, another techniques, which we employ
 * in our DWPT and DWT functions, is to simply preserve the last
 * coefficient of an odd length crystal in a specal 'extra'
 * storage vector and put it back in
 * place upon synthesis. This method has the advantage of preserving
 * the energy of the original crystal with a disadvantage of having
 * to do a little bookkeeping.
 *
 * The extra coefficients and a map to those extra coefficients is
 * facilitated by the wav\_dwpt\_extra structure whose members are
 * an extra 'atoms' vector, a 'levelmap' vector, and an 'nelem'
 * member sued to keep count of the number of extra stored atoms.
 * For a J level DWPT, the levelmap will be
 * J element vector of integers such that
 * the jth element corresponds to the jth level of a DWPT.
 * If the jth element contains a negative integer it denotes
 * that the jth DWPT level does not contain any extra coefficients.
 * Conversely, if a positive integer p appears in the jth element,
 * it means that the jth DWPT level has extra coefficients (if any
 * crystal in a given level does, then they all do in that level),
 * and p denotes the index of the levelmap vector where the extra coefficients
 * for the first crystal (with an oscillation index of zero) is stored.
 * Thus, the extra coefficient for crystal W(j,n) can be found in
 * the (p + n)th element of the atoms vector. The nelem member of the structure
 * contains the total number of extra atoms, which can be zero. Since there
 * are mutils functions whihc complain about zero lenght matrices, this additional
 * member is needed.
 *
 * NOTE 1: This extraction is only needed for synthesis operations and since the
 * DWPT synthesis is based on a subset of the DWPT, this function is called from the
 * the \Ref{wavuniv_transform_packet_basis} function and both the DWPT and
 * extra atom structure is sent to the \Ref{wavuniv_transform_packet_inverse}
 * synthesis function.
 *
 * NOTE 2: By design, there are never any extra coefficients at the last level
 * of DWPT, regardless of whether it is a full decomposition or not. This is
 * why the levelmap vector is J elements long as opposed to the J+1 levels that
 * are stored in a DWPT (since the original time series at elvel 0 is also included).
 * Since this function will typically be used prior to synthesis operations,
 * an 'extra DWPT atoms' structure of type \Ref{_wav_dwpt_extra} is returned
 * as well. This information is needed to recosntruct an odd length parent crystal.
 * See the references for more details.
 *
 * @limits Only relevant for the DWPT and MODWPT.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwtc.h
 * @source wav\_dwtc.c
 * @library wavelets
 * @usage #err = wavuniv_transform_packet_basis( &dwpt, &transform_indices, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param  dwpt Pointer to a matrix set containing a discrete wavelet
 *         packet transform such as that produced by \Ref{wavuniv_transform_packet}.
 * @param  transform_indices  Pointer to a pre-allocated universal matrix
 *         of type MUTIL\_SINT32 containing the representative indices of the
 *         crystals which define the subset of a wavelet packet transform to analyze.
 *         The wavelet packet crystals are assumed to be arranged in the so-called
 *         natural order ala [W\_{0,0} , W\_{1,0}, W\_{1,1}, W\_{2,0}, W\_{2,1},
 *         W\_{2,2}, W\_{2,3}, ... , W\_{J,0}, ...  W\_{J, N\_J}]
 *         where J  is the number of decomposition levels and $N\_J = 2^J$.
 *         By definition, W\_{0,0} is the original time series. A given crystal
 *         is identified in the W\_{j,n} form above or by a flattened index.
 *         For example, the DWPT crystal in level 2 at oscillation index 1 is
 *         identified either by {j,n} = {2,1} or by its flattened index 4
 *         (with zero based indexing, 4 represents the fifth crystal of the
 *         wavelet packet transform in natural order). If the crystals
 *         are identified using flattened indices, then the transform\_indices
 *         matrix must be a single-column or single-row matrix. Conversely, if the
 *         crystals are identified ala the {j,n} style, then the transform\_indices
 *         matrix must be either a two-column matrix with columns representing
 *         the j and n indices, respectively, or a two-row matrix with rows
 *         representing j and n indices, respectively.
 * @param  intrp_ptr  Pointer for implementation of interrupt checking.
 * @param  result  Pointer to a matrix set containing universal matrices
 *               of type MUTIL\_DOUBLE. The size of each matrix
 *               will correspond to to the crystal of the original
 *               wavelet packet transform from which it was replicated.
 *               The memory for this matrix set (and its constitutive
 *               universal matrices) is automatically allocated within
 *               the function.
 * @param  extra Pointer to a structure of type \Ref{_wav_dwpt_extra}.
 * @see wavuniv_transform_packet_convert_indices
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_maximum_overlap_packet
 * @see _wav_dwpt_extra
 * @see wavuniv_transform_packet_inverse
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_packet_basis(
  const mat_set     *dwpt,
  const univ_mat    *transform_indices,
  void              *intrp_ptr,
  mat_set           *result,
  wav_dwpt_extra    *extra );

/** Discrete wavelet packet transform subset inverse.
 * Inveets a subset of discrete wavelet packet transform (DWPT)
 * crystals representing the projection of the original
 * time series onto one (of the many possible) discrete wavelet
 * transforms available in the DWPT.
 *
 * @limits Only relevant for the DWPT.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwtc.h
 * @source wav\_dwtc.c
 * @library wavelets
 * @usage #err = wavuniv_transform_packet_invert( &dwpt_basis, &transform_indices, &filters, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param  transform_indices  Pointer to a pre-allocated universal matrix
 *         of type MUTIL\_SINT32 containing the representative indices of the
 *         crystals which define the subset of a wavelet packet transform to analyze.
 *         The wavelet packet crystals are assumed to be arranged in the so-called
 *         natural order ala [W\_{0,0} , W\_{1,0}, W\_{1,1}, W\_{2,0}, W\_{2,1},
 *         W\_{2,2}, W\_{2,3}, ... , W\_{J,0}, ...  W\_{J, N\_J}]
 *         where J  is the number of decomposition levels and $N\_J = 2^J$.
 *         By definition, W\_{0,0} is the original time series. A given crystal
 *         is identified in the W\_{j,n} form above or by a flattened index.
 *         For example, the DWPT crystal in level 2 at oscillation index 1 is
 *         identified either by {j,n} = {2,1} or by its flattened index 4
 *         (with zero based indexing, 4 represents the fifth crystal of the
 *         wavelet packet transform in natural order). If the crystals
 *         are identified using flattened indices, then the transform\_indices
 *         matrix must be a single-column or single-row matrix. Conversely, if the
 *         crystals are identified ala the {j,n} style, then the transform\_indices
 *         matrix must be either a two-column matrix with columns representing
 *         the j and n indices, respectively, or a two-row matrix with rows
 *         representing j and n indices, respectively.
 * @param  intrp_ptr  Pointer for implementation of interrupt checking.
 * @param  result  Pointer to a matrix set containing universal matrices
 *               of type MUTIL\_DOUBLE. The size of each matrix
 *               will correspond to to the crystal of the original
 *               wavelet packet transform from which it was replicated.
 *               The memory for this matrix set (and its constitutive
 *               universal matrices) is automatically allocated within
 *               the function.
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_packet_convert_indices
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_packet_inverse(
  const mat_set  *dwpt_basis,
  const wav_dwpt_extra *extra,
  const univ_mat *transform_indices,
  const mat_set  *filters,
  void           *intrp_ptr,
  univ_mat       *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_DWTC_H_*/
