
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_modw.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_WAV_MODW_H_
#define IN_WAV_MODW_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "wav_type.h"

/* This file contains function declarations the wavelet mutils library. */

#ifdef __cplusplus
extern "C" {
#endif

/** The maximum overlap discrete wavelet transform (MODWT).
 * Given j, t are the decomposition level,
 * and time index, respectively, and V(0,t) = X(t) where X(t) is an N point
 * time series, the MODWT is given by
 * \begin{center}
 * \TEX{
 * \begin{eqnarray}
 * \widetilde{W}_{j,t} &\equiv& \sum_{l=0}^{L-1}\widetilde{h}_l
 * \widetilde{V}_{j-1, t - 2^{j-1}\;l \mbox{ mod }N}\nonumber\\
 * \widetilde{V}_{j,t} &\equiv& \sum_{l=0}^{L-1}\widetilde{g}_l
 * \widetilde{V}_{j-1, t - 2^{j-1}\;l \mbox{ mod }N}\nonumber
 * \end{eqnarray}
 * }
 * \end{center}
 *
 * for t = 0, ..., N - 1. The variable L is the length of both the scaling
 * filter (g) and wavelet filter (h).
 * W(j,t) are the wavelet coefficients at decomposition level j and time
 * index t while V(j,t) are the scaling coefficients.
 * The MODWT is a non-decimated form of the discrete wavelet transform (DWT)
 * having many advantages over the DWT including the ability
 * to handle arbitrary length sequences and shift invariance. The cost of
 * the MODWT is in its redundancy. For
 * an  N point input sequence, there are  N wavelet
 * coefficients per scale. However, the number of multiplication operations is
 * O(N log2(N)) which is the same as the fast Fourier transform,
 * and is acceptably fast
 * for most situations. Along with a matrix containing the MODWT
 * wavelet coefficients in [row, col] = [ time,scale ] format, the scaling
 * coefficient matrix is also returned in the same format.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_modw.h
 * @source wav\_modw.c
 * @library wavelets
 * @usage #err = wavuniv_transform_maximum_overlap(&time_series, &filters, n_level, intrp_ptr, &result);#
 * @return Standard mutils error/OK code.
 * @param  time_series Pointer to a pre-allocated universal matrix
 *                     containing the time series to analyze.
 *                     This input must a single-row or single-column
 *                     matrix and can be of any type with the
 *                     exception of MUTIL\_DCOMPLEX.
 * @param  filters     Pointer to a matrix set containing two pre-allocated
 *                     universal matrices of type MUTIL\_DOUBLE. The matrices
 *                     must contain (in order) valid wavelet and scaling
 *                     filters. The filter matrices must be a single-column
 *                     or single-row and each must contain the same number
 *                     of elements.
 * @param  n_level     The number of decomposition levels. This input should
 *                     not exceed $2^{\lfloor \log_2{N}\rfloor}$.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  result      Pointer to a matrix set to store the result.
 *                     The matrix set header and matrices are automatically
 *                     allocated by this function and (upon return) will
 *                     contain n\_level + 1 single-row universal matrices of
 *                     type MUTIL\_DOUBLE. The order of the matrices is
 *                     given by [W\_1 , W\_2 , ... , W\_J,  V\_J] where
 *                     J  = n\_level.
 * @see wavuniv_filters_daubechies
 * @see wavuniv_transform_maximum_overlap_inverse
 * @see wavuniv_transform_maximum_overlap_packet
 * @see wavuniv_transform_packet_detail
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_maximum_overlap(
   const univ_mat *time_series,
   const mat_set  *filters,
   sint32          n_level,
   void           *intrp_ptr,
   mat_set        *result );

/** The maximum overlap discrete wavelet packet transform (MODWPT).
 * Given j, n, t are the decomposition level, local node index,
 * and time index, respectively, the MODWPT is given by
 * \begin{center}
 * $\widetilde{W}_{j,n,t} \equiv \sum_{l=0}^{L-1}\widetilde{u}_{n,l}
 * \widetilde{W}_{j-1,\lfloor n/2 \rfloor, t - 2^{j-1}\;l \mbox{ mod }N}$
 * \end{center}
 *
 * L is the length of the filters defined by
 *
 * \begin{center}
 * \TEX{
 * \begin{eqnarray}
 * \widetilde{u}_{n,l} &=& g_l/\sqrt{2} \mbox{ if }n \mbox{ mod }4 = 0
 * \mbox{ or }3 \nonumber \\
 *                     &=& h_l/\sqrt{2} \mbox{ if }n \mbox{ mod }4 = 1
 * \mbox{ or }2 \mbox{ \qquad  for } l = 0,\ldots, L-1\nonumber
 * \end{eqnarray}
 * }
 * \end{center}
 *
 * where g and h are the scaling filter and wavelet filter, respectively.
 * By definition, $\widetilde{W}_{0,0,t} \equiv X_t$ where X is the original
 * time series.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_modw.h
 * @source wav\_modw.c
 * @library wavelets
 * @usage  #err = wavuniv_transform_maximum_overlap_packet(&time_series, &filters, n_level, intrp_ptr, &result);#
 * @return Standard mutils error/OK code.
 * @param  time_series Pointer to a pre-allocated universal matrix
 *                     containing the time series to analyze.
 *                     This input must a single-row or single-column
 *                     matrix and can be of any type with the
 *                     exception of MUTIL\_DCOMPLEX.
 * @param  filters     Pointer to a matrix set containing two pre-allocated
 *                     universal matrices of type MUTIL\_DOUBLE. The matrices
 *                     must contain (in order) valid wavelet and scaling
 *                     filters. The filter matrices must be a single-column
 *                     or single-row and each must contain the same number
 *                     of elements.
 * @param  n_level     The number of decomposition levels. This input should
 *                     not exceed $2^{\lfloor \log_2{N}\rfloor}$.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
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
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_transform_maximum_overlap_inverse
 * @see wavuniv_transform_packet_detail
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_maximum_overlap_packet(
   const univ_mat *time_series,
   const mat_set  *filters,
   sint32          n_level,
   void           *intrp_ptr,
   mat_set        *result );


/** Calculate the detail sequence for a specified crystal of a 1-D wavelet transform.
 * Given the location of a single crystal in a maximum wavelet (packet) transform, 
 * this function calculates the corresponding "detail" sequence. 
 * The detail sequence is formed by reconstructing the transform
 * using only the specified crystal and with all other coefficients (found in other
 * crystals) set to zero. A set of detail sequences whose collective frequency content
 * spans the normalized frequency range [0, 1/2] forms an additive decomposition that
 * can be summed to (re)form the original time series. There are many such
 * collections that can be formed using a subset of all crystals in a 
 * wavelet packet transform (MO)DWPT, including the (MO)DWT.
 *
 * The location of a particular crystal is given by the so-called "natural order" of a
 * wavelet packet tree (j,n), where where j is the decomposition level and n is
 * the local node index. By definition, the original time series is located at
 * crystal (0,0) in the wavelet packet tree. The first level (MO)DWPT contains
 * crystals (1,0) and (1,1), the second level contains crystals (2,0),(2,1),(2,2)
 * and (2,3), and so on. Each level of a wavelet packet transform
 * contains exactly $2^j$ crystals with the local node index range given by
 * $0 \le n \le 2^j - 1$. The wavelet transform coefficients at crystal (j,n)
 * are nominally associated with the normalized frequency band
 * $f \in [n/2^{j+1}, (n+1)/2^{j+1}]$. For example, crystal (2,3)
 * nominally corresponds to frequencies [3/8, 1/2]. Thus, an increase in the
 * local node index (n) corresponds to an increase in frequency content.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_modw.h
 * @source wav\_modw.c
 * @library wavelets
 * @usage  #err = wavuniv_transform_packet_detail( &transform, &filters, level, node, type, intrp_ptr, &result);#
 * @return Standard mutils error/OK code.
 * @param  transform Pointer to a pre-allocated matrix set
 *                     containing the MODWPT, MODWT. DWPT, or DWT. Each 
 *                     crystal of the transform must be stored as a single-row
 *                     universal matrix of type MUTIL\_DOUBLE.
 *                     The collection of crystals for the (MO)DWPT
 *                     must be ordered using the "natural" order such
 *                     that (for a level 2 (MO)DWPT) the first 4 crystals
 *                     would be (0,0),(1,0),(1,1),(2,0). The (MO)DWT
 *                     should be stored from high to low frequency content,
 *                     so that in a level 2 (MO)DWT, for example, the
 *                     crystals should be stored in the order (1,1),(2,1),(2,0).
 * @param  filters Pointer to a matrix set containing two pre-allocated
 *                     universal matrices of type MUTIL\_DOUBLE. The matrices
 *                     must contain (in order) valid wavelet and scaling
 *                     filters. The filter matrices must be a single-column
 *                     or single-row and each must contain the same number
 *                     of elements.
 * @param  level The decomposition level of the crystal to decompose.
 * @param  node The local node index of the crystal to decompose.
 * @param  type An enumerated type denoting the type of 
 *              wavelet transform used in developing the "transform" argument.
 *              Acceptable types are WAV\_TRANSFORM\_DWT, WAV\_TRANSFORM\_DWPT,
 *              WAV\_TRANSFORM\_MODWT, and WAV\_TRANSFORM\_MODWPT,
 *              See \Ref{_wav_transform} for details.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  result Pointer to a single-column or single-row universal
 *                     matrix of type MUTIL\_DOUBLE to store the result.
 *
 * @see wavuniv_transform_discrete_wavelet_convolution
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_transform_maximum_overlap_packet
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_packet_detail(
   const mat_set  *transform,
   const mat_set  *filters,
   sint32          level,
   sint32          node,
   wav_transform   xformtype,
   void           *intrp_ptr,
   univ_mat       *result );


/** The inverse maximum overlap discrete wavelet transform (IMODWT).
 * Given j, t are the decomposition level, and time index, respectively,
 * V(J,t) are the MODWT scaling coefficients at the largest scale (J),
 * and W(j,t) are the MODWT wavelet coefficients, the inverse MODWT is
 * given by recursively applying
 * \begin{center}
 * \[ \widetilde{V}_{j-1,t} \equiv \sum_{l=0}^{L-1}
 *  \widetilde{h}_l\widetilde{W}_{j, t + 2^{j-1}\;l \mbox{ mod }N}
 *  + \sum_{l=0}^{L-1}\widetilde{g}_l\widetilde{V}_{j, t +
 *  2^{j-1}\;l \mbox{ mod }N} \]
 * \end{center}
 *
 * for t = 0, ..., N - 1 and j = J, J-1, ... , 1. The variable N is the
 * length of the original time series and
 * L is the length of both the scaling filter (g) and wavelet filter (h).
 * V(j,t) are the reconstructed scaling coefficients.
 * This routine calculates the inverse MODWT using largest scale
 * scaling coefficients and the all of the wavelet coefficients
 * from the smallest to the largest scale.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_modw.h
 * @source wav\_modw.c
 * @library wavelets
 * @usage #err = wavuniv_transform_maximum_overlap_inverse(&modwt, &filters, intrp_ptr, &result);#
 * @return Standard mutils error/OK code.
 * @param  modwt       Pointer to a matrix set containing the MODWT
 *                     in n\_level + 1 pre-allocated single-row universal
 *                     matrices of type MUTIL\_DOUBLE. The order of the
 *                     matrices must be [W\_1 , W\_2 , ... , W\_J,  V\_J]
 *                     where J  = n\_level.
 * @param  filters     Pointer to a matrix set containing two pre-allocated
 *                     universal matrices of type MUTIL\_DOUBLE. The matrices
 *                     must contain (in order) valid wavelet and scaling
 *                     filters. The filter matrices must be a single-column
 *                     or single-row and each must contain the same number
 *                     of elements.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  result      Pointer to universal matrix to store the synthesis.
 *                     The memory for this matrix is automatically allocated
 *                     by this function and (upon return) will contain
 *                     a single-row universal matrix of type MUTIL\_DOUBLE
 *                     containing the synthesis.
 *
 * @see wavuniv_filters_daubechies
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_transform_maximum_overlap_packet
 * @see wavuniv_transform_packet_detail
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_transform_maximum_overlap_inverse(
   const mat_set  *modwt,
   const mat_set  *filters,
   void           *intrp_ptr,
   univ_mat       *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_MODW_H_*/
