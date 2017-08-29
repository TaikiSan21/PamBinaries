
/* $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_mac.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file */

#ifndef IN_WAV_MACRO_H_
#define IN_WAV_MACRO_H_


/* This file contains wavelet-related macro definitions
 */

/* Headers to define types used in this file */
#include "ut_type.h"
#include "wav_type.h"

#ifdef __cplusplus
extern "C" {
#endif


/** Extract a desired two-dimensional discrete wavelet transform subband.
 * Given a matrix set of the kind returned by
 * \Ref{wavuniv_transform_discrete_wavelet_2d}, a desired subband, and
 * a desired decomposition level, it returns the pointer to
 * the universal matrix containing the corresponding wavelet coefficients.
 *
 * This macro does not check that the requested decomposition level is
 * acceptable. If the requested subband is MUTIL\_SUBBAND\_SMOOTH, this
 * macro returns the last matrix in the matrix set, and ignores LEVEL.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_mac.h
 * @source wav\_mac.h
 * @library wavelet
 * @usage #univ_mat_ptr = WAVDWT_EXTRACT_SUBBAND_2D(&matset, MUTIL_SUBBAND_DIAGONAL, 3);#
 * @return Pointer to universal matrix corresponding to desired subband
 *   and decomposition level.
 * @param MATSET_PTR Pointer to matrix set containing the result of a
 *   two-dimensional wavelet decomposition.
 * @param SUBBAND    Enum for the type of subband desired. Must be of type
 *   wav\_subband\_dwt\_2d.
 * @param LEVEL Desired decomposition level.
 * @see _wav_subband_dwt_2d
 * @see wavuniv_transform_discrete_wavelet_2d
 */
#define WAVDWT_EXTRACT_SUBBAND_2D( MATSET_PTR, SUBBAND, LEVEL ) \
  ( (SUBBAND) == MUTIL_SUBBAND_SMOOTH ) ? \
    &((MATSET_PTR)->mats[ (MATSET_PTR)->nelem - 1 ]) : \
    &((MATSET_PTR)->mats[ 3 * ( (LEVEL) - 1 ) + (sint32) (SUBBAND) ] )


/** Index for the probability mass values, means, and variances of an
 * HMT model structure.  Compose an index for the flat data array of
 * the matrix from a set of dimensional indices.  This macro expands
 * to a C statement that computes an index to pick out a single data
 * value from the matrix.  No checking is done to determine the
 * validity of the indices.  Zero corresponds to the first element for
 * each index.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_mac.h
 * @source wav\_mac.h
 * @library wavelet
 * @usage #index = WAVHMT_MODEL_INDEX_PMF_MEAN_VAR(&model, state, sc_lev, subband ); temp_pmf = *( model.state_pmfs.data + index);#
 * @param MODEL    Pointer to HMT model.
 * @param STATE    State number of the data.
 * @param SCALE    Scale level of the data.
 * @param SUBBAND  Data subband; 0, 1, and 2 correspond respectively to
 *    the HL, LH, and HH subbands.
 * @see _wav_hmt_model
 * @see WAVHMT_MODEL_INDEX_TRANS
 */
#define WAVHMT_MODEL_INDEX_PMF_MEAN_VAR( MODEL, STATE, SCALE, SUBBAND ) \
  ( 3 * ( STATE * (MODEL)->num_scales + SCALE ) + SUBBAND )


/** Index for the transition probabilities of an HMT model structure.
 * Compose an index for the flat data array of the matrix, from a set
 * of dimensional indices.  This macro expands to a C statement that
 * computes an index to pick out a single data value from the matrix.
 * No checking is done to determine the validity of the indices.  Zero
 * corresponds to the first element for each index.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_mac.h
 * @source wav\_mac.h
 * @library wavelet
 * @usage #index = WAVHMT_MODEL_INDEX_TRANS(&model, ch_state, par_state, sc_lev, subband); temp_tprob = *(model.trans_probs.data + index);#
 * @param MODEL    Pointer to HMT model.
 * @param CSTATE   Child state number of the desired data.
 * @param PSTATE   Parent State number of the desired data.
 * @param SCALE    Scale level of the desired data.
 * @param SUBBAND  Data subband; 0, 1, and 2 correspond respectively to
 *   the HL, LH, and HH subbands.
 * @see _wav_hmt_model
 * @see WAVHMT_MODEL_INDEX_PMF_MEAN_VAR
 */
#define WAVHMT_MODEL_INDEX_TRANS( MODEL, CSTATE, PSTATE, SCALE, SUBBAND ) \
  ( 3 * ( ((MODEL)->num_scales - 1 ) * ( CSTATE * (MODEL)->num_states \
   + PSTATE ) + SCALE ) + SUBBAND )


/** Initialize the members of an HMT model structure to zero or NULL.
 * This macro is needed to facilitate re-training of a previously trained
 * HMT model, i.e., using previous model parameters to initialize the EM
 * algorithm. This macro sets all members of a \Ref{_wav_hmt_model}
 * structure to zero or NULL, including those of its \Ref{_double_mat}
 * members and all enumerated types.
 *
 * The function \Ref{wavuniv_hidden_markov_tree_2d} checks members of
 * the \Ref{_wav_hmt_model} structure to determine whether re-training
 * is to occur, then chooses the appropriate initialization scheme for
 * the EM algorithm.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_mac.h
 * @source wav\_mac.h
 * @library wavelet
 * @usage #WAVHMT_ZERO_MODEL( &model );#
 * @param MODEL    Pointer to HMT model.
 * @see _wav_hmt_model
 */
#define WAVHMT_ZERO_MODEL( MODEL ) \
  ( MODEL )->num_train = 0; \
  ( MODEL )->num_states = 0; \
  ( MODEL )->num_scales = 0; \
  ( MODEL )->boundary = ( mutil_boundary_type ) 0; \
  ( MODEL )->state_pmfs.nelem = 0; \
  ( MODEL )->state_pmfs.nrow = 0; \
  ( MODEL )->state_pmfs.ncol = 0; \
  ( MODEL )->state_pmfs.data = ( double* ) NULL; \
  ( MODEL )->means.nelem = 0; \
  ( MODEL )->means.nrow = 0; \
  ( MODEL )->means.ncol = 0; \
  ( MODEL )->means.data = ( double* ) NULL; \
  ( MODEL )->vars.nelem = 0; \
  ( MODEL )->vars.nrow = 0; \
  ( MODEL )->vars.ncol = 0; \
  ( MODEL )->vars.data = ( double* ) NULL; \
  ( MODEL )->trans_probs.nelem = 0; \
  ( MODEL )->trans_probs.nrow = 0; \
  ( MODEL )->trans_probs.ncol = 0; \
  ( MODEL )->trans_probs.data = ( double* ) NULL; \
  ( MODEL )->final_like = 0.0; \
  ( MODEL )->iters =0; \
  ( MODEL )->filt_type = ( wav_filter_type ) 0; \
  ( MODEL )->filt_taps = 0; \
  ( MODEL )->dwt_row_dims.nelem = 0; \
  ( MODEL )->dwt_row_dims.nrow = 0; \
  ( MODEL )->dwt_row_dims.ncol = 0; \
  ( MODEL )->dwt_row_dims.data = ( sint32* ) NULL; \
  ( MODEL )->dwt_col_dims.nelem = 0; \
  ( MODEL )->dwt_col_dims.nrow = 0; \
  ( MODEL )->dwt_col_dims.ncol = 0; \
  ( MODEL )->dwt_col_dims.data = ( sint32* ) NULL; \
  ( MODEL )->scl_mean = 0.0; \
  ( MODEL )->scl_var = 0.0

#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_MACRO_H_ */



