
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/wav_hmt.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file */

#ifndef IN_WAV_HMT_H_
#define IN_WAV_HMT_H_

#include "mat_type.h"
#include "wav_type.h"

/* This file contains function declarations for the wavelet
  Hidden Markov Tree modeling of images.
 */


#ifdef __cplusplus
extern "C" {
#endif


/** Hidden Markov tree image modeling.
 * Function to compute a wavelet transform-based hidden Markov tree
 * (HMT) model from a training set of images.
 *
 * Given a matrix set of training images, the function calls
 * \Ref{wavuniv_transform_discrete_wavelet_2d} to compute the 2-D
 * discrete wavelet transform (DWT) of each image using the
 * specified boundary condition.  The tree-like structure of the
 * wavelet coefficients is exploited to compute an HMT model via
 * the iterative estimation-maximization (EM) algorithm.
 *
 * The mean and variance parameters for the EM algorithm are
 * initialized using the sample means and variances in
 * each subband.  The reciprocal of num\_states is used
 * to initialize all probability mass function (PMF)
 * and transition probabilities.
 * During each iteration of the EM alogirthm
 * the total likelihood of the image set,
 * the multiplication of the individual
 * image likelihoods, is computed, and the
 * percentage difference from the previous iteration is calculated.
 * Training will stop if either the maximum number of iterations is
 * reached or the percentage difference in likelihood is smaller than
 * the specified tolerance.
 *
 * The resulting HMT model consists of means, variances and
 * PMF values corresponding to each hidden
 * state, scale and subband, along with transition probabilities to
 * the next to last scale.  {\em Tying} is used for all wavelet
 * coefficients within a particular scale/subband, i.e., within a
 * particular scale/subband, all wavelet coefficients share a single
 * state variable.  Specifically, for an M-state model trained from
 * 2-D DWTs with L scale levels, the means, variances, and PMFs of the
 * model will each be 1 by M x (L x 3) matrices that may be indexed using
 * the \Ref{WAVHMT_MODEL_INDEX_PMF_MEAN_VAR} macro.  The
 * \Ref{WAVHMT_MODEL_INDEX_TRANS} macro may be used to index the 1 by M x
 * M x (L-1) x 3 matrix of transition probabilities.
 *
 * The input hmt\_model may point to a \Ref{_wav_hmt_model} structure that has
 * been intialized via the \Ref{WAVHMT_ZERO_MODEL} macro, or a valid HMT
 * model structure obtained from a previous call to the function. In the
 * latter case, the initial PMFs, means, variances and transition
 * probabilities are taken from the model; in addition, values passed in
 * for num\_scales, num\_states, filter, taps, and boundary are ignored and
 * determined by members of the model structure.
 *
 * For more information see:
 * M. Crouse, R. Nowak and R. Baraniuk, ``Wavelet-Based Statistical
 * Signal Processing Using Hidden Markov Models,''
 * {\em IEEE Transactions on Signal Processing,} v. 46, n. 4, 1998.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_hmt.h
 * @source wav\_hmt.c
 * @library wavelets
 * @usage #err =  wavuniv_hidden_markov_tree_2d( &images, scales, states, iters, tol, WAV_FILTER_LEAST_ASYMMETRIC, 2, MUTIL_BOUNDARY_CONTINUE, TRUE, interrupt, &model );#
 * @return Standard mutils error/OK code.
 *
 * @param imset         Pointer to matrix set of training images, all
 *    of identical dimensions and of any real type.
 *    In addition, the matrices may be non-square and non-dyadic with row
 *    and columns sizes greater than 32 and less than 4096.
 * @param num_scales    Number of scales used for the 2-D DWT,
 *    must be 1 or more, but cannot exceed log2(min(N,M)), where N and M are
 *    the row and column size of the images in imset.
 * @param num_states    Number of possible states for each hidden state
 *    variable, must be greater than 1.
 * @param iterations    Maximum number of iterations for the training
 *    (EM) algorithm, must be greater than 0.
 * @param tolerance     Tolerance level used for algorithm stopping criteria.
 *     Taken as a percentage difference in likelihood between iterations,
 *     typically in the range (0,1).
 * @param filter        Type of wavelet filter used for 2-D DWT.
 * @param taps          Length of the wavelet filter.
 * @param boundary      Boundary extension method used in the 2-D DWT.
 * @param free_images   Flag to free image memory after 2-D DWT computation.
 *     If  TRUE then all memory for imset will be freed by the function.
 * @param intrp_ptr     Pointer to implement interrupt checking.
 * @param hmt_model     Pointer to unallocated model to
 *     hold the resulting model parameters.
 *
 * @see wavuniv_hmt_image_likelihoods
 * @see wavuniv_validate_hmt_model
 * @see _wav_hmt_model
 * @see WAVHMT_MODEL_INDEX_PMF_MEAN_VAR
 * @see WAVHMT_MODEL_INDEX_TRANS
 * @see WAVHMT_ZERO_MODEL
 * @see wavuniv_transform_discrete_wavelet_2d
 * @see _wav_filter_type
 * @see _mutil_boundary_type
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_hidden_markov_tree_2d(
  mat_set*              imset,
  sint32                num_scales,
  sint32                num_states,
  sint32                iters,
  double                tolerance,
  wav_filter_type       filter,
  sint32                taps,
  mutil_boundary_type   boundary,
  const boolean         free_images,
  void*                 intrp_ptr,
  wav_hmt_model*        hmt_model );


/** Simulate a set of images from an HMT model.
 * This function simulates a set of images from a given HMT model
 * assumed to be trained using the function
 * \Ref{wavuniv_hidden_markov_tree_2d}.
 *
 * A Gaussian distribution with mean and variance equal to model->scl\_mean
 * and model->scl\_var model structure members, respectively, will be used to
 * generate a scaling coefficient matrix for each
 * desired image. The mixture density defined by the model parameters
 * will be used to generate wavelet coefficients at each scale level. The
 * dimensions of all wavelet and scaling coefficient matrices are determined
 * by the model->dwt\_row\_dims and model->dwt\_col\_dims.
 * Dimensions of the simulated images(s) will be determined by the
 * \Ref{wavuniv_transform_discrete_wavelet_inverse_2d} function.
 *
 * Input num\_sim determines the number of images to be simulated.
 * The input levs will be passed to
 * \Ref{wavuniv_transform_discrete_wavelet_inverse_2d} and
 * determines the level of detail (wavelet) coefficients used for each
 * image reconstruction. A value of zero forces all wavelet coefficints to
 * be used. See
 * \Ref{wavuniv_transform_discrete_wavelet_inverse_2d} for
 * further details. The boolean
 * input zero\_scale gives the user the option to produce sample images
 * reconstructed from simulated wavelet coefficients only. If TRUE the scaling
 * coefficients are all set to zero prior to reconstruction. Parameter imset
 * is a pointer to a \Ref{_mat_set} structure that will hold the simulated
 * images. All memory for imset is allocated by this function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_hmt.h
 * @source wav\_hmt.c
 * @library wavelets
 * @usage #err =  wavuniv_hmt_simulate_2d( &model, nsim, lev, zs, intrp, images );#
 * @return Standard mutils error/OK code.
 * @param model     Pointer to input \Ref{_wav_hmt_model} model structure.
 * @param num_sim   Number of images to simulate.
 * @param levs      Reconstruct images to this scale level. A value of zero
 *                  forces complete reconstruction using simulated wavelet
 *                  coefficients at all levels.
 * @param zero_scale Boolean flag, it TRUE then the scaling coefficients in
 *                   in the simulates 2D DWTs are all set to zero and the
 *                   sample image is reconstructed only from the wavelet
 *                   coefficients. If FALSE, the scaling coefficients
 * @param intrp_ptr Pointer to implement interrupt checking.
 * @param imset     Pointer to matrix set of simulated images.
 * @see wavuniv_hidden_markov_tree_2d
 * @see _wav_hmt_model
 * @see wavuniv_transform_discrete_wavelet_inverse_2d
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_hmt_simulate_2d(
  const wav_hmt_model*    model,
  const sint32            num_sim,
  const sint32            levs,
  const boolean           zero_scale,
  void*                   intrp_ptr,
  mat_set*                imset );


/** Compute image log-likelihoods given a hidden Markov tree (HMT)
 * model.
 *
 * Given a set of images and HMT model, the function computes the
 * log-likelihood of each image in the set. For each image, the 2-D
 * discrete wavelet transform (DWT) is computed and the means, variances,
 * probability mass functions values and transition probabilites of
 * the given HMT model are then used to perform one iteration of the
 * expectation-maximization (EM) algorithm.
 * The conditional likelihoods and joint density function values
 * produced by the EM algorithm are then used to compute the final
 * image likelihoods.
 *
 * If free\_images is TRUE then all memory for the imset
 * matrix set will be freed within the function.  The results are
 * returned in the 1 x N likelihoods matrix, where N is the
 * number of images in imset.  The data in likelihoods is
 * ordered in correspondence with the images in imset.
 *
 * For more information see:
 * M. Crouse, R. Nowak and R. Baraniuk, ``Wavelet-Based Statistical Signal
 * Processing Using Hidden Markov Models,'' {\em IEEE Trans. Sig. Proc.,}
 * v. 46, n. 4, 1998.
 *
 * @limits Only images of dyadic dimensions up
 * to 4096 x 4096 and type MUTIL\_DOUBLE are supported.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_hmt.h
 * @source wav\_hmt.c
 * @library wavelets
 * @usage #err = wavuniv_hmt_image_likelihoods( &imset, &model, TRUE, intrp_ptr, &lhoods );#
 * @return Standard mutils error/OK code.
 * @param imset         Pointer to \Ref{_mat_set} structure of training
 *                        images. Univeral matrices in the matrix set may
 *                      any real data type.
 * @param model         Pointer to \Ref{_wav_hmt_model} structure with
 *                        trained HMT model.
 * @param free_images   Flag to determine whether memory for the input
 *                        image set is freed.
 * @param intrp_ptr     Pointer to implement user interrupts.
 * @param likelihoods   Pointer to an unallocated \Ref{_double_mat} structure
 *                        to store
 *                      resulting likelihoods.
 * @see wavuniv_hidden_markov_tree_2d
 * @see _wav_hmt_model
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_hmt_image_likelihoods(
  mat_set*               imset,
  const wav_hmt_model*   model,
  const boolean          free_images,
  void*                  intrp_ptr,
  double_mat*            likelihoods );


/** Validate a given Hidden Markov Tree (HMT) model.
 * This function validates the data of a given \Ref{_wav_hmt_model} structure.
 * The structure is considered valid if the following conditions on its data
 * members are satisfied:
 * \begin{itemize}
 * \item num\_images is greater than zero
 * \item num\_scales is greater than one
 * \item num\_states is greater than one
 * \item boundary is a valid \Ref{_mutil_boundary_type}
 * \item filt\_type is a valid \Ref{_wav_filter_type}
 * \item filt\_taps is at least 2
 * \item the state\_pmfs, means, vars and trans\_probs
 * matrices are all valid according to \Ref{matuniv_validate}
 * \item the number of elements in these matrices are consistent with
 * respect to num\_states and num\_scales.
 * \end{itemize}
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_hmt.h
 * @source wav\_hmt.c
 * @library wavelets
 * @usage #err =  wavuniv_validate_hmt_model( &model );#
 * @return Standard mutils error/OK code.
 * @param model     Pointer to input \Ref{_wav_hmt_model} model structure.
 * @see wavuniv_hidden_markov_tree_2d
 * @see _wav_hmt_model
 * @see _wav_filter_type
 * @see _mutil_boundary_type
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_validate_hmt_model(
  const wav_hmt_model*   model );


/** Free allocated memory from a given Hidden Markov Tree (HMT) model structure.
 * This function frees memory previously allocated for the parameters of an
 * HMT model structure of type \Ref{_wav_hmt_model}.
 *
 * Memory for the means, variances, probability mass function values and
 * transition probabilities will be freed. The model structure is otherwise
 * unaltered.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_hmt.h
 * @source wav\_hmt.c
 * @library wavelets
 * @usage #err =  wavuniv_free_hmt_model( &model );#
 * @return Standard mutils error/OK code.
 * @param model     Pointer to model structure to be freed.
 * @see wavuniv_hidden_markov_tree_2d
 * @see _wav_hmt_model
 */
MUTIL_LIBEXPORT mutil_errcode wavuniv_hmt_model_free(
  wav_hmt_model*   model );

#ifdef __cplusplus
}
#endif

#endif /* IN_WAV_HMT_H_ */
