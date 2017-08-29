
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_wav_wtmm.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS wavelets library.

   Functions wrapped:

   wavuniv_transform_continuous_wavelet_modulus_maxima()
   wavuniv_transform_continuous_wavelet_modulus_maxima_tree()
*/

#include "wav_wtmm.h"
#include "wav_type.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** The modulus maxima of a continuous wavelet transform.
 * @source RRS\_wav\_wtmm.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_continuous_wavelet_modulus_maxima", cwt, tolerance))#
 * @return           An R ... containing ...
 * @param cwt        Pointer to an R object containing ... cwt
 * @param tolerance  Pointer to an R object containing ... tolerance
 * @see wavuniv_transform_continuous_wavelet
 * @see wavuniv_transform_continuous_wavelet_modulus_maxima_tree
*/
EXTERN_R SEXP RS_wavelets_transform_continuous_wavelet_modulus_maxima(
 SEXP pr_cwt,
 SEXP pr_tolerance,
 SEXP pr_peak_type )
{
  SEXP             pr_ret_iscale;
  SEXP             pr_ret_itime;
  SEXP             pr_ret_obj;
  mutil_data_type  type;
  mutil_errcode    err;
  univ_mat         iscale;
  univ_mat         itime;
  univ_mat         cwt;
  univ_mat         tolerance;
  void             *VPNULL = NULL;
  memlist          list;
  wav_transform_peak peak_type;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_cwt to cwt */
  READ_MATRIX_REGISTER( pr_cwt, &cwt );

  /* ... pr_tolerance to tolerance */
  READ_MATRIX_REGISTER( pr_tolerance, &tolerance );

  /* ... pr_peak_type to peak_type */
  WAV_TRANSFORM_PEAK_FROM_R( pr_peak_type, &peak_type );

  /* Call the function */
  err = wavuniv_transform_continuous_wavelet_modulus_maxima(
    &cwt,
    &tolerance,
    peak_type,
    VPNULL,
    &(itime.mat.s32mat),
    &(iscale.mat.s32mat) );
  MUTIL_FREE_WARN( matuniv, &cwt );
  MUTIL_FREE_WARN( matuniv, &tolerance );
  if ( err ){
    PROBLEM "Problem calling wavuniv_transform_continuous_wavelet_modulus_maxima() function" ERROR;
  }

  /* free registered local memory */
  MUTIL_FREE_WARN( memlist, &list );

  /* wrap the sint32 output matrices into a universal matrix */

  err = matuniv_wrap_matrix( &itime, &(itime.mat.s32mat), MUTIL_SINT32 );
  if (err){
    MUTIL_FREE_WARN( matuniv, &itime );
    PROBLEM "Problem wrapping WTMM time index matrix into a universal matrix" ERROR;
  }

  err = matuniv_wrap_matrix( &iscale, &(iscale.mat.s32mat), MUTIL_SINT32 );
  if (err){
    MUTIL_FREE_WARN( matuniv, &iscale );
    PROBLEM "Problem wrapping WTMM scale index matrix into a universal matrix" ERROR;
  }

  /* create the output R object */

  err = matuniv_to_R( &itime, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_itime );
  MUTIL_FREE_WARN( matuniv, &itime );
  if ( err ) {
      PROBLEM "Unable to convert output data to Splus format" ERROR;
  }

  err = matuniv_to_R( &iscale, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_iscale );
  MUTIL_FREE_WARN( matuniv, &iscale );
  if ( err ) {
    PROBLEM "Unable to convert output data to Splus format" ERROR;
  }

  PROTECT( pr_ret_obj = allocVector( VECSXP, 2 ) );
  SET_VECTOR_ELT( pr_ret_obj, 0, pr_ret_itime );
  SET_VECTOR_ELT( pr_ret_obj, 1, pr_ret_iscale );
  UNPROTECT(1);

  return pr_ret_obj;
}

/** The modulus maxima tree of a continuous wavelet transform.
 * @source RRS\_wav\_wtmm.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_continuous_wavelet_modulus_maxima_tree", wtmm.time.index, wtmm.scale.index, cwt, cwt.time, cwt.scale, bridge.gaps, n.octave.min, wtmm.strength.min))#
 * @return                   An R ... containing ...
 * @param wtmm.time.index    Pointer to an R object containing ... wtmm.time.index
 * @param wtmm.scale.index   Pointer to an R object containing ... wtmm.scale.index
 * @param cwt                Pointer to an R object containing ... cwt
 * @param cwt.time           Pointer to an R object containing ... cwt.time
 * @param cwt.scale          Pointer to an R object containing ... cwt.scale
 * @param bridge.gaps        Pointer to an R object containing ... bridge.gaps
 * @param n.octave.min       Pointer to an R object containing ... n.octave.min
 * @param wtmm.strength.min  Pointer to an R object containing ... wtmm.strength.min
 * @see wavuniv_transform_continuous_wavelet
 * @see wavuniv_transform_continuous_wavelet_modulus_maxima
*/
EXTERN_R SEXP RS_wavelets_transform_continuous_wavelet_modulus_maxima_tree(
 SEXP pr_wtmm_time_index,
 SEXP pr_wtmm_scale_index,
 SEXP pr_cwt,
 SEXP pr_cwt_time,
 SEXP pr_cwt_scale,
 SEXP pr_bridge_gaps,
 SEXP pr_n_octave_min,
 SEXP pr_wtmm_strength_min )
{
  boolean            bridge_gaps;
  double             n_octave_min;
  double             wtmm_strength_min;
  mat_set            result;
  memlist            list;
  mutil_data_type    type;
  mutil_errcode      err;
  SEXP               pr_ret_result;
  univ_mat           cwt;
  univ_mat           cwt_scale;
  univ_mat           cwt_time;
  univ_mat           wtmm_scale_index;
  univ_mat           wtmm_time_index;
  void              *VPNULL = NULL;

  /* Avoid lint warning */

  (void) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  READ_MATRIX_REGISTER( pr_wtmm_time_index, &wtmm_time_index );
  READ_MATRIX_REGISTER( pr_wtmm_scale_index, &wtmm_scale_index );
  READ_MATRIX_REGISTER( pr_cwt, &cwt );
  READ_MATRIX_REGISTER( pr_cwt_scale, &cwt_scale );
  READ_MATRIX_REGISTER( pr_cwt_time, &cwt_time );

  /* ... pr_bridge_gaps to bridge_gaps */
  BOOLEAN_FROM_R( pr_bridge_gaps, &bridge_gaps );

  /* ... pr_n_octave_min to n_octave_min */
  DOUBLE_FROM_R( pr_n_octave_min, &n_octave_min );

  /* ... pr_wtmm_strength_min to wtmm_strength_min */
  DOUBLE_FROM_R( pr_wtmm_strength_min, &wtmm_strength_min );

  /* Call the function */
  err = wavuniv_transform_continuous_wavelet_modulus_maxima_tree(
    &(wtmm_time_index.mat.s32mat),
    &(wtmm_scale_index.mat.s32mat),
    &(cwt.mat.cpxmat),
    &(cwt_time.mat.dblmat),
    &(cwt_scale.mat.dblmat),
    bridge_gaps,
    n_octave_min,
    wtmm_strength_min,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( wavuniv_transform_continuous_wavelet_modulus_maxima_tree, &result, &pr_ret_result );
}

