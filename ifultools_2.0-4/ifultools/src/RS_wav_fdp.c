
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_wav_fdp.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS wavelets library.

   Functions wrapped:

   wavuniv_fdp_estimator_instantaneous()
   wavuniv_fdp_estimator_block()
   wavuniv_fdp_bandpass_variance()
*/

#include "wav_fdp.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

#include "R_ext/Print.h"

/** Block-independent (instantaneous) estimation of fractionally differenced (FD) model parameters.
 * @source RS\_wav\_fdp.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_fdp_estimator_instantaneous", time.series, levels, filter.type, filter.length, estimator, biased, dof.order, delta.range))#
 * @return               An R ... containing ...
 * @param time.series    Pointer to an R object containing ... time.series
 * @param levels         Pointer to an R object containing ... levels
 * @param filter.type    Pointer to an R object containing ... filter.type
 * @param filter.length  Pointer to an R object containing ... filter.length
 * @param estimator      Pointer to an R object containing ... estimator
 * @param biased         Pointer to an R object containing ... biased
 * @param dof.order      Pointer to an R object containing ... dof.order
 * @param delta.range    Pointer to an R object containing ... delta.range
 * @see _wav_fdp_estimator
 * @see _wav_filter_type
 * @see wavuniv_fdp_estimator_block
 * @see wavuniv_fdp_bandpass_variance
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_filters_daubechies
*/
EXTERN_R SEXP RS_wavelets_fdp_estimator_instantaneous(
 SEXP pr_time_series,
 SEXP pr_levels,
 SEXP pr_filter_type,
 SEXP pr_filter_length,
 SEXP pr_estimator,
 SEXP pr_biased,
 SEXP pr_dof_order,
 SEXP pr_delta_range )
{
  SEXP              pr_ret_delta;                
  SEXP              pr_ret_innovation_variance;  
  SEXP              pr_ret_obj;                  
  SEXP              pr_ret_variance_delta;       
  boolean           biased;                      
  mutil_data_type   type;                        
  mutil_errcode     err;                         
  sint32            dof_order;                   
  sint32            filter_length;               
  univ_mat          delta;                       
  univ_mat          delta_range;                 
  univ_mat          innovation_variance;         
  univ_mat          levels;                      
  univ_mat          time_series;                 
  univ_mat          variance_delta;              
  void              *VPNULL = NULL;              
  wav_fdp_estimator estimator;                   
  wav_filter_type   filter_type;                 
  memlist           list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_levels to levels */
  READ_MATRIX_REGISTER( pr_levels, &levels );

  /* ... pr_filter_type to filter_type */
  WAV_FILTER_TYPE_FROM_R( pr_filter_type, &filter_type );

  /* ... pr_filter_length to filter_length */
  SINT32_FROM_R( pr_filter_length, &filter_length );

  /* ... pr_estimator to estimator */
  WAV_FDP_ESTIMATOR_FROM_R( pr_estimator, &estimator );

  /* ... pr_biased to biased */
  BOOLEAN_FROM_R( pr_biased, &biased );

  /* ... pr_dof_order to dof_order */
  SINT32_FROM_R( pr_dof_order, &dof_order );

  /* ... pr_delta_range to delta_range */
  READ_MATRIX_REGISTER( pr_delta_range, &delta_range );

  /* Call the function */
  err = wavuniv_fdp_estimator_instantaneous(
    &time_series,
    &levels,
    filter_type,
    filter_length,
    estimator,
    biased,
    dof_order,
    &delta_range,
    VPNULL,
    &delta,
    &variance_delta,
    &innovation_variance );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Problem calling wavuniv_fdp_estimator_instantaneous() function" );
  err = memlist_member_register( &list, &delta, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
  err = memlist_member_register( &list, &variance_delta, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
  err = memlist_member_register( &list, &innovation_variance, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );

  /* create the output R object */
  err = matuniv_to_R( &delta, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_delta );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );
  err = matuniv_to_R( &variance_delta, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_variance_delta );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );
  err = matuniv_to_R( &innovation_variance, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_innovation_variance );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  PROTECT( pr_ret_obj = allocVector( VECSXP, 3 ) );
  SET_VECTOR_ELT( pr_ret_obj, 0, pr_ret_delta );
  SET_VECTOR_ELT( pr_ret_obj, 1, pr_ret_variance_delta );
  SET_VECTOR_ELT( pr_ret_obj, 2, pr_ret_innovation_variance );
  UNPROTECT(1);

  /* free registered local memory */
  MUTIL_FREE_WARN( memlist, &list );

  return pr_ret_obj;
}

/** Block-dependent estimation of fractionally differenced (FD) model parameters.
 * @source RS\_wav\_fdp.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_fdp_estimator_block", time.series, levels, filter.type, filter.length, estimator, boundary.mode, edof.mode, sdf, delta.range))#
 * @return               An R ... containing ...
 * @param time.series    Pointer to an R object containing ... time.series
 * @param levels         Pointer to an R object containing ... levels
 * @param filter.type    Pointer to an R object containing ... filter.type
 * @param filter.length  Pointer to an R object containing ... filter.length
 * @param estimator      Pointer to an R object containing ... estimator
 * @param boundary.mode  Pointer to an R object containing ... boundary.mode
 * @param edof.mode      Pointer to an R object containing ... edof.mode
 * @param sdf            Pointer to an R object containing ... sdf
 * @param delta.range    Pointer to an R object containing ... delta.range
 * @see _wav_fdp_estimator
 * @see _wav_filter_type
 * @see wavuniv_fdp_estimator_instantaneous
 * @see wavuniv_variance_edof
 * @see wavuniv_fdp_bandpass_variance
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_filters_daubechies
 * @see wavuniv_variance
*/
EXTERN_R SEXP RS_wavelets_fdp_estimator_block(
 SEXP pr_time_series,
 SEXP pr_levels,
 SEXP pr_filter_type,
 SEXP pr_filter_length,
 SEXP pr_estimator,
 SEXP pr_boundary_mode,
 SEXP pr_edof_mode,
 SEXP pr_sdf,
 SEXP pr_delta_range )
{
  SEXP              pr_ret_delta;                
  SEXP              pr_ret_innovation_variance;  
  SEXP              pr_ret_obj;                  
  SEXP              pr_ret_variance_delta;       
  boolean           boundary_mode;               
  double            delta;                       
  double            innovation_variance;         
  double            variance_delta;              
  mutil_data_type   type;                        
  mutil_errcode     err;                         
  sint32            edof_mode;                   
  sint32            filter_length;               
  univ_mat          delta_range;                 
  univ_mat          levels;                      
  univ_mat          sdf;                         
  univ_mat          time_series;                 
  void              *VPNULL = NULL;              
  wav_fdp_estimator estimator;                   
  wav_filter_type   filter_type;                 
  memlist           list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_levels to levels */
  READ_MATRIX_REGISTER( pr_levels, &levels );

  /* ... pr_filter_type to filter_type */
  WAV_FILTER_TYPE_FROM_R( pr_filter_type, &filter_type );

  /* ... pr_filter_length to filter_length */
  SINT32_FROM_R( pr_filter_length, &filter_length );

  /* ... pr_estimator to estimator */
  WAV_FDP_ESTIMATOR_FROM_R( pr_estimator, &estimator );

  /* ... pr_boundary_mode to boundary_mode */
  BOOLEAN_FROM_R( pr_boundary_mode, &boundary_mode );

  /* ... pr_edof_mode to edof_mode */
  SINT32_FROM_R( pr_edof_mode, &edof_mode );

  /* ... pr_sdf to sdf */
  READ_MATRIX_REGISTER( pr_sdf, &sdf );

  /* ... pr_delta_range to delta_range */
  READ_MATRIX_REGISTER( pr_delta_range, &delta_range );

  /* Call the function */
  err = wavuniv_fdp_estimator_block(
    &time_series,
    &levels,
    filter_type,
    filter_length,
    estimator,
    boundary_mode,
    edof_mode,
    &sdf,
    &delta_range,
    VPNULL,
    &delta,
    &variance_delta,
    &innovation_variance );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Problem calling wavuniv_fdp_estimator_block() function" );

  /* create the output R object */
  err = double_to_R( delta, &pr_ret_delta );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );
  err = double_to_R( variance_delta, &pr_ret_variance_delta );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );
  err = double_to_R( innovation_variance, &pr_ret_innovation_variance );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  PROTECT( pr_ret_obj = allocVector( VECSXP, 3 ) );
  SET_VECTOR_ELT( pr_ret_obj, 0, pr_ret_delta );
  SET_VECTOR_ELT( pr_ret_obj, 1, pr_ret_variance_delta );
  SET_VECTOR_ELT( pr_ret_obj, 2, pr_ret_innovation_variance );
  UNPROTECT(1);

  /* free registered local memory */
  MUTIL_FREE_WARN( memlist, &list );

  return pr_ret_obj;
}

/** Mid-octave spectral density function (SDF) estimation.
 * @source RS\_wav\_fdp.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_fdp_bandpass_variance", levels, delta, n.sample, filter.length, estimator, boundary.mode, edof.mode, sdf, delta.range))#
 * @return          An R ... containing ...
 * @param levels    Pointer to an R object containing ... levels
 * @param delta     Pointer to an R object containing ... delta
 * @param n.sample  Pointer to an R object containing ... n.sample
 * @see wavuniv_fdp_estimator_instantaneous
 * @see wavuniv_variance_confidence
 * @see wavuniv_variance_edof
*/
EXTERN_R SEXP RS_wavelets_fdp_bandpass_variance(
 SEXP pr_levels,
 SEXP pr_delta,
 SEXP pr_n_sample )
{
  SEXP             pr_ret_result;   
  double           delta;           
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           n_sample;        
  univ_mat         levels;          
  univ_mat         result;          
  void             *VPNULL = NULL;  

  /* Avoid lint warning */
  (void) whatssi;

  /* Conversion of input data ... */

  /* ... pr_levels to levels */
  err = mutil_R_type( pr_levels, &type );
  if ( err ){
      PROBLEM "Unable to read pr_levels type" ERROR;
  }

  err = matuniv_from_R( pr_levels, type, &levels );
  if ( err ){
      PROBLEM "Unable to read pr_levels" ERROR;
  }

  /* ... pr_delta to delta */
  err = double_from_R( pr_delta, &delta );
  if ( err ){
    PROBLEM "Unable to convert double type argument pr_delta to delta" ERROR;
  }

  /* ... pr_n_sample to n_sample */
  err = sint32_from_R( pr_n_sample, &n_sample );
  if ( err ){
    PROBLEM "Unable to convert sint32 type argument pr_n_sample to n_sample" ERROR;
  }

  /* Call the function */

  err = wavuniv_fdp_bandpass_variance(
    &levels,
    delta,
    n_sample,
    VPNULL,
    &result );
  if ( err ){
    PROBLEM "Problem calling wavuniv_fdp_bandpass_variance() function" ERROR;
  }

  /* create the output R object */

  err = matuniv_to_R( &result, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_result );
  MUTIL_FREE_WARN( matuniv, &result );
  if ( err ) {
      PROBLEM "Unable to convert output data to R format" ERROR;
  }

  return pr_ret_result;
}

