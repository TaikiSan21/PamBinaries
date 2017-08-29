
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_wav_var.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS wavelets library.

   Functions wrapped:

   wavuniv_variance()
   wavuniv_variance_confidence()
   wavuniv_variance_edof()
*/

#include "wav_var.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Discrete wavelet variance estimation.
 * @source RS\_wav\_var.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_variance", time.series, transform.type, filter.type, filter.length, n.level, sdf))#
 * @return                An R ... containing ...
 * @param time.series     Pointer to an R object containing ... time.series
 * @param transform.type  Pointer to an R object containing ... transform.type
 * @param filter.type     Pointer to an R object containing ... filter.type
 * @param filter.length   Pointer to an R object containing ... filter.length
 * @param n.level         Pointer to an R object containing ... n.level
 * @param sdf             Pointer to an R object containing ... sdf
 * @see wavuniv_variance_confidence
 * @see wavuniv_variance_edof
*/
EXTERN_R SEXP RS_wavelets_variance(
 SEXP pr_time_series,
 SEXP pr_transform_type,
 SEXP pr_filter_type,
 SEXP pr_filter_length,
 SEXP pr_n_level,
 SEXP pr_sdf )
{
  SEXP              pr_ret_confidence;      
  SEXP              pr_ret_edof;            
  SEXP              pr_ret_obj;             
  SEXP              pr_ret_variance_block;  
  SEXP              pr_ret_variance_time;   
  mat_set           confidence;             
  mat_set           edof;                   
  mat_set           variance_block;         
  mat_set           variance_time;          
  mutil_data_type   type;                   
  mutil_errcode     err;                    
  sint32            filter_length;          
  sint32            n_level;                
  univ_mat          sdf;                    
  univ_mat          time_series;            
  void              *VPNULL = NULL;         
  wav_filter_type   filter_type;            
  wav_transform     transform_type;         
  memlist           list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_transform_type to transform_type */
  WAV_TRANSFORM_FROM_R( pr_transform_type, &transform_type );

  /* ... pr_filter_type to filter_type */
  WAV_FILTER_TYPE_FROM_R( pr_filter_type, &filter_type );

  /* ... pr_filter_length to filter_length */
  SINT32_FROM_R( pr_filter_length, &filter_length );

  /* ... pr_n_level to n_level */
  SINT32_FROM_R( pr_n_level, &n_level );

  /* ... pr_sdf to sdf */
  READ_MATRIX_REGISTER( pr_sdf, &sdf );

  /* Call the function confidence intervals are not available for DWT-based
  wavelet variance estimates */

  if ( transform_type == WAV_TRANSFORM_MODWT ){

    err = wavuniv_variance(
      &time_series,
      transform_type,
      filter_type,
      filter_length,
      n_level,
      &sdf,
      VPNULL,
      &variance_time,
      &variance_block,
      &confidence,
      &edof );
    MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Problem calling wavuniv_variance() function" );
    err = memlist_member_register( &list, &variance_time, MEMTYPE_MATSET);
    MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
    err = memlist_member_register( &list, &variance_block, MEMTYPE_MATSET);
    MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
    err = memlist_member_register( &list, &confidence, MEMTYPE_MATSET);
    MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
    err = memlist_member_register( &list, &edof, MEMTYPE_MATSET);
    MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );

    /* create the output R object */
    err = matset_to_R_list( &variance_time, &pr_ret_variance_time );
    MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );
    err = matset_to_R_list( &variance_block, &pr_ret_variance_block );
    MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );
    err = matset_to_R_list( &confidence, &pr_ret_confidence );
    MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );
    err = matset_to_R_list( &edof, &pr_ret_edof );
    MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

    PROTECT( pr_ret_obj = allocVector( VECSXP, 4 ) );
    SET_VECTOR_ELT( pr_ret_obj, 0, pr_ret_variance_time );
    SET_VECTOR_ELT( pr_ret_obj, 1, pr_ret_variance_block );
    SET_VECTOR_ELT( pr_ret_obj, 2, pr_ret_confidence );
    SET_VECTOR_ELT( pr_ret_obj, 3, pr_ret_edof );
    UNPROTECT(1);
  }
  else{
    err = wavuniv_variance(
      &time_series,
      transform_type,
      filter_type,
      filter_length,
      n_level,
      &sdf,
      VPNULL,
      &variance_time,
      &variance_block,
      (mat_set *) NULL,
      (mat_set *) NULL);
    MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Problem calling wavuniv_variance() function" );
    err = memlist_member_register( &list, &variance_time, MEMTYPE_MATSET);
    MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
    err = memlist_member_register( &list, &variance_block, MEMTYPE_MATSET);
    MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );

    /* create the output R object */
    err = matset_to_R_list( &variance_time, &pr_ret_variance_time );
    MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );
    err = matset_to_R_list( &variance_block, &pr_ret_variance_block );
    MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

    PROTECT( pr_ret_obj = allocVector( VECSXP, 2 ) );
    SET_VECTOR_ELT( pr_ret_obj, 0, pr_ret_variance_time );
    SET_VECTOR_ELT( pr_ret_obj, 1, pr_ret_variance_block );
    UNPROTECT(1);
  }

  /* free registered local memory */
  MUTIL_FREE_WARN( memlist, &list )
;
  
  return pr_ret_obj;
}

/** Confidence intervals for the unbiased and blocked averaged discrete wavelet variance estimates.
 * @source RS\_wav\_var.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_variance_confidence", variance, edof, probability, filter.length, n.level, sdf))#
 * @return             An R ... containing ...
 * @param variance     Pointer to an R object containing ... variance
 * @param edof         Pointer to an R object containing ... edof
 * @param probability  Pointer to an R object containing ... probability
 * @see wavuniv_variance
 * @see wavuniv_variance_edof
*/
EXTERN_R SEXP RS_wavelets_variance_confidence(
 SEXP pr_variance,
 SEXP pr_edof,
 SEXP pr_probability )
{
  SEXP             pr_ret_result;   
  double           probability;     
  mat_set          result;          
  mutil_data_type  type;            
  mutil_errcode    err;             
  univ_mat         edof;            
  univ_mat         variance;        
  void             *VPNULL = NULL;  
  memlist           list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_variance to variance */
  READ_MATRIX_REGISTER( pr_variance, &variance );

  /* ... pr_edof to edof */
  READ_MATRIX_REGISTER( pr_edof, &edof );

  /* ... pr_probability to probability */
  DOUBLE_FROM_R( pr_probability, &probability );

  /* Call the function */
  err = wavuniv_variance_confidence(
    &variance,
    &edof,
    probability,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( wavuniv_variance_confidence, &result, &pr_ret_result );
}

/** Equivalent degrees of freedom (EDOF) estimates for a chi-squared distribution assumption on the interior wavelet coefficients.
 * @source RS\_wav\_var.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_variance_edof", interior, num.coefs, variance, level, sdf, filter.type, filter.length))#
 * @return               An R ... containing ...
 * @param interior       Pointer to an R object containing ... interior
 * @param num.coefs      Pointer to an R object containing ... num.coefs
 * @param variance       Pointer to an R object containing ... variance
 * @param level          Pointer to an R object containing ... level
 * @param sdf            Pointer to an R object containing ... sdf
 * @param filter.type    Pointer to an R object containing ... filter.type
 * @param filter.length  Pointer to an R object containing ... filter.length
 * @see wavuniv_variance
 * @see wavuniv_variance_confidence
 * @see wavuniv_filters_daubechies
*/
EXTERN_R SEXP RS_wavelets_variance_edof(
 SEXP pr_interior,
 SEXP pr_num_coefs,
 SEXP pr_variance,
 SEXP pr_level,
 SEXP pr_sdf,
 SEXP pr_filter_type,
 SEXP pr_filter_length )
{
  SEXP             pr_ret_result;   
  mat_set          result;          
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           filter_length;   
  univ_mat         interior;        
  univ_mat         level;           
  univ_mat         num_coefs;       
  univ_mat         sdf;             
  univ_mat         variance;        
  void             *VPNULL = NULL;  
  wav_filter_type  filter_type;     
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_interior to interior */
  READ_MATRIX_REGISTER( pr_interior, &interior );

  /* ... pr_num_coefs to num_coefs */
  READ_MATRIX_REGISTER( pr_num_coefs, &num_coefs );

  /* ... pr_variance to variance */
  READ_MATRIX_REGISTER( pr_variance, &variance );

  /* ... pr_level to level */
  READ_MATRIX_REGISTER( pr_level, &level );

  /* ... pr_sdf to sdf */
  READ_MATRIX_REGISTER( pr_sdf, &sdf );

  /* ... pr_filter_type to filter_type */
  WAV_FILTER_TYPE_FROM_R( pr_filter_type, &filter_type );

  /* ... pr_filter_length to filter_length */
  SINT32_FROM_R( pr_filter_length, &filter_length );

  /* Call the function */
  err = wavuniv_variance_edof(
    &interior,
    &num_coefs,
    &variance,
    &level,
    &sdf,
    filter_type,
    filter_length,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( wavuniv_variance_confidence, &result, &pr_ret_result );
}

