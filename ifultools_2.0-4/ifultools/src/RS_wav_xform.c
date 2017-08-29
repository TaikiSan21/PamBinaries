
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_wav_xform.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS wavelets library.

   Functions wrapped:

   wavuniv_transform_continuous_wavelet()
   wavuniv_transform_discrete_wavelet_convolution()
   wavuniv_transform_packet()
   wavuniv_transform_discrete_wavelet_convolution_inverse()
   wavuniv_transform_packet_convert_indices()
   wavuniv_transform_packet_basis()
   wavuniv_transform_packet_inverse()
   wavuniv_transform_maximum_overlap()
   wavuniv_transform_maximum_overlap_packet()
   wavuniv_transform_maximum_overlap_inverse()
 */

#include "wav_modw.h"
#include "wav_dwtc.h"
#include "wav_cwt.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** The continuous wavelet transform.
 * @source RS\_wav\_xform.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_continuous_wavelet", time.series, sampling.interval, filter.type, filter.arg, scale))#
 * @return                   An R ... containing ...
 * @param time.series        Pointer to an R object containing ... time.series
 * @param sampling.interval  Pointer to an R object containing ... sampling.interval
 * @param filter.type        Pointer to an R object containing ... filter.type
 * @param filter.arg         Pointer to an R object containing ... filter.arg
 * @param scale              Pointer to an R object containing ... scale
 * @see _wav_filter_type
 * @see wavuniv_filters_continuous
*/
EXTERN_R SEXP RS_wavelets_transform_continuous_wavelet(
 SEXP pr_time_series,
 SEXP pr_sampling_interval,
 SEXP pr_filter_type,
 SEXP pr_filter_arg,
 SEXP pr_scale )
{
  SEXP              pr_ret_result;         
  double            filter_arg;         
  double            sampling_interval;  
  mutil_data_type   type;               
  mutil_errcode     err;                
  univ_mat          result;                
  univ_mat          scale;              
  univ_mat          time_series;        
  void              *VPNULL = NULL;     
  wav_filter_type   filter_type;        
  memlist           list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_sampling_interval to sampling_interval */
  DOUBLE_FROM_R( pr_sampling_interval, &sampling_interval );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert double type argument pr_sampling_interval to sampling_interval" );

  /* ... pr_filter_type to filter_type */
  WAV_FILTER_TYPE_FROM_R( pr_filter_type, &filter_type );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert wav_filter_type type argument pr_filter_type to filter_type" );

  /* ... pr_filter_arg to filter_arg */
  DOUBLE_FROM_R( pr_filter_arg, &filter_arg );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert double type argument pr_filter_arg to filter_arg" );

  /* ... pr_scale to scale */
  READ_MATRIX_REGISTER( pr_scale, &scale );


  /* Call the function */
  err = wavuniv_transform_continuous_wavelet(
    &time_series,
    sampling_interval,
    filter_type,
    filter_arg,
    &scale,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( wavuniv_transform_continuous_wavelet, &result, &pr_ret_result );
}

/** The discrete wavelet transform using convolution style filtering.
 * @source RS\_wav\_dwtc.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_discrete_wavelet_convolution", time.series, filters, n.level))#
 * @return             An R ... containing ...
 * @param time.series  Pointer to an R object containing ... time.series
 * @param filters      Pointer to an R object containing ... filters
 * @param n.level      Pointer to an R object containing ... n.level
 * @see _wav_filter_type
 * @see wavuniv_transform_discrete_wavelet
 * @see wavuniv_transform_discrete_wavelet_inverse
 * @see wavuniv_transform_discrete_wavelet_convolution_inverse
 * @see wavuniv_coefficient_zero_phase
 * @see wavuniv_coefficient_boundaries
 * @see wavuniv_transform_discrete_wavelet
 * @see wavuniv_transform_maximum_overlap
*/
EXTERN_R SEXP RS_wavelets_transform_discrete_wavelet_convolution(
 SEXP pr_time_series,
 SEXP pr_filters,
 SEXP pr_n_level )
{
  SEXP             pr_ret_result;   
  mat_set          filters;         
  mat_set          result;          
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           n_level;         
  univ_mat         time_series;     
  void             *VPNULL = NULL;  
  memlist           list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_filters to filters */
  READ_MATSET_REGISTER( pr_filters, MUTIL_DOUBLE, &filters );

  /* ... pr_n_level to n_level */
  SINT32_FROM_R( pr_n_level, &n_level );

  /* Call the function */
  err = wavuniv_transform_discrete_wavelet_convolution(
    &time_series,
    &filters,
    n_level,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( wavuniv_transform_discrete_wavelet_convolution, &result, &pr_ret_result );
}

/** The discrete wavelet packet transform using convolution style filtering.
 * @source RS\_wav\_dwtc.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_packet", time.series, filters, n.level))#
 * @return             An R ... containing ...
 * @param time.series  Pointer to an R object containing ... time.series
 * @param filters      Pointer to an R object containing ... filters
 * @param n.level      Pointer to an R object containing ... n.level
 * @see _wav_filter_type
 * @see wavuniv_transform_discrete_wavelet
 * @see wavuniv_transform_maximum_overlap_packet
 * @see wavuniv_transform_packet_detail
*/
EXTERN_R SEXP RS_wavelets_transform_packet(
 SEXP pr_time_series,
 SEXP pr_filters,
 SEXP pr_n_level )
{
  SEXP             pr_ret_result;   
  mat_set          filters;         
  mat_set          result;          
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           n_level;         
  univ_mat         time_series;     
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_filters to filters */
  READ_MATSET_REGISTER( pr_filters, MUTIL_DOUBLE, &filters );

  /* ... pr_n_level to n_level */
  SINT32_FROM_R( pr_n_level, &n_level );

  /* Call the function */
  err = wavuniv_transform_packet(
    &time_series,
    &filters,
    n_level,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( wavuniv_transform_packet, &result, &pr_ret_result );
}

/** The inverse discrete wavelet transform using convolution style filtering.
 * @source RS\_wav\_dwtc.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_discrete_wavelet_convolution_inverse", dwt, filters, n.level))#
 * @return         An R ... containing ...
 * @param dwt      Pointer to an R object containing ... dwt
 * @param filters  Pointer to an R object containing ... filters
 * @see _wav_filter_type
 * @see wavuniv_transform_discrete_wavelet
 * @see wavuniv_transform_discrete_wavelet_inverse
 * @see wavuniv_transform_discrete_wavelet_convolution
*/
EXTERN_R SEXP RS_wavelets_transform_discrete_wavelet_convolution_inverse(
 SEXP pr_dwt,
 SEXP pr_filters )
{
  SEXP             pr_ret_result;   
  mat_set          dwt;             
  mat_set          filters;         
  mutil_errcode    err;             
  univ_mat         result;          
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_dwt to dwt */
  READ_MATSET_REGISTER( pr_dwt, MUTIL_DOUBLE, &dwt );

  /* ... pr_filters to filters */
  READ_MATSET_REGISTER( pr_filters, MUTIL_DOUBLE, &filters );

  /* Call the function */
  err = wavuniv_transform_discrete_wavelet_convolution_inverse(
    &dwt,
    &filters,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( wavuniv_transform_discrete_wavelet_convolution_inverse, &result, &pr_ret_result );
}

/** Conversion and validation of wavelet packet indices.
 * @source RS\_wav\_dwtc.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_packet_convert_indices", transform.indices, filters, n.level))#
 * @return                   An R ... containing ...
 * @param transform.indices  Pointer to an R object containing ... transform.indices
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_maximum_overlap_packet
*/
EXTERN_R SEXP RS_wavelets_transform_packet_convert_indices(
 SEXP pr_transform_indices )
{
  SEXP                pr_ret_flat;        
  SEXP                pr_ret_level;       
  SEXP                pr_ret_obj;         
  SEXP                pr_ret_osc;         
  mutil_data_type     type;               
  mutil_errcode       err;                
  univ_mat            flat;               
  univ_mat            level;              
  univ_mat            osc;                
  univ_mat            transform_indices;  
  void                *VPNULL = NULL;     
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_transform_indices to transform_indices */
  READ_MATRIX_REGISTER( pr_transform_indices, &transform_indices );

  /* Call the function */
  err = wavuniv_transform_packet_convert_indices(
    &transform_indices,
    VPNULL,
    &flat,
    &level,
    &osc );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Problem calling wavuniv_transform_packet_convert_indices() function" );
  err = memlist_member_register( &list, &flat, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
  err = memlist_member_register( &list, &level, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
  err = memlist_member_register( &list, &osc, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );

  /* create the output R object */

  err = matuniv_to_R( &flat, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_flat );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  err = matuniv_to_R( &level, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_level );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  err = matuniv_to_R( &osc, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_osc );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  PROTECT( pr_ret_obj = allocVector( VECSXP, 3 ) );
  SET_VECTOR_ELT( pr_ret_obj, 0, pr_ret_flat );
  SET_VECTOR_ELT( pr_ret_obj, 1, pr_ret_level );
  SET_VECTOR_ELT( pr_ret_obj, 2, pr_ret_osc );
  UNPROTECT(1);

  /* free registered local memory */
  MUTIL_FREE_WARN( memlist, &list );

  return pr_ret_obj;
}

/** Extracts a discrete wavelet packet transform subset.
 * @source RS\_wav\_dwtc.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_packet_basis", dwpt, transform.indices, n.level))#
 * @return                   An R ... containing ...
 * @param dwpt               Pointer to an R object containing ... dwpt
 * @param transform.indices  Pointer to an R object containing ... transform.indices
 * @see wavuniv_transform_packet_convert_indices
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_maximum_overlap_packet
 * @see _wav_dwpt_extra
 * @see wavuniv_transform_packet_inverse
*/
EXTERN_R SEXP RS_wavelets_transform_packet_basis(
 SEXP pr_dwpt,
 SEXP pr_transform_indices )
{
  SEXP             pr_ret_obj;         
  SEXP             tmpobj;
  mat_set          dwpt;               
  mat_set          result;             
  mutil_data_type  type;               
  mutil_errcode    err;                
  univ_mat         transform_indices;  
  void            *VPNULL = NULL;     
  wav_dwpt_extra   extra;  
  univ_mat         extra_atoms;
  univ_mat         extra_levelmap;
  sint32           i;
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_dwpt to dwpt */
  READ_MATSET_REGISTER( pr_dwpt, MUTIL_DOUBLE, &dwpt );

  /* ... pr_transform_indices to transform_indices */
  READ_MATRIX_REGISTER( pr_transform_indices, &transform_indices );

  /* Call the function */
  err = wavuniv_transform_packet_basis(
    &dwpt,
    &transform_indices,
    VPNULL,
    &result,
    &extra );
  MUTIL_FREE_WARN( memlist, &list );
  if ( err ){
    PROBLEM "Problem calling wavuniv_transform_packet_basis() function" ERROR;
  }

  /* in the case where there are extra dwpt atoms that need to
     be stored, there will be memory allocated for the
     'extra' structure above. the easiest way to get this
     back to R is to pack everything into a list */

  if ( extra.nelem > 0 ){

    /* wrap universal matrix headers around extra matrices */
    err = matuniv_wrap_matrix( &extra_atoms, &(extra.atoms), MUTIL_DOUBLE);
    if ( err ){
      MUTIL_ERROR( "Unable to wrap universal matrix around extra atoms matrix" );
      PROBLEM "Problem creating R list" ERROR;
    }

    err = matuniv_wrap_matrix( &extra_levelmap, &(extra.levelmap), MUTIL_SINT32 );
    if ( err ){
      MUTIL_ERROR( "Unable to wrap universal matrix around extra levelmap matrix" );
      PROBLEM "Problem creating R list" ERROR;
    }

    /* create an R list object */
    PROTECT( pr_ret_obj = allocVector( VECSXP, result.nelem + 2 ) );

    /* create R matrix objects and pack the extra crystals into the list */

    for( i = 0; i < result.nelem; i++ ) {

      err = matuniv_to_R( &(result.mats[i]), (mutil_R_class_type) MUTIL_R_MATRIX, &tmpobj );
      if( err ){
	      PROBLEM "Problem adding DWPT crystal to R list" ERROR;
      }

      /*LINTED: cast OK, checked range */
      SET_VECTOR_ELT( pr_ret_obj, (int) i, tmpobj);
    }

    /* pack the extra atoms crystal into the R list */
    err = matuniv_to_R( &extra_atoms, (mutil_R_class_type) MUTIL_R_MATRIX, &tmpobj );
    if( err ){
      PROBLEM "Problem adding extra DWPT crystal to R list" ERROR;
    }

    /*LINTED: cast OK, checked range */
    SET_VECTOR_ELT( pr_ret_obj, (int) i++, tmpobj );

    /* finally, pack the levelmap vector into the list */
    err = matuniv_to_R( &extra_levelmap, (mutil_R_class_type) MUTIL_R_MATRIX, &tmpobj );
    if( err ){
      PROBLEM "Problem adding extra DWPT levelmap crystal to R list" ERROR;
    }

    /*LINTED: cast OK, checked range */
    SET_VECTOR_ELT( pr_ret_obj, (int) i, tmpobj );

    /* free the memory */

    MUTIL_FREE_WARN( matuniv, &extra_levelmap );
    MUTIL_FREE_WARN( matuniv, &extra_atoms );
    MUTIL_FREEALL_MATSET_WARN( &result );

    UNPROTECT(1);
  }
  else{

    err = matset_to_R_list( &result, &pr_ret_obj );
    MUTIL_FREEALL_MATSET_WARN( &result );
    if ( err ) {
      PROBLEM "Unable to convert output data to Splus format" ERROR;
    }
  }

  return pr_ret_obj;
}

/** Discrete wavelet packet transform subset inverse.
 * @source RS\_wav\_dwtc.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_packet_inverse", dwpt.basis, extra, transform.indices, filters))#
 * @return                   An R ... containing ...
 * @param dwpt.basis         Pointer to an R object containing ... dwpt.basis
 * @param extra              Pointer to an R object containing ... extra
 * @param atoms              Any extra DWPT atoms.
 * @param levelmap           The level map for extra DWPT atoms.
 * @param transform.indices  Pointer to an R object containing ... transform.indices
 * @param filters            Pointer to an R object containing ... filters
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_packet_convert_indices
*/
EXTERN_R SEXP RS_wavelets_transform_packet_inverse(
 SEXP pr_dwpt_basis,
 SEXP pr_nextra,
 SEXP pr_atoms,
 SEXP pr_levelmap,
 SEXP pr_transform_indices,
 SEXP pr_filters )
{
  SEXP                pr_ret_result;      
  mat_set             dwpt_basis;         
  mat_set             filters;            
  mutil_data_type     type;               
  mutil_errcode       err;                
  univ_mat            result;             
  univ_mat            transform_indices;  
  void                *VPNULL = NULL;     
  wav_dwpt_extra      extra;              
  univ_mat            um_atoms;
  univ_mat            um_levelmap;
  boolean             any_extra;

  /* Avoid lint warning */

  (void) whatssi;

  /* Conversion of input data ... */

  /* ... pr_extra to extra */
  err = sint32_from_R( pr_nextra, &(extra.nelem) );
  if ( err ) {
    PROBLEM "Unable to convert pr_nextra to sint32 value" ERROR;
  }

  any_extra = (boolean) ( extra.nelem > 0 );

  if ( any_extra ){

    err = matuniv_from_R( pr_atoms, MUTIL_DOUBLE, &um_atoms );
    if ( err ){
      PROBLEM "Unable to read pr_atoms" ERROR;
    }

    err = matuniv_from_R( pr_levelmap, MUTIL_SINT32, &um_levelmap );
    if ( err ){
      PROBLEM "Unable to read pr_levelmap" ERROR;
    }

    extra.nelem = MATUNIV_NELEM( &um_atoms );
    extra.atoms = um_atoms.mat.dblmat;
    extra.levelmap = um_levelmap.mat.s32mat;
  }

  /* ... pr_dwpt_basis to dwpt_basis */
  err = matset_from_R( pr_dwpt_basis, MUTIL_DOUBLE, &dwpt_basis );
  if ( err ){
      PROBLEM "Unable to read pr_dwpt_basis" ERROR;
  }

  /* ... pr_transform_indices to transform_indices */
  err = mutil_R_type( pr_transform_indices, &type );
  if ( err ){
      PROBLEM "Unable to read pr_transform_indices type" ERROR;
  }

  err = matuniv_from_R( pr_transform_indices, type, &transform_indices );
  if ( err ){
      PROBLEM "Unable to read pr_transform_indices" ERROR;
  }

  /* ... pr_filters to filters */
  err = matset_from_R( pr_filters, MUTIL_DOUBLE, &filters );
  if ( err ){
      PROBLEM "Unable to read pr_filters" ERROR;
  }

  /* Call the function */
  err = wavuniv_transform_packet_inverse(
    &dwpt_basis,
    &extra,
    &transform_indices,
    &filters,
    VPNULL,
    &result );
  if ( err ){
    PROBLEM "Problem calling wavuniv_transform_packet_inverse() function" ERROR;
  }

  if ( any_extra ){
    MUTIL_FREE_WARN( matuniv, &um_atoms );
    MUTIL_FREE_WARN( matuniv, &um_levelmap );
  }

  /* create the output R object */
  err = matuniv_to_R( &result, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_result );
  MUTIL_FREE_WARN( matuniv, &result );
  if ( err ) {
      PROBLEM "Unable to convert output data to R format" ERROR;
  }

  return pr_ret_result;
}

/** The maximum overlap discrete wavelet transform (MODWT).
 * @source RS\_wav\_modw.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_maximum_overlap", time.series, filters, n.level))#
 * @return             An R ... containing ...
 * @param time.series  Pointer to an R object containing ... time.series
 * @param filters      Pointer to an R object containing ... filters
 * @param n.level      Pointer to an R object containing ... n.level
 * @see wavuniv_filters_daubechies
 * @see wavuniv_transform_maximum_overlap_inverse
 * @see wavuniv_transform_maximum_overlap_packet
 * @see wavuniv_transform_packet_detail
*/
EXTERN_R SEXP RS_wavelets_transform_maximum_overlap(
 SEXP pr_time_series,
 SEXP pr_filters,
 SEXP pr_n_level )
{
  SEXP             pr_ret_result;   
  mat_set          filters;         
  mat_set          result;          
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           n_level;         
  univ_mat         time_series;     
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_filters to filters */
  READ_MATSET_REGISTER( pr_filters, MUTIL_DOUBLE, &filters );

  /* ... pr_n_level to n_level */
  SINT32_FROM_R( pr_n_level, &n_level );

  /* Call the function */
  err = wavuniv_transform_maximum_overlap(
    &time_series,
    &filters,
    n_level,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( wavuniv_transform_maximum_overlap, &result, &pr_ret_result );
}

/** The maximum overlap discrete wavelet packet transform (MODWPT).
 * @source RS\_wav\_modw.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_maximum_overlap_packet", time.series, filters, n.level))#
 * @return             An R ... containing ...
 * @param time.series  Pointer to an R object containing ... time.series
 * @param filters      Pointer to an R object containing ... filters
 * @param n.level      Pointer to an R object containing ... n.level
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_transform_maximum_overlap_inverse
 * @see wavuniv_transform_packet_detail
*/
EXTERN_R SEXP RS_wavelets_transform_maximum_overlap_packet(
 SEXP pr_time_series,
 SEXP pr_filters,
 SEXP pr_n_level )
{
  SEXP             pr_ret_result;   
  mat_set          filters;         
  mat_set          result;          
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           n_level;         
  univ_mat         time_series;     
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_filters to filters */
  READ_MATSET_REGISTER( pr_filters, MUTIL_DOUBLE, &filters );

  /* ... pr_n_level to n_level */
  SINT32_FROM_R( pr_n_level, &n_level );

  /* Call the function */
  err = wavuniv_transform_maximum_overlap_packet(
    &time_series,
    &filters,
    n_level,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( wavuniv_transform_maximum_overlap_packet, &result, &pr_ret_result );
}

/** The inverse maximum overlap discrete wavelet transform (IMODWT).
 * @source RS\_wav\_modw.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_maximum_overlap_inverse", modwt, filters, level, node, xformtype))#
 * @return         An R ... containing ...
 * @param modwt    Pointer to an R object containing ... modwt
 * @param filters  Pointer to an R object containing ... filters
 * @see wavuniv_filters_daubechies
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_transform_maximum_overlap_packet
 * @see wavuniv_transform_packet_detail
*/
EXTERN_R SEXP RS_wavelets_transform_maximum_overlap_inverse(
 SEXP pr_modwt,
 SEXP pr_filters )
{
  SEXP             pr_ret_result;   
  mat_set          filters;         
  mat_set          modwt;           
  mutil_errcode    err;             
  univ_mat         result;          
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_modwt to modwt */
  READ_MATSET_REGISTER( pr_modwt, MUTIL_DOUBLE, &modwt );

  /* ... pr_filters to filters */
  READ_MATSET_REGISTER( pr_filters, MUTIL_DOUBLE, &filters );

  /* Call the function */
  err = wavuniv_transform_maximum_overlap_inverse(
    &modwt,
    &filters,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( wavuniv_transform_maximum_overlap_inverse, &result, &pr_ret_result );
}

