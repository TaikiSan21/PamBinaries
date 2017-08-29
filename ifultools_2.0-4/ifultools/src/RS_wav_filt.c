
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_wav_filt.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS wavelets library.

   Functions wrapped:

   wavuniv_filters_daubechies
   wavuniv_filters_daubechies_verify
   wavuniv_filters_daubechies_gain
   wavuniv_filters_zero_phase
   wavuniv_filters_continuous
   wavuniv_transform_coefficient_boundaries
*/

#include "wav_filt.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "RS_mac.h"

/** Daubechies wavelet and scaling filters.
 * @source RS\_wav\_filt.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_filters_daubechies", filter.length, filter.type, normalize))#
 * @return               An R ... containing ...
 * @param filter.length  Pointer to an R object containing ... filter.length
 * @param filter.type    Pointer to an R object containing ... filter.type
 * @param normalize      Pointer to an R object containing ... normalize
 * @see _wav_filter_type
 * @see wavuniv_filters_daubechies_verify
 * @see wavuniv_transform_discrete_wavelet_convolution
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_filters_daubechies_gain
 * @see wavuniv_filters_zero_phase
 * @see wavuniv_coefficient_boundaries
*/
EXTERN_R SEXP RS_wavelets_filters_daubechies(
 SEXP pr_filter_length,
 SEXP pr_filter_type,
 SEXP pr_normalize )
{
  SEXP             pr_ret_result;
  boolean          normalize;
  mat_set          result;
  mutil_errcode    err;
  sint32           filter_length;
  void             *VPNULL = NULL;
  wav_filter_type  filter_type;

  /* Avoid lint warning */
  (void) whatssi;

  /* Conversion of input data ... */

  /* ... pr_filter_length to filter_length */
  err = sint32_from_R( pr_filter_length, &filter_length );
  if ( err ){
    PROBLEM "Unable to convert sint32 type argument pr_filter_length to filter_length" ERROR;
  }

  /* ... pr_filter_type to filter_type */
  err = wav_filter_type_from_R( pr_filter_type, &filter_type );
  if ( err ){
    PROBLEM "Unable to convert wav_filter_type type argument pr_filter_type to filter_type" ERROR;
  }

  /* ... pr_normalize to normalize */
  err = boolean_from_R( pr_normalize, &normalize );
  if ( err ){
    PROBLEM "Unable to convert boolean type argument pr_normalize to normalize" ERROR;
  }

  /* Call the function */
  err = wavuniv_filters_daubechies(
    filter_length,
    filter_type,
    normalize,
    VPNULL,
    &result );
  if ( err ){
    PROBLEM "Problem calling wavuniv_filters_daubechies() function" ERROR;
  }

  /* create the output R object */
  err = matset_to_R_list( &result, &pr_ret_result );
  MUTIL_FREEALL_MATSET_WARN( &result );
  if ( err ) {
      PROBLEM "Unable to convert output data to R format" ERROR;
  }

  return pr_ret_result;
}

/** The gain functions for Daubechies wavelet and scaling filters.
 * @source RS\_wav\_filt.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_filters_daubechies_gain", filter.type, filter.length, num.levels, num.fft, normalize))#
 * @return               An R ... containing ...
 * @param filter.type    Pointer to an R object containing ... filter.type
 * @param filter.length  Pointer to an R object containing ... filter.length
 * @param num.levels     Pointer to an R object containing ... num.levels
 * @param num.fft        Pointer to an R object containing ... num.fft
 * @param normalize      Pointer to an R object containing ... normalize
 * @see _wav_filter_type
 * @see wavuniv_filters_daubechies
 * @see wavuniv_filters_zero_phase
*/
EXTERN_R SEXP RS_wavelets_filters_daubechies_gain(
 SEXP pr_filter_type,
 SEXP pr_filter_length,
 SEXP pr_num_levels,
 SEXP pr_num_fft,
 SEXP pr_normalize )
{
  SEXP             pr_ret_gain_frequency;
  SEXP             pr_ret_gain_scaling;
  SEXP             pr_ret_gain_wavelet;
  SEXP             pr_ret_obj;
  boolean          normalize;
  mutil_errcode    err;
  sint32           filter_length;
  sint32           num_fft;
  sint32           num_levels;
  univ_mat         gain_frequency;
  univ_mat         gain_scaling;
  univ_mat         gain_wavelet;
  void             *VPNULL = NULL;
  wav_filter_type  filter_type;
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_filter_type to filter_type */
  WAV_FILTER_TYPE_FROM_R( pr_filter_type, &filter_type );

  /* ... pr_filter_length to filter_length */
  SINT32_FROM_R( pr_filter_length, &filter_length );

  /* ... pr_num_levels to num_levels */
  SINT32_FROM_R( pr_num_levels, &num_levels );

  /* ... pr_num_fft to num_fft */
  SINT32_FROM_R( pr_num_fft, &num_fft );

  /* ... pr_normalize to normalize */
  BOOLEAN_FROM_R( pr_normalize, &normalize );

  /* Call the function */
  err = wavuniv_filters_daubechies_gain(
    filter_type,
    filter_length,
    num_levels,
    num_fft,
    normalize,
    VPNULL,
    &gain_frequency,
    &gain_wavelet,
    &gain_scaling );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Problem calling wavuniv_filters_daubechies_gain() function" );
  err = memlist_member_register( &list, &gain_frequency, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
  err = memlist_member_register( &list, &gain_wavelet, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
  err = memlist_member_register( &list, &gain_scaling, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );

  /* create the output R object */

  err = matuniv_to_R( &gain_frequency, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_gain_frequency );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  err = matuniv_to_R( &gain_wavelet, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_gain_wavelet );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  err = matuniv_to_R( &gain_scaling, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_gain_scaling );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  PROTECT( pr_ret_obj = allocVector( VECSXP, 3 ) );
  SET_VECTOR_ELT( pr_ret_obj, 0, pr_ret_gain_frequency );
  SET_VECTOR_ELT( pr_ret_obj, 1, pr_ret_gain_wavelet );
  SET_VECTOR_ELT( pr_ret_obj, 2, pr_ret_gain_scaling );
  UNPROTECT(1);

  /* free registered local memory */
  MUTIL_FREE_WARN( memlist, &list );

  return pr_ret_obj;
}

/** Zero phase shift factors for Daubechies symmlet and Coiflet filters.
 * @source RS\_wav\_filt.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_filter_zero_phase", filter.type, filter.length, n.level, transform.type, normalize))#
 * @return               An R ... containing ...
 * @param filter.type    Pointer to an R object containing ... filter.type
 * @param filter.length  Pointer to an R object containing ... filter.length
 * @param n.level        Pointer to an R object containing ... n.level
 * @see _wav_filter_type
 * @see wavuniv_coefficient_boundaries
 * @see wavuniv_filters_daubechies
 * @see wavuniv_transform_maximum_overlap
*/
EXTERN_R SEXP RS_wavelets_filter_zero_phase(
 SEXP pr_filter_type,
 SEXP pr_filter_length,
 SEXP pr_n_level )
{
  SEXP             pr_ret_result;
  mat_set          result;
  mutil_errcode    err;
  sint32           filter_length;
  sint32           n_level;
  void             *VPNULL = NULL;
  wav_filter_type  filter_type;

  /* Avoid lint warning */
  (void) whatssi;

  /* Conversion of input data ... */

  /* ... pr_filter_type to filter_type */
  err = wav_filter_type_from_R( pr_filter_type, &filter_type );
  if ( err ){
    PROBLEM "Unable to convert wav_filter_type type argument pr_filter_type to filter_type" ERROR;
  }

  /* ... pr_filter_length to filter_length */
  err = sint32_from_R( pr_filter_length, &filter_length );
  if ( err ){
    PROBLEM "Unable to convert sint32 type argument pr_filter_length to filter_length" ERROR;
  }

  /* ... pr_n_level to n_level */
  err = sint32_from_R( pr_n_level, &n_level );
  if ( err ){
    PROBLEM "Unable to convert sint32 type argument pr_n_level to n_level" ERROR;
  }

  /* Call the function */
  err = wavuniv_filters_zero_phase(
    filter_type,
    filter_length,
    n_level,
    VPNULL,
    &result );
  if ( err ){
    PROBLEM "Problem calling wavuniv_filters_zero_phase() function" ERROR;
  }

  /* create the output R object */

  err = matset_to_R_list( &result, &pr_ret_result );
  MUTIL_FREEALL_MATSET_WARN( &result );
  if ( err ) {
    PROBLEM "Unable to convert output data to R format" ERROR;
  }

  return pr_ret_result;
}

/** Creates frequency domain filters for the continuous wavelet transform.
 * @source RS\_wav\_filt.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_filters_continuous", filter.type, filter.arg, frequency, transform.type, normalize))#
 * @return             An R ... containing ...
 * @param filter.type  Pointer to an R object containing ... filter.type
 * @param filter.arg   Pointer to an R object containing ... filter.arg
 * @param frequency    Pointer to an R object containing ... frequency
 * @see _wav_filter_type
 * @see wavuniv_transform_continuous_wavelet
*/
EXTERN_R SEXP RS_wavelets_filters_continuous(
 SEXP pr_filter_type,
 SEXP pr_filter_arg,
 SEXP pr_frequency )
{
  SEXP             pr_ret_result;
  double           filter_arg;
  mutil_data_type  type;
  mutil_errcode    err;
  univ_mat         frequency;
  univ_mat         result;
  void             *VPNULL = NULL;
  wav_filter_type  filter_type;

  /* Avoid lint warning */
  (void) whatssi;

  /* Conversion of input data ... */

  /* ... pr_filter_type to filter_type */
  err = wav_filter_type_from_R( pr_filter_type, &filter_type );
  if ( err ){
    PROBLEM "Unable to convert wav_filter_type type argument pr_filter_type to filter_type" ERROR;
  }

  /* ... pr_filter_arg to filter_arg */
  err = double_from_R( pr_filter_arg, &filter_arg );
  if ( err ){
    PROBLEM "Unable to convert double type argument pr_filter_arg to filter_arg" ERROR;
  }

  /* ... pr_frequency to frequency */
  err = mutil_R_type( pr_frequency, &type );
  if ( err ){
    PROBLEM "Unable to read pr_frequency type" ERROR;
  }

  err = matuniv_from_R( pr_frequency, type, &frequency );
  if ( err ){
    PROBLEM "Unable to read pr_frequency" ERROR;
  }

  err = matuniv_malloc( &result,
    MATUNIV_NROW( &frequency ),
    MATUNIV_NCOL( &frequency ),
    MUTIL_DCOMPLEX );
  if ( err ) {
    MUTIL_FREE_WARN( matuniv, &frequency );
    PROBLEM "Unable to allocate memory for frequency response matrix" ERROR;
  }

  /* Call the function */
  err = wavuniv_filters_continuous(
    filter_type,
    filter_arg,
    &frequency,
    VPNULL,
    &result );
  MUTIL_FREE_WARN( matuniv, &frequency );
  if ( err ){
    PROBLEM "Problem calling wavuniv_filters_continuous() function" ERROR;
  }

  /* create the output R object */

  err = matuniv_to_R( &result, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_result );
  MUTIL_FREE_WARN( matuniv, &result );
  if ( err ) {
      PROBLEM "Unable to convert output data to R format" ERROR;
  }

  return pr_ret_result;
}

/** Boundary and interior wavelet coefficient identification for the DWT and MODWT.
 * @source RS\\_wav\_coef.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_coefficient_boundaries", n.level, filter.length, n.sample, transform.type, normalize))#
 * @return                An R ... containing ...
 * @param n.level         Pointer to an R object containing ... n.level
 * @param filter.length   Pointer to an R object containing ... filter.length
 * @param n.sample        Pointer to an R object containing ... n.sample
 * @param transform.type  Pointer to an R object containing ... transform.type
 * @see wavuniv_filters_zero_phase
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_filters_daubechies
*/
EXTERN_R SEXP RS_wavelets_transform_coefficient_boundaries(
 SEXP pr_n_level,
 SEXP pr_filter_length,
 SEXP pr_n_sample,
 SEXP pr_transform_type )
{
  SEXP             pr_ret_result;
  mat_set          result;
  mutil_errcode    err;
  sint32           filter_length;
  sint32           n_level;
  sint32           n_sample;
  void             *VPNULL = NULL;
  wav_transform    transform_type;

  /* Avoid lint warning */
  (void) whatssi;

  /* Conversion of input data ... */

  /* ... pr_n_level to n_level */
  err = sint32_from_R( pr_n_level, &n_level );
  if ( err ){
    PROBLEM "Unable to convert sint32 type argument pr_n_level to n_level" ERROR;
  }

  /* ... pr_filter_length to filter_length */
  err = sint32_from_R( pr_filter_length, &filter_length );
  if ( err ){
    PROBLEM "Unable to convert sint32 type argument pr_filter_length to filter_length" ERROR;
  }

  /* ... pr_n_sample to n_sample */
  err = sint32_from_R( pr_n_sample, &n_sample );
  if ( err ){
    PROBLEM "Unable to convert sint32 type argument pr_n_sample to n_sample" ERROR;
  }

  /* ... pr_transform_type to transform_type */
  err = wav_transform_from_R( pr_transform_type, &transform_type );
  if ( err ){
    PROBLEM "Unable to convert wav_transform type argument pr_transform_type to transform_type" ERROR;
  }

  /* Call the function */
  err = wavuniv_transform_coefficient_boundaries(
    n_level,
    filter_length,
    n_sample,
    transform_type,
    VPNULL,
    &result );
  if ( err ){
    PROBLEM "Problem calling wavelets_transform_coefficient_boundaries() function" ERROR;
  }

  /* create the output R object */
  err = matset_to_R_list( &result, &pr_ret_result );
  MUTIL_FREEALL_MATSET_WARN( &result );
  if ( err ) {
    PROBLEM "Unable to convert output data to R format" ERROR;
  }

  return pr_ret_result;
}
