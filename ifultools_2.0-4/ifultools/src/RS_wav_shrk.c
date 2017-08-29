
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_wav_shrk.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS wavelets library.

   Functions wrapped:

   wavuniv_shrink()
   wavuniv_shrink_threshold()
*/

#include "wav_shrk.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Nonlinear noise reduction in a time series via wavelet shrinkage.
 * @source RS\_wav\_shrk.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_shrink", time.series, filters, threshold, threshold.function, threshold.scale, noise.variance, shrink.function, n.level, decimated))#
 * @return                    An R ... containing ...
 * @param time.series         Pointer to an R object containing ... time.series
 * @param filters             Pointer to an R object containing ... filters
 * @param threshold           Pointer to an R object containing ... threshold
 * @param threshold.function  Pointer to an R object containing ... threshold.function
 * @param threshold.scale     Pointer to an R object containing ... threshold.scale
 * @param noise.variance      Pointer to an R object containing ... noise.variance
 * @param shrink.function     Pointer to an R object containing ... shrink.function
 * @param n.level             Pointer to an R object containing ... n.level
 * @param decimated           Pointer to an R object containing ... decimated
 * @see _wav_shrink_function
 * @see _wav_shrink_threshold
 * @see wavuniv_shrink_threshold
 * @see wavuniv_transform_discrete_wavelet_convolution
 * @see wavuniv_transform_discrete_wavelet_convolution_inverse
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_transform_maximum_overlap_inverse
 * @see wavuniv_filters_daubechies
*/
EXTERN_R SEXP RS_wavelets_shrink(
 SEXP pr_time_series,
 SEXP pr_filters,
 SEXP pr_threshold,
 SEXP pr_threshold_function,
 SEXP pr_threshold_scale,
 SEXP pr_noise_variance,
 SEXP pr_shrink_function,
 SEXP pr_n_level,
 SEXP pr_decimated )
{
  SEXP                  pr_ret_result;        
  boolean               decimated;  
  boolean               threshold_defined;
  double                noise_variance;       
  double                threshold_scale;      
  mat_set               filters;              
  mutil_data_type       type;                 
  mutil_errcode         err;                  
  sint32                n_level;
  sint32                i;
  univ_mat              result;               
  univ_mat              threshold;            
  univ_mat              time_series;          
  void                  *VPNULL = NULL;       
  wav_shrink_function   shrink_function;      
  wav_shrink_threshold  threshold_function;   
  memlist            list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_filters to filters */
  READ_MATSET_REGISTER( pr_filters, MUTIL_DOUBLE, &filters );

  /* ... pr_threshold to threshold */
  READ_MATRIX_REGISTER( pr_threshold, &threshold );

  /* ... pr_threshold_function to threshold_function */
  WAV_SHRINK_THRESHOLD_FROM_R( pr_threshold_function, &threshold_function );

  /* ... pr_threshold_scale to threshold_scale */
  DOUBLE_FROM_R( pr_threshold_scale, &threshold_scale );

  /* ... pr_noise_variance to noise_variance */
  DOUBLE_FROM_R( pr_noise_variance, &noise_variance );

  /* ... pr_shrink_function to shrink_function */
  WAV_SHRINK_FUNCTION_FROM_R( pr_shrink_function, &shrink_function );

  /* ... pr_n_level to n_level */
  SINT32_FROM_R( pr_n_level, &n_level );

  /* ... pr_decimated to decimated */
  BOOLEAN_FROM_R( pr_decimated, &decimated );

  /* check for flag to indicate whether the thresholds
     are sent in explicitly or not */

  threshold_defined = (boolean) TRUE;

  for ( i = 0; i < MATUNIV_NELEM( &threshold ); i++ ){

    if ( threshold.mat.dblmat.data[i] < (double) 0.0 ){

      threshold_defined = (boolean) FALSE;
      break;
    }
  }

  /* Call the function */
  err = wavuniv_shrink(
    &time_series,
    &filters,
    (threshold_defined ? (&threshold) : ((univ_mat *) NULL)),
    threshold_function,
    threshold_scale,
    noise_variance,
    shrink_function,
    n_level,
    decimated,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( wavuniv_shrink, &result, &pr_ret_result );
}

/** Waveshrink thresholds.
 * @source RS\_wav\_shrk.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_shrink_threshold", wavelet.transform, decimated, threshold.type, shrink.fun, noise.variance, noise.variance, shrink.function, n.level, decimated))#
 * @return                   An R ... containing ...
 * @param wavelet.transform  Pointer to an R object containing ... wavelet.transform
 * @param decimated          Pointer to an R object containing ... decimated
 * @param threshold.type     Pointer to an R object containing ... threshold.type
 * @param shrink.fun         Pointer to an R object containing ... shrink.fun
 * @param noise.variance     Pointer to an R object containing ... noise.variance
 * @see _wav_shrink_function
 * @see _wav_shrink_threshold
 * @see wavuniv_transform_discrete_wavelet_convolution
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_shrink
*/
EXTERN_R SEXP RS_wavelets_shrink_threshold(
 SEXP pr_wavelet_transform,
 SEXP pr_decimated,
 SEXP pr_threshold_type,
 SEXP pr_shrink_fun,
 SEXP pr_noise_variance )
{
  SEXP                  pr_ret_result;        
  boolean               decimated;            
  double                noise_variance;       
  mat_set               wavelet_transform;    
  mutil_errcode         err;                  
  univ_mat              result;               
  void                  *VPNULL = NULL;       
  wav_shrink_function   shrink_fun;           
  wav_shrink_threshold  threshold_type;       
  memlist               list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_wavelet_transform to wavelet_transform */
  READ_MATSET_REGISTER( pr_wavelet_transform, MUTIL_DOUBLE, &wavelet_transform );

  /* ... pr_decimated to decimated */
  BOOLEAN_FROM_R( pr_decimated, &decimated );

  /* ... pr_threshold_type to threshold_type */
  WAV_SHRINK_THRESHOLD_FROM_R( pr_threshold_type, &threshold_type );

  /* ... pr_shrink_fun to shrink_fun */
  WAV_SHRINK_FUNCTION_FROM_R( pr_shrink_fun, &shrink_fun );

  /* ... pr_noise_variance to noise_variance */
  DOUBLE_FROM_R( pr_noise_variance, &noise_variance );

  /* Call the function */
  err = wavuniv_shrink_threshold(
    &wavelet_transform,
    decimated,
    threshold_type,
    shrink_fun,
    noise_variance,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( wavuniv_shrink_threshold, &result, &pr_ret_result );
}

