
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_surr.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS fractal library.

   Functions wrapped:

   frauniv_bootstrap_theiler()
   frauniv_bootstrap_davison_hinkley()
   frauniv_bootstrap_circulant_embedding()
*/

#include "fra_surr.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Generate surrogate time series via Theiler's methods.
 * @source RS\_fra\_surr.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_bootstrap_theiler", time.series, method, seed))#
 * @return             An R ... containing ...
 * @param time.series  Pointer to an R object containing ... time.series
 * @param method       Pointer to an R object containing ... method
 * @param seed         Pointer to an R object containing ... seed
 * @see frauniv_bootstrap_davison_hinkley
 * @see frauniv_bootstrap_circulant_embedding
 * @see _fra_surrogate
 * @see frauniv_determinism_delta_epsilon
*/
EXTERN_R SEXP RS_fractal_bootstrap_theiler(
 SEXP pr_time_series,
 SEXP pr_method,
 SEXP pr_seed )
{
  SEXP             pr_ret_result;   
  fra_surrogate    method;          
  mutil_data_type  type;            
  mutil_errcode    err;             
  uint32           seed;            
  univ_mat         result;          
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

  /* ... pr_method to method */
  err = fra_surrogate_from_R( pr_method, &method );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert fra_surrogate type argument pr_method to method" );

  /* ... pr_seed to seed */
  err = sint32_from_R( pr_seed, (sint32 *) &seed );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert uint32 type argument pr_seed to seed" );

  /* Call the function */
  err = frauniv_bootstrap_theiler(
    &time_series,
    method,
    (uint32) seed,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_bootstrap_theiler, &result, &pr_ret_result );
}

/** Generate surrogate time series via Davison-Hinkley method.
 * @source RS\_fra\_surr.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_bootstrap_davison_hinkley", time.series, seed, seed))#
 * @return             An R ... containing ...
 * @param time.series  Pointer to an R object containing ... time.series
 * @param seed         Pointer to an R object containing ... seed
 * @see frauniv_bootstrap_theiler
 * @see frauniv_bootstrap_circulant_embedding
 * @see _fra_surrogate
 * @see frauniv_determinism_delta_epsilon
*/
EXTERN_R SEXP RS_fractal_bootstrap_davison_hinkley(
 SEXP pr_time_series,
 SEXP pr_seed )
{
  SEXP             pr_ret_result;   
  mutil_data_type  type;            
  mutil_errcode    err;             
  uint32           seed;            
  univ_mat         result;          
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

  /* ... pr_seed to seed */
  err = sint32_from_R( pr_seed, (sint32 *) &seed );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert uint32 type argument pr_seed to seed" );

  /* Call the function */
  err = frauniv_bootstrap_davison_hinkley(
    &time_series,
    (uint32) seed,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_bootstrap_davison_hinkley, &result, &pr_ret_result );
}

/** Generate surrogate time series via circulant embedding.
 * @source RS\_fra\_surr.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_bootstrap_circulant_embedding", sdf, seed, seed))#
 * @return      An R ... containing ...
 * @param sdf   Pointer to an R object containing ... sdf
 * @param seed  Pointer to an R object containing ... seed
 * @see frauniv_bootstrap_theiler
 * @see frauniv_bootstrap_davison_hinkley
 * @see _fra_surrogate
 * @see frauniv_spectral_density_function_direct
 * @see frauniv_spectral_density_function_lag_window
 * @see frauniv_spectral_density_function_wosa
 * @see frauniv_spectral_density_function_multitaper
 * @see frauniv_determinism_delta_epsilon
*/
EXTERN_R SEXP RS_fractal_bootstrap_circulant_embedding(
 SEXP pr_sdf,
 SEXP pr_seed )
{
  SEXP             pr_ret_result;   
  mutil_data_type  type;            
  mutil_errcode    err;             
  uint32           seed;            
  univ_mat         result;          
  univ_mat         sdf;             
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_sdf to sdf */
  READ_MATRIX_REGISTER( pr_sdf, &sdf );

  /* ... pr_seed to seed */
  err = sint32_from_R( pr_seed, (sint32 *) &seed );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert uint32 type argument pr_seed to seed" );

  /* Call the function */
  err = frauniv_bootstrap_circulant_embedding(
    &sdf,
    (uint32) seed,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_bootstrap_circulant_embedding, &result, &pr_ret_result );
}

