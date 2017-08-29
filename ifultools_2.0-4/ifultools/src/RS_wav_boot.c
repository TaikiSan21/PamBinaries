
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_wav_boot.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS wavelets library.

   Functions wrapped:

   wavuniv_transform_packet_whitest()
   wavuniv_bootstrap()
*/

#include "wav_boot.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Find the whitest set of DWPT crystals.
 * @source RS\_wav\_boot.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_packet_whitest", dwpt, significance, white.noise.test))#
 * @return                  An R ... containing ...
 * @param dwpt              Pointer to an R object containing ... dwpt
 * @param significance      Pointer to an R object containing ... significance
 * @param white.noise.test  Pointer to an R object containing ... white.noise.test
 * @see _wav_white_test
 * @see wavuniv_bootstrap
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_packet_inverse
*/
EXTERN_R SEXP RS_wavelets_transform_packet_whitest(
 SEXP pr_dwpt,
 SEXP pr_significance,
 SEXP pr_white_noise_test )
{
  SEXP               pr_ret_result;     
  double             significance;      
  mat_set            dwpt;              
  mutil_errcode      err;               
  univ_mat           result;            
  void               *VPNULL = NULL;    
  wav_white_test     white_noise_test;
  memlist            list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_dwpt to dwpt */
  READ_MATSET_REGISTER( pr_dwpt, MUTIL_DOUBLE, &dwpt );

  /* ... pr_significance to significance */
  DOUBLE_FROM_R( pr_significance, &significance );

  /* ... pr_white_noise_test to white_noise_test */
  err = wav_white_test_from_R( pr_white_noise_test, &white_noise_test );
  CHECK_CONVERSION( wav_white_test, pr_white_noise_test, &white_noise_test );

  /* Call the function */
  err = wavuniv_transform_packet_whitest(
    &dwpt,
    significance,
    white_noise_test,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( wavuniv_transform_packet_whitest, &result, &pr_ret_result );
}

/** Adaptive wavelet-based bootstrapping.
 * @source RS\_wav\_boot.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_bootstrap", dwpt, filters, white.indices, n.realization))#
 * @return               An R ... containing ...
 * @param dwpt           Pointer to an R object containing ... dwpt
 * @param filters        Pointer to an R object containing ... filters
 * @param white.indices  Pointer to an R object containing ... white.indices
 * @param n.realization  Pointer to an R object containing ... n.realization
 * @see wavuniv_transform_packet_whitest
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_packet_inverse
*/
EXTERN_R SEXP RS_wavelets_bootstrap(
 SEXP pr_dwpt,
 SEXP pr_filters,
 SEXP pr_white_indices,
 SEXP pr_n_realization )
{
  SEXP             pr_ret_result;   
  mat_set          dwpt;            
  mat_set          filters;         
  mat_set          result;          
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           n_realization;   
  univ_mat         white_indices;   
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_dwpt to dwpt */
  READ_MATSET_REGISTER( pr_dwpt, MUTIL_DOUBLE, &dwpt );

  /* ... pr_filters to filters */
  READ_MATSET_REGISTER( pr_filters, MUTIL_DOUBLE, &filters );

  /* ... pr_white_indices to white_indices */
  READ_MATRIX_REGISTER( pr_white_indices, &white_indices );

  /* ... pr_n_realization to n_realization */
  SINT32_FROM_R( pr_n_realization, &n_realization );

  /* Call the function */
  err = wavuniv_bootstrap(
    &dwpt,
    &filters,
    &white_indices,
    n_realization,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( wavuniv_bootstrap, &result, &pr_ret_result );
}

