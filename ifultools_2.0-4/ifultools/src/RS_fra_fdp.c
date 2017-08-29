
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_fdp.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS wavelets library.

   Functions wrapped:

   frauniv_fdp_simulate()
   frauniv_fdp_simulate_weights()
*/

#include "fra_fdp.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Function to simulate an FD process with time varying model parameters.
 * @source RS\_fra\_fdp.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_fdp_simulate", delta, innovation.variance))#
 * @return                     An R ... containing ...
 * @param delta                Pointer to an R object containing ... delta
 * @param innovation.variance  Pointer to an R object containing ... innovation.variance
 * @see frauniv_fdp_simulate_weights
*/
EXTERN_R SEXP RS_wavelets_fdp_simulate(
 SEXP pr_delta,
 SEXP pr_innovation_variance )
{
  SEXP                  pr_ret_result;        
  mutil_data_type       type;                 
  mutil_errcode         err;                  
  univ_mat              delta;                
  univ_mat              innovation_variance;  
  univ_mat              result;               
  void                  *VPNULL = NULL;       
  memlist               list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_delta to delta */
  READ_MATRIX_REGISTER( pr_delta, &delta );

  /* ... pr_innovation_variance to innovation_variance */
  READ_MATRIX_REGISTER( pr_innovation_variance, &innovation_variance );

  /* Call the function */
  err = frauniv_fdp_simulate(
    &delta,
    &innovation_variance,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_fdp_simulate, &result, &pr_ret_result );
}

/** Function to generate the weights for a time varying FD process simulation.
 * @source RS\_fra\_fdp.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_fdp_simulate_weights", delta, innovation.variance))#
 * @return                     An R ... containing ...
 * @param delta                Pointer to an R object containing ... delta
 * @param innovation.variance  Pointer to an R object containing ... innovation.variance
 * @see frauniv_fdp_simulate
*/
EXTERN_R SEXP RS_wavelets_fdp_simulate_weights(
 SEXP pr_delta,
 SEXP pr_innovation_variance )
{
  SEXP                  pr_ret_result;        
  mutil_data_type       type;                 
  mutil_errcode         err;                  
  univ_mat              delta;                
  univ_mat              innovation_variance;  
  univ_mat              result;               
  void                  *VPNULL = NULL;       
  memlist               list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_delta to delta */
  READ_MATRIX_REGISTER( pr_delta, &delta );

  /* ... pr_innovation_variance to innovation_variance */
  READ_MATRIX_REGISTER( pr_innovation_variance, &innovation_variance );

  /* Call the function */

  err = frauniv_fdp_simulate_weights(
    &delta,
    &innovation_variance,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_fdp_simulate_weights, &result, &pr_ret_result );
}

