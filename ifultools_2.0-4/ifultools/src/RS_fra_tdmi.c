
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_tdmi.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS fractal library.

   Functions wrapped:

   frauniv_time_delayed_mutual_information()
*/

#include "fra_tdmi.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Time delayed mutual information.
 * @source RS\_fra\_tdmi.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_time_delayed_mutual_information", time.series, lags))#
 * @return             An R ... containing ...
 * @param time.series  Pointer to an R object containing ... time.series
 * @param lags         Pointer to an R object containing ... lags
 * @see frauniv_kernel_density_estimate
*/
EXTERN_R SEXP RS_fractal_time_delayed_mutual_information(
 SEXP pr_time_series,
 SEXP pr_lags )
{
  SEXP             pr_ret_result;     
  mutil_data_type  type;            
  mutil_errcode    err;             
  univ_mat         lags;            
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

  /* ... pr_lags to lags */
  READ_MATRIX_REGISTER( pr_lags, &lags );

  /* Call the function */
  err = frauniv_time_delayed_mutual_information(
    &time_series,
    &lags,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_time_delayed_mutual_information, &result, &pr_ret_result );
}

