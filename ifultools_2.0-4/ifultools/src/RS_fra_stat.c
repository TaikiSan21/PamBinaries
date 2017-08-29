
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_stat.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS fractal library.

   Functions wrapped:

   frauniv_stationarity_priestley_subba_rao()
*/

#include "fra_stat.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Stationarity tests for a time series.
 * @source RS\_fra\_stat.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_stationarity_priestley_subba_rao", time.series, sampling.interval, n.taper, n.block, significance, center, recenter))#
 * @return                   An R ... containing ...
 * @param time.series        Pointer to an R object containing ... time.series
 * @param sampling.interval  Pointer to an R object containing ... sampling.interval
 * @param n.taper            Pointer to an R object containing ... n.taper
 * @param n.block            Pointer to an R object containing ... n.block
 * @param significance       Pointer to an R object containing ... significance
 * @param center             Pointer to an R object containing ... center
 * @param recenter           Pointer to an R object containing ... recenter
 * @see frauniv_determinism_delta_epsilon
*/
EXTERN_R SEXP RS_fractal_stationarity_priestley_subba_rao(
 SEXP pr_time_series,
 SEXP pr_sampling_interval,
 SEXP pr_n_taper,
 SEXP pr_n_block,
 SEXP pr_significance,
 SEXP pr_center,
 SEXP pr_recenter )
{
  SEXP             pr_ret_result;      
  boolean          center;             
  boolean          recenter;           
  double           sampling_interval;  
  double           significance;       
  mat_set          result;             
  mutil_data_type  type;               
  mutil_errcode    err;                
  sint32           n_block;            
  sint32           n_taper;            
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

  /* ... pr_sampling_interval to sampling_interval */
  DOUBLE_FROM_R( pr_sampling_interval, &sampling_interval );

  /* ... pr_n_taper to n_taper */
  SINT32_FROM_R( pr_n_taper, &n_taper );

  /* ... pr_n_block to n_block */
  SINT32_FROM_R( pr_n_block, &n_block );

  /* ... pr_significance to significance */
  DOUBLE_FROM_R( pr_significance, &significance );

  /* ... pr_center to center */
  BOOLEAN_FROM_R( pr_center, &center );

  /* ... pr_recenter to recenter */
  BOOLEAN_FROM_R( pr_recenter, &recenter );

  /* Call the function */
  err = frauniv_stationarity_priestley_subba_rao(
    &time_series,
    sampling_interval,
    n_taper,
    n_block,
    significance,
    center,
    recenter,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( frauniv_stationarity_priestley_subba_rao, &result, &pr_ret_result );
}

