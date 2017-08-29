
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_filt.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS fractal library.

   Functions wrapped:

   frauniv_filter_median()
   frauniv_filter_nonlinear_local_projection()
*/

#include "fra_filt.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Median filtering of a time series.
 * @source RS\_fra\_filt.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_filter_median", time.series, order))#
 * @return             An R ... containing ...
 * @param time.series  Pointer to an R object containing ... time.series
 * @param order        Pointer to an R object containing ... order
 * @see frauniv_filter_local_projection
*/
EXTERN_R SEXP RS_fractal_filter_median(
 SEXP pr_time_series,
 SEXP pr_order )
{
  SEXP             pr_ret_result;   
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           order;           
  univ_mat         result;          
  univ_mat         time_series;     
  void             *VPNULL = NULL;  
  memlist               list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_order to order */
  SINT32_FROM_R( pr_order, &order );

  /* Call the function */
  err = frauniv_filter_median(
    &time_series,
    order,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_filter_median, &result, &pr_ret_result );
}

/** Nonlinear denoising via local projection.
 * @source RS\_fra\_filt.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_filter_nonlinear_local_projection", time.series, dim, delay, min.nneig, max.radius, distance.metric, noise.dim, correct.curvature))#
 * @return                   An R ... containing ...
 * @param time.series        Pointer to an R object containing ... time.series
 * @param dim                Pointer to an R object containing ... dim
 * @param delay              Pointer to an R object containing ... delay
 * @param min.nneig          Pointer to an R object containing ... min.nneig
 * @param max.radius         Pointer to an R object containing ... max.radius
 * @param distance.metric    Pointer to an R object containing ... distance.metric
 * @param noise.dim          Pointer to an R object containing ... noise.dim
 * @param correct.curvature  Pointer to an R object containing ... correct.curvature
 * @see _fra_distance_metric
*/
EXTERN_R SEXP RS_fractal_filter_nonlinear_local_projection(
 SEXP pr_time_series,
 SEXP pr_dim,
 SEXP pr_delay,
 SEXP pr_min_nneig,
 SEXP pr_max_radius,
 SEXP pr_distance_metric,
 SEXP pr_noise_dim,
 SEXP pr_correct_curvature )
{
  SEXP                 pr_ret_result;       
  boolean              correct_curvature;   
  double               max_radius;          
  fra_distance_metric  distance_metric;     
  mutil_data_type      type;                
  mutil_errcode        err;                 
  sint32               delay;               
  sint32               dim;                 
  sint32               min_nneig;           
  sint32               noise_dim;           
  univ_mat             result;              
  univ_mat             time_series;         
  void                 *VPNULL = NULL;      
  memlist               list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_dim to dim */
  SINT32_FROM_R( pr_dim, &dim );

  /* ... pr_delay to delay */
  SINT32_FROM_R( pr_delay, &delay );

  /* ... pr_min_nneig to min_nneig */
  SINT32_FROM_R( pr_min_nneig, &min_nneig );

  /* ... pr_max_radius to max_radius */
  DOUBLE_FROM_R( pr_max_radius, &max_radius );

  /* ... pr_distance_metric to distance_metric */
  DISTANCE_METRIC_FROM_R( pr_distance_metric, &distance_metric );

  /* ... pr_noise_dim to noise_dim */
  SINT32_FROM_R( pr_noise_dim, &noise_dim );

  /* ... pr_correct_curvature to correct_curvature */
  BOOLEAN_FROM_R( pr_correct_curvature, &correct_curvature );

  /* Call the function */
  err = frauniv_filter_nonlinear_local_projection(
    &time_series,
    dim,
    delay,
    min_nneig,
    max_radius,
    distance_metric,
    noise_dim,
    correct_curvature,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_filter_nonlinear_local_projection, &result, &pr_ret_result );
}

