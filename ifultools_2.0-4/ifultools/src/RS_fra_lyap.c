
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_lyap.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS fractal library.

   Functions wrapped:

   frauniv_local_lyapunov_spectrum()
*/

#include "fra_lyap.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Local Lyapunov spectrum estimation.
 * @source RS\_fra\_lyap.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_local_lyapunov_spectrum", time.series, embedding.dimension, time.lag, orbital.lag, sampling.interval, local.dimension, polynomial.order, global.reference, n.reference.local, metric, scale))#
 * @return                     An R ... containing ...
 * @param time.series          Pointer to an R object containing ... time.series
 * @param embedding.dimension  Pointer to an R object containing ... embedding.dimension
 * @param time.lag             Pointer to an R object containing ... time.lag
 * @param orbital.lag          Pointer to an R object containing ... orbital.lag
 * @param sampling.interval    Pointer to an R object containing ... sampling.interval
 * @param local.dimension      Pointer to an R object containing ... local.dimension
 * @param polynomial.order     Pointer to an R object containing ... polynomial.order
 * @param global.reference     Pointer to an R object containing ... global.reference
 * @param n.reference.local    Pointer to an R object containing ... n.reference.local
 * @param metric               Pointer to an R object containing ... metric
 * @param scale                Pointer to an R object containing ... scale
 * @see frauniv_embed
 * @see frauniv_dimension_correlation_summation
 * @see frauniv_dimension_information
*/
EXTERN_R SEXP RS_fractal_local_lyapunov_spectrum(
 SEXP pr_time_series,
 SEXP pr_embedding_dimension,
 SEXP pr_time_lag,
 SEXP pr_orbital_lag,
 SEXP pr_sampling_interval,
 SEXP pr_local_dimension,
 SEXP pr_polynomial_order,
 SEXP pr_global_reference,
 SEXP pr_n_reference_local,
 SEXP pr_metric,
 SEXP pr_scale )
{
  SEXP                  pr_ret_result;        
  double                sampling_interval;    
  fra_distance_metric   metric;               
  mat_set               result;               
  mutil_data_type       type;                 
  mutil_errcode         err;                  
  sint32                embedding_dimension;  
  sint32                local_dimension;      
  sint32                n_reference_local;    
  sint32                orbital_lag;          
  sint32                polynomial_order;     
  sint32                time_lag;             
  univ_mat              global_reference;     
  univ_mat              scale;                
  univ_mat              time_series;          
  void                  *VPNULL = NULL;       
  memlist               list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_embedding_dimension to embedding_dimension */
  SINT32_FROM_R( pr_embedding_dimension, &embedding_dimension );

  /* ... pr_time_lag to time_lag */
  SINT32_FROM_R( pr_time_lag, &time_lag );

  /* ... pr_orbital_lag to orbital_lag */
  SINT32_FROM_R( pr_orbital_lag, &orbital_lag );

  /* ... pr_sampling_interval to sampling_interval */
  DOUBLE_FROM_R( pr_sampling_interval, &sampling_interval );

  /* ... pr_local_dimension to local_dimension */
  SINT32_FROM_R( pr_local_dimension, &local_dimension );

  /* ... pr_polynomial_order to polynomial_order */
  SINT32_FROM_R( pr_polynomial_order, &polynomial_order );

  /* ... pr_global_reference to global_reference */
  READ_MATRIX_REGISTER( pr_global_reference, &global_reference );

  /* ... pr_n_reference_local to n_reference_local */
  SINT32_FROM_R( pr_n_reference_local, &n_reference_local );

  /* ... pr_metric to metric */
  DISTANCE_METRIC_FROM_R( pr_metric, &metric );

  /* ... pr_scale to scale */
  READ_MATRIX_REGISTER( pr_scale, &scale );

  /* Call the function */
  err = frauniv_local_lyapunov_spectrum(
    &time_series,
    embedding_dimension,
    time_lag,
    orbital_lag,
    sampling_interval,
    local_dimension,
    polynomial_order,
    &global_reference,
    n_reference_local,
    metric,
    &scale,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( frauniv_local_lyapunov_spectrum, &result, &pr_ret_result );
}

