
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_modl.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS fractal library.

   Functions wrapped:

   frauniv_determinism_delta_epsilon()
*/

#include "fra_modl.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Detecting determinism in a time series.
 * @source RS\_fra\_modl.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_determinism_delta_epsilon", time.series, embedding.dimension, time.lag, orbital.lag, image.lag, scale.min, scale.max, resolution, minimize))#
 * @return                     An R ... containing ...
 * @param time.series          Pointer to an R object containing ... time.series
 * @param embedding.dimension  Pointer to an R object containing ... embedding.dimension
 * @param time.lag             Pointer to an R object containing ... time.lag
 * @param orbital.lag          Pointer to an R object containing ... orbital.lag
 * @param image.lag            Pointer to an R object containing ... image.lag
 * @param scale.min            Pointer to an R object containing ... scale.min
 * @param scale.max            Pointer to an R object containing ... scale.max
 * @param resolution           Pointer to an R object containing ... resolution
 * @param minimize             Pointer to an R object containing ... minimize
 * @see frauniv_surrogate
*/
EXTERN_R SEXP RS_fractal_determinism_delta_epsilon(
 SEXP pr_time_series,
 SEXP pr_embedding_dimension,
 SEXP pr_time_lag,
 SEXP pr_orbital_lag,
 SEXP pr_image_lag,
 SEXP pr_scale_min,
 SEXP pr_scale_max,
 SEXP pr_resolution,
 SEXP pr_minimize )
{
  SEXP                  pr_ret_result;        
  boolean               minimize;             
  double                resolution;           
  double                scale_max;            
  double                scale_min;            
  mat_set               result;               
  mutil_data_type       type;                 
  mutil_errcode         err;                  
  sint32                embedding_dimension;  
  sint32                image_lag;            
  sint32                orbital_lag;          
  sint32                time_lag;             
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

  /* ... pr_image_lag to image_lag */
  SINT32_FROM_R( pr_image_lag, &image_lag );

  /* ... pr_scale_min to scale_min */
  DOUBLE_FROM_R( pr_scale_min, &scale_min );

  /* ... pr_scale_max to scale_max */
  DOUBLE_FROM_R( pr_scale_max, &scale_max );

  /* ... pr_resolution to resolution */
  DOUBLE_FROM_R( pr_resolution, &resolution );

  /* ... pr_minimize to minimize */
  BOOLEAN_FROM_R( pr_minimize, &minimize );

  /* Call the function */
  err = frauniv_determinism_delta_epsilon(
    &time_series,
    embedding_dimension,
    time_lag,
    orbital_lag,
    image_lag,
    scale_min,
    scale_max,
    resolution,
    minimize,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( frauniv_determinism_delta_epsilon, &result, &pr_ret_result );
}

