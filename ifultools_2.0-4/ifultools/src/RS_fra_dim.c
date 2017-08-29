
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_dim.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS fractal library.

   Functions wrapped:

   frauniv_embed()
   frauniv_dimension_information()
   frauniv_dimension_correlation_summation()
   frauniv_dimension_false_nearest_neighbors()
   frauniv_dimension_false_nearest_strands()
   frauniv_poincare_map()
   frauniv_space_time_separation_plot()
*/

#include "fra_dim.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Delay embedding of a univariate time series.
 * @source RS\_fra\_dim.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_embed", time.series, embedding.dimension, time.lag))#
 * @return                     An R ... containing ...
 * @param time.series          Pointer to an R object containing ... time.series
 * @param embedding.dimension  Pointer to an R object containing ... embedding.dimension
 * @param time.lag             Pointer to an R object containing ... time.lag
 * @see frauniv_dimension_correlation_summation
 * @see frauniv_dimension_information
*/
EXTERN_R SEXP RS_fractal_embed(
 SEXP pr_time_series,
 SEXP pr_embedding_dimension,
 SEXP pr_time_lag )
{
  SEXP                  pr_ret_result;
  mutil_data_type       type;
  mutil_errcode         err;
  sint32                embedding_dimension;
  sint32                time_lag;
  univ_mat              result;
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

  /* Call the function */
  err = frauniv_embed(
    &time_series,
    embedding_dimension,
    time_lag,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_embed, &result, &pr_ret_result );
}

/** Information dimension statistics.
 * @source RS\_fra\_dim.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_dimension_information", time.series, embedding.dimension, time.lag, orbital.lag, n.density, distance.metric, max.neighbors, n.reference))#
 * @return                     An R ... containing ...
 * @param time.series          Pointer to an R object containing ... time.series
 * @param embedding.dimension  Pointer to an R object containing ... embedding.dimension
 * @param time.lag             Pointer to an R object containing ... time.lag
 * @param orbital.lag          Pointer to an R object containing ... orbital.lag
 * @param n.density            Pointer to an R object containing ... n.density
 * @param distance.metric      Pointer to an R object containing ... distance.metric
 * @param max.neighbors        Pointer to an R object containing ... max.neighbors
 * @param n.reference          Pointer to an R object containing ... n.reference
 * @see frauniv_embed
 * @see frauniv_dimension_correlation_summation
*/
EXTERN_R SEXP RS_fractal_dimension_information(
 SEXP pr_time_series,
 SEXP pr_embedding_dimension,
 SEXP pr_time_lag,
 SEXP pr_orbital_lag,
 SEXP pr_n_density,
 SEXP pr_distance_metric,
 SEXP pr_max_neighbors,
 SEXP pr_n_reference )
{
  SEXP                  pr_ret_result;
  fra_distance_metric   distance_metric;
  mat_set               result;
  mutil_data_type       type;
  mutil_errcode         err;
  sint32                embedding_dimension;
  sint32                max_neighbors;
  sint32                n_density;
  sint32                n_reference;
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

  /* ... pr_n_density to n_density */
  SINT32_FROM_R( pr_n_density, &n_density );

  /* ... pr_distance_metric to distance_metric */
  DISTANCE_METRIC_FROM_R( pr_distance_metric, &distance_metric );

  /* ... pr_max_neighbors to max_neighbors */
  SINT32_FROM_R( pr_max_neighbors, &max_neighbors );

  /* ... pr_n_reference to n_reference */
  SINT32_FROM_R( pr_n_reference, &n_reference );

  /* Call the function */
  err = frauniv_dimension_information(
    &time_series,
    embedding_dimension,
    time_lag,
    orbital_lag,
    n_density,
    distance_metric,
    max_neighbors,
    n_reference,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( frauniv_dimension_information, &result, &pr_ret_result );
}

/** Discrete correlation integral estimation using a delay embedding of a single variable time series.
 * @source RS\_fra\_dim.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_dimension_correlation_summation", time.series, emb.dim, time.lag, orbital.lag, n.scale.octave, distance.metric, max.neighbors, n.reference))#
 * @return                An R ... containing ...
 * @param time.series     Pointer to an R object containing ... time.series
 * @param emb.dim         Pointer to an R object containing ... emb.dim
 * @param time.lag        Pointer to an R object containing ... time.lag
 * @param orbital.lag     Pointer to an R object containing ... orbital.lag
 * @param n.scale.octave  Pointer to an R object containing ... n.scale.octave
 * @see frauniv_embed
 * @see frauniv_dimension_information
*/
EXTERN_R SEXP RS_fractal_dimension_correlation_summation(
 SEXP pr_time_series,
 SEXP pr_embedding_dimension,
 SEXP pr_time_lag,
 SEXP pr_orbital_lag,
 SEXP pr_n_scale_octave )
{
  SEXP             pr_ret_result;
  mat_set          result;
  mutil_data_type  type;
  mutil_errcode    err;
  sint32           embedding_dimension;
  sint32           n_scale_octave;
  sint32           orbital_lag;
  sint32           time_lag;
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

  /* ... pr_embedding_dimension to embedding_dimension */
  SINT32_FROM_R( pr_embedding_dimension, &embedding_dimension );

  /* ... pr_time_lag to time_lag */
  SINT32_FROM_R( pr_time_lag, &time_lag );

  /* ... pr_orbital_lag to orbital_lag */
  SINT32_FROM_R( pr_orbital_lag, &orbital_lag );

  /* ... pr_n_scale_octave to n_scale_octave */
  SINT32_FROM_R( pr_n_scale_octave, &n_scale_octave );

  /* Call the function */
  err = frauniv_dimension_correlation_summation(
    &time_series,
    embedding_dimension,
    time_lag,
    orbital_lag,
    n_scale_octave,
    VPNULL,
    &result );
  CONVERT_MATSET_AND_RETURN( frauniv_dimension_correlation_summation, &result, &pr_ret_result );
}

/** Estimation of the proper embedding dimension for a single-variable time series.
 * @source RS\_fra\_dim.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_dimension_false_nearest_neighbors", time.series, embedding.dimension, time.lag, orbital.lag, rtol, atol, max.neighbors, n.reference))#
 * @return                     An R ... containing ...
 * @param time.series          Pointer to an R object containing ... time.series
 * @param embedding.dimension  Pointer to an R object containing ... embedding.dimension
 * @param time.lag             Pointer to an R object containing ... time.lag
 * @param orbital.lag          Pointer to an R object containing ... orbital.lag
 * @param rtol                 Pointer to an R object containing ... rtol
 * @param atol                 Pointer to an R object containing ... atol
 * @see frauniv_embed
 * @see frauniv_dimension_information
 * @see frauniv_dimension_correlation_summation
*/
EXTERN_R SEXP RS_fractal_dimension_false_nearest_neighbors(
 SEXP pr_time_series,
 SEXP pr_embedding_dimension,
 SEXP pr_time_lag,
 SEXP pr_orbital_lag,
 SEXP pr_rtol,
 SEXP pr_atol )
{
  SEXP                  pr_ret_result;
  double                atol;
  double                rtol;
  mutil_data_type       type;
  mutil_errcode         err;
  sint32                embedding_dimension;
  sint32                orbital_lag;
  sint32                time_lag;
  univ_mat              result;
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

  /* ... pr_rtol to rtol */
  DOUBLE_FROM_R( pr_rtol, &rtol );

  /* ... pr_atol to atol */
  DOUBLE_FROM_R( pr_atol, &atol );

  /* Call the function */
  err = frauniv_dimension_false_nearest_neighbors(
    &time_series,
    embedding_dimension,
    time_lag,
    orbital_lag,
    rtol,
    atol,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_dimension_false_nearest_neighbors, &result, &pr_ret_result );
}

/** Estimation of the proper embedding dimension for a single-variable time series.
 * @source RS\_fra\_dim.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_dimension_false_nearest_strands", time.series, embedding.dimension, time.lag, orbital.lag, iterate.tolerance, atol, max.neighbors, n.reference))#
 * @return                     An R ... containing ...
 * @param time.series          Pointer to an R object containing ... time.series
 * @param embedding.dimension  Pointer to an R object containing ... embedding.dimension
 * @param time.lag             Pointer to an R object containing ... time.lag
 * @param orbital.lag          Pointer to an R object containing ... orbital.lag
 * @param iterate.tolerance    Pointer to an R object containing ... iterate.tolerance
 * @param atol                 Pointer to an R object containing ... atol
 * @see frauniv_embed
 * @see frauniv_false_nearest_neighbors
 * @see frauniv_dimension_information
 * @see frauniv_dimension_correlation_summation
*/
EXTERN_R SEXP RS_fractal_dimension_false_nearest_strands(
 SEXP pr_time_series,
 SEXP pr_embedding_dimension,
 SEXP pr_time_lag,
 SEXP pr_orbital_lag,
 SEXP pr_iterate_tolerance,
 SEXP pr_atol )
{
  SEXP                  pr_ret_result;
  double                atol;
  mutil_data_type       type;
  mutil_errcode         err;
  sint32                embedding_dimension;
  sint32                iterate_tolerance;
  sint32                orbital_lag;
  sint32                time_lag;
  univ_mat              result;
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

  /* ... pr_iterate_tolerance to iterate_tolerance */
  SINT32_FROM_R( pr_iterate_tolerance, &iterate_tolerance );

  /* ... pr_atol to atol */
  DOUBLE_FROM_R( pr_atol, &atol );

  /* Call the function */
  err = frauniv_dimension_false_nearest_strands(
    &time_series,
    embedding_dimension,
    time_lag,
    orbital_lag,
    iterate_tolerance,
    atol,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_dimension_false_nearest_strands, &result, &pr_ret_result );
}

/** Poincare map.
 * @source RS\_fra\_dim.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_poincare_map", time.series, extrema.type, denoise, orbital.lag, iterate.tolerance, atol, max.neighbors, n.reference))#
 * @return              An R ... containing ...
 * @param time.series   Pointer to an R object containing ... time.series
 * @param extrema.type  Pointer to an R object containing ... extrema.type
 * @param denoise       Pointer to an R object containing ... denoise
 * @see _fra_extrema_type
 * @see wavuniv_shrink
*/
EXTERN_R SEXP RS_fractal_poincare_map(
 SEXP pr_time_series,
 SEXP pr_extrema_type,
 SEXP pr_denoise )
{
  SEXP                  pr_ret_extrema_amplitude;
  SEXP                  pr_ret_extrema_location;
  SEXP                  pr_ret_obj;
  boolean               denoise;
  fra_extrema_type      extrema_type;
  mutil_data_type       type;
  mutil_errcode         err;
  univ_mat              extrema_amplitude;
  univ_mat              extrema_location;
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

  /* ... pr_extrema_type to extrema_type */
  err = fra_extrema_type_from_R( pr_extrema_type, &extrema_type );
  CHECK_CONVERSION( fra_extrema_type, pr_extrema_type, &extrema_type );

  /* ... pr_denoise to denoise */
  BOOLEAN_FROM_R( pr_denoise, &denoise );

  /* Call the function */
  err = frauniv_poincare_map(
    &time_series,
    extrema_type,
    denoise,
    VPNULL,
    &extrema_location,
    &extrema_amplitude );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Problem calling frauniv_poincare_map() function" );
  err = memlist_member_register( &list, &extrema_location, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
  err = memlist_member_register( &list, &extrema_amplitude, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );

  /* create the output R object */
  err = matuniv_to_R( &extrema_location, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_extrema_location );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  err = matuniv_to_R( &extrema_amplitude, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_extrema_amplitude );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  PROTECT( pr_ret_obj = allocVector( VECSXP, 2 ) );
  SET_VECTOR_ELT( pr_ret_obj, 0, pr_ret_extrema_location );
  SET_VECTOR_ELT( pr_ret_obj, 1, pr_ret_extrema_amplitude );
  UNPROTECT(1);

  /* free registered local memory */
  MUTIL_FREE_WARN( memlist, &list );

  return pr_ret_obj;
}

/** Space time separation plot.
 * @source RS\_fra\_dim.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_space_time_separation_plot", time.series, dim, delay, orbital.lags, frac.prob, atol, max.neighbors, n.reference))#
 * @return              An R ... containing ...
 * @param time.series   Pointer to an R object containing ... time.series
 * @param dim           Pointer to an R object containing ... dim
 * @param delay         Pointer to an R object containing ... delay
 * @param orbital.lags  Pointer to an R object containing ... orbital.lags
 * @param frac.prob     Pointer to an R object containing ... frac.prob
 * @see frauniv_embed
*/
EXTERN_R SEXP RS_fractal_space_time_separation_plot(
 SEXP pr_time_series,
 SEXP pr_embedding_dimension,
 SEXP pr_time_lag,
 SEXP pr_orbital_lags,
 SEXP pr_frac_prob )
{
  SEXP             pr_ret_max_eps;
  SEXP             pr_ret_obj;
  SEXP             pr_ret_result;
  double           frac_prob;
  double           max_eps;
  mutil_data_type  type;
  mutil_errcode    err;
  sint32           time_lag;
  sint32           embedding_dimension;
  univ_mat         orbital_lags;
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

  /* ... pr_embedding_dimension to embedding_dimension */
  SINT32_FROM_R( pr_embedding_dimension, &embedding_dimension );

  /* ... pr_time_lag to time_lag */
  SINT32_FROM_R( pr_time_lag, &time_lag );

  /* ... pr_orbital_lag to orbital_lag */
  READ_MATRIX_REGISTER( pr_orbital_lags, &orbital_lags );

  /* ... pr_frac_prob to frac_prob */
  DOUBLE_FROM_R( pr_frac_prob, &frac_prob );

  /* Call the function */
  err = frauniv_space_time_separation_plot(
    &time_series,
    embedding_dimension,
    time_lag,
    &orbital_lags,
    frac_prob,
    VPNULL,
    &max_eps,
    &result );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Problem calling frauniv_space_time_separation_plot() function" );
  err = memlist_member_register( &list, &result, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );

  /* create the output R object */
  err = double_to_R( max_eps, &pr_ret_max_eps );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  err = matuniv_to_R( &result, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_result );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  PROTECT( pr_ret_obj = allocVector( VECSXP, 2 ) );
  SET_VECTOR_ELT( pr_ret_obj, 0, pr_ret_max_eps );
  SET_VECTOR_ELT( pr_ret_obj, 1, pr_ret_result );
  UNPROTECT(1);

  /* free registered local memory */
  MUTIL_FREE_WARN( memlist, &list );

  return pr_ret_obj;
}

