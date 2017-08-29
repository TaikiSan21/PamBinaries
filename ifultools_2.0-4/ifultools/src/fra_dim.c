
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_dim.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */

#include "fra_dim.h"
#include "fra_neig.h"
#include "fra_type.h"
#include "fra_util.h"

#include "mat_assn.h"
#include "mat_cast.h"
#include "mat_set.h"
#include "mat_rand.h"
#include "mat_stat.h"
#include "mat_type.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "str_type.h"

#include "ut_debug.h"
#include "ut_intrp.h"
#include "ut_limit.h"
#include "ut_math.h"
#include "ut_mem.h"
#include "ut_intrn.h"
#include "ut_rand.h"

#include "wav_cwt.h"
#include "wav_filt.h"
#include "wav_shrk.h"

#include <math.h>
#include <stdio.h>


/*
  This file contains function definitions for
  calculating phase space fractal dimensions.
  The functions are declared in fra_dim.h
*/

/* Static macro definitions */

/* Minumum time series input length */
#define FRA_POIN_MIN_N_SAMPLE        4

/* WaveShrink defaults */
#define FRA_POIN_WAVSHRK_DEF_FILT     WAV_FILTER_LEAST_ASYMMETRIC
#define FRA_POIN_WAVSHRK_DEF_FILT_LEN 8
#define FRA_POIN_WAVSHRK_DEF_FUNC     WAV_SHRINK_FUNCTION_MID
#define FRA_POIN_WAVSHRK_DEF_TSH      WAV_SHRINK_THRESHOLD_UNIVERSAL

/* Continuous wavelet transform function (fixed) arguments */
#define FRA_POIN_CWT_FILT_ARG         1.0
#define FRA_POIN_CWT_SCALE            0.25

#define FRAUNIV_STSP_MIN_FRAC_PROB   0.0001
#define FRAUNIV_STSP_NUM_EPSILON     1000

#define IS_EVEN( n ) (((n) % 2) == 0)

#undef LOCALDEF_CHECK_NULL_POINTER_DIM
#define LOCALDEF_CHECK_NULL_POINTER_DIM( DATA_PTR,          \
  DATA_TYPE, TYPE_PREFIX )                                  \
  err = TYPE_PREFIX ## _validate( DATA_PTR );               \
  if ( err ) return err;                                    \
  if ( DATA_PTR == ( DATA_TYPE * ) NULL ) {                 \
  MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
  return MUTIL_ERR_NULL_POINTER;                            \
  }

/* Static functions declared here and defined at end of file */

static void localdef_ipower( float *x, sint32 n );

/* static uint32 localdef_exponent_biased( float num ); */

static mutil_errcode localfn_average_neighborhood_radius_constant_density(
  univ_mat           *embedding,
  const mutil_kdtree *kdtree,
  const sint32         orbital_lag,
  const sint32         n_neighbor,
  const fra_distance_metric distance_metric,
  const sint32         n_average,
  void                *intrp_ptr,
  double              *average_log2_radius );

static mutil_errcode localfn_create_waveshrink_defaults(
  const sint32            n_sample,
  void                   *intrp_ptr,
  mat_set                *filters,
  sint32                 *n_level,
  wav_shrink_function    *ws_fun,
  wav_shrink_threshold   *ws_tsh );


static mutil_errcode localfn_extrema_locations(
  const univ_mat          *time_series,
  const fra_extrema_type   extrema_type,
  void                    *intrp_ptr,
  univ_mat                *extrema,
  sint8                   *concavity );

static mutil_errcode localfn_extrema_amplitudes(
  const univ_mat          *time_series,
  const univ_mat          *ex_locs,
  const fra_extrema_type   ex_type,
  const sint8             *concavity,
  univ_mat                *ex_vals );


static mutil_errcode localfn_stsp_error_checks(
  const univ_mat     *time_series,
  const sint32        delay,
  const sint32        dim,
  const univ_mat     *orbital_lags,
  const double        frac_prob,
  sint32             *ts_length,
  sint32             *nfracs,
  sint32             *nembed,
  sint32             *nolags );

static mutil_errcode localfn_kdtree(
  const sint32    n_neighbor,
  const sint32    max_neighbors,
  const univ_mat *embedding,
  const sint32    orbital_lag,
  const sint32    bucket_size,
  void           *intrp_ptr,
  mutil_kdtree   *kdtree,
  memlist        *list );


/* static boolean printit=FALSE; */

/* Delay embedding of a univariate time series   */
/*                                               */
/* Documented in fra_dim.h                       */
/* Written by William Constantine                */

mutil_errcode frauniv_embed(
 const univ_mat *time_series,
 sint32          embedding_dimension,
 const sint32    time_delay,
 void           *intrp_ptr,
 univ_mat       *result )
{
  sint32         k;
  sint32         i;
  sint32         d;
  sint32         n_embed;
  sint32         n_sample;
  sint32         max_dim;
  mutil_errcode  err;
  double        *pd_base;
  memlist        list;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_embed()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /*** check input data ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_DIM( time_series, univ_mat, matuniv );

  /* ... if number of elements is positive */

  if ( MATUNIV_NELEM( time_series ) < 1 ){
    MUTIL_ERROR( "Number of elements in input time series matrix must be positive." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if the type is double */

  if ( time_series->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Input time series matrix must be of type double." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( &( time_series->mat.dblmat ) ) ){
    MUTIL_ERROR( "Input time series must be a single-row or single-column matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /*** check embedding dimension ... ***/

  if ( embedding_dimension <= 0 ){
    MUTIL_ERROR( "Embedding dimension must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  /*** check time_delay ... ***/


  if ( time_delay < 1 ){
    MUTIL_ERROR( "Time delay must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  if ( time_delay >= MATUNIV_NELEM( time_series ) ){
    MUTIL_ERROR( "Time delay is too large. It must be less than the number of samples in time series." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  /* initialize dimension variables */

  n_sample = MATUNIV_NELEM( time_series );
  max_dim  = (sint32) floor( (double) ( n_sample - 1 ) / (double) time_delay ) + 1;

  if ( embedding_dimension > max_dim ){
    embedding_dimension = max_dim;
  }

  n_embed = n_sample - ( embedding_dimension - 1 ) * time_delay;

  /* allocate space for embedding */

  err = matuniv_malloc_register( result, n_embed, embedding_dimension, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create a delay embedding of the time series */

  pd_base = time_series->mat.dblmat.data;

  for ( k = 0, i = 0; i < n_embed; i++ ){

    for ( d = 0; d < embedding_dimension; d++ ){
      result->mat.dblmat.data[ k++ ] = pd_base[ d * time_delay ];
    }

    pd_base++;

    if ( MUTIL_INTERRUPT( 3.0 * n_embed, intrp_ptr ) ) {
         MUTIL_ERROR( "user interrupt" );
         MUTIL_FREE_WARN( memlist, &list );
         return MUTIL_ERR_INTERRUPT;
    }
  }

  /* free nodes corresponding to registered
  memory for the result, but do not free
  the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
  corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with frauniv_embed()" );

  return MUTIL_ERR_OK;
}

/* Information dimension estimation of a delay */
/* coordinate embedding                        */
/*                                             */
/* Documented in fra_dim.h                     */
/* Written by William Constantine              */

mutil_errcode frauniv_dimension_information(
  const univ_mat *time_series,
  const sint32    embedding_dimension,
  const sint32    time_delay,
  const sint32    orbital_lag,
  const sint32    n_density,
  const fra_distance_metric distance_metric,
  const sint32    max_neighbors,
  const sint32    n_reference,
  void           *intrp_ptr,
  mat_set        *result )
{
  boolean        is_delay_embedding;
  double         exponent;
  double         logmin;
  mutil_kdtree   kdtree;
  memlist        list;
  mutil_errcode  err;
  sint32         bucket_size = 1;
  sint32         dim;
  sint32         i;
  sint32         n_neighbor;
  sint32         n_embed;
  univ_mat       embedding;
  sint32         n_return = 2;
  double         density;
  sint32_mat     nrow;
  sint32_mat     ncol;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_dimension_information()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  err = frauniv_check_embedding_inputs(
    time_series, embedding_dimension, time_delay, orbital_lag,
    intrp_ptr, &is_delay_embedding, &n_embed );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( MATUNIV_NELEM( time_series ) < 1000 ){
    MUTIL_ERROR( "Time series must be at least 1000 points long" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( n_density < 1 ){
    MUTIL_ERROR( "Number of densities must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  switch( distance_metric ){
  case FRA_DISTANCE_L1:
  case FRA_DISTANCE_L2:
  case FRA_DISTANCE_LINFINITY:
    break;
  default:
    MUTIL_ERROR( "Unsupported distance metric" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( max_neighbors < 100 ){
    MUTIL_ERROR( "Maximum neighbors must be at least 100" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( max_neighbors > n_embed ){
    MUTIL_ERROR( "Maximum neighbors greater than the number of points in the embedding" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( n_reference < 10 ){
    MUTIL_ERROR( "Maximum neighbors must be at least 10" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( n_reference > n_embed ){
    MUTIL_ERROR( "Maximum neighbors greater than the number of points in the embedding" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( MUTIL_INTERRUPT( 3.0 * n_embed, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }

  /* allocate memory */

  err = mats32_malloc_register( &nrow, 1, 2, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &ncol, 1, 2, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_register( result, 1, &n_return, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* ... the D1 matrix */

  nrow.data[0] = n_density;
  ncol.data[0] = embedding_dimension - 1;

  /* ... the density vector */

  nrow.data[1] = n_density;
  ncol.data[1] = 1;

  /* ... allocate space for the matrix set and matrices */

  err = matset_malloc_matrices_arbitrary_size(
    result, &nrow, &ncol, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* form the log2 of the density vector (the values are linear
  in log2 space) */

  logmin = MUTIL_LOG2( 1.0 / (double) n_density );

  for ( i = 0; i < n_density; i++ ){

    exponent = logmin * ( 1.0 - (double) i / ( (double) n_density - 1.0 ) );

    result->mats[ 1 ].mat.dblmat.data[ i ] = exponent;
  }


  /* form the information dimension statistics for each embedding dimension  */

  for ( dim = 2; dim <= embedding_dimension; dim++ ){

    /* create embedding and register with the memory manager */

    err = frauniv_embed( time_series, dim, time_delay, intrp_ptr, &embedding );
    MEMLIST_FREE_ON_ERROR( err, &list );
    err = memlist_member_register( &list, &embedding, MEMTYPE_MATUNIV );
    if ( err ){
      MUTIL_FREE_WARN( matuniv, &embedding );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }


    /* loop through each density and calculate the average neighborhood radius */

    for ( i = 0; i < n_density; i++ ){

      density = MUTIL_POW( 2.0, result->mats[ 1 ].mat.dblmat.data[ i ] );

      /* as density is on (0,1], the variable n_neighbor is the
      associated number of neighbors we need to find
      for for each reference point in the analysis */

      n_neighbor = (sint32) ( density * (double) n_embed );

      /* form the kd-tree for fast nearest neighbor searching
      and register the associated memory with the memory manager
      (this latter part is taken care of in the called local function). */

      err = localfn_kdtree(
        n_neighbor,
        max_neighbors,
        &embedding,
        orbital_lag,
        bucket_size,
        intrp_ptr,
        &kdtree,
        &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* find the average radius of neighborhoods with a constant density */

      err = localfn_average_neighborhood_radius_constant_density(
        &embedding,
        &kdtree,
        orbital_lag,
        MUTIL_MIN( n_neighbor, max_neighbors ),
        distance_metric,
        n_reference,
        intrp_ptr,
        &( result->mats[ 0 ].mat.dblmat.data[ i * ( embedding_dimension - 1 ) + dim - 2 ] ) );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* free the kdtree */

      err = memlist_member_free( &kdtree, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

    } /* end loop for each density */


    /* free the memory for the current embedding */
    err = memlist_member_free( &embedding, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

  } /* end loop for each dimension */

  /* free nodes corresponding to registered
  memory for the result, but do not free
  the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
  corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with frauniv_dimension_information()" );

  return MUTIL_ERR_OK;
}

/* Discrete approximation to the correlation     */
/* integrals defined over a finite set of points */
/* embedded in the phase space via a single      */
/* variable time delay embedding scheme.         */
/*                                               */
/* Documented in fra_dim.h                       */
/* Written by William Constantine                */

mutil_errcode frauniv_dimension_correlation_summation(
 const univ_mat *data,
 sint32          dimension,
 sint32          time_lag,
 sint32          orbital_lag,
 sint32          resolution,
 void           *intrp_ptr,
 mat_set        *result )
{
  boolean        is_delay_embedding;
  double         cum_sum;
  double         dfac;
  double         norminv;
  double        *pd_C;
  double        *pd_corr;
  double        *pd_scale;
  double        *pd_data;
  double        *pd_x;
  double        *pd_y;
  float          distance;
  float          infinity_norm;
  memlist        list;
  mutil_errcode  err;
  sint32         E;
  sint32         dims;
  sint32         i;
  sint32         iscale;
  sint32         iscale_max;
  sint32         iscale_min;
  sint32         j;
  sint32         k;
  sint32         max_lag;
  sint32         exp_offset = 1022;
  sint32         n_embed;
  sint32         n_sample;
  sint32         n_scale;
  sint32_mat     ijscale_min;
  sint32_mat     ncol;
  sint32_mat     nrow;
  univ_mat       C2;
  univ_mat       scale;
  boolean        is_embedded;
  sint32         maxexp = 1024;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_dimension_correlation_summation()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  err = frauniv_check_embedding_inputs(
    data, dimension, time_lag, orbital_lag,
    intrp_ptr, &is_delay_embedding, &n_embed );
  MEMLIST_FREE_ON_ERROR( err, &list );

  is_embedded = (boolean) ( ( MATUNIV_NCOL( data ) > 1 ) &&
    ( MATUNIV_NROW( data ) > 1 ) );

  if ( is_embedded ){
    n_embed = MATUNIV_NROW( data );
    dimension = MATUNIV_NCOL( data );
  }
  else{

    n_sample = MATUNIV_NELEM( data );
    max_lag  = (sint32) ( (dimension - 1) * time_lag );
    n_embed = n_sample - max_lag;
  }

  if ( dimension < 2 ){
    MUTIL_ERROR( "Embedding dimension must exceed one." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( orbital_lag < 0 ){
    MUTIL_ERROR( "Trajectory lag must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( orbital_lag > n_embed ){
    MUTIL_ERROR( "Trajectory lag is too large. "
      "Decrease time lag and/or embedding dimension" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( resolution < 1 ){
    MUTIL_ERROR( "Number of scales per octave must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* define scale-based variables */

  n_scale    = maxexp + exp_offset + 1;
  iscale_min = n_scale - 1;
  iscale_max = (sint32) 0;

  /* allocate space and register with memory manager */

  err = matuniv_malloc_register( &C2, n_scale, dimension, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &scale, n_scale, 1, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create pointers */

  pd_data  = data->mat.dblmat.data;
  pd_corr  = C2.mat.dblmat.data;
  pd_scale = scale.mat.dblmat.data;

  /* initialize variables */

  for ( j = 0, i = -exp_offset; i <= maxexp; i++ ){
    pd_scale[ j++ ] = (double) i / (double) resolution;
  }

  for ( i = 0; i < ( n_scale * dimension ); i++ ){
    pd_corr[ i ] = (double) 0.0;
  }

  /* calculate the correlation sum */

  for ( i = 0; i < n_embed; i++ ){

    for ( j = i + 1 + orbital_lag; j < n_embed; j++ ){

      /* set pointers to the first coordinates of two
	 points in the phase space */

      if ( is_embedded ){

	pd_x = pd_data + i * dimension;
	pd_y = pd_data + j * dimension;
      }
      else{

	pd_x = pd_data + i;
	pd_y = pd_data + j;
      }

      /* distance_sum = (double) 0.0; */

      infinity_norm = (double) 0.0;

      for ( E = 0; E < dimension; E++){

	/* calculate coordinate separation */

	distance = (float) fabs( *pd_x - *pd_y );

        /* raise distance to the resolution power
	   to (effectively) increase the number of
	   scales of nonzero correlation summations
           by a factor of the resolution */

	localdef_ipower( &distance, resolution );

	/* use L-infinity norm for point separation */

	if ( distance > infinity_norm ) infinity_norm = distance;

	/* NOTE: THIS SECTION IS UNSTABLE FOR SOME PLATFORMS AND
           COMPILER COMBINATIONS. THEREFORE, THE DISTANCE
           EXPONENT IS CALCULATED USING THE STANDARD MUTIL_LOG2 MACRO.

	   calcuate scale index based on fast bitwise
	   masking and shift operations to reveal
	   base-2 exponent of point separation.
	   these values are biased by +127 so that
	   the index of scale can be found directly
	   without exceeding the bounds of the vector: */

        /* 	iscale = (sint32) localdef_exponent_biased( (float) infinity_norm ); */

	iscale = (sint32) ( MUTIL_LOG2( infinity_norm ) ) + exp_offset;

	/* perform index cheks just in case offset and maximum exponent are incorrect */

        if ( iscale < 0 ) iscale = 0;

        if ( iscale > n_scale - 1) iscale = n_scale - 1;

	/* record maximum and minimum scale index */

        if ( iscale < iscale_min ) iscale_min = iscale;

        if ( iscale > iscale_max ) iscale_max = iscale;

	/* increment the correlation sum at the appropriate scale */

	pd_corr[ iscale * dimension + E ] += (double) 1.0;

	/* increment pointers to next dimension */

	if ( is_embedded ){

	  pd_x++;
	  pd_y++;
	}
	else{

	  pd_x += time_lag;
	  pd_y += time_lag;
	}
      }
    }
  }

  /* further refine the iscale_min index
     by going through each dimenions of the
     correlation summation and finding the minimum
     index where there are no "gaps" (zeros) */

  err = mats32_malloc_register( &ijscale_min, 1, dimension, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  for ( j = 0; j < dimension; j++ ){
    ijscale_min.data[ j ] = iscale_min;
  }

  for ( i = iscale_min; i < iscale_max; i++ ){

    for ( j = 0; j < dimension; j++ ){

      if ( pd_corr[ i * dimension + j ] < 1.0 ){
		ijscale_min.data[ j ] = i + 1;
      }
    }
  }

  /* set iscale_min to that corresponding
     to the first dimension */

  iscale_min = ijscale_min.data[0];

  /* allocate space for the output matrix
     set, register with the memory manager,
     and copy data corresponding to usable
     scales into the nre matrix set matrices */

  n_scale = iscale_max - iscale_min + 1;

  err = mats32_malloc_register( &nrow, 1, 2, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &ncol, 1, 2, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* ... the correlation matrix */

  nrow.data[0] = n_scale;
  ncol.data[0] = dimension;

  /* ... the scale vector */

  nrow.data[1] = n_scale;
  ncol.data[1] = (sint32) 1;

  /* ... allocate space for the matrix set and matrices */

  dims = 2;
  err = matset_malloc_register( result, 1, &dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices_arbitrary_size(
    result, &nrow, &ncol, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* ... extract usable portions of each matrix
     and put into matrix set */

  err = matuniv_extract( &C2, iscale_min, 0, intrp_ptr, &result->mats[0] );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_extract( &scale, iscale_min, 0, intrp_ptr, &result->mats[1] );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* ... free old memory */

  err = memlist_member_free( &C2, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_free( &scale, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* reassign pointers */

  pd_corr = result->mats[0].mat.dblmat.data;

  /* normalize the results by doing a cumulative
     summation over each column in the correlation
     summation matrix followed by the standard
     normalization */

  dfac = (double) ( n_embed - orbital_lag );

  norminv = 2.0 / ( ( dfac - 1.0 ) *  dfac );

  for ( j = 0; j < dimension; j++ ){

    if ( MUTIL_INTERRUPT( 3.0 * n_embed, intrp_ptr ) ) {
         MUTIL_ERROR( "user interrupt" );
         MUTIL_FREE_WARN( memlist, &list );
         return MUTIL_ERR_INTERRUPT;
    }

    cum_sum = (double) 0.0;

    pd_C = pd_corr + j;

    k = ijscale_min.data[ j ] - iscale_min;

    for ( i = 0; i < n_scale; i++ ){

      if ( i >= k ){

	cum_sum += *pd_C;

	*pd_C = (double) MUTIL_LOG2( cum_sum * norminv );

      }
      else{

	/* set the current log2(C(r,E)) value to
	   be greater than unity indicating a
	   bad (impossible) value to the caller */

	*pd_C = MUTIL_DOUBLE_MAX;
      }

      pd_C += dimension;
    }
  }

  /* free nodes corresponding to registered
  memory for the result, but do not free
  the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
  corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with frauniv_dimension_correlation_summation()" );

  return MUTIL_ERR_OK;
}



/* Finding the proper embedding dimesion for     */
/* single variable time series ala the method    */
/* of false nearest neighbors                    */
/*                                               */
/* Documented in fra_dim.h                       */
/* Written by William Constantine                */

mutil_errcode frauniv_dimension_false_nearest_neighbors(
  const univ_mat *time_series,
  const sint32    embedding_dimension,
  const sint32    time_lag,
  const sint32    orbital_lag,
  const double    rtol,
  const double    atol,
  void           *intrp_ptr,
  univ_mat       *result )
{
  boolean        fail;
  boolean        is_delay_embedding;
  boolean        sort_distances = FALSE;
  double         attractor_size;
  double         extra_distance;
  double         mean;
  double         total_distance;
  double         variance;
  double        *pd_result;
  fra_distance_metric distance_metric = FRA_DISTANCE_L2;
  memlist        list;
  mutil_errcode  err;
  sint32         dim;
  sint32         fnn_atol;
  sint32         fnn_either;
  sint32         fnn_rtol;
  sint32         i;
  sint32         lag;
  sint32         n_embed;
  univ_mat       embedding;
  univ_mat       neighbor_distances;
  univ_mat       neighbor_indices;
  univ_mat       original_indices;
  sint32         ioriginal;
  sint32         ineighbor;
  sint32         n_sample;
  sint32         n_usable;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_dimension_false_nearest_neighbors()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  err = frauniv_check_embedding_inputs(
    time_series, embedding_dimension, time_lag, orbital_lag,
    intrp_ptr, &is_delay_embedding, &n_embed );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( orbital_lag < 1 ){
    MUTIL_ERROR("Orbital lag in FNN must be at least unity");
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( rtol < 1 || atol < 1 ){
    MUTIL_ERROR("FNN tolerances rtol and atol must be at least unity");
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* allocate memory */

  err = matuniv_malloc_register( result, 3, embedding_dimension, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* estimate attractor size */

  err = matuniv_mean_variance( time_series, (boolean) FALSE, intrp_ptr, &mean, &variance );
  MEMLIST_FREE_ON_ERROR( err, &list );

  attractor_size = sqrt( variance );

  /* set pointers */

  pd_result = result->mat.dblmat.data;

  /* calculate false nearest neighbors */

  for ( dim = 1; dim <= embedding_dimension; dim++ ){

    /* initialize variables */

    fnn_rtol   = 0;
    fnn_atol   = 0;
    fnn_either = 0;

    n_sample = MATUNIV_NELEM( time_series );

    /*
      embed the data in the current embedding dimension and register the
      result with the memory manager
    */

    if ( memlist_member_exist( &embedding, &list ) ){
      err = memlist_member_free( &embedding, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    err = frauniv_embed(
      time_series,
      dim,
      time_lag,
      intrp_ptr,
      &embedding );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_register( &list, &embedding, MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, &embedding );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    /* find nearest neighbor for each point in the embedding */

    if ( memlist_member_exist( &neighbor_indices, &list ) ){
      err = memlist_member_free( &neighbor_indices, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    if ( memlist_member_exist( &neighbor_distances, &list ) ){
      err = memlist_member_free( &neighbor_distances, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    err = frauniv_neighbor_find(
      &embedding,
      (sint32) 1,
      (double) 0.0,
      distance_metric,
      (univ_mat*) NULL,
      sort_distances,
      orbital_lag,
      intrp_ptr,
      &original_indices,
      &neighbor_indices,
      &neighbor_distances );
    MEMLIST_FREE_ON_ERROR( err, &list );

   /* free superfluous matrices */

    MUTIL_FREE_WARN( matuniv, &original_indices );

    /* register the nearest neighbor matrices with the memory manager */

    err = memlist_member_register( &list, &neighbor_indices, MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, &neighbor_indices );
      MUTIL_FREE_WARN( matuniv, &neighbor_distances );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    err =  memlist_member_register( &list, &neighbor_distances, MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, &neighbor_distances );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    /* perform false nearest neighbor tests */

    n_embed = MATUNIV_NELEM( &neighbor_indices );

    n_usable = 0;

    for ( i = 0; i < n_embed; i++ ){

      fail = FALSE;

      /*
	calculate indices of the time series which
	correspond to the neighbors projected up
	 one dimension
      */

      lag       = dim * time_lag;
      ioriginal = i + lag;
      ineighbor = neighbor_indices.mat.s32mat.data[ i ] + lag;

      /*
	as we increase the embedding dimension, a delay
	embedding will decrease the number of available
	points. check here to ensure that the indices do not
	exceeded the number of available points in the next
	dimension. this is equivalent to making sure that
	the corresponding indices of the original time series
	do not exceed N - 1, where N is the length of the
	time series.
      */

      if ( ioriginal < n_sample && ineighbor < n_sample ){

	extra_distance = MUTIL_ABS( time_series->mat.dblmat.data[ ioriginal ] -
	  time_series->mat.dblmat.data[ ineighbor ] );

	/* perform extra distance test */

	if ( ( extra_distance / neighbor_distances.mat.dblmat.data[ i ] ) > rtol ){

	  fnn_rtol++;
	  fail = TRUE;
	}

	/* perform attractor distance test */

	total_distance = sqrt( MUTIL_SQR( neighbor_distances.mat.dblmat.data[ i ] ) +
	  MUTIL_SQR( extra_distance ) );

	if ( ( total_distance / attractor_size ) > atol ){
	  fnn_atol++;
	  fail = TRUE;
	}

	if ( fail ) fnn_either++;

	n_usable++;
      }

    } /* end testing loop */

    /* record test results */

    if ( n_usable > 0 ){

      pd_result[ dim - 1 ] = (double) fnn_rtol / (double) n_usable * 100.0;
      pd_result[ dim - 1 + embedding_dimension ] = (double) fnn_atol / (double) n_usable * 100.0;
      pd_result[ dim - 1 + 2 * embedding_dimension ] = (double) fnn_either / (double) n_usable * 100.0;
    }
    else{
      pd_result[ dim - 1 ]                           = (double) - 1.0;
      pd_result[ dim - 1 + embedding_dimension ]     = (double) - 1.0;
      pd_result[ dim - 1 + 2 * embedding_dimension ] = (double) - 1.0;
    }

    /* check for interrupts */

    if ( MUTIL_INTERRUPT( 3.0 * n_embed, intrp_ptr ) ) {
         MUTIL_ERROR( "user interrupt" );
         MUTIL_FREE_WARN( memlist, &list );
         return MUTIL_ERR_INTERRUPT;
    }

  } /* end loop over each dimension */

  /* free nodes corresponding to registered
  memory for the result, but do not free
  the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
  corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with frauniv_dimension_false_nearest_neighbors()" );

  return MUTIL_ERR_OK;
}


/* Poincare map via local extrema */
/*                                */
/* Documented in fra_dim.h        */
/* Written by Keith L. Davidson   */
mutil_errcode frauniv_poincare_map(
  const univ_mat          *time_series,
  const fra_extrema_type   extrema_type,
  const boolean            denoise,
  void                    *intrp_ptr,
  univ_mat                *extrema_location,
  univ_mat                *extrema_amplitude )
{
  mat_set                filters;
  memlist                mlist;
  mutil_errcode          err;
  sint32                 n_level;
  sint32                 n_sample;
  sint8                 *concavity;
  univ_mat               ts_denoised;
  univ_mat              *ts = (univ_mat*) time_series;
  wav_shrink_function    ws_function;
  wav_shrink_threshold   ws_treshold;

  MUTIL_TRACE( "Start in frauniv_poincare_map()" );

  /* Avoid lint warning */

  (void) whatssi;

  /* Initialize memory management list */

  MEMLIST_INIT( mlist );

  /* Error checks ... */

  /* ... time series */

  err = matuniv_validate( time_series );
  if ( err ) {
    MUTIL_ERROR( "Input time_series is an invalid matrix" );
    return err;
  }
  if ( time_series->type != MUTIL_DOUBLE ) {
    MUTIL_ERROR( "Input time_series must be of type MUTIL_DOUBLE" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  n_sample = MATUNIV_NELEM( time_series );

  if ( !MATANY_IS_VEC( &( time_series->mat.dblmat ) ) ){
    MUTIL_ERROR( "Input time_series must be a row or column vector" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if ( n_sample < FRA_POIN_MIN_N_SAMPLE ) {
    MUTIL_ERROR( "Not enough elements in the input time series" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* extrema type */

  switch( extrema_type ) {
  case FRA_EXTREMA_MINIMA:
  case FRA_EXTREMA_MAXIMA:
  case FRA_EXTREMA_ALL:
    break;
  default:
    MUTIL_ERROR( "Extrema type unsupported" );
    return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* If we're denoising then call the WaveShrink function */

  if ( denoise ) {

    /* get WaveShrink default parameters */

    err = localfn_create_waveshrink_defaults(
      n_sample,
      intrp_ptr,
      &filters,
      &n_level,
      &ws_function,
      &ws_treshold );
    if ( err ) return err;

    /* perform denoising with WaveShrink algorithm */

    err = wavuniv_shrink(
      time_series,
      &filters,
	    (univ_mat *) NULL,
      ws_treshold,
      (double) 1.0,
      (double) -1.0,
      ws_function,
      n_level,
	    (boolean) TRUE,
      intrp_ptr,
      &ts_denoised );

    /* free filters from memory */

    MUTIL_FREE_WARN( matset, &filters );

    /* check for error during WaveShrink function */

    if ( err ) {
      MUTIL_ERROR( "Function wavuniv_shrink() failed" );
      return err;
    }

    /* update time series pointer */

    ts = &ts_denoised;

    /* add denoised time series to the memory management list */

    err = memlist_member_register( &mlist, (void*) &ts_denoised, MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, &ts_denoised );
      MUTIL_ERROR( "Could not add denoised time series to memory manager" );
      return err;
    }

  } /* if ( denoise ) */

  /* If ALL extrema are to be found then allocate memory to hold their
     concavity ( 1 = maxima, -1 = minima */
  sint8 concval = 1;
  concavity = &concval; /* BC: SET DEFAULT TO AVOID COMPILER WARNING */

  if ( extrema_type == FRA_EXTREMA_ALL ) {
    err = mutil_malloc_register(
      n_sample * sizeof( sint8 ),
      (void**) &concavity,
      &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );
  }

  /* Find the minima (or maxima) locations of the time series */

  err = localfn_extrema_locations(
    ts,
    extrema_type,
    intrp_ptr,
    extrema_location,
    concavity );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  /* Add extrema locations to memory management list */

  err = memlist_member_register( &mlist, (void*) extrema_location, MEMTYPE_MATUNIV );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  /* Allocate memory for output map */

  err = matuniv_malloc_register(
    extrema_amplitude,
    MATUNIV_NROW( extrema_location ),
    MATUNIV_NCOL( extrema_location ),
    MUTIL_DOUBLE,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  /* Compute interpolated time series values at the extrema points */

  err = localfn_extrema_amplitudes(
    time_series,
    extrema_location,
    extrema_type,
    concavity,
    extrema_amplitude );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  /* Remove output matrices from memory management list */

  err = memlist_member_unregister( (void*) extrema_location, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = memlist_member_unregister( (void*) extrema_amplitude, &mlist );
  if ( err ) {
    MUTIL_FREE_WARN( matuniv, extrema_location );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }

  /* Free all other allocated memory */

  MUTIL_FREE_WARN( memlist, &mlist );

  MUTIL_TRACE( "Done with frauniv_poincare_map()" );

  return MUTIL_ERR_OK;
}

/* Space time separation plot.  */
/*                              */
/* Documented in fra_visu.h     */
/* Written by Keith L. Davidson */
mutil_errcode frauniv_space_time_separation_plot(
  const univ_mat   *time_series,
  const sint32      dim,
  const sint32      delay,
  const univ_mat   *orbital_lags,
  const double      frac_prob,
  void             *intrp_ptr,
  double           *max_eps,
  univ_mat         *result )
{
  double         *hist;
  double         *pt1_ptr;
  double         *pt2_ptr;
  double         *res_ptr;

  double          eps_max;
  double          max  = -MUTIL_DOUBLE_MAX;
  double          min  = MUTIL_DOUBLE_MAX;

  sint32          cumsum;
  sint32          hist_idx;
  sint32          i;
  sint32          limit;
  sint32          m;
  sint32          n;
  sint32          nembed;
  sint32          nfracs;
  sint32          nolags;
  sint32          target;
  sint32          tmp_idx;
  sint32          ts_length;

  const sint32    neps = FRAUNIV_STSP_NUM_EPSILON;

  memlist         mlist;
  mutil_errcode   err;

  MUTIL_TRACE( "Start in frauniv_space_time_separation()" );

  /* Avoid lint warning */
  (void) whatssi;
  (void) intrp_ptr;

  /* Initialize memory management list */
  MEMLIST_INIT( mlist );


  /* Error checks */
  err = localfn_stsp_error_checks(
    time_series,
    delay,
    dim,
    orbital_lags,
    frac_prob,
    &ts_length,
    &nfracs,
    &nembed,
    &nolags );
  if ( err ) return err;


  /* Find the maximum epsilon. Here, this is the maximum possible infinity
     norm distance between any pair, over all possible embeddings. */
  pt1_ptr = (double*) MATUNIV_DATA( time_series );
  for ( n = 0; n < ts_length; n++, pt1_ptr++ ) {
    min = MUTIL_MIN( min, *pt1_ptr );
    max = MUTIL_MAX( max, *pt1_ptr );
  }
  eps_max  = max - min;
  *max_eps = eps_max;


  /* Check for constant time series (shouldn't happen) */
  if ( MUTIL_EQUAL_TO( eps_max + 1, 1.0, 30 ) ) {
    MUTIL_ERROR( "Input time series appears to be of constant values" );
    return MUTIL_ERR_SINGULARITY;
  }


  /* Allocate memory for the epsilon-histogram */
  err = mutil_malloc_register( neps * sizeof(double), (void**) &hist, &mlist );
  if ( err ) {
    MUTIL_ERROR( "Memory allocation error" );
    return err;
  }


  /* Allocate memory for the output */
  err = matuniv_malloc_register( result, nfracs, nolags, MUTIL_DOUBLE, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );


  /* Set a pointer to the first embedding point */
  pt1_ptr = (double*) MATUNIV_DATA( time_series );


  /* Begin looping through the orbital lags */
  for ( m = 0; m < nolags; m++ ) {

    /* zero the histogram bins */
    for ( i = 0; i < neps; i++ )
      hist[ i ] = (double) 0;

    /* set a second to the point that's orbital_lags(m) time steps beyond
       the first embedding point */
    pt2_ptr = pt1_ptr + orbital_lags->mat.s32mat.data[ m ];

    /* now we move through the embedding in temporal order, and always
       staying orbital_lags(m) time steps apart */
    limit = nembed - orbital_lags->mat.s32mat.data[ m ] - 1;

    for ( n = 0; n < limit; n++ ) {

      /* calculate the infinity norm of the current embedding vector pair */
      max = (double) 0;
      for ( tmp_idx = 0, i = 0; i < dim; i++, tmp_idx += delay ) {
        max = MUTIL_MAX(
          max,
          MUTIL_ABS( pt1_ptr[ tmp_idx ] - pt2_ptr[ tmp_idx ] ) );
      }
      pt1_ptr++;
      pt2_ptr++;

      /* calculate the histogram index */
      hist_idx = (sint32) ( ( neps - 1 ) * max / eps_max );

      /* add 1 to the appropriate histogram bin */
      hist[ hist_idx ]++;

    } /*  for ( n = 0; ... ) */

    /* reset pointer to first embedding point */
    pt1_ptr = (double*) MATUNIV_DATA( time_series );

    /* run a cumulative sum of the histogram, along the way filling in the
    contour line values for the current orbital lag */
    cumsum   = (sint32) 0;
    hist_idx = (sint32) 0;
    res_ptr  = (double*) MATUNIV_DATA( result ) + m;

    for ( n = 0; n < nfracs; n++ ) {

      target = (sint32) ( limit * ( n + 1 ) * frac_prob );

      for ( i = hist_idx; i < neps; i++ ) {

		cumsum += (sint32) hist[ hist_idx ];

        if ( cumsum >= target ) {

		  /* fills column of result matrix, each column corresponds
		  to the contour values for a particular orbital lag */
          *( res_ptr + n * nolags ) = eps_max * ( ( i + 1 ) / (double) neps );

          break;
        }

/*        cumsum += (sint32) hist[ hist_idx ]; */

        hist_idx++;

      } /* for ( i = hist_idx; ... ) */

    } /* for ( n = 0; ... ) */


  } /* for ( m = 0; ... ) */


  /* Remove the output matrix from the memory management list */
  err = memlist_member_unregister( result, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );


  /* Free all other memory */
  MUTIL_FREE_WARN( memlist, &mlist );

  MUTIL_TRACE( "Done with frauniv_space_time_separation()" );

  return MUTIL_ERR_OK;
}

/*******************************/
/* STATIC FUNCTION DEFINITIONS */
/*******************************/

/** Average radius of phase space neighborhoods with constant density.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_dim.c
 * @library fractal
 * @usage #err = localfn_average_neighborhood_radius_constant_density( &embedding, &kdtree, orbital_lag, k, distance_metric, n_reference, intrp_ptr, &radius );#
 * @return Standard mutils error code.
 * @param embedding A pointer to a universal matrix of type
 *        MUTIL\_DOUBLE containing the time delayed embedding
 *        matrix with each column a different lagged coordinate..
 * @param kdtree Pointer to a kd-tree object corresponding to
 *        the embedding.
 * @param orbital_lag The number of points along the trajectory of the
 *                    current point that must be exceeded in order for
 *                    another point in the phase space to be considered
 *                    a neighbor candidate. This argument is used
 *                    to help attenuate temporal correlation in the
 *                    the embedding which can lead to spuriously low
 *                    information dimension estimates. The orbital lag
 *                    must be positive or zero.
 * @param n_neighbor  The number of neighbors used to form the kd-tree.
 *                    This parameter is adjusted in the caller to achieve
 *                    a given density by density =
 *                    n\_neighbor / n_reference.
 * @param distance_metric The metric used to define the distance between
 *                    points in the embedding. The argument is of
 *                    type \Ref{_fra_distance_metric}.
 * @param n_reference The number of points in the embedding over which the
 *                    average neighborhood radius should be estimated
 *                    given a specified density.
 * @param intrp_ptr   Interrupt pointer.
 * @param radius      A double value representing the estimated radius of
 *                    all neighborhoods with e specified density..
 * @see frauniv_dimension_information
 * @private
 */
static mutil_errcode localfn_average_neighborhood_radius_constant_density(
  univ_mat            *embedding,
  const mutil_kdtree  *kdtree,
  const sint32         orbital_lag,
  const sint32         n_neighbor,
  const fra_distance_metric distance_metric,
  const sint32         n_reference,
  void                *intrp_ptr,
  double              *average_log2_radius )
{
  double        *pd_dist;
  double        *pd_embed;
  double        *pd_sampled;
  sint32         nrow  = MATUNIV_NROW( embedding );
  sint32         ncol  = MATUNIV_NCOL( embedding );
  //sint32         nelem = MATUNIV_NELEM( embedding );
  univ_mat       original_indices;
  univ_mat       neighbor_indices;
  univ_mat       neighbor_distances;
  univ_mat       sampled_embedding;
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         j;
  sint32         iskip;
  boolean        sort_distances = TRUE;
  sint32         max_index;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start localfn_average_neighborhood_radius_constant_density()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* allocate memory */

  err = matuniv_malloc_register( &sampled_embedding, n_reference, ncol,
    MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* uniformly sample the embedding to obtain n_reference reference points */

  max_index  = nrow - 1 - 2 * orbital_lag - n_neighbor;
  iskip      = (sint32) max_index / ( n_reference - 1 );
  if ( iskip == 0 ) iskip = 1;
  iskip *= ncol;

  pd_embed   = embedding->mat.dblmat.data;
  pd_sampled = sampled_embedding.mat.dblmat.data;

  for ( i = 0; i < n_reference; i++ ){

     for ( j = 0; j < ncol; j++ ){

       *pd_sampled = *pd_embed;
       pd_sampled++;
       pd_embed++;
     }

     pd_embed += iskip;
  }

  /* find n_neighbor nearest neighbors for the uniformly
  sampled reference points */

  err = frauniv_neighbor_find_arbitrary(
    &sampled_embedding,
    kdtree,
    n_neighbor,
    (double) 0.0,
    distance_metric,
    (univ_mat*) NULL,
    sort_distances,
    orbital_lag,
    (univ_mat *) NULL,
    intrp_ptr,
    &original_indices,
    &neighbor_indices,
    &neighbor_distances );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free superfluous data */

  MUTIL_FREE_WARN( matuniv, &original_indices );
  MUTIL_FREE_WARN( matuniv, &neighbor_indices );

  err = memlist_member_free( &sampled_embedding, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /*
     calculate the average neighborhood size.
     since the distances are sorted from smallest
     to largest, and since the number of neighbors found
     for each reference point is constant, we can simply
     uniformly subsample the neighbor_distances vector
     to obtain each neighborhood size.
  */

  *average_log2_radius = 0.0;

  /* ... point to the first distance */

  pd_dist = &( neighbor_distances.mat.dblmat.data[ n_neighbor - 1 ] );

  for ( i = 0; i < n_reference; i++ ){

    if ( *pd_dist > (double) 0.0 ){

      *average_log2_radius += MUTIL_LOG2( *pd_dist );
    }

    pd_dist += n_neighbor;

    if ( MUTIL_INTERRUPT( 3.0 * nrow, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( matuniv, &neighbor_distances );
      return MUTIL_ERR_INTERRUPT;
    }
  }

  *average_log2_radius /= (double) n_reference;

  MUTIL_FREE_WARN( matuniv, &neighbor_distances );

  /* Free all other memory */
  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_average_neighborhood_radius_constant_density()" );

  return MUTIL_ERR_OK;
}

/** In-place fast algorithm for raising a float to a positive integer power.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_dim.c
 * @library fractal
 * @usage #x = localdef_ipower( y, n );#
 * @return Nothing (void)
 * @param  y Pointer to a float value.
 *           Upon return y will contain y to the nth power.
 * @param  n Positive integer power.
 *
 * @see frauniv_dimension_correlation_summation
 * @private
 */
static void localdef_ipower( float *x, sint32 n )
{
  float aux = 1.0;

  if ( n == (sint32) 1 ) return;

  while ( n > 1 ){
    if ( IS_EVEN(n) ){
      *x *= *x;
      n = (sint32) ( (double) n * 0.5 );
    }
    else{
      aux *= *x;
      *x  *= *x;
      n = (sint32) ( (double) n * 0.5 );
    }
  }
  *x *= aux;

  return;
}

/** Bitwise masking and shift operations to extract exponent of floating point values
 * based on IEEE Standard 754 Floating Point Numbers.
 *
 * NOTE: After bitwise masking and shifting, the result (representing the exponent)
 *       is still offset by +127. So the true base-2 exponent can be obtained
 *       by subtracting 127 from the result of this function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_dim.c
 * @library fractal
 * @usage #x = localdef_exponent_biased( y );#
 * @return sint32 value representing the BIASED exponent.
 * @param  y Float value.
 *
 * @see frauniv_dimension_correlation_summation
 * @private
 */
/* static uint32 localdef_exponent_biased( float num ){ */

/*   /\* sint32 offset  = -127; *\/ */
/*   /\* Mask 0x7f800000: 0111 1111 1000 0000 0000 0000 0000 0000 *\/ */

/*   uint32 mask  = 0x7F800000; */
/*   /\*   uint32 shift = 23; *\/ */

/*   float  *pfNum = &num; */
/*   uint32 *plNum; */

/*   plNum = (uint32 *) pfNum; */

/*   return ( (*plNum & mask) >> 23 ); */
/* } */

/** Establishes the default values for for waveshrink options.
 * Waveshrink is used by the \Ref{frauniv_poincare_map}.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_dim.c
 * @library fractal
 * @usage #err = localfn_create_waveshrink_defaults( n_sample, intrp_ptr, &filters, &n_level, &ws_function, &ws_treshold );#
 * @return Standard mutils error code.
 * @param n_sample The numbe rof samples in the time series.
 * @param intrp_ptr Interrupt pointer.
 * @param filters Pointer to a matrix set to hold the wavelet and scaling
 *                filters.
 * @param n_level Pointer to an integer denoting the number of DWT
 *                decomposition levels.
 * @param ws_fun  Pointer to waveshrink function type.
 * @param ws_tsh  Pointer to waveshrink threshold type.
 *
 * @see frauniv_poincare_map
 * @see wavuniv_shrink
 * @private
 */
static mutil_errcode localfn_create_waveshrink_defaults(
  const sint32            n_sample,
  void                   *intrp_ptr,
  mat_set                *filters,
  sint32                 *n_level,
  wav_shrink_function    *ws_fun,
  wav_shrink_threshold   *ws_tsh )
{
  mutil_errcode   err;

  MUTIL_TRACE( "Start localfn_create_waveshrink_defaults()" );

  /* Determine the number of decomposition levels */

  *n_level = (sint32) floor( (double) MUTIL_LOG2( (double) n_sample ) );

  /* Set default threshold and shrinkage function types */

  *ws_fun = FRA_POIN_WAVSHRK_DEF_FUNC;
  *ws_tsh = FRA_POIN_WAVSHRK_DEF_TSH;

  /* Compute the default wavelet filters */

  err = wavuniv_filters_daubechies(
    FRA_POIN_WAVSHRK_DEF_FILT_LEN,
    FRA_POIN_WAVSHRK_DEF_FILT,
    FALSE,
    intrp_ptr,
    filters );
  if ( err ) {
    MUTIL_ERROR( "Could not create default wavelet filters" );
    return err;
  }

  MUTIL_TRACE( "Done with localfn_create_waveshrink_defaults()" );

  return MUTIL_ERR_OK;
}


/** Estimate local extrema locations of a time series.
 * A first derivative approximation is obtained
 * by extracting the continuous wavelet transform (CWT)
 * coefficients at an appropriate scale.
 * A second derivative is also estimated in
 * the same way by taking the CWT of the first derivative estimate.
 * The first derivative is used to find
 * zero crossings (identifying local extrema in the original
 * series) while the second derivative is used to infer
 * concavity (separating minima from maxima).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_dim.c
 * @library fractal
 * @usage #err = localfn_extrema_locations(ts,extrema_type,intrp_ptr,extrema_location,concavity );#
 * @return Standard mutils error code.
 * @param time_series Pointer to single-column or single-row
 *        universal matrix of type MUTIL\_DOUBLE containing
 *        the time series.
 * @param extrema_type The extrema type (minima, maxima, or both).
 * @param intrp_ptr Interrupt pointer.
 * @param extrema Pointer to a MUTIL\_DOUBLE universal matrix containing the
 *        estimated extrema locations. The memory for this matrix is
 *        automatically allocated within the function.
 * @param concavity Pointer to an integer array whose elements are either
 *        1 or -1 indicating that the corresponding local extrema are
 *        maxima or minima, respectively.
 * @see frauniv_poincare_map
 * @private
 */
static mutil_errcode localfn_extrema_locations(
  const univ_mat          *time_series,
  const fra_extrema_type   extrema_type,
  void                    *intrp_ptr,
  univ_mat                *extrema,
  sint8                   *concavity )
{
  double             scale_value = FRA_POIN_CWT_SCALE;
  memlist            mlist;
  mutil_errcode      err;
  register double    tmp_dbl;
  register double   *d1_ptr;
  register double   *d2_ptr;
  register double   *extrema_ptr;
  register sint32    extrema_count;
  register sint32    limit;
  register sint32    n;
  univ_mat           cpx_mat;
  univ_mat           deriv1;
  univ_mat           deriv2;
  univ_mat           scale;

  MUTIL_TRACE( "Start localfn_extrema_locations()" );

  /* Initialize memory management list */

  MEMLIST_INIT( mlist );

  /* Set up the scale universal matrix (a 1 x 1 matrix) */

  scale.type             = MUTIL_DOUBLE;
  scale.mat.dblmat.nrow  = (sint32) 1;
  scale.mat.dblmat.ncol  = (sint32) 1;
  scale.mat.dblmat.nelem = (sint32) 1;
  scale.mat.dblmat.data  = &scale_value;

  /* Avoid lint warning */

  (void) scale_value;

  /* Compute the first derivative of the time series via the CWT */

  err = wavuniv_transform_continuous_wavelet(
    time_series,
    (double) 1,
    WAV_FILTER_GAUSSIAN_I,
    FRA_POIN_CWT_FILT_ARG,
    &scale,
    intrp_ptr,
    &cpx_mat );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = memlist_member_register( &mlist, (void*) &cpx_mat, MEMTYPE_MATUNIV );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = matuniv_malloc_register(
    &deriv1,
    MATUNIV_NROW( &cpx_mat ),
    MATUNIV_NCOL( &cpx_mat ),
    MUTIL_DOUBLE,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = matuniv_cast( &cpx_mat, TRUE, intrp_ptr, &deriv1 );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = memlist_member_free( (void*) &cpx_mat, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  /* Compute the second derivative of the time series via the CWT */

  err = wavuniv_transform_continuous_wavelet(
    &deriv1,
    (double) 1,
    WAV_FILTER_GAUSSIAN_I,
    FRA_POIN_CWT_FILT_ARG,
    &scale,
    intrp_ptr,
    &cpx_mat );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = memlist_member_register( &mlist, (void*) &cpx_mat, MEMTYPE_MATUNIV );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = matuniv_malloc_register(
    &deriv2,
    MATUNIV_NROW( &cpx_mat ),
    MATUNIV_NCOL( &cpx_mat ),
    MUTIL_DOUBLE,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = matuniv_cast( &cpx_mat, TRUE, intrp_ptr, &deriv2 );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = memlist_member_free( (void*) &cpx_mat, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  /* Allocate memory for the extrema */

  err = matuniv_malloc_register(
    extrema,
    MATUNIV_NROW( time_series ),
    MATUNIV_NCOL( time_series ),
    MUTIL_DOUBLE,
    &mlist );
  if ( err ) return err;

  /* Begin looping to find extrema */

  limit         = MATUNIV_NELEM( time_series );
  d1_ptr        = ( (double*) MATUNIV_DATA( &deriv1 ) ) + 1;
  d2_ptr        = ( (double*) MATUNIV_DATA( &deriv2 ) ) + 1;
  extrema_ptr   = (double*) MATUNIV_DATA( extrema );
  extrema_count = (sint32) 0;

  for ( n = 1; n < limit; n++, d1_ptr++, d2_ptr++ ){

    if ( *d1_ptr == (double) 0 ) {

      /* the first derivative is EXACTLY zero (unlikely with doubles) */

      switch ( extrema_type ) {

	case FRA_EXTREMA_MINIMA:

	  if ( *d2_ptr > (double) 0 ) {

	    *extrema_ptr = (double) n;

	    extrema_ptr++;
	    extrema_count++;
	  }

	  break;

	case FRA_EXTREMA_MAXIMA:

	  if ( *d2_ptr < (double) 0 ) {

	    *extrema_ptr = (double) n;

	    extrema_ptr++;
	    extrema_count++;
	  }

	  break;

	case FRA_EXTREMA_ALL:

	  *extrema_ptr = (double) n;

	  concavity[ extrema_count ] = ( *d2_ptr < (double) 0 ) ? (sint8) 1 : (sint8) -1;

	  extrema_ptr++;
	  extrema_count++;

	  break;

	default:

	  /* should never get here, check enum definition in fra_type.h */

	  MUTIL_ERROR( "Extrema type unsupported" );
	  MUTIL_FREE_WARN(memlist, &mlist );

	  return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;

      } /* end switch on extrema type */

    }
    else{

      if ( ( *d1_ptr * *( d1_ptr - 1 ) ) < (double) 0 ){

	/* we have a sign change. use linear interpolation to find
	   first derivative zero crossing */

	switch( extrema_type ) {

	  case FRA_EXTREMA_MINIMA:

	    if ( *d2_ptr > (double) 0 ) {

	      /* tmp_dbl cannot be zero because we're subtracting two
		 numbers with different signs */

	      tmp_dbl = *d1_ptr - *( d1_ptr - 1 );

	      *extrema_ptr = (double) n - *d1_ptr / tmp_dbl;

	      extrema_ptr++;
	      extrema_count++;
	    }

	    break;

	  case FRA_EXTREMA_MAXIMA:

	    if ( *d2_ptr < (double) 0 ) {

	      /* tmp_dbl cannot be zero because we're subtracting two
		 numbers with different signs */

	      tmp_dbl = *d1_ptr - *( d1_ptr - 1 );

	      *extrema_ptr = (double) n - *d1_ptr / tmp_dbl;

	      extrema_ptr++;
	      extrema_count++;
	    }

	    break;

	  case FRA_EXTREMA_ALL:

	    /* tmp_dbl cannot be zero because we're subtracting two
	       numbers with different signs */

	    tmp_dbl = *d1_ptr - *( d1_ptr - 1 );

	    *extrema_ptr = (double) n - *d1_ptr / tmp_dbl;

	    concavity[ extrema_count ] = ( *d2_ptr < (double) 0 ) ? (sint8) 1 : (sint8) -1;

	    extrema_ptr++;
	    extrema_count++;

	    break;

	  default:

	    /* should never get here, check fra_extrema_type enum in fra_type.h */

	    MUTIL_ERROR( "Extrema type unsupported" );
	    MUTIL_FREE_WARN( memlist, &mlist );

	    return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;

	} /* end switch on extrema_type */

      } /* if ( ( *ch1 & FRA_POIN_SIGN_MASK ) != ... ) */

    } /* if ( *d1_ptr == (double) 0 ) */

  } /* for ( n = 1; ... ) */

  /* Free the first and second derivative vectors */

  err = memlist_member_free( (void*) &deriv1, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = memlist_member_free( (void*) &deriv2, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  /* Check to see if any extrema were found */

  if ( extrema_count == (sint32) 0 ) {
    MUTIL_ERROR( "No extrema were found in time series" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return MUTIL_ERR_SINGULARITY;
  }

  /* Resize the extrema vector */

  err = matuniv_realloc_register( extrema, (sint32) 1, extrema_count, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  /* Make sure extrema is a row (column) vector if the input time series is
     a row (column) vector */

  if ( MATUNIV_NROW( time_series ) != (sint32) 1 ) {
    extrema->mat.dblmat.nrow = extrema_count;
    extrema->mat.dblmat.ncol = (sint32) 1;
  }

  /* Remove the output vector of extrema from the memory management list */

  err = memlist_member_unregister( (void*) extrema, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  /* Free all remaining memory */

  MUTIL_FREE_WARN( memlist, &mlist );

  MUTIL_TRACE( "Done with localfn_extrema_locations()" );

  return MUTIL_ERR_OK;
}


/** Estimate local extrema amplitudes of a time series.
 * Using a supplied vectors of local extrema locations,
 * this function performs a second-order polynomial
 * fit to estimate the amplitudes of the extrema.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_dim.c
 * @library fractal
 * @usage #err = localfn_extrema_amplitudes(time_series,extrema_location,extrema_type,concavity,extrema_amplitude );#
 * @return Standard mutils error code.
 * @param time_series Pointer to single-column or single-row
 *        universal matrix of type MUTIL\_DOUBLE containing
 *        the time series.
 * @param extrema_location The extrema location output of
 *        \Ref{localfn_extrema_locations}.
 * @param extrema_type The extrema type (minima, maxima, or both).
 * @param concavity Pointer to an integer array whose elements are either
 *        1 or -1 indicating that the corresponding local extrema are
 *        maxima or minima, respectively. This also is an output of
 *        \Ref{localfn_extrema_locations}.
 * @param extrema_amplitude Pointer to a pre-allocated
 *        MUTIL\_DOUBLE universal matrix that upon return will
 *        containing the estimated extrema amplitudes.
 * @see localfn_extrema_locations
 * @see frauniv_poincare_map
 * @private
 */
static mutil_errcode localfn_extrema_amplitudes(
  const univ_mat          *time_series,
  const univ_mat          *extrema_location,
  const fra_extrema_type   extrema_type,
  const sint8             *concavity,
  univ_mat                *extrema_amplitude )
{
  register double   *loc_ptr = (double*) MATUNIV_DATA( extrema_location );
  register double   *ts_ptr  = (double*) MATUNIV_DATA( time_series );
  register double   *val_ptr = (double*) MATUNIV_DATA( extrema_amplitude );
  register sint8    *con_ptr = (sint8*) concavity;

  register double    a;
  register double    b;
  register double    f_0;
  register double    tmp_dbl;

  register sint32    n;
  register sint32    nlocs   = MATUNIV_NELEM( extrema_location );
  register sint32    after;

  MUTIL_TRACE( "Start in localfn_extrema_amplitudes()" );

  /* Begin looping through extrema */

  switch ( extrema_type ) {

    case FRA_EXTREMA_MINIMA:

      for ( n = 0; n < nlocs; n++, loc_ptr++, val_ptr++ ) {

	if ( *loc_ptr == (double) ( (sint32) *loc_ptr ) ) {

	  /* the extrema location is an integer, i.e., it falls on one of the
	     original sample points */

	  *val_ptr = ts_ptr[ (sint32) *loc_ptr ];

	}
	else{

	  /* the extrema location is not an integer, i.e., it falls between
	     two of the original sampling points */

	  after = ( (sint32) *loc_ptr ) + 1;

	  /* fit a parabola around the two time series points on either side
	     of the extrema location */

	  f_0 = -MUTIL_ABS( ts_ptr[ (sint32) *loc_ptr ] -
	    ts_ptr[ after ] ) / 2.0;

	  f_0 += MUTIL_MIN( ts_ptr[ (sint32) *loc_ptr ],
	    ts_ptr[ after ] );

	  /* coefficients for y = a x^2 + b x + c (c is the time series value
	     to the left of the extrema location) */

	  a = 2.0 * ( ts_ptr[ after ] + ts_ptr[ (sint32) *loc_ptr ] ) -
	    4.0 * f_0;

	  b = 4.0 * f_0 - ts_ptr[ after ] - 3.0 * ts_ptr[ (sint32) *loc_ptr ];

	  tmp_dbl = *loc_ptr - (double) ( (sint32) *loc_ptr );

	  *val_ptr = a * tmp_dbl * tmp_dbl + b * tmp_dbl +
	    ts_ptr[ (sint32) *loc_ptr ];

	} /* if ( *loc_ptr == ... ) */

      } /* for ( n = 0; ... ) */

      break;

    case FRA_EXTREMA_MAXIMA:

      for ( n = 0; n < nlocs; n++, loc_ptr++, val_ptr++ ) {

	if ( *loc_ptr == (double) ( (sint32) *loc_ptr ) ) {

	  /* the extrema location is an integer, i.e., it falls on one of the
	     original sample points */

	  *val_ptr = ts_ptr[ (sint32) *loc_ptr ];

	} else {

	  /* the extrema location is not an integer, i.e., it falls between
	     two of the original sampling points */

	  after = ( (sint32) *loc_ptr ) + 1;

	  /* fit a parabola around the two time series points on either side
	     of the extrema location */

	  f_0 = MUTIL_ABS( ts_ptr[ (sint32) *loc_ptr ] -
	    ts_ptr[ after ] ) / 2.0;

	  f_0 += MUTIL_MAX( ts_ptr[ (sint32) *loc_ptr ],
	    ts_ptr[ after ] );

          /* coefficients for y = a x^2 + b x + c (c is the time series value
	     to the left of the extrema location) */
	  a = 2.0 * ( ts_ptr[ after ] + ts_ptr[ (sint32) *loc_ptr ] ) -
	    4.0 * f_0;

	  b = 4.0 * f_0 - ts_ptr[ after ] - 3.0 * ts_ptr[ (sint32) *loc_ptr ];

	  tmp_dbl = *loc_ptr - (double) ( (sint32) *loc_ptr );

	  *val_ptr = a * tmp_dbl * tmp_dbl + b * tmp_dbl +
	    ts_ptr[ (sint32) *loc_ptr ];

	} /* if ( *loc_ptr == ... ) */

      } /* for ( n = 0; ... ) */

      break;

    case FRA_EXTREMA_ALL:

      for ( n = 0; n < nlocs; n++, loc_ptr++, val_ptr++, con_ptr++ ) {

	if ( *loc_ptr == (double) ( (sint32) *loc_ptr ) ) {

	  /* the extrema location is an integer, i.e., it falls on one of the
	     original sample points */

	  *val_ptr = ts_ptr[ (sint32) *loc_ptr ];

	}
	else{

	  /* the extrema location is not an integer, i.e., it falls between
	     two of the original sampling points */

	  after = ( (sint32) *loc_ptr ) + 1;

	  /* fit a parabola around the two time series points on either side
	     of the extrema location */

	  if ( *con_ptr > 0 ) {

	    /* we're at a maxima */
	    f_0 = MUTIL_ABS( ts_ptr[ (sint32) *loc_ptr ] -
	      ts_ptr[ after ] ) / 2.0;

	    f_0 += MUTIL_MAX( ts_ptr[ (sint32) *loc_ptr ],
	      ts_ptr[ after ] );

	  } else {

	    /* we're at a minima */
	    f_0 = -MUTIL_ABS( ts_ptr[ (sint32) *loc_ptr ] -
	      ts_ptr[ after ] ) / 2.0;

	    f_0 += MUTIL_MIN( ts_ptr[ (sint32) *loc_ptr ],
	      ts_ptr[ after ] );
	  }

	  /* coefficients for y = a x^2 + b x + c (c is the time series value
	     to the left of the extrema location) */

	  a = 2.0 * ( ts_ptr[ after ] + ts_ptr[ (sint32) *loc_ptr ] ) -
	    4.0 * f_0;

	  b = 4.0 * f_0 - ts_ptr[ after ] - 3.0 * ts_ptr[ (sint32) *loc_ptr ];

	  tmp_dbl = *loc_ptr - (double) ( (sint32) *loc_ptr );

	  *val_ptr = a * tmp_dbl * tmp_dbl + b * tmp_dbl +
	    ts_ptr[ (sint32) *loc_ptr ];

	} /* if ( *loc_ptr == ... ) */

      } /* for ( n = 0; ... ) */

      break;

    default:

      /* should never get here, check fra_extrema_type enum in fra_type.h */

      MUTIL_ERROR( "Extrema type unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;

  } /* switch ( extrema_type ) */

  MUTIL_TRACE( "Done with localfn_extrema_amplitudes()" );

  return MUTIL_ERR_OK;
}

static mutil_errcode localfn_stsp_error_checks(
  const univ_mat     *time_series,
  const sint32        delay,
  const sint32        dim,
  const univ_mat     *orbital_lags,
  const double        frac_prob,
  sint32             *ts_length,
  sint32             *nfracs,
  sint32             *nembed,
  sint32             *nolags )
{
  sint32         *s32_ptr;

  sint32          n;

  mutil_errcode   err;

  MUTIL_TRACE( "Start in localfn_stsp_error_checks()" );


  /* time_series */
  err = matuniv_validate( time_series );
  if ( err ) {
    MUTIL_ERROR( "Input time_series is invalid" );
    return err;
  }

  if ( time_series->type != MUTIL_DOUBLE ) {
    MUTIL_ERROR( "Input time_series must be of type MUTIL_DOUBLE" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  *ts_length = MATUNIV_NELEM( time_series );
  if ( *ts_length != MATUNIV_NROW( time_series ) &&
       *ts_length != MATUNIV_NCOL( time_series ) ) {
    MUTIL_ERROR( "Input time_series must be a row or column vector" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }


  /* delay */
  if ( delay <= (sint32) 0 ) {
    MUTIL_ERROR( "Input delay must be a positive integer" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }


  /* dim */
  if ( dim <= (sint32) 0 ) {
    MUTIL_ERROR( "Input dim must be a positive integer" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }


  /* Calculate the number of embedding points */
  *nembed = *ts_length - ( dim - 1 ) * delay;
  if ( *nembed < (sint32) 2 ) {
    MUTIL_ERROR( "Not enough embedding points" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }


  /* orbital_lags */
  err = matuniv_validate( orbital_lags );
  if ( err ) {
    MUTIL_ERROR( "Input orbital_lags is invalid" );
    return err;
  }

  if ( orbital_lags->type != MUTIL_SINT32 ) {
    MUTIL_ERROR( "Input orbital_lags must be of type MUTIL_DOUBLE" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  *nolags = MATUNIV_NELEM( orbital_lags );
  if ( *nolags != MATUNIV_NROW( orbital_lags ) &&
       *nolags != MATUNIV_NCOL( orbital_lags ) ) {
    MUTIL_ERROR( "Input orbital_lags must be a row or column vector" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  s32_ptr = (sint32*) MATUNIV_DATA( orbital_lags );
  for ( n = 0; n < *nolags; n++, s32_ptr++ ) {

    if ( *s32_ptr < (sint32) 1 ) {
      MUTIL_ERROR( "Input orbital_lags must have positive elements" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

    if ( *s32_ptr > ( *nembed - 1 ) ) {
      MUTIL_ERROR( "Too few embedding points for given orbital lags" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
  }

  /* frac_prob */
  if ( MUTIL_EQUAL_TO(frac_prob + 1,1.0,30) || frac_prob >= (double) 1 ) {
    MUTIL_ERROR( "Input frac_prob must be in the interval (0,1)" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }


  /* Calculate the number of probability (contour) lines to be returned */
  *nfracs = 100 / (sint32) ( frac_prob * 100 );

  MUTIL_TRACE( "Done with localfn_stsp_error_checks()" );

  return MUTIL_ERR_OK;
}



static mutil_errcode localfn_kdtree(
  const sint32    n_neighbor,
  const sint32    max_neighbors,
  const univ_mat *embedding,
  const sint32    orbital_lag,
  const sint32    bucket_size,
  void           *intrp_ptr,
  mutil_kdtree   *kdtree,
  memlist        *list )
{
  sint32         n_embed_reduced;
  sint32         n_embed;
  sint32         n_embed_min;
  sint32         row;
  sint32         col;
  sint32         k;
  sint32         dim;
  univ_mat       sampled_embedding;
  univ_mat       random_index;
  double        *pd_sampled;
  double        *pd_embed;
  sint32        *ps_index;
  memlist        mlist;
  void          *rand_ptr;
  mutil_errcode  err;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start localfn_kdtree()" );

  /* initialize local memory list (this differs from the
  one passed into the function that primarily is used to
  register the resulting kd tree */

  MEMLIST_INIT( mlist );

  /* initialize variables */

  k       = n_neighbor;
  n_embed = MATUNIV_NROW( embedding );
  dim     = MATUNIV_NCOL( embedding );

  /* if the number of requested neighbors exceeds specified maximum,
  then adjust the number of points in the embedding to achieve the
  desired density level */

  if ( k > max_neighbors ){

    /* calculate the reduced size of the embedding in order to
    accommodate the current density given a fixed number of neighbors (k).
    we have to be careful here to also accommodate a restriction on the
    size of the embedding in the nearest neighbors program which states that
    (for a fixed number of neighbors, k) the number of points in the embedding
    must exceed k + 2 * orbital_lag - 1. */

    n_embed_reduced = (sint32) ( (double) max_neighbors / (double) k * (double) n_embed );
    k               = max_neighbors;
    n_embed_min     = k + 2 * orbital_lag - 1;

    if ( n_embed_reduced < n_embed_min ){

      n_embed_reduced = n_embed_min;
    }

    /* randomly select n_embed_reduced points from the original
    embedding for the sampled_embedding matrix. use this ersatz
    topology in forming the kd-tree. do this by creating a
    randomized index vector (without replacement)
    and register with the memory manager. the range of the values
    is on [0, N - 1] where N is the number of points in the emebdding.
    we will use the first n_embed_reduced of these indices to identify
    the points (rows) to copy from the original embedding matrix to the
    sampled_embedding matrix */

    /* initiate random number generation */

    err = mutil_rand_begin( &rand_ptr );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    err = mutil_rand_set_seed( (void*) NULL, rand_ptr );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    err = matuniv_random_uniform_indices(
      n_embed, 1, FALSE, rand_ptr, intrp_ptr, &random_index );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    err = matuniv_malloc_register( &sampled_embedding, n_embed_reduced, dim,
      MUTIL_DOUBLE, &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    pd_sampled = sampled_embedding.mat.dblmat.data;
    ps_index   = random_index.mat.s32mat.data;

    for ( row = 0; row < n_embed_reduced; row++ ){

      pd_embed = embedding->mat.dblmat.data + (*ps_index) * dim;

      for ( col = 0; col < dim; col++ ){

        *pd_sampled = *pd_embed;
        pd_sampled++;
        pd_embed++;
      }

      ps_index++;

    }

    /* free the random index vector and
    end random number generation */

    MUTIL_FREE_WARN( matuniv, &random_index );

    err = mutil_rand_end( rand_ptr );
    if ( err ) {
      MUTIL_ERROR( "Problem ending random number generator" );
      return err;
    }

    /* create the kd-tree with the sampled embedding */

    err = mutil_kdtree_malloc_register(
      kdtree,
      &sampled_embedding,
      bucket_size,
      list );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    /* free the sampled embedding */

    err = memlist_member_free( &sampled_embedding, &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

  }
  else{  /* requested n_neighbor does not exceed specified max_neighbors */

    err = mutil_kdtree_malloc_register(
      kdtree,
      embedding,
      bucket_size,
      list );
    MEMLIST_FREE_ON_ERROR( err, &mlist );
  }

  /* Free all remaining memory */

  MUTIL_FREE_WARN( memlist, &mlist );

  MUTIL_TRACE( "Done with localfn_kdtree()" );

  return MUTIL_ERR_OK;

}

/* Finding the proper embedding dimesion for     */
/* single variable time series ala the method    */
/* of false nearest strands                      */
/*                                               */
/* Documented in fra_dim.h                       */
/* Written by William Constantine                */

mutil_errcode frauniv_dimension_false_nearest_strands(
  const univ_mat *time_series,
  const sint32    embedding_dimension,
  const sint32    time_lag,
  const sint32    orbital_lag,
  const sint32    iterate_tolerance,
  const double    atol,
  void           *intrp_ptr,
  univ_mat       *result )
{
  boolean        is_delay_embedding;
  boolean        sort_distances = TRUE;
  double         attractor_size;
  double         mean;
  double         variance;
  double        *pd_result;
  double        *pd_series;
  fra_distance_metric distance_metric = FRA_DISTANCE_L2;
  memlist        list;
  mutil_errcode  err;
  sint32         dim;
  sint32         i;
  sint32         j;
  sint32         n_embed;
  univ_mat       embedding;
  univ_mat       neighbor_distances;
  univ_mat       neighbor_indices;
  univ_mat       original_indices;
  /* sint32         n_sample; */
  sint32_mat     strand_class;
  sint32_mat     dindex;
  sint32_mat     n_pairs;
  sint32         n_strands;
  sint32        *ps_original;
  sint32        *ps_neighbor;
  boolean        is_image;
  double         strand_statistic;
  sint32         lag;
  double         projected_distance;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_dimension_false_nearest_strands()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  err = frauniv_check_embedding_inputs(
    time_series, embedding_dimension, time_lag, orbital_lag,
    intrp_ptr, &is_delay_embedding, &n_embed );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( orbital_lag < 1 ){
    MUTIL_ERROR("Orbital lag in FNS must be at least unity");
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( atol <= 0.0 ){
    MUTIL_ERROR("FNS tolerances atol must be positve");
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( iterate_tolerance < 0 ){
    MUTIL_ERROR("Iterate tolerance in FNS must be non-negative");
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* allocate memory for FNS percentage vector and zero it out */

  err = matuniv_malloc_register( result, embedding_dimension, 1, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_assign_scalar( 0.0, intrp_ptr, &(result->mat.dblmat) );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* estimate attractor size */

  err = matuniv_mean_variance( time_series, (boolean) FALSE, intrp_ptr, &mean, &variance );
  MEMLIST_FREE_ON_ERROR( err, &list );

  attractor_size = sqrt( variance );

  /* set pointers */

  pd_result = result->mat.dblmat.data;
  pd_series = time_series->mat.dblmat.data;

  /* initialize variables */

  /* n_sample = MATUNIV_NELEM( time_series ); */

  /* calculate false nearest neighbors */

  for ( dim = 1; dim <= embedding_dimension; dim++ ){

    /*
      embed the data in the current embedding dimension and register the
      result with the memory manager
    */

    if ( memlist_member_exist( &embedding, &list ) ){
      err = memlist_member_free( &embedding, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    err = frauniv_embed( time_series, dim, time_lag, intrp_ptr, &embedding );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_register( &list, &embedding, MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, &embedding );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    n_embed = MATUNIV_NROW( &embedding );

    /* find nearest neighbor for each point in the embedding */

    if ( memlist_member_exist( &original_indices, &list ) ){
      err = memlist_member_free( &original_indices, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }
    if ( memlist_member_exist( &neighbor_indices, &list ) ){
      err = memlist_member_free( &neighbor_indices, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }
    if ( memlist_member_exist( &neighbor_distances, &list ) ){
      err = memlist_member_free( &neighbor_distances, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    err = frauniv_neighbor_find(
      &embedding,
      (sint32) 1,
      (double) 0.0,
      distance_metric,
      (univ_mat*) NULL,
      sort_distances,
      orbital_lag,
      intrp_ptr,
      &original_indices,
      &neighbor_indices,
      &neighbor_distances );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* register the nearest neighbor matrices with the memory manager */

    err = memlist_member_register( &list, &original_indices, MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, &original_indices );
      MUTIL_FREE_WARN( matuniv, &neighbor_indices );
      MUTIL_FREE_WARN( matuniv, &neighbor_distances );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }
    err = memlist_member_register( &list, &neighbor_indices, MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, &neighbor_indices );
      MUTIL_FREE_WARN( matuniv, &neighbor_distances );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }
    err = memlist_member_register( &list, &neighbor_distances, MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, &neighbor_distances );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    /* agglomerate the neighbor pairs into temporally correlated strands ... */

    /* ... allocate memory */

    if ( memlist_member_exist( &strand_class, &list ) ){
      err = memlist_member_free( &strand_class, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }
    err = mats32_malloc_register( &strand_class, n_embed, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    if ( memlist_member_exist( &n_pairs, &list ) ){
      err = memlist_member_free( &n_pairs, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }
    err = mats32_malloc_register( &n_pairs, n_embed, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    if ( memlist_member_exist( &dindex, &list ) ){
      err = memlist_member_free( &dindex, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }
    err = mats32_malloc_register( &dindex, n_embed, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* ... initialize variables */

    n_strands = 1;

    err = mats32_assign_scalar( 0, intrp_ptr, &strand_class );
    MEMLIST_FREE_ON_ERROR( err, &list );
    err = mats32_assign_scalar( 0, intrp_ptr, &n_pairs );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* assign first pair to the first strand */

    strand_class.data[0] = 1;
    n_pairs.data[0] = 1;

    /* assign pointers */

    ps_original = original_indices.mat.s32mat.data;
    ps_neighbor = neighbor_indices.mat.s32mat.data;

    /* for each nearest neighbor pair, calculate the difference
    between their respective indices */

    for ( i = 0; i < n_embed; i++ ){

      dindex.data[i] = (sint32) MUTIL_ABS( *(ps_original + i) - *(ps_neighbor + i) );
    }

    /* initiate strand statistic. lag is the number of points that we must
    skip in the time series to find the (dim + 1)th projected distance (distance
    in next dimension) for the reference point. */

    lag = dim * time_lag;

    strand_statistic = MUTIL_ABS( *(pd_series + *ps_original + lag) - *(pd_series + *ps_neighbor + lag) );

    /* form strands starting with the second point */

    for ( i = 1; i < n_embed; i++ ){

      /* increment pointers */

      ps_original++;
      ps_neighbor++;

	    /* determine if the current point is already registered
      as an element of strand pair. if so, break and move onto next point */

      if ( strand_class.data[i] < 0){
        continue;
      }

      /* determine whether the current [index,neighbor] pair is a temporal
      image of the last olag pairs registered in the current strand.
      if so, then add the current pair to that strand. otherwise,
      create a new strand with the current pair as its starting member */

      is_image = FALSE;

      for ( j = i - 1; j >= MUTIL_MAX( 0, i - orbital_lag ); j-- ){

      /* test to see if any of the preceding olag points belonging to the current
        strand are temporal iterates */

        if ( ( MUTIL_ABS( dindex.data[i] - dindex.data[j] ) <= iterate_tolerance	) &&
          strand_class.data[j] == n_strands ){

          is_image = TRUE;
          break;
        }
      }

      /* calculated the projected (dim + 1)th distance between the current nearest neighbor pair */

      projected_distance = MUTIL_ABS( *(pd_series + *ps_original + lag) - *(pd_series + *ps_neighbor + lag) );

      if ( is_image ){

        /* add pair to current strand and nullify neighbor index
        to avoid symmetric pair formation in another strand */

        strand_class.data[ *ps_neighbor ] = (sint32) -1;

        /* increment the number of pairs in current strand */

        n_pairs.data[ n_strands - 1 ] = n_pairs.data[ n_strands - 1 ] + 1;

        /* update strand statistic */

        strand_statistic += projected_distance;

      }
      else{

        /* perform the FNS test */

        /*printf("S(%ld) = %10.4f, sum(d+1)=%10.4f, n_pairs=%ld, n_strands=%ld, Ra=%10.4f, atol=%10.4f\n",
       dim, strand_statistic / (double) n_pairs.data[ n_strands - 1] / attractor_size ,
       strand_statistic, n_pairs.data[ n_strands - 1 ], n_strands, attractor_size, atol ); */

        if ( ( strand_statistic / (double) n_pairs.data[ n_strands - 1 ] / attractor_size ) > atol ){

          (*pd_result) += 1.0;
        }

        /* create a new strand */

        n_strands++;
        n_pairs.data[ n_strands - 1 ] = (sint32) 1;

        /* reset strand statistic */

        strand_statistic = projected_distance;

      }

      /* assign strand class for current (original, neighbor) pair. we need this mainly
      as a check above to ensure that pair is a candidate for analysis (if not, then a negative
      index is assigned). here, we are just assigning a positive index, which happens to be
      the current number of counted strands (equivalent to the label or assigned class of the
      current strand). */

      strand_class.data[i] = n_strands;

    } /* end strand formation */


    /* normalize the current FNS tally by the number of recorded strands and multiply
    by 100 to form FNS percentage for current dimension */

    (*pd_result) *= 100.0 / (double) n_strands;

    /* increment pointers */

    pd_result++;

    /* check for interrupts */

    if ( MUTIL_INTERRUPT( 3.0 * n_embed, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }

  } /* end dimension loop */

  /* free nodes corresponding to registered
  memory for the result, but do not free
  the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
  corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with frauniv_dimension_false_nearest_strands()" );

  return MUTIL_ERR_OK;
}









