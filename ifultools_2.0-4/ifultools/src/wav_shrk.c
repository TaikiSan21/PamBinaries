
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/wav_shrk.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "wav_shrk.h"

#include "mat_arit.h"
#include "mat_assn.h"
#include "mat_set.h"
#include "mat_sort.h"
#include "mat_stat.h"
#include "mat_summ.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "mat_io.h"

#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_limit.h"
#include "ut_math.h"
#include "ut_mem.h"

#include "wav_dwtc.h"
#include "wav_modw.h"
#include "wav_type.h"
#include <math.h>

/* Static functions declared here and defined at end of file */

static mutil_errcode localdef_sure(
  const mat_set  *W,
  const sint32    n_level,
  const double    noise_scale,
  void           *intrp_ptr,
  double_mat     *threshold );

static mutil_errcode localfn_shrink_input_check(
  const univ_mat             *time_series,
  const mat_set              *filters,
  const univ_mat             *threshold,
  const wav_shrink_threshold  threshold_function,
  const double                threshold_scale,
  const wav_shrink_function   shrink_function,
  const sint32                n_level,
  const boolean               decimated );

static mutil_errcode localfn_filters_check(
   const mat_set *filters );

static void localdef_hard_threshold( double *x, double threshold );
static void localdef_soft_threshold( double *x, double threshold );
static void localdef_mid_threshold( double  *x, double threshold );
static double localdef_sign( double x );

static mutil_errcode localdef_mad_noise(
  const univ_mat *W,
  const boolean   decimated,
  void           *intrp_ptr,
  double         *scale);

static void localdef_waveshrink(
  mat_set    *wavelet_transform,
  sint32      n_level,
  double_mat *threshold,
  void        ( *threshold_function ) ( double *, double ) );

/* Static macro definitions */

#define LOCALDEF_CHECK_NULL_POINTER_SHRK( DATA_PTR, DATA_TYPE, \
                                     TYPE_PREFIX )             \
   err = TYPE_PREFIX ## _validate( DATA_PTR );                 \
   if ( err ) return err;                                      \
   if ( DATA_PTR == ( DATA_TYPE * ) NULL ){                    \
     MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
     return MUTIL_ERR_NULL_POINTER;                            \
   }

#undef LOCALDEF_ILOG2
#define LOCALDEF_ILOG2( VALUE ) \
(sint32) floor( MUTIL_LOG2( (double) ( VALUE ) + MUTIL_DOUBLE_EPSILON ) )

#undef LOCALDEF_IS_ODD
#define LOCALDEF_IS_ODD(N) ( (N % 2) == 1 ? 1 : 0 )

#ifdef __cplusplus
extern "C" {
#endif

  mutil_errcode ( *forward_transform )( const univ_mat *, const mat_set *, sint32, void *, mat_set * );
  mutil_errcode ( *inverse_transform )( const mat_set *, const mat_set *, void *, univ_mat * );

#ifdef __cplusplus
}
#endif

/* Nonlinear denoising via wavelet shrinkage. */
/* Documented in wav_shrk.h                   */
/* Written by William Constantine             */

mutil_errcode wavuniv_shrink(
  const univ_mat             *time_series,
  const mat_set              *filters,
  const univ_mat             *threshold,
  const wav_shrink_threshold  threshold_function,
  const double                threshold_scale,
  const double                noise_variance,
  const wav_shrink_function   shrink_function,
  const sint32                n_level,
  const boolean               decimated,
  void                       *intrp_ptr,
  univ_mat                   *result )
{
  boolean        threshold_defined;
  memlist        list;
  mutil_errcode  err;
  void          ( *threshfun ) ( double *, double );
  univ_mat       threshold_values;
	double         sqrt2 = (double) sqrt(2.0);
  mat_set        wavelet_transform;
  sint32         j;
  sint32         n_sample;
  sint32         n_level_threshold;
  /* sint32         iminimax; */

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_shrink()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check input arguments */

  err = localfn_shrink_input_check( time_series, filters,
    threshold, threshold_function, threshold_scale,
    shrink_function, n_level, decimated );
  if ( err ) return err;

  /* define threshold variables */

  threshold_defined = (boolean) ( threshold != (univ_mat *) NULL );

  if ( threshold_defined ){

	  LOCALDEF_CHECK_NULL_POINTER_SHRK( threshold, univ_mat, matuniv );

	  if ( threshold->type != MUTIL_DOUBLE ){
		  MUTIL_ERROR( "Threshold matrix must be of type MUTIL_DOUBLE." );
		  return MUTIL_ERR_ILLEGAL_TYPE;
	  }

	  if ( !MATANY_IS_VEC( &(threshold->mat.dblmat) ) ){
		  MUTIL_ERROR("Threshold matrix must contain a single-column or single-row");
		  return MUTIL_ERR_ILLEGAL_SIZE;
	  }

	  n_level_threshold = MATUNIV_NELEM( threshold );

	  if ( n_level_threshold > 1 && n_level_threshold != n_level ){

		MUTIL_ERROR("The number of specified threshold levels "
			"must either be unity or equal to the number of specified "
			"wavelet transform decomposition levels" );

		return MUTIL_ERR_ILLEGAL_SIZE;

	  }

  }

  /* initialize local variables */

  n_sample = MATUNIV_NELEM( time_series );
  /* iminimax = MUTIL_MIN( LOCALDEF_ILOG2( n_sample ), 15 ); */

  /* switch on threshold function */

  switch( shrink_function ){

    case WAV_SHRINK_FUNCTION_HARD:
      threshfun = localdef_hard_threshold;
      break;
    case WAV_SHRINK_FUNCTION_SOFT:
      threshfun = localdef_soft_threshold;
      break;
    case WAV_SHRINK_FUNCTION_MID:
      threshfun = localdef_mid_threshold;
      break;
    default:
      MUTIL_ERROR( "Wavelet shrinkage function not supported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* force adaptive thresholding to use soft thresholding */

  if ( threshold_function == WAV_SHRINK_THRESHOLD_ADAPTIVE ){

    threshfun = localdef_soft_threshold;
  }

  /* choose either a decimated or undecimated DWT */

  if ( decimated ){
		forward_transform = wavuniv_transform_discrete_wavelet_convolution;
		inverse_transform = wavuniv_transform_discrete_wavelet_convolution_inverse;
  }
  else{
		forward_transform = wavuniv_transform_maximum_overlap;
		inverse_transform = wavuniv_transform_maximum_overlap_inverse;
  }

  /* calculate the (MO)DWT and register result with the memory manager */

  err = forward_transform( time_series, filters, n_level, intrp_ptr,
	  &wavelet_transform );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &wavelet_transform, MEMTYPE_MATSET );
  if ( err ){
      MUTIL_FREE_WARN( matset, &wavelet_transform );
      MEMLIST_FREE_ON_ERROR( err, &list );
  }

  if ( !threshold_defined ){

    err = wavuniv_shrink_threshold(
      &wavelet_transform,
      decimated,
      threshold_function,
      shrink_function,
      noise_variance,
      intrp_ptr,
      &threshold_values );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_register( &list, &threshold_values, MEMTYPE_MATUNIV );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* amplify the thresholds if desired */

    if ( threshold_scale != (double) 1.0 ){

      err = matdbl_multiply_scalar( &(threshold_values.mat.dblmat),
        threshold_scale, intrp_ptr, &(threshold_values.mat.dblmat));
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

  }
  else{

    /* allocate space for thresholds */

    err = matuniv_malloc_register( &threshold_values, n_level, 1,
      MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

		if ( n_level_threshold == 1 ){

			*(threshold_values.mat.dblmat.data) = *(threshold->mat.dblmat.data);

			/* adjust for MODWT if selected */

			if ( !decimated ){

				for ( j = 1; j < n_level; j++ ){

					threshold_values.mat.dblmat.data[j] =
            threshold_values.mat.dblmat.data[j - 1] / sqrt2;
				}
			}

		}
		else{

		/* all thresholds have been specified explicitly, so just
			point the threshold_values to the correct location */

			for ( j = 0; j < n_level; j++ ){

				threshold_values.mat.dblmat.data[j] = threshold->mat.dblmat.data[j];
			}
		}

  }

  /* perform the waveshrink operations */

  localdef_waveshrink( &wavelet_transform, n_level, &(threshold_values.mat.dblmat),
	  threshfun );

  /* perform an inverse (MO)DWT using the shrunken transform */

  err = inverse_transform( &wavelet_transform, filters, intrp_ptr, result );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* make sure input and output have same dimensions */

  result->mat.dblmat.nrow  = MATUNIV_NROW( time_series );
  result->mat.dblmat.ncol  = MATUNIV_NCOL( time_series );
  result->mat.dblmat.nelem = MATUNIV_NELEM( time_series );

  if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }

  /* free memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_shrink()" );

  return MUTIL_ERR_OK;
}

/***************************************/
/* STATIC FUNCTION DEFINITIONS         */
/***************************************/

/* Checks the inputs for DWT functions */
/* Written by William Constantine */

static mutil_errcode localfn_shrink_input_check(
  const univ_mat             *time_series,
  const mat_set              *filters,
  const univ_mat             *threshold,
  const wav_shrink_threshold  threshold_function,
  const double                threshold_scale,
  const wav_shrink_function   shrink_function,
  const sint32                n_level,
  const boolean               decimated )
{
  mutil_errcode  err;
  sint32         n_level_max;
  sint32         n_sample = MATUNIV_NELEM( time_series );

  MUTIL_TRACE( "Start localfn_shrink_input_check()" );

  /*** time_series ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_SHRK( time_series, univ_mat, matuniv );

  /* ... to see if it is a vector */

  if ( ( MATUNIV_NCOL( time_series ) != 1 ) &&
    ( MATUNIV_NROW( time_series ) != 1 ) ){
    MUTIL_ERROR( "Time series matrix must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if number of elements is positive */

  if ( MATUNIV_NELEM( time_series ) < 1 ){
    MUTIL_ERROR( "Number of elements in time series must be positive." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if the type is not dcomplex */

  if ( time_series->type == MUTIL_DCOMPLEX ){
    MUTIL_ERROR( "Complex time series are currently unsupported." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /*** filters ***/

  err = localfn_filters_check( filters );
  if ( err ) return err;

  /*** threshold ***/

  /* threshold CAN be a NULL pointer */

  if ( threshold != (univ_mat *) NULL ){

    if ( threshold->type != MUTIL_DOUBLE ){
      MUTIL_ERROR( "Threshold matrix must be of type MUTIL_DOUBLE." );
      return MUTIL_ERR_ILLEGAL_TYPE;
    }

    if ( !MATANY_IS_VEC( &( threshold->mat.dblmat ) ) ){
      MUTIL_ERROR( "Threshold matrix must be a single column or row." );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
  }

  /*** threshold_function ***/

  switch( threshold_function ){
    case WAV_SHRINK_THRESHOLD_UNIVERSAL:
    case WAV_SHRINK_THRESHOLD_MINIMAX:
    case WAV_SHRINK_THRESHOLD_ADAPTIVE:
      break;
    default:
       MUTIL_ERROR( "Wavelet threshold estimator not supported" );
       return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /*** threshold_scale ***/

  if ( threshold_scale <= 0.0 ){
    MUTIL_ERROR( "Threshold scale factor must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  /*** shrink_function ***/

  switch( shrink_function ){
    case WAV_SHRINK_FUNCTION_HARD:
    case WAV_SHRINK_FUNCTION_SOFT:
    case WAV_SHRINK_FUNCTION_MID:
      break;
    default:
      MUTIL_ERROR( "Wavelet shrinkage function not supported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /*** n_level ***/

  if ( n_level <= 0 ){
    MUTIL_ERROR( "Number of decomposition levels must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  if ( decimated ){

	  /* calculate maximum possible number of scales */

	  n_level_max = LOCALDEF_ILOG2( n_sample );

	  /* verify sizes */

	  if ( n_level_max < n_level ) {
		  MUTIL_ERROR( "Number of decomposition levels exceeds maximum." );
		  return MUTIL_ERR_ILLEGAL_SIZE;
	  }
  }

  MUTIL_TRACE( "Done with localfn_shrink_input_check()" );

  return MUTIL_ERR_OK;
}

/* Checks the filters for DWT functions */
/* Written by William Constantine */

static mutil_errcode localfn_filters_check(
  const mat_set *filters )
{
  mutil_errcode err;
  sint32        filter_length_scaling;
  sint32        filter_length_wavelet;

  /*** check wavelet filter ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_SHRK( filters, mat_set, matset );

  /* ... for type MUTIL_DOUBLE */

  if ( filters->mats[0].type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Wavelet filter must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( &(filters->mats[ 0 ].mat.dblmat) ) ){
    MUTIL_ERROR( "Wavelet filter matrix must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if number of elements is positive */

  filter_length_wavelet = MATUNIV_NELEM( &filters->mats[ 0 ] );

  if ( MATUNIV_NELEM( &filters->mats[ 0 ] ) < 1 ){
    MUTIL_ERROR( "Wavelet filter length must be greater than one." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if the length of the filter is even */

  if ( ( filter_length_wavelet % 2 ) != 0 ){
    MUTIL_ERROR( "Wavelet filter length must be even." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*** check scaling filter ... ***/

  /* ... for type MUTIL_DOUBLE */

  if ( filters->mats[0].type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Scaling filter must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( &(filters->mats[ 1 ].mat.dblmat) ) ){
    MUTIL_ERROR( "Scaling filter matrix must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if number of elements is positive */

  filter_length_scaling = MATUNIV_NELEM( &filters->mats[ 1 ] );

  if ( MATUNIV_NELEM( &filters->mats[ 1 ] ) < 1 ){
    MUTIL_ERROR( "Scaling filter length must be greater than one." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if the length of the filter is even */

  if ( ( filter_length_scaling % 2 ) != 0 ){
    MUTIL_ERROR( "Scaling filter length must be even." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*** check consistent filter lengths */

  if ( filter_length_scaling != filter_length_wavelet ){
    MUTIL_ERROR( "Wavelet filter and scaling filter must be the same length." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  return MUTIL_ERR_OK;
}

/** Hard thresholding for waveshrink.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_shrk.c
 * @library fractal
 * @usage #localdef_hard_threshold( &x, threshold );#
 * @return Void.
 * @param  x Pointer to a double value representing discrete wavelet
 *           transform coefficient. Upon return, this value is
 *           replaced by the thresholded value.
 * @param  threshold Double value representing the threshold.
 *
 * @see localdef_soft_threshold
 * @see localdef_mid_threshold
 * @private
 */
static void localdef_hard_threshold( double *x, double threshold ){
  *x = ( ( MUTIL_ABS( *x ) <= threshold ) ? 0.0 : *x );
}

/** Soft thresholding for waveshrink.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_shrk.c
 * @library fractal
 * @usage #localdef_soft_threshold( &x, threshold );#
 * @return Void.
 * @param  x Pointer to a double value representing discrete wavelet
 *           transform coefficient. Upon return, this value is
 *           replaced by the thresholded value.
 * @param  threshold Double value representing the threshold.
 *
 * @see localdef_hard_threshold
 * @see localdef_mid_threshold
 * @private
 */
static void localdef_soft_threshold( double *x, double threshold ){

  double fac = MUTIL_ABS( *x ) - threshold;

  *x = localdef_sign( *x ) * ( ( fac < 0.0 ) ? 0.0 : fac );
}

/** Mid thresholding for waveshrink.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_shrk.c
 * @library fractal
 * @usage #localdef_mid_threshold( &x, threshold );#
 * @return Void.
 * @param  x Pointer to a double value representing discrete wavelet
 *           transform coefficient. Upon return, this value is
 *           replaced by the thresholded value.
 * @param  threshold Double value representing the threshold.
 *
 * @see localdef_soft_threshold
 * @see localdef_hard_threshold
 * @private
 */
static void localdef_mid_threshold( double *x, double threshold ){

  double absx = MUTIL_ABS( *x );
  double fac = absx - threshold;

  *x = localdef_sign( *x ) * ( ( absx < ( 2.0 * threshold ) ) ? 2.0 * ( ( fac < 0.0 ) ? 0.0 : fac ) : absx );
}

/** Sign of a number including zero.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_shrk.c
 * @library fractal
 * @usage #y = localdef_sign( x );#
 * @return 1.0 for x > 0.0, - 1.0 for x < 0.0, 0.0 for x = 0.0.
 * @param  x Double value.
 * @private
 */
static double localdef_sign( double x )
{
  if ( x > 0.0 ){
    return 1.0;
  }
  else if ( x == 0.0 ){
    return 0.0;
  }
  else{
    return - 1.0;
  }
}

/** Mean absolute deviation estimate of noise scale (stdev) using
 * discrete wavelet coefficients and with Gaussian normalization.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_shrk.c
 * @library fractal
 * @usage #err = localdef_mad_noise( &W, decimated, intrp_ptr, &mad );#
 * @return Standard mutils error/OK code.
 * @param  W Universal matrix containing the current scale's
 *           wavelet coefficients.
 * @param decimated Boolean. TRUE means that a decimated DWT is being used.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param mad         Pointer to a double value containing the result.
 * @private
 */
static mutil_errcode localdef_mad_noise(
  const univ_mat *W,
  const boolean   decimated,
  void           *intrp_ptr,
  double         *scale )
{
  univ_mat       absW;
  memlist        list;
  mutil_errcode  err;
  sint32         n_sample = MATUNIV_NELEM( W );

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start localdef_mad_noise()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* allocate memory */

  err = matuniv_malloc_register( &absW, MATUNIV_NROW( W ), MATUNIV_NCOL( W ),
	  MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* take absolute value of wavelet coefficients and store */

  err = matuniv_abs( W, intrp_ptr, &absW );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate median */

  err = matuniv_median( &absW, intrp_ptr, scale );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* normalize for assumed Gaussian distribution */

  *scale /= 0.6745;

  if ( !decimated ) *scale *= sqrt(2.0);

  /* free memory list */

  MUTIL_FREE_WARN( memlist, &list );

  if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE( "Done with localdef_mad_noise()" );

  return MUTIL_ERR_OK;
}

/** Waveshrink.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_shrk.c
 * @library fractal
 * @usage #localdef_waveshrink( &wavelet_transform, n_level, &threshold, &threshold_function );#
 * @return Void.
 * @param wavelet_transform Pointer to a matrix set containing the original (MO)DWT.
 *   The wavelet shrinkage is done in
 *   place so that upon return W will contain the shrunken (MO)DWT.
 * @param n_level The number of decomposition levels (to shrink).
 * @param A double matrix containing the thresholds, one per decomposition level.
 * @private
 */
static void localdef_waveshrink(
  mat_set    *wavelet_transform,
  sint32      n_level,
  double_mat *threshold,
  void        ( *threshold_function ) ( double *, double ) )
{
  sint32  j;
  sint32  t;
  double *pd_coeff;
  double  delta;

  for ( j = 0; j < n_level; j++ ){

    /* assign pointer */

    pd_coeff = wavelet_transform->mats[ j ].mat.dblmat.data;

    /* obtain threshold */

    delta = threshold->data[ j ];

    /* threshold the wavelet coefficients */

    for ( t = 0; t < wavelet_transform->mats[ j ].mat.dblmat.nelem; t++ ){

      threshold_function( pd_coeff, delta );

      /* advance pointer */

      pd_coeff++;
    }
  }
}

/** Waveshrink threshold based on Stein's unbiased risk estimator (SURE).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_shrk.c
 * @library fractal
 * @usage #err = localdef_sure( &W, noise_scale, intrp_ptr, threshold );#
 * @return Void.
 * @param W Matrix set containing the original DWT. The wavelet shrinkage is done in
 *          place so that upon return W will contain the shrunken DWT.
 * @param n_level The number of decomposition levels (to shrink).
 * @param noise_scale The estimated noise scale.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param A double matrix containing the thresholds, one per decomposition level.
 *
 * @see localdefmad_noise
 * @private
 */
static mutil_errcode localdef_sure(
  const mat_set  *W,
  const sint32    n_level,
  const double    noise_scale,
  void           *intrp_ptr,
  double_mat     *threshold )
{
  double         cumsum;
  double         risk;
  double         risk_min;
  double         s;
  double        *pd_sxi;
  memlist        list;
  mutil_errcode  err;
  sint32         Nj;
  sint32         irisk_min;
  sint32         j;
  sint32         ncol;
  sint32         nrow;
  sint32         t;
  univ_mat       Wsquared;
  univ_mat       sxi;
  univ_mat      *pd_W;
  double        *pd_thresh;
  double         mean;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start localdef_sure()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  pd_thresh = threshold->data;

  /* calculate the SURE threshold for each decomposition level */

  for ( j = 0; j < n_level; j++ ){

    pd_W = &( W->mats[ j ] );

    ncol = MATUNIV_NCOL( pd_W );
    nrow = MATUNIV_NROW( pd_W );
    Nj   = MATUNIV_NELEM( pd_W );

    /* allocate memory */

    err = matuniv_malloc_register( &Wsquared, nrow, ncol, MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matuniv_malloc_register( &sxi, nrow, ncol, MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* form cumulative sum of squares of wavelet coefficients vector */

    err = matuniv_multiply_elem( pd_W, pd_W, intrp_ptr, &Wsquared );
    MEMLIST_FREE_ON_ERROR( err, &list );


    /* if sparse, use a universal threshold at the current level. otherwise
       use the SURE estimator */

    err = matdbl_sum( &( Wsquared.mat.dblmat ), intrp_ptr, &mean );
    MEMLIST_FREE_ON_ERROR( err, &list );

    mean /= (double) Nj;

    if ( mean - 1.0 <= MUTIL_POW( MUTIL_LOG2( (double) Nj ), 1.5 ) / sqrt( (double) Nj ) ){

      *pd_thresh = sqrt( 2.0 * log( (double) Nj ) ) * noise_scale;

    }
    else{

      err = matdbl_divide_scalar( &( Wsquared.mat.dblmat ), noise_scale * noise_scale,
	TRUE, intrp_ptr, &( Wsquared.mat.dblmat ) );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = matuniv_sort( &Wsquared, intrp_ptr, &sxi );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* form risk elements and record the location of the minimum risk */

      pd_sxi = sxi.mat.dblmat.data;

      cumsum = *pd_sxi;

      for ( t = 0; t < Nj; t++ ){

	s = cumsum + (double) ( Nj - 1 - t ) * *pd_sxi;

	risk = ( (double) ( Nj - 2 * ( t + 1 ) ) + s ) / (double) Nj;

	/* increment pointers */

	pd_sxi++;

	/* add to cumsum */

	cumsum += *pd_sxi;

	if ( t == 0 ){

	  risk_min  = risk;
	  irisk_min = 0;
	}
	else{

	  if ( risk < risk_min ){

	    risk_min  = risk;
	    irisk_min = t;
	  }
	}
      }

      *pd_thresh = sqrt( sxi.mat.dblmat.data[ irisk_min ] ) * noise_scale;
    }

    pd_thresh++;

    /* free memory list */

    MUTIL_FREE_WARN( memlist, &list );

    if ( MUTIL_INTERRUPT( 3.0 * Nj, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }
  }

  MUTIL_TRACE( "Done with localdef_sure()" );

  return MUTIL_ERR_OK;
}

mutil_errcode wavuniv_shrink_threshold(
  const mat_set              *wavelet_transform,
  const boolean               decimated,
  const wav_shrink_threshold  threshold_function,
  const wav_shrink_function   shrink_function,
  const double                noise_variance,
  void                       *intrp_ptr,
  univ_mat                   *result )
{
  sint32         iminimax;
  const double   minimax_soft[] = { 0.0, 0.0, 0.0, 1.20, 1.27, 1.474, 1.669, 1.859, 2.045,
    2.226, 2.403, 2.575, 2.743, 2.906, 3.066, 3.221 };
  const double   minimax_hard[] = { 0.0, 0.0, 0.0, 2.23, 2.47, 2.697, 2.913, 3.117, 3.312,
    3.497, 3.674, 3.844, 4.008, 4.166, 4.319, 4.467 };
  double         noise_scale;
  double         delta;
  sint32         n_sample;
  sint32         n_level;
  sint32         n_sample_j;
  boolean        extra_counted;
  sint32         j;
  memlist        list;
  mutil_errcode  err;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_shrink()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* perform argument checks */

  LOCALDEF_CHECK_NULL_POINTER_SHRK( wavelet_transform, mat_set, matset );

  /* ... for type MUTIL_DOUBLE */

  if ( wavelet_transform->mats[ 0 ].type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Transform matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /*** check shrinkage function ***/

  switch( shrink_function ){
    case WAV_SHRINK_FUNCTION_HARD:
    case WAV_SHRINK_FUNCTION_SOFT:
    case WAV_SHRINK_FUNCTION_MID:
      break;
    default:
      MUTIL_ERROR( "Wavelet shrinkage function not supported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /*** check threshold estimator ***/

  switch( threshold_function ){
    case WAV_SHRINK_THRESHOLD_UNIVERSAL:
    case WAV_SHRINK_THRESHOLD_MINIMAX:
    case WAV_SHRINK_THRESHOLD_ADAPTIVE:
      break;
    default:
       MUTIL_ERROR( "Wavelet threshold estimator not supported" );
       return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  if ( !decimated && threshold_function != WAV_SHRINK_THRESHOLD_UNIVERSAL ){
    MUTIL_ERROR("Only universal thresholding is currently supported for MODWT-based waveshrink" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* calculate the number of decomposition levels:
  for the DWT, this is slightly tricky because of the storage system
  for extra scaling coefficients. */

  if ( decimated ){

    /* calculate original sample size */

    for ( n_sample = 0, j = 0; j < wavelet_transform->nelem; j++ ){

      n_sample += MATUNIV_NELEM( &wavelet_transform->mats[ j ] );
    }

    n_sample_j = n_sample;

    extra_counted = (boolean) FALSE;

    for ( j = 1; j < wavelet_transform->nelem; j++ ){

      if ( j == wavelet_transform->nelem - 1 && !extra_counted ) break;

      if ( LOCALDEF_IS_ODD( n_sample_j ) && !extra_counted ){

        extra_counted = (boolean) TRUE;
        break;
      }

      n_sample_j /= 2; /* integer division intended */
    }

    n_level = extra_counted ? wavelet_transform->nelem - 2 : wavelet_transform->nelem - 1;
  }
  else{

    n_sample = MATUNIV_NELEM( &(wavelet_transform->mats[0]) );
    n_level  = wavelet_transform->nelem - 1;

  }

  iminimax = MUTIL_MIN( LOCALDEF_ILOG2( n_sample ), 15 );

  /* estimate noise scale using level 1 wavelet coefficients */

  if ( noise_variance <= (double) 0.0 ){

    err = localdef_mad_noise( wavelet_transform->mats, decimated, intrp_ptr,
      &noise_scale );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }
  else{

    noise_scale = sqrt( noise_variance );
  }

  /* allocate memory */

  err = matuniv_malloc_register( result, n_level, 1, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* form threshold vector (one threshold per decomposition level) */

  switch( threshold_function ){

  case WAV_SHRINK_THRESHOLD_UNIVERSAL:

    /* form universal threshold */

    delta = sqrt( 2.0 * noise_scale * noise_scale * log( (double) n_sample ) );

    /* fill the threshold vector with a constant threshold */

    err = matdbl_assign_scalar( delta, intrp_ptr, &(result->mat.dblmat) );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* adjust for MODWT if selected */

    if ( !decimated ){

      for ( j = 1; j <= n_level; j++ ){

        result->mat.dblmat.data[j - 1] /=
          (double) MUTIL_POW( 2.0, (double) j / 2.0 );
      }
    }

    break;

  case WAV_SHRINK_THRESHOLD_MINIMAX:

    /* form minimax threshold */

    switch( shrink_function ){

    case WAV_SHRINK_FUNCTION_HARD:
      delta = minimax_hard[ iminimax ] * noise_scale;
      break;
    case WAV_SHRINK_FUNCTION_SOFT:
      delta = minimax_soft[ iminimax ] * noise_scale;
      break;
    default:
      MUTIL_ERROR( "Minimax threshold is only available for hard and soft thresholding functions." );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
    }

    /* fill the threshold vector with a constant threshold */

    err = matdbl_assign_scalar( delta, intrp_ptr, &(result->mat.dblmat) );
    MEMLIST_FREE_ON_ERROR( err, &list );

    break;

    case WAV_SHRINK_THRESHOLD_ADAPTIVE:

      /* form SURE thresholds, one per decomposition level */

      err = localdef_sure( wavelet_transform, n_level, noise_scale,
        intrp_ptr, &(result->mat.dblmat) );
      MEMLIST_FREE_ON_ERROR( err, &list );
      break;

    default:
      MUTIL_ERROR( "Wavelet threshold estimator not supported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }

  /* free node corresponding to registered
     output memory, but do not free the
     memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_shrink_threshold()" );

  return MUTIL_ERR_OK;
}
