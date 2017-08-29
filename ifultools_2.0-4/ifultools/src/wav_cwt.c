
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/wav_cwt.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "wav_cwt.h"
#include "wav_type.h"

#include "mat_assn.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "mth_mac.h"
#include "mth_dist.h"
#include "sig_tran.h"
#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_math.h"
#include "ut_mem.h"
#include <math.h>
#include <stdio.h>

/* Static functions declared here and defined at end of file */

static mutil_errcode localfn_cwt(
  const wav_filter_type  filter_type,
  const double           filter_arg,
  const double           scale,
  const double           sampling_interval,
  const dcomplex_mat    *dft_time_series,
  void                  *intrp_ptr,
  dcomplex_mat          *result );

static mutil_errcode localfn_wavelet_Euler_Maclaurin_expansion(
  const wav_filter_type  filter_type,
  const double           filter_arg,
  const double           tau,
  const sint32           K,
  void                  *intrp_ptr,
  dcomplex_mat          *result );

/* Static macro definitions */

#define LOCALDEF_CHECK_NULL_POINTER_CWT( DATA_PTR, DATA_TYPE,  \
                                     TYPE_PREFIX )             \
   err = TYPE_PREFIX ## _validate( DATA_PTR );                 \
   if ( err ) return err;                                      \
   if ( DATA_PTR == ( DATA_TYPE * ) NULL ){                    \
     MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
     return MUTIL_ERR_NULL_POINTER;                            \
   }

#define SQRT2 1.41421356237309510488016887242096980785697
#define TWOPI 6.2831853071795862

/* The continuous wavelet transform function  */
/* Documented in wav_cwt.h                    */
/* Written by William Constantine             */

mutil_errcode wavuniv_transform_continuous_wavelet(
  const univ_mat        *time_series,
  const double           sampling_interval,
  const wav_filter_type  filter_type,
  const double           filter_arg,
  const univ_mat        *scale,
  void                  *intrp_ptr,
  univ_mat              *cwt )
{
  memlist       list;
  mutil_errcode err;
  sint32        iscale;
  sint32        n_sample;
  sint32        n_scale;
  univ_mat      dftproduct;
  univ_mat      dftx;
  univ_mat      temp;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_transform_continuous_wavelet()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  /*** check time series ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_CWT( time_series, univ_mat, matuniv );

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

  /*** check sampling interval ... ***/

  if ( sampling_interval <= 0.0 ){
    MUTIL_ERROR( "Sampling interval of time series must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*** check filter type ... ***/

  /** To avoid redundancy and to automatically
      coordinate with any changes in wavuniv_filters_continuous()
      the checking of the filter arguments is performed
      in wavuniv_filters_continuous() **/

  /*** check scale vector ... ***/

  if ( ( MATUNIV_NCOL( scale ) != 1 ) &&
    ( MATUNIV_NROW( scale ) != 1 ) ){
    MUTIL_ERROR( "scales matrix must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if number of elements is positive */

  if ( MATUNIV_NELEM( scale ) < 1 ){
    MUTIL_ERROR( "Number of elements in scales must be positive." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if the type is MUTIL_DOUBLE */

  if ( scale->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Scale vector must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* assign dimensions to variables */

  n_sample = MATUNIV_NELEM( time_series );
  n_scale  =  MATUNIV_NELEM( scale );

  /* allocate memory */

  err = matuniv_malloc_register( cwt, n_sample, n_scale, MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &dftx, n_sample, 1, MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &dftproduct, n_sample, 1, MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* form Fourier transform of input sequence.

  NOTE: The MUTILS DFT function only takes data
  stored in columns. so trick it since we know
  we have a vector . */

  temp.type = MUTIL_DOUBLE;
  temp.mat.dblmat.nelem = n_sample;
  temp.mat.dblmat.ncol  = 1;
  temp.mat.dblmat.nrow  = n_sample;
  temp.mat.dblmat.data  = time_series->mat.dblmat.data;

  err = siguniv_transform_discrete_fourier(
    &temp,
    (boolean) FALSE,
    intrp_ptr,
    &dftx );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate the CWT */

  for ( iscale = 0; iscale < n_scale; iscale++ ){

    /* calculate the CWT coefficients at the current scale */

    err = localfn_cwt( filter_type, filter_arg,
      scale->mat.dblmat.data[ iscale ], sampling_interval,
      &( dftx.mat.cpxmat ), intrp_ptr, &( dftproduct.mat.cpxmat ) );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* copy inverse DFT into cwt matrix */

    err = matuniv_assign_submat( &dftproduct, 0, iscale, intrp_ptr, cwt );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* Check for interrupts */

    if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }

  }

  /* free nodes corresponding to output
     in memory list, but not the memory itself as
     it must remain to send back to the caller */

  err = memlist_member_unregister( cwt, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_transform_continuous_wavelet()" );

  return MUTIL_ERR_OK;
}

/** Function to calculate the CWT via the DFT at the given scale.
 * Using a Riemann sum approximation over a grid of appropriate
 * Fourier frequencies, this function calculates the frequency
 * response fo the wavelet filter at the specified scale modulated
 * by the DFT of the original time series. It then performs an inverse
 * DFT of the product to form an approximation to the CWT at the
 * specified scale.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_cwt.h
 * @source wav\_cwt.h
 * @library wavelet
 * @param  filter_type Filter type. Only the wav\_filter\_type
 *                     WAV\_FILTER\_HAAR, WAV\_FILTER\_GAUSSIAN\_I,
 *                     WAV\_FILTER\_GAUSSIAN\_II, and WAV\_FILTER\_MORLET
 *                     are supported.
 * @param  filter_arg  A double value representing a secondary argument to be
 *                     passed to the filter function. If the filter
 *                     is of type WAV\_FILTER\_GAUSSIAN\_I or
 *                     WAV\_FILTER\_GAUSSIAN\_II then this parameter
 *                     represents the standard deviation $\sigma$ of a Guassian PDF.
 *                     If the filter is of type WAV\_FILTER\_MORLET, then
 * @param  scale       The current CWT scale.
 * @param  sampling_interval The sampling interval of the time series.
 * @param  dft_time_series Pointer to a pre-allocated complex matrix
 *                     containing the discrete Fourier transform of
 *                     the time series. This matrix must be a vector, i.e.,
 *                     contain a single row or column.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  result      Pointer to a pre-allocated complex matrix containing a single-row
 *                     or single-column which (upon return) will contain the CWT
 *                     estimation at the current scale. This vector must have the
 *                     same number of elements as does the dft\_time\_series
 *                     input argument.
 * @private
 */
static mutil_errcode localfn_cwt(
  const wav_filter_type  filter_type,
  const double           filter_arg,
  const double           scale,
  const double           sampling_interval,
  const dcomplex_mat    *dft_time_series,
  void                  *intrp_ptr,
  dcomplex_mat          *result )
{
  dcomplex_mat   dft_wavelet;
  dcomplex_mat   filtered_series;
  double         tau;
  memlist        list;
  mutil_errcode  err;
  sint32         j;
  sint32         n_sample;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start localfn_cwt()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check arguments:

  The filter_type, filter_arg, scale, and sampling_interval
  arguments are checked in the main function and are constant
  so they need not be checked here.        */

  /* dft_time_series ... */

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_CWT( dft_time_series, dcomplex_mat, matcpx );

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( dft_time_series ) ){
    MUTIL_ERROR( "DFT of time series matrix must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if number of elements is positive */

  if ( dft_time_series->nelem < 1 ){
    MUTIL_ERROR( "Number of elements in DFT of time series must be positive." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* result ... */

  LOCALDEF_CHECK_NULL_POINTER_CWT( result, dcomplex_mat, matcpx );

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( result ) ){
    MUTIL_ERROR( "DFT of time series matrix must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if number of elements is the same as that of dft_time_series */

  if ( result->nelem != dft_time_series->nelem ){
    MUTIL_ERROR( "Number of elements in result must be equal " \
      "to that of the DFT of the time series matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* initialize variables */

  n_sample = dft_time_series->nelem;
  tau      = scale / sampling_interval;

  /* allocate memory for matrices */

  err = matcpx_malloc_register( &dft_wavelet, n_sample, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matcpx_malloc_register( &filtered_series, n_sample, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate scaled frequency response of the mother wavelet */

  err = localfn_wavelet_Euler_Maclaurin_expansion(
    filter_type,
    filter_arg,
    tau,
    (sint32) 5,
    intrp_ptr,
    &dft_wavelet );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* form the product of the scaled wavelet DFT approximation
     and the DFT of the time series */

  for ( j = 0; j < n_sample; j++ ){

    MUTIL_CPX_MULT( dft_time_series->data[j], dft_wavelet.data[j], filtered_series.data[j] );

    /* Check for interrupts */

    if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }
  }

  /* calculate the inverse DFT of the filtered time series */

  err = sigcpx_transform_discrete_fourier(
    &filtered_series,
    (boolean) TRUE,
    intrp_ptr,
    result );
  MEMLIST_FREE_ON_ERROR( err, &list );

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_cwt()" );

  return MUTIL_ERR_OK;
}

/** Function to estimate the DFT of a continuous scaled wavelet filter.
 * An Euler-Maclaurin series expansion is used to estimate the DFT of
 * of a scaled continuous wavelet.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_cwt.h
 * @source wav\_cwt.h
 * @library wavelet
 * @param  filter_type Filter type. Only the wav\_filter\_type
 *                     WAV\_FILTER\_HAAR, WAV\_FILTER\_GAUSSIAN\_I,
 *                     WAV\_FILTER\_GAUSSIAN\_II, and WAV\_FILTER\_MORLET
 *                     are supported.
 * @param  filter_arg  A double value representing a secondary argument to be
 *                     passed to the filter function. If the filter
 *                     is of type WAV\_FILTER\_GAUSSIAN\_I or
 *                     WAV\_FILTER\_GAUSSIAN\_II then this parameter
 *                     represents the standard deviation $\sigma$ of a Guassian PDF.
 *                     If the filter is of type WAV\_FILTER\_MORLET, then
 * @param  tau         Defined as the scale divided by the sampling interval of the
 *                     original time series.
 * @param  K           Summation index bound for the Maclaurin series approximation.
 *                     K = 5 is usually sufficient.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  result      Pointer to a pre-allocated complex matrix containing a single-row
 *                     or single-column which (upon return) will contain the result.
 * @private
 */
static mutil_errcode localfn_wavelet_Euler_Maclaurin_expansion(
  const wav_filter_type  filter_type,
  const double           filter_arg,
  const double           tau,
  const sint32           K,
  void                  *intrp_ptr,
  dcomplex_mat          *result )
{
  double         cjK        ;
  double         fj         ;
  double         fj2        ;
  double         fjK        ;
  double         fjK2       ;
  double         fjK4       ;
  double         hots       ;
  double         sum        ;

  double         Csigma     ;
  double         b          ;
  double         b3         ;
  double         a          ;
  double         aj         ;
  double         cjplus     ;
  double         cjminus    ;
  double         cjplus2    ;
  double         cjminus2   ;
  double         ajbk       ;
  double         expcjplus  ;
  double         expcjminus ;

  sint32         j;
  sint32         k;
  sint32         n_sample;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start localfn_wavelet_Euler_Maclaurin_expansion()" );

  n_sample = result->nelem;

  switch( filter_type ){

    case WAV_FILTER_GAUSSIAN_I:

      Csigma = 2 * sqrt( filter_arg ) * MUTIL_POW( MUTIL_PI, 0.25 );
      b      = TWOPI * tau * filter_arg;
      a      = b / (double) n_sample;
      b3     = b * b * b;

      for ( j = 0; j < n_sample; j++ ){

	aj         = a * (double) j;
	cjplus     = aj + b * K;
	cjminus    = aj - b * K;
	cjplus2    = cjplus * cjplus;
	cjminus2   = cjminus * cjminus;
	expcjplus  = exp( - cjplus2 / 2.0 );
	expcjminus = exp( - cjminus2 / 2.0 );

	sum = 0.0;

	for ( k = -K; k <= K; k++ ){

	  ajbk  = aj + b * (double) k;
	  sum  += ajbk * exp( - ajbk * ajbk / 2.0  );
	}

	hots = expcjplus / 2.0 * (
	  2.0 / b - cjplus - b * ( 1.0 - cjplus2 ) / 6.0 + b3 * expcjplus * (
	    cjplus2 * cjplus2 + 5.0 * cjplus2 - 3.0 ) / 360.0 )
	  - expcjminus / 2.0 * (
	    2.0 / b + cjminus - b * ( 1.0 - cjminus2 ) / 6.0 - b3 * expcjminus * (
	      cjminus2 * cjminus2 + 5.0 * cjminus2 - 3.0 ) / 360.0 );

	result->data[ j ].re = 0.0;
	result->data[ j ].im = Csigma * sqrt( tau ) * ( sum + hots );

	/* Check for interrupts */

	if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
	  MUTIL_ERROR( "user interrupt" );
	  return MUTIL_ERR_INTERRUPT;
	}
      }

      break;

    case WAV_FILTER_GAUSSIAN_II:

      Csigma = sqrt( 8.0 / 3.0 ) * sqrt( filter_arg ) * MUTIL_POW( MUTIL_PI, 0.25 );
      b      = TWOPI * tau * filter_arg;
      a      = b / (double) n_sample;
      b3     = b * b * b;

      for ( j = 0; j < n_sample; j++ ){

	aj         = a * (double) j;
	cjplus     = aj + b * K;
	cjminus    = aj - b * K;
	cjplus2    = cjplus * cjplus;
	cjminus2   = cjminus * cjminus;
	expcjplus  = exp( - cjplus2 / 2.0 );
	expcjminus = exp( - cjminus2 / 2.0 );

	sum = 0.0;

	for ( k = -K; k <= K; k++ ){

	  ajbk  = aj + b * (double) k;
	  sum  += ajbk * ajbk * exp( - ajbk * ajbk / 2.0  );
	}

	hots = sqrt( TWOPI ) / b * ( 1.0 + MUTIL_PNORM( cjminus ) - MUTIL_PNORM( cjplus ) )
	  + expcjplus / 2.0 * (
	  2.0 / b - cjplus2 - b * cjplus * ( 2.0 - cjplus2 ) / 6.0 - b3 * cjplus * expcjplus * (
	    12.0 - 9.0 * cjplus2 + cjplus2 * cjplus2 ) / 360.0 )
	  - expcjminus / 2.0 * (
	    2.0 / b - cjminus2 - b * cjminus * ( 2.0 - cjminus2 ) / 6.0 + b3 * cjminus * expcjminus * (
	      12.0 - 9.0 * cjminus2 + cjminus2 * cjminus2 ) / 360.0 );

	result->data[ j ].re = Csigma * sqrt( tau ) * ( sum + hots );
	result->data[ j ].im = 0.0;

	/* Check for interrupts */

	if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
	  MUTIL_ERROR( "user interrupt" );
	  return MUTIL_ERR_INTERRUPT;
	}
      }

      break;

    case WAV_FILTER_MORLET:

      Csigma = MUTIL_POW( MUTIL_PI, - 1.0 / 4.0 ) /
	sqrt( 1.0 - 4.0 / sqrt( 3.0 ) *
	  exp( - filter_arg * filter_arg / 4.0 ) + sqrt( 2.0 ) *
	  exp( - filter_arg * filter_arg / 2.0 ) ) ;

      b      = TWOPI * tau;
      a      = b / (double) n_sample;
      b3     = b * b * b;

      for ( j = 0; j < n_sample; j++ ){

	aj         = a * (double) j + filter_arg;
	cjplus     = aj + b * K;
	cjminus    = aj - b * K;
	cjplus2    = cjplus * cjplus;
	cjminus2   = cjminus * cjminus;
	expcjplus  = exp( - cjplus2 / 2.0 );
	expcjminus = exp( - cjminus2 / 2.0 );

	sum = 0.0;

	for ( k = -K; k <= K; k++ ){

	  ajbk  = aj + b * (double) k;
	  sum  += exp( - ajbk * ajbk / 2.0 );
	}

	hots = sqrt( TWOPI ) / b * ( 1.0 + MUTIL_PNORM( cjminus ) - MUTIL_PNORM( cjplus ) )
	  - expcjplus / 2.0 * ( b3 * cjplus * ( cjplus2 - 3.0 ) / 360.0 - b * cjplus / 6.0 + 1.0 )
	  + expcjminus / 2.0 * ( b3 * cjminus * ( cjminus2 - 3.0 ) / 360.0 - b * cjminus / 6.0 - 1.0 );

	result->data[ j ].re = Csigma * sqrt( tau * TWOPI ) * ( sum + hots );
	result->data[ j ].im = 0.0;
      }

      break;

    case WAV_FILTER_HAAR:

      /* calculate weighted frequency response */

      result->data[ 0 ].re = 0.0;
      result->data[ 0 ].im = 0.0;

      for ( j = 1; j < n_sample; j++ ){

	fj   = (double) j / (double) n_sample;
	fj2  = fj * fj;
	fjK  = fj2 - (double) ( K * K );
	fjK2 = fjK * fjK;
	fjK4 = fjK2 * fjK2;
	cjK  = fj2 + (double) ( K * K );

	hots = 1.0 / ( 2.0 * fj ) * log( MUTIL_ABS( ( (double) K - fj ) / ( (double) K + fj ) ) )
	  - 1.0 / ( 2.0 * fjK ) - (double) K / ( 6.0 * fjK2 ) + (double) K * cjK / ( 15.0 * fjK4 );

	sum = 0.0;

	for ( k = 1; k <= K; k++ ){

	  sum += 1.0 / ( fj2 - (double) ( k * k ) );
	}

	result->data[ j ].re = 0.0;
	result->data[ j ].im = SQRT2 * MUTIL_POW( sin(  tau * MUTIL_PI * fj ), 2.0 )
	  / MUTIL_PI / sqrt( tau ) * ( 1.0 / fj  + 2.0 * fj * ( sum + hots ) );

	/* Check for interrupts */

	if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
	  MUTIL_ERROR( "user interrupt" );
	  return MUTIL_ERR_INTERRUPT;
	}
      }

      break;

    default:
      MUTIL_ERROR( "Filter type is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  MUTIL_TRACE( "Done with localfn_wavelet_Euler_Maclaurin_expansion()" );

  return MUTIL_ERR_OK;
}
