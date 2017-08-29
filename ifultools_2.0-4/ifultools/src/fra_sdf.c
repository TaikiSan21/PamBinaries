
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_sdf.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */

#include "fra_sdf.h"
#include "fra_type.h"

#include "mat_arit.h"
#include "mat_assn.h"
#include "mat_stat.h"
#include "mat_summ.h"
#include "mat_type.h"
#include "mat_univ.h"
#include "mat_umat.h"

#include "sig_tran.h"

#include "str_type.h"

#include "ut_math.h"
#include "ut_mem.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_mem.h"
#include "mat_io.h"

#include <math.h>
#include <stdio.h>

/*
  This file contains function definitions for
  inferring deterministic structure in time series.
  The functions are declared in fra_sdf.h
*/

/* Static macro definitions */

#undef LOCALDEF_CHECK_NULL_POINTER_SDF
#define LOCALDEF_CHECK_NULL_POINTER_SDF( DATA_PTR,          \
  DATA_TYPE, TYPE_PREFIX )                                  \
  err = TYPE_PREFIX ## _validate( DATA_PTR );               \
  if ( err ) return err;                                    \
  if ( DATA_PTR == ( DATA_TYPE * ) NULL ) {                 \
  MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
  return MUTIL_ERR_NULL_POINTER;                            \
  }

#undef LOCALDEF_COMPLEX_CONJUGATE_PRODUCT
#define LOCALDEF_COMPLEX_CONJUGATE_PRODUCT(a, b, c) \
 (c).re = (a).re * (b).re + (a).im * (b).im; \
 (c).im = (a).re * (b).im - (a).im * (b).re


/* static function declarations */

static mutil_errcode localfn_cross_sdf(
  const univ_mat *mat,
  const sint32    n_frequency,
  void           *intrp_ptr,
  univ_mat       *result );

static mutil_errcode localfn_column_means(
  const double_mat *matrix,
  const boolean     center,
  void             *intrp_ptr,
  double_mat       *mean );

static mutil_errcode localfn_add_to_complex_matrix(
  const univ_mat *x,
  univ_mat       *y );

static mutil_errcode localfn_complex_matrix_divide_by_scalar(
  const univ_mat *x,
  double          norm );

static mutil_errcode localfn_recenter_columns(
  double_mat *matrix,
  void       *intrp_ptr );

/* Direct SDF estimator           */
/*                                */
/* Documented in fra_sdf.h        */
/* Written by William Constantine */

mutil_errcode frauniv_spectral_density_function_direct(
  const univ_mat *time_series,
  const univ_mat *taper,
  const boolean   center,
  const boolean   recenter,
  const boolean   single_sided,
  const sint32    npad,
  void           *intrp_ptr,
  univ_mat       *result )
{
  double        *pd_series;
  double        *pd_taper;
  double        *pd_hx;
  memlist        list;
  mutil_errcode  err;
  sint32         n_frequency;
  sint32         n_series;
  sint32         j;
  sint32         t;
  univ_mat       tapered_series;
  sint32         nsdf;
  sint32         n_sample;
  double_mat     mean;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_spectral_density_function_direct()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  /* ... time series */

  LOCALDEF_CHECK_NULL_POINTER_SDF( time_series, univ_mat, matuniv );

  if ( time_series->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Time series matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... taper */

  LOCALDEF_CHECK_NULL_POINTER_SDF( taper, univ_mat, matuniv );

  if ( !MATANY_IS_VEC( &(taper->mat.dblmat) ) ){
    MUTIL_ERROR( "Taper matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( taper->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Taper matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( MATUNIV_NROW( time_series ) != MATUNIV_NELEM( taper ) ){
    MUTIL_ERROR( "Time series and taper matrices must have the same number of elements." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... npad */

  if ( npad < MATUNIV_NROW( time_series ) ){
    MUTIL_ERROR( "The input npad must exceed or equal the number of samples in the time series." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* initialize variables */

  n_sample    = MATUNIV_NROW( time_series );
  n_series    = MATUNIV_NCOL( time_series );
  n_frequency = single_sided ? npad / 2 + 1 : npad;
  nsdf       = n_series * ( n_series + 1 ) / 2;

  /* allocate memory */

  err = matuniv_malloc_register( &tapered_series, npad, n_series,
		MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate mean of time series if center option is true */

  err = matdbl_malloc_register( &mean, 1, n_series, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = localfn_column_means( &( time_series->mat.dblmat ), center,
    intrp_ptr, &mean );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* taper the time series, then zero pad to n_frequency */

  for ( j = 0; j < n_series; j++ ){

    pd_taper  = taper->mat.dblmat.data;
    pd_series = time_series->mat.dblmat.data + j;
    pd_hx     = tapered_series.mat.dblmat.data + j;

    for ( t = 0; t < npad; t++ ){

      if ( t < n_sample ){

	*pd_hx = ( *pd_series - mean.data[j] ) * (*pd_taper);

	pd_taper++;
	pd_series += n_series;
      }
      else{

	*pd_hx = 0.0;
      }

      pd_hx += n_series;
    }
  }

  /* after tapering, the mean may not be zero even if
  the data was originally centered. recentered if requested */

  if (recenter){
    err = localfn_recenter_columns( &(tapered_series.mat.dblmat), intrp_ptr );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  if ( MUTIL_INTERRUPT( 10.0 * (double) n_frequency, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }

  /* calculate the direct SDF estimator */

  err = localfn_cross_sdf( &tapered_series, n_frequency, intrp_ptr, result );
  MEMLIST_FREE_ON_ERROR( err, &list );

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with frauniv_spectral_density_function_direct()" );

  return MUTIL_ERR_OK;
}

/* Lag window SDF estimator       */
/*                                */
/* Documented in fra_sdf.h        */
/* Written by William Constantine */

mutil_errcode frauniv_spectral_density_function_lag_window(
  const univ_mat *time_series,
  const univ_mat *lag_window,
  const univ_mat *taper,
  const boolean   center,
  const boolean   recenter,
  const boolean   single_sided,
  const sint32    npad,
  void           *intrp_ptr,
  univ_mat       *result )
{
  dcomplex      *pz_idft;
  dcomplex      *pz_dft;
  dcomplex      *pz_result;
  dcomplex_mat   dft;
  dcomplex_mat   idft;
  double        *pd_acvs;
  double        *pd_window;
  double_mat     acvs;
  memlist        list;
  mutil_errcode  err;
  sint32         j;
  sint32         N;
  sint32         nelem;
  sint32         n_series;
  sint32         nsdf;
  sint32         t;
  univ_mat       sdf_direct;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_spectral_density_function_lag_window()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  /* ... lag window */

  LOCALDEF_CHECK_NULL_POINTER_SDF( lag_window, univ_mat, matuniv );

  if ( !MATANY_IS_VEC( &(lag_window->mat.dblmat) ) ){
    MUTIL_ERROR( "Lag window matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( lag_window->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Lag window matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( MATUNIV_NELEM( lag_window ) != MATUNIV_NROW( time_series ) ){
    MUTIL_ERROR( "Time series and lag window matrices must have the same number of elements." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( npad < 2 * MATUNIV_NROW( time_series ) ){
    MUTIL_ERROR( "The input npad must be at least twice the length of "
      "each time series in the input matrix for lag window SDF estimates." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* calculate double-sided direct SDF estimator and register the
     result with the memory manager. in addition,
     use that program for error checking on most inputs
     to avoid redundant checks */

  err = frauniv_spectral_density_function_direct(
    time_series,
    taper,
    center,
	recenter,
    FALSE,
    npad,
    intrp_ptr,
    &sdf_direct );
  if ( err ) return err;

  err = memlist_member_register( &list, &sdf_direct, MEMTYPE_MATUNIV );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* initialize variables */

  N = MATUNIV_NROW( time_series );

  nelem    = single_sided ? npad / 2 + 1 : npad;
  n_series = MATUNIV_NCOL( time_series );
  nsdf     = MATUNIV_NCOL( &sdf_direct );

  /* allocate memory */

  err = matcpx_malloc_register( &idft, npad, nsdf, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &acvs, npad, nsdf, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matcpx_malloc_register( &dft, npad, nsdf, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( result, nelem, nsdf,
		MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate the autocovariance sequence corresponding to the
     direct SDF estimate via an inverse DFT and register the result
     with the memory manager */

  err = sigcpx_transform_discrete_fourier(
    &( sdf_direct.mat.cpxmat ), TRUE, intrp_ptr, &idft );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* form smoothed acvs(s) */
  
  for ( j = 0; j < nsdf; j++ ){
    
    pd_acvs   = acvs.data + j;
    pz_idft   = idft.data + j;
    pd_window = lag_window->mat.dblmat.data;
    
    for ( t = 0; t < npad; t++ ){
      
      if ( t < N ){
        
        *pd_acvs = (*pd_window) * pz_idft->re;
        
        pd_window++;
        pd_acvs += nsdf;
      }
      else if ( t >= N && t <= npad - N ){
        
        *pd_acvs = 0.0;
        
        pd_acvs += nsdf;
        
        /* place lag window pointer at end of lag window */
        
        if ( t == npad - N ){
          pd_window = lag_window->mat.dblmat.data + N - 1;
        }
      }
      else{
        
        *pd_acvs = (*pd_window) * pz_idft->re;
        
        pd_window--;
        pd_acvs += nsdf;
      }

      /* increment IDFT pointer to next coefficient */

      pz_idft += nsdf;

    } /* end loop over time */  
  
  } /* end loop over column */

  /* perform a forward DFT of the smoothed acvs(s) to obtain
     lag window SDF(s) */

  err = sigdbl_transform_discrete_fourier( &acvs, FALSE, intrp_ptr, &dft );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* extract the DFT coefficients (nelem is adjusted for
     single-sided if requested). */

  for ( j = 0; j < nsdf; j++ ){

    pz_dft    = dft.data + j;
    pz_result = result->mat.cpxmat.data + j;

    for ( t = 0; t < nelem; t++ ){

      (*pz_result) = (*pz_dft);

      pz_result += nsdf;
      pz_dft    += nsdf;
    }
  }

  if ( MUTIL_INTERRUPT( 10.0 * (double) N, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }

  /* free nodes corresponding to registered
     memory for the result, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
     corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with frauniv_spectral_density_function_lag_window()" );

  return MUTIL_ERR_OK;
}

/* Welch's Overlapped Segment     */
/* Averaging (WOSA) SDF estimator */
/*                                */
/* Documented in fra_sdf.h        */
/* Written by William Constantine */

mutil_errcode frauniv_spectral_density_function_wosa(
  const univ_mat *time_series,
  const univ_mat *taper,
  const double    overlap,
  const boolean   center,
  const boolean   recenter,
  const boolean   single_sided,
  const sint32    npad,
  void           *intrp_ptr,
  univ_mat       *result )
{
  double         fac;
  double        *pd_segment;
  double        *pd_series;
  double        *pd_taper;
  double        *pd_tappad;
  double_mat     mean;
  memlist        list;
  mutil_errcode  err;
  sint32         n_sample;
  sint32         n_block;
  sint32         block_size;
  sint32         i;
  sint32         j;
  sint32         nelem;
  sint32         n_series;
  sint32         nsdf;
  sint32         n;
  sint32         row_start;
  univ_mat       sdf_direct;
  univ_mat       segment;
  univ_mat       taper_padded;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_spectral_density_function_wosa()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  /* ... time series */

  LOCALDEF_CHECK_NULL_POINTER_SDF( time_series, univ_mat, matuniv );

  if ( time_series->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Time series matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... taper */

  LOCALDEF_CHECK_NULL_POINTER_SDF( taper, univ_mat, matuniv );

  if ( !MATANY_IS_VEC( &(taper->mat.dblmat) ) ){
    MUTIL_ERROR( "Taper matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( taper->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Taper matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( overlap < 0.0 || overlap >= 1.0 ){
    MUTIL_ERROR( "Overlap must be 0 <= overlap < 1." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* initialize variables */

  n_sample   = MATUNIV_NROW( time_series );
  block_size = MATUNIV_NELEM( taper );
  n_block    = 1 + (sint32) floor( (double) ( n_sample - block_size ) /
    ( (double) block_size * ( 1.0 - overlap ) ) );

  nelem    = single_sided ? npad / 2 + 1 : npad;
  n_series = MATUNIV_NCOL( time_series );
  nsdf    = n_series * ( n_series + 1 ) / 2;

  fac = ( n_block == 1 ) ? 0.0 : (double) ( n_sample - block_size ) /
    (double) ( n_block - 1 );

  if ( block_size >=  n_sample ){
    MUTIL_ERROR( "Block width must be less than number of samples in original time series." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* allocate memory */

  err = matuniv_malloc_register( &segment, n_sample, n_series, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &taper_padded, n_sample, 1, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* initialize matrices */

  pd_taper   = taper->mat.dblmat.data;
  pd_tappad  = taper_padded.mat.dblmat.data;

  /* pad block_size point taper out to n_sample points */

  for ( i = 0; i < n_sample; i++ ){

    *pd_tappad  = ( i < block_size ) ? *pd_taper : 0.0;

    pd_taper++;
    pd_tappad++;
  }

  err = matdbl_assign_scalar( 0.0, intrp_ptr, &( segment.mat.dblmat ) );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate mean of time series if center option is true */

  err = matdbl_malloc_register( &mean, 1, n_series, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = localfn_column_means( &( time_series->mat.dblmat ), center,
    intrp_ptr, &mean );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* for each block, calculate the direct SDF using the same taper
     add the result to the 'result' vector */

  for ( n = 0; n < n_block; n++ ){
    
  /* form the segment: the first block_size points are extracted
  from the appropriate block of the original time series.
  the remainder has already been zero padded out to n_sample points.
  subtract mean while we are at it (mean vector may also be all zeros
    so subtracting it in this case won't hurt anything) */
    
    row_start = (sint32) floor( (double) n * fac );
    
    pd_segment = segment.mat.dblmat.data;
    pd_series  = time_series->mat.dblmat.data  + row_start * n_series;
    
    for ( i = 0; i < block_size; i++ ){
      
      for ( j = 0; j < n_series; j++ ){
        
        *pd_segment = *pd_series - mean.data[j];
        
        pd_segment++;
        pd_series++;
      }
    }
    
    /* calculate the direct SDF with a frequency
       resolution of 1/(2*npad). do not center the segment:
       centering is accomplished above by removing the
       mean from the original time series. also, don't
       check for npad above as it will be checked by
       the following function. */

    if ( n == 0 ){

      err = frauniv_spectral_density_function_direct(
        &segment, &taper_padded, FALSE,	recenter, single_sided,
        npad, intrp_ptr, result );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_register( &list, result, MEMTYPE_MATUNIV );
      MEMLIST_FREE_ON_ERROR( err, &list );

    }
    else{

      err = frauniv_spectral_density_function_direct(
        &segment, &taper_padded, FALSE,	recenter, single_sided,
        npad, intrp_ptr, &sdf_direct );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_register( &list, &sdf_direct, MEMTYPE_MATUNIV );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = localfn_add_to_complex_matrix( &sdf_direct, result );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_free( &sdf_direct, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }
  }

  /* average the direct SDF estimates */

  err = localfn_complex_matrix_divide_by_scalar( result, (double) n_block );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* average the direct SDF estimates */

  if ( MUTIL_INTERRUPT( 3.0 * (double) npad, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }

  /* free nodes corresponding to registered
     memory for the result, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
     corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with frauniv_spectral_density_function_wosa()" );

  return MUTIL_ERR_OK;
}


/* Multitaper SDF estimator       */
/*                                */
/* Documented in fra_sdf.h        */
/* Written by William Constantine */

mutil_errcode frauniv_spectral_density_function_multitaper(
  const univ_mat *time_series,
  const univ_mat *taper,
  const boolean   center,
  const boolean   recenter,
  const boolean   single_sided,
  const sint32    npad,
  void           *intrp_ptr,
  univ_mat       *result )
{
  memlist        list;
  mutil_errcode  err;
  sint32         N;
  sint32         i;
  sint32         nelem;
  sint32         ntaper;
  univ_mat       sdf_direct;
  univ_mat       subtaper;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_spectral_density_function_multitaper()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* most inputs are checked via direct
     spectral estimator to avoid redundancy */

  /* initialize variables */

  N     = MATUNIV_NROW( time_series );
  nelem = single_sided ? npad / 2 + 1 : npad;
  ntaper = MATUNIV_NROW( taper );

  if ( MATUNIV_NCOL( taper ) != N ){
    MUTIL_ERROR( "Number of columns in taper matrix must equal "
      "number of elements in time series matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( MATUNIV_NROW( taper ) > N ){
    MUTIL_ERROR( "Number of rows in taper matrix must not exceed "
      "number of elements in time series matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* for each taper, calculate the direct SDF and
     add to the result vector */

  subtaper.type = MUTIL_DOUBLE;
  subtaper.mat.dblmat.nelem = N;
  subtaper.mat.dblmat.ncol  = 1;
  subtaper.mat.dblmat.nrow  = N;
  subtaper.mat.dblmat.data  = taper->mat.dblmat.data;

  /* form the first of the multitaper estimates and register
     the result with the memory manager */

  err = frauniv_spectral_density_function_direct(
    time_series,
    &subtaper,
    center,
	recenter,
    single_sided,
    npad,
    intrp_ptr,
    result );
  if ( err ) return err;

  err = memlist_member_register( &list, result, MEMTYPE_MATUNIV );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* update the window */

  subtaper.mat.dblmat.data += N;

  /* add on the remaining SDF estimates */

  for ( i = 1; i < ntaper; i++ ){

    err = frauniv_spectral_density_function_direct(
      time_series,
      &subtaper,
      center,
	  recenter,
      single_sided,
      npad,
      intrp_ptr,
      &sdf_direct );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_register( &list, &sdf_direct, MEMTYPE_MATUNIV );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = localfn_add_to_complex_matrix( &sdf_direct, result );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_free( &sdf_direct, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* update the window */

    subtaper.mat.dblmat.data += N;
  }

  /* average the direct SDF estimates */

  err = localfn_complex_matrix_divide_by_scalar( result, (double) ntaper );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( MUTIL_INTERRUPT( 3.0 * (double) N, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }

  /* free nodes corresponding to registered
     memory for the result, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
     corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with frauniv_spectral_density_function_multitaper()" );

  return MUTIL_ERR_OK;
}


/** Nonparametric cross-spectral density function estimation.
 * Estimates the spectral density function for a real-valued
 * uniformly sampled time series via a multitaper
 * spectral estimation scheme.
 *
 * References:
 *
 * 1. D.B. Percival and A.T. Walden,
 * ``Spectral Analysis for Physical Applications: Multitaper and
 * Conventional Univariate Techniques'', Cambridge University Press,
 * 1993..
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_sdf.h
 * @source fra\_sdf.c
 * @library fractal
 * @usage #err = localfn_cross_sdf( &mat, n_frequency, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param mat Pointer to a pre-allocated universal matrix of type MUTIL\_DOUBLE.
 *   Each column contains a time series that may have already been centered,
 *   tapered, or other.
 * @param n_frequency The total number of frequencies to extract from the spectra.
 *   Must be less than the number of rows of the input matrix.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param result      Pointer to a universal matrix of type
 *   MUTIL\_DCOMPLEX containing the cross-spectral density function estimates.
 *   For an M-column input matrix, the number columns in this result
 *   matrix will be M*(M+1)/2, which corresponds to all unique cross-spectral
 *   combinations that can be formed. Let Sij = conj(Xi(f)) * Xj(f)
 *   be the cross-spectral density estimate of the ith and jth time series
 *   where Xk(f) is the DFT of the kth time series. If M = 3,
 *   the result will contain the following columns: [S00 S01 S02 S11 S12 S22].
 *   The memory for this matrix is automatically allocated by the function.
 *
 * @see frauniv_spectral_density_function_direct
 * @private
 */
static mutil_errcode localfn_cross_sdf(
  const univ_mat *mat,
  const sint32    n_frequency,
  void           *intrp_ptr,
  univ_mat       *result )
{
  dcomplex      *pz_dft1;
  dcomplex      *pz_dft2;
  dcomplex      *pz_result;
  dcomplex_mat   dft;
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         j;
  sint32         k;
  sint32         t;
  sint32         n_sample;
  sint32         n_series;
  sint32         nsdf;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start localfn_cross_sdf()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* perform basic input checks */

  if ( mat->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Input data matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( n_frequency > MATUNIV_NROW( mat ) ){
    MUTIL_ERROR( "Number of SDF frequencies cannot exceed the number of time series samples" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* initialize variables */

  n_sample    = MATUNIV_NROW( mat );
  n_series    = MATUNIV_NCOL( mat );
  nsdf       = n_series * ( n_series + 1 ) / 2;

  /* allocate memory and register it with the memory manager */

  err = matcpx_malloc_register( &dft, n_sample, n_series, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( result, n_frequency, nsdf, MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate the DFT of the input data matrix */

  err = sigdbl_transform_discrete_fourier( &(mat->mat.dblmat),
    FALSE, intrp_ptr, &dft );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* form Sij = conj(Xi(f)) * Xj(f) for all unique
     combinations of i and j. if the number of columns
     of the original time series matrix is ncol = 3,
     the result will contain the following columns:

     [S00 S01 S02 S11 S12 S22]
  */

  pz_result = result->mat.cpxmat.data;

  for ( k = 0, i = 0; i < n_series; i++ ){
    
    for ( j = i; j < n_series; j++, k++ ){
      
      pz_dft1   = dft.data + i;
      pz_dft2   = dft.data + j;
      pz_result = result->mat.cpxmat.data + k;
      
      for ( t = 0; t < n_frequency; t++ ){
        
        LOCALDEF_COMPLEX_CONJUGATE_PRODUCT( (*pz_dft1), (*pz_dft2), (*pz_result) );
        
        pz_dft1   += n_series;
        pz_dft2   += n_series;
        pz_result += nsdf;
      }
    }
  }

  if ( MUTIL_INTERRUPT( 10.0 * (double) n_sample, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }

  /* free nodes corresponding to registered
     memory for the result, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
     corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_cross_sdf()" );

  return MUTIL_ERR_OK;
}

/** Column means of a matrix.
 * Calculates the column means of a matrix if requested.
 * If not requested, a vector of zeros is returned.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_sdf.h
 * @source fra\_sdf.c
 * @library fractal
 * @usage #err = localfn_column_means( &matrix, center, intrp_ptr, &mean );#
 * @return Standard mutils error/OK code.
 * @param matrix Pointer to a pre-allocated double matrix.
 * @param center A boolean flag. If TRUE, the means of each column are
 *   calculated and returned. Otherwise, a vector of zeros is returned.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param mean Pointer to a pre-allocated single-row double matrix
 *   containing the column means of the input matrix or zeros.
 *   The length of the vector is equal to the number of
 *   columns of the input matrix.
 * @see frauniv_spectral_density_function_direct
 * @see frauniv_spectral_density_function_wosa
 * @private
 */
static mutil_errcode localfn_column_means(
  const double_mat *matrix,
  const boolean     center,
  void             *intrp_ptr,
  double_mat       *mean )
{
  mutil_errcode  err;

  MUTIL_TRACE( "Start localfn_column_means()" );

  if ( center ){

    err = matdbl_sum_cols( matrix, intrp_ptr, mean );
    if ( err ) return err;

    err = matdbl_divide_scalar( mean, (double) matrix->nrow,
      (boolean) TRUE, intrp_ptr, mean );
    if ( err ) return err;
  }
  else{
    err = matdbl_assign_scalar( (double) 0.0, intrp_ptr, mean );
    if ( err ) return err;
  }

  MUTIL_TRACE( "Done with localfn_column_means()" );

  return MUTIL_ERR_OK;
}

/** Columnwise centering of a matrix.
 * Subtracts column means from a matrix. The recentering is
 * done in-place.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_sdf.h
 * @source fra\_sdf.c
 * @library fractal
 * @usage #err = localfn_recenter_columns(&matrix, intrp_ptr);#
 * @return Standard mutils error/OK code.
 * @param matrix Pointer to a pre-allocated double matrix.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @see frauniv_spectral_density_function_direct
 * @see frauniv_spectral_density_function_wosa
 * @private
 */
static mutil_errcode localfn_recenter_columns(
  double_mat *matrix,
  void       *intrp_ptr )
{
  double         *pd_x;
  double_mat     mean;
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         j;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  /* initialize memory list */

  MEMLIST_INIT( list );

  MUTIL_TRACE( "Start localfn_recenter_columns()" );

  /* allocate memory */

  err = matdbl_malloc_register( &mean, 1, matrix->ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );


  /* obtain column means */
  err = matdbl_sum_cols( matrix, intrp_ptr, &mean );
  MEMLIST_FREE_ON_ERROR( err, &list );


  err = matdbl_divide_scalar( &mean, (double) matrix->nrow,
    (boolean) TRUE, intrp_ptr, &mean );
  MEMLIST_FREE_ON_ERROR( err, &list );
  
  /* assign pointers */
  pd_x = matrix->data;

  /* recenter each columns */

  for ( i = 0; i < matrix->nrow; i++ ){

	  for ( j = 0; j < matrix->ncol; j++ ){

		  (*pd_x) -= mean.data[j];
		  pd_x++;
	  }
  }

  /* free local memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_recenter_columns()" );

  return MUTIL_ERR_OK;
}



/** Addition of complex matrices.
 * Adds the contents of one complex universal matrix to
 * the existing complex values of another universal matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_sdf.h
 * @source fra\_sdf.c
 * @library fractal
 * @usage #err = localfn_add_to_complex_matrix( &x, &y );#
 * @return Standard mutils error/OK code.
 * @param x Pointer to a pre-allocated universal matrix of type
 *          MUTIL\_DCOMPLEX.
 * @param y Pointer to a pre-allocated universal matrix of type
 *          MUTIL\_DCOMPLEX. The contents of x are add to this
 *          matrix and returned.
 * @see frauniv_spectral_density_function_wosa
 * @see frauniv_spectral_density_function_multitaper
 * @private
 */
static mutil_errcode localfn_add_to_complex_matrix(
  const univ_mat *x,
  univ_mat       *y )
{
  mutil_errcode   err;
  dcomplex       *pz_x;
  dcomplex       *pz_y;
  sint32          i;

  err = matuniv_verify_aresame( x, y );
  if ( err ) return err;

  if ( x->type != MUTIL_DCOMPLEX ){
    MUTIL_ERROR("Input and output matrices must be complex universal matrices");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* add x to y */

  pz_x = x->mat.cpxmat.data;
  pz_y = y->mat.cpxmat.data;

  for ( i = 0; i < MATUNIV_NELEM( x ); i++ ){

    (*pz_y).re += (*pz_x).re;
    (*pz_y).im += (*pz_x).im;

    pz_y++;
    pz_x++;
  }

  return MUTIL_ERR_OK;
}


/** Division of complex matrix by a scalar.
 * Divides the values of a complex universal matrix
 * by a scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_sdf.h
 * @source fra\_sdf.c
 * @library fractal
 * @usage #err = localfn_complex_matrix_divide_by_scalar( &x, norm );#
 * @return Standard mutils error/OK code.
 * @param x Pointer to a pre-allocated universal matrix of type
 *          MUTIL\_DCOMPLEX.
 * @param norm Double value used to divide the elements of
 *             the input matrix.
 * @see frauniv_spectral_density_function_wosa
 * @see frauniv_spectral_density_function_multitaper
 * @private
 */
static mutil_errcode localfn_complex_matrix_divide_by_scalar(
  const univ_mat *x,
  double          norm )
{

  dcomplex       *pz_x;
  mutil_errcode   err;
  sint32          i;

  LOCALDEF_CHECK_NULL_POINTER_SDF( x, univ_mat, matuniv );

  if ( x->type != MUTIL_DCOMPLEX ){
    MUTIL_ERROR("Input matrix must be complex universal matrix");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* divide each element of x by a double numeric scalar */

  pz_x = x->mat.cpxmat.data;

  for ( i = 0; i < MATUNIV_NELEM( x ); i++ ){

    (*pz_x).re /= norm;
    (*pz_x).im /= norm;

    pz_x++;
  }

  return MUTIL_ERR_OK;
}

