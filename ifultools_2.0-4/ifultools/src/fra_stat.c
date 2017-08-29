
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_stat.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */

#include "fra_stat.h"
#include "fra_sdf.h"

#include "mat_assn.h"
#include "mat_stat.h"
#include "mat_set.h"
#include "mat_summ.h"
#include "mat_type.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "mth_mac.h"
#include "mth_stat.h"
#include "mth_dist.h"

#include "sig_win.h"

#include "str_type.h"

#include "ut_math.h"
#include "ut_mem.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"

#include "ut_mem.h"
#include <math.h>
#include <stdio.h>

/*
This file contains function definitions for
testing for stationarity in a time series.
The functions are declared in fra_stat.h
*/

/* Static macro definitions */

#undef LOCALDEF_CHECK_NULL_POINTER_STAT
#define LOCALDEF_CHECK_NULL_POINTER_STAT( DATA_PTR,         \
  DATA_TYPE, TYPE_PREFIX )                                  \
  err = TYPE_PREFIX ## _validate( DATA_PTR );               \
  if ( err ) return err;                                    \
  if ( DATA_PTR == ( DATA_TYPE * ) NULL ) {                 \
  MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
  return MUTIL_ERR_NULL_POINTER;                            \
  }

/* Static functions declared here and defined at end of file */


/* Stationarity tests for           */
/* univariate time series.          */
/*                                  */
/* Documented in fra_stat.h         */
/* Written by William Constantine   */

mutil_errcode frauniv_stationarity_priestley_subba_rao(
  const univ_mat      *time_series,
  const double         sampling_interval,
  const sint32         n_taper,
  const sint32         n_block,
  const double         significance,
  const boolean        center,
  const boolean        recenter,
  void                *intrp_ptr,
  mat_set             *result )
{
  double         bandwidth;
  double         deltaf;
  double         fac;
  double         mean;
  double         nyquist;
  double        *pd_anova;
  double        *pd_col_mean;
  double        *pd_freq;
  double        *pd_grand_mean;
  double        *pd_row_mean;
  double        *pd_test;
  univ_mat       block;
  univ_mat       sdf_block;
  univ_mat       taper;
  double_mat     time_series_recentered;
  memlist        list;
  mutil_errcode  err;
  sint32         block_size;
  sint32         dims = 3;
  sint32         i;
  sint32         l;
  sint32         ideltaf;
  sint32         ifhigh;
  sint32         iflow;
  sint32         j;
  sint32         n_freq;
  sint32         n_sample;
  sint32         n_block_max;
  sint32_mat     ncol;
  sint32_mat     nrow;
  dcomplex      *pz_sdf;
  
  MUTIL_INTERRUPT_INIT( intrp_ptr );
  
  MUTIL_TRACE( "Start frauniv_stationarity_priestley_subba_rao()" );
  
  /* avoid lint warning */
  
  ( void ) whatssi;
  
  /* initialize memory list */
  
  MEMLIST_INIT( list );
  
  /*** check input data ... ***/
  
  /* ... time series */

  /* ... for valid matrix structure and NULL pointer */
  
  LOCALDEF_CHECK_NULL_POINTER_STAT( time_series, univ_mat, matuniv );
  
  /* ... if number of elements is positive */

  if ( MATUNIV_NELEM( time_series ) < 1 ){
    MUTIL_ERROR( "Number of elements in input time series matrix must be positive." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  Rprintf("1\n");    
  /* ... if it is long enough, based on the blocking
  scheme and the number of tapers to ensure Gaussianity,
  we need at least 2 * (12 * sampling_interval - 1) elements */

  if ( MATUNIV_NELEM( time_series ) < 2 * (sint32) ( 12.0 * sampling_interval - 1.0 ) ){
    MUTIL_ERROR( "Number of elements in input time series matrix must be "
      "at least 2 * ( 12 * sampling_interval - 1 )" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  Rprintf("2\n"); 
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

  /* ... positive sampling interval */

  if ( sampling_interval <= 0.0 ){
    MUTIL_ERROR( "Sampling interval must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  Rprintf("3\n"); 
  /* ... number of tapers */

  if ( n_taper < 5 ){
     MUTIL_ERROR( "Number of tapers must be >= 5" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  Rprintf("4\n");   
  if ( n_taper > MATUNIV_NELEM( time_series ) ){
     MUTIL_ERROR( "Number of tapers cannot exceed number of samples in time series" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  Rprintf("5\n");    
  /* ... number of blocks not too small */

  if ( n_block < 2 ){
    MUTIL_ERROR( "Number of blocks must be at least two" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  Rprintf("6\n");   
  /* check that maximum block size is not exceeded */

  n_block_max = (sint32) ceil( (double) MATUNIV_NELEM( time_series ) / 
     ( 12.0 * sampling_interval - 1.0 ) );

  Rprintf("N = %ld, nblock = %ld, n_block_max = %ld, dt = %10.4f\n", MATUNIV_NELEM( time_series ), n_block, n_block_max, sampling_interval);
  
  if ( n_block > n_block_max ){
     MUTIL_ERROR( "Number of blocks too large" );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }
  

  
  Rprintf("7\n"); 
  /* ... significance */

  if ( significance <= 0.0 || significance >= 1.0 ){	
    MUTIL_ERROR( "Significance must be on the interval (0,1)" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  Rprintf("8\n"); 
  /* initialize variables:
    
      deltaf is the spectral resolution returned by the frauniv_sdf()
      function. for an N-point time series, the frauniv_sdf() function
      returns a SDF estimate with a (normalized) spectral resolution of 1/(2N). 
  */
  
  nyquist    = 1.0 / ( 2.0 * sampling_interval );
  block_size = MATUNIV_NELEM( time_series ) / n_block;
  bandwidth  = (double) ( n_taper + 1 ) / (double) ( block_size + 1 );
  deltaf     = 1.0 / ( (double) ( 2 * block_size ) * sampling_interval );
  ideltaf    = (sint32) ceil( bandwidth / deltaf );
  iflow      = (sint32) ceil( bandwidth / 2.0 / deltaf );
  ifhigh     = (sint32) floor( ( nyquist - bandwidth / 2.0 ) / deltaf );
  
  n_freq     = ( ifhigh - iflow ) / ideltaf + 1;
  n_sample   = n_block * block_size;
  
  fac = log( (double) n_taper ) - MUTIL_DIGAMMA( n_taper );

  /* allocate memory ... */
  
  err = mats32_malloc_register( &nrow, 1, 3, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );
  Rprintf("9\n"); 
  err = mats32_malloc_register( &ncol, 1, 3, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );
  
  /* ... the test vector ( 'between time', 'interaction + residual',
  and 'combination' sum of squares. the next row contains the
  comparative chi-square distribution statistics */
  
  nrow.data[ 0 ] = 2;
  ncol.data[ 0 ] = 3;
  
  /* ... the Priestley-Subba Rao ANOVA table */
  
  nrow.data[ 1 ] = n_block + 1;
  ncol.data[ 1 ] = n_freq + 1;
  
  /* ... the Fourier frequency vector */
  
  nrow.data[ 2 ] = n_freq;
  ncol.data[ 2 ] = 1;
  
  /* ... allocate space for the matrix set and matrices */
  
  err = matset_malloc_register( result, 1, &dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );
  Rprintf("10\n"); 
  err = matset_malloc_matrices_arbitrary_size(
    result, &nrow, &ncol, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );
  Rprintf("11\n"); 
  err = matdbl_malloc_register( &time_series_recentered, n_sample, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );
  Rprintf("12\n"); 
  /* initialize ANOVA matrix */
  
  err = matdbl_assign_scalar( 0.0, intrp_ptr, &( result->mats[ 1 ].mat.dblmat ) );
  MEMLIST_FREE_ON_ERROR( err, &list );
  Rprintf("13\n"); 
  /* set pointers */
  
  pd_test  = result->mats[ 0 ].mat.dblmat.data;
  pd_anova = result->mats[ 1 ].mat.dblmat.data;
  pd_freq  = result->mats[ 2 ].mat.dblmat.data;
  
  /* create the Fourier frequency index vector */
  
  for ( j = 0; j < n_freq; j++ ){
    
    pd_freq[ j ] = ( iflow + ideltaf * j ) * deltaf;
  }
  
  /* recenter the time series */
  
  err = matdbl_sum( &( time_series->mat.dblmat ), intrp_ptr, &mean );
  MEMLIST_FREE_ON_ERROR( err, &list );
  Rprintf("14\n"); 
  mean /= (double) n_sample;
  
  for ( i = 0; i < n_sample; i++ ){
    
    time_series_recentered.data[ i ] = time_series->mat.dblmat.data[ i ] - mean;
  }
  
  /* develop sinusoidal multitapers and register with the memory manager */
  
  err = siguniv_taper( SIG_TAPER_SINUSOIDAL, n_taper, block_size, 
    0.0, TRUE, &intrp_ptr, &taper );
    MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &taper, MEMTYPE_MATUNIV );
  MEMLIST_FREE_ON_ERROR( err, &list );
  Rprintf("15\n"); 
  /* initialize block parameters: impose single-column structure */

  block.type             = MUTIL_DOUBLE;
  block.mat.dblmat.nelem = block_size;
  block.mat.dblmat.nrow  = block_size;
  block.mat.dblmat.ncol  = 1;
  block.mat.dblmat.data  = time_series_recentered.data;

 /* form Priestley-Subba Rao ANOVA matrix */
  
  for ( l = 0; l < n_block; l++ ){
    
    /* form the single-sided multitapered SDF estimate 
    for the current block. since the time series has already
    been centered, do not do so again. remember that these
    estimates are based on a unit sampling interval assumption, 
    so adjust accordingly */

    err = frauniv_spectral_density_function_multitaper(
      &block,
      &taper,
      center,
	  recenter,
      TRUE,
      2 * block_size,
      intrp_ptr,
      &sdf_block );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_register( &list, &sdf_block, MEMTYPE_MATUNIV );
    MEMLIST_FREE_ON_ERROR( err, &list );


    /* subsample the SDF estimator at the appropriate frequencies and
       adjust for sampling interval. form corresponding ANOVA statistic */

     pz_sdf = sdf_block.mat.cpxmat.data + iflow;

     for ( j = 0; j < n_freq; j++ ){
        
        *pd_anova = log( (*pz_sdf).re * sampling_interval ) + fac;

        pz_sdf += ideltaf;
        pd_anova++;
     }
      
     /* increment ANOVA pointer to skip row-mean column */
     
     pd_anova++;
    
     /* update block pointer */

     block.mat.dblmat.data += block_size;

     /* free the SDF matrix */

     err = memlist_member_free( &sdf_block, &list );
     MEMLIST_FREE_ON_ERROR( err, &list );
  }
  Rprintf("16\n"); 
  /* calculate ANOVA row, column and grand means */
  
  pd_anova      = result->mats[ 1 ].mat.dblmat.data;
  pd_row_mean   = pd_anova + n_freq;
  pd_col_mean   = pd_anova + n_block * ( n_freq + 1 );
  pd_grand_mean = pd_anova + ( n_block + 1 ) * ( n_freq + 1 ) - 1;
  
  for ( l = 0; l < n_block; l++ ){
    
   for ( j = 0; j < n_freq; j++ ){
      
      *pd_row_mean   += *pd_anova;
      *pd_col_mean   += *pd_anova / (double) n_block;
      *pd_grand_mean += *pd_anova;

       pd_anova++;
       pd_col_mean++; 
   }

   *pd_row_mean /= (double) n_freq;

   /* reset pointers */

   pd_anova++; 
   pd_row_mean += n_freq + 1;
   pd_col_mean -= n_freq;
  }

  *pd_grand_mean /= (double) ( n_block * n_freq );
    
  if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }
  
  /* form test statistics */
  
  pd_test[ 0 ] = 0.0;
  pd_test[ 1 ] = 0.0;
  
  pd_row_mean = &( result->mats[ 1 ].mat.dblmat.data[ n_freq ] );
  pd_anova    = result->mats[ 1 ].mat.dblmat.data;
  
  for ( l = 0; l < n_block; l++ ){
    
    pd_test[ 0 ] += MUTIL_SQR( *pd_row_mean - *pd_grand_mean );
    
    for ( j = 0; j < n_freq; j++ ){
      
      pd_test[ 1 ] += MUTIL_SQR( *pd_anova - *pd_row_mean - pd_col_mean[ j ] + *pd_grand_mean );
      
      pd_anova++;
    }
    
    pd_anova++;  /* this is needed to skip over row mean entry
                 (in the last column of ANOVA matrix and onto
    the first entry of the next row */
    
    pd_row_mean += n_freq + 1;
  }
  
  pd_test[ 0 ] *= (double) n_freq;
  
  /* normalize tests with trigamma( n_taper ) */
  
  fac = MUTIL_TRIGAMMA( n_taper );
  
  pd_test[ 0 ] /= fac;
  pd_test[ 1 ] /= fac;
  pd_test[ 2 ] = pd_test[ 0 ] + pd_test[ 1 ];
  
  /* form comparative chi-square distribution statistics */
  
  fac = 1.0 - significance;
  
  pd_test[ 3 ] = mth_qchisq( fac, (double) ( n_block - 1 ) );
  pd_test[ 4 ] = mth_qchisq( fac, (double) ( ( n_block - 1 ) * ( n_freq - 1 ) ) );
  pd_test[ 5 ] = mth_qchisq( fac, (double) ( ( n_block - 1 ) * n_freq ) );
  
  /* free nodes corresponding to registered
  memory for the result, but do not free
  the memory itself */
  
  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );
  Rprintf("17\n"); 
  /* free other malloced space and
  corresponding nodes in memory list */
  
  MUTIL_FREE_WARN( memlist, &list );
  
  MUTIL_TRACE( "Done with frauniv_embed()" );
  Rprintf("18\n"); 
  return MUTIL_ERR_OK;
}
