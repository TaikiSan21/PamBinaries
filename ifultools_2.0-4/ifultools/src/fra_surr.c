
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_surr.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */


#include "fra_surr.h"
#include "fra_type.h"

#include "mat_arit.h"
#include "mat_sort.h"
#include "mat_stat.h"
#include "mat_type.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "sig_tran.h"

#include "str_type.h"

#include "ut_math.h"
#include "ut_mem.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"

#include "ut_mem.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

/*
  This file contains function definitions for
  creating surrogate data for time series.
  The functions are declared in fra_surr.h
*/

/* Static macro definitions */

#undef LOCALDEF_IS_EVEN
#define LOCALDEF_IS_EVEN( n ) ( (boolean) (((n) % 2) == 0) )

#undef LOCALDEF_CHECK_NULL_POINTER_SURR
#define LOCALDEF_CHECK_NULL_POINTER_SURR( DATA_PTR,         \
  DATA_TYPE, TYPE_PREFIX )                                  \
  err = TYPE_PREFIX ## _validate( DATA_PTR );               \
  if ( err ) return err;                                    \
  if ( DATA_PTR == ( DATA_TYPE * ) NULL ) {                 \
  MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
  return MUTIL_ERR_NULL_POINTER;                            \
  }

/* Static functions declared here and defined at end of file */

static mutil_errcode localfn_recenter_data(
  const univ_mat *x, void *intrp_ptr, double *mean,
  double *variance, univ_mat *result );

static mutil_errcode localfn_random_normal( double_mat *x, uint32_t seed);
static mutil_errcode localfn_random_uniform_phase( dcomplex_mat *x, uint32_t seed );
static mutil_errcode localfn_circulant_embedding_weights( dcomplex_mat *x, uint32_t seed );
static mutil_errcode localfn_random_phase_weights( dcomplex_mat *x, uint32_t seed );
static mutil_errcode localfn_davison_hinkley_weights( dcomplex_mat *x, uint32_t seed );

static mutil_errcode localfn_rank_and_sort( double_mat *x,
  void *intrp_ptr, double_mat *sorted, sint32_mat *rank );

static int32_t localfn_time_seed(void);

static double localfn_random_uniform_deviate_marsaglia(
  uint32_t *pSeed );

static double localfn_random_normal_deviate_marsaglia(
  uint32_t *pSeed );

/* static variables for the random number generator */

/* static short mother1[10]; */
/* static short mother2[10]; */
static int16_t mother1[10];
static int16_t mother2[10];
static boolean INITIALIZE_RANDOM_NUMBER_GENERATOR = FALSE;

#define m16Long    65536L        /* 2^16                   */
#define m16Mask    0xFFFF        /* mask for lower 16 bits */
#define m15Mask    0x7FFF        /* mask for lower 15 bits */
#define m31Mask    0x7FFFFFFF    /* mask for 31 bits       */
#define m32Double  4294967295.0  /* 2^32-1                 */
#define E22        1.7155277699214135929603792825575449562416 /*2*sqrt(2/e)*/

/* Platform dependent code to get somewhat random numbers that
   could be used for a random seed. */

#ifdef WIN32 /* ( */
#include<windows.h>
static int32_t localfn_time_seed(void)
{
  int32_t val;
  int32_t oldval;

  val = oldval = GetTickCount();

  /* guarantee at least one millisecond has elapsed
  to generate unique seeds between successive
  calls which occur within one millisecond */

  while( val == oldval ){
    val =  GetTickCount();
  }
  return( (int32_t) MUTIL_ABS( val ) );
}
#else /* )( */
#include <sys/time.h>
static int32_t localfn_time_seed(void)
{
  struct timeval tv ;
  (void)gettimeofday(&tv, (void*)NULL);
  return( (int32_t) MUTIL_ABS( tv.tv_usec ) );

}
#endif


/* Surrogate data generation for a  */
/* univariate time series (Theiler) */
/*                                  */
/* Documented in fra_surr.h         */
/* Written by William Constantine   */

mutil_errcode frauniv_bootstrap_theiler(
  const univ_mat      *time_series,
  const fra_surrogate  method,
  const uint32         seed,
  void                *intrp_ptr,
  univ_mat            *result )
{
  dcomplex_mat   dft_inverse;
  dcomplex_mat   dft_random;
  dcomplex_mat   random_weight;
  double         mean;
  double         variance;
  double        *pd_result;
  double_mat     random;
  double_mat     sorted;
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         N;
  sint32         ncol;
  sint32         nrow;
  sint32_mat     order;
  sint32_mat     rank;
  univ_mat       random_phase_univ;
  univ_mat       random_univ;
  univ_mat       time_series_recentered;
  double_mat     time_series_sorted;
  double_mat     temp;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_bootstrap_theiler()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /*** check input data ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_SURR( time_series, univ_mat, matuniv );

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

  /* ... method */

  switch( method ){

    case FRA_SURROGATE_AAFT:
    case FRA_SURROGATE_RANDOM_PHASE:
      break;

    default:

      MUTIL_ERROR( "Surrogate data method is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* initialize variables */

  N    = MATUNIV_NELEM( time_series );
  nrow = MATUNIV_NROW( time_series );
  ncol = MATUNIV_NCOL( time_series );

  /* allocate common memory */

  err = matuniv_malloc_register( result, N, 1, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  switch( method ){

    case FRA_SURROGATE_AAFT:

      /* allocate memory */

      err = matdbl_malloc_register( &random, nrow, ncol, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = matdbl_malloc_register( &sorted, nrow, ncol, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = mats32_malloc_register( &order, nrow, ncol, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* recenter the data */

      err = localfn_recenter_data( time_series, intrp_ptr, &mean, &variance,
				    &time_series_recentered );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_register( &list, &time_series_recentered, MEMTYPE_MATUNIV );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /*
	      create a zero mean unit variance Gaussian white
	      noise sequence and sort it based upon the rank of the
	      indices of the original time series
      */

      err = localfn_random_normal( &random, (uint32_t) seed );
      MEMLIST_FREE_ON_ERROR( err, &list );
      
      err = localfn_rank_and_sort( &( time_series_recentered.mat.dblmat ), intrp_ptr,
        &time_series_sorted, &rank );
      MEMLIST_FREE_ON_ERROR( err, &list );
      
      err = memlist_member_register( &list, &time_series_sorted, MEMTYPE_MATDBL );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_register( &list, &rank, MEMTYPE_MATS32 );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* rank sort the random Gaussian realization */

      err = matdbl_sort_index_partial( &random, (sint32_mat *) NULL, intrp_ptr, &order );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = matdbl_permute( &random, &order, intrp_ptr, &sorted );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = matdbl_permute( &sorted, &rank, intrp_ptr, &random );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_free( &sorted, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_free( &rank, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_free( &order, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* phase-randomize the rank-sorted Gaussian sequence */

      err = matuniv_wrap_matrix_register( &random_univ, &random, MUTIL_DOUBLE, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = frauniv_bootstrap_theiler(
        &random_univ,
        (fra_surrogate) FRA_SURROGATE_RANDOM_PHASE,
        seed,
        intrp_ptr,
        &random_phase_univ );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_register( &list, &random_phase_univ, MEMTYPE_MATUNIV );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_free( &random_univ, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* return a surrogate series by sorting the original
      time series based upon the rank indexing of the
      randomized phase result */
      
      err = localfn_rank_and_sort( &( random_phase_univ.mat.dblmat ),
        intrp_ptr, &sorted, &rank );
      MEMLIST_FREE_ON_ERROR( err, &list );
      
      MUTIL_FREE_WARN( matdbl, &sorted );
      
      err = memlist_member_register( &list, &rank, MEMTYPE_MATDBL );
      MEMLIST_FREE_ON_ERROR( err, &list );
      
      err = matdbl_permute( &time_series_sorted, &rank,
        intrp_ptr, &( result->mat.dblmat ) );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* restore the mean */

      err = matdbl_add_scalar( &( result->mat.dblmat ), mean, intrp_ptr, &( result->mat.dblmat ) );
      MEMLIST_FREE_ON_ERROR( err, &list );

      break;

    case FRA_SURROGATE_RANDOM_PHASE:

      /* allocate memory */

      err = matcpx_malloc_register( &dft_random, N, 1, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = matcpx_malloc_register( &dft_inverse, N, 1, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = matcpx_malloc_register( &random_weight, N, 1, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* calculate the forward DFT of the series: MUTILS
       DFT function only takes dat stored in columns. so trick
       it since we knwo we have a vector . */

      temp.nelem = N;
      temp.ncol  = 1;
      temp.nrow  = N;
      temp.data  = time_series->mat.dblmat.data;

      err = sigdbl_transform_discrete_fourier( &temp,
        (boolean) FALSE, intrp_ptr, &dft_random );
      MEMLIST_FREE_ON_ERROR( err, &list );
      
      /* create the random phase weights */
      
      err = localfn_random_phase_weights( &dft_random, seed );
      MEMLIST_FREE_ON_ERROR( err, &list );
      
      /* invert with the inverse DFT to create surrogate */
      /* calculate the inverse transform */
      
      err = sigcpx_transform_discrete_fourier( &dft_random,
        (boolean) TRUE, intrp_ptr, &dft_inverse );
      MEMLIST_FREE_ON_ERROR( err, &list );
      
      pd_result = result->mat.dblmat.data;

      for ( i = 0; i < N; i++ ){
        
        *pd_result = dft_inverse.data[ i ].re;
        pd_result++;
      }

      break;

  default:

      MUTIL_ERROR( "Surrogate data method is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;

  } /* end switch */

  if ( MUTIL_INTERRUPT( 3.0 * N, intrp_ptr ) ) {
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

  MUTIL_TRACE( "Done with frauniv_bootstrap_theiler()" );

  return MUTIL_ERR_OK;
}

/* Surrogate data generation for a  */
/* univariate time series.          */
/* (Davison-Hinkley)                */
/*                                  */
/* Documented in fra_surr.h         */
/* Written by William Constantine   */

mutil_errcode frauniv_bootstrap_davison_hinkley(
  const univ_mat      *time_series,
  const uint32         seed,
  void                *intrp_ptr,
  univ_mat            *result )
{
  dcomplex_mat   dft_inverse;
  dcomplex_mat   dft_random;
  dcomplex_mat   random_weight;
  double        *pd_result;
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         N;
  double_mat     temp;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_bootstrap_davison_hinkley()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /*** check input data ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_SURR( time_series, univ_mat, matuniv );

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

  /* initialize variables */

  N = MATUNIV_NELEM( time_series );

  /* allocate common memory */

  err = matuniv_malloc_register( result, N, 1, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* allocate memory */

  err = matcpx_malloc_register( &dft_random, N, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matcpx_malloc_register( &dft_inverse, N, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matcpx_malloc_register( &random_weight, N, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate the forward DFT of the series: MUTILS
     DFT function only takes dat stored in columns. so trick
     it since we knwo we have a vector . */

  temp.nelem = N;
  temp.ncol  = 1;
  temp.nrow  = N;
  temp.data  = time_series->mat.dblmat.data;

  err = sigdbl_transform_discrete_fourier( &temp,
    (boolean) FALSE, intrp_ptr, &dft_random );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create the random phase and amplitude weights */

  err = localfn_davison_hinkley_weights( &dft_random, seed );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* invert with the inverse DFT to create surrogate */
  /* calculate the inverse transform */

  err = sigcpx_transform_discrete_fourier( &dft_random,
    (boolean) TRUE, intrp_ptr, &dft_inverse );
  MEMLIST_FREE_ON_ERROR( err, &list );

  pd_result = result->mat.dblmat.data;

  for ( i = 0; i < N; i++ ){

    *pd_result = dft_inverse.data[ i ].re;
    pd_result++;
  }

  if ( MUTIL_INTERRUPT( 3.0 * N, intrp_ptr ) ) {
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

  MUTIL_TRACE( "Done with frauniv_bootstrap_davison_hinkley()" );

  return MUTIL_ERR_OK;
}


/* Surrogate data generation for a   */
/* univariate time series (circulant */
/* embedding)                        */
/* Documented in fra_surr.h          */
/* Written by William Constantine    */

mutil_errcode frauniv_bootstrap_circulant_embedding(
  const univ_mat      *sdf,
  const uint32         seed,
  void                *intrp_ptr,
  univ_mat            *result )
{
  dcomplex_mat   dft_inverse;
  dcomplex_mat   dft_random;
  dcomplex_mat   random_weight;
  double         weight;
  double        *pd_result;
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         N;
  sint32         NN;
  dcomplex      *pz_random;
  dcomplex      *pz_random2;
  dcomplex      *pz_weight;
  double        *pd_sdf;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_bootstrap_circulant_embedding()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /*** check input data ... ***/

  /* ... SDF */

  LOCALDEF_CHECK_NULL_POINTER_SURR( sdf, univ_mat, matuniv );

  /* ... if the type is double */

  if ( sdf->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Input SDF matrix must be of type double." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( &( sdf->mat.dblmat ) ) ){
    MUTIL_ERROR( "Input SDF matrix must be a single-row or single-column matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* initialize variables */

  N    = MATUNIV_NELEM( sdf ) - 1;
  NN   = N * 2;

  /* allocate common memory */

  err = matuniv_malloc_register( result, N, 1, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* allocate memory */

  err = matcpx_malloc_register( &dft_random, NN, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matcpx_malloc_register( &dft_inverse, NN, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matcpx_malloc_register( &random_weight, N + 1, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate the circulant embedding weights */

  err = localfn_circulant_embedding_weights( &random_weight, (uint32_t) seed );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* set pointers */

  pd_result = result->mat.dblmat.data;
  pd_sdf    = sdf->mat.dblmat.data;
  pz_random = dft_random.data;
  pz_weight = random_weight.data;

  /* create complex-valued symmetric and randomized sequence
     based on the SDF and weights ... */

  /* ... upper half coefficients (including purely reals) */

  for ( i = 0; i <= N; i++ ){

    weight = sqrt( *pd_sdf );

    pz_random->re = pz_weight->re * weight;
    pz_random->im = pz_weight->im * weight;

    pz_weight++;
    pz_random++;
    pd_sdf++;
  }

  /* ... form complex conjugate symmetry in lower half */

  pz_random2 = pz_random - 2;

  for ( i = N + 1; i < NN; i++ ){

    pz_random->re = pz_random2->re;
    pz_random->im = - pz_random2->im;

    pz_random++;
    pz_random2--;
  }

  /* calculate the FORWADRD DFT */

  err = sigcpx_transform_discrete_fourier( &dft_random,
	(boolean) FALSE, intrp_ptr, &dft_inverse );
  MEMLIST_FREE_ON_ERROR( err, &list );

  for ( i = 0; i < N; i++ ){

    *pd_result = dft_inverse.data[ i ].re;
    pd_result++;
  }

  if ( MUTIL_INTERRUPT( 3.0 * N, intrp_ptr ) ) {
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

  MUTIL_TRACE( "Done with frauniv_bootstrap_circulant_embedding()" );

  return MUTIL_ERR_OK;
}


/*******************************/
/* STATIC FUNCTION DEFINITIONS */
/*******************************/

/** Realization of a Gaussian distributed random white noise process.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_surr.c
 * @library fractal
 * @usage #err = localfn_random_normal( &x, seed );#
 * @return Standard mutils error/OK code.
 * @param  x Pointer to a pre-allocated double matrix
 *           which upon return will contain the realization.
 * @param seed uint32_t integer representing the initial
 *             random seed value. If zero, the seed is reset
 *             based on the current time.
 * @see localfn_random_uniform_phase
 * @private
 */
static mutil_errcode localfn_random_normal( double_mat *x, uint32_t seed)
{
  sint32         i;
  double        *pd_data;

  MUTIL_TRACE( "Start localfn_random_normal()" );

  /* initiate random number generator and set the seed */

  INITIALIZE_RANDOM_NUMBER_GENERATOR = TRUE;
  if ( seed == 0 ) seed = localfn_time_seed();

  pd_data = x->data;

  for ( i = 0; i < x->nelem; i++ ){

    /* generate a normally distributed random deviate */

    *pd_data = localfn_random_normal_deviate_marsaglia( &seed );

    pd_data++;
  }

  MUTIL_TRACE( "Done with localfn_random_normal()" );

  return MUTIL_ERR_OK;
}

/** Complex circulant embedding weights, N(0,1) Gaussian.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_surr.c
 * @library fractal
 * @usage #err = localfn_circulant_embedding_weights( &x, seed );#
 * @return Standard mutils error/OK code.
 * @param  x Pointer to a pre-allocated complex matrix
 *           which upon return will contain the complex weights
 *           needed for the circulant embedding technique.
 * @param seed uint32_t integer representing the initial
 *             random seed value. If zero, the seed is reset
 *             based on the current time.
 * @see localfn_random_normal
 * @private
 */
static mutil_errcode localfn_circulant_embedding_weights( dcomplex_mat *x, uint32_t seed )
{
  dcomplex      *pz_weight;
  double         norm_real;
  double         norm_upper_half;
  double        *pd_deviate;
  double_mat     deviates;
  memlist        list;
  mutil_errcode  err;
  sint32         N;
  sint32         Nw;
  sint32         NN;
  sint32         k;

  MUTIL_TRACE( "Start localfn_circulant_embedding_weights()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* initialize variables */

  Nw = x->nelem;
  N  = Nw - 1;
  NN = N * 2;

  norm_real       = 1.0 / sqrt( 2.0 * (double) N );
  norm_upper_half = 1.0 / sqrt( 4.0 * (double) N );

  /* allocate memory for the Gaussian deviates */

  err = matdbl_malloc_register( &deviates, NN, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = localfn_random_normal( &deviates, seed );
  MEMLIST_FREE_ON_ERROR( err, &list );

  pd_deviate = deviates.data;
  pz_weight  = x->data;

  for ( k = 0; k <= N; k++ ){

    if ( k == 0 || k == N ){

      pz_weight->re = *pd_deviate * norm_real;
      pz_weight->im = 0.0;

      pd_deviate++;
    }
    else{

      pz_weight->re = *pd_deviate * norm_upper_half;
      pd_deviate++;

      pz_weight->im = *pd_deviate * norm_upper_half;
      pd_deviate++;
    }

    pz_weight++;
  }

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_circulant_embedding_weights()" );

  return MUTIL_ERR_OK;
}

/** Complex random phase weights for Theiler's routine.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_surr.c
 * @library fractal
 * @usage #err = localfn_random_phase_weights( &x, seed );#
 * @return Standard mutils error/OK code.
 * @param  x Pointer to a pre-allocated complex matrix
 *           containing the DFT of a time series.
 *           Upon return, the DFT series is replaced by one modulated
 *           with random phase.
 * @param seed uint32_t integer representing the initial
 *             random seed value. If zero, the seed is reset
 *             based on the current time.
 * @see localfn_random_normal
 * @private
 */
static mutil_errcode localfn_random_phase_weights( dcomplex_mat *x, uint32_t seed )
{
  dcomplex       product;
  dcomplex      *pz_phase;
  dcomplex      *pz_weight2;
  dcomplex      *pz_weight;
  dcomplex_mat   phase;
  memlist        list;
  mutil_errcode  err;
  sint32         N;
  sint32         Ntop;
  sint32         i;

  MUTIL_TRACE( "Start localfn_random_phase_weights()" );

  MEMLIST_INIT( list );

  /* initialize variables */

  N    = x->nelem;
  Ntop = ( N - 1 ) / 2;

  /* create complex-valued random phase vector
   uniformly distributed on [0, 2*PI] */

  err = matcpx_malloc_register( &phase, Ntop, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = localfn_random_uniform_phase( &phase, seed );
  if ( err ) return err;

  /* set pointer to second coefficient
   in the upper z-plane */

  pz_weight = x->data + 1;
  pz_phase  = phase.data;

  for ( i = 1; i <= Ntop; i++ ){

    MUTIL_CPX_MULT( *pz_phase, *pz_weight, product );

    pz_weight->re = product.re;
    pz_weight->im = product.im;

    pz_weight++;
    pz_phase++;
  }

  /* set second pointer to last complex
     DFT coefficient in upper z-plane */

  pz_weight2 = pz_weight - 1;

  /* if the number of coefficients is even, leave
     the DFT coefficient at normalized frequency
     1/2 alone (skip over it). */

  if ( LOCALDEF_IS_EVEN( N ) ){

    pz_weight++;
  }

  /* form complex conjugate symmetry */

  for ( i = N / 2 + 1; i < N; i++ ){

    pz_weight->re = pz_weight2->re;
    pz_weight->im = -pz_weight2->im;

    pz_weight++;
    pz_weight2--;
  }

  /* free the memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_random_phase_weights()" );

  return MUTIL_ERR_OK;
}

/** Weights for Davison-Hinkley.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_surr.c
 * @library fractal
 * @usage #err = localfn_davison_hinkley_weights( &x, seed );#
 * @return Standard mutils error/OK code.
 * @param  x Pointer to a pre-allocated complex matrix
 *           containing the DFT of a time series.
 *           Upon return, the DFT series is replaced by one modulated
 *           with random phase and amplitude.
 * @param seed uint32_t integer representing the initial
 *             random seed value. If zero, the seed is reset
 *             based on the current time.
 * @see localfn_random_normal
 * @private
 */
static mutil_errcode localfn_davison_hinkley_weights( dcomplex_mat *x, uint32_t seed )
{
  dcomplex       z;
  dcomplex      *pz_phase;
  dcomplex      *pz_weight2;
  dcomplex      *pz_weight;
  dcomplex_mat   phase;
  double         norm = 0.70710678118654746; /* 2^{-0.5} */
  memlist        list;
  mutil_errcode  err;
  sint32         N;
  sint32         Ntop;
  sint32         i;

  MUTIL_TRACE( "Start localfn_davison_hinkley_weights()" );

  MEMLIST_INIT( list );

  /* initialize variables */

  N    = x->nelem;
  Ntop = ( N - 1 ) / 2;

  /* create complex-valued random phase vector
   uniformly distributed on [0, 2*PI] */

  err = matcpx_malloc_register( &phase, N - 1, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = localfn_random_uniform_phase( &phase, seed );
  if ( err ) return err;

  /* set pointer to second coefficient
     in the upper z-plane. randomize phase
     for all coefficients except the first */

  pz_weight = x->data + 1;
  pz_phase  = phase.data;

  for ( i = 1; i <= N - 1; i++ ){

    MUTIL_CPX_MULT( *pz_phase, *pz_weight, z );

    pz_weight->re = z.re;
    pz_weight->im = z.im;

    pz_weight++;
    pz_phase++;
  }

  /* average DFT coefficients that are opposite
     each other in the unit circle of the z-plane.
     store the result in the upper half and the conjugate
     of the result in the lower half coefficient */

  pz_weight  = x->data + 1;
  pz_weight2 = x->data + N - 1;

  for ( i = 1; i < Ntop; i++ ){

    MUTIL_CPX_ADD( *pz_weight, *pz_weight2, z );

    z.re *= norm;
    z.im *= norm;

    pz_weight->re = z.re;
    pz_weight->im = z.im;

    pz_weight2->re = z.re;
    pz_weight2->im = -z.im;

    pz_weight++;
    pz_weight2--;
  }

  /* address the purely real DFT coefficient at
     at normalized frequency 1/2 if N is even */

  if ( LOCALDEF_IS_EVEN( N ) ){

    pz_weight++;

    MUTIL_CPX_ADD( *pz_weight, *pz_weight, z );

    z.re *= norm;
    z.im *= norm;

    pz_weight->re = z.re;
    pz_weight->im = z.im;
  }

  /* free the memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_davison_hinkley_weights()" );

  return MUTIL_ERR_OK;
}


/** Creates a uniformly distributed random phase vector on [0, 2*PI].
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_surr.c
 * @library fractal
 * @usage #err = localfn_random_uniform_phase( &x, seed );#
 * @return Standard mutils error/OK code.
 * @param  x Pointer to a pre-allocated complex matrix
 *           which upon return will contain the realization.
 * @param seed uint32_t integer representing the initial
 *             random seed value. If zero, the seed is reset
 *             based on the current time.
 * @see localfn_random_normal
 * @private
 */
static mutil_errcode localfn_random_uniform_phase( dcomplex_mat *x, uint32_t seed )
{
  sint32         i;
  double         deviate;

  MUTIL_TRACE( "Start localfn_random_uniform_phase()" );

  /* initialize the random number generator and set seed */

  INITIALIZE_RANDOM_NUMBER_GENERATOR = TRUE;
  if ( seed == 0 ) seed = localfn_time_seed();

  for ( i = 0; i < x->nelem; i++ ){

    /* generate a uniformly distributed random deviate */

    deviate = localfn_random_uniform_deviate_marsaglia( &seed );

    deviate *= 2.0 * MUTIL_PI;

    x->data[ i ].re = cos( deviate );
    x->data[ i ].im = sin( deviate );
  }

  MUTIL_TRACE( "Done with localfn_random_uniform_phase()" );

  return MUTIL_ERR_OK;
}


/** Index and rank tables.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_surr.c
 * @library fractal
 * @usage #err = localfn_rank_and_sort( &x, intrp_ptr, &index, &rank );#
 * @return Standard mutils error/OK code.
 * @param  x Pointer to a pre-allocated double matrix
 *         to form the index and rank tables.
 * @param  intrp_ptr Pointer for implementation of interrupt checking.
 * @param  index Pointer to a pre-allocated sint32 matrix
 *         of the same size as the input. Upon return,
 *         this matrix will contain the flattened index table
 *         corresponding to x such that x.data[index.data[i]]
 *         is the ith smallest value of x.
 * @param  rank Pointer to a pre-allocated sint32 matrix
 *         of the same size as the input. Upon return,
 *         this matrix will contain the flattened rank table
 *         corresponding to x such that rank.data[i]
 *         is the rank of the ith value in x (in increasing
 *         algebraic order). e.g. a rank of 5 means that
 *         the current element ranks 5th out of all the data
 *         in x.
 * @private
 */
static mutil_errcode localfn_rank_and_sort(
  double_mat *x, void *intrp_ptr, double_mat *sorted, sint32_mat *rank )
{
  mutil_errcode  err;
  sint32         i;
  sint32        *ps_index;
  double_mat     temp;
  sint32_mat     index;
  sint32         nrow = x->nrow;
  sint32         ncol = x->ncol;
  memlist        list;

  /* initialize memory list */

  MEMLIST_INIT( list );

  MUTIL_TRACE( "Start localfn_rank_and_sort()" );

  err = matdbl_malloc_register( &temp, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( sorted, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( rank, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &index, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create index table */

  err = matdbl_sort_index_partial( x, (sint32_mat *) NULL, intrp_ptr, &index );
  if ( err ) {
    MUTIL_ERROR( "Problem encountered in creating index table." );
    return err;
  }

  /* create rank table */

  ps_index = index.data;

  for ( i = 0; i < x->nelem; i++ ){

    rank->data[ *ps_index ] = i;

    ps_index++;
  }

  /* sort the input */

  err = matdbl_permute( x, &index, intrp_ptr, sorted );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free nodes corresponding to registered
  memory for the result, but do not free
  the memory itself */

  err = memlist_member_unregister( sorted, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_unregister( rank, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
  corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_rank_and_sort()" );

  return MUTIL_ERR_OK;
}

static mutil_errcode localfn_recenter_data(
  const univ_mat *x, void *intrp_ptr, double *mean, double *variance, univ_mat *result )
{
  mutil_errcode  err;
  memlist        list;

  /* initialize memory list */

  MEMLIST_INIT( list );

  MUTIL_TRACE( "Start localfn_recenter_data()" );

  /* allocate space for the recentered data */

  err = matuniv_malloc_register( result, MATUNIV_NROW( x ), MATUNIV_NCOL( x ), x->type, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate the mean and sample variance */

  err = matuniv_mean_variance( x, FALSE, intrp_ptr, mean, variance );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* form the result */

  err = matdbl_add_scalar( &( x->mat.dblmat ), - *mean, intrp_ptr, &( result->mat.dblmat ) );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free nodes corresponding to registered
     memory for the result, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
  corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_recenter_data()" );

  return MUTIL_ERR_OK;
}



static double localfn_random_uniform_deviate_marsaglia(
  uint32_t *pSeed )
{
  //uint32_t number, number1, number2;
  uint32_t number, number1, number2;
    
  int16_t n, *p;

  uint16_t sNumber;

  /* Initialize motheri with 9 random values the first time */

  if ( INITIALIZE_RANDOM_NUMBER_GENERATOR ){

    sNumber = (uint16_t) (*pSeed & m16Mask);   /* The low 16 bits */
    number  = *pSeed & m31Mask;   /* Only want 31 bits */

    p = mother1;

    for ( n = 18; n--; ){

      number = 30903 * sNumber + ( number >> 16 );   /* One line multiply-with-carry */

      sNumber = (uint16_t) ( number & m16Mask );

      *p++ = (int16_t) ( number & m16Mask );

      if ( n == 9 ){
	      p = mother2;
      }
    }

    /* make cary 15 bits */

    mother1[0] &= m15Mask;
    mother2[0] &= m15Mask;
    INITIALIZE_RANDOM_NUMBER_GENERATOR = FALSE;
  }

  /* Move elements 1 to 8 to 2 to 9 */

  /* (void *) memmove( mother1 + 2, mother1 + 1, 8 * sizeof(short) ); */
  /* (void *) memmove( mother2 + 2, mother2 + 1, 8 * sizeof(short) ); */
  memmove( mother1 + 2, mother1 + 1, 8 * sizeof(int16_t) );
  memmove( mother2 + 2, mother2 + 1, 8 * sizeof(int16_t) );

  /* Put the carry values in numberi */

  number1 = mother1[ 0 ];
  number2 = mother2[ 0 ];
  
  /* printf("\nnumber1: %" PRIu32 ", mother1[0]: %" PRIi16 "\n", number1, mother1[0]); */

  /* Form the linear combinations */

  number1 += 1941 * mother1[2] + 1860 * mother1[3] + 1812 * mother1[4] + 1776*mother1[5]+

    1492*mother1[6]+1215*mother1[7]+1066*mother1[8]+12013*mother1[9];

  number2+=1111*mother2[2]+2222*mother2[3]+3333*mother2[4]+4444*mother2[5]+

    5555*mother2[6]+6666*mother2[7]+7777*mother2[8]+9272*mother2[9];

  /* Save the high bits of numberi as the new carry */

  mother1[ 0 ] = (int16_t) ( number1 / m16Long );
  mother2[ 0 ] = (int16_t) ( number2 / m16Long );

  /* Put the low bits of numberi into motheri[1] */

  mother1[ 1 ] = (int16_t) ( m16Mask & number1 );
  mother2[ 1 ] = (int16_t) ( m16Mask & number2 );

  /* Combine the two 16 bit random numbers into one 32 bit */

  *pSeed = ( ( (int32_t) mother1[ 1 ] ) << 16 ) + (int32_t) mother2[ 1 ];

  /* Return a double value between 0 and 1 */

  return ( (double) *pSeed ) / m32Double;
}

/*
 * Gaussian random variable.
 * Reference: Devroye, 194-199.  This is a ratio-of-uniforms method with
 * quick acceptance and quick rejection conditions.  The constant E22
 * must be not less than sqrt(8/e), and the expected number of uniforms
 * per gaussian is the constant times 4/sqrt(2*pi); taking the constant to
 * be exactly sqrt(8/e), this gives 8/sqrt(pi*e) or about 2.738 uniforms
 * per gaussian.
 */
static double localfn_random_normal_deviate_marsaglia(
  uint32_t *pSeed )
{

  double rnormk, u, x2;

  do {
    u = localfn_random_uniform_deviate_marsaglia( pSeed );
    rnormk = E22 * ( localfn_random_uniform_deviate_marsaglia( pSeed )-0.5) / u;
    x2 = rnormk * rnormk / 4;
  } while(x2 > 1-u && x2 > -log(u));

  return( rnormk );
}

