
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/wav_boot.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "wav_boot.h"
#include "wav_type.h"
#include "wav_dwtc.h"

#include "mat_io.h"
#include "mat_arit.h"
#include "mat_assn.h"
#include "mat_rand.h"
#include "mat_set.h"
#include "mat_sort.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "mth_dist.h"
#include "mth_var.h"

#include "sig_tran.h"

#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_limit.h"
#include "ut_math.h"
#include "ut_mem.h"
#include "ut_rand.h"

#include <math.h>


/* Static functions declared here and defined at end of file */

static mutil_errcode localfn_seek_whitest_basis(
  const mat_set *dwpt,
  const double   significance,
  const sint32   level,
  const sint32   osc,
  const sint32   max_level,
  const wav_white_test  white_noise_test,
  sint32        *n_white_crystal,
  void          *intrp_ptr,
  univ_mat      *result );

static mutil_errcode localdef_acf(
  const univ_mat *time_series,
  const boolean   recenter,
  void           *intrp_ptr,
  univ_mat       *result );

static double localdef_box_pierce_statistic(
  const univ_mat *acf,
  const sint32    n_term );

static double localdef_ljung_box_pierce_statistic(
  const univ_mat *acf,
  const sint32    n_term );

static mutil_errcode localdef_normalized_cumulative_periodogram_statistic(
  const univ_mat *dft,
  double         *result );

static mutil_errcode localfn_shuffle_dwpt_basis(
  const mat_set *dwpt_basis,
  void          *rand_ptr,
  void          *intrp_ptr,
  mat_set       *dwpt_basis_random );

static mutil_errcode localdef_test_dwpt_crystal_whiteness(
  const univ_mat       *crystal,
  const wav_white_test  white_noise_test,
  const double          significance,
  void                 *intrp_ptr,
  boolean              *is_white );

static mutil_errcode localfn_filters_check(
  const mat_set *filters );

/* Static macro definitions */

#undef LOCALDEF_ILOG2
#define LOCALDEF_ILOG2( VALUE ) \
(sint32) floor( MUTIL_LOG2( (double) ( VALUE ) + MUTIL_DOUBLE_EPSILON ) )

#undef LOCALDEF_IS_ODD
#define LOCALDEF_IS_ODD(N) ( (N % 2) == 1 ? 1 : 0 )

#undef LOCALDEF_SEQUENCY2FLAT
#define LOCALDEF_SEQUENCY2FLAT( LEVEL, OSCILLATION_INDEX ) \
  ( 1 << ( LEVEL ) ) - 1 + ( OSCILLATION_INDEX )

#undef LOCALDEF_FLAT2SEQUENCY
#define LOCALDEF_FLAT2SEQUENCY( FLAT_INDEX,  LEVEL, OSCILLATION_INDEX )  \
( LEVEL ) = LOCALDEF_ILOG2( ( FLAT_INDEX ) + 1 ) ); \
( OSCILLATION_INDEX ) = ( FLAT_INDEX ) + 1 - ( 1 << ( LEVEL ) )

#undef LOCALDEF_CHECK_NULL_MATRIX_POINTER
#define LOCALDEF_CHECK_NULL_MATRIX_POINTER( DATA_PTR, DATA_TYPE,  \
                                     TYPE_PREFIX )                \
   err = TYPE_PREFIX ## _validate( DATA_PTR );                    \
   if ( err ) return err;                                         \
   if ( DATA_PTR == ( DATA_TYPE * ) NULL ){                       \
     MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" );    \
     return MUTIL_ERR_NULL_POINTER;                               \
   }

/* Adaptive bootstrapping using the DWPT */
/* Documented in wav_boot.h              */
/* Written by William Constantine        */
mutil_errcode wavuniv_bootstrap(
  const mat_set    *dwpt,
  const mat_set    *filters,
  const univ_mat   *white_indices,
  const sint32      n_realization,
  void             *intrp_ptr,
  mat_set          *result )
{
  double_mat     *pd_dwpt;
  mat_set         white_basis;
  mat_set         white_basis_random;
  memlist         list;
  mutil_errcode   err;
  sint32          dims;
  sint32          i;
  sint32_mat      ncol;
  sint32_mat      nrow;
  univ_mat        synthesis;
  void           *rand_ptr;
  wav_dwpt_extra  extra;

  MUTIL_TRACE( "Start wavuniv_bootstrap()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs ... */

  /*** check dwpt matrix set */

  LOCALDEF_CHECK_NULL_MATRIX_POINTER( dwpt, mat_set, matset );

  if ( dwpt->mats[ 0 ].type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "DWPT matrix set must contain universal matrices of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( !LOCALDEF_IS_ODD( dwpt->nelem ) ){
    MUTIL_ERROR( "DWPT matrix set must contain an odd number of crystals" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /*** check filters ***/

  err = localfn_filters_check( filters );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /*** check white_indices ***/

  LOCALDEF_CHECK_NULL_MATRIX_POINTER( white_indices, univ_mat, matuniv );

  if ( white_indices->type != MUTIL_SINT32 ){
    MUTIL_ERROR( "Transform indices matrix must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /*** check n_realization */

  if ( n_realization < 1 ){
    MUTIL_ERROR("The number of realizations must be positive");
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* obtain the specified white transform
     (which is a subset of the original
     dwpt matrix set) and register the result
     with the memory manager */

  err = wavuniv_transform_packet_basis( dwpt, white_indices, intrp_ptr, &white_basis, &extra );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &white_basis, MEMTYPE_MATSET );
  if ( err ){
    MUTIL_FREEALL_MATSET_WARN( &white_basis );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* allocate memory for a copy of the white basis */

  dims = white_basis.nelem;

  err = mats32_malloc_register( &nrow, 1, dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &ncol, 1, dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  for ( i = 0; i < dims; i++ ){

    pd_dwpt = &( white_basis.mats[ i ].mat.dblmat );

    nrow.data[ i ] = pd_dwpt->nrow;
    ncol.data[ i ] = pd_dwpt->ncol;
  }

  err = matset_malloc_register( &white_basis_random, 1, &dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices_arbitrary_size(
    &white_basis_random, &nrow, &ncol, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* allocate memory for the result */

  err = matset_malloc_register( result, 1, &n_realization, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* initiate random number generation */

  err = mutil_rand_begin( &rand_ptr );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mutil_rand_set_seed( (void*) NULL, rand_ptr );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create the surrogate data */

  for ( i = 0; i < n_realization; i++ ){

    /* randomize the DWPT basis by randomly
       shuffling the coefficients within each
       crystal with replacement */

    err = localfn_shuffle_dwpt_basis( &white_basis, rand_ptr, intrp_ptr,
      &white_basis_random );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* now perform synthesis and store result */

    err = wavuniv_transform_packet_inverse( &white_basis_random, &extra,
      white_indices, filters, intrp_ptr, &synthesis );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matuniv_wrap_univ_matrix( &( result->mats[ i ] ), &synthesis );
    if ( err ){
      MUTIL_FREE_WARN( matuniv, &synthesis );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }
  }

  /* end random number generation */

  err = mutil_rand_end( rand_ptr );
  if ( err ) {
    MUTIL_ERROR( "Problem ending random number generator" );
    return err;
  }

  /* free extra DWPT atom structure */

  if ( extra.nelem > 0 ){

    MUTIL_FREE_WARN( matdbl, &( extra.atoms ) );
    MUTIL_FREE_WARN( mats32, &( extra.levelmap ) );
  }

  /* free nodes corresponding to registered
     modwpt matrix memory, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free memory registered with the memory manager */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_bootstrap()" );

  return MUTIL_ERR_OK;
}


/* Finding the whitest set of DWPT crystals */
/* Documented in wav_boot.h                 */
/* Written by William Constantine           */

mutil_errcode wavuniv_transform_packet_whitest(
  const mat_set  *dwpt,
  const double    significance,
  const wav_white_test  white_noise_test,
  void           *intrp_ptr,
  univ_mat       *result )
{
  mutil_errcode  err;
  memlist        list;
  sint32         n_white_crystal;
  sint32         n_level;
  univ_mat       temp;
  univ_mat       flat;
  univ_mat       level;
  univ_mat       osc;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_transform_packet_whitest()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check input ... */

  /* ... dwpt matrix set */

  LOCALDEF_CHECK_NULL_MATRIX_POINTER( dwpt, mat_set, matset );

  if ( dwpt->mats[ 0 ].type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "DWPT matrix set must contain universal matrices of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( !LOCALDEF_IS_ODD( dwpt->nelem ) ){
    MUTIL_ERROR( "DWPT matrix set must contain an odd number of crystals" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... significance */

  if ( significance <= 0.0 || significance >= 1 ){
    MUTIL_ERROR( "significance must be on (0,1)" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* ... white noise test type */

  switch( white_noise_test ){
    case WAV_WHITE_TEST_PORTMANTEAU_I:
    case WAV_WHITE_TEST_PORTMANTEAU_II:
    case WAV_WHITE_TEST_PORTMANTEAU_III:
    case WAV_WHITE_TEST_CUMULATIVE_PERIODOGRAM:
      break;
    default:
      MUTIL_ERROR( "White noise test is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* initialize variables */

  n_level = LOCALDEF_ILOG2( dwpt->nelem );

  n_white_crystal = 0;

  /* allocate memory */

  err = matuniv_malloc_register( &temp, 2, (sint32) MUTIL_POW( 2, n_level ),
    MUTIL_SINT32, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* find the whitest basis */

  err = localfn_seek_whitest_basis(
    dwpt,
    significance,
    0, 0,
    n_level,
    white_noise_test,
    &n_white_crystal,
    intrp_ptr,
    &temp );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( MUTIL_INTERRUPT( 3.0 * n_white_crystal, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }

  /* reduce the number of columns in the result matrix
     to the number of crystals which passed the whiteness test ... */

  err = matuniv_malloc_register( result, 2, n_white_crystal, MUTIL_SINT32, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_extract( &temp, 0, 0, intrp_ptr, result );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* sort the indices */

  err = wavuniv_transform_packet_convert_indices(
    result, intrp_ptr, &flat, &level, &osc );
  MEMLIST_FREE_ON_ERROR( err, &list );

  MUTIL_FREE_WARN( matuniv, &flat );

  err = memlist_member_unregister( result, &list );
  if ( err ){
    MUTIL_FREE_WARN( matuniv, &level );
    MUTIL_FREE_WARN( matuniv, &osc );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  err = matuniv_assign_submat( &level, 0, 0, intrp_ptr, result );
  if ( err ){
    MUTIL_FREE_WARN( matuniv, &level );
    MUTIL_FREE_WARN( matuniv, &osc );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  err = matuniv_assign_submat( &osc, 1, 0, intrp_ptr, result );
  if ( err ){
    MUTIL_FREE_WARN( matuniv, &level );
    MUTIL_FREE_WARN( matuniv, &osc );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  MUTIL_FREE_WARN( matuniv, &level );
  MUTIL_FREE_WARN( matuniv, &osc );

  /* free malloced space and corresponding linked list structures */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_transform_packet_whitest()" );

  return MUTIL_ERR_OK;
}

/** Tests a single discrete wavelet packet crystal for whiteness
 * based on one of a variety of tests slected by the user.
 *
 * References:
 *
 * 1. D. B. Percival, S. Sardy and A. C. Davison, {\it Wavestrapping Time Series:
 * Adaptive Wavelet-Based Bootstrapping}, In W. J. Fitzgerald, R. L. Smith,
 * A. T. Walden and P. C. Young (Eds.), Nonlinear and Nonstationary Signal
 * Processing, Cambridge, England: Cambridge University Press.
 *
 * 2. Stephens, M.A. (1974), 'EDF Statistics for Goodness of Fit and
 * Some Comparisons', Journal of the American Statistical Association,
 * 69, 730-737.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_boot.h
 * @source wav\_boot.c
 * @library wavelets
 * @usage #err = localdef_test_dwpt_crystal_whiteness( &crystal, white_noise_test, significance, intrp_ptr, &is_white );#
 * @return Standard mutils error/OK code.
 * @param  crystal    Pointer to a pre-allocated universal matrix
 *                    of type MUTIL\_DOUBLE containing a DWPT crystal.
 * @param  white_noise_test  An enum of type \Ref{_wav_white_test}
 *                    denoting the white noise test to use.
 * @param  significance The significance. The white noise test involves the calculation
 *                    of a white noise test statistic which is compared
 *                    to to a chi-square (1 - significance) x 100 perctage point
 *                    with a (pre-calculated) number of degrees of
 *                    freedom..
 * @param  intrp_ptr  Pointer for implementation of interrupt checking.
 * @param  is_white   Pointer to a boolean value to contain the result.
 * @see localfn_seek_whitest_basis
 * @see wavuniv_transform_packet
 * @see localdef_box_pierce_statistic
 * @see localdef_ljung_box_pierce_statistic
 * @see localdef_normalized_cumulative_periodogram_statistic
 * @private
 */
static mutil_errcode localdef_test_dwpt_crystal_whiteness(
  const univ_mat       *crystal,
  const wav_white_test  white_noise_test,
  const double          significance,
  void                 *intrp_ptr,
  boolean              *is_white )
{
  memlist        list;
  mutil_errcode  err;
  univ_mat       statistic;
  univ_mat       crystal_squared;
  univ_mat       dft_series;
  sint32         nrow;
  sint32         ncol;
  sint32         n_sample;
  double         test_value;
  double         hypothetical_value;
  sint32         n_lag;
  double         Calpha;
  double         fac;
  univ_mat       crystal_column;

  MUTIL_TRACE( "Start localdef_test_dwpt_crystal_whiteness()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs */

  /* ... if series is a valid universal matrix */

  LOCALDEF_CHECK_NULL_MATRIX_POINTER( crystal, univ_mat, matuniv  );

  /* ... if the series type is double */

  if ( crystal->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "DWPT crystal matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... if it is a vector */

  if ( !MATANY_IS_VEC( &(crystal->mat.dblmat) ) ){
    MUTIL_ERROR( "DWPT crystal matrix must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* check whiteness test type */

  switch( white_noise_test ){
    case WAV_WHITE_TEST_PORTMANTEAU_I:
    case WAV_WHITE_TEST_PORTMANTEAU_II:
    case WAV_WHITE_TEST_PORTMANTEAU_III:
    case WAV_WHITE_TEST_CUMULATIVE_PERIODOGRAM:
      break;
    default:
      MUTIL_ERROR( "Whiteness test type is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* check significance */

  if ( significance < 0 || significance > 1 ){
    MUTIL_ERROR( "Significnace must be on the interval [0,1]" );
  }

  /* initialize varaibles */

  nrow     = MATUNIV_NROW( crystal );
  ncol     = MATUNIV_NCOL( crystal );
  n_sample = MATUNIV_NELEM( crystal );
  n_lag    = MUTIL_MAX( 2, MUTIL_MIN( 20, (sint32) ( n_sample / 10 ) ) );

  /* form the test statistic */

  switch( white_noise_test ){

    case WAV_WHITE_TEST_PORTMANTEAU_I:

      err = localdef_acf( crystal, FALSE, intrp_ptr, &statistic );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_register( &list, &statistic, MEMTYPE_MATUNIV );
      if ( err ){
	MUTIL_FREE_WARN( matuniv, &statistic );
	MEMLIST_FREE_ON_ERROR( err, &list );
      }

      test_value = localdef_box_pierce_statistic( &statistic, n_lag );

      break;

    case WAV_WHITE_TEST_PORTMANTEAU_II:

      err = localdef_acf( crystal, FALSE, intrp_ptr, &statistic );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_register( &list, &statistic, MEMTYPE_MATUNIV );
      if ( err ){
	MUTIL_FREE_WARN( matuniv, &statistic );
	MEMLIST_FREE_ON_ERROR( err, &list );
      }

      test_value = localdef_ljung_box_pierce_statistic( &statistic, n_lag );

      break;

    case WAV_WHITE_TEST_PORTMANTEAU_III:

      err = matuniv_malloc_register( &crystal_squared, nrow, ncol, MUTIL_DOUBLE, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = matuniv_multiply_elem( crystal, crystal, intrp_ptr, &crystal_squared );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = localdef_acf( &crystal_squared, TRUE, intrp_ptr, &statistic );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_register( &list, &statistic, MEMTYPE_MATUNIV );
      if ( err ){
	MUTIL_FREE_WARN( matuniv, &statistic );
	MEMLIST_FREE_ON_ERROR( err, &list );
      }

      err = memlist_member_free( &crystal_squared, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      test_value = localdef_ljung_box_pierce_statistic( &statistic, n_lag );

      break;

    case WAV_WHITE_TEST_CUMULATIVE_PERIODOGRAM:

      err = matuniv_malloc_register( &statistic, nrow, ncol, MUTIL_DOUBLE, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = matuniv_malloc_register( &dft_series, n_sample, 1, MUTIL_DCOMPLEX, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* trick the DFT into thinking the current crystal is a column-vector,
	 regardless of whether it is or not. the DFT expects all data to be
	 in column form */

      crystal_column.type = MUTIL_DOUBLE;
      crystal_column.mat.dblmat.nelem = n_sample;
      crystal_column.mat.dblmat.nrow  = n_sample;
      crystal_column.mat.dblmat.ncol  = (sint32) 1;
      crystal_column.mat.dblmat.data  = crystal->mat.dblmat.data;

      err = siguniv_transform_discrete_fourier( &crystal_column, FALSE, intrp_ptr, &dft_series );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = localdef_normalized_cumulative_periodogram_statistic( &dft_series, &test_value );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_free( &dft_series, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      break;

    default:
      MUTIL_ERROR( "Whiteness test type is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* form hypothetical value aligned with the NULL hypothesis of white noise */

  if ( white_noise_test == WAV_WHITE_TEST_CUMULATIVE_PERIODOGRAM ){

    /* use Stephens' Monte Carlo approximations for the distribution
       of D(alpha) statistic */

    if ( significance == 0.1 ){
      Calpha = 1.224;
    }
    else if ( significance == 0.05 ){
      Calpha = 1.358;
    }
    else if ( significance == 0.01 ){
      Calpha = 1.628;
    }
    else{
      MUTIL_ERROR( "Significance must be 0.1, 0.05, or 0.01 for cumulative periodogram test" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
    }

    fac = sqrt( (double) ( n_sample / 2 - 2 ) ); /* integer division okay here */

    hypothetical_value = Calpha / ( fac + 0.12 + 0.11 / fac );

  }
  else{

    hypothetical_value = mth_qchisq( 1.0 - significance, (double) n_lag );
  }

  *is_white = ( test_value < hypothetical_value ) ? TRUE : FALSE;

  /* free malloced space and corresponding linked list structures */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localdef_test_dwpt_crystal_whiteness()" );

  return MUTIL_ERR_OK;
}


/** Calculates the single-sided non-recentered biased autocorrelation sequence
 * for an input time series.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_boot.h
 * @source wav\_boot.c
 * @library wavelets
 * @usage #err = localdef_acf( &time_series, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param  time_series Pointer to a pre-allocated single-column or single-row
 *                    universal matrix of type MUTIL\_DOUBLE containing a time
 *                    series (or DWPT crystal).
 * @param  recenter   Boolean. If TRUE, the sample mean is first subtracted from the series.
 * @param  intrp_ptr  Pointer for implementation of interrupt checking.
 * @param  result     Pointer to a universal matrix of type MUTIL\_DOUBLE
 *                    containing the resulting acf. The memory for the
 *                    data is automatically allocated wihtin the function.
 * @see localdef_test_dwpt_crystal_whiteness
 * @private
 */
static mutil_errcode localdef_acf(
  const univ_mat *time_series,
  const boolean   recenter,
  void           *intrp_ptr,
  univ_mat       *result )
{
  mutil_errcode  err;
  boolean        biased = TRUE;

  /* calculate the acvs */

  err = mthuniv_acvs( time_series, biased, recenter, intrp_ptr, result );
  if ( err ) return err;

  /* normalize the acvs by its first element (representing the
     acvs at lag zero) to form the sample autocorrelation sequence */

  err = matdbl_divide_scalar( &( result->mat.dblmat ), result->mat.dblmat.data[ 0 ],
    TRUE, intrp_ptr, &( result->mat.dblmat ) );
  if ( err ) return err;

  return MUTIL_ERR_OK;
}


/** The Box-Pierce (Portmanteau I) whiteness test.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_boot.h
 * @source wav\_boot.c
 * @library wavelets
 * @usage #result = localdef_box_pierce_statistic( &acf, n_term );#
 * @return The whiteness statistic.
 * @param  acf Pointer to a pre-allocated single-column or single-row
 *             universal matrix of type MUTIL\_DOUBLE containing an
 *             autocorrelation sequence for a DWPT crystal.
 * @param  n_term The number of terms of the acf sequence to use
 *                in calculating the statistic.No check is made to
 *                ensure that n\_term is less than the number of
 *                elements in the acf matrix.
 * @see localdef_test_dwpt_crystal_whiteness
 * @see localdef_ljung_box_pierce_statistic
 * @see localdef_normalized_cumulative_periodogram_statistic
 * @private
 */
static double localdef_box_pierce_statistic(
  const univ_mat *acf,
  const sint32    n_term )
{
  double  result   = 0.0;
  double *pd_acf   = acf->mat.dblmat.data + 1;
  sint32  n_sample = MATUNIV_NELEM( acf );
  sint32  t;

  for ( t = 1; t <= n_term; t++ ){

    result += MUTIL_SQR( *pd_acf );

    pd_acf++;
  }

  result *= (double) n_sample;

  return result;
}

/** The Ljung-Box-Pierce (Portmanteau II, III) whiteness test.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_boot.h
 * @source wav\_boot.c
 * @library wavelets
 * @usage #result = localdef_ljung_box_pierce_statistic( &acf, n_term );#
 * @return The whiteness statistic.
 * @param  acf Pointer to a pre-allocated single-column or single-row
 *             universal matrix of type MUTIL\_DOUBLE containing an
 *             autocorrelation sequence for a DWPT crystal.
 * @param  n_term The number of terms of the acf sequence to use
 *                in calculating the statistic.No check is made to
 *                ensure that n\_term is less than the number of
 *                elements in the acf matrix.
 * @see localdef_test_dwpt_crystal_whiteness
 * @see localdef_box_pierce_statistic
 * @see localdef_normalized_cumulative_periodogram_statistic
 * @private
 */
static double localdef_ljung_box_pierce_statistic(
  const univ_mat *acf,
  const sint32    n_term )
{
  double  result   = 0.0;
  double *pd_acf   = acf->mat.dblmat.data + 1;
  sint32  n_sample = MATUNIV_NELEM( acf );
  sint32  t;

  for ( t = 1; t <= MUTIL_MIN( n_sample - 1, n_term ); t++ ){

    result += MUTIL_SQR( *pd_acf ) / (double) ( n_sample - t );

    pd_acf++;
  }

  result *= (double) ( n_sample * ( n_sample + 2 ) );

  return result;
}

/** The normalized cumulative periodogram whiteness statistic.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_boot.h
 * @source wav\_boot.c
 * @library wavelets
 * @usage #err = localdef_normalized_cumulative_periodogram_statistic( &dft, &result );#
 * @return Standard mutils error/OK code.
 * @param  dft Pointer to a pre-allocated single-column or single-row
 *             universal matrix of type MUTIL\_DCOMPLEX containing the
 *             DFT coefficients for a DWPT crystal.
 * @param  result Pointer to a double value containing the whiteness test statistic.
 * @see localdef_test_dwpt_crystal_whiteness
 * @see localdef_box_pierce_statistic
 * @see localdef_ljung_box_pierce_statistic
 * @private
 */
static mutil_errcode localdef_normalized_cumulative_periodogram_statistic(
  const univ_mat *dft,
  double         *result )
{
  dcomplex      *pc_dft;
  double         fac;
  double         norm;
  double         val;
  double        *pd_dft2;
  double_mat     dft_squared;
  memlist        list;
  mutil_errcode  err;
  sint32         Mj;
  sint32         i;
  sint32         n_sample;
  double         cumsum;
  double         P;

  /* initialize memory list */

  MEMLIST_INIT( list );

  MUTIL_TRACE( "Start localdef_normalized_cumulative_periodogram_statistic()" );

  /* initialize variables */

  n_sample = MATUNIV_NELEM( dft );
  Mj       = n_sample / 2 - 1;

  norm    = 0.0;
  *result = 0.0;
  fac     = (double) ( Mj - 1 );

  /* allocate  memory */

  err = matdbl_malloc_register( &dft_squared, Mj, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* assign pointers */

  pd_dft2 = dft_squared.data;
  pc_dft  = dft->mat.cpxmat.data + 1;

  /* form squared modulus of dft coefficients */

  for ( i = 0; i < Mj; i++ ){

    val = MUTIL_CPX_ABS( *pc_dft );

    *pd_dft2 = val * val;

    norm += *pd_dft2;

    pd_dft2++;
    pc_dft++;
  }

  pd_dft2 = dft_squared.data;
  cumsum = 0.0;

  for ( i = 1; i <= Mj - 1; i++ ){

    cumsum += *pd_dft2;

    P = cumsum / norm;

    /* form D+ */

    val = (double) i / fac - P;
    if ( val > *result ) *result = val;

    /* form D- */

    val = P - (double) ( i - 1 ) / fac;
    if ( val > *result ) *result = val;

    /* increment pointer */

    pd_dft2++;
  }

  /* free malloced space and corresponding linked list structures */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localdef_normalized_cumulative_periodogram_statistic()" );

  return MUTIL_ERR_OK;
}


/** Seeks the whitest transform of many possible transforms
 * in a discrete wavelet packet transform. This function
 * is recursive and performs the following algorithm:
 * (1) a parent crystal is tested for whiteness.
 * (2) if the parent is not white, then the two children
 * are tested (here, the localfn\_seek\_whitest\_basis
 * function is called recursively, once for each child).
 * in this sceanario, the recursive structure will allow the
 * program to seek the whitest basis via a binary tree. If
 * a crystal is at the last level of the transform, then it is
 * automatically added as a member of the DWPT subset comprising
 * the whitest transform. The idea here is that children in the
 * the last level are at least as white as the parent, i.e.,
 * dividing the corresponding band in the SDF will result in
 * two bands which are at least as white in each as the parent.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_boot.h
 * @source wav\_boot.c
 * @library wavelets
 * @usage #err = localfn_seek_whitest_basis( &dwpt, significance, level, osc, max_level, white_noise_test, n_white_crystal, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param  dwpt Pointer to a pre-allocated matrix set of type MUTIL\_DOUBLE
 *              containing a discrete wavelet packet transform
 *              such as that retuned by \Ref{wavuniv_transform_packet}.
 * @param  significance  The significance level to use in calculating comparative chi-square distribution
 *                     p x 100 percentage points where p = 1 - significance (the chi-square degrees of
 *                     freedom are estimated automatically within the specified white noise test).
 * @param level The level of the current crystal being analyzed.
 * @param osc The oscillation index of the current crystal being analyzed.
 * @param max_level The maximum level to consider in seeking the whitest transform.
 * @param  white_noise_test An enum of type \Ref{_wav_white_test} denoting the white noise
 *                     test to use in testing wavelet packet crystals.
 * @param n_white_crystal The current count of the number of white crystals
 * found in the transform.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  result      Pointer to a pre-allocated two-row universal matrix of
 *                     type MUTIL\_SINT32 which upon return will contain
 *                     the sequency ordered coordinates of the whitest set of
 *                     crystals in the DWPT bssed on the supplied
 *                     white\_noise\_test.test. The rows of this output
 *                     contain the decomposition levels and the oscillation
 *                     indices for the whitest set. Each column thus defines
 *                     the location of one such crystal in the DWPT in sequency
 *                     order. The result matrix must be large enough to store the
 *                     largest possible transform, which is equivalent to the
 *                     number of crystals in the max\_level of the DWPT.
 * @see wavuniv_transform_packet_whitest
 * @see wavuniv_transform_packet
 * @see localdef_box_pierce_statistic
 * @see localdef_ljung_box_pierce_statistic
 * @see localdef_normalized_cumulative_periodogram_statistic
 * @private
 */
static mutil_errcode localfn_seek_whitest_basis(
  const mat_set *dwpt,
  const double   significance,
  const sint32   level,
  const sint32   osc,
  const sint32   max_level,
  const wav_white_test  white_noise_test,
  sint32        *n_white_crystal,
  void          *intrp_ptr,
  univ_mat      *result )
{
  boolean        is_white = FALSE;
  mutil_errcode  err;
  sint32         i;
  sint32         index;
  sint32         ncol = MATUNIV_NCOL( result );
  sint32        *white_level;
  sint32        *white_osc;

  MUTIL_TRACE( "Start localfn_seek_whitest_basis()" );

  /* perform perfunctary error checks */

  if ( *n_white_crystal > MATUNIV_NCOL( result ) ){

    MUTIL_ERROR("The current count of white crystals exceeds the number " \
      "of columns in result matrix" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* perform whiteness test on specified crystal */

  if ( level < max_level ){

    index = LOCALDEF_SEQUENCY2FLAT( level, osc );

    err = localdef_test_dwpt_crystal_whiteness(
      &( dwpt->mats[ index ] ), white_noise_test,  significance,
      intrp_ptr, &is_white );
    if ( err ){
      MUTIL_ERROR("Problem in testing DWPT crystal for whiteness");
      return err;
    }
  }

  if ( is_white || level == max_level ){

    /* assign pointers */

    white_level = result->mat.s32mat.data;
    white_osc   = white_level + ncol;

    /* calculate column of result matrix to store test result */

    white_level[ *n_white_crystal ] = level;
    white_osc[ *n_white_crystal ]   = osc;

    /* update the count of crystals that have is_whiteed the test */

    (*n_white_crystal)++;

    return MUTIL_ERR_OK;
  }

  if ( !is_white && level < max_level ){

    for ( i = 0; i <= 1; i++ ){

      err = localfn_seek_whitest_basis(
	dwpt,
	significance,
	level + 1,
	2 * osc + i,
	max_level,
	white_noise_test,
	n_white_crystal,
	intrp_ptr,
	result );
      if ( err ) return err;
    }
  }

  MUTIL_TRACE( "Done with localfn_seek_whitest_basis()" );

  return MUTIL_ERR_OK;
}


/** Shuffles each crystal's coefficients (with replacement)
 * given an input DWPT basis.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_boot.h
 * @source wav\_boot.c
 * @library wavelets
 * @usage #err = localfn_shuffle_dwpt_basis( &dwpt_basis, rand_ptr, intrp_ptr, &dwpt_basis_random );#
 * @return Standard mutils error/OK code.
 * @param  dwpt Pointer to a pre-allocated matrix set of type MUTIL\_DOUBLE
 *              containing one transform (of the many possible) in a discrete
 *              wavelet packet transform such as that retuned by
 *              \Ref{wavuniv_transform_packet}. The transform must be a projection
 *              onto a legitimate basis, i.e., one in which the union of all
 *              passbands corresponding to each crystal of the DWPT subset
 *              must span the normalized frequency range on [ 0, 1/2].
 * @param  rand_ptr  Initialized pointer for implementing random number generator.
 * @param  intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  dwpt_basis_random  Pointer to a pre-allocated matrix set whose
 *                     set of universal matrices are of the same type and dimension
 *                     as that of the input matrix set. Upon return, this
 *                     matrix set will contain the shuffled coefficients
 *                     of the input DWPT basis.
 * @see wavuniv_bootstrap
 * @private
 */
static mutil_errcode localfn_shuffle_dwpt_basis(
  const mat_set *dwpt_basis,
  void          *rand_ptr,
  void          *intrp_ptr,
  mat_set       *dwpt_basis_random )
{
  sint32         j;
  mutil_errcode  err;
  memlist        list;
  univ_mat      *pum_dwpt_basis;
  univ_mat      *pum_dwpt_basis_random;
  boolean        replacement = TRUE;
  univ_mat       random_index;

  MUTIL_TRACE( "Start localfn_shuffle_dwpt_basis()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  for ( j = 0; j < dwpt_basis->nelem; j++ ){

    /* assign pointers */

    pum_dwpt_basis = &( dwpt_basis->mats[ j ] );
    pum_dwpt_basis_random = &( dwpt_basis_random->mats[ j ] );

    /* create randomized index vector (with replacement)
       and register with the memory manager */

    err = matuniv_random_uniform_indices(
      MATUNIV_NROW( pum_dwpt_basis ), MATUNIV_NCOL( pum_dwpt_basis ),
      replacement, rand_ptr, intrp_ptr, &random_index );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_register( &list, &random_index, MEMTYPE_MATUNIV );
    if ( err ){
      MUTIL_FREE_WARN( matuniv, &random_index );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    /* permute the original crystal based on the random index
       and the store the result in the corresponding location
       of the matrix set copy */

    err = matuniv_permute( pum_dwpt_basis, &random_index, intrp_ptr,
      pum_dwpt_basis_random );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* free the random index memory */

    err = memlist_member_free( &random_index, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* free malloced space and corresponding linked list structures */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_shuffle_dwpt_basis()" );

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

   LOCALDEF_CHECK_NULL_MATRIX_POINTER( filters, mat_set, matset );

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
