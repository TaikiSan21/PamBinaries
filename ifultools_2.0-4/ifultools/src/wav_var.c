
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/wav_var.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "wav_var.h"

#include "mat_assn.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "mat_io.h"

#include "mth_dist.h"
#include "mth_var.h"

#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_math.h"
#include "ut_mem.h"

#include "wav_dwtc.h"
#include "wav_filt.h"
#include "wav_math.h"
#include "wav_modw.h"
#include "wav_type.h"
#include <math.h>

/* Static functions declared here and defined at end of file */

static mutil_errcode localfn_wavuniv_variance_confidence_input_check(
  const univ_mat *variance,
  const univ_mat *edof,
  double          probability);

static mutil_errcode localfn_wavuniv_variance_edof_input_check(
  const univ_mat *interior,
  const univ_mat *n_coeff,
  const univ_mat *variance,
  const univ_mat *level,
  const univ_mat *sdf);

static mutil_errcode localfn_wavuniv_variance_input_check(
 const univ_mat  *time_series,
 wav_transform    transform_type,
 wav_filter_type  filter_type,
 sint32           filter_length,
 sint32           n_level,
 const univ_mat  *sdf );

static mutil_errcode localfn_check_sdf(
  const univ_mat *sdf );

/* Static macro definitions */

#define LOCALDEF_CHECK_NULL_POINTER_VAR( DATA_PTR, DATA_TYPE,  \
                                     TYPE_PREFIX )             \
   err = TYPE_PREFIX ## _validate( DATA_PTR );                 \
   if ( err ) return err;                                      \
   if ( DATA_PTR == ( DATA_TYPE * ) NULL ){                    \
     MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
     return MUTIL_ERR_NULL_POINTER;                            \
   }

#define LOCALDEF_CHECK_UNIVMAT_VECTOR( MAT_PTR, MAT_TYPE, TYPE_PREFIX )        \
                                                                               \
 /* validate matrix structure and check for NULL pointer */                    \
                                                                               \
 LOCALDEF_CHECK_NULL_POINTER_VAR( MAT_PTR, univ_mat, matuniv );                \
                                                                               \
 /* check matrix type */                                                       \
                                                                               \
 if ( ( MAT_PTR )->type != MAT_TYPE ){                                         \
   MUTIL_ERROR( #MAT_PTR " matrix must be of type " #MAT_TYPE );               \
   return MUTIL_ERR_ILLEGAL_TYPE;                                              \
 }                                                                             \
                                                                               \
 /* ... to see if it is a vector */                                            \
                                                                               \
 if ( !MATANY_IS_VEC( &( ( MAT_PTR )->mat. TYPE_PREFIX ) ) ){                  \
   MUTIL_ERROR( #MAT_PTR " matrix must be a single column or row." );          \
   return MUTIL_ERR_ILLEGAL_SIZE;                                              \
 }                                                                             \
                                                                               \
 /* ... if number of elements is positive */                                   \
                                                                               \
 if ( MATUNIV_NELEM( MAT_PTR ) < 1 ){                                          \
   MUTIL_ERROR("Number of elements in " #MAT_PTR " matrix must be positive."); \
   return MUTIL_ERR_ILLEGAL_SIZE;                                              \
 }

#define CALCULATE_CONFIDENCE( EDOF_ROW )                         \
  err = wavuniv_variance_confidence(                             \
    &variance_block->mats[ 1 ],                                  \
    &edof->mats[ EDOF_ROW ],                                     \
    ( double ) 0.975,                                            \
    intrp_ptr,                                                   \
    &conf );                                                     \
  MEMLIST_FREE_ON_ERROR( err, &list );                           \
                                                                 \
  /* register the new matset with the memory manager */          \
                                                                 \
  err = memlist_member_register( &list, &conf, MEMTYPE_MATSET ); \
  MEMLIST_FREE_ON_ERROR( err, &list );                           \
                                                                 \
  /* map the current confidence interval matrix      */          \
  /* set of 2 [1 x J] vectors into a single [2 x J]  */          \
  /* matrix. Store the result in the appropriate     */          \
  /* location of the global matrix set for the	     */          \
  /* confidence intervals ...                        */          \
                                                                 \
  /* ... assign lower limits of confidence intervals */          \
                                                                 \
  err = matuniv_assign_submat(                                   \
    &conf.mats[ 0 ], 0, 0, intrp_ptr,                            \
    &confidence->mats[ EDOF_ROW ] );                             \
  MEMLIST_FREE_ON_ERROR( err, &list );                           \
                                                                 \
  /* ... assign upper limits of confidence intervals */          \
                                                                 \
  err = matuniv_assign_submat(                                   \
    &conf.mats[ 1 ], 1, 0, intrp_ptr,                            \
    &confidence->mats[ EDOF_ROW ] );                             \
  MEMLIST_FREE_ON_ERROR( err, &list );                           \
                                                                 \
  /* free the two vectors in the current confidence   */         \
  /* interval matrix set                              */         \
                                                                 \
  err = memlist_member_free( &conf, &list );                     \
  MEMLIST_FREE_ON_ERROR( err, &list )


/* Time (in)dependent (un)biased wavelet variance    */
/* and corresponding confidence interval estimation. */
/* Documented in wav_var.h                           */
/* Written by William Constantine                    */

mutil_errcode wavuniv_variance(
  const univ_mat  *time_series,
  wav_transform    transform_type,
  wav_filter_type  filter_type,
  sint32           filter_length,
  sint32           n_level,
  const univ_mat  *sdf,
  void            *intrp_ptr,
  mat_set         *variance_time,
  mat_set         *variance_block,
  mat_set         *confidence,
  mat_set         *edof )
{
  boolean        normalize_filters;
  boolean        calculate_edof;
  boolean        calculate_confidence;
  double         norm;
  double         num_ops = 0.0;
  double         sum_var_biased;
  double         sum_var_unbiased;
  double         var;
  double        *Wj;
  double        *pd_var_block_biased;
  double        *pd_var_block_unbiased;
  double        *pd_var_time_biased;
  double        *pd_var_time_unbiased;
  mat_set        boundary;
  mat_set        conf;
  mat_set        filters;
  memlist        list;
  mutil_errcode  err;
  sint32         dims_confidence = 3;
  sint32         dims_variance_block = 2;
  sint32         dims_variance_time = 2;
  sint32         ii;
  sint32         j;
  sint32         n_sample;
  sint32         ndim = 1;
  sint32         sum_Lj;
  sint32         t;
  sint32         t_interior;
  sint32        *Lj;
  sint32        *Nj;
  univ_mat       interior;
  univ_mat       levels;
  mat_set        transform;
  univ_mat       var_block_biased;
  univ_mat       var_block_unbiased;
  univ_mat       var_time_biased;
  univ_mat       var_time_unbiased;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_variance()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check input arguments */
  err = localfn_wavuniv_variance_input_check(
    time_series,
    transform_type,
    filter_type,
    filter_length,
    n_level,
    sdf );
  if ( err ) return err;

  /* set confidence interval related booleans */

  calculate_edof = ( edof != (mat_set *) NULL );
  calculate_confidence = ( confidence != (mat_set *) NULL );

  /* force the calculation of the EDOF if confidence intervals
  are desired */

  if ( calculate_confidence && !calculate_edof ){
    MUTIL_ERROR( "Confidence intervals cannot be calculated unless "
      "EDOF values are also allowed to be calculated" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* ensure that confidence intervals and EDOF are
  not calculated if the transform is not a MODWT */

  if ( ( transform_type != (wav_transform) WAV_TRANSFORM_MODWT ) &&
    ( calculate_confidence || calculate_edof ) ){
    MUTIL_ERROR( "Confidence intervals and corresponding EDOF "
      "values can only be calculated using the MODWT" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* initialize local variables */

  n_sample = MATUNIV_NELEM( time_series );

  /* allocate space for the output ... */

  /* ... create the matset headers */

  err = matset_malloc_register(
    variance_time,
    ndim,
    &dims_variance_time,
    &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_register(
    variance_block,
    ndim,
    &dims_variance_block,
    &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( calculate_confidence ){

    err = matset_malloc_register(
      confidence,
      ndim,
      &dims_confidence,
      &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matset_malloc_matrices(
      confidence,
      2,
      n_level,
      MUTIL_DOUBLE );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* ... allocate space for matset universal matrices */

  err = matset_malloc_matrices(
    variance_time,
    n_level,
    n_sample,
    MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices(
    variance_block,
    1,
    n_level,
    MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* ... time dependent wavelet variance matrices */

  /* obtain the wavelet and scaling filters
     and register them with the memory manager */


  if ( transform_type == (wav_transform) WAV_TRANSFORM_MODWT ){
    normalize_filters = (boolean) TRUE;
  }
  else{
    normalize_filters = (boolean) FALSE;
  }

  err = wavuniv_filters_daubechies(
    filter_length,
    filter_type,
    normalize_filters,
    intrp_ptr,
    &filters );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &filters, MEMTYPE_MATSET );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate a wavelet transform of the time series */

  if ( transform_type == WAV_TRANSFORM_MODWT ){

    err = wavuniv_transform_maximum_overlap(
      time_series,
      &filters,
      n_level,
      intrp_ptr,
      &transform);
    MEMLIST_FREE_ON_ERROR( err, &list );
  }
  else if ( transform_type == (wav_transform) WAV_TRANSFORM_DWT ){

    err = wavuniv_transform_discrete_wavelet_convolution(
      time_series,
      &filters,
      n_level,
      intrp_ptr,
      &transform );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }
  else{
    MUTIL_ERROR( "Wavelet transform type is currently unsupported" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* register the wavelet transform with the memory manager */

  err = memlist_member_register( &list, &transform, MEMTYPE_MATSET );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* Calculate the location of the boundary and interior
     wavelet coefficients. the call to the wavuniv_transform_coefficient_boundaries()
     function returns the following information in a matrix set:

     Matrix 0: smallest index of interior coefficients (index base 1)
     Matrix 1: largest index of interior coefficients (index base 1)
     Matrix 2: number of interior coefficients
     Matrix 3: number of boundary coefficients
     Matrix 4: number of all coefficients       */

  err = wavuniv_transform_coefficient_boundaries(
    n_level,
    filter_length,
    n_sample,
    transform_type,
    intrp_ptr,
    &boundary );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* register boundary matrix with memory manager */

  err = memlist_member_register( &list, &boundary, MEMTYPE_MATSET );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* allocate memory for confidence interval estimation */

  err = matuniv_malloc_register( &levels, 1, n_level,
    MUTIL_SINT32, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* allocate output memory */

  err = matuniv_malloc_register( &var_time_biased, n_level, n_sample,
    MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &var_time_unbiased, n_level, n_sample,
    MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &var_block_biased, 1, n_level,
    MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &var_block_unbiased, 1, n_level,
    MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* assign pointers */

  Lj = boundary.mats[ 2 ].mat.s32mat.data;
  Nj = boundary.mats[ 4 ].mat.s32mat.data;

  pd_var_block_biased   = variance_block->mats[ 0 ].mat.dblmat.data;
  pd_var_block_unbiased = variance_block->mats[ 1 ].mat.dblmat.data;

  /* initialize index variables */

  ii     = 0;
  sum_Lj = 0;

  /* calculate the total number of interior
     wavelet coefficients and fill in the
     levels vector. */

  for ( j = 0; j < n_level; j++ ){

    sum_Lj += Lj[ j ];

    levels.mat.s32mat.data[ j ] = j + 1;
  }

  /* allocate memory to store a concatenated
     collection of interior wavelet coefficients */

  err = matuniv_malloc_register( &interior, 1, sum_Lj,
    MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate the wavelet variance estimates */

  for ( j = 0; j < n_level; j++ ){

    /* assign pointer to current level's
       wavelet transform coefficients */

    Wj = transform.mats[ j ].mat.dblmat.data;

    /* set pointers in biased and unbiased universal matrices */

    pd_var_time_biased   =
      &( variance_time->mats[ 0 ].mat.dblmat.data[ j * n_sample ] );
    pd_var_time_unbiased =
      &( variance_time->mats[ 1 ].mat.dblmat.data[ j * n_sample ] );

    /* initialize variables */

    sum_var_biased   = ( double ) 0.0;
    sum_var_unbiased = ( double ) 0.0;
    t_interior       = Nj[ j ] - Lj[ j ];

    /* define normalization factor */

    if ( transform_type == WAV_TRANSFORM_DWT ){
      norm = ( double ) MUTIL_POW( 2.0, ( double ) ( j + 1 ) );
    }
    else norm = ( double ) 1.0;

    /* form time dependent wavelet variance estimates */

    for ( t = 0; t < Nj[ j ]; t++ ){

      var = Wj[ t ] * Wj[ t ] / norm;

      pd_var_time_biased[ t ] = var / ( double ) Nj[ j ];

      sum_var_biased += pd_var_time_biased[ t ];

      if ( t >= t_interior ){

        interior.mat.dblmat.data[ ii++ ] = Wj[ t ];

        pd_var_time_unbiased[ t ] = var / ( double ) Lj[ j ];

        sum_var_unbiased += pd_var_time_unbiased[ t ];
      }
      else{
        pd_var_time_unbiased[ t ] = ( double ) -1.0;
      }

      /* check for interrupts */

      num_ops += 10.0 * ( double ) Nj[ j ];
      if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
        MUTIL_ERROR( "user interrupt" );
        MUTIL_FREE_WARN( memlist, &list );
        return MUTIL_ERR_INTERRUPT;
      }

    } /* end of loop over time index */

    /* form blocked variance estimates */

    pd_var_block_biased[ j ] = sum_var_biased;

    if ( Lj[ j ] > 0 ){

      pd_var_block_unbiased[ j ] = sum_var_unbiased;
    }
    else{
      pd_var_block_unbiased[ j ] = ( double ) -1.0;
    }

  } /* end of loop over level index */

  if ( transform_type == (wav_transform) WAV_TRANSFORM_MODWT ){

    /* calculate the equivalent degrees of freedom
       for the confidence interval estimates
       and register output with the memory manager  */

    if ( calculate_edof ){

      err = wavuniv_variance_edof(
        &interior,
        &boundary.mats[ 2 ],
        &variance_block->mats[ 1 ],
        &levels,
        sdf,
        filter_type,
        filter_length,
        intrp_ptr,
        edof );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_register( &list, edof, MEMTYPE_MATSET );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    /* calculate the confidence intervals for each EDOF mode */

    if ( calculate_confidence ){
      CALCULATE_CONFIDENCE( 0 );
      CALCULATE_CONFIDENCE( 1 );
      CALCULATE_CONFIDENCE( 2 );
    }
  }

  /* remove nodes corresponding to
     output from the memory list */

  err = memlist_member_unregister( variance_time, &list);
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_unregister( variance_block, &list);
  if ( err ){
    MUTIL_FREEALL_MATSET_WARN( variance_time );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  if ( calculate_confidence ){
    err = memlist_member_unregister( confidence, &list);
    if ( err ){
      MUTIL_FREEALL_MATSET_WARN( variance_time );
      MUTIL_FREEALL_MATSET_WARN( variance_block );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }
  }

  if ( calculate_edof ){
    err = memlist_member_unregister( edof, &list);
    if ( err ){
      if ( calculate_confidence ){
	MUTIL_FREEALL_MATSET_WARN( confidence );
      }
      MUTIL_FREEALL_MATSET_WARN( variance_time );
      MUTIL_FREEALL_MATSET_WARN( variance_block );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }
  }

  /* free memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_variance()" );

  return MUTIL_ERR_OK;
}

/* Equivalent degress of freedom estimation for a */
/* chi-squared distribution assumption on the     */
/* interior wavelet coefficients                  */
/* Documented in wav_var.h                        */
/* Written by William Constantine                 */

mutil_errcode wavuniv_variance_edof(
  const univ_mat  *interior,
  const univ_mat  *num_coefs,
  const univ_mat  *var_block_unbiased,
  const univ_mat  *level,
  const univ_mat  *sdf,
  wav_filter_type  filter_type,
  sint32           filter_length,
  void            *intrp_ptr,
  mat_set         *result)
{
  boolean        biased = TRUE;
  boolean        normalize = TRUE;
  boolean        recentered = FALSE;
  boolean        use_edof2;
  dcomplex      *Hj;
  double         Ahat;
  double         Cj;
  double         Hsqr;
  double         num_ops = 0.0;
  double         scale;
  double         sum_Cj = 0.0;
  double         sum_Cj_squared = 0.0;
  double        *Hj2;
  double        *Hsqrj;
  double        *Sx;
  double        *edof1;
  double        *edof2;
  double        *edof3;
  double        *pvar;
  double        *s;
  mutil_errcode  err;
  sint32         Hjrow;
  sint32         extract_start = 0;
  sint32         i;
  sint32         j;
  sint32         max_level = 0;
  sint32         num_fourier;
  sint32         num_gain_frequencies = 0;
  sint32         n_level;
  sint32         num_sdf_frequencies = 0;
  sint32         num_sqrgain_frequencies;
  sint32        *N;
  sint32        *levelj;
  univ_mat       Fourier;
  univ_mat       W;
  univ_mat       acvs;
  univ_mat       gain_frequency;
  univ_mat       gain_wavelet;
  univ_mat       gain_scaling;
  univ_mat       sdf_freq;
  univ_mat       sdf_interp;
  univ_mat       wavelet_sqrgain;
  univ_mat       wavelet_sqrgain_frequency;
  univ_mat       wavelet_sqrgain_j;
  univ_mat       wavelet_sqrgain_j_interp;
  memlist        list;
  sint32         ndim;
  sint32         dims;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_variance_edof()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

   /* validate inputs */

  err = localfn_wavuniv_variance_edof_input_check(
    interior, num_coefs, var_block_unbiased, level, sdf );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* obtain sizes */

  n_level = MATUNIV_NELEM( level );

  /* allocate memory for output */

  ndim = 1;
  dims = 3;

  err = matset_malloc_register( result, ndim, &dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices( result, 1, n_level, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* define size variables */

  for ( j = 0; j < n_level; j++ ){

    max_level = MUTIL_MAX( max_level,
      level->mat.s32mat.data[ j ] );

    num_gain_frequencies = MUTIL_MAX( num_gain_frequencies,
      num_coefs->mat.s32mat.data[ j ] );
  }

  num_sqrgain_frequencies =
    (sint32) ceil( ( double ) num_gain_frequencies / 2.0 ) + 1;

  /* decide whether EDOF2 is usable */

  if ( sdf != (univ_mat *) NULL ){

    num_sdf_frequencies = MATUNIV_NELEM( sdf );

    if ( ( num_sdf_frequencies == 1 ) & ( sdf->mat.dblmat.data[ 0 ] < 0.0 ) ){
      use_edof2 = FALSE;
    }
    else{
      use_edof2 = TRUE;
    }
  }
  else{
    use_edof2 = FALSE;
  }

  /* set pointers */

  edof1  = result->mats[ 0 ].mat.dblmat.data;
  edof2  = result->mats[ 1 ].mat.dblmat.data;
  edof3  = result->mats[ 2 ].mat.dblmat.data;
  N      = &( num_coefs->mat.s32mat.data[ 0 ] );
  levelj = &( level->mat.s32mat.data[ 0 ] );

  /* allocate space */

  if ( use_edof2 ){

    err = matuniv_malloc_register( &sdf_freq, num_sdf_frequencies, 1,
      MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matuniv_malloc_register( &wavelet_sqrgain, n_level,
      num_sqrgain_frequencies, MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matuniv_malloc_register( &wavelet_sqrgain_frequency,
      num_sqrgain_frequencies, 1, MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matuniv_malloc_register( &wavelet_sqrgain_j, num_sqrgain_frequencies,
      1, MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* fill in SDF frequency vector. we impose a structure on the
       incoming sdf such that the values supplied correspond to the
       frequencies f in [0, 1/P , 2/P, 3/P, ..., (M-1)/P] where P = 2*(M-1)
       where M is the number of points in the supplied sdf input vector */

    for ( i = 0; i < num_sdf_frequencies; i++ ){

      sdf_freq.mat.dblmat.data[ i ] =
        ( double ) i  / ( double ) ( 2.0 * ( num_sdf_frequencies - 1 ) );
    }

    /* create the squared gain function matrix */

    err = wavuniv_filters_daubechies_gain(
      filter_type,
      filter_length,
      max_level,
      num_gain_frequencies,
      normalize,
      intrp_ptr,
      &gain_frequency,
      &gain_wavelet,
      &gain_scaling);
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* fill in the gain frequency vector and squared gain matrix */

    for ( j = 0; j < n_level; j++ ){

      Hjrow = levelj[ j ] - 1;

      Hj   = &( gain_wavelet.mat.cpxmat.data[ Hjrow * num_gain_frequencies ] );
      Hj2  =
        &( wavelet_sqrgain.mat.dblmat.data[ j * num_sqrgain_frequencies ] );

      for ( i = 0; i < num_sqrgain_frequencies; i++ ){

        if ( j == 0 ){
          wavelet_sqrgain_frequency.mat.dblmat.data[ i ] =
            ( double ) i / ( double ) num_gain_frequencies;
        }

        Hj2[ i ] = ( double ) MUTIL_POW( MUTIL_CPX_ABS( Hj[ i ] ), 2.0 );

      }
    }

    /* free the gain matrix */

    MUTIL_FREE_WARN( matuniv, &gain_frequency );
    MUTIL_FREE_WARN( matuniv, &gain_wavelet );
    MUTIL_FREE_WARN( matuniv, &gain_scaling );
 }


  /* begin loop over each level, obtaining the EDOFs for
     each along the way */

  for ( j = 0; j < n_level; j++ ){

    num_fourier = (sint32) floor( ( double ) ( N[ j ] - 1 ) / 2.0 );

    if (j == 0){

      /* allocate space for temporary matrices */

      err = matuniv_malloc_register( &W, N[ j ], 1, MUTIL_DOUBLE, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      if ( use_edof2 ){

        err = matuniv_malloc_register( &Fourier, num_fourier, 1,
          MUTIL_DOUBLE, &list );
        MEMLIST_FREE_ON_ERROR( err, &list );

        err = matuniv_malloc_register( &sdf_interp, num_fourier, 1,
          MUTIL_DOUBLE, &list );
        MEMLIST_FREE_ON_ERROR( err, &list );

        err = matuniv_malloc_register( &wavelet_sqrgain_j_interp,
          num_fourier, 1, MUTIL_DOUBLE, &list);
        MEMLIST_FREE_ON_ERROR( err, &list );
      }
    }
    else{

      /* reallocate space for temporary matrices */

      err = matuniv_realloc_register( &W, N[ j ], 1, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      if ( memlist_member_exist( &acvs, &list ) ){

	err = memlist_member_free( &acvs, &list );
	MEMLIST_FREE_ON_ERROR( err, &list );
      }

      if ( use_edof2 ){

        err = matuniv_realloc_register( &Fourier, num_fourier,
          1, &list );
        MEMLIST_FREE_ON_ERROR( err, &list );

        err = matuniv_realloc_register( &sdf_interp, num_fourier,
          1, &list );
        MEMLIST_FREE_ON_ERROR( err, &list );

        err = matuniv_realloc_register( &wavelet_sqrgain_j_interp,
          num_fourier, 1, &list );
        MEMLIST_FREE_ON_ERROR( err, &list );
      }
    }

    /* set pointers */

    pvar = var_block_unbiased->mat.dblmat.data;

    /* copy interior coefficients at level j into temporary vector */

    for ( i = 0; i < N[ j ]; i++ ){

      W.mat.dblmat.data[ i ] =
        interior->mat.dblmat.data[ i + extract_start ];
    }

    /* increment extraction starting point */

    extract_start += N[ j ];

    /* develop various statistics used to calculate the EDOFs ... */

    /* ... the biased, not recentered, acvs of the current
       level's interior wavelet coefficients     */

    err = mthuniv_acvs( &W, biased, recentered, intrp_ptr, &acvs );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* register the acvs with the memory manager */

    err = memlist_member_register( &list, &acvs, MEMTYPE_MATUNIV );
    MEMLIST_FREE_ON_ERROR( err, &list );

    if ( use_edof2 ){

      /* set pointers */

      Hsqrj = &( wavelet_sqrgain_j.mat.dblmat.data[ 0 ]);
      Sx    = &( sdf_interp.mat.dblmat.data[ 0 ] );

      /* ... the Fourier frequency vector for the SDF interpolation,
             i.e. the frequencies over which the appximatin is made */

      for ( i = 0; i < num_fourier; i++ ){
        Fourier.mat.dblmat.data[ i ] = ( ( double ) ( i + 1 ) ) / N[ j ];
      }

      /* ... linearly interpolate the squared gain function */

      Hj2  = &( wavelet_sqrgain.mat.dblmat.data[ j * num_sqrgain_frequencies ] );

      for ( i = 0; i < num_sqrgain_frequencies; i++ ){

        Hsqrj[ i ] = Hj2[ i ];
      }

      err = wavuniv_statistic_interpolation_linear(
        &wavelet_sqrgain_j,
        &wavelet_sqrgain_frequency,
        &Fourier,
        intrp_ptr,
        &wavelet_sqrgain_j_interp);
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* ... linearly interpolated lookup of supplied SDF over
      freqeuncies f in [1/N, 2/N, ..., M/N] where M is the number
      of Fourier frequencies and N is the number of wavelet coefficients
      at the current level */

      err = wavuniv_statistic_interpolation_linear(
        sdf,
        &sdf_freq,
        &Fourier,
        intrp_ptr,
        &sdf_interp);
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* ... band limited variance estimates based on the product of the
             supplied sdf and the wavelet gain function. here we have to
             be careful to skip the proper number of coefficients to
             obtain the correct Fourier frequency gain value. */

      sum_Cj = 0.0;
      sum_Cj_squared = 0.0;

      for ( i = 0; i < num_fourier; i++ ){

        Hsqr = wavelet_sqrgain_j_interp.mat.dblmat.data[ i ];
        Cj   =  Hsqr * Sx[ i ];

        sum_Cj += Cj;
        sum_Cj_squared += ( Cj * Cj );
      }

    }

    /* ... Ahat = ( s[0]^2 ) / 2 + sum ( s[k]^2 ) for k = 1, ..., N[j]
       where N[j] is the number of wavelet coefficients in the current level
       and s is the acvs of the interior wavelet coefficients in the current level */

    s = acvs.mat.dblmat.data;

    Ahat = s[ 0 ] * s[ 0 ] / 2.0;

    for ( i = 1; i < MATUNIV_NELEM( &acvs ); i++ ){
      Ahat += s[ i ] * s[ i ];
    }

    /* ... the scale at the current level */

    scale = MUTIL_POW( 2.0, ( double ) ( level->mat.s32mat.data[ j ] - 1 ) );

    /* caluclate the EDOFs ... */

    edof1[ j ] = ( double ) N[ j ] * pvar[ j ] * pvar[ j ] / Ahat;

    if ( use_edof2 ){
      edof2[ j ] = 2.0 * sum_Cj * sum_Cj / sum_Cj_squared;
    }
    else{
      edof2[ j ] = -1.0;
    }

    edof3[ j ] = MUTIL_MAX( N[ j ] / ( scale * 2.0 ), 1.0 );

    /* check for interrupts */

    num_ops += 10.0 * ( double ) N[ j ];
    if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }

  }  /* end loop over each level */

  /* free node corresponding to registered
     ouput memory, but do not free the
     memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_variance_edof()" );

  return MUTIL_ERR_OK;
}

/* Confidence limits for the unbiased blocked */
/* discrete wavelet variance estimates        */
/* Documented in wav_var.h                    */
/* Written by William Constantine             */

mutil_errcode wavuniv_variance_confidence(
  const univ_mat *variance,
  const univ_mat *edof,
  double          probability,
  void           *intrp_ptr,
  mat_set        *result )
{
  double            qchisq_low;
  double            qchisq_high;
  double           *pedof;
  double           *pvar;
  double           *low;
  double           *high;
  double            num_ops = 0.0;
  /* double            temp; */
  mutil_errcode     err;
  sint32            j;
  sint32            n_dim = 1;
  sint32            dims = 2;
  sint32            n_level;
  memlist           list;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_variance_confidence()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* validate input arguments */

  err = localfn_wavuniv_variance_confidence_input_check(
   variance, edof, probability );
  if ( err ) return err;

  /* obtain matrix sizes */

  n_level = MATUNIV_NELEM( edof );

  /* allocate memory for output and
     register with memory manager */

  err = matset_malloc_register( result, n_dim, &dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices( result, 1, n_level, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* set pointers */

  low   = result->mats[ 0 ].mat.dblmat.data;
  high  = result->mats[ 1 ].mat.dblmat.data;
  pedof = &( edof->mat.dblmat.data[ 0 ] );
  pvar  = &( variance->mat.dblmat.data[ 0 ] );

  /* loop through each level and calculate
     the low and high confidence limits */

  for ( j = 0; j < n_level; j++ ){

    /* calculate the chi-square quantiles */

    qchisq_low  = mth_qchisq( probability,  pedof[ j ] );
    qchisq_high = mth_qchisq( 1.0 - probability,  pedof[ j ] );

    low[ j ]  = pedof[ j ] * pvar[ j ] / qchisq_low;
    high[ j ] = pedof[ j ] * pvar[ j ] / qchisq_high;

    /* ensure low < high. if not, swap */

   /* if ( low[ j ] > high[ j ] ){

      temp      = low[ j ];
      low[ j ]  = high[ j ];
      high[ j ] = temp;
    }*/

    /* Check for interrupt */

    num_ops += 3.0 * n_level;
    if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }
  }

  /* free node corresponding to registered
     ouput memory, but do not free the
     memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_variance_confidence()" );

  return MUTIL_ERR_OK;
}

/***************************************/
/* STATIC FUNCTION DEFINITIONS         */
/***************************************/

/** Validate the inputs to the wavelet variance confidence function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_var.c
 * @usage #err = localfn_wavuniv_variance_confidence_input_check( &variance, &edof, probability );#
 * @return              Standard mutils error/OK code.
 * @param  variance     Pointer to pre-allocated universal matrix of type
 *                      MUTIL\_DOUBLE containing the blocked biased
 *                      discrete wavelet variance estimates (one for each
 *                      decomposition level). This matrix must be a
 *                      row or column vector with J elements where J is the
 *                      number of decomposition levels.
 * @param  edof         Pointer to pre-allocated universal matrix of type
 *                      MUTIL\_DOUBLE containing the equivalent
 *                      degrees of freedom estimates (one for each
 *                      decomposition level). This matrix must be a
 *                      a row or column vector with J elements where J is
 *                      the number of decomposition levels.
 * @param  probability  The probability desired for the confidence intervals.
 *                      Allowable values are 0.005, 0.025, 0.05, 0.95, 0.975,
 *                      and 0.995.
 * @see wavuniv_variance_confidence
 */
static mutil_errcode localfn_wavuniv_variance_confidence_input_check(
  const univ_mat *variance,
  const univ_mat *edof,
  double          probability)
{
  boolean       found = FALSE;
  double        eps = 0.000001;
  double        p[] = {0.005, 0.025, 0.05, 0.95, 0.975, 0.995};
  mutil_errcode err;
  sint32        i;

  MUTIL_TRACE( "Start localfn_wavuniv_variance_confidence_input_check()" );

  /* check matrix inputs */

  LOCALDEF_CHECK_UNIVMAT_VECTOR( variance, MUTIL_DOUBLE, dblmat );
  LOCALDEF_CHECK_UNIVMAT_VECTOR( edof, MUTIL_DOUBLE, dblmat );

  /* match desired probability with acceptable probabilities */

  for ( i = 0; i < 6; i++ ){

    if ( MUTIL_ABS( p[ i ] - probability ) < eps ){
      found = TRUE;
      break;
    }
  }

  if ( !found ){
     MUTIL_ERROR( "Unsupported probability" );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }

  MUTIL_TRACE( "Done with localfn_wavuniv_variance_confidence_input_check()" );

  return MUTIL_ERR_OK;
}

/** Validate the inputs to the wavelet variance confidence function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_var.c
 * @usage #err = localfn_wavuniv_variance_edof_input_check( &interior, &num_coefs, &variance, &level, &sdf );#
 * @return            Standard mutils error/OK code.
 * @param interior    Pointer to a single-column or single-row
 *                    universal matrix of type MUTIL\_DOUBLE
 *                    containing a concatenated collection of
 *                    interior wavelet coefficients in order
 *                    of wavelet transform decomposition level.
 * @param num_coefs   Pointer to a single-column or single-row
 *                    universal matrix of type MUTIL\_SINT32
 *                    containing the number of wavelet transform
 *                    coefficents in each level.
 * @param variance    Pointer to a single-column or single-row
 *                    universal matrix of type MUTIL\_DOUBLE
 *                    containing the biased wavelet variance estimates.
 * @param level       Pointer to a single-column or single-row
 *                    universal matrix of type MUTIL\_SINT32
 *                    containing the decomposition levels.
 *                    This argument coordinates
 *                    with the interior, num\_coefs, and variance
 *                    arguments. It therefore does not need to contain
 *                    a monotonic continuous progression of levels.
 * @param sdf         Pointer to a single-column or single-row
 *                    universal matrix of type MUTIL\_DOUBLE
 *                    containing a discretized approximation of the process
 *                    spectral density function. The coefficients
 *                    of this argument should correspond exactly
 *                    with the normalized Fourier frequencies
 *                    f = [0, 1/P , 2/P, 3/P, ..., (M-1)/P] where
 *                    P = 2*(M-1) and M is the number of points
 *                    in the sdf vector. For example, if the sdf
 *                    vector contains 5 elements,
 *                    the corresponding frequencies will be
 *                    f = [0, 1/8, 1/4, 3/8, 1/2]. This input
 *                    can also be a NULL pointer, and (in that
 *                    case signifies that an EDOF mode 2 is
 *                    not to be calculated.
 *
 * @see wavuniv_variance_edof
 */
static mutil_errcode localfn_wavuniv_variance_edof_input_check(
  const univ_mat *interior,
  const univ_mat *n_coeff,
  const univ_mat *variance,
  const univ_mat *level,
  const univ_mat *sdf)
{
  mutil_errcode err;

  MUTIL_TRACE( "Start localfn_wavuniv_variance_edof_input_check()" );

  /* check matrix inputs */

  LOCALDEF_CHECK_UNIVMAT_VECTOR( interior, MUTIL_DOUBLE, dblmat );
  LOCALDEF_CHECK_UNIVMAT_VECTOR( n_coeff, MUTIL_SINT32, s32mat );
  LOCALDEF_CHECK_UNIVMAT_VECTOR( variance, MUTIL_DOUBLE, dblmat );
  LOCALDEF_CHECK_UNIVMAT_VECTOR( level, MUTIL_SINT32, s32mat );

  /*** check sdf argument */

  err = localfn_check_sdf( sdf );
  if ( err ) return err;

  MUTIL_TRACE( "Done with localfn_wavuniv_variance_edof_input_check()" );

  return MUTIL_ERR_OK;
}

/** Validate the inputs to the wavelet variance function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_var.c
 * @usage #err = localfn_wavuniv_variance_input_check( &time_series, transform_type, filter_type, filter_length, n_level, &sdf );#
 * @return                  Standard mutils error/OK code.
 * @param time_series       Pointer to a pre-allocated universal matrix
 *                          containing the time series to analyze.
 *                          This input must a single-row or single-column
 *                          matrix and can be of any type with the
 *                          exception of MUTIL\_DCOMPLEX.
 * @param transform_type    The type of wavelet transform to perform.
 *                          Valid choices are WAV\_TRANSFORM\_DWT and
 *                          WAV\_TRANSFORM\_MODWT.
 * @param filter_type       The type of Daubechies filter.
 * @param filter_length     The length of the wavelet filter.
 * @param n_level           The number of wavelet transform decomposition
 *                          levels.
 * @param sdf               Pointer to a single-column or single-row
 *                          universal matrix of type MUTIL\_DOUBLE
 *                          containing a discretized approximation of
 *                          the process spectral density function.
 *                          The coefficients of this argument should
 *                          correspond exactly with the normalized
 *                          Fourier frequencies
 *                          f = [0, 1/P , 2/P, 3/P, ..., (M-1)/P] where
 *                          P = 2*(M-1) and M is the number of points
 *                          in the sdf vector. For example, if the sdf
 *                          vector contains 5 elements,
 *                          the corresponding frequencies will be
 *                          f = [0, 1/8, 1/4, 3/8, 1/2]. This input
 *                          can also be a NULL pointer, and (in that
 *                          case signifies that an EDOF mode 2 is
 *                          not to be calculated.
 * @see wavuniv_variance
 */
static mutil_errcode localfn_wavuniv_variance_input_check(
 const univ_mat  *time_series,
 wav_transform    transform_type,
 wav_filter_type  filter_type,
 sint32           filter_length,
 sint32           n_level,
 const univ_mat  *sdf )
{
  mutil_errcode err;

  MUTIL_TRACE( "Start localfn_wavuniv_variance_input_check()" );

  /*** check time series ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_VAR( time_series, univ_mat, matuniv );

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( &(time_series->mat.dblmat) ) ){
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

  /*** check the transform type ... ***/

  if ( ( transform_type != (wav_transform) WAV_TRANSFORM_DWT ) &&
       ( transform_type != (wav_transform) WAV_TRANSFORM_MODWT ) ){
    MUTIL_ERROR( "Transofrm type must be either "
		 "WAV_TRANSFORM_DWT or WAV_TRANSFORM_MODWT." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  /*** filter_length argument ... */

  /* ... if filter_length is positive */

  if ( filter_length <= 1 ){
    MUTIL_ERROR( "Filter length must be greater than one." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if filter_length is even */

  if ( ( filter_length % 2 ) != 0 ){
    MUTIL_ERROR( "Filter length must be even." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*** check filter_type argument */

  switch( filter_type ){
  case WAV_FILTER_EXTREMAL_PHASE:
  case WAV_FILTER_LEAST_ASYMMETRIC:
  case WAV_FILTER_BEST_LOCALIZED:
  case WAV_FILTER_COIFLET:
  case WAV_FILTER_HAAR:
    break;
  default:
    MUTIL_ERROR( "Filter type is unsupported" );
    return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /*** check n_level argument */

  if ( n_level <= 0 ){
    MUTIL_ERROR( "Number of decomposition levels must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  /*** check sdf argument */

  err = localfn_check_sdf( sdf );
  if ( err ) return err;

  MUTIL_TRACE( "Done with localfn_wavuniv_variance_input_check()" );

  return MUTIL_ERR_OK;
}

/** Validate the SDF input for wavelet variance functions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_var.c
 * @usage #err = localfn_check_sdf( &sdf );#
 * @return                  Standard mutils error/OK code.
 * @param sdf               Pointer to a single-column or single-row
 *                          universal matrix of type MUTIL\_DOUBLE
 *                          containing a discretized approximation of
 *                          the process spectral density function.
 *                          The coefficients of this argument should
 *                          correspond exactly with the normalized
 *                          Fourier frequencies
 *                          f = [0, 1/P , 2/P, 3/P, ..., (M-1)/P] where
 *                          P = 2*(M-1) and M is the number of points
 *                          in the sdf vector. For example, if the sdf
 *                          vector contains 5 elements,
 *                          the corresponding frequencies will be
 *                          f = [0, 1/8, 1/4, 3/8, 1/2]. This input
 *                          can also be a NULL pointer, and (in that
 *                          case signifies that an EDOF mode 2 is
 *                          not to be calculated.
 *
 * @see localfn_wavuniv_variance_input_check
 * @see localfn_wavuniv_variance_edof_input_check
 */
static mutil_errcode localfn_check_sdf(
  const univ_mat *sdf )
{
  mutil_errcode err;

  /*** check SDF argument ...

       NOTE: THIS ARGUMENT CAN BE A NULL POINTER.
       IF IT IS, IT SIGNIFIES THAT THE USER DOES
       NOT WISH TO CALCULATE AN EDOF MODE 2   */

  if ( sdf != ( univ_mat * ) NULL ){

    /* ... for valid matrix structure */

    err = matuniv_validate( sdf );
    if ( err ) return err;

    /* ... for the proper type */

    if ( sdf->type != MUTIL_DOUBLE ){
      MUTIL_ERROR( "SDF matrix must be of type MUTIL_DOUBLE." );
      return MUTIL_ERR_ILLEGAL_TYPE;
    }

    /* ... to see if it is a vector */

    if ( !MATANY_IS_VEC( &(sdf->mat.dblmat) ) ){
      MUTIL_ERROR( "SDF matrix must be a single column or row." );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* ... if number of elements is positive */

    if ( MATUNIV_NELEM( sdf ) < 1 ){
      MUTIL_ERROR( "Number of elements in SDF must be positive." );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
  }

  return MUTIL_ERR_OK;
}
