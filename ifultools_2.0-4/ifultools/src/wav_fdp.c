
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/wav_fdp.c $: $Revision: #1 $, $Date: 2008/03/21 $";

#include "wav_dwtc.h"
#include "wav_fdp.h"
#include "wav_filt.h"
#include "wav_look.h"
#include "wav_modw.h"
#include "wav_type.h"
#include "wav_var.h"

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_limit.h"
#include "ut_math.h"
#include "ut_mem.h"
#include "ut_rand.h"

#include "mat_any.h"
#include "mat_arit.h"
#include "mat_assn.h"
#include "mat_set.h"
#include "mat_sort.h"
#include "mat_summ.h"
#include "mat_type.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "mth_mac.h"

/***********************************************/
/* STATIC FDP ESTIMATION FUNCTION DECLARATIONS */
/***********************************************/

static mutil_errcode localfn_parabola_minimum_abscissa(
  const double_mat  *x,
  const double_mat  *y,
  double            *result);

static mutil_errcode localfn_Cprime_lookup_table(
  const univ_mat    *levels,
  const univ_mat    *delta_range,
  void              *intrp_ptr,
  double_mat        *result);

static mutil_errcode localfn_minimum_bracket(
  const double_mat *Cprime,
  wave_energy      *wave,
  const double_mat *sum_NjlogCprime,
  sint32            delta_start,
  double_mat       *xmin,
  double_mat       *ymin);

static mutil_errcode localfn_wavelet_energy_instantaneous(
  const mat_set     *modwt,
  const mat_set     *shift,
  const univ_mat    *levels,
  const univ_mat    *interior,
  sint32             n_level_usable,
  boolean            biased,
  sint32             time,
  sint32             dof_order,
  wav_fdp_estimator  estimator,
  sint32_mat        *n_coeff,
  double_mat        *Wtmodsqr );

static mutil_errcode localfn_wavuniv_transform_energy_block(
  const univ_mat    *time_series,
  const univ_mat    *levels,
  wav_filter_type    filter_type,
  sint32             filter_length,
  boolean            biased,
  wav_fdp_estimator  estimator,
  sint32             edof_mode,
  const univ_mat    *sdf,
  void              *intrp_ptr,
  double_mat        *edof,
  sint32_mat        *n_coeff,
  double_mat        *energy );

static mutil_errcode localfn_innovation_variance(
  const double_mat *Cprime,
  wave_energy      *wave,
  sint32            row_delta,
  sint32           *n_coeff_used,
  double           *innovation_variance);

static mutil_errcode localfn_delta_indices(
  const univ_mat    *delta_range,
  sint32            *start,
  sint32            *end,
  sint32            *n_delta);

static mutil_errcode localfn_reduced_log_likelihood(
  const double_mat *Cprime,
  wave_energy      *wave,
  const double_mat *sum_NjlogCprime,
  sint32            delta_index,
  double           *result);

static mutil_errcode localfn_interior_coefficient_index_range(
  const mat_set     *shift,
  const univ_mat    *levels,
  sint32             n_level_max,
  sint32             filter_length,
  sint32             n_sample,
  void              *intrp_ptr,
  univ_mat          *result);

static mutil_errcode localfn_number_usable_levels(
  const univ_mat    *interior,
  boolean            biased,
  sint32             time,
  sint32            *result);

static mutil_errcode localfn_fdp_estimator_input_check(
  const univ_mat    *time_series,
  const univ_mat    *levels,
  wav_filter_type    filter_type,
  sint32             filter_length,
  boolean            biased,
  wav_fdp_estimator  estimator,
  sint32             dof,
  const univ_mat    *delta_range,
  void              *intrp_ptr,
  boolean            block );

/***********************************************/
/* STATIC FDP SDF FUNCTION DECLARATIONS        */
/***********************************************/

static mutil_errcode localfn_taylorseries_u(
  const univ_mat *x,
  double          beta,
  univ_mat       *result );

static mutil_errcode localfn_taylorseries_du(
  const univ_mat *x,
  double          beta,
  univ_mat       *result );

static mutil_errcode localfn_taylorseries_ddu(
  const univ_mat *x,
  double          beta,
  univ_mat       *result );

static mutil_errcode localfn_taylorseries_coeff_weights(
  const univ_mat *levels,
  sint32          n,
  double          beta,
  univ_mat       *result );

static mutil_errcode localfn_taylorseries_coeff(
  const univ_mat *levels,
  double          delta,
  univ_mat       *result );

static mutil_errcode localfn_fdp_bandpass_levels12(
  double    delta,
  double    innovation_variance,
  void     *intrp_ptr,
  univ_mat *result );

static mutil_errcode localfn_fdp_bandpass_levels3J(
  const univ_mat *levels,
  double          delta,
  double          innovation_variance,
  univ_mat       *result );

static mutil_errcode localfn_fdp_bandpass_scaling(
  univ_mat   *Cprime,
  double      delta,
  double      innovation_variance,
  sint32      n_sample);

static mutil_errcode localfn_fdp_bandpass_arg_check(
  const univ_mat *levels,
  sint32          n_sample,
  univ_mat       *result );

/***********************************************/
/* STATIC FDP ESTIMATION MACRO DECLARATIONS    */
/***********************************************/

#define SHIFT2(a,b,c) (a)=(b);(b)=(c);
#define SHIFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define MINIMUM_USABLE_LEVELS 2
#define COERCE(x) ( ( (x) >= 0 ) ? ceil(x) : floor(x) )

#define LOCALDEF_CHECK_NULL_POINTER_FDP( DATA_PTR, DATA_TYPE,  \
                                     TYPE_PREFIX )             \
   err = TYPE_PREFIX ## _validate( DATA_PTR );                 \
   if ( err ) return err;                                      \
   if ( DATA_PTR == ( DATA_TYPE * ) NULL ){              \
     MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
     return MUTIL_ERR_NULL_POINTER;                            \
   }

/***********************************************/
/* STATIC FDP SDF MACRO DEFINITIONS            */
/***********************************************/

/* NONE */


/***********************************************/
/* FDP ESTIMATION FUNCTION DEFINITIONS         */
/***********************************************/

/* Function to estimate instantaneous */
/* FD model parameters                */
/* Function documented in wav_fdp.h   */
/* Written by William Constantine     */

mutil_errcode wavuniv_fdp_estimator_instantaneous(
  const univ_mat    *time_series,
  const univ_mat    *levels,
  wav_filter_type    filter_type,
  sint32             filter_length,
  wav_fdp_estimator  estimator,
  boolean            biased,
  sint32             dof_order,
  const univ_mat    *delta_range,
  void              *intrp_ptr,
  univ_mat          *delta,
  univ_mat          *variance_delta,
  univ_mat          *innovation_variance)
{
  boolean       normalize_filters = TRUE;
  double        Yt;
  double        den;
  double        error_constant = 0.0;
  double        half_dof;
  double        num;
  double        sqrsumlogtau;
  double        sum;
  double        sumYt;
  double        sumlogtau;
  double        sumlogtauYt;
  double        sumlogtausqr;
  double        var_weight = 0.0;
  double_mat    Cprime;
  double_mat    Wtmodsqr;
  double_mat    logtau;
  double_mat    logtausqr;
  double_mat    mle_delta_bracket;
  double_mat    mle_rll_bracket;
  double_mat    sum_NjlogCprime;
  mat_set       filters;
  mat_set       modwt;
  memlist       list;
  mutil_errcode err;
  sint32        col;
  sint32        delta_end;
  sint32        delta_start = 0;
  sint32        j;
  sint32        level;
  sint32        n_delta = 0;
  sint32        n_level_max = 0;
  sint32        n_level_min;
  sint32        n_level;
  sint32        n_level_usable = 0;
  sint32        n_sample;
  sint32        row;
  sint32        row_global;
  sint32        time;
  sint32        n_coeff_used_innovation_variance;
  sint32_mat    n_coeff;
  univ_mat      Cprime_mle;
  univ_mat      interior;
  univ_mat      levels_sorted;
  mat_set       shift;
  wave_energy   wave;

  /* avoid lint warning */

  (void) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  MUTIL_TRACE( "Start wavuniv_fdp_estimator_instantaneous() ::" );

  /* check inputs for errors */

  err = localfn_fdp_estimator_input_check(
    time_series,
    levels,
    filter_type,
    filter_length,
    biased,
    estimator,
    dof_order,
    delta_range,
    intrp_ptr,
    FALSE );
  if ( err ) return err;

  /* obtain or specify matrix sizes */

  n_sample = MATUNIV_NELEM( time_series );
  n_level = MATUNIV_NELEM( levels );

  /* allocate space for output */

  err = matuniv_malloc_register( delta, n_sample, 1,
    MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( variance_delta, n_sample, 1,
    MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( innovation_variance, n_sample, 1,
    MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* obtain the range of desired decomposition levels */

  err = mats32_range( &(levels->mat.s32mat), intrp_ptr,
                      &n_level_min, &n_level_max );
  if ( err ) return err;

  /* obtain the number of delta values to
     access in the Cprime table based on the
     range of delta supplied by the user */

  err = localfn_delta_indices(
    delta_range, &delta_start, &delta_end, &n_delta );
  if ( err ) return err;

  /* malloc space for common objects */

  err = matdbl_malloc_register( &Wtmodsqr, n_level, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &n_coeff, n_level, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &levels_sorted, n_level, 1,
                             MUTIL_SINT32, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &interior, 2, n_level, MUTIL_SINT32, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( estimator == WAV_FDP_MAXIMUM_LIKELIHOOD ){

    /* malloc space for MLE objects */

    err = matdbl_malloc_register( &mle_rll_bracket, 3, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( &mle_delta_bracket, 3, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( &Cprime, n_delta, n_level, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( &sum_NjlogCprime, n_delta, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }
  else if ( estimator == WAV_FDP_LEAST_SQUARES ){

    /* malloc space for least squares variables */

    err = matdbl_malloc_register( &logtau, n_level, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( &logtausqr, n_level, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* obtain the wavelet and scaling filters
     and register them with the memory manager */

  err = wavuniv_filters_daubechies(
    filter_length,
    filter_type,
    normalize_filters,
    intrp_ptr,
    &filters );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &filters, MEMTYPE_MATSET );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* sort the levels vector to make things tidier */

  err = matuniv_sort(
    levels,
    intrp_ptr,
    &levels_sorted );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate the MODWT and register
     the memory with the memory manager */

  err = wavuniv_transform_maximum_overlap(
      time_series,
      &filters,
      n_level_max,
      intrp_ptr,
      &modwt);
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &modwt, MEMTYPE_MATSET );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* obtain MODWT and DWT zero phase shift factors
     and register with the memory manager */

  err = wavuniv_filters_zero_phase(
      filter_type,
      filter_length,
      n_level_max,
      intrp_ptr,
      &shift );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &shift, MEMTYPE_MATSET );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* obtain the range of indices for the (zero phase)
     shifted interior wavelet coefficients */

  err = localfn_interior_coefficient_index_range(
    &shift,
    &levels_sorted,
    n_level_max,
    filter_length,
    n_sample,
    intrp_ptr,
    &interior);
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( estimator == WAV_FDP_MAXIMUM_LIKELIHOOD ){

    /* obtain the mid-octave SDF values for the specified
       levels and range of delta via a lookup table */

    err = localfn_Cprime_lookup_table(
      &levels_sorted,
      delta_range,
      intrp_ptr,
      &Cprime);
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* form the sum of the log of the Cprime values
       for each delta (row) in the Cprime matrix */

    for ( row = 0; row < n_delta; row++ ){

      sum = 0.0;

      row_global = row * n_level;

      for ( col = 0; col < n_level; col++ ){
        sum += log( Cprime.data[ row_global + col ] );
      }

      sum_NjlogCprime.data[ row ] = sum;
    }
  }
  else if ( estimator == WAV_FDP_LEAST_SQUARES ){

    /* set factors relating to the chi-squared degrees of freedom */

    half_dof       = (double) dof_order + 0.5;
    var_weight     = MUTIL_TRIGAMMA( half_dof );
    error_constant = - MUTIL_DIGAMMA( half_dof ) + log( half_dof );

    /* set level dependent least squares variables */

    for ( j = 0; j < n_level; j++ ){

      level               = levels_sorted.mat.s32mat.data[ j ];
      logtau.data[ j ]    = log( MUTIL_POW( 2.0, (double) ( level - 1 ) ) );
      logtausqr.data[ j ] = logtau.data[ j ] * logtau.data[ j ];
    }
  }

  /* set level independent information for
     the wavelet energy structure */

  wave.level      = &levels_sorted.mat.s32mat;
  wave.n_level    = n_level;
  wave.n_sample   = n_sample;
  wave.biased     = biased;
  wave.stationary = (boolean) FALSE;

  /* loop through each point in time and calculate the
     MLE or LSE of delta, the corresponding innovation
     variance, and the variance of the delta estimate */

  for ( time = 0; time < n_sample; time++ ){

    /* make sure that at the minimum
       number of usable levels is achieved */

    err = localfn_number_usable_levels(
      &interior,
      biased,
      time,
      &n_level_usable);
    MEMLIST_FREE_ON_ERROR( err, &list );

    if ( n_level_usable >= MINIMUM_USABLE_LEVELS ){

      /* calculate the square of the
         (zero phase) shifted MODWT
         coefficients for each specified
         level at the current time */

      err = localfn_wavelet_energy_instantaneous(
        &modwt,
        &shift,
        &levels_sorted,
        &interior,
        n_level_usable,
        biased,
        time,
        dof_order,
        estimator,
        &n_coeff,
        &Wtmodsqr );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* set the values in the energy structure
      that change with time */

      wave.energy         = &Wtmodsqr;
      wave.n_coeff        = &n_coeff;
      wave.n_level_usable = n_level_usable;

      /* calculate the output */

      if ( estimator == WAV_FDP_MAXIMUM_LIKELIHOOD ){

        /* bracket the minimum reduced log-likelihood value */

        err = localfn_minimum_bracket(
          &Cprime,
          &wave,
          &sum_NjlogCprime,
          delta_start,
          &mle_delta_bracket,
          &mle_rll_bracket);
        MEMLIST_FREE_ON_ERROR( err, &list );

        /* calculate the MLE of delta ... */

        /* ... if no bracket of a minimum
           values was found, then the first
           vlaue in the rll_delta_bracket vector
           will be set to MUTIL_DOUBLE_EPSILON and the
           second value will be set to the boundary
           encountered during bracketing (this will be either
           FDP_SDF_TABLE_DELTA_MIN or FDP_SDF_TABLE_DELTA_MAX
           in the FDP_SDF_TABLE )  */

        if ( mle_delta_bracket.data[ 0 ] == MUTIL_DOUBLE_EPSILON ){
          delta->mat.dblmat.data[ time ] =  mle_delta_bracket.data[ 1 ];
        }
        else{

          err = localfn_parabola_minimum_abscissa(
            &mle_delta_bracket,
            &mle_rll_bracket,
            &(delta->mat.dblmat.data[ time ]) );
          MEMLIST_FREE_ON_ERROR( err, &list );
        }

        /* calculate corresponding MLE of the innovation variance ... */

        /* ... first calculate the mid-octave SDF value for
           an FD process with FD parameter equal to the MLE
           of delta. here we do not use the lookup table but
           call the appropriate function directly. */

        err = wavuniv_fdp_bandpass_variance(
          &levels_sorted,
          delta->mat.dblmat.data[ time ],
          -1,
          intrp_ptr,
          &Cprime_mle );
        MEMLIST_FREE_ON_ERROR( err, &list );

        err = localfn_innovation_variance(
          &(Cprime_mle.mat.dblmat),
          &wave,
          0,
          &n_coeff_used_innovation_variance,
          innovation_variance->mat.dblmat.data + time );
        MEMLIST_FREE_ON_ERROR( err, &list );

        /* free the Cprime_mle vector */

        MUTIL_FREE_WARN( matuniv, &Cprime_mle );

      } /* end of MLE FDP estimation */

        else if ( estimator == WAV_FDP_LEAST_SQUARES ){

          sumYt        = 0.0;
          sumlogtau    = 0.0;
          sumlogtauYt  = 0.0;
          sumlogtausqr = 0.0;

          for ( j = 0; j < n_level_usable; j++ ){

            Yt = log ( Wtmodsqr.data[ j ] ) + error_constant;

            sumYt        += Yt;
            sumlogtau    += logtau.data[ j ];
            sumlogtauYt  += logtau.data[ j ] * Yt;
            sumlogtausqr += logtausqr.data[ j ];
          }

          sqrsumlogtau = sumlogtau * sumlogtau;

          num = ( (double) n_level_usable * sumlogtauYt ) - sumlogtau * sumYt;
          den = ( (double) n_level_usable * sumlogtausqr ) - sqrsumlogtau;

          delta->mat.dblmat.data[ time ] = ( ( num / den ) + 1.0 ) / 2.0;

          variance_delta->mat.dblmat.data[ time ] =
            (double) n_level_usable * var_weight / den / 4.0;

      } /* end of LSE FDP estimation */

    }
    else{

      /* there are not enough usable levels to make an MLE */

      innovation_variance->mat.dblmat.data[ time ] = MUTIL_DOUBLE_MAX;
      delta->mat.dblmat.data[ time ] = MUTIL_DOUBLE_MAX;

    } /* end of minimum usable levels check */
  } /* end of loop through time */


  /* free the node corresponding to the
     outputs in the memory list but not the
     memory itself */

  err = memlist_member_unregister( delta, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_unregister( variance_delta, &list );
  if ( err ){
    MUTIL_FREE_WARN( matuniv, delta );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  err = memlist_member_unregister( innovation_variance, &list );
  if ( err ){
    MUTIL_FREE_WARN( matuniv, delta );
    MUTIL_FREE_WARN( matuniv, variance_delta );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* free the remaining registered memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( ":: Done with wavuniv_fdp_estimator_instantaneous()" );

  return MUTIL_ERR_OK;
}

/* Function to estimate block */
/* FD model parameters                */
/* Function documented in wav_fdp.h   */
/* Written by William Constantine     */

mutil_errcode wavuniv_fdp_estimator_block(
  const univ_mat    *time_series,
  const univ_mat    *levels,
  wav_filter_type    filter_type,
  sint32             filter_length,
  wav_fdp_estimator  estimator,
  boolean            biased,
  sint32             edof_mode,
  const univ_mat    *sdf,
  const univ_mat    *delta_range,
  void              *intrp_ptr,
  double            *delta,
  double            *variance_delta,
  double            *innovation_variance )
{
  boolean       normalize_filters = TRUE;
  double        Y;
  double        beta;
  double        den;
  double        half_edof;
  double        logtau;
  double        num;
  double        sum;
  double        sumlogtau;
  double        sumwj;
  double        sumwjY;
  double        sumwjlogtau;
  double        sumwjlogtauY;
  double        sumwjlogtausqr;
  double        tau;
  double        var;
  double        wj;
  double        wjlogtau;
  double_mat    Cprime;
  double_mat    edof;
  double_mat    energy;
  double_mat    mle_delta_bracket;
  double_mat    mle_rll_bracket;
  double_mat    sum_NjlogCprime;
  mat_set       filters;
  memlist       list;
  mutil_errcode err;
  sint32        col;
  sint32        delta_end;
  sint32        delta_start = 0;
  sint32        j;
  sint32        jj;
  sint32        n_coeff_used_innovation_variance;
  sint32        n_delta = 0;
  sint32        n_level;
  sint32        n_level_max = 0;
  sint32        n_level_min;
  sint32        n_sample;
  sint32        row;
  sint32        row_global;
  sint32_mat    n_coeff;
  univ_mat      Cprime_mle;
  univ_mat      levels_sorted;
  wave_energy   wave;

  /* avoid lint warning */

  (void) whatssi;

  MUTIL_TRACE( "Start wavuniv_fdp_estimator_block() ::" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  err = localfn_fdp_estimator_input_check(
    time_series,
    levels,
    filter_type,
    filter_length,
    biased,
    estimator,
    edof_mode,
    delta_range,
    intrp_ptr,
    TRUE );
  if ( err ) return err;

  /* make sure that the pointer to the sdf matrix
     is not NULL if the EDOF mode is equal to 2 */

  if ( ( edof_mode == 2 ) && ( sdf == (univ_mat *) NULL ) ){

     MUTIL_ERROR( "An SDF vector must be supplied for "
       "EDOF mode 2 calculations." );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* obtain or specify matrix sizes */

  n_level  = MATUNIV_NELEM( levels );
  n_sample = MATUNIV_NELEM( time_series );

  /* obtain the range of desired decomposition levels */

  err = mats32_range( &(levels->mat.s32mat), intrp_ptr,
                      &n_level_min, &n_level_max );
  if ( err ) return err;

   /* obtain the number of delta values to
     access in the Cprime table based on the
     range of delta supplied by the user */

  err = localfn_delta_indices(
    delta_range, &delta_start, &delta_end, &n_delta );
  if ( err ) return err;

  /* malloc space for common objects */

  err = matuniv_malloc_register( &levels_sorted, n_level, 1,
                             MUTIL_SINT32, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( estimator == WAV_FDP_MAXIMUM_LIKELIHOOD ){

    /* malloc space for MLE objects */

    err = matdbl_malloc_register( &mle_rll_bracket, 3, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( &mle_delta_bracket, 3, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( &Cprime, n_delta, n_level, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( &sum_NjlogCprime, n_delta, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* obtain the wavelet and scaling filters
     and register them with the memory manager */

  err = wavuniv_filters_daubechies(
    filter_length,
    filter_type,
    normalize_filters,
    intrp_ptr,
    &filters );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &filters, MEMTYPE_MATSET );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* sort the levels vector to make things tidier */

  err = matuniv_sort( levels, intrp_ptr, &levels_sorted );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate energy of the wavelet coefficients */

  err = localfn_wavuniv_transform_energy_block(
    time_series,
    &levels_sorted,
    filter_type,
    filter_length,
    biased,
    estimator,
    edof_mode,
    sdf,
    intrp_ptr,
    &edof,
    &n_coeff,
    &energy );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( estimator == WAV_FDP_LEAST_SQUARES ){

    err = memlist_member_register( &list, &edof, MEMTYPE_MATDBL );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  err = memlist_member_register( &list, &energy, MEMTYPE_MATDBL );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &n_coeff, MEMTYPE_MATS32 );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( estimator == WAV_FDP_LEAST_SQUARES ){

    /* initialize summation variables */

    sumwj          = 0.0;
    sumwjY         = 0.0;
    sumlogtau      = 0.0;
    sumwjlogtau    = 0.0;
    sumwjlogtauY   = 0.0;
    sumwjlogtausqr = 0.0;

    for ( j = 0; j < n_level; j++ ){

      jj = levels_sorted.mat.s32mat.data[ j ] - 1;

      var = energy.data[ j ] / n_coeff.data[ j ];
      half_edof = edof.data[ j ] / 2.0;

      Y        = log( var ) - MUTIL_DIGAMMA( half_edof ) + log( half_edof );
      wj       = 1.0 / MUTIL_TRIGAMMA( half_edof );
      tau      = (double) MUTIL_POW( 2.0, (double) jj );
      logtau   = log( tau );
      wjlogtau = wj * logtau;

      sumwj          += wj;
      sumlogtau      += logtau;
      sumwjlogtau    += wjlogtau;
      sumwjlogtauY   += wjlogtau * Y;
      sumwjlogtausqr += wjlogtau * logtau;
      sumwjY         += wj * Y;
    }

    num  = sumwj * sumwjlogtauY - sumwjlogtau * sumwjY;
    den  = sumwj * sumwjlogtausqr - sumwjlogtau * sumwjlogtau;
    beta = num / den;

    *delta          = ( beta + 1.0 ) / 2.0;
    *variance_delta = sumwj / den / 4.0;

    /* we do not have a method to calculate the innovations variance for
       WLSE, so just set it to some unrealistic number here */

    *innovation_variance = -1.0;
  }
  else if ( estimator == WAV_FDP_MAXIMUM_LIKELIHOOD ){

    /* obtain the mid-octave SDF values for the specified
       levels and range of delta via a lookup table */

    err = localfn_Cprime_lookup_table(
      &levels_sorted,
      delta_range,
      intrp_ptr,
      &Cprime);
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* form the sum of the log of the Cprime values
       for each delta (row) in the Cprime matrix */

    for ( row = 0; row < n_delta; row++ ){

      sum = 0.0;

      row_global = row * n_level;

      for ( col = 0; col < n_level; col++ ){
        sum += (double) n_coeff.data[ col ] *
          log( Cprime.data[ row_global + col ] );
      }

      sum_NjlogCprime.data[ row ] = sum;
    }

    /* set information in the wavelet energy structure */

    wave.n_sample       = n_sample;
    wave.biased         = biased;
    wave.energy         = &energy;
    wave.level          = &levels_sorted.mat.s32mat;
    wave.n_coeff        = &n_coeff;
    wave.n_level        = n_level;
    wave.n_level_usable = n_level;

    wave.delta_min =
      FDP_SDF_TABLE_DELTA_MIN + delta_start * FDP_SDF_TABLE_dDELTA;

    wave.stationary =
      ( estimator == (wav_fdp_estimator) WAV_FDP_MAXIMUM_LIKELIHOOD ) &&
      biased;

    /* bracket a three-point minimum of the
       reduced log-likelihood function */

    err = localfn_minimum_bracket(
      &Cprime,
      &wave,
      &sum_NjlogCprime,
      delta_start,
      &mle_delta_bracket,
      &mle_rll_bracket);
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* calculate the MLE of delta ... */

    /* ... if no bracket of a minimum
       values was found, then the first
       value in the rll_delta_bracket vector
       will be set to MUTIL_DOUBLE_EPSILON and the
       second value will be set to the boundary
       encountered during bracketing (this will be either
       FDP_SDF_TABLE_DELTA_MIN or FDP_SDF_TABLE_DELTA_MAX
       in the FDP_SDF_TABLE )  */

    if ( mle_delta_bracket.data[ 0 ] == MUTIL_DOUBLE_EPSILON ){
      *delta = mle_delta_bracket.data[ 1 ];
    }
    else{

      err = localfn_parabola_minimum_abscissa(
        &mle_delta_bracket,
        &mle_rll_bracket,
        delta );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    /* calculate corresponding MLE of the innovation variance ... */

    /* ... first calculate the mid-octave SDF value for
       an FD process with FD parameter equal to the MLE
       of delta. here we do not use the lookup table but
       call the appropriate function directly. */

    err = wavuniv_fdp_bandpass_variance(
      &levels_sorted,
      *delta,
      -1,
      intrp_ptr,
      &Cprime_mle );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* ... register the Cprime_mle matrix with
       the memory manager */

    err = memlist_member_register( &list, &Cprime_mle, MEMTYPE_MATUNIV );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* calculate the innovation variance */

    err = localfn_innovation_variance(
      &(Cprime_mle.mat.dblmat),
      &wave,
      0,
      &n_coeff_used_innovation_variance,
      innovation_variance );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* we do not have a method to calculate the variance of delta for
       MLE, so just set it to some unrealistic number here */

    *variance_delta = -1.0;
  }

  /* free the remaining registered memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( ":: Done with wavuniv_fdp_estimator_block()" );

  return MUTIL_ERR_OK;
}


/***********************************************/
/* FDP SDF FUNCTION DEFINITIONS                */
/***********************************************/

/* Mid-octave spectral density       */
/* function (SDF) estimation.        */
/* Function documented in wav_fdp.h  */
/* Written by William Constantine    */
mutil_errcode wavuniv_fdp_bandpass_variance(
  const univ_mat *levels,
  double          delta,
  sint32          n_sample,
  void           *intrp_ptr,
  univ_mat       *result )
{
  double         innovation_variance = 1.0;
  double        *pd_result;
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         n_count;
  sint32         n_level_max;
  sint32         n_level_min = 0;
  sint32         n_level3J = 0;
  sint32         n_level;
  sint32        *ps_levels;
  univ_mat       um_Cprime12;
  univ_mat       um_Cprime3J;
  univ_mat       um_levels3J;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_fdp_bandpass_variance()" );

  /* avoid lint warning */

  (void) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* initialize variables */

  n_level = MATUNIV_NELEM( levels );

  /* allocate space for the output */

  if ( n_sample > 0 ){
    err = matuniv_malloc_register( result, 1, n_level + 1,
      MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }
  else{
    err = matuniv_malloc_register( result, 1, n_level,
      MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* check the argument list */

  err = localfn_fdp_bandpass_arg_check(
    levels,
    n_sample,
    result );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* obtain the range of levels */

  err = mats32_range( &(levels->mat.s32mat), intrp_ptr,
                      &n_level_min, &n_level_max );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create pointers */

  ps_levels = levels->mat.s32mat.data;
  pd_result = result->mat.dblmat.data;

  /* in the case where the scaling band estimates
     are to be calculated (prompted by n_sample > 0),
     the levels vector must contain the elements
     1,2,3, ..., J where J is the maximum wavelet
     transform decomposition level */

  for ( i = 0; i < n_level; i++ ){

    if ( n_sample > 0 ){
      if ( ps_levels[ i ] != ( i+1 ) ){
        MUTIL_ERROR( "Levels vector must contain the elements "
                     "1,2,3, ..., J where J is the maximum wavelet "
                     "transform decomposition level");
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
    }

    /* count the number of levels > 2 */

    if ( ps_levels[ i ] > 2 )  n_level3J++;
  }

  /* if the minimum level is less than or equal to 2
     then (1) allocate a double vector with 2 elements,
     (2) calculate the mid-octave SDF values for
     levels 1 and 2, and (2) store the result in the
     return vector */

  if ( n_level_min <= 2 ){

    err = matuniv_malloc_register( &um_Cprime12, 2, 1,
                               MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = localfn_fdp_bandpass_levels12(
      delta,
      innovation_variance,
      intrp_ptr,
      &um_Cprime12);
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* fill in the results vector in the order specified by the levels
       vector. */

    for ( i = 0; i < n_level; i++ ){

      if ( ps_levels[ i ] < 3 ){
        pd_result[ i ] =
          um_Cprime12.mat.dblmat.data[ ps_levels[ i ] - 1 ];
      }

      /* Check for interrupts */

      if ( MUTIL_INTERRUPT( 3.0 * n_level, intrp_ptr ) ) {
        MUTIL_FREE_WARN ( memlist, &list );
        MUTIL_ERROR( "user interrupt" );
        return MUTIL_ERR_INTERRUPT;
      }
    }
  }

  /* if there are (n_level3J) levels greater than zero,
     then (1) allocate a double vector with n_level3J elements,
     (2) allocate a large_levels vector which contains
     those levels, (3) calculate the mid-octave SDF values for
     these levels and (2) store the result in the
     return vector */

  if ( n_level3J > 0 ){

    /* malloc space */

    err = matuniv_malloc_register( &um_Cprime3J, n_level3J, 1,
                               MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matuniv_malloc_register( &um_levels3J, n_level3J, 1,
                               MUTIL_SINT32, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* fill in the um_levels3J vector */

    n_count = 0;

    for ( i = 0; i < n_level; i++ ){

      if (ps_levels[ i ] > 2){
        um_levels3J.mat.s32mat.data[ n_count++ ] = ps_levels[ i ];
      }
    }

    /* calculate the mid-octave SDF values for wavelet decompoition
       levels greater than 2 */

    err = localfn_fdp_bandpass_levels3J( &um_levels3J,
                                         delta,
                                         innovation_variance,
                                         &um_Cprime3J );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* fill in the result vector */

    n_count = 0;

    for ( i = 0; i < n_level; i++ ){

      if ( ps_levels[ i ] > 2 ){
        pd_result[ i ] = um_Cprime3J.mat.dblmat.data[ n_count++ ];
      }
    }
  }

  /* append the average mid-octave SDF value for the
     octave corresponding to the scaling coefficients */

  if ( n_sample > 0 ){

    err = localfn_fdp_bandpass_scaling( result,
                                        delta,
                                        innovation_variance,
                                        n_sample );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* free the node corresponding to the
     outputs in the memory list but not the
     memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free allocated memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_fdp_bandpass_variance()" );

  return MUTIL_ERR_OK;
}


/***********************************************/
/* STATIC FDP ESTIMATION FUNCTION DEFINITIONS  */
/***********************************************/

/** Calculates the number of "usable" levels in a maximum overlap wavelet
 * transform for instantaneous FD estimates at a given point in time.
 * For unbiased instantaneous FD parameter estimates, the the number
 * of usable levels is defined as the number of wavelet decomposition
 * levels which contain an interior (zero phase shifted) wavelet coefficient
 * at a given location in time. For biased estimates, the number of usable
 * levels is simply the number of levels over which the FD estimates
 * are calculated.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_number_usable_levels( &interior, TRUE, time, &result );#
 * @return Standard mutils error/OK code.
 * @param interior   Pointer to a pre-allocated universal matrix of type
 *                   MUTIL\_SINT32 containing
 *                   the range of indices for the interior wavelet
 *                   coefficients at specific decomposition levels.
 *                   This interior matrix is 2xJ where J is the the total
 *                   number of specified levels. The first and second rows
 *                   contain the low and high indices, respectively.
 * @param biased     Boolean flag denoting biased or unbiased estimates.
 * @param time       Current time index.
 * @param result     Pointer to a signed integer to contain the result.
 *
 * @see wavuniv_transform_coefficient_boundaries
 * @private
 */
static mutil_errcode localfn_number_usable_levels(
  const univ_mat *interior,
  boolean         biased,
  sint32          time,
  sint32         *result)
{
  sint32  j;
  sint32  n_level;
  sint32 *ps_interior_high;
  sint32 *ps_interior_low;
  mutil_errcode err;

  MUTIL_TRACE( "Start localfn_number_usable_levels() ::" );

  /*** check interior ... */

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( interior, univ_mat, matuniv );

  /* ... for type MUTIL_SINT32 */

  if ( interior->type != MUTIL_SINT32 ){
    MUTIL_ERROR( "Interior matrix must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... for size */

  if ( MATUNIV_NROW( interior ) != 2 ){
    MUTIL_ERROR( "Interior matrix must have two rows." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* obtain the number of decomposition levels */

  n_level = MATUNIV_NCOL( interior );

  if ( !biased ){

    /* initial result */

    *result = 0;

    /* create pointers */

    ps_interior_low  = interior->mat.s32mat.data;
    ps_interior_high = interior->mat.s32mat.data + n_level;

    for ( j = 0; j < n_level; j++ ){

      if ( ( time >= ps_interior_low[ j ] ) &&
        ( time <= ps_interior_high[ j ] ) ) (*result)++;
    }
  }
  else *result = n_level;

  return MUTIL_ERR_OK;
}

/** Returns the range of indices for the interior wavelet
 * coefficients at the user specified decomposition levels.
 *
 * NOTE: levels can be skipped ( e.g. j = 1,3,5 ).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_interior_coefficient_index_range( &shift, &levels, 5, 8, 1024, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param shift          Pointer to a universal matrix set which
 *                       must contain two matrices, each of
 *                       type MUTIL\_SINT32 and of size
 *                       [1 x ( num\_levels * 2 )].
 *                       The contents of the matrix set should be (in order)
 *                       the zero phase shifts for the DWT and the MODWT.
 *                       MUTIL\_SINT32. The columns in each matrix correspond to
 *                       [W1 | ... | WJ | V1 | ... | VJ] where W and V are
 *                       the wavelet and scaling coefficient vectors,
 *                       respectively, and J = num\_levels is the number
 *                       of decomposition levels. A negative shift factor
 *                       implies an advance (circular shift to the left)
 *                       of the wavelet transform coefficient vectors W or V.
 * @param  levels        Pointer to a pre-allocated single-column or
 *                       single-row universal matrix of type MUTIL\_SINT32.
 *                       This argument contains the decomposition
 *                       levels. The levels can be given in any order, but
 *                       must be positive.
 * @param n_level_max    The maximum value in the levels vector.
 * @param filter_length  The length of the wavelet filter.
 * @param n_sample       The number of samples in the time series.
 * @param intrp_ptr      Interrupt pointer.
 * @param result         Pointer to a pre-allocated universal matrix of type
 *                       MUTIL\_SINT32. The result is a 2xJ matrix where J
 *                       is the the total number of specified levels. The
 *                       first and second rows contain the low and high
 *                       indices, respectively.
 *
 * @see wavuniv_filters_zero_phase
 * @private
 */
static mutil_errcode localfn_interior_coefficient_index_range(
  const mat_set  *shift,
  const univ_mat *levels,
  sint32          n_level_max,
  sint32          filter_length,
  sint32          n_sample,
  void           *intrp_ptr,
  univ_mat       *result)
{
  mat_set        boundary;
  memlist        list;
  mutil_errcode  err;
  sint32         j;
  sint32         level;
  sint32         n_level;
  sint32        *ps_boundary_start;
  sint32        *ps_interior_high;
  sint32        *ps_interior_length;
  sint32        *ps_interior_low;
  sint32        *ps_levels;
  sint32        *ps_shift;

  /* initialize memory list */

  MEMLIST_INIT( list );

  MUTIL_TRACE( "Start localfn_interior_coefficient_index_range() ::" );

  /* begin I/O checking */
  /* validate matrices  */

  /*** check shift ... */

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( shift, mat_set, matset );

  /* ...for type */

  if ( shift->mats[ 0 ].type != MUTIL_SINT32 ){
    MUTIL_ERROR( "The matrices in the zero phase shift "
      "matrix set must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /*** check levels ... */

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( levels, univ_mat, matuniv );

  /* ...for type */

  if ( levels->type != MUTIL_SINT32  ){
    MUTIL_ERROR( "Levels matrix must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  err = matuniv_validate( result );
  if ( err ) return err;

  if ( result->type != MUTIL_SINT32  ){
    MUTIL_ERROR( "Result matrix must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  if ( MATUNIV_NROW( result ) != 2 ){
    MUTIL_ERROR( "Result matrix must have two rows." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* obtain number of decomposition levels */

  n_level = MATUNIV_NELEM( levels );

  /* check sizes */

  if (MATUNIV_NCOL( &shift->mats[ 0 ] ) != ( 2 * n_level_max ) ){
    MUTIL_ERROR( "The matrices in the zero phase shift matrix set "
      "must have twice the number of columns "
      "as does the boundary matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if (MATUNIV_NCOL( result ) != n_level ){
    MUTIL_ERROR( "Result matrix must have the same number of columns "
      "as does elements in the levels vector." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* Calculate the location of the boundary and interior
     wavelet coefficients. */

  err = wavuniv_transform_coefficient_boundaries(
    n_level_max,
    filter_length,
    n_sample,
    ( wav_transform ) WAV_TRANSFORM_MODWT,
    intrp_ptr,
    &boundary );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* register boundary matrix with memory manager */

  err = memlist_member_register( &list, &boundary, MEMTYPE_MATSET );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create pointers */

  ps_boundary_start  = boundary.mats[ 0 ].mat.s32mat.data;
  ps_levels          = levels->mat.s32mat.data;
  ps_shift           = shift->mats[ 1 ].mat.s32mat.data;
  ps_interior_length = boundary.mats[ 2 ].mat.s32mat.data;
  ps_interior_low    = result->mat.s32mat.data;
  ps_interior_high   = result->mat.s32mat.data + n_level;

  for ( j = 0; j < n_level; j++ ){

    /* obtain level. decrement by one to coordinate
       with zero based indexing */

    level = ps_levels[ j ] - 1;

    ps_interior_low[ j ]  =
      (ps_boundary_start[level] + ps_shift[level] - 1) % n_sample;

    ps_interior_high[ j ] =
      ps_interior_low[ j ] + ps_interior_length[level] - 1;
  }

  /* free malloced space and corresponding linked list structures */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( ":: Done with localfn_interior_coefficient_index_range()" );

  return MUTIL_ERR_OK;
}

/** Returns the abscissa of the estimated minimum value
 * of (a set of three) points using a simple parabolic fit.
 * The three points are input as (x,y) rectangular coordinates.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_parabola_minimum_abscissa( &x, &y, &result );#
 * @return Standard mutils error/OK code.
 * @param x       Pointer to a pre-allocated double matrix containing
 *                the abscissa values. This matrix must be a single-column
 *                or single-row and must contain three elements.
 * @param y       Pointer to a pre-allocated double matrix containing
 *                the ordinate values. This matrix must be a single-column
 *                or single-row and must contain three elements.
 * @param result  Pointer to double value containing the result.
 *
 * @see wavuniv_fdp_estimator_instantaneous
 * @private
 */
static mutil_errcode localfn_parabola_minimum_abscissa(
  const double_mat *x,
  const double_mat *y,
  double           *result)
{
  double *pd_x;
  double *pd_y;
  double  denominator;
  double  dx10;
  double  dx12;
  double  dy10;
  double  dy12;
  mutil_errcode err;

  MUTIL_TRACE( "Start localfn_parabola_minimum_abscissa() ::" );

  /*** check x ... */

  /* ... for valid matrix structure and NULL pointer */
  LOCALDEF_CHECK_NULL_POINTER_FDP( x, double_mat, matdbl );

  /* ... for size */
  if ( !MATANY_IS_VEC( x ) ){
    MUTIL_ERROR( "Input matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if ( x->nelem != 3 ){
    MUTIL_ERROR( "Input matrix must have 3 elements" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /*** check y ... */

  /* ... for valid matrix structure and NULL pointer */
  LOCALDEF_CHECK_NULL_POINTER_FDP( y, double_mat, matdbl );

  /* ... for size */
  if ( !MATANY_IS_VEC( y ) ){
    MUTIL_ERROR( "Input matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if ( y->nelem != 3 ){
    MUTIL_ERROR( "Input matrix must have 3 elements" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* create pointers */

  pd_x = x->data;
  pd_y = y->data;

  /* calculate local differences */

  dx10 = pd_x[ 1 ] - pd_x[ 0 ];
  dx12 = pd_x[ 1 ] - pd_x[ 2 ];
  dy10 = pd_y[ 1 ] - pd_y[ 0 ];
  dy12 = pd_y[ 1 ] - pd_y[ 2 ];

  /* ensure that the three points are not colinear */

  denominator = dx10 * dy12 - dx12 * dy10;

  if ( (double) MUTIL_ABS( denominator ) < MUTIL_DOUBLE_EPSILON ){
    MUTIL_ERROR( "The three points given are colinear and " \
                 "thus have no minimum");
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  else *result = pd_x[ 1 ] -
    ( dx10 * dx10 * dy12 - dx12 * dx12 * dy10 ) /
    denominator / 2.0;

  MUTIL_TRACE( ":: Done with localfn_parabola_minimum_abscissa()" );

  return MUTIL_ERR_OK;
}

/** Brackets the minimum of a reduced log-likelihood
 * function used for instantaneous FD parameter estimation.
 * This function returns the abscissa and ordinate of a three point
 * bracketed minimum of the reduced log-likelihood function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_minimum_bracket( &Cprime, &wave, &sum_NjlogCprime, -10.0, &xmin, &ymin );#
 * @return Standard mutils error/OK code.
 * @param Cprime       Pointer to a pre-allocated double matrix containing
 *                     containing the average mid-octave SDF
 *                     values. The rows of this matrix correspond to
 *                     the FD parameter (delta) and the columns to
 *                     the wavelet decomposition levels. The size of this
 *                     matrix is D x J where D is the number of deltas
 *                     in the table and J is the number of wavelet
 *                     decomposition levels. It is acceptable for levels
 *                     to be skipped.
 * @param wave         Pointer to a \Ref{_wave_energy} structure containing
 *                     a wavelet coefficient energy vector and (other)
 *                     relevant information. All vectors in this structure
 *                     must correpond to the levels (columns) of the
 *                     Cprime matrix.
 * @param sum_NjlogCprime Pointer to a pre-allocated double matrix containing
 *                     the sum of the log of the Cprime matrix. The sum
 *                     is performed across the columns for each row in Cprime.
 *                     The size of this matrix is D x 1.
 * @param delta_start  The value of delta which corresponds to the first row
 *                     in Cprime.
 * @param xmin         A pointer to a pre-allocated double matrix with a
 *                     single-column or single-row containing exactly
 *                     three elements. This vector contains the abscissa
 *                     values of the bracketed minimum.
 * @param ymin         A pointer to a pre-allocated double matrix with a
 *                     single-column or single-row containing exactly
 *                     three elements. This vector contains the ordinate
 *                     values of the bracketed minimum.
 *
 * @see _wave_energy
 * @see wavuniv_fdp_estimator_instantaneous
 * @private
 */
static mutil_errcode localfn_minimum_bracket(
  const double_mat *Cprime,
  wave_energy      *wave,
  const double_mat *sum_NjlogCprime,
  sint32            delta_start,
  double_mat       *xmin,
  double_mat       *ymin)
{
  double R = 0.61803399;
  double C = 1.0 - R;

  double f1 = 0.0;
  double f2 = 0.0;
  double ftemp;
  double fxmin;
  mutil_errcode err;
  sint32 ax = 0;
  sint32 bx;
  sint32 cx;
  sint32 dx = 1;
  sint32 j;
  sint32 n_delta;
  sint32 x0;
  sint32 x1;
  sint32 x2;
  sint32 x3;
  sint32 xtemp;

  MUTIL_TRACE( "Start localfn_minimum_bracket() ::" );

  /* validate input matrices and check for NULL */

  LOCALDEF_CHECK_NULL_POINTER_FDP( Cprime, double_mat, matdbl );
  LOCALDEF_CHECK_NULL_POINTER_FDP( wave->energy, double_mat, matdbl );
  LOCALDEF_CHECK_NULL_POINTER_FDP( sum_NjlogCprime, double_mat, matdbl );

  /* check sizes */

  if ( xmin->nelem != 3 ){
    MUTIL_ERROR( "xmin matrix must have 3 elements" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( ymin->nelem != 3 ){
    MUTIL_ERROR( "xmin matrix must have 3 elements" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  err = matdbl_validate( xmin );
  if ( err ) return err;

  err = matdbl_validate( ymin );
  if ( err ) return err;

  /* obtain number of deltas in Cprime table */

  n_delta = Cprime->nrow;

  if ( sum_NjlogCprime->nelem != n_delta ){
    MUTIL_ERROR( "Number of elements in sum_NjlogCprime vector must " \
                 "be equal to the number of rows in the Cprime matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* initialize indices */

  bx = n_delta / 2;
  cx = n_delta - 1;
  x0 = ax;
  x3 = cx;

  /* find the side that is larger (the left
     or the right), and calculate the next
     abscissa value to compare with ... */

  if ( MUTIL_ABS( cx - bx ) > MUTIL_ABS( bx - ax ) ){
    x1  =  bx;
    x2  =  (sint32) ( MUTIL_ROUND( (double) bx + C * (double) ( cx - bx ) ) );
  }
  else {
    x2  =  bx;
    x1  =  (sint32) ( MUTIL_ROUND( (double) bx - C * (double) ( bx - ax ) ) );
  }

  /* evaluate function in the new interval */

  err = localfn_reduced_log_likelihood(
    Cprime,
    wave,
    sum_NjlogCprime,
    x1,
    &f1 );
  if (err) return err;

  err = localfn_reduced_log_likelihood(
    Cprime,
    wave,
    sum_NjlogCprime,
    x2,
    &f2 );
  if ( err ) return err;

  /* shift points to bracket towards a minimum */

  while ( MUTIL_ABS( x3 - x0 ) > dx ) {

    if ( f2 < f1 ) {

      SHIFT3( x0, x1, x2,
        (sint32) ( MUTIL_ROUND( R * (double) x1 + C * (double) x3 ) ) );

      err = localfn_reduced_log_likelihood(
        Cprime,
        wave,
        sum_NjlogCprime,
        x2,
        &ftemp );
      if ( err ) return err;

      SHIFT2( f1, f2, ftemp );
    }
    else{
      SHIFT3( x3, x2, x1,
            (sint32) ( MUTIL_ROUND( R * (double) x2 + C * (double) x0) ) );

      err = localfn_reduced_log_likelihood(
        Cprime,
        wave,
        sum_NjlogCprime,
        x1,
        &ftemp );
      if ( err ) return err;

      SHIFT2( f2, f1, ftemp );
    }
  }
  if ( f1 < f2 ) {
    xtemp = x1;
    fxmin = f1;
  }
  else{
    xtemp  = x2;
    fxmin = f2;
  }

  /* because we are using a lookup table version of golden section
     bracketing, we may be off to the left or right by one point
     in seeking a minimum. therefore, we should evaluate the points
     on each side of xtemp to establish a true abscissa value which
     locates the minimum                                       */

  if ( ( xtemp == 0 ) || ( xtemp == n_delta - 1 ) ){

    /* in this case a minimum was never bracketed.
       return a ridiculous value for the first
       abscissa value as a flag to the caller.
       set the second value as the value of delta
       where we ran off the Cprime table */

    xmin->data[ 0 ] = MUTIL_DOUBLE_EPSILON;

    if ( xtemp == 0 ) xmin->data[ 1 ] = FDP_SDF_TABLE_DELTA_MIN;
    else xmin->data[ 1 ] = FDP_SDF_TABLE_DELTA_MAX;
  }
  else{

    for ( j = -1; j <= 1; j++ ){

      xmin->data[ j + 1 ] =
        ( xtemp + j + delta_start ) * FDP_SDF_TABLE_dDELTA +
        FDP_SDF_TABLE_DELTA_MIN;

      if ( j == 0 ){
        ymin->data[ j + 1 ] = fxmin;
      }
      else{
        err = localfn_reduced_log_likelihood(
          Cprime,
          wave,
          sum_NjlogCprime,
          xtemp + j,
          &ftemp );
        if ( err ) return err;

        ymin->data[ j + 1 ] = ftemp;
      }
    }
  }

  MUTIL_TRACE( ":: Done with localfn_minimum_bracket()" );

  return MUTIL_ERR_OK;
}

/** Returns the square of the zero phase shifted MODWT wavelet
 * coefficients for levels j in [1,J] corresponding
 * to a particular time. Here, J is the number of wavelet
 * transform decomposition levels.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_wavelet_energy_instantaneous( &modwt, &shift, &levels, &interior, n_level_usable, biased, time, dof_order, WAV_FDP_LEAST_SQUARES, &n_coeff, &Wtmodsqr );#
 * @return Standard mutils error/OK code.
 * @param  modwt       Pointer to a matrix set containing the MODWT
 *                     in n\_level + 1 pre-allocated single-row universal
 *                     matrices of type MUTIL\_DOUBLE. The order of the
 *                     matrices must be [W1 , W2 , ... , WJ,  VJ].
 * @param shift        Pointer to a universal matrix set which
 *                     must contain two matrices, each of
 *                     type MUTIL\_SINT32 and of size
 *                     [1 x ( num\_levels * 2 )].
 *                     The contents of the matrix set should be (in order)
 *                     the zero phase shifts for the DWT and the MODWT.
 *                     MUTIL\_SINT32. The columns in each matrix correspond to
 *                     [W1 | ... | WJ | V1 | ... | VJ] where W and V are
 *                     the wavelet and scaling coefficient vectors,
 *                     respectively, and J = num\_levels is the number
 *                     of decomposition levels. A negative shift factor
 *                     implies an advance (circular shift to the left)
 *                     of the wavelet transform coefficient vectors W or V.
 * @param  levels      Pointer to a pre-allocated single-column or
 *                     single-row universal matrix of type MUTIL\_SINT32.
 *                     This argument contains the decomposition
 *                     levels. The levels can be given in any order, but
 *                     must be positive.
 * @param interior     Pointer to a pre-allocated universal matrix of type
 *                     MUTIL\_SINT32 containing
 *                     the range of indices for the interior wavelet
 *                     coefficients at specific decomposition levels.
 *                     This interior matrix is 2xJ where J is the the total
 *                     number of specified levels. The first and second rows
 *                     contain the low and high indices, respectively.
 * @param n_level_usable Integer denoting the number of modwt coefficients
 *                     to use in forming the result.
 * @param biased       Boolean denoting biased or unbiased estimation of the
 *                     instantaneous FD parameters.
 * @param time         The current time index.
 * @param dof_order    The order (K) of the the degrees of freedom used to
 *                     estimate the FD parameters. The actual degrees
 *                     of freedom is equal to 2*(K+1).
 * @param estimator    The type of instantantaneous estimator used
 *                     to form the FD parameter estimates. Acceptable
 *                     types are WAV\_FDP\_MAXIMUM\_LIKELIHOOD and
 *                     WAV\_FDP\_LEAST\_SQUARES.
 * @param n_coeff      A pointer to a pre-allocated sint32 matrix with a
 *                     single-column or single-row containing the number
 *                     of wavelet coefficients used to for the energies
 *                     returned in the Wtmodsqr vector.
 * @param Wtmodsqr     A pointer to a pre-allocated double matrix with a
 *                     single-column or single-row containing the square
 *                     of the zero phase shifted MODWT coefficients
 *                     corresponding to the given time.
 *
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_filters_zero_phase
 * @see wavuniv_transform_coefficient_boundaries
 * @see wavuniv_fdp_estimator_instantaneous
 * @private
 */
static mutil_errcode localfn_wavelet_energy_instantaneous(
  const mat_set     *modwt,
  const mat_set     *shift,
  const univ_mat    *levels,
  const univ_mat    *interior,
  sint32             n_level_usable,
  boolean            biased,
  sint32             time,
  sint32             dof_order,
  wav_fdp_estimator  estimator,
  sint32_mat        *n_coeff,
  double_mat        *Wtmodsqr )
{
  double         Wshift;
  double         scale_factor;
  double         sumWtmodsqr;
  double        *pd_modwt;
  mutil_errcode  err;
  sint32         dofs_used;
  sint32         j;
  sint32         k;
  sint32         level;
  sint32         lower_bound;
  sint32         n_level;
  sint32         n_level_modwt;
  sint32         n_sample;
  sint32         scale;
  sint32         t;
  sint32         tdof;
  sint32         upper_bound;
  sint32        *ps_interior_high;
  sint32        *ps_interior_low;
  sint32        *ps_levels;
  sint32        *ps_shift;

  MUTIL_TRACE( "Start localfn_wavelet_energy_instantaneous() ::" );

  /*** check time series ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( modwt, mat_set, matset );

  /* ... for type MUTIL_DOUBLE */

  if ( modwt->mats[ 0 ].type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "MODWT matrix set must contain universal "
      "matrices of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /*** check shift ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( shift, mat_set, matset );

  /* ... for type MUTIL_SINT32 */

  if ( shift->mats[ 0 ].type != MUTIL_SINT32 ){
    MUTIL_ERROR( "Matrices in the zero phase shift matrix set "
      "must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /*** check levels ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( levels, univ_mat, matuniv );

  /* ... for type MUTIL_SINT32 */

  if ( levels->type != MUTIL_SINT32 ){
    MUTIL_ERROR( "Levels matrix must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /*** check interior ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( interior, univ_mat, matuniv );

  /* ... for type MUTIL_SINT32 */

  if ( interior->type != MUTIL_SINT32 ){
    MUTIL_ERROR( "Interior matrix must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* validate output */

  err = matdbl_validate( Wtmodsqr );
  if ( err ) return err;

  if ( !MATANY_IS_VEC( Wtmodsqr ) ){
    MUTIL_ERROR( "Energy matrix must be a single column or a single row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  err = mats32_validate( n_coeff );
  if ( err ) return err;

  if ( !MATANY_IS_VEC( n_coeff ) ){
    MUTIL_ERROR( "Number of coefficients matrix must be a "
      "single column or a single row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* check sizes:
     n_level_modwt contains the total number of decomposition
     levels used to develop a full MODWT. n_level is the
     total number of levels that are to be extracted from the
     MODWT matrix. accordingly, n_level <= n_level_modwt */

  n_level       = MATUNIV_NELEM( levels );
  n_level_modwt = modwt->nelem - 1;
  n_sample      = MATUNIV_NELEM( &modwt->mats[ 0 ] );

  if ( n_level > n_level_modwt ){
    MUTIL_ERROR( "Number of levels exceeds that available in " \
                 "the full MODWT." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if( MATUNIV_NCOL( &shift->mats[ 0 ] ) != ( 2 * n_level_modwt ) ){
    MUTIL_ERROR( "Number of columns in matrices in zero phase shift "
      "matrix set must be twice the number of "
      "decomposition levels." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

 if ( Wtmodsqr->nelem != n_level ){
   MUTIL_ERROR( "Number of elements in Wtmodsqr vector must be " \
                "equal to the number of decomposition levels." );
   return MUTIL_ERR_ILLEGAL_SIZE;
 }
 if ( n_coeff->nelem != n_level ){
   MUTIL_ERROR( "Number of elements in n_coeff vector must be " \
                "equal to the number of decomposition levels." );
   return MUTIL_ERR_ILLEGAL_SIZE;
 }

  if ( MATUNIV_NCOL( interior ) != n_level ){
    MUTIL_ERROR( "Number of columns in interior matrix must be " \
                 "equal to the number of decomposition levels." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* create pointers */

  ps_levels        = levels->mat.s32mat.data;
  ps_shift         = shift->mats[ 1 ].mat.s32mat.data;
  ps_interior_low  = interior->mat.s32mat.data;
  ps_interior_high = interior->mat.s32mat.data + n_level;

  for( j = 0; j < n_level; j++ ){

    /* obtain the desired level from which to extract
       the wavelet coefficients. we decrement by one
       because the levels are unit-based while C
       indexing is zero-based, so level 1 should
       be accessing the zeroth element, and so on. */

    level = ps_levels[ j ] - 1;

    pd_modwt = modwt->mats[ level ].mat.dblmat.data;

    /* calculate the scale of the current level */

    /* NOTE: Bill, could/should this be 1<<level ??? */
    scale = (sint32) MUTIL_POW( 2, level );

    /* add shift to approximate zero phase */

    t = time - ps_shift[ level ];

    /* impose periodic extension
       through mod operations */

    while ( t < 0 ) t+= n_sample;
    while ( t >= n_sample ) t-= n_sample;

    /* form lower and upper bounds of interior wavelet coefficients.
       the "interior" matrix contains the index range of the
       interior wavelet coefficient for each level specified
       by the "levels" vector, i.e. the range given by
       [lower_bound, upperbound] is post-shifted. however, the
       time variable "t" has not shifted for */

    if ( !biased ){
      lower_bound = ps_interior_low[ j ] - ps_shift[ level ];
      upper_bound = ps_interior_high[ j ] - ps_shift[ level ];
    }
    else{
      lower_bound = 0;
      upper_bound = n_sample - 1;
    }

    /* initialize summation variables */

    sumWtmodsqr = 0.0;
    dofs_used   = 0;

    /* form a dof-point average in time of the
       (zero phase shifted) squared wavelet
       coefficients */

    if ( j < n_level_usable ){

      for ( tdof = - dof_order; tdof <= dof_order; tdof++ ){

        k = t + ( tdof * scale );

        if ( ( k  >= lower_bound ) && ( k <= upper_bound ) ){

          Wshift       = pd_modwt[ k ];
          sumWtmodsqr += Wshift * Wshift;

          dofs_used++;
        }
      }
    }

    /* normalize by the number of dofs actually used */

    if ( dofs_used > 0 ) sumWtmodsqr /= (double) dofs_used;

    /* the Wtmodsqr is scaled by 2^level for the MLE case */

    if ( estimator == WAV_FDP_MAXIMUM_LIKELIHOOD ){
      scale_factor = MUTIL_POW( 2.0, (double) (level + 1) );
    }
    else scale_factor = 1.0;

    /* set level-dependent data in output */

    Wtmodsqr->data[ j ]  = sumWtmodsqr * scale_factor;
    n_coeff->data[ j ]   = 1;
  }

  MUTIL_TRACE( ":: Done with localfn_wavelet_energy_instantaneous()" );

  return MUTIL_ERR_OK;
}

/** Energy of the wavelet coefficients at each decomposition level.
 * In addition, the number of coefficients in each level and (in
 * the case where the estimator is a WLSE) the EDOF of a specified
 * mode for each level is returned.
 *
 * NOTE: The length of all vectors returned by this function is equal
 * to that of the levels input vector. The levels vector may
 * contain skipped level values (e.g. levels = {2, 4, 5}) and
 * all outputs from this function coordinate with those levels
 * specified by the user.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_wavuniv_transform_energy_block( &time_series, &levels, filter_type, filter_length, biased, estimator, edof_mode, &sdf, , intrp_ptr, &edof, &n_coeff, &energy );#
 * @return Standard mutils error/OK code.
 * @param time_series    Pointer to a single column universal matrix
 *                       of type MUTIL\_DOUBLE containing the time series
 *                       to analyze.
 * @param levels         Pointer to a pre-allocated single-column
 *                       universal matrix of type MUTIL\_SINT32.
 *                       This argument contains the decomposition
 *                       levels. The levels can be given in any order, but
 *                       must be positive.
 * @param filter_type    The type of Daubechies filter.
 * @param filter_length  The length of the wavelet filter.
 * @param biased         Boolean flag denoting biased or unbiased
 *                       estimates. Biased estimates are those which
 *                       use all available levels in calculating
 *                       the FD model parameters. Unbiased estimates
 *                       are calculated with only those wavelet
 *                       coefficients not subject to circular filter
 *                       operations, i.e. only the interior wavelet
 *                       coefficients are used in calculating unbiased
 *                       estimates.
 * @param estimator      The method to estimate the FD model parameters.
 *                       This argument is an enumerated type
 *                       \Ref{_wav_fdp_estimator}.
 * @param edof_mode      The equivalent degrees of freedom (EDOF) mode.
 * @param intrp_ptr      Pointer for implementation of interrupt checking.
 * @param n_coeff        A pointer to a sint32 matrix with a
 *                       single-column or single-row containing the number
 *                       of wavelet coefficients used to for the energies
 *                       returned in the energy vector. This vector
 *                       is automatically allocated by this function.
 *                       In the case where the estimator is
 *                       WAV\_FDP\_MAXIMUM\_LIKELIHHOD and the biased
 *                       argument is TRUE
 * @param energy         A pointer to a double matrix with a
 *                       single-column or single-row containing the
 *                       energy of all coefficients in each level of
 *                       a wavelettransform. This vector
 *                       is automatically allocated by this function.
 *
 * @see wavuniv_variance
 * @see wavuniv_transform_coefficient_boundaries
 * @see wavuniv_fdp_estimator_block
 * @private
 */
static mutil_errcode localfn_wavuniv_transform_energy_block(
  const univ_mat    *time_series,
  const univ_mat    *levels,
  wav_filter_type    filter_type,
  sint32             filter_length,
  boolean            biased,
  wav_fdp_estimator  estimator,
  sint32             edof_mode,
  const univ_mat    *sdf,
  void              *intrp_ptr,
  double_mat        *edof,
  sint32_mat        *n_coeff,
  double_mat        *energy )
{
  boolean        calculate_edof;
  boolean        normalize_filters;
  boolean        stationary;
  double         Wtsqr;
  double         num_ops = 0.0;
  double         sum_energy;
  double        *Wj;
  double        *pd_edof;
  double        *pd_var;
  mat_set        boundary;
  mat_set        edof_all;
  mat_set        filters;
  mat_set        transform;
  memlist        list;
  mutil_errcode  err;
  sint32         ii = 0;
  sint32         j;
  sint32         jj;
  sint32         n_level;
  sint32         n_level_max = 0;
  sint32         n_level_min;
  sint32         n_sample;
  sint32         sum_Lj;
  sint32         t;
  sint32         t_interior;
  sint32        *Lj;
  sint32        *Nj;
  univ_mat       interior;
  univ_mat       n_coeff_block_unbiased;
  univ_mat       var_block_unbiased;
  wav_transform  transform_type;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start localfn_wavuniv_transform_energy_block()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check the input arguments */

  /*** check time series ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( time_series, univ_mat, matuniv );

  /*** check levels ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( levels, univ_mat, matuniv );

  /* ... for type MUTIL_SINT32 */

  if ( levels->type != MUTIL_SINT32 ){
    MUTIL_ERROR( "Levels matrix must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* initialize variables */

  n_sample = MATUNIV_NELEM( time_series );
  n_level  = MATUNIV_NELEM( levels );

  err = mats32_range( &(levels->mat.s32mat), intrp_ptr,
                      &n_level_min, &n_level_max );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* set EDOF boolean:

     if the estimator is an MLE, then
     we do not need to calculate the EDOF since they are
     only used in WLSE of delta calculations */

  calculate_edof =
    ( estimator == (wav_fdp_estimator) WAV_FDP_LEAST_SQUARES );

  /* set stationary boolean:

     if we are calculating an MLE of delta and the biased
     boolean is TRUE, it means that we are calculating
     a blocked stationary MLE of delta using a full DWT  */

  stationary =
    ( estimator == (wav_fdp_estimator) WAV_FDP_MAXIMUM_LIKELIHOOD ) &&
    biased;

  if ( calculate_edof ){

    err = matdbl_malloc_register( edof, n_level, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    transform_type = (wav_transform) WAV_TRANSFORM_MODWT;
    normalize_filters = (boolean) TRUE;
  }
  else{

    transform_type = (wav_transform) WAV_TRANSFORM_DWT;
    normalize_filters = (boolean) FALSE;
  }

  if ( stationary ){

    err = mats32_malloc_register( n_coeff, n_level + 1, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( energy, n_level + 1, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }
  else{

    err = mats32_malloc_register( n_coeff, n_level, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( energy, n_level, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* obtain the wavelet and scaling filters
     and register them with the memory manager */

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
      time_series, &filters, n_level_max, intrp_ptr, &transform );
  }
  else if ( transform_type == (wav_transform) WAV_TRANSFORM_DWT ){

    err = wavuniv_transform_discrete_wavelet_convolution(
      time_series, &filters, n_level_max, intrp_ptr, &transform );
  }
  MEMLIST_FREE_ON_ERROR( err, &list );

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
    n_level_max, filter_length, n_sample,
    transform_type, intrp_ptr, &boundary );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &boundary, MEMTYPE_MATSET );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* assign pointers */

  Lj = boundary.mats[ 2 ].mat.s32mat.data;
  Nj = boundary.mats[ 4 ].mat.s32mat.data;

  /* allocate space and fill variables for the edof output */

  if ( calculate_edof ){

    err = matuniv_malloc_register(
      &var_block_unbiased, 1, n_level, MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matuniv_malloc_register(
      &n_coeff_block_unbiased, 1, n_level, MUTIL_SINT32, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* initialize index variables */

    ii     = 0;
    sum_Lj = 0;

    /* calculate the total number of interior */

    for ( j = 0; j < n_level; j++ ){

      jj = levels->mat.s32mat.data[ j ] - 1;

      /* here we allow the setting of
      "interior" wavelet coefficients
      even though (for the biased case)
      not all coefficients at a given level
      are interior wavelet coefficients. this
      is meant as a work-around for calculating
      biased WLSE, since the statistics for
      them have not been worked out as of yet */

      if ( biased ){

        sum_Lj += Nj[ jj ];
      }
      else{
        sum_Lj += Lj[ jj ];
      }
    }

    /* allocate memory to store a concatenated
       collection of interior wavelet coefficients */

    err = matuniv_malloc_register(
      &interior, 1, sum_Lj, MUTIL_DOUBLE, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* calculate the energy of the
     wavelet transform at each level. */

  for ( j = 0; j < n_level; j++ ){

    /* extract the desired decompostion level */

    jj = levels->mat.s32mat.data[ j ] - 1;

    /* assign pointer to current level's
       wavelet transform coefficients */

    Wj = transform.mats[ jj ].mat.dblmat.data;

    /* initialize variables */

    sum_energy = ( double ) 0.0;

    t_interior = Nj[ jj ] - Lj[ jj ];

    /* form time dependent wavelet variance estimates */

    for ( t = 0; t < Nj[ jj ]; t++ ){

      Wtsqr = Wj[ t ] * Wj[ t ];

      if ( biased ){
        sum_energy += Wtsqr;

        /* here we allow the setting of
        "interior" wavelet coefficients
        even though (for the biased case)
        not all coefficients at a given level
        are interior wavelet coefficients. this
        is meant as a work-around for calcualting
        biased WLSE, since the statistics for
        them have not been worked out as of yet */

        if ( calculate_edof ){
          interior.mat.dblmat.data[ ii++ ] = Wj[ t ];
        }
      }
      else{

        if ( t >= t_interior ){

          sum_energy += Wtsqr;

          /* store the interior wavelet coefficients
          for edof calcluation */

          if ( calculate_edof ){
            interior.mat.dblmat.data[ ii++ ] = Wj[ t ];
          }
        }
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

    energy->data[ j ] = sum_energy;

    /* set number of wavelet coefficients
       at current decompostion level */

    if ( biased ){

      n_coeff->data[ j ] = Nj[ jj ];
    }
    else{

      n_coeff->data[ j ] = Lj[ jj ];
    }

  } /* end of loop over level index */

  /* set scaling coefficient energy and
     number of coefficients if a stationary
     MLE is being performed. in this case,
     a full DWT has been performed and there
     can be only one scaling coefficient */

  if ( stationary ){

    Wj = transform.mats[ n_level ].mat.dblmat.data;

    energy->data[ n_level ]  = Wj[ 0 ] * Wj[ 0 ];
    n_coeff->data[ n_level ] = (sint32) 1;
  }

  /* calculate the edof */

  if ( calculate_edof ){

    /* form biased blocked wavelet variance estimates */

    for ( j = 0; j < n_level; j++ ){

      /* we can only get here if we are
      using the MODWT AND we have
      calculated UNBIASED wavelet
      energies. so, we develop the
      unbiased wavelet variance
      estimates accordingly */

      /* set pointers */

      pd_var = var_block_unbiased.mat.dblmat.data;

      n_coeff_block_unbiased.mat.s32mat.data[ j ] =
        n_coeff->data[ j ];

      /* define normalization factor */

      if ( n_coeff->data[ j ] > 0 ){

        pd_var[ j ] =
          energy->data[ j ] / n_coeff->data[ j ];
      }
      else{

        MUTIL_WARN( "There are no interior MODWT wavelet "
          "coefficients to use in calculating "
          "an unbiased wavelet variance estimate "
          "for the current decomposition level." );
        pd_var[ j ] = (double) -1.0;
      }
    }

    /* call the edof function */

    err = wavuniv_variance_edof(
      &interior,
      &n_coeff_block_unbiased,
      &var_block_unbiased,
      levels,
      sdf,
      filter_type,
      filter_length,
      intrp_ptr,
      &edof_all );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_register( &list, &edof_all, MEMTYPE_MATSET );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* set pointer to propoer edof matrix and fill */

    pd_edof = edof_all.mats[ edof_mode - 1 ].mat.dblmat.data;

    for ( j = 0; j < n_level; j++ ){

      /* obtain the current level and
      convert to base zero index */

      edof->data[ j ] = pd_edof[ j ];
    }
  }

  /* remove nodes corresponding to
     output from the memory list */

  err = memlist_member_unregister( energy, &list);
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_unregister( n_coeff, &list);
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( calculate_edof ){

    err = memlist_member_unregister( edof, &list);
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* free memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_wavuniv_transform_energy_block()" );

  return MUTIL_ERR_OK;
}


/** Estimates the innovation variance for an FD process.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_innovation_variance( &Cprime, &wave, row_delta, n_energy, n_level_usable, &result );#
 * @return Standard mutils error/OK code.
 * @param Cprime       Pointer to a pre-allocated double matrix containing
 *                     containing the average mid-octave SDF
 *                     values. The rows of this matrix correspond to
 *                     the FD parameter (delta) and the columns to
 *                     the wavelet decomposition levels. The size of this
 *                     matrix is D x J where D is the number of deltas
 *                     in the table and J is the number of wavelet
 *                     decomposition levels. It is acceptable for levels
 *                     to be skipped.
 * @param wave         Pointer to a \Ref{_wave_energy} structure containing
 *                     a wavelet coefficient energy vector and (other)
 *                     relevant information. All vectors in this structure
 *                     must correpond to the levels (columns) of the
 *                     Cprime matrix.
 * @param row_delta    An integer denoting the row of the Cprime matrix.
 *                     Each row of this table corresponds to specific
 *                     value of the FD parameter (delta).
 * @param n_coeff_used Pointer to an sint32 value which (upon return) will
 *                     contain the total number of wavelet coefficients used
 *                     to calculate the total energy of the wavelet
 *                     coefficients. This factor is used within the code to
 *                     normalize the innovation variance.
 * @param innovation_variance Pointer to a double value which (upon return)
 *                     contain the estimated innovation variance.
 *
 * @see _wave_energy
 * @see localfn_reduced_log_likelihood
 * @see wavuniv_fdp_estimator_instantaneous
 * @private
 */
static mutil_errcode localfn_innovation_variance(
  const double_mat *Cprime,
  wave_energy      *wave,
  sint32            row_delta,
  sint32           *n_coeff_used,
  double           *innovation_variance)
{
  double        *pd_Cprime;
  double         variance_fdp;
  double         Cprime_scaling;
  double         delta;
  sint32         j;
  sint32         n_level;
  sint32         n_level_usable;
  mutil_errcode  err;

  MUTIL_TRACE( "Start localfn_innovation_variance() ::" );

  /*** check Cprime ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( Cprime, double_mat, matdbl );

  /*** check wave energy matrices ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( wave->energy, double_mat, matdbl );
  LOCALDEF_CHECK_NULL_POINTER_FDP( wave->n_coeff, sint32_mat, mats32 );

  /* obtain level information:

     NOTE:

     n_level represents the number of
     levels set by the user (outside
     of this function). It does NOT
     necessarily represent the total
     number of levels used in the original
     wavelet decomposition. If levels
     were skipped, e.g. levels = {2,4,5}
     the n_level = 3 here. The value set
     for the n_level_usable variable is
     also relative to the levels vector so
     that if n_level_usable = 2 and
     levels = {2,4,5}, only the information
     set for levels 2 and 4 is used to calculate
     the innovation variance. The Cprime matrix
     is also made relative to the original
     levels vector so that the number of columns
     in the Cprime matrix should equal n_level.
  */

  n_level        = wave->n_level;
  n_level_usable = wave->n_level_usable;

  /* validate sizes of input */

  if ( Cprime->ncol != n_level){
    MUTIL_ERROR( "Number of columns in Cprime matrix must be equal " \
                 "to the number of decomposition levels." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( n_level_usable > n_level ){
    MUTIL_ERROR( "The number of usable levels cannot exceed the " \
                 "number of levels in the wavelet energy vector." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* validate components of the wave_energy structure
     in the case where a blocked stationary MLE of
     delta method is being used. here, the lengths of
     the energy and n_coeff vectors (packed into the
     wave_energy structure as an input to this function)
     need to be n_level + 1 so that information regarding the
     scaling coefficients may be included. By consequence
     of design, the length of the scaling coefficient
     vector should always be unity since it is expected
     that (in the case of a blocked stationary MLE of delta)
     a full wavelet decomposition has been performed.
  */

  if ( wave->stationary ){

    if ( ( ( wave->energy )->nelem != ( n_level + 1 ) ) ||
      ( ( wave->n_coeff )->nelem != ( n_level + 1 ) ) ){

      MUTIL_ERROR( "The length of the energy and n_coeff "
        "vectors in the wave_energy structure "
        "must be equal to J + 1 (where J is the "
        "number of decomposition levels) for "
        "blocked stationary MLEs of delta." );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    if ( ( wave->n_coeff )->data[ ( wave->n_coeff )->nelem - 1 ] !=
      (sint32) 1 ){

      MUTIL_ERROR( "The final value of the n_coeff "
        "vector in the wave_energy structure "
        "must be equal to J + 1 (where J is the "
        "number of decomposition levels) for "
        "blocked stationary MLEs of delta." );
      return MUTIL_ERR_ILLEGAL_SIZE;

    }
  }

  /* create pointers */

  pd_Cprime = &( Cprime->data[ row_delta * n_level ] );

  /* initialize variables */

  *innovation_variance = 0.0;
  *n_coeff_used = 0;

  /* calculate innovation variance */

  for( j = 0; j < n_level_usable; j++){

    *innovation_variance += ( wave->energy )->data[ j ] / pd_Cprime[ j ];
    *n_coeff_used        += ( wave->n_coeff )->data[ j ];
  }

  /* include the energy form the (single) scaling coefficient
     if a stationary blocked MLE of delta method is used. */

  if ( wave->stationary ){

  /* first disclose the value of the
  FD parameter. if it is not in the
  stationary regime, then we cannot
  use the scaling coefficient energy
    in estimating the innovation variance */

    delta = wave->delta_min +
      row_delta * FDP_SDF_TABLE_dDELTA;

    if ( delta < 0.5 ){

    /* develop the theoretical variance of an
       FD process (the following result
       represents an FD process variance
       implicitly divided by the FD
       innovation variance)  */

      variance_fdp = wave->n_sample *
        MUTIL_GAMMA( 1.0 - 2.0 * delta ) /
        MUTIL_POW( MUTIL_GAMMA( 1.0 - delta ), 2.0 );

      /* calculate the mid-octave SDF value
        of an FD process for the octave
        associated with the scaling
        coefficients */

      Cprime_scaling = variance_fdp - *innovation_variance;

      /* the mid-octave value for the
        scaling coefficients cannot be negative.
        if it is, reset it to unity so it
        will have no impact on the innvotion
        variance or reduced log-likelihood */

      if ( Cprime_scaling > 0.0 ){

	/* add the scaling energy to the current
	   innovation variance summation */

	*innovation_variance +=
	  ( wave->energy )->data[ n_level ] / Cprime_scaling;

        /* increment the number of coefficients used
           in the calculation of the innovation variance */

	(*n_coeff_used)++;

	/* pack the Cprime_scaling value into the wave_energy
	   structure for later use in calculating the
	   resuced log-likelihood */
      }

      wave->Cprime_scaling = Cprime_scaling;
    }
  }

  /* calculate the innovation variance */

  *innovation_variance /= (double) *n_coeff_used;

  MUTIL_TRACE( ":: Done with localfn_innovation_variance()" );

  return MUTIL_ERR_OK;
}

/** Given a range of the FD parameter (delta), this function calculates
 * the corresponding range of rows in the mid-octave
 * SDF table (Cprime) for a FD process. Given the range of delta
 * (say (DMIN, DMAX)) specified by the user,
 * this function returns (1) the starting row of
 * the Cprime table which most closely matches DMIN,
 * (2) the ending row of the Cprime table which
 * most closely matches DMAX and (3) the number of
 * rows (different deltas) in the Cprime table spanned
 * by DMIN through DMAX.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_delta_indices( &delta_range, &start, &end, &n_delta );#
 * @return Standard mutils error/OK code.
 * @param delta_range  Pointer to single-row or single-column
 *                     universal matrix of type MUTIL\_DOUBLE
 *                     with two elements and containing the desired
 *                     range of delta in (delta\_low, delta\_high) format.
 * @param start        The row of Cprime corresponding to delta\_low.
 * @param end          The row of Cprime corresponding to delta\_high.
 * @param n_delta      The number of rows in Cprime from delta\_low to
 *                     delta\_high.
 *
 * @see localfn_Cprime_lookup_table
 * @see wavuniv_fdp_estimator_instantaneous
 * @private
 */
static mutil_errcode localfn_delta_indices(
  const univ_mat *delta_range,
  sint32         *start,
  sint32         *end,
  sint32         *n_delta)
{
  double  delta_min = 0.0;
  double  delta_max = 0.0;
  sint32  offset;
  void   *intrp_ptr = NULL;
  mutil_errcode err;

  MUTIL_TRACE( "Start localfn_delta_indices() ::" );

  /*** check delta_range ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( delta_range, univ_mat, matuniv );

  err = matdbl_range( &(delta_range->mat.dblmat),
    intrp_ptr, &delta_min, &delta_max );
  if ( err ) return err;

  /* make sure values are within table range */

  if ( ( delta_min < FDP_SDF_TABLE_DELTA_MIN ) ||
    ( delta_min > FDP_SDF_TABLE_DELTA_MAX ) ||
    ( delta_max < FDP_SDF_TABLE_DELTA_MIN ) ||
    ( delta_max > FDP_SDF_TABLE_DELTA_MAX ) ){

    MUTIL_ERROR( "Specified delta range exceeds table values." );
    return MUTIL_ERR_OUT_OF_BOUNDS;
  }

  offset = (sint32) COERCE( FDP_SDF_TABLE_DELTA_MIN / FDP_SDF_TABLE_dDELTA );
  *start = (sint32) ( COERCE( delta_min / FDP_SDF_TABLE_dDELTA ) - offset );
  *end   = (sint32) ( COERCE( delta_max / FDP_SDF_TABLE_dDELTA ) - offset );

  *n_delta = *end - *start + 1;

  MUTIL_TRACE( ":: Done with localfn_delta_indices()" );

  return MUTIL_ERR_OK;
}

/** Reduced log-likelihood function evaluation for an MLE of delta.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_reduced_log_likelihood( &Cprime, &wave, &sum_NjlogCprime, delta_index, &result );#
 * @return Standard mutils error/OK code.
 * @param Cprime       Pointer to a pre-allocated double matrix containing
 *                     containing the average mid-octave SDF
 *                     values. The rows of this matrix correspond to
 *                     the FD parameter (delta) and the columns to
 *                     the wavelet decomposition levels. The size of this
 *                     matrix is D x J where D is the number of deltas
 *                     in the table and J is the number of wavelet
 *                     decomposition levels. It is acceptable for levels
 *                     to be skipped.
 * @param wave         Pointer to a \Ref{_wave_energy} structure containing
 *                     a wavelet coefficient energy vector and (other)
 *                     relevant information. All vectors in this structure
 *                     must correpond to the levels (columns) of the
 *                     Cprime matrix.
 * @param sum_NjlogCprime Pointer to a pre-allocated double matrix containing
 *                     the sum of the log of the Cprime matrix. The sum
 *                     is performed across the columns for each row in Cprime.
 *                     The size of this matrix is D x 1.
 * @param delta_index  The row of the Cprime table which corresponds to the
 *                     current estimate of delta.
 * @param result       A pointer to a pre-allocated double matrix with a
 *                     single-column or single-row containing exactly
 *                     three elements. This vector contains the ordinate
 *                     values of the bracketed minimum.
 *
 * @see _wave_energy
 * @see localfn_minimum_bracket
 * @private
 */
static mutil_errcode localfn_reduced_log_likelihood(
  const double_mat *Cprime,
  wave_energy      *wave,
  const double_mat *sum_NjlogCprime,
  sint32            delta_index,
  double           *result)
{
  double         innovation_variance = 0.0;
  double         sum_NjlogCprime_correction;
  mutil_errcode  err;
  sint32         j;
  sint32         n_level = wave->n_level;
  sint32         n_level_usable;
  sint32         n_coeff_used = 0;
  sint32         row_global;

  MUTIL_TRACE( "Start localfn_reduced_log_likelihood() ::" );

  /* validate input matrices and check for NULL */

  LOCALDEF_CHECK_NULL_POINTER_FDP( Cprime, double_mat, matdbl );
  LOCALDEF_CHECK_NULL_POINTER_FDP( wave->energy, double_mat, matdbl );
  LOCALDEF_CHECK_NULL_POINTER_FDP( sum_NjlogCprime, double_mat, matdbl );

  /* initialize variables */

  sum_NjlogCprime_correction = 0.0;
  n_level_usable             = wave->n_level_usable;

  /* calculate the correction for sum_NjlogCprime
     (which is obtained by using ALL available
     levels specified by the user in the levels
     vector) by subtracting the sum_NjlogCprime
     over those levels not included in the
     unbiased case  */

  if ( !wave->biased ){

    /* define (global) row to access in Cprime table */

    row_global = delta_index * n_level;

    for ( j = n_level_usable; j < n_level; j++ ){

      sum_NjlogCprime_correction +=
        ( wave->n_coeff )->data[ j ] *
        log( Cprime->data[ row_global + j ] );
    }
  }

  /* calculate the innovation variance */

  err = localfn_innovation_variance(
    Cprime,
    wave,
    delta_index,
    &n_coeff_used,
    &innovation_variance);
  if ( err ) return err;


  if ( n_level_usable >= MINIMUM_USABLE_LEVELS ){


    *result = (double) ( n_coeff_used ) * log( innovation_variance ) +
      sum_NjlogCprime->data[ delta_index ] - sum_NjlogCprime_correction;

    /* add contribution from the scaling coefficients if a
       blocked MLE for a stationary FD process model is
       being used to estimate the FD model parameters */

    if ( wave->stationary ){

      if ( wave->Cprime_scaling > (double) 0.0 ){

        *result += (double) log( wave->Cprime_scaling );
      }
    }
  }
  else *result = 0.0;

  MUTIL_TRACE( ":: Done with localfn_reduced_log_likelihood()" );

  return MUTIL_ERR_OK;
}

/** Extraction of mid-octave SDF values for an FD process
 * from a pre-defined table.
 * Returns a portion of a pre-defined table (Cprime) containing
 * mid-octave average SDF values for a fractionally differenced
 * process.
 *
 * @limits The range of delta (over which the values in the pre-defined
 *         Cprime table are extracted) is limited to that defined by the macros
 *         (DELTA\_MIN, DELTA\_MAX). These macros are defined in wav\_fdp.c.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_Cprime_lookup_table( &levels, &delta_range, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param levels      Pointer to a pre-allocated single-column or single-row
 *                    universal matrix of type MUTIL\_SINT32 containing the
 *                    decomposition levels used in the MODWT.
 * @param delta_range Pointer to a pre-allocated single-column
 *                    or single-row universal matrix of type
 *                    MUTIL\_DOUBLE containing exactly two elements.
 *                    This vector holds the range of delta over which
 *                    the FD estimates are restricted.
 * @param intrp_ptr   Interrupt pointer.
 * @param result      A pointer to a pre-allocated double matrix containing
 *                    the result. The number of rows in this matrix
 *                    can be calculated by the localfn\_delta\_indices()
 *                    function. The number of columns should be equal to
 *                    the length of the levels vector.
 *
 * @see wavuniv_fdp_estimator_instantaneous
 * @private
 */
static mutil_errcode localfn_Cprime_lookup_table(
  const univ_mat *levels,
  const univ_mat *delta_range,
  void           *intrp_ptr,
  double_mat     *result )
{
  double *pd_delta_range;
  mutil_errcode err;
  sint32  col;
  sint32  i;
  sint32  j;
  sint32  n_delta;
  sint32  n_level_max = 0;
  sint32  n_level_min = 0;
  sint32  row;
  sint32  row_end = 0;
  sint32  row_start;
  sint32 *ps_levels;

  /* the following macro is defined in wav_look.h
     and loads the Cprime lookup table. see that header
     file for details */

  MUTIL_TRACE( "Start localfn_Cprime_lookup_table() ::" );

  /*** check levels ... ***/

  /* ... for valid matrix structure and NULL pointer */
  LOCALDEF_CHECK_NULL_POINTER_FDP( levels, univ_mat, matuniv );

  /* ... for type MUTIL_SINT32 */
  if ( levels->type != MUTIL_SINT32 ){
    MUTIL_ERROR( "Levels matrix must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... for vector */
  if ( !MATANY_IS_VEC( &(levels->mat.s32mat) ) ){
    MUTIL_ERROR( "Levels matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /*** check delta_range ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( delta_range, univ_mat, matuniv );

  /* ... for type MUTIL_DOUBLE */

  if ( delta_range->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "delta_range matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... for vector */

  if ( !MATANY_IS_VEC( &(delta_range->mat.dblmat) ) ){
    MUTIL_ERROR( "delta_range matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... for size */

  if ( MATUNIV_NELEM( delta_range ) != 2 ){
    MUTIL_ERROR( "delta_range matrix must have 2 elements" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* validate output */

  err = matdbl_validate( result );
  if ( err ) return err;

  /* create pointers */

  pd_delta_range = delta_range->mat.dblmat.data;
  ps_levels      = levels->mat.s32mat.data;

  /* check to make sure that the user specifed search range for delta is
     in the range of the Cprime table */

  if ( ( pd_delta_range[ 0 ] < FDP_SDF_TABLE_DELTA_MIN ) ||
       ( pd_delta_range[ 1 ] > FDP_SDF_TABLE_DELTA_MAX ) ){

    MUTIL_ERROR( "Specified delta range exceeds table values. " \
                 "Reduce range to -15 <= delta <= 15" );
    return MUTIL_ERR_OUT_OF_BOUNDS;
  }

  /* check to make sure that the user specified decomposition levels do
     not exceed Cprime table  */

  err = mats32_range( &(levels->mat.s32mat),
                      intrp_ptr,
                      &n_level_min,
                      &n_level_max );
  if ( err ) return err;

  if ( ( n_level_min < 1 ) || ( n_level_max > 20 ) ){

    MUTIL_ERROR( "Specified levels exceed table values. " \
                 "Coerce levels such that 1 <= levels <= 20" );
    return MUTIL_ERR_OUT_OF_BOUNDS;
  }

  /* calculate the rows to access in Cprime */

  err =  localfn_delta_indices(
    delta_range,
    &row_start,
    &row_end,
    &n_delta);
  if ( err ) return err;

  /* extract table values and return result */

  for ( i = 0, row = row_start; row <= row_end; row++ ){

    for( j = 0; j < MATUNIV_NELEM( levels ); j++ ){

      col = ps_levels[ j ] - 1;

      result->data[ i++ ] = FDP_SDF_TABLE[ row ][ col ];
    }
  }

  MUTIL_TRACE( ":: Done with localfn_Cprime_lookup_table()" );

  return MUTIL_ERR_OK;
}

/** Validate the inputs to the FD parameter
 * estimation functions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @usage #err = localfn_fdp_estimator_input_check( &time_series, &levels, filter_type, filter_length, biased, estimator, dof, &delta_range, intrp_ptr, FALSE );#
 * @return       Standard mutils error/OK code.
 * @param  time_series    Pointer to a single column universal matrix
 *                        of type MUTIL\_DOUBLE containing the time series
 *                        to analyze.
 * @param  levels         Pointer to a pre-allocated single-column
 *                        universal matrix of type MUTIL\_SINT32.
 *                        This argument contains the decomposition
 *                        levels. The levels can be given in any order, but
 *                        must be positive.
 * @param  filter_type    The type of Daubechies filter.
 * @param  filter_length  The length of the wavelet filter.
 * @param  estimator      The method to estimate the FD model parameters.
 *                        This argument is an enumerated type
 *                        \Ref{_wav_fdp_estimator}.
 * @param  biased         Boolean flag denoting biased or unbiased
 *                        estimates. Biased estimates are those which
 *                        use all available levels in calculating
 *                        the FD model parameters. Unbiased estimates
 *                        are calculated with only those wavelet
 *                        coefficients not subject to circular filter
 *                        operations, i.e. only the interior wavelet
 *                        coefficients are used in calculating unbiased
 *                        estimates.
 * @param  dof            FOR INSTANTANEOUS ESIMTATORS:
 *
 *                        The degree of freedom (dof) order. The number of
 *                        chi-squared dofs used in estimating the FD
 *                        parameters is equal to 2K + 1, where K is
 *                        the dof\_order such that (K > 0). As the order
 *                        increases, the estimates will become smoother
 *                        but less localized in time.
 *
 *                        FOR BLOCK ESTIMATORS:
 *
 *                        The equivalent degrees of freedom (EDOF) mode.
 * @param  delta_range    Pointer to a two element universal matrix of
 *                        type MUTIL\_DOUBLE.  This vector contains the
 *                        minimum and maximum search range
 *                        to use in estimating the FD parameter delta.
 * @param  intrp_ptr      Pointer for implementation of interrupt checking.
 * @param  block          Boolean value indicating whether the caller is
 *                        a blocked estimator or not (instantaneous).
 *
 * @see wavuniv_fdp_estimator_block
 * @see wavuniv_fdp_estimator_instantaneous
 * @private
 */
static mutil_errcode localfn_fdp_estimator_input_check(
  const univ_mat    *time_series,
  const univ_mat    *levels,
  wav_filter_type    filter_type,
  sint32             filter_length,
  boolean            biased,
  wav_fdp_estimator  estimator,
  sint32             dof,
  const univ_mat    *delta_range,
  void              *intrp_ptr,
  boolean            block )
{
  mutil_errcode err;
  sint32        n_level;
  sint32        n_level_boundary;
  sint32        n_level_max = 0;
  sint32        n_level_min = 0;
  sint32        n_sample;
  univ_mat      unique_levels;

  MUTIL_TRACE( "Start localfn_fdp_estimator_input_check()" );

  /*** check time series ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( time_series, univ_mat, matuniv );

  /* ... for type MUTIL_DOUBLE */

  if ( time_series->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Time series must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

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

  /*** check levels ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( levels, univ_mat, matuniv );

  /* ... for type MUTIL_SINT32 */

  if ( levels->type != MUTIL_SINT32 ){
    MUTIL_ERROR( "Levels matrix must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( &(levels->mat.s32mat) ) ){
    MUTIL_ERROR( "Levels matrix must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if number of elements is greater than 1 */

  n_level = MATUNIV_NELEM( levels );

  if ( n_level < 2 ){
    MUTIL_ERROR( "Number of levels must be greater than one." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... to ensure that all levels are positive */

  err = mats32_range( &(levels->mat.s32mat),
                      intrp_ptr,
                      &n_level_min,
                      &n_level_max );
  if ( err ) return err;

  if ( n_level_min <= 0 ){
    MUTIL_ERROR( "All values in the levels matrix must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  /* ... to ensure that there are no redundant levels specified */

  err = matuniv_unique( levels, FALSE, intrp_ptr, &unique_levels );
  if ( err ) return err;

  if ( MATUNIV_NCOL( &unique_levels ) != n_level ){
    MUTIL_FREE_WARN( matuniv, &unique_levels );
    MUTIL_ERROR( "Redundant values in the levels matrix are not allowed." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  MUTIL_FREE_WARN( matuniv, &unique_levels );

  n_sample = MATUNIV_NELEM( time_series );

  if ( !block ){

    /* ... for the unbiased case:
       that the maximum level does not exceed
       that at which there exists at least one
       interior (non-boundary) wavelet coefficient.
       along the way check to ensure that the filter
       length is less than the number of samples in
       the time series. */

    if ( !biased ){

      if ( n_sample > filter_length ){

        n_level_boundary =
          (sint32) floor( MUTIL_LOG2( ( (double) ( n_sample - 1 ) /
          (double) ( filter_length - 1 ) +
          1.0 ) ) );
      }
      else{
        MUTIL_ERROR( "The maximum value in the levels matrix should not " \
          "exceed the maximum level in which there exists at " \
          "least one interior wavelet coefficient." );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }

      if ( n_level_max > n_level_boundary ){
        MUTIL_ERROR( "The maximum value in the levels matrix should not " \
          "exceed the maximum level in which there exists at " \
          "least one interior wavelet coefficient." );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
    }
  }
  else{

    /* in the case we are calculating a block MLE of delta,
       we have to make sure that the maximum level does not
       exceed that of the DWT */

    if ( estimator == WAV_FDP_MAXIMUM_LIKELIHOOD ){

      if ( n_level_max > (sint32)
        floor( log( ( double ) n_sample ) / log( 2.0 ) ) ) {

        MUTIL_ERROR( "The maximum value in the levels matrix should not " \
          "exceed the maximum allowable level in a DWT." );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
    }
  }

  /*** check filter length ***/

  if ( filter_length < 2 ){
    MUTIL_ERROR( "Filter length must be greater than or equal to 2." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( ( filter_length % 2 ) != 0 ){
    MUTIL_ERROR( "Filter length must be even." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*** check filter type ***/

  switch( filter_type ){
  case WAV_FILTER_HAAR:
  case WAV_FILTER_EXTREMAL_PHASE:
  case WAV_FILTER_LEAST_ASYMMETRIC:
  case WAV_FILTER_BEST_LOCALIZED:
  case WAV_FILTER_COIFLET:
    break;
  default:
    MUTIL_ERROR( "Model filter type is invalid" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*** check estimator for correct type ***/

  switch( estimator ) {
  case WAV_FDP_MAXIMUM_LIKELIHOOD:
  case WAV_FDP_LEAST_SQUARES:
    break;
  default:
    MUTIL_ERROR( "FD process parameter estimation scheme " \
                 "currently unsupported" );
    return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }


  /*** check the dof variable ... */

  if ( !block ){

    /* ... for positive value or zero value ***/

    if ( dof < 0 ){
      MUTIL_ERROR( "Chi-squared degree of freedom order must be " \
        "greater than or equal to zero." );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
  }
  else{

    if ( estimator == (wav_fdp_estimator) WAV_FDP_LEAST_SQUARES ){

      if ( dof < 1 || dof > 3 ){
        MUTIL_ERROR( "EDOF mode must be a 1, 2, or 3." );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
    }
  }

  /*** check delta_range ... */

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( delta_range, univ_mat, matuniv );

  /* ... for type MUTIL_DOUBLE */

  if ( delta_range->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Delta range matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( &(delta_range->mat.dblmat) ) ){
    MUTIL_ERROR( "Delta range matrix must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if number of elements is positive */

  if ( MATUNIV_NELEM( delta_range ) != 2 ){
    MUTIL_ERROR( "Number of elements in delta range matrix must equal two." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* .. if the range of delta is supported */

  if (  delta_range->mat.dblmat.data[ 0 ] >=
        delta_range->mat.dblmat.data[ 1 ] ){
    MUTIL_ERROR( "The first value of the delta range vector must be " \
                 "less than the second value." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( ( delta_range->mat.dblmat.data[ 0 ] < FDP_SDF_TABLE_DELTA_MIN ) ||
       ( delta_range->mat.dblmat.data[ 1 ] > FDP_SDF_TABLE_DELTA_MAX ) ){

    MUTIL_ERROR( "Specified delta range exceeds table values. " \
                 "Reduce range to -15 <= delta <= 15" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  MUTIL_TRACE( "End localfn_fdp_estimator_input_check()" );
  return MUTIL_ERR_OK;
}


/***********************************************/
/* STATIC FDP SDF FUNCTION DEFINITIONS         */
/***********************************************/

/** Calculates a function value used in a Taylor series
 * approximation of the mid-octave SDF values for a
 * FD process over a frequency range corresponding to
 * wavelet decomposition levels j = 1,2.
 * The function is u(x) = pow(x,b) / pow(sin(x), b)
 * where
 *
 *     x = pi*f
 *     b = 2*delta
 *
 * The variable f is normalized frequency and delta is the FD parameter.
 *
 * Reference:
 *
 * D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press,
 * 2000, pp. 353.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @usage #err = localfn_taylorseries_u( &x, beta, &result );#
 * @return       Standard mutils error/OK code.
 * @param x      Pointer to universal matrix containing mid-octave
 *               normalized frequencies corresponding to the wavelet
 *               subbands over which the SDF for an FD process is
 *               being approximated.
 * @param beta   A double value obtained by multiplying the
 *               FD parameter (delta) by 2.0.
 * @param result A pointer to a universal matrix containing u(x).
 *
 * @see wavuniv_fdp_bandpass_variance
 * @see localfn_taylor_series_du
 * @see localfn_taylor_series_ddu
 * @see localfn_taylorseries_coeff_weights
 * @see localfn_taylorseries_coeff
 * @see localfn_fdp_bandpass_levels12
 * @private
 */
static mutil_errcode localfn_taylorseries_u(
  const univ_mat *x,
  double          beta,
  univ_mat       *result )
{
  double        *pd_result;
  double        *pd_x;
  mutil_errcode  err;
  sint32         i;

  MUTIL_TRACE( "Start localfn_taylorseries_u() ::" );

  /* begin I/O checking */
  /* validate matrices  */

  err = matuniv_validate( x );
  if ( err ) return err;

  if ( x->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Input matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  err = matuniv_validate( result );
  if ( err ) return err;

  if ( result->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Output matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* create pointer to data */

  pd_x      = x->mat.dblmat.data;
  pd_result = result->mat.dblmat.data;

  /* calculate results */

  for ( i = 0; i < MATUNIV_NELEM( x ); i++ ){
    pd_result[ i ] = MUTIL_POW( pd_x[ i ] / sin( pd_x[ i ] ), beta );
  }

  MUTIL_TRACE( ":: Done with localfn_taylorseries_u()" );

  return MUTIL_ERR_OK;
}

/** Calculates the first derivative of a function u(x)
 * which is used in a Taylor series
 * approximation of the mid-octave SDF values for a
 * FD process over a frequency range corresponding to
 * wavelet decomposition levels j = 1,2.
 * The function u(x) is defined as
 *
 *     u(x) = pow(x,b) / pow(sin(x), b)
 *
 * where
 *
 *     x = pi*f
 *     b = 2*delta
 *
 * The variable f is normalized frequency and delta is the FD parameter.
 *
 * Reference:
 *
 * D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press,
 * 2000, pp. 354.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @usage #err = localfn_taylorseries_du( &x, beta, &result );#
 * @return  Standard mutils error/OK code.
 * @param x      Pointer to universal matrix containing mid-octave
 *               normalized frequencies corresponding to the wavelet
 *               subbands over which the SDF for an FD process is
 *               being approximated.
 * @param beta   A double value obtained by multiplying the
 *               FD parameter (delta) by 2.0.
 * @param result A pointer to a universal matrix containing u'(x).
 *
 * @see wavuniv_fdp_bandpass_variance
 * @see localfn_taylor_series_u
 * @see localfn_taylor_series_ddu
 * @see localfn_taylorseries_coeff_weights
 * @see localfn_taylorseries_coeff
 * @see localfn_fdp_bandpass_levels12
 * @private
 */
static mutil_errcode localfn_taylorseries_du(
  const univ_mat *x,
  double          beta,
  univ_mat       *result )
{
  double         b   = beta;
  double         bm1 = beta - 1;
  double         bp1 = beta + 1;
  double        *pd_result;
  double        *pd_x;
  mutil_errcode  err;
  sint32         i;

  MUTIL_TRACE( "Start localfn_taylorseries_du() ::" );

  /* begin I/O checking */
  /* validate matrices  */

  err = matuniv_validate( x );
  if ( err ) return err;

  if ( x->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Input matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  err = matuniv_validate( result );
  if ( err ) return err;

  if ( result->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Output matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* create pointer to data */

  pd_x      = x->mat.dblmat.data;
  pd_result = result->mat.dblmat.data;

  /* calculate results */

  for ( i = 0; i < MATUNIV_NELEM( x ); i++ ){
    pd_result[ i ] = b * MUTIL_POW( pd_x[ i ], bm1 ) /
      MUTIL_POW( sin( pd_x[ i ]), b ) -
      ( b * MUTIL_POW( pd_x[ i ], b ) * cos( pd_x[ i ] ) ) /
      MUTIL_POW( sin( pd_x[ i ] ), bp1 );
  }

  MUTIL_TRACE( ":: Done with localfn_taylorseries_du()" );

  return MUTIL_ERR_OK;
}

/** Calculates the second derivative of a function u(x)
 * which is used in a Taylor series
 * approximation of the mid-octave SDF values for a
 * FD process over a frequency range corresponding to
 * wavelet decomposition levels j = 1,2.
 * The function u(x) is defined as
 *
 *     u(x) = pow(x,b) / pow(sin(x), b)
 *
 * where
 *
 *     x = pi*f
 *     b = 2*delta
 *
 * The variable f is normalized frequency and delta is the FD parameter.
 *
 * Reference:
 *
 * D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press,
 * 2000, pp. 354.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @usage #err = localfn_taylorseries_ddu( &x, beta, &result );#
 * @return  Standard mutils error/OK code.
 * @param x      Pointer to universal matrix containing mid-octave
 *               normalized frequencies corresponding to the wavelet
 *               subbands over which the SDF for an FD process is
 *               being approximated.
 * @param beta   A double value obtained by multiplying the
 *               FD parameter (delta) by 2.0.
 * @param result A pointer to a universal matrix containing u''(x).
 *
 * @see wavuniv_fdp_bandpass_variance
 * @see localfn_taylor_series_u
 * @see localfn_taylor_series_du
 * @see localfn_taylorseries_coeff_weights
 * @see localfn_taylorseries_coeff
 * @see localfn_fdp_bandpass_levels12
 * @private
 */
static mutil_errcode localfn_taylorseries_ddu(
  const univ_mat *x,
  double          beta,
  univ_mat       *result )
{
  double         b   = beta;
  double         b2  = beta * beta;
  double         bm1 = beta - 1;
  double         bm2 = beta - 2;
  double         bp1 = beta + 1;
  double         bp2 = beta + 2;
  double        *pd_result;
  double        *pd_x;
  mutil_errcode  err;
  sint32         i;

  MUTIL_TRACE( "Start localfn_taylorseries_ddu() ::" );

  /* begin I/O checking */
  /* validate matrices  */

  err = matuniv_validate( x );
  if ( err ) return err;

  if ( x->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Input matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  err = matuniv_validate( result );
  if ( err ) return err;

  if ( result->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Output matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* create pointer to data */

  pd_x      = x->mat.dblmat.data;
  pd_result = result->mat.dblmat.data;

  /* calculate results */

  for ( i = 0; i < MATUNIV_NELEM( x ); i++ ){

    pd_result[ i ] = ( b * MUTIL_POW( pd_x[ i ], b ) +
                       b * bm1 * MUTIL_POW( pd_x[ i ], bm2 ) ) /
      MUTIL_POW( sin( pd_x[ i ] ), b ) -
      ( 2 * b2 * MUTIL_POW( pd_x[ i ], bm1 ) * cos( pd_x[ i ] ) ) /
      MUTIL_POW( sin( pd_x[ i ] ), bp1 ) +
      ( b * bp1 * MUTIL_POW( pd_x[ i ], b ) *
      MUTIL_POW( cos( pd_x[ i ] ), 2) ) /
      MUTIL_POW( sin( pd_x[ i ] ), bp2 );
  }

  MUTIL_TRACE( ":: Done with localfn_taylorseries_ddu()" );

  return MUTIL_ERR_OK;
}

/** Calculates the weights applied to the Taylor series
 * coefficients in approximating the mid-octave SDF values for a
 * FD process over a frequency range corresponding to
 * wavelet decomposition levels j = 1,2.
 *
 * Reference:
 *
 * D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press,
 * 2000, pp. 354.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @usage #err = localfn_taylorseries_coeff_weights( &levels, 0, beta, &result );#
 * @return       Standard mutils error/OK code.
 * @param levels Pointer to universal matrix containing
 *               the indices of the wavelet decomposition
 *               levels over which the mid-octave FD process
 *               SDF is being calculated. Typically this
 *               vector = {1, 2}.
 * @param n      Order of Taylor series coefficient.
 * @param beta   A double value obtained by multiplying the
 *               FD parameter (delta) by 2.0.
 * @param result A pointer to a universal matrix containing the weights.
 *
 * @see wavuniv_fdp_bandpass_variance
 * @see localfn_taylor_series_u
 * @see localfn_taylor_series_du
 * @see localfn_taylorseries_coeff_weights
 * @see localfn_taylorseries_coeff
 * @see localfn_fdp_bandpass_levels12
 *
 * @see wavuniv_fdp_bandpass_variance
 * @private
 */
static mutil_errcode localfn_taylorseries_coeff_weights(
  const univ_mat *levels,
  sint32          n,
  double          beta,
  univ_mat       *result )
{
  sint32        *ps_levels;
  double        *pd_result;
  double         fac = (double) n + 1.0 - beta;
  mutil_errcode  err;
  sint32         i;

  MUTIL_TRACE( "Start localfn_taylorseries_coeff_weights() ::" );

  /* begin I/O checking */
  /* validate matrices  */

  err = matuniv_validate( levels );
  if ( err ) return err;

  if ( levels->type != MUTIL_SINT32  ){
    MUTIL_ERROR( "Input matrix must be of type SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  err = matuniv_validate( result );
  if ( err ) return err;

  if ( result->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Output matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* create pointer to data */

  ps_levels = levels->mat.s32mat.data;
  pd_result = result->mat.dblmat.data;

  /* calculate results */

  for ( i = 0; i < MATUNIV_NELEM( levels ); i++ ){

    if ( MUTIL_ABS( beta - (double) (n+1) ) < MUTIL_DOUBLE_EPSILON )
      pd_result[ i ] = log( 2.0 );
    else pd_result[ i ] =
           MUTIL_POW( MUTIL_PI /
           MUTIL_POW( 2, ps_levels[ i ] + 1 ), fac ) *
           ( MUTIL_POW( 2, fac ) - 1 ) / fac;
  }

  MUTIL_TRACE( ":: Done with localfn_taylorseries_coeff_weights()" );

  return MUTIL_ERR_OK;
}

/** Calculates the Taylor series coefficients used to
 * approximate the mid-octave SDF values for a
 * FD process over a frequency range corresponding to
 * wavelet decomposition levels j = 1,2.
 *
 * Reference:
 *
 * D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press,
 * 2000, pp. 354.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @usage #err = localfn_taylorseries_coeff( &levels, delta, &result );#
 * @return       Standard mutils error/OK code.
 * @param levels Pointer to universal matrix containing
 *               the indices of the wavelet decomposition
 *               levels over which the mid-octave FD process
 *               SDF is being calculated. Typically this
 *               vector = {1, 2}.
 * @param delta  The FD parameter.
 * @param result A pointer to a universal matrix containing
 *               the Taylor series coefficients.
 *
 * @see wavuniv_fdp_bandpass_variance
 * @see localfn_taylorseries_coeff_weights
 * @private
 */
static mutil_errcode localfn_taylorseries_coeff(
  const univ_mat *levels,
  double          delta,
  univ_mat       *result )
{
  double         beta = 2.0 * delta;
  double        *pd_ddu;
  double        *pd_du;
  double        *pd_fj;
  double        *pd_result;
  double        *pd_u;
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         n_level;
  sint32        *ps_levels;
  univ_mat       um_mid_octave_frequency;
  univ_mat       um_taylorseries_ddu;
  univ_mat       um_taylorseries_du;
  univ_mat       um_taylorseries_u;

  MUTIL_TRACE( "Start localfn_taylorseries_coeff() ::" );

  MEMLIST_INIT( list );

  /* obtain the number of levels */

  n_level = MATUNIV_NELEM( levels );

  /* begin I/O checking */
  /* validate matrices  */

  err = matuniv_validate( levels );
  if ( err ) return err;

  if ( levels->type != MUTIL_SINT32  ){
    MUTIL_ERROR( "Input matrix must be of type SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  err = matuniv_validate( result );
  if ( err ) return err;

  if ( MATUNIV_NROW( result ) != (sint32) 3 ){
    MUTIL_ERROR( "Result matrix must have three rows." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if ( MATUNIV_NCOL( result ) != n_level ){
    MUTIL_ERROR( "Result matrix must have the same number of columns as "
      "there are levels." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if ( result->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Output matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* allocate space for vectors used to calculate the
     Taylor series coefficients */

  err = matuniv_malloc_register( &um_mid_octave_frequency, n_level, 1,
                             MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &um_taylorseries_u, n_level, 1,
                             MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &um_taylorseries_du, n_level, 1,
                             MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &um_taylorseries_ddu, n_level, 1,
                             MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create pointer to data */

  ps_levels = levels->mat.s32mat.data;
  pd_result = result->mat.dblmat.data;
  pd_u      = um_taylorseries_u.mat.dblmat.data;
  pd_du     = um_taylorseries_du.mat.dblmat.data;
  pd_ddu    = um_taylorseries_ddu.mat.dblmat.data;
  pd_fj     = um_mid_octave_frequency.mat.dblmat.data;

  /* fill up the mid-octave frequency vector */

  for ( i = 0; i < n_level; i++ ){
    um_mid_octave_frequency.mat.dblmat.data[ i ] = 3.0 * MUTIL_PI /
      (double) MUTIL_POW( 2, ps_levels[ i ] + 2 );
  }

  /* calculate the Taylor series coefficients out to order 2 */

  err = localfn_taylorseries_u(
    &um_mid_octave_frequency,
    beta,
    &um_taylorseries_u );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = localfn_taylorseries_du(
    &um_mid_octave_frequency,
    beta,
    &um_taylorseries_du );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = localfn_taylorseries_ddu(
    &um_mid_octave_frequency,
    beta,
    &um_taylorseries_ddu );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate results */

  for ( i = 0; i < n_level; i++ ){

    /* a0 coefficients */

    pd_result[ i ] = pd_u[ i ] - pd_fj[ i ] * pd_du[ i ] +
      MUTIL_POW( pd_fj[ i ], 2 ) * pd_ddu[ i ] / 2.0;

    /* a1 coefficients */

    pd_result[ i + n_level ] = pd_du[ i ] - pd_fj[ i ] * pd_ddu[ i ];

    /* a2 coefficients */

    pd_result[ i + 2 * n_level ] = pd_ddu[ i ] / 2.0;
  }

  /* clean up and return */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( ":: Done with localfn_taylorseries_coeff()" );

  return MUTIL_ERR_OK;
}

/** Estimates the mid-octave SDF values for a
 * FD process over a frequency range corresponding to
 * wavelet decomposition levels j = 1,2.
 *
 * Reference:
 *
 * D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press,
 * 2000, pp. 354.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @usage #err = localfn_fdp_bandpass_levels12( delta, 1.0, intrp_ptr, &result );#
 * @return                    Standard mutils error/OK code.
 * @param delta               The FD parameter.
 * @param innovation_variance The FD innovation variance.
 * @param intrp_ptr           Pointer for implementation of interrupt checking.
 * @param result              A pointer to a universal matrix containing the result.
 *
 * @see wavuniv_fdp_bandpass_variance
 * @see localfn_fdp_bandpass_levels3J
 * @see localfn_fdp_bandpass_scaling
 * @private
 */
static mutil_errcode localfn_fdp_bandpass_levels12(
  double    delta,
  double    innovation_variance,
  void     *intrp_ptr,
  univ_mat *result )
{
  double        beta = 2.0 * delta;
  double        sum_A;
  memlist       list;
  mutil_errcode err;
  sint32        i;
  sint32        j;
  univ_mat      um_A_weight;
  univ_mat      um_A_weight_all;
  univ_mat      um_levels;
  univ_mat      um_taylor_coefs;

  MUTIL_TRACE( "Start localfn_fdp_bandpass_levels12() ::" );

  MEMLIST_INIT( list );

  /* begin I/O checking */
  /* validate matrices  */

  err = matuniv_validate( result );
  if ( err ) return err;

  if ( result->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Output matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  if( MATUNIV_NELEM(result) != 2){
    MUTIL_ERROR( "Output matrix must be contain 2 elements." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* allocate space */

  err = matuniv_malloc_register( &um_levels, 2, 1,
                             MUTIL_SINT32, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &um_taylor_coefs, 3, 2,
                             MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &um_A_weight, 1, 2,
                             MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &um_A_weight_all, 3, 2,
                             MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* define the levels vector */
  /* j = 1: f in [1/4, 1/2]   */
  /* j = 2: f in [1/8, 1/4]   */

  um_levels.mat.s32mat.data[0] = 1;
  um_levels.mat.s32mat.data[1] = 2;

  /* calcuate the Taylor series coefficients ... */

  /* ... obtain the integral over each octave
     corresponding to the wavelet decomposition levels */

  for ( i = 0; i <= 2; i++ ){

    /* TS coefficient for order i */

    err = localfn_taylorseries_coeff_weights(
      &um_levels, i, beta, &um_A_weight );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* assign the TS of order i to corresponding column
       in the um_A_weight_all matrix                  */

    err = matuniv_assign_submat(&um_A_weight, i, 0, intrp_ptr,
                                &um_A_weight_all);
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* ... calculate the Taylor series weights */

  err = localfn_taylorseries_coeff( &um_levels, delta, &um_taylor_coefs );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* perform an element by element multiplication of the
     um_taylor_coefs and um_A_weight_all matrices. assign
     the result to the former matrix                   */

  err = matuniv_multiply_elem(
    &um_taylor_coefs,
    &um_A_weight_all,
    intrp_ptr,
    &um_taylor_coefs);
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* finally, add the columns of the um_taylor_coefs
     matrix and weight the result                 */

  for ( i = 0; i < 2; i++ ){

    sum_A = 0.0;

    for ( j = 0; j < 3; j++ ){
      sum_A += um_taylor_coefs.mat.dblmat.data[ i + 2 * j ];
    }

    result->mat.dblmat.data[ i ] =
      MUTIL_POW( 2, um_levels.mat.s32mat.data[ i ] + 1 ) *
      innovation_variance * sum_A /
      ( MUTIL_PI * MUTIL_POW( 2, beta ) );
  }

  /* clean up */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( ":: Done with localfn_fdp_bandpass_levels12()" );

  return MUTIL_ERR_OK;
}

/** Estimates the mid-octave SDF values for a
 * FD process over a frequency range corresponding to
 * wavelet decomposition levels j = 3,...,J where
 * J is the number of discrete wavelet transform levels.
 *
 * Reference:
 *
 * D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press,
 * 2000, pp. 354.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @usage #err = localfn_fdp_bandpass_levels3J( delta, 1.0, intrp_ptr, &result );#
 * @return                    Standard mutils error/OK code.
 * @param levels              Pointer to universal matrix containing
 *                            the indices of the wavelet decomposition
 *                            levels over which the mid-octave FD process
 *                            SDF is being calculated. The values in this vector
 *                            must be greater than 2.
 * @param delta               The FD parameter.
 * @param innovation_variance The FD innovation variance.
 * @param result              A pointer to a universal matrix containing the result.
 *
 * @see wavuniv_fdp_bandpass_variance
 * @see localfn_fdp_bandpass_levels12
 * @see localfn_fdp_bandpass_scaling
 * @private
 */
static mutil_errcode localfn_fdp_bandpass_levels3J(
  const univ_mat *levels,
  double          delta,
  double          innovation_variance,
  univ_mat       *result )
{
  double         beta = 2.0 * delta;
  double        *pd_result;
  mutil_errcode  err;
  sint32         i;
  sint32        *ps_levels;

  MUTIL_TRACE( "Start localfn_fdp_bandpass_levels3J() ::" );

  /* validate matrices  */

  err = matuniv_validate( levels );
  if ( err ) return err;

  if ( levels->type != MUTIL_SINT32  ){
    MUTIL_ERROR( "Input matrix must be of type SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  err = matuniv_validate( result );
  if ( err ) return err;

  if ( result->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Output matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( MATUNIV_NELEM( result ) != MATUNIV_NELEM( levels ) ){
    MUTIL_ERROR( "Result and levels matrix must have same "
                 "number of elements." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* create pointers */

  ps_levels = levels->mat.s32mat.data;
  pd_result = result->mat.dblmat.data;

  /* calculate bandpass approximation of mid-octave SDF
     values under an FDP model and assuming a small
     angle approximation                             */

  if ( MUTIL_ABS( delta - 0.5 ) < MUTIL_DOUBLE_EPSILON ){

    for ( i = 0; i < MATUNIV_NELEM(levels); i++ ){

      pd_result[ i ] = innovation_variance *
        MUTIL_POW( 2.0, (double) ps_levels[ i ] + 0.5 ) *
        log(2.0) / MUTIL_PI / sqrt(2.0);
    }
  }
  else{
    for ( i = 0; i < MATUNIV_NELEM( levels ); i++ ){

      pd_result[ i ] =
        MUTIL_POW( 2.0, ( (double) ps_levels[ i ] + 0.5 ) * beta ) *
        MUTIL_POW( 2.0, -3.0 * delta ) * MUTIL_POW( MUTIL_PI, -beta ) *
        innovation_variance * ( 2.0 - MUTIL_POW( 2.0, beta ) ) /
        ( 1.0 - beta );
    }
  }

  MUTIL_TRACE( ":: Done with localfn_fdp_bandpass_levels3J()" );

  return MUTIL_ERR_OK;
}

/** Estimates the mid-octave SDF values for a
 * FD process over a frequency range corresponding to
 * scaling coefficients of a discrete wavelet transform.
 *
 * Reference:
 *
 * D.B. Percival and  A.T. Walden,
 * ``Wavelet Methods for Time Series Analysis'',
 * Cambridge University Press,
 * 2000, pp. 354.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @usage #err = localfn_fdp_bandpass_scaling( delta, 1.0, intrp_ptr, &result );#
 * @return                    Standard mutils error/OK code.
 * @param Cprime              Pointer to universal matrix containing
 *                            the mid-octave SDF estimates corresponding
 *                            to wavelet decomposition levels 1,...,J
 *                            where J is the number of wavelet transform
 *                            decomposition levels. The number of elements
 *                            must be equal to J+1. The result is returned
 *                            in this vector in the (J+1)th element.
 * @param delta               The FD parameter.
 * @param innovation_variance The FD innovation variance.
 * @param n_sample            The number of samples in the time series.
 *
 * @see wavuniv_fdp_bandpass_variance
 * @see localfn_fdp_bandpass_levels12
 * @see localfn_fdp_bandpass_levels3J
 * @private
 */
static mutil_errcode localfn_fdp_bandpass_scaling(
  univ_mat *Cprime,
  double    delta,
  double    innovation_variance,
  sint32    n_sample )
{
  double        sum = 0.0;
  mutil_errcode err;
  sint32        i;
  sint32        n_level;

  MUTIL_TRACE( "Start localfn_fdp_bandpass_scaling() ::" );

  /* validate matrices  */

  err = matuniv_validate( Cprime );
  if ( err ) return err;

  if ( Cprime->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Input matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* obtain the number of decomposition levels. this value be one less
     than the number of elements in the Cprime vector since the
     last element of Cprime will be used to store the mid-octave
     SDF of the octave correpsonding to the scaling coefficients */

  n_level = MATUNIV_NELEM( Cprime ) - 1;

  /* estimate mid-octave SDF value for scaling coefficients band */

  if ( delta < 0.5 ){

    for ( i = 0; i < n_level; i++ ){

      sum += (double) n_sample /
        MUTIL_POW( 2.0, (double) ( i + 1 ) ) *
        Cprime->mat.dblmat.data[ i ];
    }

    Cprime->mat.dblmat.data[ n_level ] =
      innovation_variance *
      ( (double) n_sample *
        MUTIL_GAMMA( 1.0 - 2.0 * delta ) /
        MUTIL_POW( MUTIL_GAMMA( 1.0 - delta ), 2.0 ) - sum );

    MUTIL_TRACE( ":: Done with localfn_fdp_bandpass_scaling()" );

    return MUTIL_ERR_OK;
  }
  else return MUTIL_ERR_ILLEGAL_VALUE;
}

/** Argument check for the wavuniv\_fdp\_bandpass\_variance() function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_fdp_bandpass_arg_check( &levels, n_sample, &result);#
 * @return Standard mutils error/OK code.
 * @param  levels            Pointer to a pre-allocated single-column or
 *                           single-row universal matrix of type MUTIL\_SINT32.
 *                           This argument contains the decomposition
 *                           levels. The levels can be given in any order,
 *                           but must be positive. However, for n\_sample > 0,
 *                           the levels vector must contain the values
 *                           1, 2, 3, ..., J where J is the maximum
 *                           wavelet transform decomposition level.
 * @param  n_sample          The number of samples in the time series.
 * @param  result            Pointer to pre-allocated single-column universal
 *                           matrix of type MUTIL\_DOUBLE.  In the case where
 *                           n\_sample is negative or equal to zero, the
 *                           result matrix must have the same number of
 *                           elements as does the levels matrix.  For the
 *                           case where n\_sample is positive, the length of
 *                           the results vector must be equal to J+1 where
 *                           J is the wavelet transform decomposition level.
 *
 * @see wavuniv_fdp_bandpass_variance
 * @private
 */
static mutil_errcode localfn_fdp_bandpass_arg_check(
  const univ_mat *levels,
  sint32          n_sample,
  univ_mat       *result )
{
  mutil_errcode err;
  sint32        n_level;

  err = matuniv_validate( levels );
  if ( err ) return err;

  if ( levels->type != MUTIL_SINT32  ){
    MUTIL_ERROR( "Input matrix must be of type SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  if ( MATUNIV_NCOL( levels ) != 1 ){
    MUTIL_ERROR( "Levels matrix must be a one column vector." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  err = matuniv_validate( result );
  if ( err ) return err;

  if ( result->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Output matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  if ( ( MATUNIV_NCOL( result ) > 1 ) & ( MATUNIV_NROW( result ) > 1 ) ){
    MUTIL_ERROR( "Result matrix must be a single coulmn or row vector." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* make sure that the result vector is the correct size */

  n_level = MATUNIV_NELEM( levels );

  if ( n_sample > 0 ){

    if ( MATUNIV_NELEM( result ) != n_level + 1 ){
      MUTIL_ERROR( "Result must have one more element than the "
                   "levels matrix so the result for the scaling "
                   "coefficient band can also be returned." );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

  }
  else{
    if ( MATUNIV_NELEM( result ) != n_level ){
      MUTIL_ERROR( "Result and levels matrices must have same "
                   "number of elements." );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
  }

  return MUTIL_ERR_OK;
}
