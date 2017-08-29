
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_filt.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */

/*
 *****************************************************************
 *                                                               *
 * This file contains routines implementing non-linear filtering *
 * techniques for time series analysis.                          *
 *                                                               *
 *****************************************************************
*/


#include "fra_filt.h"

#include "fra_dim.h"
#include "fra_neig.h"
#include "fra_type.h"

#include "mat_num.h"
#include "mat_stat.h"
#include "mat_type.h"
#include "mat_univ.h"
#include "mat_umat.h"

#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_math.h"
#include "ut_mem.h"


/* Static macro definitions */
#define FRA_DENOISE_MINIMUM_DIM   3


#undef LOCALDEF_CHECK_NULL_POINTER_NOISE
#define LOCALDEF_CHECK_NULL_POINTER_NOISE( DATA_PTR,        \
  DATA_TYPE, TYPE_PREFIX )                                  \
  err = TYPE_PREFIX ## _validate( DATA_PTR );               \
  if ( err ) return err;                                    \
  if ( DATA_PTR == ( DATA_TYPE * ) NULL ) {                 \
  MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
  return MUTIL_ERR_NULL_POINTER;                            \
  }

#undef LOCALDEF_RESIZE_VECTOR
#define LOCALDEF_RESIZE_VECTOR( MATRIX, NEW_LENGTH )        \
  switch( (MATRIX).type ){				    \
    case MUTIL_UINT8:					    \
      (MATRIX).mat.u8mat.nelem = NEW_LENGTH;	            \
      (MATRIX).mat.u8mat.nrow  = NEW_LENGTH;	            \
      (MATRIX).mat.u8mat.ncol  = 1;		            \
      break;						    \
    case MUTIL_SINT8:					    \
      (MATRIX).mat.s8mat.nelem = NEW_LENGTH;	            \
      (MATRIX).mat.s8mat.nrow  = NEW_LENGTH;	            \
      (MATRIX).mat.s8mat.ncol  = 1;		            \
      break;						    \
    case MUTIL_UINT16:					    \
      (MATRIX).mat.u16mat.nelem = NEW_LENGTH;	            \
      (MATRIX).mat.u16mat.nrow  = NEW_LENGTH;	            \
      (MATRIX).mat.u16mat.ncol  = 1;		            \
      break;						    \
    case MUTIL_SINT16:					    \
      (MATRIX).mat.s16mat.nelem = NEW_LENGTH;	            \
      (MATRIX).mat.s16mat.nrow  = NEW_LENGTH;	            \
      (MATRIX).mat.s16mat.ncol  = 1;		            \
      break;						    \
    case MUTIL_UINT32:					    \
      (MATRIX).mat.u32mat.nelem = NEW_LENGTH;	            \
      (MATRIX).mat.u32mat.nrow  = NEW_LENGTH;	            \
      (MATRIX).mat.u32mat.ncol  = 1;		            \
      break;						    \
    case MUTIL_SINT32:					    \
      (MATRIX).mat.s32mat.nelem = NEW_LENGTH;	            \
      (MATRIX).mat.s32mat.nrow  = NEW_LENGTH;	            \
      (MATRIX).mat.s32mat.ncol  = 1;		            \
      break;						    \
    case MUTIL_FLOAT:					    \
      (MATRIX).mat.fltmat.nelem = NEW_LENGTH;	            \
      (MATRIX).mat.fltmat.nrow  = NEW_LENGTH;	            \
      (MATRIX).mat.fltmat.ncol  = 1;		            \
      break;						    \
    case MUTIL_DOUBLE:					    \
      (MATRIX).mat.dblmat.nelem = NEW_LENGTH;	            \
      (MATRIX).mat.dblmat.nrow  = NEW_LENGTH;	            \
      (MATRIX).mat.dblmat.ncol  = 1;		            \
      break;						    \
    default:                                                \
      /* Will never get here due to prior error checking */ \
      break;                                                \
  }

#undef LOCALDEF_INCREMENT_UNIVMAT_DATA_POINTER
#define LOCALDEF_INCREMENT_UNIVMAT_DATA_POINTER( MATRIX,    \
  INCREMENT )                                               \
  switch( (MATRIX).type ){				    \
    case MUTIL_UINT8:					    \
      (MATRIX).mat.u8mat.data  += INCREMENT;	            \
      break;						    \
    case MUTIL_SINT8:					    \
      (MATRIX).mat.s8mat.data += INCREMENT; 	            \
      break;						    \
    case MUTIL_UINT16:					    \
      (MATRIX).mat.u16mat.data += INCREMENT;	            \
      break;						    \
    case MUTIL_SINT16:					    \
      (MATRIX).mat.s16mat.data += INCREMENT;	            \
      break;						    \
    case MUTIL_UINT32:					    \
      (MATRIX).mat.u32mat.data += INCREMENT;	            \
      break;						    \
    case MUTIL_SINT32:					    \
      (MATRIX).mat.s32mat.data += INCREMENT;	            \
      break;						    \
    case MUTIL_FLOAT:					    \
      (MATRIX).mat.fltmat.data += INCREMENT;	            \
      break;						    \
    case MUTIL_DOUBLE:					    \
       (MATRIX).mat.dblmat.data += INCREMENT;	            \
      break;						    \
    default:                                                \
      /* Will never get here due to prior error checking */ \
      break;                                                \
  }

/* Median filtering of a time series. */
/*                                    */
/* Documented in fra_filt.h           */
/* Written by William Constantine     */

mutil_errcode frauniv_filter_median(
  const univ_mat *time_series,
  const sint32    order,
  void           *intrp_ptr,
  univ_mat       *result )
{
  memlist        list;
  mutil_errcode  err;
  sint32         k;
  sint32         ilo;
  sint32         ihi;
  sint32         offset;
  sint32         dhi;
  sint32         n_sample;
  boolean        is_even = (boolean) ( ( order % 2 ) == 0 );
  univ_mat       ts;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_filter_median()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /*** check input data ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_NOISE( time_series, univ_mat, matuniv );

  /* ... to make sure it is not a complex matrix */

  if ( time_series->type == MUTIL_DCOMPLEX ){
    MUTIL_ERROR( "Time series must not be of type dcomplex." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... to see if it is a vector */

  if ( MATUNIV_NROW( time_series ) > 1 && MATUNIV_NCOL( time_series ) > 1 ){
    MUTIL_ERROR( "Time_Series matrix must be a single-row or single-column vector." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  n_sample = MATUNIV_NELEM( time_series );

  /* check order argument */

  if ( order < 1 ){
    MUTIL_ERROR( "Order must be a positive integer." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( order > 2 * n_sample ){
    MUTIL_ERROR( "Order must be less than or equal to 2 * N where N is the number "
      " of samples in the original time series." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* allocate space and register with memory manager */

  err = matuniv_malloc_register( result, MATUNIV_NROW( time_series ),
    MATUNIV_NCOL( time_series ), MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* use integer division to form half-window length */

  offset = order / 2;

  dhi = ( is_even ? offset - 1 : offset );

  /* copy the universal matrix header for the time series
     into another univ_mat header so that we can mess
     with the nelem, nrow, and ncol elements */

  err = matuniv_wrap_univ_matrix( &ts, (univ_mat *) time_series );
  if ( err ) return( err );

  /* set time series matrix components */

  LOCALDEF_RESIZE_VECTOR( ts, order );

  for ( k = 0; k < n_sample; k++ ){

    /* set index range of current window */

    ilo = k - offset;
    ihi = k + dhi;

    /* adjust for boundaries */

    if ( ilo < 0 ){
      ilo = 0;
    }

    if ( ihi > n_sample - 1 ){
      ihi = n_sample - 1;
    }

    /* trick the matuniv_median() function
       by resetting the time series length
       parameters via the universal matrix
       components */

    LOCALDEF_RESIZE_VECTOR( ts, ihi - ilo + 1 );

    err = matuniv_median(
      &ts,
      intrp_ptr,
      &( result->mat.dblmat.data[ k ] ) );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /*
      update pointer of faux time series.
      don't increment pointer until left hand side
      of window has reached the zero index of the
      time series
    */

    if ( k - offset >= 0 ){

      LOCALDEF_INCREMENT_UNIVMAT_DATA_POINTER( ts, 1 );
    }

    /* check for interrupts */

    if ( MUTIL_INTERRUPT( 3.0 * order, intrp_ptr ) ) {
	 MEMLIST_FREE_ON_ERROR( err, &list );
         MUTIL_ERROR( "user interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }
  }

  /* free nodes corresponding to registered
     memory for the result, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  MUTIL_TRACE( "Done with frauniv_filter_median()" );

  return MUTIL_ERR_OK;
}



/* Denoising of a time series using a local projections. */
/*                                                       */
/* Documented in fra_filt.h                              */
/* Written by Keith L. Davidson                          */
mutil_errcode frauniv_filter_nonlinear_local_projection(
  const univ_mat            *time_series,
  const sint32               dim,
  const sint32               delay,
  const sint32               min_nneig,
  const double               max_radius,
  const fra_distance_metric  distance_metric,
  const sint32               noise_dim,
  const boolean              correct_curvature,
  void                      *intrp_ptr,
  univ_mat                  *result )
{
  register double   *cm_ptr;
  register double   *dbl_ptr;
  register double   *eig_ptr;
  register double   *ts_ptr;
  register double    eig_tmp;
  register double    sum;

  register sint32   *s32_ptr;
  register sint32    rdelay = delay;
  register sint32    rdim   = dim;
  register sint32    idx;
  register sint32    nneig;
  register sint32    d;
  register sint32    i;
  register sint32    j;
  register sint32    n;
  register sint32    skip;

  univ_mat           cov;
  univ_mat           eigval;
  univ_mat           eigvec;
  univ_mat           embedding;
  univ_mat           orig_idx;
  univ_mat           neig_idx;
  univ_mat           neig_dist;
  univ_mat           point;
  double_mat         cent_mass;
  double_mat         temp_cm;
  sint32_mat         neig_found;
  sint32_mat         neighbors;

  double            *res_ptr;
  double             num_ops = (double) 0;
  sint32             nembed;
  sint32             neig_cntr = (sint32) 0;
  sint32             neig_count;
  sint32             run;
  sint32             ts_length;

  boolean            first_run;
  mutil_kdtree       kdtree;
  memlist            mlist;
  mutil_errcode      err;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start in frauniv_filter_nonlinear_local_projection()" );

  /* Initialize memory management list */
  MEMLIST_INIT( mlist );


  /* Check inputs for errors */

  /* dim */
  if ( dim < FRA_DENOISE_MINIMUM_DIM ) {
    MUTIL_ERROR( "Embedding dimension too small" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* delay */
  if ( delay < (sint32) 1 ) {
    MUTIL_ERROR( "Embedding delay must be at least 1" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* time_series (if OK, create an embedding and add that to the memory
     manager) */
  err = matuniv_validate( time_series );
  if ( err ) {
    MUTIL_ERROR( "Input time series is an invalid matrix" );
    return err;
  }
  if ( time_series->type != MUTIL_DOUBLE ) {
    MUTIL_ERROR( "Input time series matrix has illegal type" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  ts_length = MUTIL_MAX( MATUNIV_NROW( time_series ),
    MATUNIV_NCOL( time_series ) );
  if ( ts_length !=  MATUNIV_NELEM( time_series ) ) {
    MUTIL_ERROR( "Input time series must be a row or column vector" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  err = frauniv_embed( time_series, dim, delay, intrp_ptr, &embedding );
  if ( err ) {
    MUTIL_ERROR( "Could not create embedding for given time series" );
    return err;
  }
  err = memlist_member_register( &mlist, (void*) &embedding, MEMTYPE_MATUNIV );
  if ( err ) {
    MUTIL_FREE_WARN( matuniv, &embedding );
    return err;
  }
  nembed = MATUNIV_NROW( &embedding );

  /* min_nneig */
  if ( min_nneig <= dim ) {
    MUTIL_ERROR( "Minimum number of neighbors must be greater "
      "than the embedding dimension" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  if ( min_nneig > ( nembed - 1 ) ) {
    MUTIL_ERROR( "Input min_nneig too large, not enough points in the delay "
      "embedding" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* max_radius */
  if ( max_radius <= (double) 0 ) {
    MUTIL_ERROR( "Maximum search radius must be a positive number" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* distance_metric */
  switch ( distance_metric ) {

  case FRA_DISTANCE_L1:
  case FRA_DISTANCE_L2:
  case FRA_DISTANCE_LINFINITY:
    break;

  default:
    MUTIL_ERROR( "Input distance_metric is invalid" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* noise_dim */
  if ( noise_dim < (sint32) 1 ) {
    MUTIL_ERROR( "Noise dimension must be at least 1" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  if ( noise_dim >= dim ) {
    MUTIL_ERROR( "Noise dimension must be less than embedding dimension" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }


  /* Create a kd-tree from the embedding for subsequent neighbor searching.
     Then remove the embedding matrix from the memory management list */
  err = mutil_kdtree_malloc_register( &kdtree, &embedding, (sint32) 1, &mlist );
  if ( err ) {
    MUTIL_FREE_WARN( memlist, &mlist );
    MUTIL_ERROR( "Could not create kd-tree from embedding matrix" );
    return err;
  }
  num_ops += (double) nembed * log( (double) nembed );


  /* Allocate memory to hold the neighbor indices, the number of neighbors
     for each embedding point, the center-of-mass vectors, and, if needed,
     temporary storage used for center-of-mass vector correction */
  err = mats32_malloc_register( &neig_found, nembed, (sint32) 1, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );
  err = mats32_malloc_register( &neighbors, (sint32) 1, (sint32) 1, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );
  err = matdbl_malloc_register( &cent_mass, nembed, dim, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );
  if ( correct_curvature ) {
    err = matdbl_malloc_register( &temp_cm, nembed, dim, &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );
  }


  /* Initializations */
  neig_count             = (sint32) 0;
  first_run              = TRUE;
  point.type             = MUTIL_DOUBLE;
  point.mat.dblmat.nrow  = (sint32) 1;
  point.mat.dblmat.ncol  = (sint32) dim;
  point.mat.dblmat.nelem = (sint32) dim;
  point.mat.dblmat.data  = embedding.mat.dblmat.data;


  /* Run through the embedding twice. In the first run we search for neighbors
     and compute the center-of-mass vectors. In the second run we compute the
     corrected center-of-mass vectors, if the user desires, and the final noise
     component vectors to be subtracted from the embedding. */
  for ( run = 0; run < (sint32) 2; run++, first_run = FALSE ) {

    /* Loop through the embedding points */
    for ( n = 0; n < nembed; n++ ) {

      /* If we're in the first run through the embedding points then search for
         neighbors and record the neighbor indices */
      if ( first_run ) {

        /* set dummy matrix to the n-th point of the embedding matrix */

        /* search for neighbors */
        err = frauniv_neighbor_find_arbitrary(
          &point,
          &kdtree,
          (sint32) 0,
          max_radius,
          distance_metric,
          (univ_mat*) NULL,
          FALSE,
          (sint32) 0,
          (univ_mat*) NULL,
          (void*) NULL,
          &orig_idx,
          &neig_idx,
          &neig_dist );
        if ( err && ( err != MUTIL_ERR_ZERO_NEIGHBORS_FOUND ) ) {
          MUTIL_FREE_WARN( memlist, &mlist );
          return err;
        }

        /* free unneeded memory */
        if ( err != MUTIL_ERR_ZERO_NEIGHBORS_FOUND ) {
          MUTIL_FREE_WARN( matuniv, &orig_idx );
          MUTIL_FREE_WARN( matuniv, &neig_dist );
        }

        /* check for interrupts */
        num_ops += MUTIL_MIN( 0, MATUNIV_NELEM( &neig_idx ) );
        if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
          if ( err != MUTIL_ERR_ZERO_NEIGHBORS_FOUND )
            MUTIL_FREE_WARN( matuniv, &neig_idx );
          MEMLIST_FREE_ON_ERROR( err, &mlist );
          MUTIL_ERROR( "user interrupt" );
          return MUTIL_ERR_INTERRUPT;
        }

        /* if an insufficient number of neighbors were found using the
           distance criteria then grab the minimum number of neighbors
           requested by the user. */
        if ( neig_idx.mat.s32mat.nelem < min_nneig ||
          err == MUTIL_ERR_ZERO_NEIGHBORS_FOUND ) {

          /* free the neighbor indices from previous search */
          if ( err != MUTIL_ERR_ZERO_NEIGHBORS_FOUND )
            MUTIL_FREE_WARN( matuniv, &neig_idx );

          /* find min_nneig neigbors */
          err = frauniv_neighbor_find_arbitrary(
            &point,
            &kdtree,
            min_nneig,
            (double) 0,
            distance_metric,
            (univ_mat*) NULL,
            FALSE,
            (sint32) 0,
            (univ_mat*) NULL,
            (void*) NULL,
            &orig_idx,
            &neig_idx,
            &neig_dist );
          if ( err ) {
            MUTIL_FREE_WARN( memlist, &mlist );
            return err;
          }

          /* free unneeded memory */
          MUTIL_FREE_WARN( matuniv, &orig_idx );
          MUTIL_FREE_WARN( matuniv, &neig_dist );

          /* check for interrupts */
          num_ops += min_nneig;
          if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
            if ( err != MUTIL_ERR_ZERO_NEIGHBORS_FOUND )
              MUTIL_FREE_WARN( matuniv, &neig_idx );
            MEMLIST_FREE_ON_ERROR( err, &mlist );
            MUTIL_ERROR( "user interrupt" );
            return MUTIL_ERR_INTERRUPT;
          }
        } /* if ( neig_idx.mat.s32mat.nelem < min_nneig || ... ) */

        /* update the neighbor count */
        neig_count += neig_idx.mat.s32mat.nelem;

        /* record the number of neighbors found */
        neig_found.data[ n ] = neig_idx.mat.s32mat.nelem;

        /* reallocate memory for the neighbor indices */
        err = mats32_realloc_register(
          &neighbors,
          neig_count,
          (sint32) 1,
          &mlist );
        if ( err ) {
          MUTIL_FREE_WARN( matuniv, &neig_idx );
          MUTIL_FREE_WARN( memlist, &mlist );
          return err;
        }

        /* initialize the neighbor indices pointer */
        s32_ptr = neighbors.data + neig_count - neig_idx.mat.s32mat.nelem;
        cm_ptr  = cent_mass.data + n * dim;
        dbl_ptr = temp_cm.data + n * dim;
        ts_ptr  = embedding.mat.dblmat.data;

        /* record the neighbor indices and calculate the n-th center-of-mass
           vector */
        for ( d = 0; d < rdim; d++, ts_ptr++ ) {

          sum    = (double) 0;

          for ( i = 0; i < neig_idx.mat.s32mat.nelem; i++ ) {

            if ( d == (sint32) 0 )
              s32_ptr[ i ] = neig_idx.mat.s32mat.data[ i ];

            sum += *( ts_ptr + s32_ptr[ i ] * rdim );
          }

          if ( correct_curvature ) {

            *cm_ptr = *dbl_ptr =
              sum / (double) neig_idx.mat.s32mat.nelem;

            dbl_ptr++;

          } else {

            *cm_ptr = sum / (double) neig_idx.mat.s32mat.nelem;
          }

          cm_ptr++;

        } /* for ( d = 0, idx = n * rdim; ... ) */

        /* free memory */
        MUTIL_FREE_WARN( matuniv, &neig_idx );

        /* check for interrupts */
        num_ops += 2.0 * dim * ( 1.0 + MATUNIV_NELEM( &neig_idx ) );
        if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
          MEMLIST_FREE_ON_ERROR( err, &mlist );
          MUTIL_ERROR( "user interrupt" );
          return MUTIL_ERR_INTERRUPT;
        }

        /* update the pointer into the embedding matrix */
        point.mat.dblmat.data += dim;

      } else { /* if ( first_run ) */

        /* if we're at the first embedding point in the second run through
           the embedding points then we need to do a few things ... */
        if ( n == (sint32) 0 ) {

          /* remove the embedding matrix and kd-tree from the memory
             management list, won't need to search for neighbors during
             the second run through the embedding points */
          err = memlist_member_free( (void*) &kdtree, &mlist );
          MEMLIST_FREE_ON_ERROR( err, &mlist );
          err = memlist_member_free( (void*) &embedding, &mlist );
          MEMLIST_FREE_ON_ERROR( err, &mlist );

          /* Compute the z vectors, z_n = (s_n-\bar{s}), corresponding to
             each center-of-mass vector, first correcting the center-of-mass
             vector for curvature effects if the user desires */
          if ( correct_curvature ) {

            s32_ptr = neighbors.data;

            for ( i = 0; i < nembed; i++ ) {

              /* grab the number of neighbors for the nn-th embedding point */
              nneig = neig_found.data[ i ];

              /* initialize the center-of-mass, time series, and weight vector
                 data pointers */
              dbl_ptr = cent_mass.data + i * rdim;
              ts_ptr  = time_series->mat.dblmat.data + i;

              for ( d = 0; d < rdim; d++, dbl_ptr++, ts_ptr += rdelay ) {

                sum = (double) 0;

                for ( j = 0; j < nneig; j++, s32_ptr++ )
                  sum += temp_cm.data[ *s32_ptr * dim + d ];

                (*dbl_ptr) *= 2.0;
                (*dbl_ptr) -= sum / (double) nneig;

                *dbl_ptr = ( *ts_ptr - *dbl_ptr );

                if ( d != (rdim - 1) ) s32_ptr -= nneig;
              }

              /* check for interrupts */
              num_ops += 4.0 + 5.0 * dim * ( nneig + 1 );
              if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
                MEMLIST_FREE_ON_ERROR( err, &mlist );
                MUTIL_ERROR( "user interrupt" );
                return MUTIL_ERR_INTERRUPT;
              }

            } /* for ( i = 0, idx = 0; ... ) */

          } else { /* if ( correct_curvature ) */

            for ( i = 0; i < nembed; i++ ) {

              /* initialize the center-of-mass, time series, and weight vector
                 data pointers */
              dbl_ptr = cent_mass.data + i * rdim;
              ts_ptr  = time_series->mat.dblmat.data + i;

              for ( d = 0; d < rdim; d++, dbl_ptr++, ts_ptr += rdelay )
                *dbl_ptr = ( *ts_ptr - *dbl_ptr );

            } /* for ( i = 0; ... ) */

            /* check for interrupts */
            num_ops += 4.0 * nembed * (dim + 1);
            if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
              MEMLIST_FREE_ON_ERROR( err, &mlist );
              MUTIL_ERROR( "user interrupt" );
              return MUTIL_ERR_INTERRUPT;
            }

          } /* if ( correct_curvature ) */

          /* allocate memory for the covariance matrix, eigenvalues, and
             eignenvectors. */
          err = matuniv_malloc_register( &cov, dim, dim, MUTIL_DOUBLE, &mlist );
          MEMLIST_FREE_ON_ERROR( err, &mlist );

          err = matuniv_malloc_register(
            &eigvec,
            dim,
            dim,
            MUTIL_DOUBLE,
            &mlist );
          MEMLIST_FREE_ON_ERROR( err, &mlist );

          err = matuniv_malloc_register(
            &eigval,
            dim,
            (sint32) 1,
            MUTIL_DOUBLE,
            &mlist );
          MEMLIST_FREE_ON_ERROR( err, &mlist );

          /* allocate memory for the correction vectors */
          if ( !correct_curvature ) {
            err = matdbl_malloc_register(
              &temp_cm,
              nembed,
              dim,
              &mlist );
            MEMLIST_FREE_ON_ERROR( err, &mlist );
          }
        } /* if ( n == (sint32) 0 ) */

        /* grab the number of neighbors for the n-th embedding point */
        nneig = neig_found.data[ n ];

        /* compute the covariance matrix corresponding to the n-th
           embedding point */
        for ( i = 0; i < rdim; i++ ) {

          idx = i * rdim;

          for ( j = i; j < rdim; j++ ) {

            sum = (double) 0;

            for ( d = 0; d < nneig; d++ ) {

              dbl_ptr = cent_mass.data +
                neighbors.data[ neig_cntr + d ] * rdim;

              sum += dbl_ptr[ i ] * dbl_ptr[ j ];
            }

            cov.mat.dblmat.data[ idx + j ] =
              cov.mat.dblmat.data[ j * dim + i ] =
                sum / (double) nneig;

          } /* for ( j = i; ... ) */

        } /* for ( i = 0; ... ) */

        /* check for interrupts */
        num_ops += (double) dim * ( dim * ( 6.0 * nneig + 4.0 ) + 1.0 );
        if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
          MEMLIST_FREE_ON_ERROR( err, &mlist );
          MUTIL_ERROR( "user interrupt" );
          return MUTIL_ERR_INTERRUPT;
        }

        /* update the neighbor counter */
        neig_cntr += nneig;

        /* compute the eignvector of the covariance matrix */
        err = matuniv_eigen_jacobi(
          &cov,
          intrp_ptr,
          &eigval,
          &eigvec );
        MEMLIST_FREE_ON_ERROR( err, &mlist );

        /* sort eigenvalues in descending order and permute eigenvectors
        accordingly */
        err = matuniv_eigen_sort(
          &eigval,
          &eigvec,
          intrp_ptr,
          &eigval,
          &eigvec );
        MEMLIST_FREE_ON_ERROR( err, &mlist );

        /* check for interrupts */
        num_ops += 10.0 * dim * dim * dim;
        if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
          MEMLIST_FREE_ON_ERROR( err, &mlist );
          MUTIL_ERROR( "user interrupt" );
          return MUTIL_ERR_INTERRUPT;
        }

        /* initialize data pointers */
        cm_ptr  = cent_mass.data + n * rdim;
        dbl_ptr = temp_cm.data + n * rdim;

        /* calculate the correction vector, \delta-s */
        for ( i = 0; i < rdim; i++ ) {

          idx = i * dim;
          sum = (double) 0;

          for ( d = rdim - noise_dim; d < rdim; d++ ) {

            eig_ptr = eigvec.mat.dblmat.data + d;
            eig_tmp = *( eig_ptr + idx );

            for ( j = 0; j < rdim; j++ )
              sum += cm_ptr[ j ] * eig_tmp * *( eig_ptr + j * dim );
          }

          dbl_ptr[ i ] = sum;
        }

        /* check for interrupts */
        num_ops += dim * ( noise_dim * ( 8.0 * dim + 2.0 ) + 2.0 );
        if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
          MEMLIST_FREE_ON_ERROR( err, &mlist );
          MUTIL_ERROR( "user interrupt" );
          return MUTIL_ERR_INTERRUPT;
        }

      } /* else --- if ( first_run ) */

    } /* for ( n = 0; n < nembed; ... ) */

  } /* for ( run = 0; run < (sint32) 2; ... ) */


  /* Remove the correction vectors from the memory management list, free
  unneeded memory, and allocate memory for the final result */
  err = memlist_member_unregister( (void*) &temp_cm, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );
  MUTIL_FREE_WARN( memlist, &mlist );

  err = matuniv_malloc(
    result,
    MATUNIV_NROW( time_series ),
    MATUNIV_NCOL( time_series ),
    MUTIL_DOUBLE );
  if ( err ) {
    MUTIL_FREE_WARN( matdbl, &temp_cm );
    return err;
  }


  /* Correct (denoise) the time series */
  res_ptr = (double*) MATUNIV_DATA( result );
  ts_ptr  = (double*) MATUNIV_DATA( time_series );
  nneig   = ts_length;
  skip    = dim * delay - 1;

  for ( n = 0; n < nneig; n++, res_ptr++, ts_ptr++ ) {

    if ( n < nembed ) {

      idx = n * rdim;
      i   = n;
      j   = (sint32) 0;

    } else {

      j = n - nembed;
      i = nembed - rdelay + ( j % rdelay );
      j = j / rdelay + 1;
    }

    dbl_ptr = temp_cm.data + idx;

    sum = (double) 0;
    d   = (sint32) 0;

    do {

      sum     += *dbl_ptr;
      dbl_ptr -= skip;

      i -= delay;
      j++;

      d++;

    } while ( (i >= 0) && (j < rdim) );

    *res_ptr = *ts_ptr - sum / (double) d;

  } /* for ( n = 0; n < nneig; ... ) */


  /* Free remaining memory */
  MUTIL_FREE_WARN( matdbl, &temp_cm );

  MUTIL_TRACE(  "Done with frauniv_filter_nonlinear_local_projection()" );

  return MUTIL_ERR_OK;
}



