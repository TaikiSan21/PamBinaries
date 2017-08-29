
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_modl.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */

#include "fra_modl.h"
#include "fra_type.h"
#include "fra_util.h"

#include "mat_assn.h"
#include "mat_set.h"
#include "mat_stat.h"
#include "mat_summ.h"
#include "mat_type.h"
#include "mat_univ.h"
#include "mat_umat.h"

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

/*
  This file contains function definitions for
  inferring deterministic structure in time series.
  The functions are declared in fra_modl.h
*/

/* Static macro definitions */

#undef LOCALDEF_EXTRACT_UMAT_ROWS
#define LOCALDEF_EXTRACT_UMAT_ROWS( MAT_PTR, START_ROW, END_ROW ) \
num_rows = END_ROW - START_ROW + 1;						          \
	if ( num_rows <= 0 ){                                         \
  MUTIL_ERROR("Number of rows to extract from result matrices must be positive" ); \
  MUTIL_FREE_WARN( memlist, &list );						      \
  return MUTIL_ERR_ILLEGAL_VALUE;                                 \
	}										                      \
num_cols = MATUNIV_NCOL( MAT_PTR );				     		      \
pd_data  = (MAT_PTR)->mat.dblmat.data;				     		  \
								     		                      \
tmp = memmove( pd_data, pd_data + START_ROW * num_cols,	          \
  num_rows * num_cols * sizeof( double ) );	                      \
MUTIL_ASSERT( tmp );						                      \
								     		                      \
err = matuniv_realloc( MAT_PTR, num_rows, num_cols );             \
MEMLIST_FREE_ON_ERROR( err, &list )

#undef LOCALDEF_CHECK_NULL_POINTER_MODL
#define LOCALDEF_CHECK_NULL_POINTER_MODL( DATA_PTR,         \
  DATA_TYPE, TYPE_PREFIX )                                  \
  err = TYPE_PREFIX ## _validate( DATA_PTR );               \
  if ( err ) return err;                                    \
  if ( DATA_PTR == ( DATA_TYPE * ) NULL ) {                 \
  MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
  return MUTIL_ERR_NULL_POINTER;                            \
  }

#undef LOCALDEF_DISTANCE2BIN
#define LOCALDEF_DISTANCE2BIN( DISTANCE, DISTANCE_MINIMUM, DISTANCE_INTERVAL ) \
 (sint32) ( DISTANCE < DISTANCE_MINIMUM ? 0 : floor( ( DISTANCE - DISTANCE_MINIMUM ) / DISTANCE_INTERVAL ) )


static mutil_errcode localfn_calculate_Dij_I(
  const sint32   i,
  const sint32   j,
  const sint32   jj,
  const sint32   embedding_dimension,
  const sint32   time_lag,
  const sint32   orbital_lag,
  const sint32   image_lag,
  const sint32   n_scale,
  const sint32   n_embed,
  const double   scale_min,
  const double   resolution,
  const double  *pd_series,
  sint32_mat    *future_bins,
  sint32_mat    *bin_count,
  sint32        *iaccess,
  sint32        *istored,
  mat_set       *result);

static mutil_errcode localfn_calculate_Dij_II(
  const sint32   i,
  const sint32   j,
  const sint32   embedding_dimension,
  const sint32   time_lag,
  const sint32   orbital_lag,
  const sint32   image_lag,
  const double   scale_min,
  const double   scale_max,
  const double   resolution,
  const double  *pd_series,
  sint32_mat    *bin_count,
  mat_set       *result);


/* Static functions declared here and defined at end of file */

/* Detecting determinism in a time series */
/*                                        */
/* Documented in fra_modl.h               */
/* Written by William Constantine         */

mutil_errcode frauniv_determinism_delta_epsilon(
 const univ_mat *time_series,
 const sint32    embedding_dimension,
 const sint32    time_lag,
 const sint32    orbital_lag,
 const sint32    image_lag,
 const double    scale_min,
 const double    scale_max,
 const double    resolution,
 const boolean   minimize,
 void           *intrp_ptr,
 mat_set        *result )
{
  boolean        is_delay_embedding;
  double         cumsum;
  double        *pd_series;
  memlist        list;
  mutil_errcode  err;
  sint32         d;
  sint32         dims = 3;
  sint32         i;
  sint32         index;
  sint32         n_epsilon;
  sint32         iaccess = 0;
  sint32         istored = 0;
  sint32         start_row;
  sint32         end_row;
  sint32         j;
  sint32         n_embed;
  sint32         P;
  sint32         Nr;
  sint32         jj;
  sint32         jmax;
  sint32         jmin;
  sint32         imax;
  sint32_mat     future_bins;
  sint32_mat     ncol;
  sint32_mat     nrow;
  sint32_mat     bin_count;
  double         xmin;
  double         xmax;
  sint32         n_scale;
  sint32         num_cols;
  sint32         num_rows;
  sint32         bin_count_sum;
  double        *pd_eps;
  double        *pd_cumeps;
  double        *pd_data;
  void          *tmp;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_determinism_delta_epsilon()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /*
    Important algorithm notes:

    The delta-epsilon algorithm works by tracking the Euclidean
    distances between neighbors in the phase space and their
    images (points in the future along the same orbits).
    The distance between two points in the phase space is denoted
    as delta, and the distance between their images is denoted as epsilon.
    Delta is used to define the "bin" that the distance metric
    epsilon is held.

    A standard algorithm used to calculate these distances would typically
    do so for every possible combination of points in the phase space.
    Thus, for an Ne point embedding (and with a unit orbital lag so that the very
    next point in the orbit is a neighbor candidate), there are Ne * (Ne - 1) / 2 unique
    pairs (deltas) that could be calculated. However, we can only use a subset
    of those pairs because for each delta we must calculate an epsilon which
    is defined as the Euclidean distance between two future points (image_lag points ahead)
    along the same orbits of the points used to form delta.
    This "look" into the future of an orbit means that we cannot use the deltas
    near the end of the orbits. 
    
    Let Dij(E) denote delta and represent the
    Euclidean distance between points i and j in an E-dimensional embedding.
    Our limitation is defined by the index j as a function of i:
    jmin = i + orbital_lag while jmax = Ne - 1 - image_lag. Thus, for each i,
    we can calculate jmax - jmin + 1 deltas (and epsilons) which is to say
    Nj = Ne - image_lag - orbital_lag - i.
    Using zero-based indexing, the largest value that i can take on is
    imax = Ne - 1 - orbital_lag - image_lag. This means that there will be
    imax + 1 "i" groups calculated, each group containing Nj deltas and epsilons.

    Doing these calculations blindly will lead to a very inefficient program
    since, after calculating i = image_lag groups, the Dij(E) will be
    redundant in that they will have been previously calculated
    as images in a previous group! This is to say that a portion
    of the epsilons calculated for the first groups can be used as deltas
    for later groups (image_lag groups ahead). For example, the number of
    redundant distances in the i = 0 group are P = Ne - orbital_lag - 2 * image_lag.
    For the i = 1 group, there are P - 1 redundant distances and, in general,
    there are P - i epsilons in group i that will become deltas image_lag groups
    later. The total number of redundant distances is based on the series summation
    1 + 2 + ... + P, which reduces to Nr = P * ( P + 1 ) / 2. Thus, in the program
    you will see a special storage system used to hold and retrieve redundant
    distances and promises to speed up the code enormously at the cost of a little
    bookkeeping.
  */

  /* check inputs for errors */

  err = frauniv_check_embedding_inputs(
    time_series, embedding_dimension, time_lag, orbital_lag,
    intrp_ptr, &is_delay_embedding, &n_embed );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( !is_delay_embedding ){
    MUTIL_ERROR( "Only delay embeddings are supported in the delta-epsilon determinism test" );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }


  if ( orbital_lag < 1 ){
    MUTIL_ERROR( "Orbital lag must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  if ( image_lag < 1 ){
    MUTIL_ERROR( "Image lag must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  if ( image_lag + orbital_lag >= n_embed ){
    MUTIL_ERROR( "The sum of the image and orbital lags "       
      "must not exceed the number of points in the embedding" );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  if ( n_embed < 2 ){
    MUTIL_ERROR( "The number of points in the embedding must exceed unity" );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  if ( scale_min < 0.0 ){
    MUTIL_ERROR( "Minimum scale must be positive or zero" );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  if ( scale_max <= 0.0 ){
    MUTIL_ERROR( "Maximum scale must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  if ( scale_max <= scale_min ){
    MUTIL_ERROR( "Maximum scale must be larger than minimum scale" );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  if ( MUTIL_INTERRUPT( 3.0 * n_embed, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_INTERRUPT;
  }

  /*
     Initialize variables:

     P is the number of redundant distances calculated in the
     first (i=0) group.

     Nr is the total number of redundant distances for all groups
     and thus is used to form the dimension of the
     redundant storage vector
  */

  P  = n_embed - orbital_lag - 2 * image_lag;
  Nr = P * ( P + 1 ) / 2;

  /* initialize scale attributes */

  err = matdbl_range( &( time_series->mat.dblmat ), intrp_ptr, &xmin, &xmax );
  MEMLIST_FREE_ON_ERROR( err, &list );


  n_scale = LOCALDEF_DISTANCE2BIN( scale_max, scale_min, resolution ) + 1;

  /* allocate memory for the matrix set: scale vector, epsilon matrix,
     and cumulative epsilon matrix */

  err = mats32_malloc_register( &nrow, 1, 3, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &ncol, 1, 3, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* ... the scale vector */

  nrow.data[0] = n_scale;
  ncol.data[0] = 1;

  /* ... the epsilon matrix */

  nrow.data[1] = n_scale;
  ncol.data[1] = embedding_dimension;

  /* ... the cumulative sum over epsilon matrix */

  nrow.data[2] = n_scale;
  ncol.data[2] = embedding_dimension;

  err = matset_malloc_register( result, 1, &dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices_arbitrary_size(
    result, &nrow, &ncol, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* ... allocate space for the future bins (storage for redundant deltas) */

  err = mats32_malloc_register( &future_bins, Nr, embedding_dimension, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* ... allocate space for the matrix used to coun the number of epsilons recorded
     in each bin */

  err = mats32_malloc_register( &bin_count, n_scale, embedding_dimension, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create pointers */

  pd_series = time_series->mat.dblmat.data;

  /* initialize matrices */

  err = matdbl_assign_scalar( 0.0, intrp_ptr, &( result->mats[1].mat.dblmat ) );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_assign_scalar( 0.0, intrp_ptr, &( result->mats[2].mat.dblmat ) );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_assign_scalar( 0, intrp_ptr, &bin_count );
  MEMLIST_FREE_ON_ERROR( err, &list );

  for ( i = 0; i < n_scale; i++ ){

    result->mats[0].mat.dblmat.data[ i ] = scale_min + (double) i * resolution;
  }

  /* Calculate the delta-epsilon statistics. */

  imax = n_embed - 1 - orbital_lag - image_lag;
  jmax = imax + orbital_lag;

  for ( i = 0; i <= imax ; i++ ){
    
    jmin = i + orbital_lag;
    
    for ( jj = 0, j = jmin; j <= jmax; j++, jj++ ){
      
      if (1){

        err = localfn_calculate_Dij_I(
          i,
          j,
          jj,
          embedding_dimension,
          time_lag,
          orbital_lag,
          image_lag,
          n_scale, 
          n_embed,
          scale_min,
          resolution,
          pd_series,
          &future_bins,
          &bin_count,
          &iaccess,
          &istored,
          result);
      }
      else{

        err = localfn_calculate_Dij_II(
          i,
          j,
          embedding_dimension,
          time_lag,
          orbital_lag,
          image_lag,
          scale_min,
          scale_max,
          resolution,
          pd_series,
          &bin_count,
          result);
      }
      MEMLIST_FREE_ON_ERROR( err, &list );
    
              
    } /* end loop over j */
    
  } /* end loop over i */
  
  /* form bin and cumulative average over each dimension (column)
     of the epsilon matrix */

  pd_eps    = result->mats[ 1 ].mat.dblmat.data;
  pd_cumeps = result->mats[ 2 ].mat.dblmat.data;
  
  for ( d = 0; d < embedding_dimension; d++ ){
    
    cumsum = 0.0;
    bin_count_sum = 0;
    
    for ( i = 0; i < n_scale ; i++ ){
      
      index = i * embedding_dimension + d;
      
      n_epsilon = bin_count.data[ index ];
      
      cumsum += pd_eps[ index ];
      
      bin_count_sum += n_epsilon;
      
      /* record bin average */
      
      if ( n_epsilon > 0 ){
        
        pd_eps[ index ] /= (double) n_epsilon;
      }
      
      /* record cumulative bin average */
      
      if ( bin_count_sum > 0 ){
        
        pd_cumeps[ index ] = cumsum / (double) bin_count_sum;
      }
    }
  }

  if ( minimize ){
    
    /* trim the resulting matrices down so
    that the data corresponding to the
    smallest scale in the first dimension
    contains a nonzero value. */
    
    start_row = 0;
    end_row = n_scale - 1;
    
    for ( i = 0; i < n_scale; i++ ){
      
      index = i * embedding_dimension;
      
      if ( bin_count.data[ index ] > 0 ){
        start_row = i;
        break;
      }
    }
    
    for ( i = n_scale - 1; i > 0; i-- ){
      
      index = ( i + 1 ) * embedding_dimension - 1;
      
      if ( bin_count.data[ index ] > 0 ){
        end_row = i;
        break;
      }
    }

    /* copy the latter part of the data to the
       location where the memory originally began.
       make sure to use memmove, so that the overlap
       is handled correctly ... */

    if ( start_row > 0 ){
      LOCALDEF_EXTRACT_UMAT_ROWS( &(result->mats[0]), start_row, end_row );
      LOCALDEF_EXTRACT_UMAT_ROWS( &(result->mats[1]), start_row, end_row );
      LOCALDEF_EXTRACT_UMAT_ROWS( &(result->mats[2]), start_row, end_row );
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

  MUTIL_TRACE( "Done with frauniv_determinism_delta_epsilon()" );

  return MUTIL_ERR_OK;
}


static mutil_errcode localfn_calculate_Dij_I(
  const sint32   i,
  const sint32   j,
  const sint32   jj,
  const sint32   embedding_dimension,
  const sint32   time_lag,
  const sint32   orbital_lag,
  const sint32   image_lag,
  const sint32   n_scale,
  const sint32   n_embed,
  const double   scale_min,
  const double   resolution,
  const double  *pd_series,
  sint32_mat    *future_bins,
  sint32_mat    *bin_count,
  sint32        *iaccess,
  sint32        *istored,
  mat_set       *result)
{
  sint32 P;
  sint32 Nr;
  sint32 d;
  sint32 bin;
  sint32 index;
  sint32 lag;
  double delta;
  double diff;
  double diffsum= 0.0;
  double epsilon;
  double imagediff;
  double imagediffsum = 0.0;

  P  = n_embed - orbital_lag - 2 * image_lag;
  Nr = P * ( P + 1 ) / 2;
    
  for ( d = 0; d < embedding_dimension; d++ ){
    
    lag = d * time_lag;
    
    /* locate the proper bin to store epsilon */
    
    if ( i < image_lag ){
      
      diff = pd_series[ i + lag ] - pd_series[ j + lag ];
      
      diffsum += diff * diff;
      
      delta = sqrt( diffsum );
      
      /*
      given delta,
      the proper bin to place epsilon is determined by
      floor( (delta - scale_min)/ resolution ) where
      delta_min is the minimum distance between any two points
      in the phase space (at a given embedding dimension),
      
        In the case where delta < delta_min, bin is set
        to bin = 0 to avoid negative indexing.
      */
      
      bin = LOCALDEF_DISTANCE2BIN( delta, scale_min, resolution );
    }
    else{
      
      bin = future_bins->data[ (*iaccess)++ ];
    }
    
    /* calculate epsilon and store corresponding
    bin if necessary for future delta bin */
    
    imagediff = pd_series[ i + lag + image_lag ] -
      pd_series[ j + lag + image_lag ];
    
    imagediffsum += imagediff * imagediff;
    
    epsilon = sqrt( imagediffsum );
    
    if ( jj < P - i ){
      
      future_bins->data[ (*istored)++ ] = LOCALDEF_DISTANCE2BIN( epsilon, scale_min, resolution );
    }
    
    /* store epsilon if the corresponding delta does not exceed the maximum scale */
    
    if ( bin < n_scale ){
      
      index = bin * embedding_dimension + d;
      
      result->mats[ 1 ].mat.dblmat.data[ index ] += epsilon;
      
      bin_count->data[ index ]++;
    }
 /*   else{
      
      printf("\nbin >= n_scale");
    }
   */ 
    
  } /* end loop thru each dimension */
  
  return MUTIL_ERR_OK;
}


static mutil_errcode localfn_calculate_Dij_II(
  const sint32   i,
  const sint32   j,
  const sint32   embedding_dimension,
  const sint32   time_lag,
  const sint32   orbital_lag,
  const sint32   image_lag,
  const double   scale_min,
  const double   scale_max,
  const double   resolution,
  const double  *pd_series,
  sint32_mat    *bin_count,
  mat_set       *result)
{
  sint32 d;
  sint32 bin;
  sint32 index;
  sint32 lag;
  double delta;
  double diff;
  double diffsum= 0.0;
  double epsilon;
  double imagediff;
  double imagediffsum = 0.0;

  for ( d = 0; d < embedding_dimension; d++ ){
    
    lag = d * time_lag;
    
    /* locate the proper bin to store epsilon */
    
    
    diff = pd_series[ i + lag ] - pd_series[ j + lag ];
    
    diffsum += diff * diff;
    
    delta = sqrt( diffsum );
    
    /* no need to continue onto higher dimensions
    if delta exceeds the maximum scale since we will
    only possibly add to the distance */

    if (delta > scale_max)
      break;
    
    /*
    given delta,
    the proper bin to place epsilon is determined by
    floor( (delta - scale_min)/ resolution ) where
    delta_min is the minimum distance between any two points
    in the phase space (at a given embedding dimension),
    
      In the case where delta < delta_min, bin is set
      to bin = 0 to avoid negative indexing.
    */
    
    bin = LOCALDEF_DISTANCE2BIN( delta, scale_min, resolution );

    /* calculate epsilon */
    
    imagediff = pd_series[ i + lag + image_lag ] -
      pd_series[ j + lag + image_lag ];
    
    imagediffsum += imagediff * imagediff;
    
    epsilon = sqrt( imagediffsum );
    
    index = bin * embedding_dimension + d;
      
    result->mats[ 1 ].mat.dblmat.data[ index ] += epsilon;
      
    bin_count->data[ index ]++;
    
  } /* end loop thru each dimension */
  
  return MUTIL_ERR_OK;
}
