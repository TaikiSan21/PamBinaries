
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_neig.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */


/* This file contains definitions for functions used to find neighbors    */
/* of points in N-dimensional Euclidian space. The functions are declared */
/* in fra_neig.h.                                                         */


#include "fra_neig.h"

#include "mat_umat.h"
#include "mat_univ.h"
#include "mat_sort.h"

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_kdtr.h"
#include "ut_math.h"
#include "ut_mem.h"

#include <stdio.h>


/* Macro definitions */
#undef LOCALDEF_CHECK_NULL_POINTER_NEIGHBOR
#define LOCALDEF_CHECK_NULL_POINTER_NEIGHBOR( DATA_PTR,     \
  DATA_TYPE, TYPE_PREFIX )                                  \
  err = TYPE_PREFIX ## _validate( DATA_PTR );               \
  if ( err ) return err;                                    \
  if ( DATA_PTR == ( DATA_TYPE * ) NULL ) {                 \
  MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
  return MUTIL_ERR_NULL_POINTER;                            \
  }


#define FRA_DISTANCE_CONST       1e20
#define FRA_WORK_CONST           1e54


/*
 ****************************
 *                          *
 * STATIC (local) VARIABLES *
 *                          *
 ****************************
 */

/* Functions frauniv_neighbor_find() and  frauniv_neighbor_find_arbitrary() */
/* share local functions which need to know from which they were called.    */
static boolean global_called_from_arbitrary;


/* Tells local functions whether the neighbor search criteria is (1) */
/* within a given distance or (2) for a certain number.              */
static boolean global_searching_by_distance;


/* These global variables used to point to memory shared by many of the */
/* functions and use fewer arguments in recursive functions.            */
static double   *global_kdtree_data;
static double   *global_current_point;
static double   *global_medians;
static double   *global_distances;
static double   *global_work1;
static double   *global_work2;
static double   *global_mahalanobis;
static double    global_maxdist;
static sint32   *global_row_index;
static sint32   *global_split_index;
static sint32   *global_indices;
static sint32    global_kdtree_npoints;
static sint32    global_dim;
static sint32    global_k;
static sint32    global_metric;
static sint32    global_bucketsize;



/*
 ****************************************
 *                                      *
 * STATIC (local) FUNCTION DECLARATIONS *
 *                                      *
 ****************************************
 */

static double localfn_distance( register sint32 j );

static boolean localfn_ballwithinbounds( void );

static boolean localfn_boundsoverlapball( void );

static boolean localfn_searchkdtree_distance(
  const sint32   i,
  const sint32   nrow );

static boolean localfn_searchkdtree_number(
  const sint32   i,
  const sint32   nrow );

static mutil_errcode localfn_get_neighbors_distance(
  const univ_mat             *points,
  const mutil_kdtree         *kdtree,
  const double                distance_max,
  const fra_distance_metric   distance_metric,
  univ_mat                   *indices,
  univ_mat                   *distances );

static mutil_errcode localfn_get_neighbors_number(
  const univ_mat             *points,
  const mutil_kdtree         *kdtree,
  const sint32                n_neighbor,
  const fra_distance_metric   distance_metric,
  univ_mat                   *indices,
  univ_mat                   *distances );

static mutil_errcode localfn_threshold_with_orbital_lag(
  univ_mat   *orig_idx,
  univ_mat   *neig_idx,
  univ_mat   *distances,
  sint32      orb_lag,
  memlist    *list );

static mutil_errcode localfn_sort_distances(
  univ_mat   *orig_idx,
  univ_mat   *neig_idx,
  univ_mat   *neig_dist );

static void localfn_correct_k_neighbors(
  univ_mat         *indices,
  const univ_mat   *distances,
  const sint32      divisor );



/*
 *********************************
 *                               *
 * LIBRARY FUNCTION DEFINITIONS  *
 *                               *
 *********************************
 */

/* Nearest-neighbor search in a multidimensional      */
/* embedding.                                         */
/*                                                    */
/* Documented in fra_neig.h                           */
/* Written by William Constantine & Keith L. Davidson */
mutil_errcode frauniv_neighbor_find(
  const univ_mat            *embedding,
  const sint32               n_neighbor,
  const double               distance_max,
  const fra_distance_metric  distance_metric,
  const univ_mat            *pdmatrix,
  const boolean              sort_distances,
  const sint32               orbital_lag,
  void                      *intrp_ptr,
  univ_mat                  *original_indices,
  univ_mat                  *neighbor_indices,
  univ_mat                  *neighbor_distances )
{
  sint32          divisor;
  sint32          index;
  sint32          npoints;
  sint32          i;
  boolean         use_radius;
  memlist         mlist;
  mutil_kdtree    kdtree;
  mutil_errcode   err;

  MUTIL_TRACE( "Start frauniv_neighbor_find()" );

  /* Avoid lint warnings */

  ( void ) whatssi;

  ( void ) intrp_ptr;

  /* Initialize memory list */

  MEMLIST_INIT( mlist );

  /* Set a global variable used by local functions */

  global_called_from_arbitrary = FALSE;

  /* initialize variables */

  use_radius = (boolean) ( n_neighbor <= 0 );


  /* Check inputs for errors */

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_NEIGHBOR( embedding, univ_mat, matuniv );

  /* ... if the type is double */

  if ( embedding->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Input embedding matrix must be of type MUTIL_DOUBLE" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  npoints = MATUNIV_NROW( embedding );

  /* ... for valid distance_metric */

  switch( distance_metric ) {

  case FRA_DISTANCE_L1:
  case FRA_DISTANCE_L2:
  case FRA_DISTANCE_LINFINITY:
    break;

  case FRA_DISTANCE_MAHALANOBIS:

    LOCALDEF_CHECK_NULL_POINTER_NEIGHBOR( pdmatrix, univ_mat, matuniv );

    if ( MATUNIV_NROW( pdmatrix ) != MATUNIV_NCOL( embedding ) ||
      MATUNIV_NCOL( pdmatrix ) != MATUNIV_NCOL( embedding ) ) {
      MUTIL_ERROR( "Input pdmatrix has incompatible dimensions" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* NOTE: There should be here a test for positive definiteness of input
       pdmatrix. The documentation warns of no such check. */

    /* set a global pointer to the matrix data */
    global_mahalanobis = (double*) MATUNIV_DATA( pdmatrix );

    break;

  default:
    MUTIL_ERROR( "Input distance_metric is invalid" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* check for enough points */

  if ( n_neighbor > embedding->mat.dblmat.nrow ) {
    MUTIL_ERROR( "Input n_neighbor can not be larger than "
      "the number of points in the embedding" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* check for valid search criteria */

  if ( ( distance_max <= (double) 0 ) && ( n_neighbor <= (sint32) 0 ) ) {
    MUTIL_ERROR( "Either distance_max or n_neighbor must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( orbital_lag >= MATUNIV_NROW( embedding ) ) {
    MUTIL_ERROR( "Input orbital_lag must be smaller than the "
      "number of available points" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( n_neighbor > MATUNIV_NROW( embedding ) ){
    MUTIL_ERROR( "n_neighbor is larger than number of available points" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }


  /* Create a kd-tree from the set of points and register it */
  /* the memory management list                              */

  err = mutil_kdtree_malloc_register( &kdtree, embedding, (sint32) 1, &mlist );
  if ( err ) {
    MUTIL_ERROR( "Could not create a kd-tree using the input embedding" );
    return err;
  }

  /* If the number of neighbors to find is specified
  then call the correpsonding routine. otherwise,
  call the routine that finds all the neighbors
  within the specified maximum distance */

  if ( !use_radius ) {

    global_searching_by_distance = FALSE;

    if ( orbital_lag > 0 ) {

      divisor = n_neighbor + 2 * orbital_lag - 1;

      if ( divisor > kdtree.npoints ) {
        MUTIL_ERROR( "Not enough points for given n_neighbor "
          "and orbital_lag" );
        MUTIL_FREE_WARN( memlist, &mlist );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }

    } else {

      divisor = n_neighbor;

    }

    /* find specified number of neighbors */

    err = localfn_get_neighbors_number(
      embedding,
      &kdtree,
      divisor,
      distance_metric,
      neighbor_indices,
      neighbor_distances );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    /* At this point, particularly for the case when the
       embedding is 1-D and there are repeated points,
       the search algorithm may not have found a point to be its
       own neighbor. We fix that here by replacing the index
       corresponding to the largest distance with the appropriate
       "self" index. We replace the index corresponding to the
       largest distance */
    localfn_correct_k_neighbors(
      neighbor_indices,
      neighbor_distances,
      divisor );

    /* free memory for the kd-tree */

    err = memlist_member_free( (void*) &kdtree, &mlist );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, neighbor_indices );
      MUTIL_FREE_WARN( matuniv, neighbor_distances );
      MEMLIST_FREE_ON_ERROR( err, &mlist );
    }

    /* register indices and distances with the memory management list */

    err = memlist_member_register( &mlist, (void*) neighbor_indices,
      MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, neighbor_indices );
      MUTIL_FREE_WARN( matuniv, neighbor_distances );
      MEMLIST_FREE_ON_ERROR( err, &mlist );
    }
    err = memlist_member_register( &mlist, (void*) neighbor_distances,
      MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, neighbor_distances );
      MEMLIST_FREE_ON_ERROR( err, &mlist );
    }

    /* fill in the vector of original indices */

    err = matuniv_malloc_register(
      original_indices,
      MATUNIV_NROW( neighbor_indices ),
      MATUNIV_NCOL( neighbor_indices ),
      MUTIL_SINT32,
      &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    for ( i = 0; i < original_indices->mat.s32mat.nelem; i++ )
      original_indices->mat.s32mat.data[ i ] = i / divisor;

    /* if a non-zero orbital lag was given then threshold the neighbors */

    if ( orbital_lag > (sint32) 0 ) {

      err = localfn_threshold_with_orbital_lag(
        original_indices,
        neighbor_indices,
        neighbor_distances,
        orbital_lag,
        &mlist );

      if ( err == MUTIL_ERR_ILLEGAL_SIZE ) { /* no neighbors survived */

        MUTIL_FREE_WARN( memlist, &mlist );
        return MUTIL_ERR_ZERO_NEIGHBORS_FOUND;

      } else if ( err ) {

        MUTIL_ERROR( "Error while thresholding neighbors" );
        MEMLIST_FREE_ON_ERROR( err, &mlist );
      }
    }
  }
  else{

    global_searching_by_distance = TRUE;

    /* find all neighbors within specified maximum distance */

    err = localfn_get_neighbors_distance(
      embedding,
      &kdtree,
      distance_max,
      distance_metric,
      neighbor_indices,
      neighbor_distances );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    /* free memory for the kd-tree */

    err = memlist_member_free( (void*) &kdtree, &mlist );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, neighbor_indices );
      MUTIL_FREE_WARN( matuniv, neighbor_distances );
      MEMLIST_FREE_ON_ERROR( err, &mlist );
    }

    /* register indices and distances with the memory management list */

    err = memlist_member_register( &mlist, (void*) neighbor_indices,
      MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, neighbor_indices );
      MUTIL_FREE_WARN( matuniv, neighbor_distances );
      MEMLIST_FREE_ON_ERROR( err, &mlist );
    }
    err =  memlist_member_register( &mlist, (void*) neighbor_distances,
      MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, neighbor_distances );
      MEMLIST_FREE_ON_ERROR( err, &mlist );
    }

    /* fill in the vector of original indices */

    err = matuniv_malloc_register(
      original_indices,
      MATUNIV_NROW( neighbor_indices ),
      MATUNIV_NCOL( neighbor_indices ),
      MUTIL_SINT32,
      &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    for ( index = 0, i = 0; i < original_indices->mat.s32mat.nelem; i++ ) {

      if ( neighbor_indices->mat.s32mat.data[ i ] == (sint32) -1 ) {

        original_indices->mat.s32mat.data[ i ] = (sint32) -1;
        index++;

      } else
        original_indices->mat.s32mat.data[ i ] = index;
    }

    /* if a non-zero orbital delay was given then threshold the neighbors */

    if ( orbital_lag > (sint32) 0 ) {

      err = localfn_threshold_with_orbital_lag(
        original_indices,
        neighbor_indices,
        neighbor_distances,
        orbital_lag,
        &mlist );

      if ( err == MUTIL_ERR_ILLEGAL_SIZE ) { /* no neighbors survived */

        MUTIL_FREE_WARN( memlist, &mlist );
        return MUTIL_ERR_ZERO_NEIGHBORS_FOUND;

      } else if ( err ) {

        MUTIL_ERROR( "Error while thresholding neighbors" );
        MEMLIST_FREE_ON_ERROR( err, &mlist );
      }
    }

    /* remove the -1 indices and corresponding elements of distances */

    { /* begin program block */

      sint32   *oidx_ptr  = (sint32*) MATUNIV_DATA( original_indices );
      sint32   *nidx_ptr  = (sint32*) MATUNIV_DATA( neighbor_indices );
      double   *ndist_ptr = (double*) MATUNIV_DATA( neighbor_distances );

      const sint32 nelem  = MATUNIV_NELEM( original_indices );

      for ( index = 0, i = 0; i < nelem; i++ ) {

        if ( nidx_ptr[ i ] != (sint32) -1 ) {

          nidx_ptr[ index ]  = nidx_ptr[ i ];
          oidx_ptr[ index ]  = oidx_ptr[ i ];
          ndist_ptr[ index ] = ndist_ptr[ i ];
          index++;
        }
      }
    } /* end program block */

    if ( matuniv_realloc_register( original_indices, index, (sint32) 1, &mlist ) ||
      matuniv_realloc_register( neighbor_indices, index, (sint32) 1, &mlist ) ||
      matuniv_realloc_register( neighbor_distances, index, (sint32) 1, &mlist ) ) {

      MUTIL_FREE_WARN( memlist, &mlist );
      return MUTIL_ERR_MEM_ALLOC;
    }

  } /* if ( !use_radius ) */


  /* Sort the distances, if the user so desires */

  if ( sort_distances ) {
    err = localfn_sort_distances(
      original_indices,
      neighbor_indices,
      neighbor_distances );
    MEMLIST_FREE_ON_ERROR( err, &mlist );
  }

  /* In the case where a specific number of neighbors were requested and
     orbital lag is greater than zero there may be more than n-neighbors
     found at this point. So, we truncate the matrices holding the
     indices and distances. */
  if ( !global_searching_by_distance &&
    MATUNIV_NELEM( neighbor_distances ) > n_neighbor * npoints ) {

    sint32    curr_idx;

    sint32    run_idx     = (sint32) 0;
    sint32    idx_counter = (sint32) 0;
    sint32   *odat        = (sint32*) MATUNIV_DATA( original_indices );
    sint32   *ndat        = (sint32*) MATUNIV_DATA( neighbor_indices );
    double   *ddat        = (double*) MATUNIV_DATA( neighbor_distances );

    const sint32 nelem    = MATUNIV_NELEM( original_indices );

    curr_idx = *odat;

    do {

      if ( original_indices->mat.s32mat.data[ run_idx ] == curr_idx ) {

        if ( idx_counter < n_neighbor ) {

          *odat = original_indices->mat.s32mat.data[ run_idx ];
          *ndat = neighbor_indices->mat.s32mat.data[ run_idx ];
          *ddat = neighbor_distances->mat.dblmat.data[ run_idx ];

          odat++;
          ndat++;
          ddat++;
        }

        idx_counter++;
      }
      else {

        curr_idx = original_indices->mat.s32mat.data[ run_idx ];
        idx_counter = 0;
        continue;
      }

      run_idx++;

    } while ( run_idx < nelem );

    index = npoints * n_neighbor;

    err = matuniv_realloc_register( original_indices, index, (sint32) 1, &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    err = matuniv_realloc_register( neighbor_indices, index, (sint32) 1, &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    err = matuniv_realloc_register( neighbor_distances, index, (sint32) 1, &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );
  }

  /* remove matrices from the memory management list */

  err = memlist_member_unregister( (void*) neighbor_indices, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = memlist_member_unregister( (void*) neighbor_distances, &mlist );
  if ( err ) {
    MUTIL_FREE_WARN( memlist, &mlist );
    MUTIL_FREE_WARN( matuniv, neighbor_indices );
    return err;
  }
  err = memlist_member_unregister( (void*) original_indices, &mlist );
  if ( err ) {
    MUTIL_FREE_WARN( memlist, &mlist );
    MUTIL_FREE_WARN( matuniv, neighbor_indices );
    MUTIL_FREE_WARN( matuniv, neighbor_distances );
    return err;
  }

  /* Free other malloced space and
  corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &mlist );

  MUTIL_TRACE( "Done with frauniv_neighbor_find()" );

  return MUTIL_ERR_OK;
}



/* Nearest-neighbor search finding neighbors from an  */
/* arbitrary set of points in the phase space.        */
/*                                                    */
/* Documented in fra_neig.h                           */
/* Written by Keith L. Davidson & William Constantine */
mutil_errcode frauniv_neighbor_find_arbitrary(
  const univ_mat            *points,
  const mutil_kdtree        *kdtree,
  const sint32               n_neighbor,
  const double               distance_max,
  const fra_distance_metric  distance_metric,
  const univ_mat            *pdmatrix,
  const boolean              sort_distances,
  const sint32               orbital_lag,
  const univ_mat            *timestamp,
  void                      *intrp_ptr,
  univ_mat                  *original_indices,
  univ_mat                  *neighbor_indices,
  univ_mat                  *neighbor_distances )
{
  sint32           index;
  sint32           i;
  sint32           divisor;
  sint32           npoints;
  boolean          use_radius;
  memlist          mlist;
  mutil_errcode    err;

  MUTIL_TRACE( "Start frauniv_neighbor_find_arbitrary()" );

  /* Avoid lint warnings */
  (void) intrp_ptr;


  /* Initialize memory list */
  MEMLIST_INIT( mlist );


  /* Set a global variable used by local functions */
  global_called_from_arbitrary = TRUE;


  /* Set search criteria */
  use_radius = (boolean) ( n_neighbor <= 0 );

  /* Check inputs for errors */

  /* points */
  err = matuniv_validate( points );
  if ( err ) {
    MUTIL_ERROR( "Input points not a valid universal matrix" );
    return err;
  }
  if ( points->type != MUTIL_DOUBLE ) {
    MUTIL_ERROR( "Input points must be of type MUTIL_DOUBLE" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  npoints = MATUNIV_NROW( points );

  /* kd-tree */
  err = mutil_kdtree_validate( kdtree );
  if ( err ) {
    MUTIL_ERROR( "Input kdtree is an invalid kd-tree structure" );
    return err;
  }

  /* make sure kdtree represents data from a space with dimension
  identical to that of embedding */
  if ( kdtree->dim != points->mat.dblmat.ncol ) {
    MUTIL_ERROR( "The number of columns in points must be "
      "identical to kdtree->dim" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* n_neighbor */
  if ( n_neighbor < (sint32) 0 ) {
    MUTIL_ERROR( "Input n_neighbor must be non-negative" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* distance_max */
  if ( distance_max < (double) 0 ) {
    MUTIL_ERROR( "Input distance_max must be non-negative" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* distance_metric */
  switch( distance_metric ) {

  case FRA_DISTANCE_L1:
  case FRA_DISTANCE_L2:
  case FRA_DISTANCE_LINFINITY:
    break;

  case FRA_DISTANCE_MAHALANOBIS:

    LOCALDEF_CHECK_NULL_POINTER_NEIGHBOR( pdmatrix, univ_mat, matuniv );

    if ( MATUNIV_NROW( pdmatrix ) != MATUNIV_NCOL( points ) ||
      MATUNIV_NCOL( pdmatrix ) != MATUNIV_NCOL( points ) ) {
      MUTIL_ERROR( "Input pdmatrix has incompatible dimensions" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* NOTE: There should be here a test for positive definiteness of input
       pdmatrix. The documentation warns of no such check. */

    /* set a global pointer to the matrix data */
    global_mahalanobis = (double*) MATUNIV_DATA( pdmatrix );

    break;

  default:
    MUTIL_ERROR( "Input distance_metric is invalid" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }


  /* check search criteria */
  if ( ( n_neighbor <= (sint32) 0 ) && ( distance_max <= (double) 0 ) ) {
    MUTIL_ERROR( "Either distance_max or n_neighbor must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* if we're searching for a specific number of neighbors then we must
  make sure that there are enough points in the kd-tree to choose from */
  if ( !use_radius && ( n_neighbor > kdtree->npoints ) ) {
    MUTIL_ERROR( "Input n_neighbor too large, not enough points in "
      "kdtree to choose from" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* orbital_lag */
  if ( ( orbital_lag < 0 ) ) {
    MUTIL_ERROR( "Input orbital_lag must be a non-negative integer" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* check for valid timestamp */
  if ( timestamp != (univ_mat*) NULL ) {

    if ( timestamp->type != MUTIL_SINT32 ) {
      MUTIL_ERROR( "Input timestamp must be of type MUTIL_SINT32" );
      return MUTIL_ERR_ILLEGAL_TYPE;
    }
    err = matuniv_validate( timestamp );
    if ( err ) {
      MUTIL_ERROR( "Input timestamp is not a valid universal matrix" );
      return err;
    }

    if ( npoints != MATUNIV_NELEM( timestamp ) ) {
      MUTIL_ERROR( "The number of elements in timestamp must equal the "
        "number of rows in points" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
    for ( i = 0; i < npoints; i++ ) {

      if ( timestamp->mat.s32mat.data[ i ] < 0 ||
        timestamp->mat.s32mat.data[ i ] > ( kdtree->npoints - 1 ) ) {

        MUTIL_ERROR( "The values in timestamp must be non-negative and "
          "less than kdtree->npoints" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
    }
  }

/* If the number of neighbors to find is specified then call the
corresponding routine. otherwise, call the routine that finds
all the neighbors within the specified maximum distance   */

if ( !use_radius ) {
  
  /* find a specified number of neighbors */
  global_searching_by_distance = FALSE;
  
  if ( orbital_lag > 0 ) {
    
    divisor = n_neighbor + 2 * orbital_lag - 1;
    
    if ( divisor > kdtree->npoints ) {
      MUTIL_ERROR( "Not enough points for given n_neighbor "
        "and orbital_lag" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
    
  } 
  else {
    
    divisor = n_neighbor;
    
  }

    err = localfn_get_neighbors_number(
      points,
      kdtree,
      divisor,
      distance_metric,
      neighbor_indices,
      neighbor_distances );
    if ( err ) {
      MUTIL_ERROR( "Error while searching for neighbors" );
      return err;
    }

    /* At this point in frauniv_neighbor_find() above we make a correction
       for the case when the neighbor search does does not find a point to
       be its own neighbor. We cannot make that correction here because
       the points for for which we are doing the neighbor search do not
       necessarily correspond to those in the kd-tree. */

    /* register indices and distances with the memory management list */
    err = memlist_member_register( &mlist, (void*) neighbor_indices, MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, neighbor_indices );
      MUTIL_FREE_WARN( matuniv, neighbor_distances );
      return err;
    }
    err =  memlist_member_register( &mlist, (void*) neighbor_distances, MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( memlist, &mlist );
      MUTIL_FREE_WARN( matuniv, neighbor_distances );
      return err;
    }

    /* allocate memory for the original indices */
    err = matuniv_malloc_register(
      original_indices,
      MATUNIV_NROW( neighbor_indices ),
      MATUNIV_NCOL( neighbor_indices ),
      MUTIL_SINT32,
      &mlist );
    if ( err ) {
      MUTIL_FREE_WARN( memlist, &mlist );
      return MUTIL_ERR_MEM_ALLOC;
    }

    /* fill in the original indices */
    if ( timestamp == (univ_mat*) NULL ) {

      /* assuming original indices are 0, 1, 2, ... */
      for ( i = 0; i < original_indices->mat.s32mat.nelem; i++ )
        original_indices->mat.s32mat.data[ i ] = i / divisor;

    } else {

      /* indices provided by user */
      for ( i = 0; i < original_indices->mat.s32mat.nelem; i++ )
        original_indices->mat.s32mat.data[ i ] =
          timestamp->mat.s32mat.data[ i / divisor ];
    }

    /* if a non-zero orbital lag was given then threshold the neighbors */
    if ( orbital_lag > (sint32) 0 ) {

      err = localfn_threshold_with_orbital_lag(
        original_indices,
        neighbor_indices,
        neighbor_distances,
        orbital_lag,
        &mlist );


      if ( err == MUTIL_ERR_ILLEGAL_SIZE ) {

        /* no neighbors survived thresholding. free all memory and
           set output matrix structure fields to zero or NULL */

        MUTIL_FREE_WARN( memlist, &mlist );
        return MUTIL_ERR_ZERO_NEIGHBORS_FOUND;

      } else if ( err ) {

        MUTIL_ERROR( "Error while thresholding neighbors" );
        MUTIL_FREE_WARN( memlist, &mlist );
        return err;
      }
    }

  } else { /* if ( !use_radius ) */

    /* find all neighbors within specified maximum distance */
    global_searching_by_distance = TRUE;

    err = localfn_get_neighbors_distance(
      points,
      kdtree,
      distance_max,
      distance_metric,
      neighbor_indices,
      neighbor_distances );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    /* register indices and distances with the memory management list */
    err = memlist_member_register( &mlist, (void*) neighbor_indices,
      MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( matuniv, neighbor_indices );
      MUTIL_FREE_WARN( matuniv, neighbor_distances );
      return err;
    }
    err =  memlist_member_register( &mlist, (void*) neighbor_distances,
      MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_FREE_WARN( memlist, &mlist );
      MUTIL_FREE_WARN( matuniv, neighbor_distances );
      return err;
    }

    /* allocate memory for the original indices */
    err = matuniv_malloc_register(
      original_indices,
      MATUNIV_NROW( neighbor_indices ),
      MATUNIV_NCOL( neighbor_indices ),
      MUTIL_SINT32,
      &mlist );
    if ( err ) {
      MUTIL_FREE_WARN( memlist, &mlist );
      return MUTIL_ERR_MEM_ALLOC;
    }

    /* fill in original indices */
    for ( index = 0, i = 0; i < original_indices->mat.s32mat.nelem; i++ ) {

      if ( neighbor_indices->mat.s32mat.data[ i ] == (sint32) -1 ) {

        original_indices->mat.s32mat.data[ i ] = (sint32) -1;
        index++;

      } else {

        if ( timestamp == (univ_mat*) NULL )
          original_indices->mat.s32mat.data[ i ] = index;
        else
          original_indices->mat.s32mat.data[ i ] =
            timestamp->mat.s32mat.data[ index ];
      }
    }

    /* if a non-zero orbital delay was given then threshold the neighbors */
    if ( orbital_lag > (sint32) 0 ) {
      
      err = localfn_threshold_with_orbital_lag(
        original_indices,
        neighbor_indices,
        neighbor_distances,
        orbital_lag,
        &mlist );
      
      if ( err == MUTIL_ERR_ILLEGAL_SIZE ) { /* no neighbors survived */
        
        MUTIL_FREE_WARN( memlist, &mlist );
        return MUTIL_ERR_ZERO_NEIGHBORS_FOUND;
        
      } else if ( err ) {
        
        MUTIL_ERROR( "Error while thresholding neighbors" );
        MUTIL_FREE_WARN( memlist, &mlist );
        return err;
      }
    }

    /* remove the -1 indices and corresponding elements of distances */
    {
      /* begin program block */
      sint32   *oidx_ptr  = (sint32*) MATUNIV_DATA( original_indices );
      sint32   *nidx_ptr  = (sint32*) MATUNIV_DATA( neighbor_indices );
      double   *ndist_ptr = (double*) MATUNIV_DATA( neighbor_distances );

      const sint32 nelem  = MATUNIV_NELEM( original_indices );

      for ( index = 0, i = 0; i < nelem; i++ ) {

        if ( nidx_ptr[ i ] != (sint32) -1 ) {

          nidx_ptr[ index ]  = nidx_ptr[ i ];
          oidx_ptr[ index ]  = oidx_ptr[ i ];
          ndist_ptr[ index ] = ndist_ptr[ i ];
          index++;
        }
      }

    } /* end program block */

      
      /* if index is still 0 then all elements of indices were -1 and no
       neighbors were found */
    if ( index == (sint32) 0 ) {
      MUTIL_FREE_WARN( memlist, &mlist );
      return MUTIL_ERR_ZERO_NEIGHBORS_FOUND;
    }

    /* resize the indices and distances now that the -1s have been removed */
    if ( matuniv_realloc_register( original_indices, index, (sint32) 1, &mlist ) ||
      matuniv_realloc_register( neighbor_indices, index, (sint32) 1, &mlist ) ||
      matuniv_realloc_register( neighbor_distances, index, (sint32) 1, &mlist ) ) {

      MUTIL_FREE_WARN( memlist, &mlist );
      return MUTIL_ERR_MEM_ALLOC;
    }

  } /* if ( !use_radius ) */


  /* Sort the distances, if the user so desires */
  if ( sort_distances || ( !global_searching_by_distance &&
    MATUNIV_NELEM( neighbor_distances ) > ( n_neighbor * npoints ) ) ) {

    err = localfn_sort_distances(
      original_indices,
      neighbor_indices,
      neighbor_distances );
    MEMLIST_FREE_ON_ERROR( err, &mlist );
  }


  /* In the case where a specific number of neighbors were requested and
     orbital lag is greater than zero there may be more than n-neighbors
     found at this point. So, we truncate the matrices holding the
      indices and distances. */
  if ( !global_searching_by_distance &&
    MATUNIV_NELEM( neighbor_distances ) > ( n_neighbor * npoints ) ) {


    sint32    curr_idx;
    sint32    run_idx     = (sint32) 0;
    sint32    idx_counter = (sint32) 0;
    sint32   *odat        = (sint32*) MATUNIV_DATA( original_indices );
    sint32   *ndat        = (sint32*) MATUNIV_DATA( neighbor_indices );
    double   *ddat        = (double*) MATUNIV_DATA( neighbor_distances );

    const sint32 nelem    = MATUNIV_NELEM( original_indices );

    curr_idx = *odat;

    do {

      if ( original_indices->mat.s32mat.data[ run_idx ] == curr_idx ) {

        if ( idx_counter < n_neighbor ) {

          *odat = original_indices->mat.s32mat.data[ run_idx ];
          *ndat = neighbor_indices->mat.s32mat.data[ run_idx ];
          *ddat = neighbor_distances->mat.dblmat.data[ run_idx ];

          odat++;
          ndat++;
          ddat++;
        }

        idx_counter++;
      }
      else {

        curr_idx = original_indices->mat.s32mat.data[ run_idx ];
        idx_counter = 0;
        continue;
      }

      run_idx++;

    } while ( run_idx < nelem );

    index = npoints * n_neighbor;

    err = matuniv_realloc_register( original_indices, index, (sint32) 1, &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    err = matuniv_realloc_register( neighbor_indices, index, (sint32) 1, &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

    err = matuniv_realloc_register( neighbor_distances, index, (sint32) 1, &mlist );
    MEMLIST_FREE_ON_ERROR( err, &mlist );

  }


  /* Remove matrices from the memory management list */
  err = memlist_member_unregister( (void*) neighbor_indices, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = memlist_member_unregister( (void*) neighbor_distances, &mlist );
  if ( err ) {
    MUTIL_FREE_WARN( memlist, &mlist );
    MUTIL_FREE_WARN( matuniv, neighbor_indices );
    return err;
  }
  err = memlist_member_unregister( (void*) original_indices, &mlist );
  if ( err ) {
    MUTIL_FREE_WARN( memlist, &mlist );
    MUTIL_FREE_WARN( matuniv, neighbor_indices );
    MUTIL_FREE_WARN( matuniv, neighbor_distances );
    return err;
  }

  /* Free all other memory */
  MUTIL_FREE_WARN( memlist, &mlist );

  MUTIL_TRACE( "Done with frauniv_neighbor_find_arbitrary()" );

  return MUTIL_ERR_OK;
}




/*
 *******************************
 *                             *
 * STATIC FUNCTION DEFINITIONS *
 *                             *
 *******************************
 */

/** Calculate the nearest neigbhors based on a maximum distance specification.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_neig.h
 * @source fra\_neig.h
 * @library fractal
 * @param embedding       Pointer to a pre-allocated universal matrix
 *                        of type MUTIL\_DOUBLE containing the embedding
 *                        matrix. Each column of the matrix respresents
 *                        a single coordinate of the embedding and each
 *                        row denotes the coordiantes of a single point
 *                        in the embedding.
 * @param distance_max    The maximum distance to
 *                        search relative to the current point in the
 *                        phase space.
 * @param distance_metric The metric used to define the distance between
 *                        points in the embedding. The argument is of
 *                        type \Ref{_fra_distance_metric}.
 * @param intrp_ptr       Pointer for implementation of interrupt checking.
 * @param indices         Pointer to a single-column matrix of type
 *                        MUTIL\_SINT32 which (upon return) will contain the
 *                        indices of the neighbors found for those listed in
 *                        the original\_indices vector. The memory for this
 *                        matrix is allocated within the function.
 * @param distances       Pointer to a single-column matrix of type
 *                        MUTIL\_DOUBLE which (upon return) will contain the
 *                        distances between the points denoted by the
 *                        original\_indices vectors and those denoted by the
 *                        correspoding neighbor\_indices vector.
 *                        The memory for this matrix is allocated within the
 *                        function. The distance is based on the metric
 *                        specified by distance\_metric.
 *
 * @see _fra_distance_metric
 * @see localfn_get_neighbors_number
 * @see frauniv_neighbor_find
 * @private
 */
static mutil_errcode localfn_get_neighbors_distance(
  const univ_mat            *points,
  const mutil_kdtree        *kdtree,
  const double               distance_max,
  const fra_distance_metric  distance_metric,
  univ_mat                  *indices,
  univ_mat                  *distances )
{
  double*         current_point;
  double*         work1;
  double*         work2;
  const sint32    two_x_nrow = (sint32) ( 2 * kdtree->npoints );
  sint32          final_length;
  sint32          tmp;
  sint32          j;
  sint32          m;
  sint32          m1;
  memlist         mlist;
  mutil_errcode   err;

  MUTIL_TRACE( "Start localfn_get_neighbors_distance()" );

  /* Initialize the memory management list */
  MEMLIST_INIT( mlist );


  /* Allocate work memory used by subroutines */
  err = mutil_malloc_register(
    (sint32) ( kdtree->dim * sizeof(double) ),
    (void**) &current_point,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = mutil_malloc_register(
    (sint32) ( kdtree->dim * sizeof(double) ),
    (void**) &work1,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = mutil_malloc_register(
    (sint32) ( kdtree->dim * sizeof(double) ),
    (void**) &work2,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );


  /* Allocate memory for output data */
  err = matuniv_malloc_register(
    distances,
    two_x_nrow,
    (sint32) 1,
    MUTIL_DOUBLE,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = matuniv_malloc_register(
    indices,
    two_x_nrow,
    (sint32) 1,
    MUTIL_SINT32,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );


  /* Initialize global variables */
  global_kdtree_data    = kdtree->points.data;
  global_medians        = kdtree->medians.data;
  global_current_point  = current_point;
  global_row_index      = kdtree->point_index.data;
  global_split_index    = kdtree->split_index.data;
  global_bucketsize     = kdtree->bucket_size;
  global_dim            = kdtree->dim;
  global_kdtree_npoints = kdtree->npoints;
  global_work1          = work1;
  global_work2          = work2;
  global_metric         = distance_metric;
  global_maxdist        = distance_max;


  /* Initialize local variables */
  final_length = (sint32) -1;
  m1           = two_x_nrow;


  /* Begin searching for neighbors */
  for ( j = 0; j < two_x_nrow; j++ )
    distances->mat.dblmat.data[ j ] = FRA_DISTANCE_CONST;

  for ( m = 0; m < points->mat.dblmat.nrow; m++ ) {

    if ( (sint32) ( m1 - final_length - 3 ) < kdtree->npoints ) {

      err = matuniv_realloc_register(
        indices,
        m1 + kdtree->npoints,
        (sint32) 1,
        &mlist );
      if ( err ) {
        MUTIL_FREE_WARN( memlist, &mlist );
        return MUTIL_ERR_MEM_ALLOC;
      }
      err = matuniv_realloc_register(
        distances,
        m1 + kdtree->npoints,
        (sint32) 1,
        &mlist );
      if ( err ) {
        MUTIL_FREE_WARN( memlist, &mlist );
        return MUTIL_ERR_MEM_ALLOC;
      }

      tmp = m1 + kdtree->npoints;
      for ( j = m1 + 1; j < tmp; j++ )
        distances->mat.dblmat.data[ j ] = FRA_DISTANCE_CONST;

      m1 += kdtree->npoints;
    }

    global_distances = distances->mat.dblmat.data + final_length + 1;
    global_indices   = indices->mat.s32mat.data + final_length + 1;
    global_k         = (sint32) -1;

    /* initialize work space */
    for ( j = 0; j < kdtree->dim; j++ ) {
      global_work1[ j ] = -FRA_WORK_CONST;
      global_work2[ j ] = FRA_WORK_CONST;
    }

    /* fill in current point */
    if ( global_called_from_arbitrary ) {

      tmp = m * kdtree->dim;
      for ( j = 0; j < kdtree->dim; j++ )
        global_current_point[ j ] = points->mat.dblmat.data[ tmp + j ];

    } else {

      for ( j = 0; j < kdtree->dim; j++ ) {
        global_current_point[ j ] =
          kdtree->points.data[ kdtree->npoints * j + m ];
      }
    }

    (void) localfn_searchkdtree_distance(
      ( kdtree->npoints - 1 ) / 2,
      kdtree->npoints );

    global_k++;
    final_length += global_k + 1;

    indices->mat.s32mat.data[ final_length ]   = (sint32) -1;
    distances->mat.dblmat.data[ final_length ] = (double) -1;
  }

  final_length++;

  /* Remove output data from memory management list and free
  all other memory */
  err = memlist_member_unregister( distances, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = memlist_member_unregister( indices, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  MUTIL_FREE_WARN( memlist, &mlist );


  /* The distance and indices matrices grew (after neighbors were found for
  for each point. Must reallocate memory */
  err = matuniv_realloc( indices, (sint32) 1, final_length );
  if ( err )  {
    MUTIL_FREE_WARN( matuniv, indices );
    MUTIL_FREE_WARN( matuniv, distances );
    return err;
  }
  err = matuniv_realloc( distances, (sint32) 1, final_length );
  if ( err )  {
    MUTIL_FREE_WARN( matuniv, indices );
    MUTIL_FREE_WARN( matuniv, distances );
    return err;
  }

  MUTIL_TRACE( "Done with localfn_get_neighbors_distance()" );

  return MUTIL_ERR_OK;
}


/** Calculate the nearest neigbhors based on a maximum distance specification.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_neig.h
 * @source fra\_neig.h
 * @library fractal
 * @param embedding       Pointer to a pre-allocated universal matrix
 *                        of type MUTIL\_DOUBLE containing the embedding
 *                        matrix. Each column of the matrix respresents
 *                        a single coordinate of the embedding and each
 *                        row denotes the coordiantes of a single point
 *                        in the embedding.
 * @param n_number        The number of neighbors to find.
 * @param distance_metric The metric used to define the distance between
 *                        points in the embedding. The argument is of
 *                        type \Ref{_fra_distance_metric}.
 * @param intrp_ptr       Pointer for implementation of interrupt checking.
 * @param indices         Pointer to a single-column matrix of type
 *                        MUTIL\_SINT32 which (upon return) will contain the
 *                        indices of the neighbors found for those listed in
 *                        the original\_indices vector. The memory for this
 *                        matrix is allocated within the function.
 * @param distances       Pointer to a single-column matrix of type
 *                        MUTIL\_DOUBLE which (upon return) will contain the
 *                        distances between the points denoted by the
 *                        original\_indices vectors and those denoted by the
 *                        correspoding neighbor\_indices vector.
 *                        The memory for this matrix is allocated within the
 *                        function. The distance is based on the metric
 *                        specified by distance\_metric.
 *
 * @see _fra_distance_metric
 * @see localfn_get_neighbors_distance
 * @see frauniv_neighbor_find
 * @private
 */
static mutil_errcode localfn_get_neighbors_number(
  const univ_mat             *points,
  const mutil_kdtree         *kdtree,
  const sint32                n_number,
  const fra_distance_metric   distance_metric,
  univ_mat                   *indices,
  univ_mat                   *distances )
{
  double         *current_point;
  double         *work1;
  double         *work2;
  sint32          tmp;
  sint32          m;
  sint32          j;
  memlist         mlist;
  mutil_errcode   err;

  MUTIL_TRACE( "Start localfn_get_neighbors_number()" );

  /* Initialize the memory management list */
  MEMLIST_INIT( mlist );


  /* Allocate work memory used by subroutines */
  err = mutil_malloc_register(
    (sint32) ( kdtree->dim * sizeof(double) ),
    (void**) &current_point,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = mutil_malloc_register(
    (sint32) ( kdtree->dim * sizeof(double) ),
    (void**) &work1,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = mutil_malloc_register(
    (sint32) ( kdtree->dim * sizeof(double) ),
    (void**) &work2,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );


  /* Allocate memory for output data */
  err = matuniv_malloc_register(
    distances,
    n_number * points->mat.dblmat.nrow,
    (sint32) 1,
    MUTIL_DOUBLE,
    &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = matuniv_malloc_register(
    indices,
    n_number * points->mat.dblmat.nrow,
    (sint32) 1,
    MUTIL_SINT32,
    &mlist);
  MEMLIST_FREE_ON_ERROR( err, &mlist );


  /* Initialize global variables */
  global_kdtree_data    = kdtree->points.data;
  global_medians        = kdtree->medians.data;
  global_row_index      = kdtree->point_index.data;
  global_split_index    = kdtree->split_index.data;
  global_bucketsize     = kdtree->bucket_size;
  global_kdtree_npoints = kdtree->npoints;
  global_dim            = kdtree->dim;
  global_current_point  = current_point;
  global_work1          = work1;
  global_work2          = work2;
  global_metric         = distance_metric;
  global_k              = n_number;


  /* Begin searching for neighbors */
  for ( m = 0; m < points->mat.dblmat.nrow; m++ ) {

    /* set global pointers to correct position */
    tmp              = m * n_number;
    global_distances = distances->mat.dblmat.data + tmp;
    global_indices   = indices->mat.s32mat.data + tmp;

    /* initialize work space */
    for ( j = 0; j < kdtree->dim; j++ ) {
      global_work1[ j ] = -FRA_DISTANCE_CONST;
      global_work2[ j ] = FRA_DISTANCE_CONST;
    }

    /* grab current point */
    if ( global_called_from_arbitrary ) {

      tmp = m * kdtree->dim;
      for ( j = 0; j < kdtree->dim; j++ )
        global_current_point[ j ] = points->mat.dblmat.data[ tmp + j ];

    } else {

      for ( j = 0; j < kdtree->dim; j++ ) {
        global_current_point[ j ] =
          kdtree->points.data[ kdtree->npoints * j + m ];
      }
    }

    /* initialize distances */
    for ( j = 0; j < n_number; j++ )
      global_distances[ j ] = FRA_DISTANCE_CONST;

    /* search the kd-tree for neighbors */
    (void) localfn_searchkdtree_number(
      ( kdtree->npoints - 1 ) / 2, kdtree->npoints );
  }


  /* Remove output data from memory management list and free
  all other memory */
  err = memlist_member_unregister( distances, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = memlist_member_unregister( indices, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  MUTIL_FREE_WARN( memlist, &mlist );


  MUTIL_TRACE( "Done with localfn_get_neighbors_number()" );

  return MUTIL_ERR_OK;
}



static boolean localfn_searchkdtree_distance(
  const sint32   i,
  const sint32   nrow )
{
  /*
     This is a recursive routine for finding the nearest neighbors in a
     kd-tree.  For leaves in the tree with "global_bucketsize" or fewer
     elements, the routine processes each element in the leaf to
     determine if the element is a neighbor.  Otherwise the node is
     split, and a check is made to see if the resulting leaves should
     be searched.

     The data is defined statically to avoid passing data in recursive
     calls.
     */

  static double   dtmp;
  double          tmp;
  double          p;
  static sint32   j;
  static sint32   itmp;
  static sint32   n;
  static sint32   tmp_idx;
  sint32          d;

  MUTIL_TRACE( "Start localfn_searchkdtree_distance()" );

  if ( nrow <= global_bucketsize ) { /* Check the leaf if appropriate */

    tmp_idx = i - ( nrow - 1 ) / 2 + nrow;
    for ( j = tmp_idx - nrow; j < tmp_idx; j++ ) {

      if ( ( dtmp = localfn_distance( global_row_index[ j ] ) ) <
        global_maxdist ) {

        global_k++;
        global_distances[ global_k ] = dtmp;
        global_indices[ global_k ]   = global_row_index[ j ];
      }
    }

    if ( localfn_ballwithinbounds() )
      return TRUE;  /* Check for additional leafs */
    else
      return FALSE;
  }

  d = global_split_index[ i ];
  p = global_medians[ i ];

  if ( global_current_point[ d ] <= p ) { /*  Process each leaf */

    itmp              = i - ( nrow + 1 ) / 4;
    n                 = ( nrow + 1 ) / 2;
    tmp               = global_work2[ d ];
    global_work2[ d ] = p;

    if ( localfn_searchkdtree_distance( itmp, n ) )
      return TRUE;

    global_work2[ d ] = tmp;

  } else {

    itmp              = i + ( nrow - 2 ) / 4 + 1;
    n                 = nrow / 2;
    tmp               = global_work1[ d ];
    global_work1[ d ] = p;

    if ( localfn_searchkdtree_distance( itmp, n ) )
      return TRUE;

    global_work1[ d ] = tmp;
  }

  if ( global_current_point[ d ] <= p ) {

    itmp              = i + ( nrow - 2 ) / 4 + 1;
    n                 = nrow / 2;
    tmp               = global_work1[ d ];
    global_work1[ d ] = p;

    /* Check if a search is required */
    if ( localfn_boundsoverlapball() )
      if ( localfn_searchkdtree_distance( itmp, n ) )
        return TRUE;

    global_work1[ d ] = tmp;

  } else {

    itmp              = i - ( nrow + 1 ) / 4;
    n                 = ( nrow + 1 ) / 2;
    tmp               = global_work2[ d ];
    global_work2[ d ] = p;

    if ( localfn_boundsoverlapball() )
      if ( localfn_searchkdtree_distance( itmp, n ) )
        return TRUE;

    global_work2[ d ] = tmp;

  }

  if ( localfn_ballwithinbounds() )
    return TRUE;

  MUTIL_TRACE( "Done with localfn_searchkdtree_distance()" );

  return FALSE;
}


static boolean localfn_searchkdtree_number(
  const sint32   i,
  const sint32   nrow )
{
        /*
          This is a recursive routine for finding the nearest neighbors in a
          kd-tree.  For leaves in the tree with "bucketSize" or fewer
          elements, the routine processes each element in the leaf to
          determine if the element is a neighbor.  Otherwise the node is
          split, and a check is made to see if the resulting leaves should
          be searched.

          The data is defined statically to avoid passing data in recursive
          calls.
        */

  static double   dtmp;
  double          tmp;
  double          p;
  static sint32   j;
  static sint32   l;
  static sint32   itmp;
  static sint32   n;
  static sint32   tmp_idx;
  sint32          d;

  MUTIL_TRACE( "Start localfn_searchkdtree_number()" );

  if ( nrow <= global_bucketsize)  { /* Check the leaf if appropriate */

    tmp_idx = i - ( nrow - 1 ) / 2 + nrow;

    for ( j = tmp_idx - nrow; j < tmp_idx; j++ ) {

      if ( ( dtmp = localfn_distance( global_row_index[ j ] ) ) <
        global_distances[ 0 ] ) {

        global_distances[ 0 ] = dtmp;
        global_indices[ 0 ]   = global_row_index[ j ];

        /* sort so largest is in position 1 */
        for ( l = 1; l < global_k; l++ ) {

          if ( global_distances[ l ] > global_distances[ 0 ] ) {

            tmp                   = global_distances[ 0 ];
            itmp                  = global_indices[ 0 ];
            global_distances[ 0 ] = global_distances[ l ];
            global_indices[ 0 ]   = global_indices[ l ];
            global_distances[ l ] = tmp;
            global_indices[ l ]   = itmp;
          }
        }
      }
    }

    if ( localfn_ballwithinbounds() )
      return TRUE;  /* Check for additional leafs */
    else
      return FALSE;
  }

  d = global_split_index[i];
  p = global_medians[i];

  if (global_current_point[d] <= p) { /*  Process each leaf */

    itmp              = i - ( nrow + 1 ) / 4;
    n                 = ( nrow + 1 ) / 2;
    tmp               = global_work2[ d ];
    global_work2[ d ] = p;

    if ( localfn_searchkdtree_number( itmp, n ) )
      return TRUE;

    global_work2[ d ] = tmp;

  } else {

    itmp              = i + ( nrow - 2 ) / 4 + 1;
    n                 = nrow / 2;
    tmp               = global_work1[ d ];
    global_work1[ d ] = p;

    if ( localfn_searchkdtree_number( itmp, n ) )
      return TRUE;

    global_work1[ d ] = tmp;
  }

  if ( global_current_point[ d ] <= p ) {

    itmp              = i + ( nrow - 2 ) / 4 + 1;
    n                 = nrow / 2;
    tmp               = global_work1[ d ];
    global_work1[ d ] = p;

    /* Check if a search is required */
    if ( localfn_boundsoverlapball() ) {

      if ( localfn_searchkdtree_number( itmp, n ) ) return TRUE;
    }

    global_work1[d] = tmp;

  } else {

    itmp              = i - ( nrow + 1 ) / 4;
    n                 = ( nrow + 1 ) / 2;
    tmp               = global_work2[ d ];
    global_work2[ d ] = p;

    if ( localfn_boundsoverlapball() ) {

      if (localfn_searchkdtree_number( itmp, n ) ) return TRUE;
    }

    global_work2[ d ] = tmp;

  }

  if ( localfn_ballwithinbounds() )
    return TRUE;

  MUTIL_TRACE( "Done with localfn_searchkdtree_number()" );

  return FALSE;
}


static boolean localfn_ballwithinbounds( void )
{
/* Checks to see if the furthest distance to a neighbor is entirely within
the bounds of the current partition defined by the kd-tree.  If so,
then we are all done.  Otherwise, we must look in neighboring partitions.
 */

  register double   local_dist = global_maxdist;
  register sint32   i;

  MUTIL_TRACE( "Start localfn_ballwithinbounds()" );

  if ( global_searching_by_distance == FALSE )
    local_dist = global_distances[ 0 ];

  if ( global_distances[ 0 ] < 1.0e10 ) {

    for ( i = 0; i < global_dim; i++ ) {

      if ( MUTIL_ABS( global_current_point[ i ] - global_work1[ i ] ) <=
        local_dist )
        return FALSE;
      else if ( MUTIL_ABS( global_current_point[ i ] - global_work2[ i ] ) <=
        local_dist )
        return FALSE;
    }

    return TRUE;

  }

  MUTIL_TRACE( "Done with localfn_ballwithinbounds()" );

  return FALSE;
}


static boolean localfn_boundsoverlapball( void )
{
/*This routine is used to checks if the bounds defining a partition
of the kd-tree could contain a nearest neighbor.  The variables
are defined statically to avoid passing arguments in recursive
calls.
  */

  register double   sum;
  register double   tmp;
  register double   local_dist = global_maxdist;
  register sint32   i;
  register sint32   j;

  MUTIL_TRACE( "Start localfn_boundsoverlapball()" );

  if ( global_searching_by_distance == FALSE )
    local_dist = global_distances[ 0 ];

  sum = (double) 0;

  if ( global_distances[ 0 ] > 1.0e10 )
    return TRUE;

  for ( i = 0; i < global_dim; i++ ) {

    if ( global_current_point[ i ] < global_work1[ i ] ) {

      switch( global_metric ) {

      case FRA_DISTANCE_L1: /* L1-Norm */

        sum += MUTIL_ABS( global_current_point[ i ] - global_work1[ i ] );

        if ( sum > local_dist )
          return FALSE;

        break;

      case FRA_DISTANCE_L2: /* L2-Norm */

        tmp = global_current_point[ i ] - global_work1[ i ];

        sum += tmp * tmp;

        if ( sum > ( local_dist * local_dist ) )
          return FALSE;

        break;

      case FRA_DISTANCE_LINFINITY: /* Infinity-Norm */

        if ( (tmp = MUTIL_ABS(global_current_point[ i ] -
          global_work1[ i ] ) ) > sum )
          sum = tmp;

        if ( sum > local_dist )
          return FALSE;

        break;

      case FRA_DISTANCE_MAHALANOBIS: /* Mahalanobis-Norm */

        tmp = 0;

        for ( j = 0; j < global_dim; j++ ) {
          tmp += global_mahalanobis[ j * global_dim + i ] *
            ( global_current_point[ i ] - global_work1[ i ] ) *
            ( global_current_point[ j ] - global_work1[ j ] );
        }

        sum += tmp;

        if ( sum > ( local_dist * local_dist ) )
          return FALSE;

        break;

      default:
        /* Should never get here due to prior error checking */
        MUTIL_ERROR( "Should not have reached default case in switch statement" );
        break;
      }

    } else if ( global_current_point[ i ] > global_work2[ i ] ) {

      switch( global_metric ) {

      case FRA_DISTANCE_L1: /* L1-Norm */

        sum += MUTIL_ABS( global_current_point[ i ] - global_work2[ i ] );

        if ( sum > local_dist )
          return FALSE;

        break;

      case FRA_DISTANCE_L2: /* L2-Norm */

        tmp = global_current_point[ i ] - global_work2[ i ];
        sum += tmp * tmp;

        if ( sum > ( local_dist * local_dist ) )
          return FALSE;

        break;

      case FRA_DISTANCE_LINFINITY: /* Infinity-Norm */

        if ( ( tmp = MUTIL_ABS( global_current_point[ i ] -
          global_work2[ i ] ) ) > sum )
          sum = tmp;

        if ( sum > local_dist )
          return FALSE;

        break;

      case FRA_DISTANCE_MAHALANOBIS: /* Mahalanobis-Norm */

        tmp = 0;

        for ( j = 0; j < global_dim; j++ ) {
          tmp += global_mahalanobis[ j * global_dim + i ] *
            ( global_current_point[ i ] - global_work2[ i ] ) *
            ( global_current_point[ j ] - global_work2[ j ] );
        }

        sum += tmp;

        if ( sum > ( local_dist * local_dist ) )
          return FALSE;

        break;

      default:
        /* Should never get here due to prior error checking */
        MUTIL_ERROR( "Should not have reached default case in switch statement" );
        break;
      }
    }
  }

  MUTIL_TRACE( "Done with localfn_boundsoverlapball()" );

  return TRUE;
}



static double localfn_distance( register sint32 j )
{
  /*
     Finds the distance between the key and a point in the kd-tree using
     one of four metrics.  The input variables are defined statically
     to avoid passing them as arguments in recursive function calls
     (distance is not recursive).
     */

  register double   sum;
  register double   tmp;
  register sint32   i;
  register sint32   p;
  register sint32   q;

  MUTIL_TRACE( "Start localfn_distance()" );

  sum = (double) 0;

  switch ( global_metric ) {

  case FRA_DISTANCE_L1: /* L1-Norm */

     for ( i = 0; i < global_dim; i++ )
        sum += MUTIL_ABS( global_current_point[ i ] -
           global_kdtree_data[ j + i * global_kdtree_npoints ] );

     break;

  case FRA_DISTANCE_L2: /* L2-Norm */

     for ( i = 0; i < global_dim; i++ ) {
        tmp = global_current_point[ i ] -
           global_kdtree_data[ j + i * global_kdtree_npoints ];
        sum += tmp * tmp;
     }
     sum = sqrt( sum );

     break;

  case FRA_DISTANCE_LINFINITY: /* Infinity-Norm */

     for ( i = 0; i < global_dim; i++ ) {
        tmp = MUTIL_ABS( global_current_point[ i ] -
           global_kdtree_data[ j + i * global_kdtree_npoints ] );
        if ( tmp > sum ) sum = tmp;
     }

     break;

  case FRA_DISTANCE_MAHALANOBIS: /* Mahalanobis-Norm */

    sum = 0.0;

     for ( p = 0; p < global_dim; p++ ) {

       for ( q = p; q < global_dim; q++ ) {

        tmp = global_mahalanobis[ p * global_dim + q ] *
          ( global_current_point[ p ] -
            global_kdtree_data[ j + p * global_kdtree_npoints ] );

        tmp *= ( global_current_point[ q ] -
           global_kdtree_data[ j + q * global_kdtree_npoints ] );

        sum += ( p == q ? tmp : 2.0 * tmp );
       }
     }

     sum = sqrt( sum );

     break;

  default:
    /* Should never get here due to prior error checking */
    MUTIL_ERROR( "Should not have reached default case in switch statement" );
    break;
  }

  MUTIL_TRACE( "Done with localfn_distance()" );

  return sum;
}


static mutil_errcode localfn_threshold_with_orbital_lag(
  univ_mat   *orig_idx,
  univ_mat   *neig_idx,
  univ_mat   *distances,
  sint32      orb_lag,
  memlist    *mlist )
{
  register sint32   *orig_data;
  register sint32   *neig_data;
  register double   *dist_data;
  register sint32    nelem;
  register sint32    index;
  register sint32    i;
  register sint32    local_ol;
  sint32             ones_count = (sint32) -1;

  MUTIL_TRACE( "Start localfn_threshold_with_orbital_lag()" );

  /* Initialize variables */
  orig_data = (sint32*) MATUNIV_DATA( orig_idx );
  neig_data = (sint32*) MATUNIV_DATA( neig_idx );
  dist_data = (double*) MATUNIV_DATA( distances );
  nelem     = MATUNIV_NELEM( orig_idx );
  local_ol  = orb_lag;


  /* Begin thresholding */
  if ( global_searching_by_distance ) {

    /* All neighbors are within a user-defined distance. The neighbor */
    /* indices are partitioned using -1's.                            */
    for ( ones_count = 0, index = 0, i = 0; i < nelem; i++ ) {

      if ( neig_data[ i ] == (sint32) -1 ) {

        orig_data[ index ] = (sint32) -1;
        neig_data[ index ] = (sint32) -1;
        dist_data[ index ] = (double) -1;
        index++;
        ones_count++;

      } else if ( (sint32) MUTIL_ABS( orig_data[ i ] - neig_data[ i ] )
               >= local_ol ) {

        orig_data[ index ] = orig_data[ i ];
        neig_data[ index ] = neig_data[ i ];
        dist_data[ index ] = dist_data[ i ];
        index++;
      }
    }

  } else {

    /* A specific number of neighbors were found (no -1 values) */
    for ( index = 0, i = 0; i < nelem; i++ ) {

      if ( (sint32) MUTIL_ABS( orig_data[ i ] - neig_data[ i ] )
        >= local_ol ) {

        orig_data[ index ] = orig_data[ i ];
        neig_data[ index ] = neig_data[ i ];
        dist_data[ index ] = dist_data[ i ];
        index++;
      }
    }
  } /* if ( global_searching_by_distance ) */


  /* Consider some special cases in which no neighbors */
  /* survived the thresholding.                        */
  if ( global_searching_by_distance && ( ones_count == index ) )
    return MUTIL_ERR_ILLEGAL_SIZE;
  else if ( index == (sint32) 0 )
    return MUTIL_ERR_ILLEGAL_SIZE;


  /* Resize the data matrices if all neighbors did not */
  /* survive the thresholding.                         */

  if ( index < nelem ) {
    if ( matuniv_realloc_register( orig_idx, index, (sint32) 1, mlist ) ||
      matuniv_realloc_register( neig_idx, index, (sint32) 1, mlist ) ||
      matuniv_realloc_register( distances, index, (sint32) 1, mlist ) )
      return MUTIL_ERR_MEM_ALLOC;
  }

  MUTIL_TRACE( "Done with localfn_threshold_with_orbital_lag()" );

  return MUTIL_ERR_OK;
}


static mutil_errcode localfn_sort_distances(
  univ_mat   *orig_idx,
  univ_mat   *neig_idx,
  univ_mat   *neig_dist )
{
  double_mat      dist_wrap_mat;
  double         *dist_data;
  sint32_mat      perm_indices;
  sint32         *nidx_data;
  sint32         *oidx_data;
  sint32          nelem;
  sint32          start;
  sint32          current_idx;
  sint32          current_length;
  sint32          idx_counter;
  boolean         in_bounds = TRUE;
  const boolean   one       = TRUE;
  mutil_errcode   err;

  MUTIL_TRACE( "Start localfn_sort_distances()" );

  /* Allocate memory for permutation indices. We will */
  /* grow this vector as needed.                      */
  err = mats32_malloc( &perm_indices, (sint32) 1, (sint32) 1 );
  if ( err ) return MUTIL_ERR_MEM_ALLOC;


  /* Grab the number of elements in all input matrices */
  nelem = MATUNIV_NELEM( orig_idx );


  /* We'll always wrap dist_wrap_mat around a column vector so */
  /* we can initialize the ncol field outside of any loops.    */
  dist_wrap_mat.ncol = (sint32) 1;
  (void) dist_wrap_mat.ncol; /* this supresses a lint warning */


  /* Initialize a starting index and current length of the */
  /* permutation indices vector.                           */
  start          = (sint32) 0;
  current_length = (sint32) 1;


  /* Initialize data pointers */
  oidx_data = (sint32*) MATUNIV_DATA( orig_idx );
  nidx_data = (sint32*) MATUNIV_DATA( neig_idx );
  dist_data = (double*) MATUNIV_DATA( neig_dist );


  /* Begin sorting */
  do {

    /* count the number of neighbor indices for current index */
    idx_counter = (sint32) 0;
    current_idx = oidx_data[ start ];

    while ( one ) {

      if ( oidx_data[ start + idx_counter ] == current_idx )
        idx_counter++;
      else
        break;

      if ( start + idx_counter >= nelem )
        break;
    }

    /* NOTE: idx_counter is at least 1 at this point */

    /* check to see if we should sort, i.e., if idx_counter > 1 */
    if ( !( idx_counter > (sint32) 1 ) )

      /* update our place in the input vectors */
      start++;

    else { /* more than one neighbor, we're sorting! */

      /* make sure there is enough memory for the permutation indices */
      if ( current_length < idx_counter ) {

        current_length = idx_counter;

        err = mats32_realloc( &perm_indices, current_length, (sint32) 1 );
        if ( err ) {
          MUTIL_FREE_WARN( mats32, &perm_indices );
          return MUTIL_ERR_MEM_ALLOC;
        }
      }

      /* manually wrap a matrix around the distances, placing the    */
      /* internal data pointer at the appropriate point of the input */
      /* distances -- the ncol field is initialized to 1 above       */
      dist_wrap_mat.nrow  = idx_counter;
      dist_wrap_mat.nelem = idx_counter;
      dist_wrap_mat.data  = dist_data + start;

      /* we must make the perm_indices matrix appear to be a */
      /* vector of length idx_counter. the actual length is  */
      /* current_length which is always >= idx_counter.      */
      perm_indices.nrow  = idx_counter;
      perm_indices.nelem = idx_counter;

      /* call MUTILS routine to obtain perutation indices */
      err = matdbl_sort_index_partial(
        &dist_wrap_mat,
        (sint32_mat*) NULL, /* full sort */
        (void*) NULL,
        &perm_indices );
      if ( err ) {

        /* make sure the length of perm_indices is correct before */
        /* memory deallocation                                    */
        perm_indices.nrow  = current_length;
        perm_indices.nelem = current_length;

        MUTIL_FREE_WARN( mats32, &perm_indices );
        return err;
      }

      /* now permute the current section of the input distance and    */
      /* neighbor indices. note that orig_idx does not need permuting */
      /* since all elements in the section are identical              */

      { /* begin program block */

        register double  *dist_ptr  = dist_data + start;
        register sint32  *nidx_ptr  = nidx_data + start;
        register sint32  *perm_ptr  = perm_indices.data;
        register sint32   ctr       = idx_counter;
        register sint32   i;
        double            *dist_buff;
        sint32            *nidx_buff;

        /* allocate work buffers */
        err = mutil_malloc(
          (sint32) ( ctr * sizeof(double) ),
          (void**) &dist_buff );
        if ( err ) {
          MUTIL_FREE_WARN( mats32, &perm_indices );
          return MUTIL_ERR_MEM_ALLOC;
        }
        err = mutil_malloc(
          (sint32) ( ctr * sizeof(sint32) ),
          (void**) &nidx_buff );
        if ( err ) {
          (void) mutil_free( dist_buff, sizeof( ctr * sizeof( double ) ) );
          MUTIL_FREE_WARN( mats32, &perm_indices );
          return MUTIL_ERR_MEM_ALLOC;
        }

        /* copy distances and neighbor indices into buffers */
        for ( i = 0; i < ctr; i++ ) {
          dist_buff[ i ] = dist_ptr[ i ];
          nidx_buff[ i ] = nidx_ptr[ i ];
        }

        /* permute distances and neighbor indices */
        for ( i = 0; i < ctr; i++ ) {
          dist_ptr[ i ] = dist_buff[ perm_ptr[ i ] ];
          nidx_ptr[ i ] = nidx_buff[ perm_ptr[ i ] ];
        }


        /* free buffer memory */
        (void) mutil_free( dist_buff, sizeof( ctr * sizeof( double ) ) );
        (void) mutil_free( nidx_buff, sizeof( ctr * sizeof( sint32 ) ) );

      } /* end of progam block */


      /* update our place in the input vectors */
      start += idx_counter;

    } /* if ( !( idx_counter > (sint32) 1 ) ) */


    /* are we finished? */
    if ( start >= nelem )
      in_bounds = FALSE;

    /* make sure the length of perm_indices is correct before */
    /* memory deallocation                                    */
    perm_indices.nrow  = current_length;
    perm_indices.nelem = current_length;

  } while ( in_bounds );

  /* Free memory for the permutation indices */
  MUTIL_FREE_WARN( mats32, &perm_indices );


  MUTIL_TRACE( "Done with localfn_sort_distances()" );

  return MUTIL_ERR_OK;
}



static void localfn_correct_k_neighbors(
  univ_mat         *indices,
  const univ_mat   *distances,
  const sint32      divisor )
{
  sint32   i;
  sint32   n;
  sint32   idx;
  sint32   npoints;
  sint32   max_dist_idx;
  sint32   tmp_idx;

  double   max_dist;
  boolean  self_neighbor;


  MUTIL_TRACE( "Start localfn_correct_k_neighbors()" );

  npoints = MATUNIV_NELEM( indices ) / divisor;

  for ( n = 0; n < npoints; n++ ) {

    tmp_idx      = n * divisor;
    max_dist_idx = tmp_idx;
    max_dist     = (double) 0;

    self_neighbor = FALSE;

    for ( i = 0; i < divisor; i++ ) {

      idx = tmp_idx + i;

      if ( distances->mat.dblmat.data[ idx ] > max_dist ) {

        max_dist     = distances->mat.dblmat.data[ idx ];
        max_dist_idx = idx;

      }

      if ( indices->mat.s32mat.data[ idx ] == n ) {
        self_neighbor = TRUE;
        break;
      }
    }


    if ( !self_neighbor ) indices->mat.s32mat.data[ max_dist_idx ] = n;
  }


  MUTIL_TRACE( "Done with localfn_correct_k_neighbors()" );

  return;
}
