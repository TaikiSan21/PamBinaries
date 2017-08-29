
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_lyap.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */


/* This file contains definitions for functions used to compute */
/* time delayed mutual information for a given time series.     */

#include "fra_lyap.h"
#include "fra_dim.h"

#include "fra_neig.h"
#include "fra_util.h"

#include "mat_arit.h"
#include "mat_assn.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "mat_summ.h"
#include "mat_io.h"

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_math.h"
#include "ut_mem.h"
#include "ut_rand.h"

#include<stdio.h>

/* Macros */

#define FRA_LOCAL_LOG(x) \
( ( (x) >  1.175494351e-38 ) ? log( (x) ) : log( 1.175494351e-38) )


#define LOCALDEF_SCALE_MAX( X, Y, SCALE ) \
  (SCALE) = (double) MUTIL_MAX( (X), (Y) ); \
  if ( (SCALE) == (double) 0.0 ){  \
  SCALE = 1.0; \
  }

#define LOCALDEF_LYAPUNOV_CHANGE_MATSET_DIMENSIONS( ML_POINTER, RL_POINTER, NMATS ) \
  *((ML_POINTER)->dims)  = (NMATS); \
  (ML_POINTER)->nelem    = (NMATS); \
  *((RL_POINTER)->dims)  = (NMATS); \
  (RL_POINTER)->nelem    = (NMATS); \
  /* the following suppresses lint warnings */ \
  (void) (ML_POINTER)->dims;  \
  (void) (ML_POINTER)->nelem; \
  (void) (RL_POINTER)->dims;  \
  (void) (RL_POINTER)->nelem

#define LOCALDEF_FREE_LYAPUNOV_MEMORY_ON_ERROR( ERROR ) \
  if ( ERROR ){ \
  LOCALDEF_LYAPUNOV_CHANGE_MATSET_DIMENSIONS( &ML, &RL, dims );\
  MUTIL_FREE_WARN( memlist, &list ); \
  return ( ERROR ); \
  }


#define LOCALDEF_FREE_LYAPUNOV_MEMORY_ON_ERROR_II( ERROR ) \
  if ( ERROR ){ \
  LOCALDEF_LYAPUNOV_CHANGE_MATSET_DIMENSIONS( ML, RL, dims ); \
  MUTIL_FREE_WARN( memlist, &list ); \
  return ( ERROR ); \
  }

#define LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( ERROR ) \
  if ( ERROR ){ \
  (void) matset_free( &ML ); \
  MUTIL_FREE_WARN( memlist, &list ); \
  return ( ERROR ); \
  }

#undef LOCALDEF_CHECK_NULL_POINTER_LYAP
#define LOCALDEF_CHECK_NULL_POINTER_LYAP( DATA_PTR,         \
  DATA_TYPE, TYPE_PREFIX )                                  \
  err = TYPE_PREFIX ## _validate( DATA_PTR );               \
  if ( err ) return err;                                    \
  if ( DATA_PTR == ( DATA_TYPE * ) NULL ) {                 \
  MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
  return MUTIL_ERR_NULL_POINTER;                            \
  }


/*
****************************************
*                                      *
* STATIC (local) FUNCTION DECLARATIONS *
*                                      *
****************************************
*/

static mutil_errcode localfn_lyap_sf_compute(
  const univ_mat   *time_series,
  const mat_set    *ref_set,
  const univ_mat   *ref_indices,
  const univ_mat   *ref_nneig,
  const sint32      dimension,
  const sint32      delay,
  sint32            nfuture,
  univ_mat         *scaling );


static mutil_errcode localfn_lyap_sf_findref(
  const univ_mat              *embedding,
  const mutil_kdtree          *kdtree,
  const double                 epsilon,
  univ_mat                    *ref_indices,
  univ_mat                    *ref_nneig,
  const sint32                 nref_points,
  const sint32                 min_nneig,
  const fra_distance_metric    metric,
  const sint32                 orbital_lag,
  sint32                      *ref_count,
  mat_set                     *ref_set );


static mutil_errcode localfn_lyap_sf_memory(
  const sint32     nref_points,
  const sint32     nfuture,
  univ_mat        *ref_indices,
  univ_mat        *ref_nneig,
  univ_mat        *scaling,
  memlist         *mlist  );

static mutil_errcode localfn_permute_ref_indices( univ_mat *indices );


static mutil_errcode localfn_lyap_maximal_error_checks(
  const univ_mat             *time_series,
  const sint32                delay,
  const univ_mat             *dims,
  const univ_mat             *epsilons,
  const fra_distance_metric   metric,
  const sint32                nref_points,
  const sint32                min_nneig,
  const sint32                orbital_lag,
  const sint32                nfuture,
  const mat_set              *result );


static mutil_errcode localfn_lyap_sf_error_checks(
  const univ_mat             *time_series,
  const sint32                delay,
  const sint32                dim,
  const double                epsilon,
  const fra_distance_metric   metric,
  const sint32                nref_points,
  const sint32                min_nneig,
  const sint32                orbital_lag,
  const sint32                nfuture,
  sint32                     *embed_npoints,
  sint32                     *npoints,
  univ_mat                   *scaling );


static mutil_errcode localfn_least_squares_polynomial_design_matrix(
  const double_mat *z,
  const double_mat *zmax,
  const sint32      order,
  double_mat       *result );


static mutil_errcode localfn_fill_polyorder(
  sint32       current_order,
  sint32       current_index,
  sint32       max_order,
  sint32       embedding_dimension,
  double     **poly,
  double_mat  *product,
  double_mat  *scale,
  double      *pd_z,
  double      *pd_zmax );

static mutil_errcode localfn_qr_decomposition(
  const double_mat *A,
  void             *intrp_ptr,
  double_mat       *Q,
  double_mat       *R );

static mutil_errcode localfn_qr_solve(
  const double_mat *Q,
  const double_mat *R,
  const double_mat *B,
  double_mat       *X );

static mutil_errcode localfn_vector_dot_product(
  double *x,
  double *y,
  sint32 n,
  double *dotprod );

static mutil_errcode localfn_matrix_product_qr_decomposition(
  const mat_set *ML,
  void          *intrp_ptr,
  univ_mat      *Q,
  mat_set       *RL );

static mutil_errcode localfn_calculate_jacobian(
  const univ_mat *embedding,
  const sint32    global_reference_point,
  const sint32    n_reference,
  const sint32    local_dimension,
  const sint32    orbital_lag,
  const sint32    scale_max,
  const sint32    polyorder,
  const boolean   forward,
  const fra_distance_metric metric,
  void           *intrp_ptr,
  mat_set        *jacobian );

static mutil_errcode localfn_matset_transpose(
  const mat_set *x,
  mat_set       *transpose_x );

static void localfn_transpose_fast(
  const double_mat *x,
  double_mat       *xt );


/*
*********************************
*                               *
* LIBRARY FUNCTION DEFINITIONS  *
*                               *
*********************************
*/


/* Lyapunov scaling function curves */
/*                                  */
/* Documented in fra_lyap.h         */
/* Written by Keith L. Davidson     */
mutil_errcode frauniv_lyapunov_maximal(
  const univ_mat             *time_series,
  const sint32                delay,
  const univ_mat             *dims,
  const univ_mat             *epsilons,
  const fra_distance_metric   metric,
  const sint32                nref_points,
  const sint32                min_nneig,
  const sint32                orbital_lag,
  const sint32                nfuture,
  void                       *intrp_ptr,
  mat_set                    *result )
{
  mat_set          ref_set;

  univ_mat         ref_indices;
  univ_mat         ref_nneig;
  univ_mat         embedding;
  univ_mat         scaling;

  sint32_mat       nrow_mat;
  sint32_mat       ncol_mat;

  const sint32     two = (sint32) 2;

  double          *sc_data;

  sint32           ref_count;
  sint32           ndims;
  sint32           neps;
  sint32           offset;
  sint32           d;
  sint32           e;
  sint32           i;

  mutil_kdtree     kdtree;
  memlist          mlist;
  mutil_errcode    err;

  MUTIL_TRACE( "Start in frauniv_lyapunov_maximal()" );

  /* Avoid lint warning */
  (void) whatssi;


  /* Initialize memory management list */
  MEMLIST_INIT( mlist );


  /* Error checks */
  err = localfn_lyap_maximal_error_checks(
    time_series,
    delay,
    dims,
    epsilons,
    metric,
    nref_points,
    min_nneig,
    orbital_lag,
    nfuture,
    result );
  if ( err ) return err;


  /* Grab dimensions */
  ndims  = MATUNIV_NELEM( dims );
  neps   = MATUNIV_NELEM( epsilons );
  offset = ndims * neps;


  /* Allocate memory for output */
  err = matset_malloc_register( result, (sint32) 1, &two, &mlist );
  if ( err ) {
    MUTIL_ERROR( "Could not allocate memory for output matrix set" );
    return MUTIL_ERR_MEM_ALLOC;
  }

  err = mats32_malloc_register( &nrow_mat, two, (sint32) 1, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  err = mats32_malloc_register( &ncol_mat, two, (sint32) 1, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );

  nrow_mat.data[ 0 ] = nfuture + 1;
  nrow_mat.data[ 1 ] = nfuture + 1;
  ncol_mat.data[ 0 ] = ndims * neps;
  ncol_mat.data[ 1 ] = (sint32) 1;

  err = matset_malloc_matrices_arbitrary_size(
    result,
    &nrow_mat,
    &ncol_mat,
    MUTIL_DOUBLE );
  if ( err ) {
    MUTIL_ERROR( "Could not allocate memory for output matrix set matrices" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return MUTIL_ERR_MEM_ALLOC;
  }


  /* Initialize output matrix to all zeros */
  err = matdbl_assign_scalar(
    (double) 0,
    intrp_ptr,
    &(result->mats[ 0 ].mat.dblmat) );
  MEMLIST_FREE_ON_ERROR( err, &mlist );


  /* De-allocate row and column dimension matrices */
  err = memlist_member_free( (void*) &nrow_mat, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );
  err = memlist_member_free( (void*) &ncol_mat, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );


  /* Allocate "work" memory */
  err = localfn_lyap_sf_memory(
    nref_points,
    nfuture,
    &ref_indices,
    &ref_nneig,
    &scaling,
    &mlist  );
  if ( err ) {
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }

  /* Begin calculating the scaling function for each dimension and each
  epsilon */
  for ( d = 0; d < ndims; d++ ) {

    /* form a delay embedding */
    err = frauniv_embed(
      time_series,
      dims->mat.s32mat.data[ d ],
      delay,
      intrp_ptr,
      &embedding );
    if ( err ) {
      MUTIL_ERROR( "Could not create delay embedding" );
      MUTIL_FREE_WARN( memlist, &mlist );
      return err;
    }


    /* add embedding matrix to memory management list */
    err = memlist_member_register( &mlist, (void*) &embedding, MEMTYPE_MATUNIV );
    if ( err ) {
      MUTIL_ERROR( "Memory management error" );
      MUTIL_FREE_WARN( matuniv, &embedding );
      MUTIL_FREE_WARN( memlist, &mlist );
      return err;
    }


    /* create a kd-tree for neighbor searching */
    err = mutil_kdtree_malloc_register( &kdtree, &embedding, (sint32) 1, &mlist );
    if ( err ) {
      MUTIL_ERROR( "Could not create kd-tree for neighbor search" );
      MUTIL_FREE_WARN( memlist, &mlist );
      return err;
    }

    /* loop over all epsilons */
    for ( e = 0; e < neps; e++ ) {

      /* find reference points in the embedding */
      err = localfn_lyap_sf_findref(
        &embedding,
        &kdtree,
        epsilons->mat.dblmat.data[ e ],
        &ref_indices,
        &ref_nneig,
        nref_points,
        min_nneig,
        metric,
        orbital_lag,
        &ref_count,
        &ref_set );
      if ( err ) {
        MUTIL_FREE_WARN( memlist, &mlist );
        return err;
      }

      /* If no reference points were found then skip to the next epsilon */
      if ( ref_count == (sint32) 0 ) continue;

      /* Compute the scaling function */
      err = localfn_lyap_sf_compute(
        time_series,
        &ref_set,
        &ref_indices,
        &ref_nneig,
        dims->mat.s32mat.data[ d ],
        delay,
        nfuture,
        &scaling );
      if ( err ) {
        MUTIL_ERROR( "Error while computing scaling function" );
        MUTIL_FREE_WARN( matset, &ref_set );
        MUTIL_FREE_WARN( memlist, &mlist );
        return err;
      }

      /* Free the set of reference point neighbor indices */
      MUTIL_FREE_WARN( matset, &ref_set );

      /* Copy the scaling function into the first matrix
      of the output matrix set. Each column is for a given (dim,epsilon)
      pair, with epsilon as the fast running index */
      sc_data = result->mats[ 0 ].mat.dblmat.data + d * neps + e;
      for ( i = 0; i <= nfuture; i++, sc_data += offset )
        *sc_data = scaling.mat.dblmat.data[ i ];

    } /* for ( e = 0; ... ) */


    /* remove the kd-tree from memory management list */
    err = memlist_member_free( (void*) &kdtree, &mlist );
    if ( err ) {
      MUTIL_ERROR( "Memory management error" );
      MUTIL_FREE_WARN( memlist, &mlist );
      return err;
    }


    /* remove the embedding matrix from memory management list */
    err = memlist_member_free( (void*) &embedding, &mlist );
    if ( err ) {
      MUTIL_ERROR( "Memory management error" );
      MUTIL_FREE_WARN( memlist, &mlist );
      return err;
    }

  } /* for ( d = 0; ... ) */


  /* Remove the output matrix set from the memory management list */
  err = memlist_member_unregister( (void*) result, &mlist );
  MEMLIST_FREE_ON_ERROR( err, &mlist );


  /* Free all other memory */
  MUTIL_FREE_WARN( memlist, &mlist );


  /* Fill the scaling function argument values (abcissa) */
  for ( i = 0; i <= nfuture; i++ )
    result->mats[ 1 ].mat.dblmat.data[ i ] = (double) i;


  MUTIL_TRACE( "Done with frauniv_lyapunov_maximal()" );

  return MUTIL_ERR_OK;
}



/* Lyapunov scaling function.   */
/*                              */
/* Documented in fra_lyap.h     */
/* Written by Keith L. Davidson */
mutil_errcode frauniv_lyapunov_scaling_function(
  const univ_mat             *time_series,
  const sint32                delay,
  const sint32                dim,
  const double                epsilon,
  const fra_distance_metric   metric,
  const sint32                nref_points,
  const sint32                min_nneig,
  const sint32                orbital_lag,
  const sint32                nfuture,
  void                       *intrp_ptr,
  sint32                     *nref_found,
  univ_mat                   *ref_indices,
  univ_mat                   *ref_nneig,
  univ_mat                   *scaling )
{
  mat_set          ref_set;
  univ_mat         embedding;

  sint32           embed_npoints;
  sint32           npoints;
  sint32           ref_count;

  mutil_kdtree     kdtree;
  memlist          mlist;
  mutil_errcode    err;

  MUTIL_TRACE( "Start in frauniv_lyapunov_scaling()" );

  /* Initialize memory management list */
  MEMLIST_INIT( mlist );


  /* Error checks */
  err = localfn_lyap_sf_error_checks(
    time_series,
    delay,
    dim,
    epsilon,
    metric,
    nref_points,
    min_nneig,
    orbital_lag,
    nfuture,
    &embed_npoints,
    &npoints,
    scaling );
  if ( err ) return err;


  /* Form a delay embedding */
  err = frauniv_embed(
    time_series,
    dim,
    delay,
    intrp_ptr,
    &embedding );
  if ( err ) {
    MUTIL_ERROR( "Could not create delay embedding" );
    return err;
  }


  /* Add embedding matrix to memory management list */
  err = memlist_member_register( &mlist, (void*) &embedding, MEMTYPE_MATUNIV );
  if ( err ) {
    MUTIL_ERROR( "Memory management error" );
    MUTIL_FREE_WARN( matuniv, &embedding );
    return err;
  }


  /* Create a kd-tree for neighbor searching */
  err = mutil_kdtree_malloc_register( &kdtree, &embedding, (sint32) 1, &mlist);
  if ( err ) {
    MUTIL_ERROR( "Could not create kd-tree for neighbor search" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }


  /* Allocate necessary memory */
  err = localfn_lyap_sf_memory(
    nref_points,
    nfuture,
    ref_indices,
    ref_nneig,
    scaling,
    &mlist  );
  if ( err ) {
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }


  /* Find reference points in the embedding */
  err = localfn_lyap_sf_findref(
    &embedding,
    &kdtree,
    epsilon,
    ref_indices,
    ref_nneig,
    nref_points,
    min_nneig,
    metric,
    orbital_lag,
    &ref_count,
    &ref_set );
  if ( err ) {
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }


  /* Record the number of reference points found */
  *nref_found = ref_count;


  /* If no reference points were found then exit */
  if ( *nref_found == (sint32) 0 ) {
    MUTIL_FREE_WARN( memlist, &mlist );
    return MUTIL_ERR_OK;

  } else {

    /* add the matrix set of neighbor indices to the memory management list */
    err = memlist_member_register( &mlist, (void*) &ref_set, MEMTYPE_MATSET );
    if ( err ) {
      MUTIL_FREE_WARN( memlist, &mlist );
      return err;
    }

    /* if necessary, reallocate memory for the reference point indices and
    number of neighbors */
    if ( ref_count < nref_points ) {
      err = matuniv_realloc_register( ref_indices, ref_count, (sint32) 1, &mlist );
      if ( err ) {
        MUTIL_FREE_WARN( memlist, &mlist );
        return err;
      }
      err = matuniv_realloc_register( ref_nneig, ref_count, (sint32) 1, &mlist );
      if ( err ) {
        MUTIL_FREE_WARN( memlist, &mlist );
        return err;
      }
    }
  }


  /* Remove the kd-tree from memory management list */
  err = memlist_member_free( (void*) &kdtree, &mlist );
  if ( err ) {
    MUTIL_ERROR( "Memory management error" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }


  /* Remove the embedding matrix from memory management list */
  err = memlist_member_free( (void*) &embedding, &mlist );
  if ( err ) {
    MUTIL_ERROR( "Memory management error" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }


  /* Compute the scaling function */
  err = localfn_lyap_sf_compute(
    time_series,
    &ref_set,
    ref_indices,
    ref_nneig,
    dim,
    delay,
    nfuture,
    scaling );
  if ( err ) {
    MUTIL_ERROR( "Error while computing scaling function" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }


  /* Remove output matrices from the memory management list */
  err = memlist_member_unregister( (void*) scaling, &mlist );
  if ( err ) {
    MUTIL_ERROR( "Memory management error" );
    MUTIL_FREE_WARN( memlist, &mlist );
  }
  err = memlist_member_unregister( (void*) ref_indices, &mlist );
  if ( err ) {
    MUTIL_ERROR( "Memory management error" );
    MUTIL_FREE_WARN( matuniv, scaling );
    MUTIL_FREE_WARN( memlist, &mlist );
  }
  err = memlist_member_unregister( (void*) ref_nneig, &mlist );
  if ( err ) {
    MUTIL_ERROR( "Memory management error" );
    MUTIL_FREE_WARN( matuniv, scaling );
    MUTIL_FREE_WARN( matuniv, ref_indices );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }


  /* Free remaining memory */
  MUTIL_FREE_WARN( memlist, &mlist );

  MUTIL_TRACE( "Done with frauniv_lyapunov_scaling()" );

  return MUTIL_ERR_OK;
}






/*
*******************************
*                             *
* STATIC FUNCTION DEFINITIONS *
*                             *
*******************************
*/

static mutil_errcode localfn_lyap_sf_compute(
  const univ_mat   *time_series,
  const mat_set    *ref_set,
  const univ_mat   *ref_indices,
  const univ_mat   *ref_nneig,
  const sint32      dimension,
  const sint32      delay,
  sint32            nfuture,
  univ_mat         *scaling )
{
  register double    sum_i;
  register double    sum_r;
  register double   *ts_data;
  register double    tmp_dbl;
  register sint32    offset;
  register sint32    i;
  register sint32    idx1;
  register sint32    idx2;
  register sint32    n;
  register sint32    nelem;
  register sint32    nref;
  register sint32    r;


  /* Add 1 to the number of future time steps to span (for-loop and indexing
  convenience), we're computing delta-n from 0 to the user-specified
  nfuture */
  nfuture++;

  /* Set the indexing offset */
  offset = ( dimension - 1 ) * delay;

  /* Initialize data pointer and number of elements */
  ts_data = (double*) MATUNIV_DATA( time_series );
  nelem   = MATUNIV_NELEM( time_series );
  nref    = ref_set->nelem;

  /* Begin computing the scaling function */
  for ( n = 0; n < nfuture; n++ ) {

    for ( sum_r = 0.0, r = 0; r < nref; r++ ) {

      idx1 = ref_indices->mat.s32mat.data[ r ] + offset + n;

      if ( idx1 < nelem ) {

        tmp_dbl = ts_data[ idx1 ];

        for ( sum_i = 0.0, i = 0; i < ref_nneig->mat.s32mat.data[ r ]; i++ ) {

          idx2 = ref_set->mats[ r ].mat.s32mat.data[ i ] + offset + n;

          if ( idx2 < nelem ) sum_i += MUTIL_ABS( tmp_dbl - ts_data[ idx2 ] );

        } /* for ( sum_i = 0.0, ... ) */

        sum_r += FRA_LOCAL_LOG( sum_i / ref_nneig->mat.s32mat.data[ r ] );

      } /* if ( idx1 < nelem ) */

    } /* for ( sum_r = 0.0, ... ) */

    *( scaling->mat.dblmat.data + n ) = sum_r / (double) nref;

  } /* for ( n = 0; ... ) */

  return MUTIL_ERR_OK;
}



static mutil_errcode localfn_lyap_sf_findref(
  const univ_mat              *embedding,
  const mutil_kdtree          *kdtree,
  const double                 epsilon,
  univ_mat                    *ref_indices,
  univ_mat                    *ref_nneig,
  const sint32                 nref_points,
  const sint32                 min_nneig,
  const fra_distance_metric    metric,
  const sint32                 orbital_lag,
  sint32                      *ref_count,
  mat_set                     *ref_set )
{
  univ_mat       *matrix_ptr;
  univ_mat        orig_idx;
  univ_mat        neig_idx;
  univ_mat        neig_dist;
  univ_mat        point;
  univ_mat        perm_idx;
  sint32          dim;
  sint32          i;
  sint32          nembed;
  sint32          nperm;
  sint32          point_counter;

  mutil_errcode   err;

  MUTIL_TRACE( "Start in localfn_lyap_sf_findref()" );

  /* Initialize variables that are needed to search for reference points */
  *ref_count          = (sint32) 0;
  point_counter       = (sint32) 0;
  nembed              = MATUNIV_NROW( embedding );
  dim                 = MATUNIV_NCOL( embedding );
  ref_set->ndim       = (sint32) 1;
  ref_set->contiguous = FALSE;
  ref_set->nelem      = (sint32) 0;


  /* Allocate memory used to manually construct the matrix set that will hold
  the neighbor indices for each found reference point */
  err = mutil_malloc( sizeof( sint32 ), (void**) &(ref_set->dims) );
  if ( err ) return MUTIL_ERR_MEM_ALLOC;
  ref_set->dims[ 0 ] = (sint32) 0;

  err = mutil_malloc(
    nref_points * sizeof( univ_mat ),
    (void**) &(ref_set->mats) );
  if ( err ) {
    (void) mutil_free( (void*) ref_set->dims, sizeof( sint32 ) );
    return MUTIL_ERR_MEM_ALLOC;
  }

  point.type             = MUTIL_DOUBLE;
  point.mat.dblmat.nrow  = (sint32) 1;
  point.mat.dblmat.ncol  = dim;
  point.mat.dblmat.nelem = dim;

  (void) point.type;              /* supresses lint warning */
  (void) point.mat.dblmat.nrow;   /* supresses lint warning */
  (void) point.mat.dblmat.ncol;   /* supresses lint warning */
  (void) point.mat.dblmat.nelem;  /* supresses lint warning */


#define FRA_LYAP_FREE_REF_SET \
  \
  (void) mutil_free( (void*) ref_set->dims, sizeof( sint32 ) ); \
  for ( i = 0; i < *ref_count; i++ ) \
  MUTIL_FREE_WARN( matuniv, ref_set->mats + i ); \
  \
  (void) mutil_free( (void*) ref_set->mats, nref_points * sizeof( univ_mat ) )


  /* Allocate memory for and assign indices */
  nperm = nembed / orbital_lag; /* orbital_lag is always > 0 */
  err = matuniv_malloc(
    &perm_idx,
    nperm,
    (sint32) 1,
    MUTIL_SINT32 );
  if ( err ) {
    FRA_LYAP_FREE_REF_SET;
    return err;
  }
  for ( i = 0; i < nperm; i++ )
    perm_idx.mat.s32mat.data[ i ] = i * orbital_lag;

  /* Permute the indices of possible reference points */
  err = localfn_permute_ref_indices( &perm_idx );
  if ( err ) {
    FRA_LYAP_FREE_REF_SET;
    MUTIL_FREE_WARN( matuniv, &perm_idx );
    return err;
  }


  /* Find reference points, manually constructing a matrix set along the way */
  while ( *ref_count < nref_points && point_counter < nperm ) {

    /* set the point matrix data pointer */
    point.mat.dblmat.data =
      embedding->mat.dblmat.data +
      perm_idx.mat.s32mat.data[ point_counter ] * dim;

    /* find neighbors of the current point */
    err = frauniv_neighbor_find_arbitrary(
      &point,
      kdtree,
      (sint32) 0,
      epsilon,
      metric,
      (univ_mat*) NULL,
      (boolean) FALSE,
      orbital_lag,
      (univ_mat*) NULL,
      (void*) NULL,
      &orig_idx,
      &neig_idx,
      &neig_dist );
    if ( err == MUTIL_ERR_ZERO_NEIGHBORS_FOUND ) {

      point_counter++;
      continue;

    } else if ( err ) {

      FRA_LYAP_FREE_REF_SET;
      MUTIL_FREE_WARN( matuniv, &perm_idx );
      MUTIL_ERROR( "Error while searching for neighbors" );
      return err;
    }

    /* Free unneeded memory */
    MUTIL_FREE_WARN( matuniv, &orig_idx );
    MUTIL_FREE_WARN( matuniv, &neig_dist );


    /* if the current point has the required number of neighbors then grab
    the matrix holding the neighbor indices and record the point index */
    if ( neig_idx.mat.s32mat.nelem >= min_nneig ) {

      ref_indices->mat.s32mat.data[ *ref_count ] =
        perm_idx.mat.s32mat.data[ point_counter ];

      ref_nneig->mat.s32mat.data[ *ref_count ] = neig_idx.mat.s32mat.nelem;

      matrix_ptr = ref_set->mats + *ref_count;

      matrix_ptr->type             = MUTIL_SINT32;
      matrix_ptr->mat.s32mat.nrow  = neig_idx.mat.s32mat.nrow;
      matrix_ptr->mat.s32mat.ncol  = neig_idx.mat.s32mat.ncol;
      matrix_ptr->mat.s32mat.nelem = neig_idx.mat.s32mat.nelem;
      matrix_ptr->mat.s32mat.data  = neig_idx.mat.s32mat.data;

      ref_set->dims[ 0 ]++;
      ref_set->nelem++;
      ref_count[ 0 ]++;

    } else {

      MUTIL_FREE_WARN( matuniv, &neig_idx );
    }

    /* increment the point counter */
    point_counter++;
  }

  /* Free memory for permutation matrices */
  MUTIL_FREE_WARN( matuniv, &perm_idx );

  /* Inspect the number of found reference points. If no reference points */
  /* are found then just return. Otherwise, if the number of points found */
  /* is less than what the user specified then reallocate (truncation)    */
  /* memory appropriately.                                                */
  if ( *ref_count == (sint32) 0 ) {

    FRA_LYAP_FREE_REF_SET;
    return MUTIL_ERR_OK;

  } else if ( *ref_count < nref_points ) {

  /* truncate the memory pointed to by the matrix set internal universal
  matrix pointer if the number of reference points found is less than
    that requested by the user. */
    err = mutil_realloc(
      (void**) &(ref_set->mats),
      *ref_count * sizeof( univ_mat ),
      nref_points * sizeof( univ_mat ) );
    if ( err ) {
      FRA_LYAP_FREE_REF_SET;
      return err;
    }

  } /* if ( *ref_count == (sint32) 0 )  */

  MUTIL_TRACE( "Done with localfn_lyap_sf_findref()" );

  return MUTIL_ERR_OK;
}


static mutil_errcode localfn_lyap_sf_memory(
  const sint32     nref_points,
  const sint32     nfuture,
  univ_mat        *ref_indices,
  univ_mat        *ref_nneig,
  univ_mat        *scaling,
  memlist         *mlist  )
{
  mutil_errcode   err;

  MUTIL_TRACE( "Start in localfn_lyap_sf_memory()" );

  /* scaling function values */
  err = matuniv_malloc_register(
    scaling,
    nfuture + 1,
    (sint32) 1,
    MUTIL_DOUBLE,
    mlist );
  if ( err ) {
    MUTIL_ERROR( "Memory allocation error" );
    return MUTIL_ERR_MEM_ALLOC;
  }


  /* reference point indices */
  err = matuniv_malloc_register(
    ref_indices,
    nref_points,
    (sint32) 1,
    MUTIL_SINT32,
    mlist );
  if ( err ) {
    MUTIL_ERROR( "Memory allocation error" );
    return MUTIL_ERR_MEM_ALLOC;
  }


  /* reference point neighbor count */
  err = matuniv_malloc_register(
    ref_nneig,
    nref_points,
    (sint32) 1,
    MUTIL_SINT32,
    mlist );
  if ( err ) {
    MUTIL_ERROR( "Memory allocation error" );
    return MUTIL_ERR_MEM_ALLOC;
  }

  MUTIL_TRACE( "Done with localfn_lyap_sf_memory()" );

  return MUTIL_ERR_OK;
}


static mutil_errcode localfn_permute_ref_indices( univ_mat *indices )
{
  register sint32     i;
  register sint32     j;
  register sint32     tmp;
  register sint32     n    = MATUNIV_NELEM( indices );
  register sint32    *data = (sint32*) MATUNIV_DATA( indices );

  double              tmp_rand;

  void               *rand_ptr;
  mutil_errcode       err;

  MUTIL_TRACE( "Start in localfn_permute_ref_indices()" );


  /* Initialize the random number generator */
  err = mutil_rand_begin( &rand_ptr );
  if ( err ) {
    MUTIL_ERROR( "Could not initialize the random number generator" );
    return err;
  }
  err = mutil_rand_set_seed( (void*) NULL, rand_ptr );
  if ( err ) {
    MUTIL_ERROR( "Could not set the random number generator seed" );
    return err;
  }


  /* Begin permutation */
  for ( i = n - 1; i > 0; i-- ) {

    /* generate uniformly-distributed random number */
    err = mutil_rand_uniform( rand_ptr, &tmp_rand );
    if ( err ) {
      MUTIL_ERROR( "Could not generate a uniformly-distributed "
        "random number" );
      return err;
    }

    /* compute random index */
    j = (sint32) floor( i * tmp_rand );

    /* swap values */
    tmp       = data[ i ];
    data[ i ] = data[ j ];
    data[ j ] = tmp;
  }

  MUTIL_TRACE( "Done with localfn_permute_ref_indices()" );

  return MUTIL_ERR_OK;
}



static mutil_errcode localfn_lyap_maximal_error_checks(
  const univ_mat             *time_series,
  const sint32                delay,
  const univ_mat             *dims,
  const univ_mat             *epsilons,
  const fra_distance_metric   metric,
  const sint32                nref_points,
  const sint32                min_nneig,
  const sint32                orbital_lag,
  const sint32                nfuture,
  const mat_set              *result )
{
  sint32          npoints;
  sint32          embed_npoints;
  sint32          i;
  boolean         is_delay_embedding;
  mutil_errcode   err;

  MUTIL_TRACE( "Start in localfn_lyap_maximal_error_checks()" );

  /* time_series */
  if ( time_series == (univ_mat*) NULL ) {
    MUTIL_ERROR( "Input time_series is NULL pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }
  if ( time_series->type != MUTIL_DOUBLE ) {
    MUTIL_ERROR( "Input time_series must be of type MUTIL_DOUBLE" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  npoints = MATUNIV_NELEM( time_series );
  if ( npoints != MATUNIV_NROW( time_series ) &&
    npoints != MATUNIV_NCOL( time_series ) ) {
    MUTIL_ERROR( "Input time_series must be a row or column vector" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* delay */
  if ( delay < (sint32) 1 ) {
    MUTIL_ERROR( "Input delay must be a positive integer" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* nref_points */
  if ( nref_points < (sint32) 1 ) {
    MUTIL_ERROR( "Input nref_points must be a positive integer" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* min_nneig */
  if ( min_nneig < (sint32) 1 ) {
    MUTIL_ERROR( "Input min_nneig must be a positive integer" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* orbital_lag */
  if ( orbital_lag < (sint32) 1 ) {
    MUTIL_ERROR( "Input orbital_lag must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* dims */
  npoints = MATUNIV_NELEM( dims );
  if ( npoints != MATUNIV_NROW( dims ) &&
    npoints != MATUNIV_NCOL( dims ) ) {
    MUTIL_ERROR( "Input dims must be a row or column vector" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  if ( dims->type != MUTIL_SINT32 ) {
    MUTIL_ERROR( "Input dims must be of type MUTIL_SINT32" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  for ( i = 0; i < npoints; i++ ) {

    /* make sure all dimensions are greater than 1 */
    if ( dims->mat.s32mat.data[ i ] < 2 ) {
      MUTIL_ERROR( "Input dims must elements greater than 1" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

    /* make sure the input values allow an embedding */
    err = frauniv_check_embedding_inputs(
      time_series,
      dims->mat.s32mat.data[ i ],
      delay,
      orbital_lag,
      (void*) NULL,
      &is_delay_embedding,
      &embed_npoints );
    if ( err ) return err;

    /* nref_points (again) */
    if ( nref_points > embed_npoints / orbital_lag ) {
      MUTIL_ERROR( "Input nref_points is too large for the number of  "
        "points in one of the delay embeddings, and the given orbital lag" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

    /* min_nneig (again) */
    if ( min_nneig > ( embed_npoints - 1 ) ) {
      MUTIL_ERROR( "Input min_nneig is too large for the number of  "
        "points in one of the delay embeddings" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

    /* orbital_lag (again) */
    if ( orbital_lag > ( embed_npoints - 1 ) ) {
      MUTIL_ERROR( "Input min_nneig is too large for the number of  "
        "points in one of the delay embeddings" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
  }

  /* epsilons */
  npoints = MATUNIV_NELEM( epsilons );
  if ( epsilons->type != MUTIL_DOUBLE ) {
    MUTIL_ERROR( "Input epsilons must be of type MUTIL_DOUBLE" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  for ( i = 0; i < npoints; i++ ) {
    if ( epsilons->mat.dblmat.data[ i ] <= (double) 0 ) {
      MUTIL_ERROR( "Input epsilons must be a positive elements" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
  }

  /* metric */
  switch ( metric ) {
  case FRA_DISTANCE_L1:
  case FRA_DISTANCE_L2:
  case FRA_DISTANCE_LINFINITY:
    break;
  default:
    MUTIL_ERROR( "Input metric is invalid" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* nfuture */
  if ( nfuture < (sint32) 0 ) {
    MUTIL_ERROR( "Input nfuture must be non-negative" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* result (output matrix) */
  if ( result == (mat_set*) NULL ) {
    MUTIL_ERROR( "Output argument result is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  MUTIL_TRACE( "Done with localfn_lyap_maximal_error_checks()" );

  return MUTIL_ERR_OK;
}





static mutil_errcode localfn_lyap_sf_error_checks(
  const univ_mat             *time_series,
  const sint32                delay,
  const sint32                dim,
  const double                epsilon,
  const fra_distance_metric   metric,
  const sint32                nref_points,
  const sint32                min_nneig,
  const sint32                orbital_lag,
  const sint32                nfuture,
  sint32                     *embed_npoints,
  sint32                     *npoints,
  univ_mat                   *scaling )
{
  mutil_errcode   err;
  boolean         is_delay_embedding;

  MUTIL_TRACE( "Start in localfn_lyap_sf_error_checks()" );

  /* time_series */
  if ( time_series == (univ_mat*) NULL ) {
    MUTIL_ERROR( "Input time_series is NULL pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }
  if ( time_series->type != MUTIL_DOUBLE ) {
    MUTIL_ERROR( "Input time_series must be of type MUTIL_DOUBLE" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  *npoints = MATUNIV_NELEM( time_series );
  if ( *npoints != MATUNIV_NROW( time_series ) &&
    *npoints != MATUNIV_NCOL( time_series ) ) {
    MUTIL_ERROR( "Input time_series must be a row or column vector" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* delay */
  if ( delay < (sint32) 1 ) {
    MUTIL_ERROR( "Input delay must be a positive integer" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* dim */
  if ( dim < 2 ) {
    MUTIL_ERROR( "Input dim must be greater than 1" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* make sure the input values allow an embedding */
  err = frauniv_check_embedding_inputs( time_series, dim, delay, orbital_lag, (void *) NULL,
    &is_delay_embedding, embed_npoints );
  if ( err ) return err;

  /* epsilon */
  if ( epsilon <= (double) 0 ) {
    MUTIL_ERROR( "Input max_dist must be a positive number" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* metric */
  switch ( metric ) {
  case FRA_DISTANCE_L1:
  case FRA_DISTANCE_L2:
  case FRA_DISTANCE_LINFINITY:
    break;
  default:
    MUTIL_ERROR( "Input metric is invalid" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* nref_points */
  if ( nref_points < (sint32) 1 ) {
    MUTIL_ERROR( "Input nref_points must be a positive integer" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  } else if ( nref_points > *embed_npoints ) {
    MUTIL_ERROR( "Input nref_points is greater than the number of "
      "points in the phase space embedding" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* min_nneig */
  if ( min_nneig < (sint32) 1 ) {
    MUTIL_ERROR( "Input min_nneig must be a positive integer" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  } else if ( min_nneig > *embed_npoints ) {
    MUTIL_ERROR( "Input min_nneig too large, exceeds number "
		 "of embedded points" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* orbital_lag */
  if ( orbital_lag < (sint32) 1 ) {
    MUTIL_ERROR( "Input orbital_lag must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  } else if ( orbital_lag > (*embed_npoints - 1) ) {
    MUTIL_ERROR( "Input orbital_lag too large, exceeds number "
		 "of lags in the embedded vector time series" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* nfuture */
  if ( nfuture < (sint32) 0 ) {
    MUTIL_ERROR( "Input nfuture must be non-negative" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* scaling (output matrix) */
  if ( scaling == (univ_mat*) NULL ) {
    MUTIL_ERROR( "Output argument scalar is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  MUTIL_TRACE( "Done with localfn_lyap_sf_error_checks()" );

  return MUTIL_ERR_OK;
}



/* Local Lyapunov spectrum.       */
/*                                */
/* Documented in fra_lyap.h       */
/* Written by William Constantine */

mutil_errcode frauniv_local_lyapunov_spectrum(
  const univ_mat *time_series,
  const sint32    embedding_dimension,
  const sint32    time_lag,
  const sint32    orbital_lag,
  const double    sampling_interval,
  const sint32    local_dimension,
  const sint32    polynomial_order,
  const univ_mat *global_reference,
  const sint32    n_reference_local,
  const fra_distance_metric metric,
  const univ_mat *scale,
  void           *intrp_ptr,
  mat_set        *result )
{
  boolean        is_delay_embedding;
  double         offdiag;
  mat_set        ML;
  mat_set        RL;
  mat_set        jacobian;
  mat_set        jacobian_transpose;
  memlist        list;
  mutil_errcode  err;
  sint32         L2;
  sint32         L;
  sint32         dims;
  sint32         index;
  sint32         iscale;
  sint32         n;
  sint32         n_global_reference_point;
  sint32         max_reference_index;
  sint32         i;
  sint32         j;
  sint32         k;
  sint32         n_scale;
  sint32         n_embed;
  sint32         iref;
  sint32         scale_max;
  sint32         scale_min;
  sint32         global_reference_max;
  sint32         global_reference_min;
  univ_mat       Q;
  univ_mat       T;
  univ_mat       embedding;

  double sum;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_embed()" );

  /* avoid lint warning */

  ( void ) whatssi;

  MUTIL_TRACE( "Start frauniv_local_lyapunov_spectrum()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  err = frauniv_check_embedding_inputs(
    time_series, embedding_dimension, time_lag, orbital_lag,
    intrp_ptr, &is_delay_embedding, &n_embed );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( !MATANY_IS_VEC( &(time_series->mat.dblmat) ) ){
    MUTIL_ERROR( "Time series matrix must contain a single-column or single-row" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( sampling_interval <= 0.0 ){
    MUTIL_ERROR( "Sampling interval must be a positive value" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( local_dimension < 1 || local_dimension > embedding_dimension ){
    MUTIL_ERROR( "Local dimension must be at least one and not greater than the embedding dimension" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  if ( polynomial_order < 1 ){
    MUTIL_ERROR( "Polynomial order must be at least unity" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  switch( metric ){
    case FRA_DISTANCE_L1:
    case FRA_DISTANCE_L2:
    case FRA_DISTANCE_LINFINITY:
      break;
    default:
      MUTIL_ERROR( "Unsupported distance metric" );
      return MUTIL_ERR_ILLEGAL_VALUE;
  }

  LOCALDEF_CHECK_NULL_POINTER_LYAP( global_reference, univ_mat, matuniv );

  if ( global_reference->type != MUTIL_SINT32 ){
    MUTIL_ERROR( "Global_Reference matrix must be of type MUTIL_SINT32" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( !MATANY_IS_VEC( &(global_reference->mat.s32mat) ) ){
    MUTIL_ERROR( "Global_Reference matrix must contain a single-column or single-row" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  err = mats32_range( &( global_reference->mat.s32mat ), intrp_ptr, 
	  &global_reference_min, &global_reference_max );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if (global_reference_min < 0 || global_reference_max < 0){
	  MUTIL_ERROR( "Global reference indices must be non-negative" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( n_reference_local < 1 ){
	  MUTIL_ERROR( "Number of local reference points must exceed unity" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  LOCALDEF_CHECK_NULL_POINTER_LYAP( scale, univ_mat, matuniv );

  if ( !MATANY_IS_VEC( &(scale->mat.s32mat) ) ){
    MUTIL_ERROR( "Scale matrix must contain a single-column or single-row" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  err = mats32_range( &( scale->mat.s32mat ), intrp_ptr, &scale_min, &scale_max );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( scale_min < 1 ){
    MUTIL_ERROR( "Minimum scale must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( n_reference_local > n_embed - 2 - scale_max - global_reference_max ){
	  MUTIL_ERROR( "Number of local reference points too large" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( n_embed  < n_reference_local + 2 ){
    MUTIL_ERROR( "Number of embedding points is too large. " \
      "Decrease the number of reference points or increase the number " \
      "of points in the embedding" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* the number of Jacobians that we will form is Nj = n_reference_local + scale_max.
     given n_embed points, we have to make sure that the starting reference index
     is positive and sufficiently small so that we don't exceed the number of
     points in the embedding when forming the Jacobians (we also have to subtract
     2 from Nj to accommodate the shortened embeddding). */

  max_reference_index = n_embed - n_reference_local - scale_max - 2;

  if ( max_reference_index < 1 ){
	  MUTIL_ERROR( "Number of local reference points too large" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  for ( iref = 0; iref < MATUNIV_NELEM( global_reference ); iref++ ){

    if ( global_reference->mat.s32mat.data[iref] > max_reference_index ){
      MUTIL_ERROR( "Initial reference point is too large." );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

    if ( global_reference->mat.s32mat.data[iref] < 0 ){
      MUTIL_ERROR( "Starting reference index must be non-negative." );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
  }

  dims        = scale_max * 2;
  n_scale     = MATUNIV_NELEM( scale );
  n_global_reference_point = MATUNIV_NELEM( global_reference );

  /* allocate memory ... */

  err = matuniv_malloc_register( &Q, local_dimension,
    local_dimension, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &T, local_dimension,
    local_dimension, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* ... for RL we need to allocate enough space
     for twice the maximum scale. RL is shorthand
     for the R matrix of a QR decomposition of L matrix products.
  */

  err = matset_malloc_register( &RL, 1, &dims, &list );
  LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );

  err = matset_malloc_matrices( &RL, local_dimension,
    local_dimension, MUTIL_DOUBLE );
  LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );

  /* ... only need to allocate the space for the
     ML matrix set headers and not the data since
     we simply point the ML matrices to existing
     Jacobian (transform) matrices. ML is shorthand
     for L matrix products
  */

  err = matset_malloc( &ML, 1, &dims );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* ... the local Lyapunov exponents */

  err = matset_malloc_register( result, 1, &local_dimension, &list );
  LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );

  err = matset_malloc_matrices( result, n_global_reference_point,
    n_scale, MUTIL_DOUBLE );
  LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );

  /* create a time delayed embedding and register
     result with the memory manager */

  err = frauniv_embed( time_series, embedding_dimension,
    time_lag, intrp_ptr, &embedding );
  LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );

  err = memlist_member_register( &list, &embedding, MEMTYPE_MATUNIV );
  LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );

  /* zero out the result matrix set */

  for ( i = 0; i < local_dimension; i++ ){

    err = matdbl_assign_scalar( 0.0, intrp_ptr, &( result->mats[i].mat.dblmat ) );
    LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );
  }

  /* calculate the Lyapunov spectrum, once for each starting
     point (global_reference_point) */

  for ( n = 0; n < n_global_reference_point; n++ ){

    /* free previously allocated and registered memory */

    if ( memlist_member_exist( &jacobian, &list ) ){
      err = memlist_member_free( &jacobian, &list );
      LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );
    }
    if ( memlist_member_exist( &jacobian_transpose, &list ) ){
      err = memlist_member_free( &jacobian_transpose, &list );
      LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );
    }

    /* estimate the Jacobians of the embedding for forward exponents.
       we will obtain n_reference_local + scale_max matrices. register the
       Jacobian matrix set with the memory manager
    */

    err = localfn_calculate_jacobian( &embedding,
      global_reference->mat.s32mat.data[n], n_reference_local,
      local_dimension, orbital_lag, scale_max, polynomial_order,
      TRUE, metric, intrp_ptr, &jacobian );
    LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );

    err = memlist_member_register( &list, &jacobian, MEMTYPE_MATSET );
    LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );

    /* calculate the transpose of the Jacobians and register with
       the memory manager */

    err = localfn_matset_transpose( &jacobian, &jacobian_transpose );
    LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );

    err = memlist_member_register( &list, &jacobian_transpose, MEMTYPE_MATSET );
    LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );

    /* calculate the local Lyapunov exponents. keep in mind that the
       first matrix of the Jacobian matrix set is relative to the
       starting reference input established by the user, e.g., if
       reference = 3, the first Jacobian is estimated about the 4th
       point in the embedding (reference is in index base 0). */

    for ( iref = 0; iref < n_reference_local; iref++ ){

      for ( iscale = 0; iscale < n_scale; iscale++ ){

	L  = scale->mat.s32mat.data[iscale];
	L2 = L * 2;

	/* trick MUTILS into thinking that the ML and RL matrix sets
	   only contain L2 universal matrices */

	LOCALDEF_LYAPUNOV_CHANGE_MATSET_DIMENSIONS( &ML, &RL, L2 );

	/* set data pointers to the appropriate locations
	   in the input Jacobian matrix sets. here we form
	   a concatenation of the forward Jacobians and the
	   reverse of the transpose of the Jacobians so that
	   if L = 3 and iref = 2 then we would have

	   ML = { J[2], J[3], J[4], J'[4], J'[3], J'[2] }

	   where J' is the transpose of the Jacbian J matrix.
	*/

	index = iref;
	/*       index = iref + L - 1; */

	for ( j = 0; j < L2; j++ ){

	  ML.mats[j] = ( j < L ) ? jacobian.mats[ index++ ] : jacobian_transpose.mats[ --index ];
	  /*         ML.mats[j] = ( j < L ) ? jacobian.mats[ index-- ] : jacobian_transpose.mats[ ++index ]; */
	}

	/* perform matrix product QR decomposition */

	err = localfn_matrix_product_qr_decomposition( &ML, intrp_ptr, &Q, &RL );
	LOCALDEF_FREE_LYAPUNOV_MEMORY_ON_ERROR( err );

	/* refine the Q matrix until it is almost exactly the identity matrix.
	   this is done by testing the sum of the off-diagonal elements against
	   a threshold. If Q is an identity matrix, then the diagonal elements
	   of the product of the RL[i] matrices are the eigenvalues. do this refinement
	   at least one time, hence the 'do' loop */

	do{

	  /* make RL[0] = RL[0] * Q */

	  err = matuniv_multiply( &( RL.mats[0] ), &Q, intrp_ptr, &T );
	  LOCALDEF_FREE_LYAPUNOV_MEMORY_ON_ERROR( err );

	  err = matuniv_assign( &T, intrp_ptr, &( RL.mats[0] ) );
	  LOCALDEF_FREE_LYAPUNOV_MEMORY_ON_ERROR( err );

	  /* set ML[j] = RL[j] */

	  for ( j = 0; j < L2; j++ ){

	    ML.mats[j] = RL.mats[j];
	  }

	  /* perform matrix product QR decomposition (again) */

	  err = localfn_matrix_product_qr_decomposition( &ML, intrp_ptr, &Q, &RL );
	  LOCALDEF_FREE_LYAPUNOV_MEMORY_ON_ERROR( err );

	  /* calculate the Q matrix off-diagonal sum. if it is sufficiently
	     small, then end 'do' loop */

	  offdiag = 0.0;

	  for ( j = 0; j < local_dimension - 1; j++ ){

	    offdiag += (double) MUTIL_ABS(
	      Q.mat.dblmat.data[ j * ( local_dimension + 1 ) + 1 ] );
	  }

	} while ( offdiag > 1.0e-5 );

	/* the Q matrix is now almost the identity matrix and thus the
	   RL[i] eigenvalues are simply the diagonal elements of
	   RL[i]. It turns out that the if you multiply upper triangular
	   matrices, then the eigenvalues of that product are also
	   the product of the eigenvalues for each matrix. since these
	   eigenvalues are simply the diagonal elements, we can use the
	   diagonal elements of the RL[i] to estimate the local lyapunov
	   exponents:

	   we can sum the log of the (j,j) elements of each RL matrix.
	   normalizing this sum by the product of the sampling
	   interval and the current scale yields an estimate of the
	   local Lyapunov exponent. equivalently, we can take the
	   log of the product of these compnents so that we only
	   use the log() function once to save computations. */

	for ( j = 0; j < local_dimension; j++ ){

	  sum = 0.0;

	  index = j * ( local_dimension + 1 );

	  for ( k = 0; k < L2; k++ ){

	    sum += log( MUTIL_ABS( RL.mats[k].mat.dblmat.data[index] ) );
	  }

	  /* The exponent matrix set contains local_dimension universal matrices
	     of type MUTIL_DOUBLE, (one local Lypaunov set per dimension). Each
	     matrix is an n_scale vector, where n_scale is the number
	     of scales being investigated. The local Lyapunov exponents at all
	     scales are calculated below and added to the running total in the
	     corresponding exponent matrix. Later, these sums will be divided
	     by the number of reference points to form an average */

	  result->mats[j].mat.dblmat.data[ n * n_scale + iscale ] +=
	    sum / ( sampling_interval * (double) ( L2 * n_reference_local ) );
	}

      } /* end loop over each scale */

      /* check for interrupts */

      if ( MUTIL_INTERRUPT( 3.0 * n_reference_local, intrp_ptr ) ) {

	MUTIL_ERROR( "user interrupt" );

	/* restore original ML and RL matrix set dimensions */

	LOCALDEF_LYAPUNOV_CHANGE_MATSET_DIMENSIONS( &ML, &RL, dims );

	/* free unregisterd local memory */

	(void) matset_free( &ML );

	MUTIL_FREE_WARN( memlist, &list );
	return MUTIL_ERR_INTERRUPT;
      }

    } /* end loop over each reference point */

  } /* end loop over each starting point (global_reference_point) */

  /* restore original ML and RL matrix set dimensions */

  LOCALDEF_LYAPUNOV_CHANGE_MATSET_DIMENSIONS( &ML, &RL, dims );

  /* unregister memory associated with output */

  err = memlist_member_unregister( result, &list );
  LOCALDEF_FREE_LYAPUNOV_SPECTRUM_MEMORY_ON_ERROR( err );

  /* free unregisterd local memory */

  (void) matset_free( &ML );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with frauniv_local_lyapunov_spectrum()" );

  return MUTIL_ERR_OK;
}



/** Calculate the number of parameters needed to fit
 * a polynomial of a given order to data embedded in
 * a given dimension.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = localfn_number_polynomial_parameters( dimension, order );#
 * @return Standard mutils error code.
 * @param dimension The embedding dimension.
 * @param order Order of the polynomial to use in fitting the data.
 * @private
 */
static sint32 localfn_number_polynomial_parameters(
  const sint32 dimension,
  const sint32 order )
{
  sint32 den;
  sint32 i;
  sint32 num;

  num = den = 1;

  for ( i = 1; i <= order; i++ ){
    num *= ( dimension + i );
    den *= i;
  }

  return ((num - den) / den);
}

/** Calculate the design matrix for a multidimensional
 * least squares regression.
 * Given the model:
 * \[ z(n + 1) =
 *   sum_{i=1}^E b_i*z_i +
 *   sum_{i=1}^E sum_{j=i}^E c_{ij}*z_i*z_j +
 *   sum_{i=1}^E sum_{j=i}^E sum_{k=j}^E d_{ijk}*z_i*z_j*z_k + ...  \]
 *
 * where the z\_i for i = 1 ... E are the coordinates we wish to
 * use to fit the observed values y. In our case, the z\_i are defined
 * (in another function) as \[ z_i(n) = y^{(r)}_i(n) -  y^{(0)}_i(n) \]
 * so that the z vector represents the difference between between the
 * rth neighbor of some reference point y{0} in the phase space at time location
 * n. There will be a collection of such neighbors, and the above model
 * is used to form the design matrix, where each row is the model for
 * a particular neigbor and the columns correspond to the multidimensional
 * polynomial fit variables
 *
 * [ z1 z2 ... zE z1*z1 z1*z2 ...z1*zE ... zE*zE ... ]
 *
 * with P total terms (the underscores are removed for readability).
 * Given an N x E matrix
 * of neightbors z, where N is the number of neighbors and E is the
 * embedding dimension, the design matrix will be an N x P matrix
 * where P is the total number of terms in the regression function above.
 * If the polynomial order is 1, then only linear terms are used and
 * the regression model is [ z1 z2 ... zE ], thus in this case there are
 * only E terms. If the polynomial order is 2, then the model is based
 * on linear and quadratic terms and the total number of terms will
 * be P = E + E(E+1)/2. In general, the total number of terms up to
 * model order p for data embedded in dimension E will be
 * P = (p + E)!/(p!E!) - 1.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = localfn_least_squares_polynomial_design_matrix( &z, &zmax, order, intrpptr, &result )#
 * @return Standard mutils error code.
 * @param z An N x E double matrix containing
 *          the difference vectors (z) in a local neighborhood around
 *          some reference point in the phase space (See above for details).
 * @param zmax A 1 x E double amtrix containing
 *             the maximum value for each coordinate (dimension) of the z matrix.
 * @param order Order of the polynomial to use in fitting the data. 5 is the
 *              maximum order allowed.
 * @param result      Pointer to a double matrix
 *                    containing the design matrix. The memory for this matrix
 *                    is automatically allocated within the function.
 * @private
 */
static mutil_errcode localfn_least_squares_polynomial_design_matrix(
  const double_mat *z,
  const double_mat *zmax,
  const sint32      order,
  double_mat       *result )
{
  double         den;
  double         num;
  double        *pd_result;
  double        *pd_z;
  double        *pd_zmax;
  double       **poly;
  double_mat     product;
  double_mat     scale;
  memlist        list;
  mutil_errcode  err;
  sint32         embedding_dimension;
  sint32         i;
  sint32         n;
  sint32         n_neighbor;
  sint32         n_polyterm;
  sint32_mat     poly_offset;

  MUTIL_TRACE( "Starting localfn_least_squares_polynomial_design_matrix()" );

  /* initialize memory management list */

  MEMLIST_INIT( list );

  /* perform small check */

  if ( order < 1 || order > 5 ){

    MUTIL_ERROR( "Polynomial order must be on interval 1,...,5." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* initialize variables */

  n_neighbor = z->nrow;
  embedding_dimension = z->ncol;

  /* calculate the number of terms for each
     order of the polynomial. simlutaneously,
     create and initialize an array of double
     pointers which will be used to point to
     polynomial terms of like order within
     the result matrix. The total number of
     terms for a given polynomial order p and
     dimension E is given by

     N(p,E) = (p + E)! / ( p!E!) - 1

     The -1 part is there because we aren't using
     the zeroth order polynomial term in our model,
     mainly because ultimately we will use these
     polynomials to calculate the Jacobians, making
     them superfluous */

  err = mutil_malloc_register( order * sizeof(double *),
    (void **) &poly, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &poly_offset, 1, order, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  num = den = (sint32) 1;

  poly_offset.data[0] = 0;

  for ( i = 1; i <= order; i++ ){

    /* (E + 1)*(E + 2)*...*(E + i) */

    num *= ( embedding_dimension + i );

    /* factorial(i) */

    den *= i;

    /* total number of terms up to current order */

    n_polyterm = (sint32) ( num / den ) - 1;

    if ( i < order ){

      poly_offset.data[i] = n_polyterm;
    }
  }

  /* allocate memory and register with the memory manager */

  err = matdbl_malloc_register( result, n_neighbor, n_polyterm, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &scale, 1, order, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &product, 1, order, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* initialize pointers */

  pd_z      = z->data;
  pd_zmax   = zmax->data;
  pd_result = result->data;

  for ( n = 0; n < n_neighbor; n++ ){

    /* initialize polynomial pointers:

    poly[0] points to the beginning of the
    first order terms corresponding to the
    nth neighbor. poly_offset.data[i] contains
    contains the offset from poly[0] needed
    to point to the beginning of the ith order
    polynomial terms. */

    poly[0] = pd_result + ( n * n_polyterm );

    for ( i = 1; i < order; i++ ){

      poly[i] = poly[0] + poly_offset.data[i];
    }

    /* now calculate the design matrix */

    err = localfn_fill_polyorder(
      1,
      0,
      order,
      embedding_dimension,
      poly,
      &product,
      &scale,
      pd_z,
      pd_zmax );

    /* update pointer pd_z to point to next point (row)
       of z matrix */

    pd_z += embedding_dimension;

  } /* end loop over each neighbor */

  /* unregister memory associated with output */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_least_squares_polynomial_design_matrix()" );

  return MUTIL_ERR_OK;
}


/** Recursive function used to calculate all the terms of
 * a multidimensional polynomial fit (sans the zeroth order term).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = localfn_fill_polyorder( 1, 0, order, embedding_dimension, poly, &product, &scale, pd_z, pd_zmax )#
 * @param  current_order The current polynomial order being evaluated.
 * @param  current_index The current loop index (j on 1, ..., embedding\_dimension)
 *                       correpsonding to the current polynomial order.
 * @param  max_order     The maximum polynomial order to evaluate. Serves as a
 *                       stopping criteria for the recursive call structure and
 *                       is equivalent to the order parameter in the
 *                       \Ref{localfn_least_squares_polynomial_design_matrix} caller.
 * @param  embedding_dimension The (embedding) dimension of the data to be fit.
 * @param  poly  A double pointer to a double pointer.
 *   poly[p] initially points to the beginning of the
 *   pth order terms in the caller's result matrix. These
 *   are updated (incrememted) as each term is calculated.
 * @param  product A single-row double matrix containing max\_order
 *   elements. The first element will simply contain the z\_1 variable.
 *   The second element will initiially contain the z\_1*z\_1 product, and
 *   in general, the pth element will initially contain z\_1 raised to the pth power.
 *   These terms are updated as needed. For example, when the z\_1*z\_1*z\_3
 *   is being calculated, the product[1] (which contains second order products)
 *   will contain the product z\_1*z\_1. This will be the same for calculating
 *   all z\_1*z\_1*z\_j for j = 1, ..., embedding\_dimension and thus saves
 *   us from calculating the parent product redundantly.
 * @param  scale A single-row double matrix containing max\_order
 *   elements. The pth element contains the scale of the current
 *   order p-1 terms being evaluated. For example, if z\_1*z\_1*z\_3
 *   is being evaluated, the scale[2] element will be set to the maximum
 *   of zmax[0] and zmax[2] (corresponding to absolute maximum values of
 *   coordinates z\_1 and z\_3 taken over neighbors defined by the z matrix
 *   of the caller.
 * @param  pd_z Pointer to a double. This pointer points to elements of the z
 *   matrix in the caller.
 * @param  pd_zmax Pointer to a double. This pointer points to elements of the
 *   zmax vector in the caller.
 * @return Standard mutils error code.
 * @private
 */
 static mutil_errcode localfn_fill_polyorder(
   sint32       current_order,
   sint32       current_index,
   sint32       max_order,
   sint32       embedding_dimension,
   double     **poly,
   double_mat  *product,
   double_mat  *scale,
   double      *pd_z,
   double      *pd_zmax )
{
   mutil_errcode err;
   sint32        j;
   sint32        p = current_order - 1;

   /* input argument checks */

   if ( current_order < 1 ){

     MUTIL_ERROR( "Polynomial order must exceed unity" );
     return MUTIL_ERR_ILLEGAL_VALUE;
   }

   if ( current_index > embedding_dimension ){

     MUTIL_ERROR( "Current polynomial loop index cannot exceed the embedding dimension" );
     return MUTIL_ERR_ILLEGAL_VALUE;
   }

   if ( current_order <= max_order ){

     for ( j = current_index; j < embedding_dimension; j++ ){

       if ( current_order > 1 ){

         /* obtain the scale of the data using L-infinity norm */

         LOCALDEF_SCALE_MAX( scale->data[ p - 1 ], *( pd_zmax + j ), scale->data[ p ] );

         /* calculate the product for the current term.

          e.g., if we are ultimately forming z_1*z_1*z_3, then product->data[ p - 1 ]
          would contain z_1*z_1, and (pd_z + j) would point to the z_3
         coordinate of the current neighbor */

         product->data[ p ] = product->data[ p - 1 ] * *( pd_z + j );

         *( poly[ p ] ) = product->data[ p ] / scale->data[ p ];

       }
       else{

         scale->data[ 0 ]   = *( pd_zmax + j );
         product->data[ 0 ] = *( pd_z + j );

         *( poly[ 0 ] ) = product->data[ 0 ];
       }

       /* update pointer array element (corresponding
       to the current order) to point to the next
       term in in the same order */

       ( poly[ p ] )++;

       err = localfn_fill_polyorder(
         current_order + 1,
         j,
         max_order,
         embedding_dimension,
         poly,
         product,
         scale,
         pd_z,
         pd_zmax );
       if ( err ) return err;
     }
   }

   return MUTIL_ERR_OK;
 }

/** QR decomposition of a (non)square matrix.
 * Calculates the QR decomposition of the matrix A
 * such that A = Q*R, where Q is an orthonormal matrix
 * such that Q\^T * Q = I, and R is an upper triangular
 * matrix. All input and output matrices are expected
 * to be allocated by the user. The matrix A need not
 * be square.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = localfn_qr_decomposition( &A, intrp_ptr, &Q, &R )#
 * @param  A Pointer to a double matrix of dimension N x M..
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  Q Pointer to the resulting Q matrix of dimension N x M.
 * @param  R Pointer to the resulting R matrix of dimension M x M.
 * @return Standard mutils error code.
 * @see localfn_qr_solve
 * @private
 */
static mutil_errcode localfn_qr_decomposition(
  const double_mat *A,
  void             *intrp_ptr,
  double_mat       *Q,
  double_mat       *R )
{
  double         norm;
  double         sqrtnorm;
  double        *pd_Q;
  double        *pd_q;
  double        *pd_v;
  double        *pd_R;
  double_mat     Atemp;
  double_mat     q;
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         j;
  sint32         k;
  sint32         ncol;
  sint32         nrow;

  MUTIL_TRACE( "Starting localfn_qr_decomposition()" );

  /* initialize memory management list */

  MEMLIST_INIT( list );

  /* initilize variables */

  nrow = A->nrow;
  ncol = A->ncol;

  /* quick dimension check */

  if ( !MATANY_EQUAL_DIM( A, Q ) ){
    MUTIL_ERROR("A and Q matrix must be the same dimension for QR decomposition" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if ( R->nrow != ncol || R->ncol != ncol ){
    MUTIL_ERROR("R matrix must be M x M where M is the number of columns in A matrix" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* allocate memory */

  err = matdbl_malloc_register( &Atemp, ncol, nrow, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &q, nrow, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create a copy of the input matrix and transpose to
    facilitate C's row major format */

  localfn_transpose_fast( A, &Atemp );

  /* zero out the Q and R matrix */

  err = matdbl_assign_scalar( 0.0, intrp_ptr, Q );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_assign_scalar( 0.0, intrp_ptr, R );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* assign pointers */

  pd_R = R->data;

  /* start main loop */

  for ( i = 0; i < ncol; i++ ){

    /* extract ith column of A and store in vector v */

    pd_v = Atemp.data + i * nrow;

    /* form dot product of v with itself */

    err = localfn_vector_dot_product( pd_v, pd_v, nrow, &norm );
    MEMLIST_FREE_ON_ERROR( err, &list );

    sqrtnorm = sqrt( norm );

    /* assign sqrt( norm ) to ith diagonal element of R matrix */

    *( pd_R + i ) = sqrtnorm;

    /* form q vector and copy into ith column of Q matrix */

    pd_Q = Q->data + i;
    pd_q = q.data;

    for ( j = 0; j < nrow; j++ ){

      *pd_q = *pd_v / sqrtnorm;

      *pd_Q = *pd_q;

      /* increment pointers */

      pd_Q += ncol;
      pd_q++;
      pd_v++;
    }

    for ( j = i + 1; j < ncol; j++ ){

      /* extract jth column of A and store in vector av */

      pd_v = Atemp.data + j * nrow;

      /* form dot product of q and av */

      err = localfn_vector_dot_product( q.data, pd_v, nrow, &norm );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* assign norm to R(i,j) */

      *( pd_R + j ) = norm;

      pd_Q = Q->data + i;

      /* A(k,j) -= Q(k,i) * norm */

      for( k = 0; k < nrow; k++ ){

        (*pd_v) -= ( *pd_Q * norm );

        pd_v++;
        pd_Q += ncol;
      }

    } /* j loop */

    /* increment pointers */

    pd_R += ncol;

  } /* i loop */

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_qr_decomposition()" );

  return MUTIL_ERR_OK;
}

/** Solves for the vector X in the matrix
 * equation A * X = B using back substitution
 * via the QR decomposition of matrix A.
 * All input and output matrices are expected
 * to be allocated by the user. The A matrix need not
 * be square.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = localfn_qr_solve( &Q, &R, &B, &X )#
 * @param  Q Pointer to double matrix Q of dimension N x M.
 *  This has the same dimension as the original A matrix.
 * @param  R Pointer to a double matrix R of dimension M x M.
 * @param  B Pointer to a double matrix B, a vector of length M.
 * @param  X Pointer to a double matrix, a vector of length N containing the result.
 * @return Standard mutils error code.
 * @see localfn_qr_decomposition
 * @private
 */
static mutil_errcode localfn_qr_solve(
  const double_mat *Q,
  const double_mat *R,
  const double_mat *B,
  double_mat       *X )
{
  double         temp = 0.0;
  double        *pd_Q;
  double        *pd_R;
  double        *pd_bb;
  double        *pd_B;
  double        *pd_X;
  double_mat     bb;
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         j;
  sint32         k;
  sint32         ncol;
  sint32         nrow;

  MUTIL_TRACE( "Starting localfn_qr_solve()" );

  /* initialize memory management list */

  MEMLIST_INIT( list );

  /* initilize variables */

  nrow = Q->nrow;
  ncol = Q->ncol;

  /* allocate memory */

  err = matdbl_malloc_register( &bb, nrow, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* assign pointers */

  pd_bb = bb.data;

  /* begin back substitution */

  for( i = 0; i < ncol; i++ ){

    *pd_bb = 0.0;

    pd_Q = Q->data + i;
    pd_B = B->data;

    for(j = 0; j < nrow; j++){

      *pd_bb += *pd_Q * *pd_B;

      pd_Q += ncol;
      pd_B++;
    }

    pd_bb++;
  }

  for ( k = ( ncol - 1 ); k >= 0; k-- ){

    i = k + 1;

    pd_R = R->data + k * ncol;

    pd_X = X->data + i;

    while( i <= ( ncol - 1 ) ){

      temp += *(pd_R + i) * *pd_X;

      pd_X++;
      i++;
    }

    X->data[k] = ( bb.data[k] - temp ) / *(pd_R + k);

    temp = 0.0;
  }

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_qr_solve()" );

  return MUTIL_ERR_OK;
}

/** Computes the dot product of two vectors.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = localfn_vector_dot_product( &x, &y, 3, &dotprod )#
 * @param  x Pointer to double value (beginning of x block).
 * @param  y Pointer to double value (beginning of y block).
 * @param  n Number of points in each vector. No error checking is done
 *  on array limits so be careful.
 * @param  dotprod Pointer to a double value containing the result.
 * @return Standard mutils error code.
 * @see localfn_qr_decomposition
 * @private
 */
static mutil_errcode localfn_vector_dot_product(
  double *x,
  double *y,
  sint32  n,
  double *dotprod )
{
  sint32 i;
  double *pd_x;
  double *pd_y;

  /* initialize variables */

  *dotprod = 0.0;

  /* initialize pointers */

  pd_x = x;
  pd_y = y;

  for ( i = 0; i < n; i++ ){

    (*dotprod) += *pd_x * *pd_y;

    pd_x++;
    pd_y++;
  }

  return MUTIL_ERR_OK;
}


/** Forms a mulitplicative QR decomposition of a
 * series of matrix products.
 * Given a matrix M formed by the product
 *
 *   M = ML(L-1) * ML(L-2) * ... * ML(0)
 *
 * where ML(j) is the jth square matrix of the product,
 * this function stably computes the QR decomposition
 * of M via its formative products to form
 *
 *   M = Q * R(L-1) * R(L-2) * ... * R(0)
 *
 * The Q matrix and the R matrices are outputs of the function.
 *
 * Note: The order of the matrices in ML is expected to be
 * ML : ML[0] ML[1] ... ML[L-1] using the above notation, i.e.,
 * do not reverse the order.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = localfn_matrix_product_qr_decomposition( &ML, intrp_ptr, &Q, &RL )#
 * @param  ML Pointer to matrix set containing universal
 *  matrices of type MUTIL\_DOUBLE forming the ML matrices as defined above.
 *  Each matrix must be square and the same size.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  Q A pointer to a previously allocated universal matrix of type MUTIL\_DOUBLE
 *  to contain the Q matrix of the QR decomposition. It dimensions
 *  are the same as the ML(j) matrices.
 * @param  RL Pointer to a previously allocated matrix set containing universal
 *  matrices of type MUTIL\_DOUBLE forming the RL matrices as defined above.
 *  Each matrix must be square and the same size.
 * @return Standard mutils error code.
 * @see localfn_qr_decomposition
 * @private
 */
static mutil_errcode localfn_matrix_product_qr_decomposition(
  const mat_set *ML,
  void          *intrp_ptr,
  univ_mat      *Q,
  mat_set       *RL )
{
  memlist        list;
  mutil_errcode  err;
  sint32         dimension;
  sint32         n;
  univ_mat       T;

  MUTIL_TRACE( "Starting localfn_matrix_product_qr_decomposition()" );

  /* initialize memory management list */

  MEMLIST_INIT( list );

  /* check inputs */

  if ( ML->nelem != RL->nelem ){
    MUTIL_ERROR( "Number of matrices in matrix set in input and output matrix sets must be equal" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  dimension = MATUNIV_NCOL( &( ML->mats[0] ) );

  /* allocate memory */

  err = matuniv_malloc_register( &T, dimension, dimension, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* form QR decomposition of the first matrix in ML */

  err = localfn_qr_decomposition( &( ML->mats[0].mat.dblmat ), intrp_ptr,
    &(Q->mat.dblmat), &( RL->mats[0].mat.dblmat ) );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* now coninute QR decomposition on subsequent matrices in ML */

  for ( n = 1; n < ML->nelem; n++ ){

    /* form T = ML[n] * Q */

    err = matuniv_multiply( &( ML->mats[n] ), Q, intrp_ptr, &T );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* calculate QR decomposition of T matrix, forming new Q matrix
    and RL[n] matrix */

    err = localfn_qr_decomposition( &(T.mat.dblmat), intrp_ptr, &(Q->mat.dblmat), &( RL->mats[n].mat.dblmat ) );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_matrix_product_qr_decomposition()" );

  return MUTIL_ERR_OK;
}

/** Estimates the Jacobian matrices for points in a
 * multidimensional embedding.
 * For a selected subset of points in an multivariate embedding,
 * the Jacobian matrix is estimated using a multidimensional
 * polynomial fit to the flow of the data in local neighborhoods.
 * These Jacobian estimates are to be used in a subsequent estimation
 * of the local Lyapunov exponents (estimated outside of this function).
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = localfn_calculate_jacobian( &embedding, global_reference_point, n_reference_local, local_dimension, orbital_lag, scale_max, polyorder, forward, metric, intrp_ptr, &jacobian )#
 * @param  embedding Pointer to universal matrix of type MUTIL\_DOUBLE
 *  containing the embedding amtrix (each column is a coordinate, each row
 *  a point in the phase space.
 * @param global_reference_point Initial reference point to begin forming Jacobians (zero
 *  based indexing..
 * @param n_reference_local Number of reference points to explore beyond a given global
 *  reference point.
 * @param local_dimension While the neighbors are found in the (full) dimension
 *  of the embedding, the Jacobians are calculated in a reduced embedding space.
 *  This local\_dimension must be less than or equal to the original embedding
 *  dimension implicitly specified by the number of rows of the embedding matrix.
 * @param orbital_lag The orbital lag (number of neighbors to exlcude
 *  in time along a trajectory when searchign for nearest neighbors).
 * @param scale_max The maximum scale over which the local Lyapunv exponents will
 *  be estimated. Scale implies the number of successive reference points along
 *  an orbit that will ultimately be used to approximate the Lyapunov exponents.
 * @param polyorder The order of the polnomial to use in fitting the local neighborhoods
 *  surround a given reference point.
 * @param forward A logical flag. If TRUE, the forward local Lyapunov exponents
 *  are estimated.
 * @param metric Nearest neighbor search metric.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param jacobian A pointer to a matrix set that, upon return, will contain
 *  n\_reference + scale\_max Jacobian matrix estimates. Each matrix
 *  in the matrix set will be a d x d universal matrix of type MUTIL\_DOUBLE
 *  where d is the local\_dimension. The memory for this output is automatically
 *  allocated within the function.
 * @return Standard mutils error code.
 * @see localfn_qr_decomposition
 * @private
 */
static mutil_errcode localfn_calculate_jacobian(
  const univ_mat *embedding,
  const sint32    global_reference_point,
  const sint32    n_reference_local,
  const sint32    local_dimension,
  const sint32    orbital_lag,
  const sint32    scale_max,
  const sint32    polyorder,
  const boolean   forward,
  const fra_distance_metric metric,
  void           *intrp_ptr,
  mat_set        *jacobian )
{
  boolean        sort_distances = FALSE;
  double        *pd_first_point;
  double        *pd_last_point;
  double        *pd_neighbor;
  double        *pd_neighbor_image;
  double        *pd_reference_image;
  double        *pd_result;
  double        *pd_zmax;
  double        *pd_znew;
  double        *pd_zold;
  double_mat     Q;
  double_mat     R;
  double_mat     design_matrix;
  double_mat     polynomial_coeffs;
  double_mat     reference;
  double_mat     solution;
  double_mat     zmax;
  double_mat     znew;
  double_mat     zold;
  memlist        list;
  mutil_errcode  err;
  mutil_kdtree   kdtree;
  sint32         dims;
  sint32         embedding_dimension;
  sint32         i;
  sint32         image_index;
  sint32         image_lag;
  sint32         index;
  sint32         j;
  sint32         k;
  sint32         n;
  sint32         n_embed;
  sint32         n_polyterm;
  sint32         n_neighbor;
  sint32        *ps_index_neighbor;
  univ_mat       embedding_usable;
  univ_mat       neighbor_distances;
  univ_mat       neighbor_indices;
  univ_mat       original_indices;
  univ_mat       time_stamp;

  MUTIL_TRACE( "Starting localfn_calculate_jacobian()" );

  /* initialize memory management list */

  MEMLIST_INIT( list );

  /* initilize variables */

  embedding_dimension = MATUNIV_NCOL( embedding );
  n_embed = MATUNIV_NROW( embedding );

  n_polyterm = localfn_number_polynomial_parameters(
    local_dimension, polyorder );

  n_neighbor = n_polyterm * 2;

  /* allocate memory */

  dims = n_reference_local + scale_max;

  err = matset_malloc_register( jacobian, 1, &dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices( jacobian, local_dimension,
    local_dimension, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &zold, n_neighbor, local_dimension, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &znew, local_dimension, n_neighbor, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &zmax, 1, local_dimension, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &Q, n_neighbor, n_polyterm, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &R, n_polyterm, n_polyterm, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &polynomial_coeffs, n_polyterm, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &time_stamp, n_reference_local + scale_max, 1,
     MUTIL_SINT32, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* fill out the solution vector header
     (dimension members only, the data pointer
     will be set later on) */

  solution.ncol  = 1;
  solution.nrow  = n_neighbor;
  solution.nelem = n_neighbor;

  /* the following supresses lint warnings */

  (void) solution.nelem;
  (void) solution.ncol;
  (void) solution.nrow;

  /* For a given reference point in the embedding, we will need to access
     its image (point either before or after it in time) as well as the images
     for the neighbors of that reference point. Inasumch, we will lop off the first
     and last point of the embedding so that we will never exceed the boundaries
     of the embedding matrix when accessing the images. For example, if we find a
     neighbor associated with index 0 of the embedding (i.e., the first point)
     and we are looking to find its image one time step back, we will essentially
     be trying to access an array at index -1, which MUTILS does not support. We
     circumvent the problem by storing the end points of the embedding and thus
     remove the complication of having to resize various matrices during the processing
     of the data vectors comprising the neighborhood of a reference point */

  embedding_usable.type = MUTIL_DOUBLE;
  embedding_usable.mat.dblmat.nrow  = n_embed - 2;
  embedding_usable.mat.dblmat.ncol  = embedding_dimension;
  embedding_usable.mat.dblmat.nelem = ( n_embed - 2 ) * embedding_dimension;
  embedding_usable.mat.dblmat.data  = embedding->mat.dblmat.data + embedding_dimension;

  pd_first_point = embedding->mat.dblmat.data;
  pd_last_point  = pd_first_point + ( n_embed - 1 ) * embedding_dimension;

  /* initialize variables */

  image_lag = ( forward ) ? 1 : - 1;

  /* create a kd-tree for neighbor searching using all usable points */

  err = mutil_kdtree_malloc_register( &kdtree, &embedding_usable, (sint32) 15, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* find the nearest n_neighbor neighbors only for the first USABLE
     n_reference_local + scale_max points of the embedding starting from
     the given global_reference_point. we point
     data to the (global_reference_point + 1)th row of the ORIGINAL embedding matrix:
     a +1 offset is needed to skip the first embedding point as
     described above. global_reference_point is assumed to be given
     in zero based indexing.

     a timestamp vector of indices is also formed to map maintain
     the original index mapping of the extracted subset we are sending
     into the frauniv_neighbor_find_arbitrary() function. this is only
     necessary because we are specifying an orbital lag and need to know
     how close in time our reference points are to the neighbor candidates.
     if we did not specify a timestamp for the subset, indices 0, 1, 2....
     would be assigned to the rows of the subset by default.

     keep in mind also that index 0 of the subset actually refers to the second point
     (index 1) of the original embedding since we have skipped the first point. we
     want to keep this setup and use it to our advantage in finding images
     without overstepping array boundaries.
   */

  embedding_usable.mat.dblmat.nrow  = n_reference_local + scale_max;
  embedding_usable.mat.dblmat.nelem = ( n_reference_local + scale_max ) * embedding_dimension;
  embedding_usable.mat.dblmat.data  = embedding->mat.dblmat.data +
    ( global_reference_point + 1 ) * embedding_dimension;

  for ( i = 0; i < MATUNIV_NELEM( &time_stamp ); i++ ){

    time_stamp.mat.s32mat.data[i] = global_reference_point + i;
  }

  err = frauniv_neighbor_find_arbitrary(
    &embedding_usable,
    &kdtree,
    n_neighbor,
    (double) 0.0,
    metric,
    (univ_mat*) NULL,
    sort_distances,
    orbital_lag,
    &time_stamp,
    intrp_ptr,
    &original_indices,
    &neighbor_indices,
    &neighbor_distances );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free superfluous matrices */

  MUTIL_FREE_WARN( matuniv, &original_indices );
  MUTIL_FREE_WARN( matuniv, &neighbor_distances );

  err = memlist_member_free( &kdtree, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* register the nearest neighbor indices with the memory manager */

  err = memlist_member_register( &list, &neighbor_indices, MEMTYPE_MATUNIV );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* form reference vector, i.e. the point that
     we are forming the Jacobian around */

  reference.nelem = local_dimension;
  reference.ncol  = local_dimension;
  reference.nrow  = 1;
  reference.data  = embedding_usable.mat.dblmat.data;

  /* the following supresses lint warnings */

  (void) reference.nelem;
  (void) reference.ncol;
  (void) reference.nrow;
  (void) reference.data;

  /* set pointer to neighbor indices */

  ps_index_neighbor = neighbor_indices.mat.s32mat.data;

  /* obtain the Jacobians for the neighborhoods
     surrounding each reference point in the phase space */

  for ( n = 0, i = global_reference_point; i < global_reference_point + n_reference_local + scale_max; i++, n++ ){

    pd_zold = zold.data;

    pd_reference_image = reference.data + image_lag * embedding_dimension;

    pd_result = jacobian->mats[n].mat.dblmat.data;

    /* zero out the zmax vector */

    err = matdbl_assign_scalar( 0.0, intrp_ptr, &zmax );
    MEMLIST_FREE_ON_ERROR( err, &list );

    for ( j = 0; j < n_neighbor; j++ ){

      index       = *ps_index_neighbor;
      image_index = index + image_lag;
      ps_index_neighbor++;

      pd_znew     = znew.data + j;
      pd_zmax     = zmax.data;
      pd_neighbor = embedding->mat.dblmat.data + ( index + 1 ) * embedding_dimension;

      if ( image_index < 0 ){

        pd_neighbor_image = pd_first_point;
      }
      else if ( image_index > n_embed - 3 ){ /* number of embedding points = n_embed - 2 */

        pd_neighbor_image = pd_last_point;
      }
      else{

        pd_neighbor_image = pd_neighbor + image_lag * embedding_dimension;
      }

      /* form the z vectors:

      zold: the vector difference between the neighbor and
      the reference. each ROW of the matrix is zold
      for a different neighbor.

      znew: the vector difference between the image of the neighbor and
      the image of the reference. each COLUMN of the matrix is znew
      for a different neighbor.

      zmax: a vector to contain the max absolute value of each coordinate
      across all neighbors.

      */

      for ( k = 0; k < local_dimension; k++ ){

        *pd_zold = pd_neighbor[k] - reference.data[k];
        *pd_znew = pd_neighbor_image[k] - pd_reference_image[k];

        *pd_zmax = MUTIL_MAX( *pd_zmax, MUTIL_ABS( *pd_zold ) );

        pd_zold++;
        pd_zmax++;
        pd_znew += n_neighbor;
      }

    } /* end loop over neighbors */

    /* calculate the design matrix for a multidimensional least squares fit */

    if ( memlist_member_exist( &design_matrix, &list ) ){
      err = memlist_member_free( &design_matrix, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    err = localfn_least_squares_polynomial_design_matrix(
      &zold, &zmax, polyorder, &design_matrix );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_register( &list, &design_matrix, MEMTYPE_MATDBL );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* calculate QR decomposition */

    err = localfn_qr_decomposition( &design_matrix, intrp_ptr, &Q, &R );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* loop through each dimension of the
       znew matrix forming a 'solution' vector to
       the matrix equation: znew[j] = design_matrix * coeffs
       where znew[j] is the jth column of znew with dimensions
       n_neighbor x 1, the design matrix has dimensions
       n_neighbor x n_polyterm, and coeffs is an
       n_polyterm x 1 vector to contain the polynomial coefficients
       which best fit the model */

    for ( j = 0; j < local_dimension; j++ ){

    /* point solution vector to the proper row (dimension) of the
      znew matrix. */

      solution.data = znew.data + j * n_neighbor;

      /* now perform backsubstitution to solve for
      polynomial coefficients in the model */

      err = localfn_qr_solve( &Q, &R, &solution, &polynomial_coeffs );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* store subset of coefficients in the proper location in
      the Jacobian matrix set */

      for ( k = 0; k < local_dimension; k++ ){

        *pd_result = polynomial_coeffs.data[ k ];

        pd_result++;
      }
    }

    /* increment variables and pointers */

    reference.data += embedding_dimension;

  } /* end loop over each reference point */


  /* unregister memory associated with output */

  err = memlist_member_unregister( jacobian, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_calculate_jacobian()" );

  return MUTIL_ERR_OK;
}


/** Tranpose of matrix set matrices.
 * Forms the transpose of each matrix in a matrix
 * set.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = localfn_matset_transpose( &x, intrp_ptr, &transpose_x )#
 * @param  x Pointer to matrix set of universal matrices of type MUTIL\_DOUBLE.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  transpose_x Pointer to matrix set that, upon return, will contain the
 *    result. The memory for the matrix set headers and universal matrices
 *    are automatically allocated within the function.
 * @return Standard mutils error code.
 * @see frauniv_local_lyapunov_spectrum
 * @private
 */
static mutil_errcode localfn_matset_transpose(
  const mat_set *x,
  mat_set       *transpose_x )
{
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         ncol;
  sint32         nmats;
  sint32         nrow;

  MUTIL_TRACE( "Starting localfn_matset_transpose()" );

  /* initialize memory management list */

  MEMLIST_INIT( list );

  /* gather dimensions */

  nmats = x->nelem;
  nrow  = MATUNIV_NROW( &( x->mats[0] ) );
  ncol  = MATUNIV_NCOL( &( x->mats[0] ) );

  err = matset_malloc_register( transpose_x, 1, &nmats, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices( transpose_x, ncol, nrow, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* loop through each matrix and form transpose */

  for ( i = 0; i < nmats; i++ ){

    localfn_transpose_fast( &( x->mats[i].mat.dblmat ),
      &( transpose_x->mats[i].mat.dblmat ) );
  }

  /* unregister memory associated with output */

  err = memlist_member_unregister( transpose_x, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_matset_transpose()" );

  return MUTIL_ERR_OK;
}


/** Tranpose of matrix set matrices.
 * Fast transpose of a double matrix.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_lyap.c
 * @library fractal
 * @usage #err = localfn_transpose_fast( &x, &transpose_x )#
 * @param  x Pointer to double matrix.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param  transpose_x Pointer to double matrix that, upon return, will contain the
 *    result.
 * @return Standard mutils error code.
 * @see frauniv_local_lyapunov_spectrum
 * @see localfn_matset_transpose
 * @private
 */
static void localfn_transpose_fast( const double_mat *x, double_mat *transpose_x )
{
  sint32 i;
  sint32 j;

  double *pd_x = x->data;
  double *pd_transpose_x;

  for ( i = 0; i < x->nrow; i++ ){

    pd_transpose_x = transpose_x->data + i;

    for ( j = 0; j < x->ncol; j++ ){

      *pd_transpose_x = *pd_x;
      pd_x++;
      pd_transpose_x += x->nrow;
    }

  }

}
