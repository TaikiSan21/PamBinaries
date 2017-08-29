
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_kde.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */


/* This file contains definitions for functions used to compute */
/* multi-dimensional kernel density estimates.                  */

#include "fra_kde.h"

#include "fra_neig.h"
#include "mat_umat.h"
#include "mat_arit.h"
#include "mat_assn.h"
#include "mat_univ.h"
#include "mth_stat.h"

#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_limit.h"
#include "ut_math.h"


/* Macros */
#define FRA_MINIMUM_PROBABILITY   1e-38


/* Global variables */
static sint32            global_dim;


/*
 ****************************************
 *                                      *
 * STATIC (local) FUNCTION DECLARATIONS *
 *                                      *
 ****************************************
 */

static mutil_errcode localfn_compute_kde(
  const univ_mat         *data,
  const univ_mat         *points,
  void                   *intr_ptr,
  univ_mat               *kde );


static mutil_errcode localfn_compute_bw(
  const univ_mat   *data,
  double           *bw,
  double           *num_ops );


/*
 *********************************
 *                               *
 * LIBRARY FUNCTION DEFINITIONS  *
 *                               *
 *********************************
 */

/* Multivariate kernel density estimation. */
/*                                         */
/* Documented in fra_kde.h                 */
/* Written by Keith L. Davidson            */
mutil_errcode frauniv_kernel_density_estimate(
  const univ_mat         *data,
  const univ_mat         *points,
  void                   *intr_ptr,
  univ_mat               *kde )
{
  const univ_mat   **points_ptr;
  mutil_errcode      err;

  MUTIL_TRACE( "Start in frauniv_kernel_density_estimate()" );

  /* Avoid lint warnings */
  ( void ) whatssi;


  /* Check input for errors */

  /* data */
  if ( data == (univ_mat*) NULL ) {
    MUTIL_ERROR( "Input data is NULL pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }
  err = matuniv_validate( data );
  if ( err ) {
    MUTIL_ERROR( "Input data points to an invalid universal matrix" );
    return err;
  }
  if ( MATUNIV_NROW( data ) < (sint32) 2 ) {
    MUTIL_ERROR( "Cannot compute density estimate using a single point" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if ( data->type != MUTIL_DOUBLE ) {
    MUTIL_ERROR( "Input data must be of type MUTIL_DOUBLE" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  global_dim = MATUNIV_NCOL( data );


  /* points */
  if ( points != (univ_mat*) NULL ) {
    err = matuniv_validate( points );
    if ( err ) {
      MUTIL_ERROR( "Input points points to an invalid universal matrix" );
      return err;
    }
    if ( MATUNIV_NCOL( points ) != global_dim ) {
      MUTIL_ERROR( "Input points has invalid number of columns" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
    if ( points->type != MUTIL_DOUBLE ) {
      MUTIL_ERROR( "Input points must be of type MUTIL_DOUBLE" );
      return MUTIL_ERR_ILLEGAL_TYPE;
    }
    points_ptr = &points;

  } else {
    points_ptr = &data;
  }

  /* Compute kernel density estimate */
  err = localfn_compute_kde(
    data,
    *points_ptr,
    intr_ptr,
    kde );
  if ( err == MUTIL_ERR_ZERO_NEIGHBORS_FOUND ) {
    return MUTIL_ERR_OK;
  } else if ( err ) {
    MUTIL_ERROR( "Could not compute kernel desity estimate" );
    return err;
  }


  MUTIL_TRACE( "Finished with frauniv_kernel_density_estimate()" );

  return MUTIL_ERR_OK;
}




/*
 *******************************
 *                             *
 * STATIC FUNCTION DEFINITIONS *
 *                             *
 *******************************
 */

/** Work function for frauniv_kernel_density_estimate(). This function
 * actually computes the kernel density estimate.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_kde.c
 * @library fractal
 * @usage #err = localfn_compute_kde(&data,&pts,TRUE,intrp_ptr,&density);#
 * @return Standard mutils error/OK code.
 * @param data   Pointer to matrix data from which to estimate the density,
 *               must be of type MUTIL\_DOUBLE. The number of columns is
 *               the dimension of the space.
 * @param points Pointer to a set of points at which the density will be
 *               estimates. Must be of type MUTIL_DOUBLE and have the same
 *               number of columns as data.
 * @param kde    A pointer to a universal matrix which, upon return, will
 *               hold the computed density estimate.
 *
 * @private
 */
static mutil_errcode localfn_compute_kde(
  const univ_mat         *data,
  const univ_mat         *points,
  void                   *intrp_ptr,
  univ_mat               *kde )
{
  register double   *dbl_ptr;
  register double   *ndist_ptr;
  register double    reg_bw;
  register double    nfactor;
  register sint32   *oidx_ptr;
  register sint32    i;
  register sint32    num_neig;

  const sint32       npoints = points->mat.dblmat.nrow;

  univ_mat           orig_idx;
  univ_mat           neig_idx;
  univ_mat           neig_dist;
  double             bw;
  double             num_ops = (double) 0;
  double             tmp_nfactor;
  memlist            mlist;
  mutil_kdtree       kdtree;
  mutil_errcode      err;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start in localfn_compute_kde()" );

  /* Initialize the memory management list */
  MEMLIST_INIT( mlist );


  /* Allocate necessary memory */
  err = matuniv_malloc_register( kde, npoints, (sint32) 1, MUTIL_DOUBLE, &mlist );
  if ( err ) {
    MUTIL_ERROR( "Could not allocate memory for the density" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }
  err = matdbl_assign_scalar(
    (double) FRA_MINIMUM_PROBABILITY,
    intrp_ptr,
    &(kde->mat.dblmat) );
  if ( err ) {
    MUTIL_ERROR( "Could not initialize density matrix to zeros" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }


  /* Compute a kd-tree for the data points and register it with the
     memory manager */
  err = mutil_kdtree_malloc( &kdtree, data, (sint32) 1 );
  if ( err ) {
    MUTIL_ERROR( "Could not compute kd-tree for input data" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }
  err = memlist_member_register( &mlist, (void*) &kdtree, MEMTYPE_KDTREE );
  if ( err ) {
    MUTIL_ERROR( "Could not allocate memory for the density" );
    MUTIL_FREE_WARN( mutil_kdtree, &kdtree );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }
  num_ops += data->mat.dblmat.nelem;

  /* Compute bandwidth used to compute pilot density */
  err = localfn_compute_bw( data, &bw, &num_ops );
  if ( err ) {
    MUTIL_ERROR( "Could not compute kernel bandwidths" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }


  /* Check for user interrupt */
  if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
    MUTIL_ERROR( "User interrupt" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return MUTIL_ERR_INTERRUPT;
  }


  /* Compute density */

  /* get neighbors within bandwidth distance (of each point) */
  err = frauniv_neighbor_find_arbitrary(
    points,
    &kdtree,
    (sint32) 0,
    bw,
    FRA_DISTANCE_L2,
    (univ_mat*) NULL,
    FALSE,
    (sint32) 0,
    (univ_mat*) NULL,
    intrp_ptr,
    &orig_idx,
    &neig_idx,
    &neig_dist );
  if ( err == MUTIL_ERR_ZERO_NEIGHBORS_FOUND ) {
    /* return pdf with all zeros */
    err = memlist_member_unregister( (void*) kde, &mlist );
    MUTIL_FREE_WARN( memlist, &mlist );
    if ( err ) return err;
    return MUTIL_ERR_ZERO_NEIGHBORS_FOUND;
  } else if ( err ) {
    MUTIL_ERROR( "Nearest neighbor search failed" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }
  MUTIL_FREE_WARN( matuniv, &neig_idx );

  /* record the total number of neighbors found and initialize
  data pointers */
  num_neig  = orig_idx.mat.s32mat.nelem;
  oidx_ptr  = orig_idx.mat.s32mat.data;
  ndist_ptr = neig_dist.mat.dblmat.data;
  dbl_ptr   = kde->mat.dblmat.data;

  /* compute Epanechnikov kernel normalization factor */
  tmp_nfactor = ( global_dim + 2.0 ) * mth_gamma( global_dim / 2.0 + 1.0 )
                 / (2.0 * MUTIL_POW( MUTIL_PI, global_dim / 2.0 ) );
  nfactor     = tmp_nfactor / MUTIL_POW( bw, (double) global_dim );

  /* begin accumulating probability */
  reg_bw = bw;
  for ( i = 0; i < num_neig; i++ ) {

    dbl_ptr[ oidx_ptr[ i ] ] +=

      nfactor * ( 1.0 - ndist_ptr[ i ] * ndist_ptr[ i ] / ( reg_bw * reg_bw ) );
  }
  num_ops += 6.0 * num_neig + 10;


  /* free unneeded memory */
  MUTIL_FREE_WARN( matuniv, &orig_idx );
  MUTIL_FREE_WARN( matuniv, &neig_dist );


  /* Check for user interrupt */
  if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
    MUTIL_ERROR( "User interrupt" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return MUTIL_ERR_INTERRUPT;
  }


  /* Divide density estimate by npoints */
  err = matdbl_multiply_scalar(
    &(kde->mat.dblmat),
    1.0 / npoints,
    intrp_ptr,
    &(kde->mat.dblmat) );
  if ( err ) {
    MUTIL_FREE_WARN( memlist, &mlist );
    return err;
  }
  num_ops += kde->mat.dblmat.nelem + 1.0;

  /* Check for user interrupt */
  if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
    MUTIL_ERROR( "User interrupt" );
    MUTIL_FREE_WARN( memlist, &mlist );
    return MUTIL_ERR_INTERRUPT;
  }

  /* Remove density estimate matrix from memory manager */

  err = memlist_member_unregister( (void*) kde, &mlist );
  MUTIL_FREE_WARN( memlist, &mlist );
  if ( err ) return err;

  MUTIL_TRACE( "Finished with localfn_compute_kde()" );

  return MUTIL_ERR_OK;
}


/** Function to compute the bandwidth used for kernel density estimate.
 * The bandwidth is a smoothness parameter. The function computes the
 * minimum variance of all dimensions (columns) and uses that in Scott's
 * rule for kernel bandwidth. Increased bandwidth decreases variance at
 * the expense of bias.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_kde.c
 * @library fractal
 * @usage #err = localfn_compute_bw(&data,&bandwidth,&ops);#
 * @return Standard mutils error/OK code.
 * @param data   Pointer to matrix data from which the density is to be
 *               estimated, must be of type MUTIL\_DOUBLE. The number of
 *               columns is the dimension of the space.
 * @param bw     Pointer to a double, which will hold the computed bandwidth.
 * @param num_ops  Pointer to a double which holds the number of operations
 *               performed by this function. Used for user interrupts.
 * @private
 */
static mutil_errcode localfn_compute_bw(
  const univ_mat   *data,
  double           *bw,
  double           *num_ops )
{
  register double  *dptr;
  register double   mn;
  register double   tmp;
  register double   var;
  register sint32   dim;
  register sint32   npts;
  register sint32   d;
  register sint32   i;

  MUTIL_TRACE( "Start in localfn_compute_bw()" );

  /* Initialize data pointers and a local copies other variables */
  dptr = (double*) MATUNIV_DATA( data );
  npts = MATUNIV_NROW( data );
  dim  = global_dim;
  *bw  = MUTIL_DOUBLE_MAX;

  /* Compute variances for each dimension */
  for ( d = 0; d < global_dim; d++ ) {

    /* compute the sample mean for the current column */
    mn = (double) 0;
    for ( i = 0; i < npts; i++ )
      mn += dptr[ i * dim + d ];
    mn /= (double) npts;

    *num_ops += 4.0 * npts + 1;

    /* compute the (unbiased) variance estimate */
    var = (double) 0;
    for ( i = 0; i < npts; i++ ) {
      tmp  = dptr[ i * dim + d ] - mn;
      var += tmp * tmp;
    }

    *num_ops += 6.0 * npts + 1;

    /* hold the smallest variance over all dimensions */
    *bw = MUTIL_MIN( *bw, var );
  }

  /* Now use Scott's rule and the smallest variance to compute
     a final bandwidth */
  *bw = sqrt( *bw / (npts-1) ) * MUTIL_POW( npts, ((double) -1) / (dim+4) );

  *num_ops += 7.0;


  MUTIL_TRACE( "Finished with localfn_compute_bw()" );

  return MUTIL_ERR_OK;
}
