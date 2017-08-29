
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_scale.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */

#include "fra_scale.h"
#include "mat_io.h"

#include "mat_assn.h"
#include "mat_stat.h"
#include "mat_type.h"
#include "mat_umat.h"
#include "mat_univ.h"

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
  finding appropriate linear scaling regions
  in a time series.

  The functions are declared in fra_scale.h
*/

/* Static macro definitions */

#undef LOCALDEF_CHECK_NULL_POINTER_SCALE
#define LOCALDEF_CHECK_NULL_POINTER_SCALE( DATA_PTR,        \
  DATA_TYPE, TYPE_PREFIX )                                  \
  err = TYPE_PREFIX ## _validate( DATA_PTR );               \
  if ( err ) return err;                                    \
  if ( DATA_PTR == ( DATA_TYPE * ) NULL ) {                 \
  MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
  return MUTIL_ERR_NULL_POINTER;                            \
  }

/* Static functions declared here and defined at end of file */

static mutil_errcode localfn_least_squares_slope(
  univ_mat *x,
  univ_mat *y,
  double   *result );

/* Finding linear trends in an arbitrary time series */
/*                                                   */
/* Documented in fra_scale.h                         */
/* Written by William Constantine                    */

mutil_errcode frauniv_piecwise_linear_segmentation(
  const univ_mat *xdata,
  const univ_mat *ydata,
  const sint32    n_fit,
  const double    angle_tolerance,
  void           *intrp_ptr,
  univ_mat       *result )
{
  double         dangle;
  double         new_slope;
  double         slope;
  double         sum_slope;
  memlist        list;
  mutil_errcode  err;
  sint32         count;
  sint32         n_avg;
  sint32         n_sample;
  sint32         right;
  univ_mat       x;
  univ_mat       y;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_piecwise_linear_segmentation()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /*** check input data ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_SCALE( xdata, univ_mat, matuniv );
  LOCALDEF_CHECK_NULL_POINTER_SCALE( ydata, univ_mat, matuniv );

  /* ... if number of elements is positive */

  if ( MATUNIV_NELEM( xdata ) < 1 ){
    MUTIL_ERROR( "Number of elements in input xdata matrix must be positive." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( MATUNIV_NELEM( ydata ) < 1 ){
    MUTIL_ERROR( "Number of elements in input ydata matrix must be positive." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if the type is double */

  if ( xdata->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Input xdata matrix must be of type double." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( ydata->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Input ydata matrix must be of type double." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( &( xdata->mat.dblmat ) ) ){
    MUTIL_ERROR( "Input xdata must be a single-row or single-column matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( !MATANY_IS_VEC( &( ydata->mat.dblmat ) ) ){
    MUTIL_ERROR( "Input ydata must be a single-row or single-column matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... that the lengths of xdata and ydata are the same */

  if ( MATUNIV_NELEM( xdata ) != MATUNIV_NELEM( ydata ) ){
    MUTIL_ERROR( "Input xdata and ydata must have the same number of elements." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /*** check number of points to fit ... ***/

  if ( n_fit <= 0 ){
    MUTIL_ERROR( "Number of points to fit must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  if ( n_fit > MATUNIV_NELEM( xdata ) ){
    MUTIL_ERROR( "Number of points to fit exceeds length of input xdata" );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  /*** check angle_tolerance */

  if ( angle_tolerance <= 0.0 || angle_tolerance >= 180.0 ){
    MUTIL_ERROR( "Angle tolerance must be 0 < angle_tolerance < 180." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  /* wrap temporary universal matrices */

  err = matuniv_wrap_univ_matrix( &x, (univ_mat *) xdata );
  if ( err ) return( err );

  err = matuniv_wrap_univ_matrix( &y, (univ_mat *) ydata );
  if ( err ) return( err );

  /* reset the dimensions of the temporary vectors */

  x.mat.dblmat.nelem = n_fit;
  x.mat.dblmat.nrow  = n_fit;
  x.mat.dblmat.ncol  = 1;

  y.mat.dblmat.nelem = n_fit;
  y.mat.dblmat.nrow  = n_fit;
  y.mat.dblmat.ncol  = 1;

  /* allocate memory for the result */

  err = matuniv_malloc_register( result, 1, 1, MUTIL_SINT32, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /*  initialize variables */

  n_sample = MATUNIV_NELEM( xdata );
  count    = 0;

  err = localfn_least_squares_slope( &x, &y, &slope );
  MEMLIST_FREE_ON_ERROR( err, &list );

  sum_slope = slope;
  n_avg     = 1;
  right     = n_fit - 1;

  while ( right < n_sample ){

    /*  increment data pointers */

    x.mat.dblmat.data++;
    y.mat.dblmat.data++;
    right++;

    /*  find new slope */

    err = localfn_least_squares_slope( &x, &y, &new_slope );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /*
      compare current slope to running average.
      if the deviation is greater than that
      specified by the user mark a break point.
      otherwise, add the current slope to the
      running average
    */

    dangle = MUTIL_ABS( atan( new_slope ) - atan( slope ) ) * 180.0 / MUTIL_PI;

    if ( dangle > angle_tolerance ){

      if ( n_sample - right >= n_fit ){

	err = matuniv_realloc_register( result, ++count, 1, &list );
	MEMLIST_FREE_ON_ERROR( err, &list );

	result->mat.s32mat.data[ count - 1 ] = right - 1;
      }

      /*  reinitialize slope variables */

      x.mat.dblmat.data = &( xdata->mat.dblmat.data[ right ] );
      y.mat.dblmat.data = &( ydata->mat.dblmat.data[ right ] );
      right += n_fit - 1;

      if ( right < n_sample ){

	err = localfn_least_squares_slope( &x, &y, &slope );
	MEMLIST_FREE_ON_ERROR( err, &list );

	sum_slope = slope;
	n_avg     = 1;
      }
    }
    else{

      /* update running average of slope */

      n_avg++;
      sum_slope += new_slope;
      slope      = sum_slope / (double) n_avg;
    }

    if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }

  } /* end loop over sliding window */

  /* if no breaks were found, return only the last index */

  if ( count == 0 ){
    result->mat.s32mat.data[ 0 ] = n_sample - 1;
  }

  /* free nodes corresponding to registered
     memory for the result, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
     corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with frauniv_piecwise_linear_segmentation()" );

  return MUTIL_ERR_OK;
}

static mutil_errcode localfn_least_squares_slope(
  univ_mat *x,
  univ_mat *y,
  double   *result )
{
  sint32 i;
  sint32 nelem  = MATUNIV_NELEM( x );
  double meanx  = 0.0;
  double meany  = 0.0;
  double sum_xy = 0.0;
  double sum_xx = 0.0;
  double xshift;
  double *pdx;
  double *pdy;

  MUTIL_TRACE( "Start with localfn_least_squares_slope()" );

  /* create pointers */

  pdx = x->mat.dblmat.data;
  pdy = y->mat.dblmat.data;

  /* find mean of xdata and ydata */

  for ( i = 0; i < nelem; i++ ){

    meanx += pdx[ i ];
    meany += pdy[ i ];
  }

  meanx /= (double) nelem;
  meany /= (double) nelem;

  /* find least squares slope */

  for ( i = 0; i < nelem; i++ ){

    xshift = pdx[ i ] - meanx;

    sum_xy += xshift * ( pdy[ i ] - meany );
    sum_xx += xshift * xshift;
  }

  *result = sum_xy / sum_xx;

  MUTIL_TRACE( "Done with localfn_least_squares_slope()" );

  return MUTIL_ERR_OK;
}
