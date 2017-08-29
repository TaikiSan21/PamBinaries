
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_tdmi.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */


/* This file contains definitions for functions used to compute */
/* time delayed mutual information for a given time series.     */

#include "fra_tdmi.h"

#include "fra_kde.h"
#include "mat_univ.h"
#include "mat_comp.h"
#include "mat_sort.h"
#include "mat_summ.h"
#include "mat_umat.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_math.h"
#include "ut_mem.h"


/* Macros */

#define FRA_TDMI_LOG_ARG_MIN   ((double) 1.175494351e-38) /* FLT_MIN */

#undef LOCALDEF_CHECK_NULL_POINTER_STAT
#define LOCALDEF_CHECK_NULL_POINTER_STAT( DATA_PTR,         \
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
static double localfn_get_arb_data(
  const univ_mat   *data,
  const sint32      index );


/*
 *********************************
 *                               *
 * LIBRARY FUNCTION DEFINITIONS  *
 *                               *
 *********************************
 */

/* Time delayed mutual information. */
/*                                  */
/* Documented in fra_tdmi.h         */
/* Written by Keith L. Davidson     */
/* Subsequently altered by William Constantine */

mutil_errcode frauniv_time_delayed_mutual_information(
  const univ_mat   *time_series,
  const univ_mat   *lags,
  void             *intrp_ptr,
  univ_mat         *tdmi )
{
  double          tmp_dbl;
  double         *pd_result;
  memlist         list;
  mutil_errcode   err;
  sint32          d;
  sint32          i;
  sint32          j;
  sint32          max_lag;
  sint32          min_lag;
  sint32          n_lag;
  sint32          n_sample;
  sint32          neg;
  sint32          tmp_n_sample;
  sint32          two_x_lag;
  sint32         *ps_lag;
  univ_mat        data_2d;
  univ_mat        pdf_1d;
  univ_mat        pdf_2d;
  univ_mat        unique_lags;
  univ_mat        time_series_column;

  MUTIL_TRACE( "Start in frauniv_time_delay_mutual_info()" );

  /* Avoid lint warning */

  (void) whatssi;

  /* Initialize the memory management list */

  MEMLIST_INIT( list );

  /* Error checks */

  /* time_series */

  LOCALDEF_CHECK_NULL_POINTER_STAT( time_series, univ_mat, matuniv );

  if ( time_series->type == MUTIL_DCOMPLEX ) {
    MUTIL_ERROR( "Input time_series must be a real data type" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  n_sample = MATUNIV_NELEM( time_series );

  if ( !MATANY_IS_VEC( &( time_series->mat.dblmat ) ) ){
    MUTIL_ERROR( "Input time series must be a single-row or single-column matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( n_sample < (sint32) 3 ) {
    MUTIL_ERROR( "Input time_series must have at least 3 elements" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* lags vector */

  LOCALDEF_CHECK_NULL_POINTER_STAT( lags, univ_mat, matuniv );

  if ( lags->type != MUTIL_SINT32 ) {
    MUTIL_ERROR( "Input lags must be of type MUTIL_SINT32" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  n_lag = MATUNIV_NELEM( lags );

  if ( !MATANY_IS_VEC( &( lags->mat.s32mat ) ) ){
    MUTIL_ERROR( "Input lag series must be a single-row or single-column matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... make sure all lags are non-negative */

  err = mats32_number_less_than_scalar( &( lags->mat.s32mat ), (sint32) 0, intrp_ptr, &neg );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( neg > 0 ){
      MUTIL_ERROR( "Input lags must have non-negative elements" );
      return MUTIL_ERR_ILLEGAL_VALUE;
  }

  err = mats32_range( &( lags->mat.s32mat ), intrp_ptr, &min_lag, &max_lag );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* Make sure delay isn't too large */

  if ( max_lag > ( n_sample - 2 ) ) {
    MUTIL_ERROR( "Not enough data points for desired delay(s)" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* make sure all lag values are unique */

  err = matuniv_unique( lags, FALSE, intrp_ptr, &unique_lags );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_validate( &unique_lags );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( !err ){

    if ( MATUNIV_NELEM( &unique_lags ) != n_lag ){
      MUTIL_ERROR( "Redundant lag values found in lag vector" );
      MUTIL_FREE_WARN( memlist, &list );
      MUTIL_FREE_WARN( matuniv, &unique_lags );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
  }

  MUTIL_FREE_WARN( matuniv, &unique_lags );

  /* Allocate memory for 2-D delay series and copy original time series */

  /* into the first column.                                             */

  err = matuniv_malloc_register(
    &data_2d,
    n_sample,
    (sint32) 2,
    MUTIL_DOUBLE,
    &list );
  if ( err ) return MUTIL_ERR_MEM_ALLOC;
  for ( i = 0; i < n_sample; i++ ) {
    data_2d.mat.dblmat.data[ 2 * i ] = localfn_get_arb_data( time_series, i );
  }

  /* Allocate memory for the mutual information */

  err = matuniv_malloc_register(
    tdmi,
    lags->mat.s32mat.nrow,
    lags->mat.s32mat.ncol,
    MUTIL_DOUBLE,
    &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* Compute the univariate PDF estimate of the original time series.
    the KDE function expects the time series to be in a single-column,
	otherwise in interprets all N points in the time series as a
    single point in N dimensions */

  time_series_column.type             = MUTIL_DOUBLE;
  time_series_column.mat.dblmat.ncol  = 1;
  time_series_column.mat.dblmat.nrow  = n_sample;
  time_series_column.mat.dblmat.nelem = n_sample;
  time_series_column.mat.dblmat.data  = time_series->mat.dblmat.data;

  err = frauniv_kernel_density_estimate(
    &time_series_column,
    (univ_mat*) NULL,
    intrp_ptr,
    &pdf_1d );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* Add the 1-D pdf to the memory management list */

  err = memlist_member_register( &list, (void*) &pdf_1d, MEMTYPE_MATUNIV );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* assign pointers */

  ps_lag    = lags->mat.s32mat.data;
  pd_result = tdmi->mat.dblmat.data;

  /* Compute mutual information for each delay */

  for ( d = 0; d < n_lag; d++ ) {

    /* temporary number of points in delay matrix data_2d */

    tmp_n_sample = n_sample - *ps_lag;

    /* used to skip through data_2d */

    two_x_lag = 2 * *ps_lag;

    /* copy the delayed time series into the second column of data_2d */

    for ( i = 0; i < tmp_n_sample; i++ ) {
      j = 2 * i;
      data_2d.mat.dblmat.data[ j + 1 ] = data_2d.mat.dblmat.data[ j + two_x_lag ];
    }

    /* make data_2d appear to have fewer rows (depends on the delay value) */

    data_2d.mat.dblmat.nrow  = tmp_n_sample;
    data_2d.mat.dblmat.nelem = tmp_n_sample * 2;

    /* compute the bivariate PDF estimate of the delay matrix */

    err = frauniv_kernel_density_estimate(
      &data_2d,
      (univ_mat*) NULL,
      intrp_ptr,
      &pdf_2d );

    /* correct the dimensions in case we exit on an error */

    data_2d.mat.dblmat.nrow  = n_sample;
    data_2d.mat.dblmat.nelem = n_sample * 2;

    /* avoid lint warnings */

    (void) data_2d.mat.dblmat.nrow;
    (void) data_2d.mat.dblmat.nelem;

    /* exit if 2-D pdf was not computed */

    MEMLIST_FREE_ON_ERROR( err, &list );

    /* compute mutual information for current delay */

    *pd_result = (double) 0.0;

    for ( i = 0; i < tmp_n_sample; i++ ) {

      tmp_dbl = pdf_1d.mat.dblmat.data[ i ] *
        pdf_1d.mat.dblmat.data[ i + *ps_lag ];

      tmp_dbl = pdf_2d.mat.dblmat.data[ i ] /
        MUTIL_MAX( FRA_TDMI_LOG_ARG_MIN, tmp_dbl );

      *pd_result += log( MUTIL_MAX( FRA_TDMI_LOG_ARG_MIN, tmp_dbl ) );
    }

    *pd_result /= tmp_n_sample;

    /* Free memory for the 2-D pdf estimate */

    MUTIL_FREE_WARN( matuniv, &pdf_2d );

    /* increment pointers */

    ps_lag++;
    pd_result++;

  } /* end loop over each lag */

  /* unregister the mutual information from the memory management list */

  err = memlist_member_unregister( (void*) tdmi, &list );
  MUTIL_FREE_WARN( memlist, &list );
  if ( err ) return err;

  MUTIL_TRACE( "Done with frauniv_time_delay_mutual_info()" );

  return MUTIL_ERR_OK;
}

/*
 *******************************
 *                             *
 * STATIC FUNCTION DEFINITIONS *
 *                             *
 *******************************
 */

/** Function to extract data from an arbitrary universal matrix.
 * The frauniv_time_delayed_mutual_information() function supports all real
 * data types and calls this function to copy the input data into
 * a universal matrix of type MUTIL_DOUBLE.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source fra\_tdmi.c
 * @library fractal
 * @usage #err = localfn_get_arb_data(&data,idx);#
 * @return Standard mutils error/OK code.
 * @param data   Pointer to matrix data from which to estimate the density,
 *               must be of type MUTIL\_DOUBLE. The number of columns is
 *               the dimension of the space.
 * @param points Pointer to a set of points at which the density will be
 *               estimates. Must be of type MUTIL_DOUBLE and have the same
 *               number of columns as data.
 *
 * @see frauniv_time_delayed_mutual_information
 * @private
 */

 static double localfn_get_arb_data(
   const univ_mat   *data,
   const sint32      index )
 {
   double   answer;

   MUTIL_TRACE( "Start in localfn_get_arb_data()" );

   switch ( data->type ) {

   case MUTIL_UINT8:

     answer = (double) data->mat.u8mat.data[ index ];
     break;

   case MUTIL_SINT8:

     answer = (double) data->mat.s8mat.data[ index ];
     break;

   case MUTIL_UINT16:

     answer = (double) data->mat.u16mat.data[ index ];
     break;

   case MUTIL_SINT16:

     answer = (double) data->mat.s16mat.data[ index ];
     break;

   case MUTIL_UINT32:

     answer = (double) data->mat.u32mat.data[ index ];
     break;

   case MUTIL_SINT32:

     answer = (double) data->mat.s32mat.data[ index ];
     break;

   case MUTIL_FLOAT:

     answer = (double) data->mat.fltmat.data[ index ];
     break;

   case MUTIL_DOUBLE:

     answer = (double) data->mat.dblmat.data[ index ];
     break;

   default:
     /* will never get here because of error checks */
     answer = 0.0; /* BC: avoid compiler warning */

     break;
   }

   MUTIL_TRACE( "Done with localfn_get_arb_data()" );

   return answer;
 }
