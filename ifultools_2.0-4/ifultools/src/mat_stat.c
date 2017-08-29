
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_stat.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include <math.h>

#include "mat_stat.h"

#include "mat_arit.h"
#include "mat_assn.h"
#include "mat_cast.h"
#include "mat_set.h"
#include "mat_sort.h"
#include "mat_summ.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "mat_usca.h"

#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_math.h"

/* This file contains implementations of the functions in the
   mat_stat.h header file, which are functions for matrix statistics.
*/


/*
***************************
 Universal matrix functions
***************************
 */


/* function documented in mat_stat.h */
/* written by Jennifer Hodgdon */
mutil_errcode matuniv_median( const univ_mat *mat, void *intrp_ptr,
  double *median )
{

  mutil_errcode trouble;
  sint32_mat    index;
  univ_mat      uindex;
  univ_mat      partials;
  sint32        inds[2];
  sint32        nelem;
  sint32        nrow;
  sint32        ncol;

  MUTIL_TRACE("Start matuniv_median()");

  /* sanity checks */

  trouble = matuniv_validate( mat );
  if( trouble ) return trouble;

  if( !median ) {
    MUTIL_ERROR( "NULL pointer for returned value" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* The median is either the middle sorted element, or the
     average of the two middle sorted elements */

  /* so figure out which 1 or 2 indices we need to sort for */

  nelem   = MATUNIV_NELEM( mat );
  nrow    = MATUNIV_NROW( mat );
  ncol    = MATUNIV_NCOL( mat );

  /* if nelem is odd inds[0] == inds[1] == (nelem-1)/2 */
  /* if nelem is even inds[1] == inds[0] + 1  */
  inds[0] = inds[1] = ( nelem - 1 ) / 2;
  if( ( nelem % 2 ) == 0 ) {
    inds[1]++;
  }

  trouble = matuniv_wrap_data( &partials, inds, 2, 1, MUTIL_SINT32 );
  if( trouble ) return trouble;

  /* call sort to find the one or two elements we need */

  trouble = mats32_malloc( &index, nrow, ncol );
  if( trouble ) return trouble;
  trouble = matuniv_wrap_matrix( &uindex, &index, MUTIL_SINT32 );
  if( trouble ) {
    MUTIL_FREE_WARN( mats32, &index );
    return trouble;
  }

  trouble = matuniv_sort_index_partial( mat, &partials, intrp_ptr, &uindex );
  if( trouble ) {
    MUTIL_FREE_WARN( mats32, &index );
    return trouble;
  }

  *median = ( MATUNIV_ELEM( mat, 0, index.data[ inds[0] ]) +
    MATUNIV_ELEM( mat, 0, index.data[ inds[1] ]) ) * 0.5;

  MUTIL_FREE_WARN( mats32, &index );

  MUTIL_TRACE("matuniv_median() done");
  return MUTIL_ERR_OK;
}


/* function documented in mat_stat.h */
/* written by Luca Cazzanti */
mutil_errcode matuniv_mean_variance(
  const univ_mat *mat,
  boolean         unbiased,
  void           *intrp_ptr,
  double         *mean,
  double         *variance )
{
  mutil_errcode    err;
  mutil_data_type  type;
  univ_scalar      usum;
  univ_mat         cast_mat;
  boolean          alloc_cast;
  double          *data;
  sint32           idx;
  sint32           nelem;

  double           num_ops = 0;
  double           tmp_variance;
  double           tmp_mean;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start  matuniv_mean_variance()" );

  if( !mean && !variance ) {
    MUTIL_ERROR( "NULL pointers for all return values" );
    return MUTIL_ERR_NULL_POINTER;
  }

  err = matuniv_validate( mat );
  if( err ) {
    MUTIL_ERROR( "Invalid input matrix argument" );
    return err;
  }

  /* complex numbers not supported yet */

  if( mat->type == MUTIL_DCOMPLEX ) {
    MUTIL_ERROR( "Unsupported matrix type" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }


  /* cast input to doubles and initialize univ scalar or mean result */
  type = MUTIL_DOUBLE;
  MUTIL_DO_STANDARD_CASTING( mat, type, cast_mat, intrp_ptr,
    alloc_cast, err );
  if( err ) {
    return err;
  }
  SCAUNIV_INIT( usum, MUTIL_DOUBLE, 0.0 );

  /* compute the mean */

  err = matuniv_sum( &cast_mat, intrp_ptr, &usum );
  if( err ) {
    MUTIL_FREE_STANDARD_CASTING( cast_mat, alloc_cast );
    return err;
  }
  nelem    = MATUNIV_NELEM( mat );
  tmp_mean = usum.num.dbl / (double) nelem;

  /* compute variance */
  data         = (double*) MATUNIV_DATA( &cast_mat );
  tmp_variance = 0.0;
  for( idx = 0; idx < nelem; idx++ ) {
    tmp_variance += ( data[idx] - tmp_mean ) *
      ( data[idx] - tmp_mean );
  }
  /* don't need cast data any more */
  MUTIL_FREE_STANDARD_CASTING( cast_mat, alloc_cast );

  num_ops += (nelem * 5.0);
  if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
    MUTIL_ERROR( "User Interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  /* finish up variance computation */

  if( unbiased == TRUE ) {
    if( nelem == 1 ) { /* avoid divide by zero */
      tmp_variance = 0;
    }
    else {
      tmp_variance /= (double) (nelem - 1);
    }
  }
  else {
    tmp_variance /= (double) (nelem);
  }

  if( mean ) {
    *mean = tmp_mean;
  }

  if( variance ) {
    *variance = tmp_variance;
  }

  MUTIL_TRACE( "matuniv_mean_variance() done" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_stat.h */
/* written by Jennifer Hodgdon */
mutil_errcode matuniv_quantiles( const univ_mat *mat,
  const univ_mat *probs, void *intrp_ptr, univ_mat *quantiles )
{
  mutil_errcode trouble;
  sint32_mat    index;
  sint32_mat    partials;
  univ_mat      uindex;
  univ_mat      upartials;
  double        trueind;
  double        weight;
  sint32        i;
  sint32        nprobs;
  sint32        mat_nelem;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE("Start matuniv_quantiles()");

  /* sanity checks */

  trouble = matuniv_validate( mat );
  if( trouble ) return trouble;
  trouble = matuniv_validate( probs );
  if( trouble ) return trouble;
  trouble = matuniv_validate( quantiles );
  if( trouble ) return trouble;

  if( probs->type != MUTIL_DOUBLE ||
    quantiles->type != MUTIL_DOUBLE ) {
    MUTIL_ERROR( "Input probability and output quantile matrices must be type double" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( !MATANY_EQUAL_DIM( &(probs->mat.dblmat), &(quantiles->mat.dblmat) )) {
    MUTIL_ERROR( "Input probability and output quantile matrices must be same size" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* Quantiles are calculated by doing a partial sort, where we
     need the indexed elements right around the probabilities */

  nprobs    = probs->mat.dblmat.nelem;
  mat_nelem = MATUNIV_NELEM( mat );

  trouble = mats32_malloc( &partials, nprobs * 2, 1 );
  if( trouble ) return trouble;

  for( i = 0; i < nprobs; i++ ) {
    if( probs->mat.dblmat.data[i] < 0.0 || probs->mat.dblmat.data[i] > 1.0 ) {
      MUTIL_FREE_WARN( mats32, &partials );
      MUTIL_ERROR( "Probabilities must be between 0 and 1" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

    partials.data[2 * i] = (sint32) floor( probs->mat.dblmat.data[i] *
      (mat_nelem - 1));
    partials.data[2 * i + 1] = partials.data[2 * i] + 1;
  }

  for( i = 0; i < partials.nelem; i++ ) {
    partials.data[i] = MUTIL_MIN( mat_nelem - 1, partials.data[i] );
    partials.data[i] = MUTIL_MAX( 0, partials.data[i] );
  }

  /* call sort to find the elements we need */

  trouble = mats32_malloc( &index, MATUNIV_NROW( mat ), MATUNIV_NCOL( mat ));
  if( trouble ) {
    MUTIL_FREE_WARN( mats32, &partials );
    return trouble;
  }

  trouble = matuniv_wrap_matrix( &uindex, &index, MUTIL_SINT32 );
  if( trouble ) {
    MUTIL_FREE_WARN( mats32, &partials );
    return trouble;
  }

  trouble = matuniv_wrap_matrix( &upartials, &partials, MUTIL_SINT32 );
  if( trouble ) {
    MUTIL_FREE_WARN( mats32, &partials );
    return trouble;
  }

  /* partial sort takes care of calling interrupt */
  trouble = matuniv_sort_index_partial( mat, &upartials, intrp_ptr, &uindex );
  if( trouble ) {
    MUTIL_FREE_WARN( mats32, &partials );
    MUTIL_FREE_WARN( mats32, &index );
    return trouble;
  }

  /* calculate the quantiles */

  for( i = 0; i < nprobs; i++ ) {
    /* each quantile is a weighted average of two data values on either side
       of the probability */

    trueind = probs->mat.dblmat.data[i] * ( mat_nelem - 1 );
    weight  = trueind - partials.data[ 2 * i ];

    quantiles->mat.dblmat.data[i] =
      ( 1.0 - weight ) *
      MATUNIV_ELEM( mat, 0, index.data[ partials.data[ 2 * i ]]) +
      (weight) *
      MATUNIV_ELEM( mat, 0, index.data[ partials.data[ 2 * i + 1 ]]);
  }

  MUTIL_FREE_WARN( mats32, &index );
  MUTIL_FREE_WARN( mats32, &partials );

  if( MUTIL_INTERRUPT( 12.0 * nprobs, intrp_ptr )) {
    MUTIL_ERROR( "User interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE("matuniv_quantiles() done");
  return MUTIL_ERR_OK;
}


/* function documented in mat_stat.h */
/* written by Andrea Borning */
mutil_errcode matset_histogram( const mat_set *matset, const double start_val,
  const double end_val, const boolean include_end, void *intrp_ptr,
  univ_mat *univhistogram )
{
  mutil_errcode   errcode;
  univ_mat        hist_tmp;
  univ_mat*       mat_tmp;
  univ_scalar     zero;
  sint32          index;

  MUTIL_TRACE( "Start matset_histogram()" );

  /* avoid lint warning */
  (void) whatssi;

  if( !matset || !univhistogram ) {
    MUTIL_ERROR( "NULL pointer for operand");
    return MUTIL_ERR_NULL_POINTER;
  }

  errcode = matset_verify_allsame( matset );
  if( errcode ) return errcode;

  if( univhistogram->type != MUTIL_UINT32 ) {
    MUTIL_ERROR( "Universal matrix for histogram must be of type "
      "MUTIL_UINT32" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* create a temporary histogram to hold the histogram of
     the matrix being considered from the matrix set */

  errcode = matuniv_malloc( &hist_tmp, univhistogram->mat.u32mat.nrow,
    univhistogram->mat.u32mat.ncol, MUTIL_UINT32 );
  if( errcode ) {
    MUTIL_FREE_WARN( matuniv, &hist_tmp );
    return errcode;
  }

  /* set zero to 0 */
  zero.type = MUTIL_UINT32;
  zero.num.u32 = ( uint32 ) 0;

  /* init output histogram to zeros */
  errcode = matuniv_assign_scalar( zero, intrp_ptr, univhistogram );
  if( errcode ) {
    MUTIL_FREE_WARN( matuniv, &hist_tmp );
    return errcode;
  }

  /* consider each matrix in the input */
  for( index = 0; index < matset->nelem; index++ ) {

    /* get a matrix from the input set */
    mat_tmp = &( matset->mats[ index ] );

    /* determine the histogram of the matrix under consideration */
    errcode = matuniv_histogram(
      mat_tmp, start_val, end_val, include_end, intrp_ptr, &hist_tmp );
    if( errcode ) {
      MUTIL_FREE_WARN( matuniv, &hist_tmp );
      return errcode;
    }

    /* add the histogram to the output histogram containing
       the cumulative sum of the individual histograms */
    errcode = matuniv_add( &hist_tmp, univhistogram, intrp_ptr,
      univhistogram );
    if( errcode ) {
      MUTIL_FREE_WARN( matuniv, &hist_tmp );
      return errcode;
    }
  }

  /* done with hist_tmp */
  MUTIL_FREE_WARN( matuniv, &hist_tmp );

  MUTIL_TRACE( "matset_histogram() done" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_stat.h */
/* written by Jennifer Hodgdon */
mutil_errcode matuniv_histogram( const univ_mat *mat, double start_val,
  double end_val, boolean include_end, void *intrp_ptr, univ_mat *histogram )
{
  mutil_errcode errcode;

  MUTIL_TRACE("Start matuniv_histogram()");

  if( !mat || !histogram ) {
    MUTIL_ERROR( "NULL pointer for operand");
    return MUTIL_ERR_NULL_POINTER;
  }

  if( histogram->type != MUTIL_UINT32 ) {
    MUTIL_ERROR( "Histogram matrix must be of type uint32" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch( mat->type ) {
    case MUTIL_DOUBLE:
      errcode = matdbl_histogram( &(mat->mat.dblmat), start_val,
        end_val, include_end, intrp_ptr, &(histogram->mat.u32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_histogram( &(mat->mat.fltmat), start_val,
        end_val, include_end, intrp_ptr, &(histogram->mat.u32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_histogram( &(mat->mat.u8mat), start_val,
        end_val, include_end, intrp_ptr, &(histogram->mat.u32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_histogram( &(mat->mat.u16mat), start_val,
        end_val, include_end, intrp_ptr, &(histogram->mat.u32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_histogram( &(mat->mat.u32mat), start_val,
        end_val, include_end, intrp_ptr, &(histogram->mat.u32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_histogram( &(mat->mat.s16mat), start_val,
        end_val, include_end, intrp_ptr, &(histogram->mat.u32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_histogram( &(mat->mat.s32mat), start_val,
        end_val, include_end, intrp_ptr, &(histogram->mat.u32mat) );
      if(errcode) return errcode;
      break;

    /* nothing else available now */
    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE("matuniv_histogram() done");
  return MUTIL_ERR_OK;
}


/** Template macro for histogram function.
 * Macro that expands to the body of a non-universal
 * histogram function, such as matdbl\_histogram.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_stat.c
 * @library matrix
 * @param FN_PREFIX      Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to take histogram of
 *    (function argument).
 * @param HIST_OUT_PTR   Matrix pointer to put histogram in
 *    (function argument).
 * @param START_VAL      Low endpoint of histogram (function argument).
 * @param END_VAL        High endpoint of histogram (function argument).
 * @param INCLUDE_END    Include or not include endpoint (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_histogram function:
 *     #TMPL_MAT_HIST( dbl, mat, histogram, start_val, end_val, include_end, intrp_ptr );#
 * @private
 */
#define TMPL_MAT_HIST( FN_PREFIX, MAT_IN_PTR, HIST_OUT_PTR, START_VAL,\
  END_VAL, INCLUDE_END, INTRP_PTR )\
  mutil_errcode errcode; \
  sint32        numbins; \
  sint32        indx; \
  sint32        ndat; \
  sint32        binnum; \
  double        binmult; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE("Start mat" #FN_PREFIX "_histogram()"); \
  \
  /* sanity checks */ \
  \
  errcode = mat ## FN_PREFIX ## _validate( MAT_IN_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = matu32_validate( HIST_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  if( (HIST_OUT_PTR)->ncol != 1 ) { \
    MUTIL_ERROR( "Histogram must have exactly one column"); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  if( START_VAL >= END_VAL ) { \
    MUTIL_ERROR( "Start value for histogram must be below end value"); \
    return MUTIL_ERR_ILLEGAL_VALUE; \
  } \
  \
  /* initialize bins to zero */ \
  \
  errcode = matu32_assign_scalar( 0, INTRP_PTR, HIST_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  /* calculate histogram */ \
  \
  numbins = (HIST_OUT_PTR)->nrow; \
  ndat    = (MAT_IN_PTR)->nelem; \
  binmult = numbins / (double) (END_VAL - START_VAL ); \
  \
  for( indx = 0; indx < ndat; indx++ ) { \
    /* check first end point */ \
    if( (MAT_IN_PTR)->data[indx] == START_VAL ) { \
      /* separate for roundoff */ \
      if( INCLUDE_END ) { \
        binnum = 0; \
      } \
      else { \
        continue; \
      } \
    } \
    \
    /* If data values are outside the start-end range, skip them. This */ \
    /* prevents subtractions of START_VAL from potentially large numbers, */ \
    /* which can cause floating point exceptions.  [botoole 3/24/99] */ \
    else if ( ( (MAT_IN_PTR)->data[indx] < START_VAL ) || \
                ( (MAT_IN_PTR)->data[indx] > END_VAL ) ) { \
      continue; \
    } \
    else { \
       /* if we are on a bin boundary, want to make the count go \
          to the lower bin */ \
      binnum =  (sint32) ceil(( (MAT_IN_PTR)->data[indx]  - START_VAL ) \
        * binmult ) - 1; \
      \
      /* roundoff error can occur when computing binmult, leading to error \
         when trying to place the maximum value */ \
      if ( binnum == numbins ) binnum -= 1; \
    } \
    \
    /* I was a bit concerned about the ++ below overflowing, but it can't */ \
    /* nelem is always s32 and so total count must be less than u32 */ \
    \
    if( binnum >= 0 && binnum < numbins ) { \
      ( (HIST_OUT_PTR)->data[binnum] )++; \
    } \
  } \
  \
  if( MUTIL_INTERRUPT( 4.0 * ndat, INTRP_PTR )) { \
    MUTIL_ERROR( "User interrupt" ); \
    return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE("mat" #FN_PREFIX "_histogram() done"); \
  \
  return MUTIL_ERR_OK


/* function documented in mat_stat.h */
/* written by Jennifer Hodgdon */
mutil_errcode matdbl_histogram( const double_mat *mat, double start_val,
  double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram )
{
  TMPL_MAT_HIST( dbl, mat, histogram, start_val, end_val, include_end,
    intrp_ptr );
}


/* function documented in mat_stat.h */
/* written by Andrea Borning */
mutil_errcode matflt_histogram( const float_mat *mat, double start_val,
  double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram )
{
  TMPL_MAT_HIST( flt, mat, histogram, start_val, end_val, include_end,
    intrp_ptr );
}


/* function documented in mat_stat.h */
/* written by Jennifer Hodgdon */
mutil_errcode matu8_histogram( const uint8_mat *mat, double start_val,
  double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram )
{
  TMPL_MAT_HIST( u8, mat, histogram, start_val, end_val, include_end,
    intrp_ptr );
}


/* function documented in mat_stat.h */
/* written by Jennifer Hodgdon */
mutil_errcode matu16_histogram( const uint16_mat *mat, double start_val,
  double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram )
{
  TMPL_MAT_HIST( u16, mat, histogram, start_val, end_val, include_end,
    intrp_ptr );
}


/* function documented in mat_stat.h */
/* written by Andrea Borning */
mutil_errcode matu32_histogram( const uint32_mat *mat, double start_val,
  double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram )
{
  TMPL_MAT_HIST( u32, mat, histogram, start_val, end_val, include_end,
    intrp_ptr );
}


/* function documented in mat_stat.h */
/* written by Andrea Borning */
mutil_errcode mats16_histogram( const sint16_mat *mat, double start_val,
  double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram )
{
  TMPL_MAT_HIST( s16, mat, histogram, start_val, end_val, include_end,
    intrp_ptr );
}


/* function documented in mat_stat.h */
/* written by Jill Goldschneider */
mutil_errcode mats32_histogram( const sint32_mat *mat, double start_val,
  double end_val, boolean include_end, void *intrp_ptr, uint32_mat *histogram )
{
  TMPL_MAT_HIST( s32, mat, histogram, start_val, end_val, include_end,
    intrp_ptr );
}
