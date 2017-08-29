
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/wav_wtmm.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "wav_wtmm.h"

#include "mat_any.h"
#include "mat_assn.h"
#include "mat_comp.h"
#include "mat_set.h"
#include "mat_summ.h"
#include "mat_type.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_math.h"
#include "ut_mem.h"
#include <math.h>
#include <stdio.h>

/* Static functions declared here and defined at end of file */

/* Static macro definitions */

#define LOCALDEF_WTMM( CWT, ROW, COL ) \
 (double) MUTIL_CPX_ABS( (CWT)->data[ (ROW) * (CWT)->ncol + (COL) ] )

#define LOCALDEF_CHECK_NULL_POINTER( DATA_PTR, DATA_TYPE,      \
                                     TYPE_PREFIX )             \
   err = TYPE_PREFIX ## _validate( DATA_PTR );                 \
   if ( err ) return err;                                      \
   if ( DATA_PTR == ( DATA_TYPE * ) NULL ){                    \
     MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
     return MUTIL_ERR_NULL_POINTER;                            \
   }

#define LOCALDEF_CWT_MODULUS( TIME, SCALE ) \
MUTIL_CPX_ABS( cwt->mat.cpxmat.data[(TIME) * n_scale + (SCALE)] )

#define LOCALDEF_CWT_REAL( TIME, SCALE ) \
cwt->mat.cpxmat.data[(TIME) * n_scale + (SCALE)].re

#define LOCALDEF_REVERSE_VECTOR( MAT_PTR, MAT_FN_TYPE ) \
 if ( MATANY_IS_VEC_COL( MAT_PTR ) ){ \
   err = mat ## MAT_FN_TYPE ## _flip_up_down( MAT_PTR, intrp_ptr, MAT_PTR ); \
 }\
 else{ \
   err = mat ## MAT_FN_TYPE ## _flip_left_right( MAT_PTR, intrp_ptr, MAT_PTR ); \
 }\
 MEMLIST_FREE_ON_ERROR( err, &list )

#define LOCALDEF_UPDATE_BRANCH( WTMM_TIME, WTMM_SCALE, INDEX_DELETE ) \
 /* ... update branch lists */\
 \
 err = mats32_realloc_register( &branch_length, n_branch, 1, &list ); \
 MEMLIST_FREE_ON_ERROR( err, &list ); \
 \
 branch_time.data[ count ] = WTMM_TIME; \
 branch_scale.data[ count ] = WTMM_SCALE; \
 branch_length.data[ n_branch - 1 ]++; \
 \
 /* ... finally, set a flag in the temporary WTMM vectors \
    indicating that we have processed the current modulus maxima */ \
 \
 wtmm_time_index_temp.data[ INDEX_DELETE ]  = -1; \
 wtmm_scale_index_temp.data[ INDEX_DELETE ] = -1; \
 \
 count++

#define LOCALDEF_UPDATE_WTMM_NEIGHBOR_CANDIDATE \
 if ( diff < min ){ \
   min = diff; \
   idelete = match_index.data[ t ]; \
   ttdata = tt.data[ t ]; \
 }

/* The modulus maxima of a continuous wavelet transform   */
/* Documented in wav_wtmm.h                                */
/* Written by William Constantine                         */

mutil_errcode wavuniv_transform_continuous_wavelet_modulus_maxima(
  const univ_mat        *cwt,
  const univ_mat        *tolerance,
  const wav_transform_peak peak_type,
  void                  *intrp_ptr,
  sint32_mat            *itime,
  sint32_mat            *iscale )
{
  memlist       list;
  mutil_errcode err;
  sint32        t;
  sint32        j;
  sint32        n_sample;
  sint32        n_scale;
  double_mat    modcwt;
  double        max = 0.0;
  double        atol;
  sint32        imax;
  boolean       up;
  sint32        plateau;
  sint32        count;
  double        cwtreal;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_transform_continuous_wavelet_modulus_maxima()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  /*** check cwt ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER( cwt, univ_mat, matuniv );

  /* ... if the type is complex double */

  if ( cwt->type != MUTIL_DCOMPLEX ){
    MUTIL_ERROR( "CWT matrix must be of type MUTIL_DCOMPLEX." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* obtain dimensions of CWT matrix */

  n_scale  = MATUNIV_NCOL( cwt );
  n_sample = MATUNIV_NROW( cwt );

  /* ... for minimum dimensions */

  if ( n_sample < 2 ){
    MUTIL_ERROR( "Time coordinates of CWT matrix must be at least 2 elements long." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /*** check tolerance ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER( tolerance, univ_mat, matuniv );

  if ( tolerance->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Tolerance matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( !MATANY_IS_VEC( &(tolerance->mat.dblmat) ) ){
    MUTIL_ERROR( "Tolerance matrix must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( MATUNIV_NELEM( tolerance ) != n_scale ){
    MUTIL_ERROR( "Tolerance matrix must contain the same number of elements " \
      "as there are columns in the CWT input matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /*** check peak type ... ***/
  switch( peak_type ){
    case WAV_TRANSFORM_PEAK_EXTREMA:
    case WAV_TRANSFORM_PEAK_MAXIMA:
    case WAV_TRANSFORM_PEAK_MINIMA:
      break;
    default:
      MUTIL_ERROR( "Wavelet transform peak type is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* allocate space for temporary matrices */

  err = matdbl_malloc_register( &modcwt, n_sample, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* allocate space for output */

  err = mats32_malloc_register( itime, 1, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( iscale, 1, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* initialize the count of the WTMM found */

  count = 0;

  /* create the WTMM, a binary mask */

  for ( j = 0; j < n_scale; j++ ){

    /* copy the modulus of the current scale's CWT coefficients into a vector */

    for ( t = 0; t < n_sample; t++ ){
      switch( peak_type ){
      case WAV_TRANSFORM_PEAK_EXTREMA:
        modcwt.data[t] = LOCALDEF_CWT_MODULUS( t, j );
        break;
      case WAV_TRANSFORM_PEAK_MAXIMA:
        cwtreal = LOCALDEF_CWT_REAL( t, j );
        modcwt.data[t] = ( ( cwtreal > 0.0 ) ? cwtreal : (double) 0.0 );
        break;
      case WAV_TRANSFORM_PEAK_MINIMA:
        cwtreal = LOCALDEF_CWT_REAL( t, j );
        modcwt.data[t] = ( ( cwtreal < 0.0 ) ? MUTIL_ABS(cwtreal) : (double) 0.0 );
        break;
      }
    }

    /* initialize variables */

    atol    = tolerance->mat.dblmat.data[j];
    up      = (boolean) ( modcwt.data[1] < ( modcwt.data[0] + atol ) );
    max     = modcwt.data[0];
    imax    = 0;
    plateau = 0;

    for ( t = 1; t < n_sample; t++ ){

      if ( ( ( modcwt.data[t] + atol ) < max )  && ( max > atol ) && up ){

	/* grow wtmm lists */

	err = mats32_realloc_register( itime, count + 1, 1, &list );
	MEMLIST_FREE_ON_ERROR( err, &list );

	err = mats32_realloc_register( iscale, count + 1, 1, &list );
	MEMLIST_FREE_ON_ERROR( err, &list );

	/* assign current WTMM data to lists */

	itime->data[ count ]  = imax - (sint32) floor( (double) plateau / 2.0 );
	iscale->data[ count ] = j;

	/* update variables */

	count++;
	plateau = 0;
	max     = -1.0;
      }

      if ( ( MUTIL_ABS( modcwt.data[t] - modcwt.data[t-1] ) < atol )  && up ){
	up = (boolean) TRUE;
	plateau++;
      }
      else{
	up = (boolean) ( modcwt.data[t] > ( modcwt.data[t - 1] + atol ) );
      }

      if ( up ){
	max = modcwt.data[t];
	imax = t;
      }

      /* Check for interrupts */

      if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
	MUTIL_ERROR( "user interrupt" );
	MUTIL_FREE_WARN( memlist, &list );
	return MUTIL_ERR_INTERRUPT;
      }

    } /* end loop over time */

  } /* end loop over scales */


  if ( count == 0 ){

    MUTIL_FREE_WARN( memlist, &list );

    itime->nelem = (sint32) 0;
    itime->nrow  = (sint32) 0;
    itime->ncol  = (sint32) 0;

    iscale->nelem = (sint32) 0;
    iscale->nrow  = (sint32) 0;
    iscale->ncol  = (sint32) 0;
  }
  else{

    /* free nodes corresponding to output
       in memory list, but not the memory itself as
       it hust remain to send back to the caller */

    err = memlist_member_unregister( itime, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_unregister( iscale, &list );
    if ( err ){
      MUTIL_FREE_WARN( mats32, itime );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    MUTIL_FREE_WARN( memlist, &list );
  }

  MUTIL_TRACE( "Done with wavuniv_transform_continuous_wavelet_modulus_maxima()" );

  return MUTIL_ERR_OK;
}

/* The tree of the modulus maxima of */
/* a continuous wavelet transform        */
/* Documented in wav_wtmm.h              */
/* Written by William Constantine        */

mutil_errcode wavuniv_transform_continuous_wavelet_modulus_maxima_tree(
  const sint32_mat   *wtmm_time_index,
  const sint32_mat   *wtmm_scale_index,
  const dcomplex_mat *cwt,
  const double_mat   *cwt_time,
  const double_mat   *cwt_scale,
  const boolean       bridge_gaps,
  const double        n_octave_min,
  const double        wtmm_strength_min,
  void               *intrp_ptr,
  mat_set            *tree )
{
  boolean       reverse;
  double        time_search_limit;
  double        wtmm;
  double        wtmm_max;
  sint32        n_scale;
  double        n_voice;
  double        n_octave;
  memlist       list;
  mutil_errcode err;
  sint32        count;
  sint32        diff;
  sint32        dims;
  sint32        ibranch;
  sint32        idelete;
  sint32        ifound;
  sint32        n_scale_min;
  sint32        iscale;
  sint32        iscale_max;
  sint32        iscale_min;
  sint32        itime_max;
  sint32        itime_min;
  sint32        itime;
  sint32        min;
  sint32        n_branch;
  sint32        n_found;
  sint32        n_maxima;
  sint32        n_sample;
  sint32        pointer;
  sint32        t;
  sint32        i;
  sint32        ttdata;
  sint32_mat    branch_length;
  sint32_mat    branch_scale;
  sint32_mat    branch_time;
  sint32_mat    match_index;
  sint32_mat    branch_length_long;
  sint32_mat    nrow;
  sint32_mat    tt;
  sint32_mat    wtmm_scale_index_temp;
  sint32_mat    wtmm_time_index_temp;
  double        sampling_interval;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_transform_continuous_wavelet_modulus_maxima_tree()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check arguments ... */

  /* ... for NULL pointers */

  LOCALDEF_CHECK_NULL_POINTER( wtmm_time_index, sint32_mat, mats32 );
  LOCALDEF_CHECK_NULL_POINTER( wtmm_scale_index, sint32_mat, mats32 );
  LOCALDEF_CHECK_NULL_POINTER( cwt_time, double_mat, matdbl );
  LOCALDEF_CHECK_NULL_POINTER( cwt_scale, double_mat, matdbl );
  LOCALDEF_CHECK_NULL_POINTER( cwt, dcomplex_mat, matcpx );

  /* ... for vectors */

  if ( !MATANY_IS_VEC( wtmm_time_index ) ||
    !MATANY_IS_VEC( wtmm_scale_index ) ||
    !MATANY_IS_VEC( cwt_scale ) ||
    !MATANY_IS_VEC( cwt_time ) ){
    MUTIL_ERROR( "WTMM time and scale index matrices as well as CWT time " \
      "and scale matrices must be vectors." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... for consistent dimensions */

  if ( wtmm_time_index->nelem != wtmm_scale_index->nelem ||
    wtmm_time_index->nelem  < 1 ||
    wtmm_scale_index->nelem < 1 ){
    MUTIL_ERROR( "Number of elements in input time index " \
      "and scale index matrices must be the same and positive." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( cwt->nrow != cwt_time->nelem ||  cwt->ncol != cwt_scale->nelem ){
    MUTIL_ERROR( "CWT matrix dimensions are inconsistent with CWT time " \
      "and/or CWT scale vectors." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ...for sufficiently sized dimensions */

  if ( cwt->nrow < 2 || cwt->ncol < 2 ){
    MUTIL_ERROR( "CWT matrix must be at least a 2 x 2." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... for consistent index ranges */

  err = mats32_range( wtmm_scale_index, intrp_ptr, &iscale_min, &iscale_max );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_range( wtmm_time_index, intrp_ptr, &itime_min, &itime_max );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( itime_max > cwt->nrow || iscale_max > cwt->ncol ){
    MUTIL_ERROR( "Index range in either the WTMM scale index or time " \
      "index vectors exceeds the corresponding dimensions of the input CWT matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* .. for acceptable pruning factors */

  if ( n_octave_min < 0.0 ){
    MUTIL_ERROR( "Minimum number of ocatves to span must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( wtmm_strength_min < 0.0 || wtmm_strength_min > 1.0 ){
    MUTIL_ERROR( "WTMM strength pruning factor must be in the range of [0,1]." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* intialize variables */

  n_branch = 0;
  count    = 0;
  n_scale  = cwt_scale->nelem;
  n_sample = cwt_time->nelem;
  n_maxima = wtmm_time_index->nrow;
  sampling_interval = (double) MUTIL_ABS( cwt_time->data[ 1 ] - cwt_time->data[ 0 ] );

  reverse  = (boolean) ( cwt_scale->data[ 0 ] < cwt_scale->data[ n_scale - 1 ] );

  /* calculate the minimum number of scales that must comprise a WTMM
     branch. lengths which are less than this value are removed
     (pruned) fomr the branch list */

  if ( reverse ){
    n_octave = (double) MUTIL_LOG2( cwt_scale->data[ n_scale - 1 ] / cwt_scale->data[ 0 ] );
  }
  else{
    n_octave = (double) MUTIL_LOG2( cwt_scale->data[ 0 ] / cwt_scale->data[ n_scale - 1 ] );
  }

  n_voice     = (double) ( n_scale - 1 ) / n_octave;
  n_scale_min = (sint32) floor( n_voice * n_octave_min );

  /* allocate memory for branch storage */

  err = mats32_malloc_register( &branch_time, 1, n_maxima, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &branch_scale, 1, n_maxima, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &branch_length, 1, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* copy WTMM input vectors into temporary vectors for manipulation */

  err = mats32_malloc_register( &wtmm_time_index_temp, wtmm_time_index->nrow, wtmm_time_index->ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_assign( wtmm_time_index, intrp_ptr, &wtmm_time_index_temp );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &wtmm_scale_index_temp, wtmm_scale_index->nrow, wtmm_scale_index->ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_assign( wtmm_scale_index, intrp_ptr, &wtmm_scale_index_temp );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* if necessary, reverse the scale and corresponding input
     vectors so that the scale is ordered from large to small */

  if ( reverse ){
    LOCALDEF_REVERSE_VECTOR( &wtmm_time_index_temp, s32 );
    LOCALDEF_REVERSE_VECTOR( &wtmm_scale_index_temp, s32 );
  }

  /* prune weak valued WTMM from the branch list */

  for ( iscale = iscale_min; iscale <= iscale_max; iscale++ ){

    /* locate the WTMM at the current scale */

    err = mats32_compare_scalar( &wtmm_scale_index_temp, MUTIL_RELATION_EQUAL, iscale,
      intrp_ptr, &match_index, (sint32_mat *) NULL );
    if ( err ){
      MUTIL_FREE_WARN( mats32, &match_index );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    /* find the maximum WTMM at the current scale */

    wtmm_max = 0.0;

    for ( i = 0; i < match_index.nelem; i++ ){

      wtmm = LOCALDEF_WTMM( cwt, wtmm_time_index_temp.data[ match_index.data[ i ] ], iscale );

      if ( wtmm > wtmm_max ){
	wtmm_max = wtmm;
      }
    }

    /* sift through each WTMM at the current scale,
       and if it is not sufficiently strong, remove it
       from the branch lists */

    for ( i = 0; i < match_index.nelem; i++ ){

      wtmm = LOCALDEF_WTMM( cwt, wtmm_time_index_temp.data[ match_index.data[ i ] ], iscale );

      if ( wtmm <= wtmm_max * wtmm_strength_min ){

	wtmm_time_index_temp.data[ match_index.data[ i ] ]  = -1;
	wtmm_scale_index_temp.data[ match_index.data[ i ] ] = -1;
      }
    }

    /* free the match index vector for the next loop */

    if ( match_index.nelem > 0 ){
      MUTIL_FREE_WARN( mats32, &match_index );
    }
  }

  /* start branching the WTMM */

  while ( count < wtmm_scale_index->nelem ){

    /* Check for interrupts */

    if ( MUTIL_INTERRUPT( 3.0 * count, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }

    /* start a new branch ... */

    /* ... find the index locations of the maxima we have not yet processed */

    err = mats32_compare_scalar( &wtmm_time_index_temp,
      MUTIL_RELATION_GREATER_THAN_OR_EQUAL, (sint32) 0,
      intrp_ptr, &match_index, (sint32_mat *) NULL );
    if ( err ){
      MUTIL_FREE_WARN( mats32, &match_index );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    if ( match_index.nelem > 0 ){

      /* ... record the first maxima found (out of the pool of unprocessed maxima)
	 as the beginning of the new branch */

      ifound = match_index.data[ 0 ];
      itime  = wtmm_time_index_temp.data[ ifound ];
      iscale = wtmm_scale_index_temp.data[ ifound ];

      /* ... destroy the memory of the match index for a fresh start
	 on the next loop */
      MUTIL_FREE_WARN( mats32, &match_index );

      /* ... update branch lists */

      n_branch++;

      LOCALDEF_UPDATE_BRANCH( itime, iscale, ifound );

			//begin Lixin Gong
			/* initialize the length of the new branch */
			branch_length.data[n_branch-1] = 1;
			//end Lixin Gong

      /* pursue the rest of branch */

      while( iscale > iscale_min ){

	/* Check for interrupts */

	if ( MUTIL_INTERRUPT( 3.0 * count, intrp_ptr ) ) {
	  MUTIL_ERROR( "user interrupt" );
	  MUTIL_FREE_WARN( memlist, &list );
	  return MUTIL_ERR_INTERRUPT;
	}

	/* advance to the next scale (coarse to fine) */

	iscale--;

	/* find all of the WTMM recorded at the current scale */

	err = mats32_compare_scalar( &wtmm_scale_index_temp, MUTIL_RELATION_EQUAL, iscale,
	  intrp_ptr, &match_index, (sint32_mat *) NULL );
	if ( err ){
	  MUTIL_FREE_WARN( mats32, &match_index );
	  MEMLIST_FREE_ON_ERROR( err, &list );
	}

	n_found = match_index.nelem;

	/* if any WTMM were found at the next level, update
	   the current branch with the candidate closest in
	   (circular) time to the last recorded member in the
	   WTMM branch */

	if ( n_found > 0 ){

	  /* grab the corresponding indices of time */

	  err = mats32_malloc_register( &tt, 1, n_found, &list );
	  MEMLIST_FREE_ON_ERROR( err, &list );

	  for ( t = 0; t < n_found; t++ ){

	    tt.data[ t ] = wtmm_time_index_temp.data[ match_index.data[ t ] ];

	    if ( t == 0 ){

	      ttdata  = tt.data[ 0 ];
	      idelete = match_index.data[ 0 ];
	      min     = (sint32) MUTIL_ABS( tt.data[ t ] - itime );
	    }
	    else{

	      diff = (sint32) MUTIL_ABS( tt.data[ t ] - itime );
	      LOCALDEF_UPDATE_WTMM_NEIGHBOR_CANDIDATE;
	    }

	    diff = (sint32) MUTIL_ABS( tt.data[ t ] - itime - n_sample );
	    LOCALDEF_UPDATE_WTMM_NEIGHBOR_CANDIDATE;

	    diff = (sint32) MUTIL_ABS( tt.data[ t ] - itime + n_sample );
	    LOCALDEF_UPDATE_WTMM_NEIGHBOR_CANDIDATE;
	  }

	  /* clear allocated memory for next loop */

	  MUTIL_FREE_WARN( mats32, &match_index );

	  err = memlist_member_free( &tt, &list);
	  MEMLIST_FREE_ON_ERROR( err, &list );

	  time_search_limit = cwt_scale->data[ iscale ];

	  if ( ( (double) MUTIL_ABS( itime - ttdata ) * sampling_interval ) < time_search_limit ||
	    (double) MUTIL_ABS( itime - ttdata - n_sample ) < time_search_limit ||
	    (double) MUTIL_ABS( itime - ttdata + n_sample ) < time_search_limit ){

	    itime = ttdata;

	    LOCALDEF_UPDATE_BRANCH( itime, iscale, idelete );
	  }
	  else if ( !bridge_gaps ){
	    break;
	  }
	}
	else if ( !bridge_gaps ){

	  /* we have reached a gap in the branch
	     but do not wish to jump over the gap.
	     in this case, we break out of the inner
	     while loop used to seek the remainder of
	     the current branch */
	  break;
	}
      }

    }
    else{
      break;
    }
  }

  /* map the branches into a matrix set ... */

  /* ... identify those branches that are sufficiently long */

  err = mats32_compare_scalar( &branch_length,
    MUTIL_RELATION_GREATER_THAN, n_scale_min,
    intrp_ptr, &match_index, &branch_length_long );
  if ( err ){
    MUTIL_FREE_WARN( mats32, &match_index );
    MUTIL_FREE_WARN( mats32, &branch_length_long );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* ... allocate memory for matrix set and
     coresponding matrices */

  dims = branch_length_long.nelem;

  if ( dims > 0 ){

    err = mats32_malloc_register( &nrow, 1, dims, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = mats32_assign_scalar( (sint32) 5, intrp_ptr, &nrow );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matset_malloc_register( tree, (sint32) 1, &dims, &list );
    if ( err ) return err;

    err = matset_malloc_matrices_arbitrary_size(
      tree, &nrow, &branch_length_long, MUTIL_DOUBLE );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* ... fill the matrices in the matrix set */

    pointer = 0;

    i = 0;

    for ( ibranch = 0; ibranch < branch_length.nelem; ibranch++  ){

      if ( i < dims && ibranch == match_index.data[ i ] ){

	for ( t = 0; t < branch_length.data[ ibranch ]; t++ ){

	  itime  = branch_time.data[ pointer + t ];
	  iscale = branch_scale.data[ pointer + t ];
	  wtmm   = LOCALDEF_WTMM( cwt, itime, iscale );

	  MATUNIV_ELEM_ASSIGN( &( tree->mats[ i ] ), 0, t, (double) itime );
	  MATUNIV_ELEM_ASSIGN( &( tree->mats[ i ] ), 1, t, (double) iscale );
	  MATUNIV_ELEM_ASSIGN( &( tree->mats[ i ] ), 2, t, cwt_time->data[ itime ] );
	  MATUNIV_ELEM_ASSIGN( &( tree->mats[ i ] ), 3, t, cwt_scale->data[ iscale ] );
	  MATUNIV_ELEM_ASSIGN( &( tree->mats[ i ] ), 4, t, wtmm );
	}

	i++;
      }

      pointer += branch_length.data[ ibranch ];
    }
  }

  /* free memory used to find sufficiently long branches */

  MUTIL_FREE_WARN( mats32, &match_index );
  MUTIL_FREE_WARN( mats32, &branch_length_long );

  /* free nodes corresponding to output
     in memory list, but not the memory itself as
     it hust remain to send back to the caller */

  err = memlist_member_unregister( tree, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_transform_continuous_wavelet_modulus_maxima_tree()" );

  return MUTIL_ERR_OK;
}
