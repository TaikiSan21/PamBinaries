
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_util.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */

#include "fra_util.h"

#include "mat_univ.h"
#include "mat_umat.h"
#include "mat_any.h"
#include "mat_set.h"
#include "mat_type.h"

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrp.h"

#include <math.h>

/*
This file contains utility functions
for fractal functions
The functions are declared in fra_util.h
*/

/* Static macro definitions */

#undef LOCALDEF_CHECK_NULL_POINTER_UTIL
#define LOCALDEF_CHECK_NULL_POINTER_UTIL( DATA_PTR,         \
  DATA_TYPE, TYPE_PREFIX )                                  \
  err = TYPE_PREFIX ## _validate( DATA_PTR );               \
  if ( err ) return err;                                    \
  if ( DATA_PTR == ( DATA_TYPE * ) NULL ) {                 \
  MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
  return MUTIL_ERR_NULL_POINTER;                            \
  }



mutil_errcode frauniv_check_embedding_inputs(
  const univ_mat     *data,
  sint32              dim,
  sint32              time_lag,
  sint32              orbital_lag,
  void               *intrp_ptr,
  boolean            *is_delay_embedding,
  sint32             *n_embed )
{
  mutil_errcode err;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_check_embedding_inputs()" );

  /* avoid lint warning */

  ( void ) whatssi;

  if ( MUTIL_INTERRUPT( 3.0, intrp_ptr ) ) {
    MUTIL_ERROR( "user interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  /*** check input data ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_UTIL( data, univ_mat, matuniv );

  /* ... if number of elements is positive */

  if ( MATUNIV_NELEM( data ) < 1 ){
    MUTIL_ERROR( "Number of elements in input data matrix must be positive." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if the type is double */

  if ( data->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Input data matrix must be of type double." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  /* ... to see if it is a vector */

  *is_delay_embedding = (boolean) ( MATANY_IS_VEC( &( data->mat.dblmat ) ) );

  /*** check embedding dimension ... ***/

  if ( dim <= 0 ){
    MUTIL_ERROR( "Embedding dimension must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  if ( !( *is_delay_embedding ) ){

    *n_embed = MATUNIV_NROW( data );

    if ( dim > MATUNIV_NCOL( data ) ){

      MUTIL_ERROR( "Specified embedding dimension is greater than "
	"number of columns in embedding matrix." );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

  }
  else{

    /*** check time lag ... ***/

    if ( time_lag <= 0 ){
      MUTIL_ERROR( "Time lag must be positive." );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

    if ( dim < 1 ){
      MUTIL_ERROR( "Embedding dimension must be positive." );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

    *n_embed  = MATUNIV_NELEM( data ) - ( dim - 1 ) * time_lag;

    if ( *n_embed < 1 ){
      MUTIL_ERROR( "Number of points in input data matrix not be at least "
	"n_sample - (dim - 1) * time_lag." );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
  }

  if ( orbital_lag < 0 ){
    MUTIL_ERROR( "Orbital lag must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( orbital_lag > (sint32) floor( (double) *n_embed / 2.0 ) ){
    MUTIL_ERROR( "Orbital lag is too large. "
      "Decrease time lag and/or embedding dimension" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  MUTIL_TRACE( "Done with frauniv_check_embedding_inputs()" );

  return MUTIL_ERR_OK;
}

mutil_errcode frauniv_check_partition(
  const mat_set *partition,
  double         scale )
{
  mutil_errcode err;

  MUTIL_TRACE( "Start frauniv_check_partition()" );

  LOCALDEF_CHECK_NULL_POINTER_UTIL( partition, mat_set, matset );

  if ( partition->nelem != 2 ){
    MUTIL_ERROR( "Partition matrix set must contain two matrices." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( !MATANY_IS_VEC( &(partition->mats[0].mat.s32mat) ) ||
    !MATANY_IS_VEC( &(partition->mats[1].mat.s32mat) ) ){

    MUTIL_ERROR( "Partition matrices in matrix set must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( partition->mats[ 0 ].type != MUTIL_SINT32 ){
    MUTIL_ERROR( "The matrices in the partition "
      "matrix set must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( scale < 0.0 ){
    MUTIL_ERROR( "Scale must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  MUTIL_TRACE( "Done with frauniv_check_partition()" );

  return MUTIL_ERR_OK;
}
