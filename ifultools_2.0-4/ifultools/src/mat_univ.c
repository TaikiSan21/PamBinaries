
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_univ.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "mat_univ.h"
#include "mat_umat.h"
#include "mat_tmpl.h"

#include "ut_debug.h"

/*
   This file contains the definitions for the functions for
   universal matrices declared in mat_univ.h.
*/


/* Function documented in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode matuniv_wrap_matrix( univ_mat *umat, void *matrix_ptr,
  mutil_data_type type )
{
  MUTIL_TRACE( "Start matuniv_wrap_matrix()" );

  /* avoid lint warning */
  (void) whatssi;

  if( !umat ) {
    MUTIL_ERROR( "NULL universal matrix pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if( !matrix_ptr ) {
    MUTIL_ERROR( "NULL matrix pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }

  umat->type = type;

  switch( type ) {
    case MUTIL_UINT8:
      umat->mat.u8mat = *( (uint8_mat *) matrix_ptr );
      break;
    case MUTIL_SINT8:
      umat->mat.s8mat = *( (sint8_mat *) matrix_ptr );
      break;
    case MUTIL_UINT16:
      umat->mat.u16mat = *( (uint16_mat *) matrix_ptr );
      break;
    case MUTIL_SINT16:
      umat->mat.s16mat = *( (sint16_mat *) matrix_ptr );
      break;
    case MUTIL_UINT32:
      umat->mat.u32mat = *( (uint32_mat *) matrix_ptr );
      break;
    case MUTIL_SINT32:
      umat->mat.s32mat = *( (sint32_mat *) matrix_ptr) ;
      break;
    case MUTIL_FLOAT:
      umat->mat.fltmat = *( (float_mat *) matrix_ptr );
      break;
    case MUTIL_DOUBLE:
      umat->mat.dblmat = *( (double_mat *) matrix_ptr );
      break;
    case MUTIL_DCOMPLEX:
      umat->mat.cpxmat = *( (dcomplex_mat *) matrix_ptr );
      break;
    default:
      MUTIL_ERROR( "Unknown value for data type" );
      return MUTIL_ERR_ILLEGAL_VALUE;
  }

  MUTIL_TRACE( "matuniv_wrap_matrix() done" );
  return MUTIL_ERR_OK;
}


/* This function is documented in mat_univ.h */
/* Written by Andrea Borning. */
mutil_errcode matuniv_wrap_univ_matrix( univ_mat *mat1, univ_mat *mat2 )
{
  mutil_errcode errcode;

  if( !mat1 || !mat2 ) {
    MUTIL_ERROR( "NULL universal matrix pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }

  MUTIL_TRACE( "Start matuniv_wrap_univ_matrix()" );
  switch ( mat2->type ) {
    case MUTIL_UINT8:
      errcode = matuniv_wrap_matrix( mat1, &(mat2->mat.u8mat ), mat2->type );
      break;
    case MUTIL_SINT8:
      errcode = matuniv_wrap_matrix( mat1, &(mat2->mat.s8mat ), mat2->type );
      break;
    case MUTIL_UINT16:
      errcode = matuniv_wrap_matrix( mat1, &(mat2->mat.u16mat), mat2->type );
      break;
   case MUTIL_SINT16:
      errcode = matuniv_wrap_matrix( mat1, &(mat2->mat.s16mat), mat2->type );
      break;
    case MUTIL_UINT32:
      errcode = matuniv_wrap_matrix( mat1, &(mat2->mat.s32mat), mat2->type );
      break;
    case MUTIL_SINT32:
      errcode = matuniv_wrap_matrix( mat1, &(mat2->mat.s32mat), mat2->type );
      break;
    case MUTIL_FLOAT:
      errcode = matuniv_wrap_matrix( mat1, &(mat2->mat.fltmat), mat2->type );
      break;
    case MUTIL_DOUBLE:
      errcode = matuniv_wrap_matrix( mat1, &(mat2->mat.dblmat), mat2->type );
      break;
    case MUTIL_DCOMPLEX:
      errcode = matuniv_wrap_matrix( mat1, &(mat2->mat.cpxmat), mat2->type );
      break;
    default:
      MUTIL_ERROR( "Unknown value for data type" );
      errcode = MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( errcode ) return errcode;

  MUTIL_TRACE( "matuniv_wrap_univ_matrix() done" );
  return MUTIL_ERR_OK;
}


/* Function documented in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matuniv_wrap_data( univ_mat *matrix, void *data,
  sint32 nrow, sint32 ncol, mutil_data_type type )
{
  mutil_errcode errcode;

  MUTIL_TRACE( "Start matuniv_wrap_data()" );

  if( !matrix ){
    MUTIL_ERROR( "NULL universal matrix pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if( !data ) {
    MUTIL_ERROR( "NULL data pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }

  matrix->type = type;

  switch( type ) {
    case MUTIL_UINT8:
      errcode = matu8_wrap_data( &(matrix->mat.u8mat), (uint8 *) data,
        nrow, ncol );
      break;

    case MUTIL_SINT8:
      errcode = MUTIL_ERR_OK;
      matrix->mat.s8mat.data  = (sint8 *) data;
      matrix->mat.s8mat.nrow  = nrow;
      matrix->mat.s8mat.ncol  = ncol;
      matrix->mat.s8mat.nelem = ncol * nrow;
      if( nrow <= 0 || ncol <= 0 ) {
        MUTIL_ERROR("Illegal matrix dimension value");
        errcode = MUTIL_ERR_ILLEGAL_SIZE;
      }
      if( nrow * ncol > MUTIL_SINT32_MAX ) {
        MUTIL_ERROR("Number of elements is too large");
        errcode = MUTIL_ERR_ILLEGAL_SIZE;
      }
      /* errcode = mats8_wrap_data( &(matrix->mat.fltmat), (float *) data,
        nrow, ncol ); */
      break;

    case MUTIL_UINT16:
      errcode = matu16_wrap_data( &(matrix->mat.u16mat), (uint16 *) data,
        nrow, ncol );
      break;

    case MUTIL_SINT16:
      errcode = mats16_wrap_data( &(matrix->mat.s16mat), (sint16 *) data,
        nrow, ncol );
      break;

    case MUTIL_UINT32:
      errcode = matu32_wrap_data( &(matrix->mat.u32mat), (uint32 *) data,
        nrow, ncol );
      break;

    case MUTIL_SINT32:
      errcode = mats32_wrap_data( &(matrix->mat.s32mat), (sint32 *) data,
        nrow, ncol );
      break;

    case MUTIL_FLOAT:
      errcode = matflt_wrap_data( &(matrix->mat.fltmat), (float *) data,
        nrow, ncol );
      break;

    case MUTIL_DOUBLE:
      errcode = matdbl_wrap_data( &(matrix->mat.dblmat), (double *) data,
        nrow, ncol );
      break;

    case MUTIL_DCOMPLEX:
      errcode = matcpx_wrap_data( &(matrix->mat.cpxmat), (dcomplex *) data,
        nrow, ncol );
      break;

    default:
      MUTIL_ERROR( "Unknown value for data type" );
      return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( errcode ) return errcode;

  MUTIL_TRACE( "matuniv_wrap_data() done" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_wrap_data in mat_univ.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matcpx_wrap_data( dcomplex_mat *matrix, dcomplex *data,
  sint32 nrow, sint32 ncol )
{
  TMPL_MAT_WRAP( matcpx, matrix, nrow, ncol, data );
}


/* This function is documented under matuniv_wrap_data in mat_univ.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_wrap_data( double_mat *matrix, double *data, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_WRAP( matdbl, matrix, nrow, ncol, data );
}


/* This function is documented under matuniv_wrap_data in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matflt_wrap_data( float_mat *matrix, float *data, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_WRAP( matflt, matrix, nrow, ncol, data );
}


/* This function is documented under matuniv_wrap_data in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode mats32_wrap_data( sint32_mat *matrix, sint32 *data,
  sint32 nrow, sint32 ncol )
{
  TMPL_MAT_WRAP( mats32, matrix, nrow, ncol, data );
}


/* This function is documented under matuniv_wrap_data in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats16_wrap_data( sint16_mat *matrix, sint16 *data,
  sint32 nrow, sint32 ncol )
{
  TMPL_MAT_WRAP( mats16, matrix, nrow, ncol, data );
}


/* This function is documented under matuniv_wrap_data in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode matu32_wrap_data( uint32_mat *matrix, uint32 *data, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_WRAP( matu32, matrix, nrow, ncol, data );
}


/* This function is documented under matuniv_wrap_data in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu16_wrap_data( uint16_mat *matrix, uint16 *data,
  sint32 nrow, sint32 ncol )
{
  TMPL_MAT_WRAP( matu16, matrix, nrow, ncol, data );
}


/* This function is documented under matuniv_wrap_data in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu8_wrap_data( uint8_mat *matrix, uint8 *data, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_WRAP( matu8, matrix, nrow, ncol, data );
}


/* Function documented in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode matuniv_malloc( univ_mat *matrix, sint32 nrow,
  sint32 ncol, mutil_data_type type )
{
  mutil_errcode errcode;

  MUTIL_TRACE( "Start matuniv_malloc()" );

  /* sanity checks... */

  if( !matrix ) {
    MUTIL_ERROR( "NULL pointer for returned matrix" );
    return MUTIL_ERR_NULL_POINTER;
  }

  matrix->type = type;

  switch( type ) {
    case MUTIL_DOUBLE:
      errcode = matdbl_malloc( &(matrix->mat.dblmat), nrow, ncol );
      break;

    case MUTIL_FLOAT:
      errcode = matflt_malloc( &(matrix->mat.fltmat), nrow, ncol );
      break;

    case MUTIL_UINT8:
      errcode = matu8_malloc( &(matrix->mat.u8mat), nrow, ncol );
      break;

    case MUTIL_UINT16:
      errcode = matu16_malloc( &(matrix->mat.u16mat), nrow, ncol );
      break;

    case MUTIL_UINT32:
      errcode = matu32_malloc( &(matrix->mat.u32mat), nrow, ncol );
      break;

    case MUTIL_SINT16:
      errcode = mats16_malloc( &(matrix->mat.s16mat), nrow, ncol );
      break;

    case MUTIL_SINT32:
      errcode = mats32_malloc( &(matrix->mat.s32mat), nrow, ncol );
      break;

    case MUTIL_DCOMPLEX:
      errcode = matcpx_malloc( &(matrix->mat.cpxmat), nrow, ncol );
      break;

    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( errcode ) return errcode;

  MUTIL_TRACE( "matuniv_malloc() done" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_malloc in mat_univ.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_malloc( double_mat *matrix, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_MALLOC( matdbl, matrix, nrow, ncol, double );
}


/* This function is documented under matuniv_malloc in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matflt_malloc( float_mat *matrix, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_MALLOC( matflt, matrix, nrow, ncol, float );
}


/* This function is documented under matuniv_malloc in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode mats32_malloc( sint32_mat *matrix, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_MALLOC( mats32, matrix, nrow, ncol, sint32 );
}


/* This function is documented under matuniv_malloc in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats16_malloc( sint16_mat *matrix, sint32 nrow, sint32 ncol )
{
  TMPL_MAT_MALLOC( mats16, matrix, nrow, ncol, sint16 );
}


/* This function is documented under matuniv_malloc in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode matu32_malloc( uint32_mat *matrix, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_MALLOC( matu32, matrix, nrow, ncol, uint32 );
}


/* This function is documented under matuniv_malloc in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu16_malloc( uint16_mat *matrix, sint32 nrow, sint32 ncol )
{
  TMPL_MAT_MALLOC( matu16, matrix, nrow, ncol, uint16 );
}


/* This function is documented under matuniv_malloc in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu8_malloc( uint8_mat *matrix, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_MALLOC( matu8, matrix, nrow, ncol, uint8 );
}


/* This function is documented under matuniv_malloc in mat_univ.h */
/* Written by Luca Cazzanti */
mutil_errcode matcpx_malloc( dcomplex_mat *matrix, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_MALLOC( matcpx, matrix, nrow, ncol, dcomplex );
}


/* Function documented in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode matuniv_realloc( univ_mat *matrix, sint32 nrow, sint32 ncol )
{
  mutil_errcode errcode;
  MUTIL_TRACE( "Start matuniv_realloc()" );

  /* sanity checks... */

  if( !matrix ) {
    MUTIL_ERROR( "Cannot reallocate NULL matrix" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do the reallocation */

  switch( matrix->type ) {
    case MUTIL_DOUBLE:
      errcode = matdbl_realloc( &(matrix->mat.dblmat), nrow, ncol );
      break;

    case MUTIL_FLOAT:
      errcode = matflt_realloc( &(matrix->mat.fltmat), nrow, ncol );
      break;

    case MUTIL_UINT8:
      errcode = matu8_realloc( &(matrix->mat.u8mat), nrow, ncol );
      break;

    case MUTIL_UINT16:
      errcode = matu16_realloc( &(matrix->mat.u16mat), nrow, ncol );
      break;

    case MUTIL_UINT32:
      errcode = matu32_realloc( &(matrix->mat.u32mat), nrow, ncol );
      break;

    case MUTIL_SINT16:
      errcode = mats16_realloc( &(matrix->mat.s16mat), nrow, ncol );
      break;

    case MUTIL_SINT32:
      errcode = mats32_realloc( &(matrix->mat.s32mat), nrow, ncol );
      break;

    case MUTIL_DCOMPLEX:
      errcode = matcpx_realloc( &(matrix->mat.cpxmat), nrow, ncol );
      break;

    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( errcode ) return errcode;

  MUTIL_TRACE( "matuniv_realloc() done" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_realloc in mat_univ.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_realloc( double_mat *matrix, sint32 nrow,
  sint32 ncol)
{
  TMPL_MAT_REALLOC( matdbl, matrix, nrow, ncol, double );
}


/* This function is documented under matuniv_realloc in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matflt_realloc( float_mat *matrix, sint32 nrow, sint32 ncol )
{
  TMPL_MAT_REALLOC( matflt, matrix, nrow, ncol, float );
}


/* This function is documented under matuniv_realloc in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode mats32_realloc( sint32_mat *matrix, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_REALLOC( mats32, matrix, nrow, ncol, sint32 );
}


/* This function is documented under matuniv_realloc in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats16_realloc( sint16_mat *matrix, sint32 nrow, sint32 ncol )
{
  TMPL_MAT_REALLOC( mats16, matrix, nrow, ncol, sint16 );
}


/* This function is documented under matuniv_realloc in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode matu32_realloc( uint32_mat *matrix, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_REALLOC( matu32, matrix, nrow, ncol, uint32 );
}


/* This function is documented under matuniv_realloc in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu16_realloc( uint16_mat *matrix, sint32 nrow, sint32 ncol )
{
  TMPL_MAT_REALLOC( matu16, matrix, nrow, ncol, uint16 );
}


/* This function is documented under matuniv_realloc in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu8_realloc( uint8_mat *matrix, sint32 nrow, sint32 ncol )
{
  TMPL_MAT_REALLOC( matu8, matrix, nrow, ncol, uint8 );
}


/* This function is documented under matuniv_realloc in mat_univ.h */
/* Written by Luca Cazzanti */
mutil_errcode matcpx_realloc( dcomplex_mat *matrix, sint32 nrow,
  sint32 ncol )
{
  TMPL_MAT_REALLOC( matcpx, matrix, nrow, ncol, dcomplex );
}


/* Function documented in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode matuniv_free( univ_mat *mat )
{
  mutil_errcode errcode;

  MUTIL_TRACE( "Start matuniv_free()" );

  if( !mat ) {
    MUTIL_WARN( "Attempt to free NULL address" );
    return MUTIL_ERR_OK;
  }

  switch( mat->type ) {
    case MUTIL_DOUBLE:
      errcode = matdbl_free( &(mat->mat.dblmat) );
      break;

    case MUTIL_FLOAT:
      errcode = matflt_free( &(mat->mat.fltmat) );
      break;

    case MUTIL_UINT8:
      errcode = matu8_free( &(mat->mat.u8mat) );
      break;

    case MUTIL_UINT16:
      errcode = matu16_free( &(mat->mat.u16mat) );
      break;

    case MUTIL_UINT32:
      errcode = matu32_free( &(mat->mat.u32mat) );
      break;

    case MUTIL_SINT16:
      errcode = mats16_free( &(mat->mat.s16mat) );
      break;

    case MUTIL_SINT32:
      errcode = mats32_free( &(mat->mat.s32mat) );
      break;

    case MUTIL_DCOMPLEX:
      errcode = matcpx_free( &(mat->mat.cpxmat) );
      break;

    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( errcode ) return errcode;

  MUTIL_TRACE( "matuniv_free() done" );
  return MUTIL_ERR_OK;
}


/* Function documented under matuniv_free in mat_univ.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_free( double_mat *mat )
{
  TMPL_MAT_FREE( matdbl, mat, sizeof( double ));
}


/* Function documented under matuniv_free in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matflt_free( float_mat *mat )
{
  TMPL_MAT_FREE( matflt, mat, sizeof( float ));
}


/* Function documented under matuniv_free in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode mats32_free( sint32_mat *mat )
{
  TMPL_MAT_FREE( mats32, mat, sizeof( sint32 ));
}


/* Function documented under matuniv_free in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats16_free( sint16_mat *mat )
{
  TMPL_MAT_FREE( mats16, mat, sizeof( sint16 ));
}


/* Function documented under matuniv_free in mat_univ.h */
/* Written by Qin Cai */
mutil_errcode matu32_free( uint32_mat *mat )
{
  TMPL_MAT_FREE( matu32, mat, sizeof( uint32 ));
}


/* Function documented under matuniv_free in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu16_free( uint16_mat *mat )
{
  TMPL_MAT_FREE( matu16, mat, sizeof( uint16 ));
}


/* Function documented under matuniv_free in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu8_free( uint8_mat *mat )
{
  TMPL_MAT_FREE( matu8, mat, sizeof( uint8 ));
}


/* Function documented under matuniv_free in mat_univ.h */
/* Written by Luca Cazzanti */
mutil_errcode matcpx_free( dcomplex_mat *matrix)
{
  TMPL_MAT_FREE( matcpx, matrix, sizeof( dcomplex ));
}


/* Function documented in mat_univ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matuniv_validate( const univ_mat *matrix )
{
  mutil_errcode errcode;

  MUTIL_TRACE( "Start matuniv_validate()" );

  if( !matrix ){
    MUTIL_ERROR( "Matrix pointer is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  switch( matrix->type ) {
    case MUTIL_DOUBLE:
      errcode = matdbl_validate(&(matrix->mat.dblmat));
      break;

    case MUTIL_FLOAT:
      errcode = matflt_validate(&(matrix->mat.fltmat));
      break;

    case MUTIL_UINT8:
      errcode = matu8_validate(&(matrix->mat.u8mat));
      break;

    case MUTIL_UINT16:
      errcode = matu16_validate(&(matrix->mat.u16mat));
      break;

    case MUTIL_UINT32:
      errcode = matu32_validate(&(matrix->mat.u32mat));
      break;

    case MUTIL_SINT16:
      errcode = mats16_validate(&(matrix->mat.s16mat));
      break;

    case MUTIL_SINT32:
      errcode = mats32_validate(&(matrix->mat.s32mat));
      break;

    case MUTIL_DCOMPLEX:
      errcode = matcpx_validate(&(matrix->mat.cpxmat));
      break;

    /* only some matrices available now */
    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( errcode ) return errcode;

  MUTIL_TRACE( "matuniv_validate() done" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_validate in mat_univ.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_validate( const double_mat *mat )
{
  TMPL_MAT_VALIDATE( matdbl, mat );
}


/* This function is documented under matuniv_validate in mat_univ.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matflt_validate( const float_mat *mat )
{
  TMPL_MAT_VALIDATE( matflt, mat );
}


/* This function is documented under matuniv_validate in mat_univ.h */
/* Written by Qin Cai, uses template */
mutil_errcode mats32_validate( const sint32_mat *mat )
{
  TMPL_MAT_VALIDATE( mats32, mat );
}


/* This function is documented under matuniv_validate in mat_univ.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode mats16_validate( const sint16_mat *mat )
{
  TMPL_MAT_VALIDATE( mats16, mat );
}


/* This function is documented under matuniv_validate in mat_univ.h */
/* Written by Qin Cai, uses template */
mutil_errcode matu32_validate( const uint32_mat *mat )
{
  TMPL_MAT_VALIDATE( matu32, mat );
}


/* This function is documented under matuniv_validate in mat_univ.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu16_validate( const uint16_mat *mat )
{
  TMPL_MAT_VALIDATE( matu16, mat );
}


/* This function is documented under matuniv_validate in mat_univ.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu8_validate( const uint8_mat *mat )
{
  TMPL_MAT_VALIDATE( matu8, mat );
}


/* This function is documented under matuniv_validate in mat_univ.h */
/* Written by Luca Cazzanti */
mutil_errcode matcpx_validate( const dcomplex_mat *mat )
{
  TMPL_MAT_VALIDATE( matcpx, mat );
}


/* This function is documented in mat_univ.h */
/* Written by Andrea Borning. */
mutil_errcode matuniv_verify_aresame(
  const univ_mat *mat1, const univ_mat *mat2 )
{
  mutil_errcode errcode;

  MUTIL_TRACE( "Start matuniv_verify_aresame()" );

  errcode = matuniv_validate( mat1 );
  if ( errcode ) return errcode;

  errcode = matuniv_validate( mat2 );
  if ( errcode ) return errcode;

  if ( mat1->type != mat2->type ) {
    MUTIL_ERROR( "Matrices have different types" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( MATUNIV_NROW( mat1 ) != MATUNIV_NROW( mat2 ) ||
       MATUNIV_NCOL( mat1 ) != MATUNIV_NCOL( mat2 ) ) {
    MUTIL_ERROR( "Matrices have different dimensions" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  MUTIL_TRACE( "matuniv_verify_aresame() done" );
  return MUTIL_ERR_OK;
}
