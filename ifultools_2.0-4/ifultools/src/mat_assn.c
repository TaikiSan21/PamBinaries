
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_assn.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "mat_assn.h"
#include "mat_tmpl.h"

#include "mat_umat.h"
#include "mat_univ.h"

#include "ut_debug.h"

/*
   This file contains the definitions for the assignment and
   extraction functions for universal matrices declared in mat_assn.h.
*/


/* Function documented in mat_assn.h */
/* Written by Qin Cai */
mutil_errcode matuniv_assign( const univ_mat *mat, void *intrp_ptr,
  univ_mat *result )
{
  mutil_errcode errcode;

  MUTIL_TRACE("Start matuniv_assign()");

  /* avoid lint warning */
  (void) whatssi;

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result");
    return MUTIL_ERR_NULL_POINTER;
  }

  if( !MATUNIV_CHECK_TYPE( mat, result )) {
    MUTIL_ERROR("Data types of operand and result are inconsistent");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch( mat->type ) {
    case MUTIL_DOUBLE:
      errcode = matdbl_assign( &(mat->mat.dblmat), intrp_ptr,
        &(result->mat.dblmat) );
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_assign( &(mat->mat.fltmat), intrp_ptr,
        &(result->mat.fltmat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_assign( &(mat->mat.u8mat), intrp_ptr,
        &(result->mat.u8mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_assign( &(mat->mat.u16mat), intrp_ptr,
        &(result->mat.u16mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_assign( &(mat->mat.u32mat), intrp_ptr,
        &(result->mat.u32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_assign( &(mat->mat.s16mat), intrp_ptr,
        &(result->mat.s16mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_assign( &(mat->mat.s32mat), intrp_ptr,
        &(result->mat.s32mat) );
      if(errcode) return errcode;
      break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "matuniv_assign() done" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_assign in mat_assn.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_assign( const double_mat *mat,
  void *intrp_ptr, double_mat *result )
{
  /*LINTED: Assignment is OK -- checked range. */
  TMPL_MAT_ASSIGN( matdbl, mat, result, double, intrp_ptr );
}


/* This function is documented under matuniv_assign in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matflt_assign( const float_mat *mat,
  void *intrp_ptr, float_mat *result )
{
  /*LINTED: Assignment is OK -- checked range. */
  TMPL_MAT_ASSIGN( matflt, mat, result, float, intrp_ptr );
}


/* This function is documented under matuniv_assign in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode mats32_assign( const sint32_mat *mat,
  void *intrp_ptr, sint32_mat *result )
{
  /*LINTED: Assignment is OK -- checked range. */
  TMPL_MAT_ASSIGN( mats32, mat, result, sint32, intrp_ptr );
}


/* This function is documented under matuniv_assign in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode mats16_assign( const sint16_mat *mat,
  void *intrp_ptr, sint16_mat *result )
{
  /*LINTED: Assignment is OK -- checked range. */
  TMPL_MAT_ASSIGN( mats16, mat, result, sint16, intrp_ptr );
}


/* This function is documented under matuniv_assign in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu32_assign( const uint32_mat *mat,
  void *intrp_ptr, uint32_mat *result )
{
  /*LINTED: Assignment is OK -- checked range. */
  TMPL_MAT_ASSIGN( matu32, mat, result, uint32, intrp_ptr );
}


/* This function is documented under matuniv_assign in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu16_assign( const uint16_mat *mat,
  void *intrp_ptr, uint16_mat *result )
{
  /*LINTED: Assignment is OK -- checked range. */
  TMPL_MAT_ASSIGN( matu16, mat, result, uint16, intrp_ptr );
}


/* This function is documented under matuniv_assign in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu8_assign( const uint8_mat *mat,
  void *intrp_ptr, uint8_mat *result )
{
  /*LINTED: Assignment is OK -- checked range. */
  TMPL_MAT_ASSIGN( matu8, mat, result, uint8, intrp_ptr );
}


/* Function documented in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matuniv_assign_scalar( univ_scalar scalar, void *intrp_ptr,
  univ_mat *mat )
{
  mutil_errcode errcode;

  MUTIL_TRACE("Start matuniv_assign_scalar()");

  if( !mat ) {
    MUTIL_ERROR( "NULL pointer for result");
    return MUTIL_ERR_NULL_POINTER;
  }

  if(!MATUNIV_CHECK_TYPE(mat, &scalar)){
    MUTIL_ERROR("Data types of operand and result are inconsistent");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch(mat->type) {

    case MUTIL_DOUBLE:
      errcode = matdbl_assign_scalar( scalar.num.dbl,
        intrp_ptr, &(mat->mat.dblmat));
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_assign_scalar( scalar.num.flt,
        intrp_ptr, &(mat->mat.fltmat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_assign_scalar( scalar.num.u8,
        intrp_ptr, &(mat->mat.u8mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_assign_scalar( scalar.num.u16,
        intrp_ptr, &(mat->mat.u16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_assign_scalar( scalar.num.u32,
        intrp_ptr, &(mat->mat.u32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_assign_scalar( scalar.num.s16,
        intrp_ptr, &(mat->mat.s16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_assign_scalar( scalar.num.s32,
        intrp_ptr, &(mat->mat.s32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_DCOMPLEX:
      errcode = matcpx_assign_scalar( scalar.num.cpx,
        intrp_ptr, &(mat->mat.cpxmat));
      if(errcode) return errcode;
      break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE("matuniv_assign_scalar() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
/* Originally written by Qin Cai.  Modified to use template by
   Jennifer Hodgdon */
mutil_errcode matdbl_assign_scalar( double scalar, void *intrp_ptr,
  double_mat *mat)
{
  TMPL_MAT_ASSIGN_SCALAR( matdbl, scalar, mat, double, intrp_ptr );
}


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matflt_assign_scalar( float scalar,
  void *intrp_ptr, float_mat *mat )
{
  TMPL_MAT_ASSIGN_SCALAR( matflt, scalar, mat, float, intrp_ptr );
}


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode mats32_assign_scalar( sint32 scalar,
  void *intrp_ptr, sint32_mat *mat )
{
  TMPL_MAT_ASSIGN_SCALAR( mats32, scalar, mat, sint32, intrp_ptr );
}


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode mats16_assign_scalar( sint16 scalar,
  void *intrp_ptr, sint16_mat *mat )
{
  TMPL_MAT_ASSIGN_SCALAR( mats16, scalar, mat, sint16, intrp_ptr );
}


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu32_assign_scalar( uint32 scalar,
  void *intrp_ptr, uint32_mat *mat )
{
  TMPL_MAT_ASSIGN_SCALAR( matu32, scalar, mat, uint32, intrp_ptr );
}


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu16_assign_scalar( uint16 scalar,
  void *intrp_ptr, uint16_mat *mat )
{
  TMPL_MAT_ASSIGN_SCALAR( matu16, scalar, mat, uint16, intrp_ptr );
}


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu8_assign_scalar( uint8 scalar,
  void *intrp_ptr, uint8_mat *mat )
{
  TMPL_MAT_ASSIGN_SCALAR( matu8, scalar, mat, uint8, intrp_ptr );
}


/* This function is documented under matuniv_assign_scalar in mat_assn.h */
/* Written by Jill Goldschneider, uses template. */
mutil_errcode matcpx_assign_scalar( dcomplex scalar,
  void *intrp_ptr, dcomplex_mat *mat )
{
  TMPL_MAT_ASSIGN_SCALAR( matcpx, scalar, mat, dcomplex, intrp_ptr );
}


/* Function documented in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matuniv_assign_zeropad(const univ_mat *smallmat,
  void *intrp_ptr, univ_mat *result)
{
  mutil_errcode errcode;

  MUTIL_TRACE("Start matuniv_assign_zeropad()");

  if( !smallmat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result");
    return MUTIL_ERR_NULL_POINTER;
  }

  if(!MATUNIV_CHECK_TYPE(smallmat, result)){
    MUTIL_ERROR("Data types of operand and result are inconsistent");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch(smallmat->type) {
    case MUTIL_DOUBLE:
      errcode = matdbl_assign_zeropad(&(smallmat->mat.dblmat),
        intrp_ptr, &(result->mat.dblmat));
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_assign_zeropad(&(smallmat->mat.fltmat),
        intrp_ptr, &(result->mat.fltmat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_assign_zeropad(&(smallmat->mat.u8mat),
        intrp_ptr, &(result->mat.u8mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_assign_zeropad(&(smallmat->mat.u16mat),
        intrp_ptr, &(result->mat.u16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_assign_zeropad(&(smallmat->mat.u32mat),
        intrp_ptr, &(result->mat.u32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_assign_zeropad(&(smallmat->mat.s16mat),
        intrp_ptr, &(result->mat.s16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_assign_zeropad(&(smallmat->mat.s32mat),
        intrp_ptr, &(result->mat.s32mat));
      if(errcode) return errcode;
      break;

    /* not all types available yet */
    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE("matuniv_assign_zeropad() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_assign_zeropad( const double_mat *smallmat,
  void *intrp_ptr, double_mat *result )
{
  TMPL_MAT_ASSIGN_PAD( matdbl, smallmat, result, intrp_ptr );
}


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matflt_assign_zeropad( const float_mat *smallmat,
  void *intrp_ptr, float_mat *result )
{
  TMPL_MAT_ASSIGN_PAD( matflt, smallmat, result, intrp_ptr );
}


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode mats32_assign_zeropad( const sint32_mat *smallmat,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_MAT_ASSIGN_PAD( mats32, smallmat, result, intrp_ptr );
}


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode mats16_assign_zeropad( const sint16_mat *smallmat,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_MAT_ASSIGN_PAD( mats16, smallmat, result, intrp_ptr );
}


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu32_assign_zeropad( const uint32_mat *smallmat,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_MAT_ASSIGN_PAD( matu32, smallmat, result, intrp_ptr );
}


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu16_assign_zeropad( const uint16_mat *smallmat,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_MAT_ASSIGN_PAD( matu16, smallmat, result, intrp_ptr );
}


/* This function is documented under matuniv_assign_zeropad in mat_assn.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu8_assign_zeropad( const uint8_mat *smallmat,
  void *intrp_ptr, uint8_mat *result )
{
  TMPL_MAT_ASSIGN_PAD( matu8, smallmat, result, intrp_ptr );
}


/* Function documented in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matuniv_assign_submat(const univ_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, univ_mat *result)
{
  mutil_errcode errcode;

  MUTIL_TRACE("Start matuniv_assign_submat()");

  if( !smallmat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result");
    return MUTIL_ERR_NULL_POINTER;
  }

  if(!MATUNIV_CHECK_TYPE(smallmat, result)) {
    MUTIL_ERROR("Data types of operand and result are inconsistent");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch(smallmat->type) {
    case MUTIL_DOUBLE:
      errcode = matdbl_assign_submat(&(smallmat->mat.dblmat),
        row, col, intrp_ptr, &(result->mat.dblmat));
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_assign_submat(&(smallmat->mat.fltmat),
        row, col, intrp_ptr, &(result->mat.fltmat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_assign_submat(&(smallmat->mat.u8mat),
        row, col, intrp_ptr, &(result->mat.u8mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_assign_submat(&(smallmat->mat.u16mat),
        row, col, intrp_ptr, &(result->mat.u16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_assign_submat(&(smallmat->mat.u32mat),
        row, col, intrp_ptr, &(result->mat.u32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_assign_submat(&(smallmat->mat.s16mat),
        row, col, intrp_ptr, &(result->mat.s16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_assign_submat(&(smallmat->mat.s32mat),
        row, col, intrp_ptr, &(result->mat.s32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_DCOMPLEX:
      errcode = matcpx_assign_submat(&(smallmat->mat.cpxmat),
        row, col, intrp_ptr, &(result->mat.cpxmat));
      if(errcode) return errcode;
      break;

    /* not all types available yet */
    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE("matuniv_assign_submat() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_assign_submat in mat_assn.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_assign_submat( const double_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, double_mat *result )
{
  TMPL_MAT_ASSIGN_SUBMAT( matdbl, smallmat, result, row, col, intrp_ptr );
}


/* This function is documented under matuniv_assign_submat in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matflt_assign_submat( const float_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, float_mat *result )
{
  TMPL_MAT_ASSIGN_SUBMAT( matflt, smallmat, result, row, col, intrp_ptr );
}


/* This function is documented under matuniv_assign_submat in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats32_assign_submat( const sint32_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, sint32_mat *result )
{
  TMPL_MAT_ASSIGN_SUBMAT( mats32, smallmat, result, row, col, intrp_ptr );
}


/* This function is documented under matuniv_assign_submat in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats16_assign_submat( const sint16_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, sint16_mat *result )
{
  TMPL_MAT_ASSIGN_SUBMAT( mats16, smallmat, result, row, col, intrp_ptr );
}


/* This function is documented under matuniv_assign_submat in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu32_assign_submat( const uint32_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, uint32_mat *result )
{
  TMPL_MAT_ASSIGN_SUBMAT( matu32, smallmat, result, row, col, intrp_ptr );
}


/* This function is documented under matuniv_assign_submat in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu16_assign_submat( const uint16_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, uint16_mat *result )
{
  TMPL_MAT_ASSIGN_SUBMAT( matu16, smallmat, result, row, col, intrp_ptr );
}


/* This function is documented under matuniv_assign_submat in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu8_assign_submat( const uint8_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, uint8_mat *result )
{
  TMPL_MAT_ASSIGN_SUBMAT( matu8, smallmat, result, row, col, intrp_ptr );
}


/* This function is documented under matuniv_assign_submat in mat_assn.h */
/* Written by Luca Cazzanti */
mutil_errcode matcpx_assign_submat( const dcomplex_mat *smallmat,
  sint32 row, sint32 col, void *intrp_ptr, dcomplex_mat *result )
{
  TMPL_MAT_ASSIGN_SUBMAT( matcpx, smallmat, result, row, col, intrp_ptr );
}


/* Function documented in mat_assn.h */
/* Written by Qin Cai */
mutil_errcode matuniv_extract( const univ_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;

  MUTIL_TRACE("Start matuniv_extract()");

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result");
    return MUTIL_ERR_NULL_POINTER;
  }

  if(!MATUNIV_CHECK_TYPE(mat, result)){
    MUTIL_ERROR("Data types of operand and result are inconsistent");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch(mat->type){
    case MUTIL_DOUBLE:
      errcode = matdbl_extract(&(mat->mat.dblmat),
        start_row, start_col, intrp_ptr, &(result->mat.dblmat));
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_extract(&(mat->mat.fltmat),
        start_row, start_col, intrp_ptr, &(result->mat.fltmat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_extract(&(mat->mat.u8mat),
        start_row, start_col, intrp_ptr, &(result->mat.u8mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_extract(&(mat->mat.u16mat),
        start_row, start_col, intrp_ptr, &(result->mat.u16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_extract(&(mat->mat.u32mat),
        start_row, start_col, intrp_ptr, &(result->mat.u32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_extract(&(mat->mat.s16mat),
        start_row, start_col, intrp_ptr, &(result->mat.s16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_extract(&(mat->mat.s32mat),
        start_row, start_col, intrp_ptr, &(result->mat.s32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_DCOMPLEX:
      errcode = matcpx_extract(&(mat->mat.cpxmat),
        start_row, start_col, intrp_ptr, &(result->mat.cpxmat));
      if(errcode) return errcode;
      break;

    /* only some matrix types available now */
    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE("matuniv_extract() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_extract in mat_assn.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_extract( const double_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, double_mat *result)
{
  TMPL_MAT_EXTRACT( matdbl, mat, result, start_row, start_col, intrp_ptr );
}


/* This function is documented under matuniv_extract in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matflt_extract( const float_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, float_mat *result )
{
  TMPL_MAT_EXTRACT( matflt, mat, result, start_row, start_col, intrp_ptr );
}


/* This function is documented under matuniv_extract in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats32_extract( const sint32_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, sint32_mat *result )
{
  TMPL_MAT_EXTRACT( mats32, mat, result, start_row, start_col, intrp_ptr );
}


/* This function is documented under matuniv_extract in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats16_extract( const sint16_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, sint16_mat *result )
{
  TMPL_MAT_EXTRACT( mats16, mat, result, start_row, start_col, intrp_ptr );
}


/* This function is documented under matuniv_extract in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu32_extract( const uint32_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, uint32_mat *result )
{
  TMPL_MAT_EXTRACT( matu32, mat, result, start_row, start_col, intrp_ptr );
}


/* This function is documented under matuniv_extract in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu16_extract( const uint16_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, uint16_mat *result )
{
  TMPL_MAT_EXTRACT( matu16, mat, result, start_row, start_col, intrp_ptr );
}


/* This function is documented under matuniv_extract in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu8_extract( const uint8_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, uint8_mat *result )
{
  TMPL_MAT_EXTRACT( matu8, mat, result, start_row, start_col, intrp_ptr );
}

/* This function is documented under matuniv_extract in mat_assn.h */
/* Written by Luca Cazzanti */
mutil_errcode matcpx_extract( const dcomplex_mat *mat,
  sint32 start_row, sint32 start_col, void *intrp_ptr, dcomplex_mat *result )
{
  TMPL_MAT_EXTRACT( matcpx, mat, result, start_row, start_col, intrp_ptr );
}


/* Function documented in mat_assn.h */
/* Written by Qin Cai */
mutil_errcode matuniv_transpose( const univ_mat *mat,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;

  MUTIL_TRACE("Start matuniv_transpose()");

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result");
    return MUTIL_ERR_NULL_POINTER;
  }

  if(!MATUNIV_CHECK_TYPE(mat, result)) {
    MUTIL_ERROR("Data types of operand and result are inconsistent");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch(mat->type) {
    case MUTIL_DOUBLE:
      errcode = matdbl_transpose(&(mat->mat.dblmat),
        intrp_ptr, &(result->mat.dblmat));
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_transpose(&(mat->mat.fltmat),
        intrp_ptr, &(result->mat.fltmat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_transpose(&(mat->mat.u8mat),
        intrp_ptr, &(result->mat.u8mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_transpose(&(mat->mat.u16mat),
        intrp_ptr, &(result->mat.u16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_transpose(&(mat->mat.u32mat),
        intrp_ptr, &(result->mat.u32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_transpose(&(mat->mat.s16mat),
        intrp_ptr, &(result->mat.s16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_transpose(&(mat->mat.s32mat),
        intrp_ptr, &(result->mat.s32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_DCOMPLEX:
      errcode = matcpx_transpose(&(mat->mat.cpxmat),
        intrp_ptr, &(result->mat.cpxmat));
      if(errcode) return errcode;
      break;

    /* other types not available now */
    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE("matuniv_transpose() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_transpose in mat_assn.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_transpose( const double_mat *mat,
  void *intrp_ptr, double_mat *result )
{
  TMPL_MAT_TRANSPOSE( matdbl, mat, result, intrp_ptr );
}


/* This function is documented under matuniv_transpose in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matflt_transpose( const float_mat *mat,
  void *intrp_ptr, float_mat *result )
{
  TMPL_MAT_TRANSPOSE( matflt, mat, result, intrp_ptr );
}


/* This function is documented under matuniv_transpose in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats32_transpose( const sint32_mat *mat,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_MAT_TRANSPOSE( mats32, mat, result, intrp_ptr );
}


/* This function is documented under matuniv_transpose in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats16_transpose( const sint16_mat *mat,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_MAT_TRANSPOSE( mats16, mat, result, intrp_ptr );
}


/* This function is documented under matuniv_transpose in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu32_transpose( const uint32_mat *mat,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_MAT_TRANSPOSE( matu32, mat, result, intrp_ptr );
}


/* This function is documented under matuniv_transpose in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu16_transpose( const uint16_mat *mat,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_MAT_TRANSPOSE( matu16, mat, result, intrp_ptr );
}


/* This function is documented under matuniv_transpose in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu8_transpose( const uint8_mat *mat,
  void *intrp_ptr, uint8_mat *result )
{
  TMPL_MAT_TRANSPOSE( matu8, mat, result, intrp_ptr );
}



/* This function is documented under matuniv_transpose in mat_assn.h */
/* Written by Luca Cazzanti */
mutil_errcode matcpx_transpose( const dcomplex_mat *cmat,
  void *intrp_ptr, dcomplex_mat *result)
{
  TMPL_MAT_TRANSPOSE(matcpx, cmat, result, intrp_ptr );
}


/* Function documented in mat_assn.h */
/* Written by Luca Cazzanti */
mutil_errcode matuniv_translate( const univ_mat *matrix,
  sint32 row_shift, sint32 col_shift, univ_scalar pad_value,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;

  MUTIL_TRACE("Start matuniv_translate()");

  if( !matrix || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result");
    return MUTIL_ERR_NULL_POINTER;
  }

  if(!MATUNIV_CHECK_TYPE(matrix, result)) {
    MUTIL_ERROR("Data types of matrices  are inconsistent");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if(!MATUNIV_CHECK_TYPE(matrix, &pad_value )) {
    MUTIL_ERROR("Data types of matrix and pad value are inconsistent");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch(matrix->type) {
    case MUTIL_DOUBLE:
      errcode = matdbl_translate(&(matrix->mat.dblmat),
        row_shift, col_shift, pad_value.num.dbl,
        intrp_ptr, &(result->mat.dblmat));
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_translate(&(matrix->mat.fltmat),
        row_shift, col_shift, pad_value.num.flt,
        intrp_ptr, &(result->mat.fltmat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_translate(&(matrix->mat.u8mat),
        row_shift, col_shift, pad_value.num.u8,
        intrp_ptr, &(result->mat.u8mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_translate(&(matrix->mat.u16mat),
        row_shift, col_shift, pad_value.num.u16,
        intrp_ptr, &(result->mat.u16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_translate(&(matrix->mat.u32mat),
        row_shift, col_shift, pad_value.num.u32,
        intrp_ptr, &(result->mat.u32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_translate(&(matrix->mat.s16mat),
        row_shift, col_shift, pad_value.num.s16,
        intrp_ptr, &(result->mat.s16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_translate(&(matrix->mat.s32mat),
        row_shift, col_shift, pad_value.num.s32,
        intrp_ptr, &(result->mat.s32mat));
      if(errcode) return errcode;
      break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE("matuniv_translate() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_translate in mat_assn.h */
/* Written by Luca Cazzanti */
mutil_errcode matdbl_translate(const double_mat *mat,
  sint32 row_shift, sint32 col_shift, double pad_value,
  void *intrp_ptr, double_mat *result )
{
  TMPL_MAT_TRANSLATE(matdbl, mat, result, row_shift, col_shift,
    pad_value, intrp_ptr );
}


/* This function is documented under matuniv_translate in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matflt_translate( const float_mat *mat,
  sint32 row_shift, sint32 col_shift, float pad_value,
  void *intrp_ptr, float_mat *result )
{
  TMPL_MAT_TRANSLATE(matflt, mat, result, row_shift, col_shift,
    pad_value, intrp_ptr );
}


/* This function is documented under matuniv_translate in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats32_translate( const sint32_mat *mat,
  sint32 row_shift, sint32 col_shift, sint32 pad_value,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_MAT_TRANSLATE(mats32, mat, result, row_shift, col_shift,
    pad_value, intrp_ptr );
}


/* This function is documented under matuniv_translate in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats16_translate( const sint16_mat *mat,
  sint32 row_shift, sint32 col_shift, sint16 pad_value,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_MAT_TRANSLATE(mats16, mat, result, row_shift, col_shift,
    pad_value, intrp_ptr );
}


/* This function is documented under matuniv_translate in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu32_translate( const uint32_mat *mat,
  sint32 row_shift, sint32 col_shift, uint32 pad_value,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_MAT_TRANSLATE(matu32, mat, result, row_shift, col_shift,
    pad_value, intrp_ptr );
}


/* This function is documented under matuniv_translate in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu16_translate( const uint16_mat *mat,
  sint32 row_shift, sint32 col_shift, uint16 pad_value,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_MAT_TRANSLATE(matu16, mat, result, row_shift, col_shift,
    pad_value, intrp_ptr );
}


/* This function is documented under matuniv_translate in mat_assn.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu8_translate( const uint8_mat *mat,
  sint32 row_shift, sint32 col_shift, uint8 pad_value,
  void *intrp_ptr, uint8_mat *result )
{
  TMPL_MAT_TRANSLATE(matu8, mat, result, row_shift, col_shift,
    pad_value, intrp_ptr );
}


/* Function documented in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matuniv_flip_left_right( const univ_mat *mat,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;

  MUTIL_TRACE("Start matuniv_flip_left_right()");

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result");
    return MUTIL_ERR_NULL_POINTER;
  }

  if(!MATUNIV_CHECK_TYPE(mat, result)) {
    MUTIL_ERROR("Data types of operand and result are inconsistent");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch(mat->type) {
    case MUTIL_DOUBLE:
      errcode = matdbl_flip_left_right(&(mat->mat.dblmat),
        intrp_ptr, &(result->mat.dblmat));
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_flip_left_right(&(mat->mat.fltmat),
        intrp_ptr, &(result->mat.fltmat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_flip_left_right(&(mat->mat.u8mat),
        intrp_ptr, &(result->mat.u8mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_flip_left_right(&(mat->mat.u16mat),
        intrp_ptr, &(result->mat.u16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_flip_left_right(&(mat->mat.u32mat),
        intrp_ptr, &(result->mat.u32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_flip_left_right(&(mat->mat.s16mat),
        intrp_ptr, &(result->mat.s16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_flip_left_right(&(mat->mat.s32mat),
        intrp_ptr, &(result->mat.s32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_DCOMPLEX:
      errcode = matcpx_flip_left_right(&(mat->mat.cpxmat),
        intrp_ptr, &(result->mat.cpxmat));
      if(errcode) return errcode;
      break;

    /* other types not available now */
    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE("matuniv_flip_left_right() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Jill Goldschneider. */
mutil_errcode matdbl_flip_left_right( const double_mat *mat,
  void *intrp_ptr, double_mat *result )
{
  TMPL_MAT_FLIP_LEFT_RIGHT( matdbl, mat, intrp_ptr, result, double );
}


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matflt_flip_left_right( const float_mat *mat,
  void *intrp_ptr, float_mat *result )
{
  TMPL_MAT_FLIP_LEFT_RIGHT( matflt, mat, intrp_ptr, result, float );
}


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode mats32_flip_left_right( const sint32_mat *mat,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_MAT_FLIP_LEFT_RIGHT( mats32, mat, intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode mats16_flip_left_right( const sint16_mat *mat,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_MAT_FLIP_LEFT_RIGHT( mats16, mat, intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matu32_flip_left_right( const uint32_mat *mat,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_MAT_FLIP_LEFT_RIGHT( matu32, mat, intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matu16_flip_left_right( const uint16_mat *mat,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_MAT_FLIP_LEFT_RIGHT( matu16, mat, intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_flip_left_right in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matu8_flip_left_right( const uint8_mat *mat,
  void *intrp_ptr, uint8_mat *result )
{
  TMPL_MAT_FLIP_LEFT_RIGHT( matu8, mat, intrp_ptr, result, uint8 );
}



/* This function is documented under matuniv_flip_left_right in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matcpx_flip_left_right( const dcomplex_mat *cmat,
  void *intrp_ptr, dcomplex_mat *result)
{
  TMPL_MAT_FLIP_LEFT_RIGHT(matcpx, cmat, intrp_ptr, result, dcomplex );
}


/* Function documented in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matuniv_flip_up_down( const univ_mat *mat,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;

  MUTIL_TRACE("Start matuniv_flip_up_down()");

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result");
    return MUTIL_ERR_NULL_POINTER;
  }

  if(!MATUNIV_CHECK_TYPE(mat, result)) {
    MUTIL_ERROR("Data types of operand and result are inconsistent");
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch(mat->type) {
    case MUTIL_DOUBLE:
      errcode = matdbl_flip_up_down(&(mat->mat.dblmat),
        intrp_ptr, &(result->mat.dblmat));
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_flip_up_down(&(mat->mat.fltmat),
        intrp_ptr, &(result->mat.fltmat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_flip_up_down(&(mat->mat.u8mat),
        intrp_ptr, &(result->mat.u8mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_flip_up_down(&(mat->mat.u16mat),
        intrp_ptr, &(result->mat.u16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_flip_up_down(&(mat->mat.u32mat),
        intrp_ptr, &(result->mat.u32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_flip_up_down(&(mat->mat.s16mat),
        intrp_ptr, &(result->mat.s16mat));
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_flip_up_down(&(mat->mat.s32mat),
        intrp_ptr, &(result->mat.s32mat));
      if(errcode) return errcode;
      break;

    case MUTIL_DCOMPLEX:
      errcode = matcpx_flip_up_down(&(mat->mat.cpxmat),
        intrp_ptr, &(result->mat.cpxmat));
      if(errcode) return errcode;
      break;

    /* other types not available now */
    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE("matuniv_flip_up_down() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Jill Goldschneider. */
mutil_errcode matdbl_flip_up_down( const double_mat *mat,
  void *intrp_ptr, double_mat *result )
{
  TMPL_MAT_FLIP_UP_DOWN( matdbl, mat, intrp_ptr, result, double );
}


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matflt_flip_up_down( const float_mat *mat,
  void *intrp_ptr, float_mat *result )
{
  TMPL_MAT_FLIP_UP_DOWN( matflt, mat, intrp_ptr, result, float );
}


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode mats32_flip_up_down( const sint32_mat *mat,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_MAT_FLIP_UP_DOWN( mats32, mat, intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode mats16_flip_up_down( const sint16_mat *mat,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_MAT_FLIP_UP_DOWN( mats16, mat, intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matu32_flip_up_down( const uint32_mat *mat,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_MAT_FLIP_UP_DOWN( matu32, mat, intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matu16_flip_up_down( const uint16_mat *mat,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_MAT_FLIP_UP_DOWN( matu16, mat, intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_flip_up_down in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matu8_flip_up_down( const uint8_mat *mat,
  void *intrp_ptr, uint8_mat *result )
{
  TMPL_MAT_FLIP_UP_DOWN( matu8, mat, intrp_ptr, result, uint8 );
}



/* This function is documented under matuniv_flip_up_down in mat_assn.h */
/* Written by Jill Goldschneider */
mutil_errcode matcpx_flip_up_down( const dcomplex_mat *cmat,
  void *intrp_ptr, dcomplex_mat *result)
{
  TMPL_MAT_FLIP_UP_DOWN(matcpx, cmat, intrp_ptr, result, dcomplex );
}


