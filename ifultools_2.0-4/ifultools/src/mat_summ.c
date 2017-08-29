
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_summ.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "mat_summ.h"
#include "mat_tmpl.h"

#include "mat_assn.h"
#include "mat_cast.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "ut_debug.h"
#include "ut_intrn.h"

/*
   This file contains the definitions for the functions for
   universal matrices declared in mat_summ.h.
   This file contains matrix summary functions such as sum and range.
*/


/* Function documented in mat_summ.h */
/* Written by Qin Cai */
/* Edited to generalize for other data types by Vikram Chalana */
mutil_errcode matuniv_sum( const univ_mat *mat, void *intrp_ptr,
  univ_scalar *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat;
  boolean       alloc_cast;

  MUTIL_TRACE("Start matuniv_sum()");

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do standard casting */

  MUTIL_DO_STANDARD_CASTING( mat, result->type, cast_mat,
    intrp_ptr, alloc_cast, errcode );
  if( errcode ) {
    return errcode;
  }


  /* call the appropriate summary function by
     switching on the type of the output data */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_sum( &(cast_mat.mat.dblmat), intrp_ptr,
         &result->num.dbl);
       break;

     case MUTIL_FLOAT:
       errcode = matflt_sum( &(cast_mat.mat.fltmat), intrp_ptr,
         &result->num.flt);
       break;

     case MUTIL_UINT8:
       errcode = matu8_sum( &(cast_mat.mat.u8mat), intrp_ptr,
         &result->num.u8);
       break;

     case MUTIL_UINT16:
       errcode = matu16_sum( &(cast_mat.mat.u16mat), intrp_ptr,
         &result->num.u16);
       break;

     case MUTIL_UINT32:
       errcode = matu32_sum( &(cast_mat.mat.u32mat), intrp_ptr,
         &result->num.u32);
       break;

     case MUTIL_SINT16:
       errcode = mats16_sum( &(cast_mat.mat.s16mat), intrp_ptr,
         &result->num.s16);
       break;

     case MUTIL_SINT32:
       errcode = mats32_sum( &(cast_mat.mat.s32mat), intrp_ptr,
         &result->num.s32);
       break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrix */
  MUTIL_FREE_STANDARD_CASTING( cast_mat, alloc_cast );

  if( errcode ) {
    return errcode;
  }

  MUTIL_TRACE("matuniv_sum() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_sum in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matdbl_sum( const double_mat *mat, void *intrp_ptr,
  double *result)
{
  TMPL_MAT_SUM( matdbl, mat, intrp_ptr, result, double );
}


/* This function is documented under matuniv_sum in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matflt_sum( const float_mat *mat, void *intrp_ptr,
  float *result)
{
  TMPL_MAT_SUM( matflt, mat, intrp_ptr, result, float );
}


/* This function is documented under matuniv_sum in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode mats32_sum( const sint32_mat *mat, void *intrp_ptr,
  sint32 *result )
{
  TMPL_MAT_SUM( mats32, mat, intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_sum in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode mats16_sum( const sint16_mat *mat, void *intrp_ptr,
  sint16 *result)
{
  TMPL_MAT_SUM( mats16, mat, intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_sum in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu32_sum( const uint32_mat *mat, void *intrp_ptr,
  uint32 *result )
{
  TMPL_MAT_SUM( matu32, mat, intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_sum in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu16_sum( const uint16_mat *mat, void *intrp_ptr,
  uint16 *result)
{
  TMPL_MAT_SUM( matu16, mat, intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_sum in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu8_sum( const uint8_mat *mat, void *intrp_ptr,
  uint8 *result)
{
  TMPL_MAT_SUM( matu8, mat, intrp_ptr, result, uint8 );
}


/* Function documented in mat_summ.h */
/* Written by Jill Goldschneider */
mutil_errcode matuniv_cumulative_sum( const univ_mat *umat_in,
  void *intrp_ptr, univ_mat *umat_out )
{
  mutil_errcode trouble;

  MUTIL_TRACE( "Start matuniv_cumulative_sum()" );

  /* avoid lint warning */
  (void) whatssi;

  if ( !umat_in || !umat_out ) {
    MUTIL_ERROR( "NULL universal matrix pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if ( !MATUNIV_CHECK_TYPE( umat_in, umat_out ) ) {
    MUTIL_ERROR( "Data types of operand and result are inconsistent" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch( umat_in->type ) {

    case MUTIL_DOUBLE:
      trouble = matdbl_cumulative_sum( &(umat_in->mat.dblmat),
        intrp_ptr, &(umat_out->mat.dblmat) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_FLOAT:
      trouble = matflt_cumulative_sum( &(umat_in->mat.fltmat),
        intrp_ptr, &(umat_out->mat.fltmat) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_UINT8:
      trouble = matu8_cumulative_sum( &(umat_in->mat.u8mat),
        intrp_ptr, &(umat_out->mat.u8mat) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_UINT16:
      trouble = matu16_cumulative_sum( &(umat_in->mat.u16mat),
        intrp_ptr, &(umat_out->mat.u16mat) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_UINT32:
      trouble = matu32_cumulative_sum( &(umat_in->mat.u32mat),
        intrp_ptr, &(umat_out->mat.u32mat) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_SINT16:
      trouble = mats16_cumulative_sum( &(umat_in->mat.s16mat),
        intrp_ptr, &(umat_out->mat.s16mat) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_SINT32:
      trouble = mats32_cumulative_sum( &(umat_in->mat.s32mat),
        intrp_ptr, &(umat_out->mat.s32mat) );
      if ( trouble ) return trouble;
      break;

    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "Done with matuniv_cumulative_sum()" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
/* Written by Jill Goldschneider */
mutil_errcode matdbl_cumulative_sum( const double_mat *mat_in,
 void *intrp_ptr, double_mat *mat_out )
{
  TMPL_MAT_CUM_SUM( matdbl, mat_in, intrp_ptr, mat_out );
}


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
/* Written by Jill Goldschneider */
mutil_errcode matflt_cumulative_sum( const float_mat *mat_in,
 void *intrp_ptr, float_mat *mat_out )
{
  TMPL_MAT_CUM_SUM( matflt, mat_in, intrp_ptr, mat_out );
}


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
/* Written by Jill Goldschneider */
mutil_errcode mats32_cumulative_sum( const sint32_mat *mat_in,
 void *intrp_ptr, sint32_mat *mat_out )
{
  TMPL_MAT_CUM_SUM( mats32, mat_in, intrp_ptr, mat_out );
}


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
/* Written by Jill Goldschneider */
mutil_errcode mats16_cumulative_sum( const sint16_mat *mat_in,
 void *intrp_ptr, sint16_mat *mat_out )
{
  TMPL_MAT_CUM_SUM( mats16, mat_in, intrp_ptr, mat_out );
}


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
/* Written by Jill Goldschneider */
mutil_errcode matu32_cumulative_sum( const uint32_mat *mat_in,
 void *intrp_ptr, uint32_mat *mat_out )
{
  TMPL_MAT_CUM_SUM( matu32, mat_in, intrp_ptr, mat_out );
}


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
/* Written by Jill Goldschneider */
mutil_errcode matu16_cumulative_sum( const uint16_mat *mat_in,
 void *intrp_ptr, uint16_mat *mat_out )
{
  TMPL_MAT_CUM_SUM( matu16, mat_in, intrp_ptr, mat_out );
}


/* This function is documented under matuniv_cumulative_sum in mat_summ.h */
/* Written by Jill Goldschneider */
mutil_errcode matu8_cumulative_sum( const uint8_mat *mat_in,
 void *intrp_ptr, uint8_mat *mat_out )
{
  TMPL_MAT_CUM_SUM( matu8, mat_in, intrp_ptr, mat_out );
}


/* Function documented in mat_summ.h */
/* Written by Qin Cai */
/* Edited to generalize for other data types by Vikram Chalana */
mutil_errcode matuniv_sum_rows( const univ_mat *mat, void *intrp_ptr,
  univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat;
  boolean       alloc_cast;

  MUTIL_TRACE("Start matuniv_sum_rows()");

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do standard casting */

  MUTIL_DO_STANDARD_CASTING( mat, result->type, cast_mat,
    intrp_ptr, alloc_cast, errcode );
  if( errcode ) {
    return errcode;
  }


  /* call the appropriate summary function by
     switching on the type of the output data */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_sum_rows( &(cast_mat.mat.dblmat), intrp_ptr,
         &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_sum_rows( &(cast_mat.mat.fltmat), intrp_ptr,
         &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_sum_rows( &(cast_mat.mat.u8mat), intrp_ptr,
         &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_sum_rows( &(cast_mat.mat.u16mat), intrp_ptr,
         &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_sum_rows( &(cast_mat.mat.u32mat), intrp_ptr,
         &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_sum_rows( &(cast_mat.mat.s16mat), intrp_ptr,
         &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_sum_rows( &(cast_mat.mat.s32mat), intrp_ptr,
         &(result->mat.s32mat) );
       break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrix */
  MUTIL_FREE_STANDARD_CASTING( cast_mat, alloc_cast );

  if( errcode ) {
    return errcode;
  }

  MUTIL_TRACE("matuniv_sum_rows() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_sum_rows in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matdbl_sum_rows( const double_mat *mat, void *intrp_ptr,
  double_mat *result)
{
  TMPL_MAT_SUM_ROWS( matdbl, mat, intrp_ptr, result, double );
}


/* This function is documented under matuniv_sum_rows in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matflt_sum_rows( const float_mat *mat, void *intrp_ptr,
  float_mat *result)
{
  TMPL_MAT_SUM_ROWS( matflt, mat, intrp_ptr, result, float );
}


/* This function is documented under matuniv_sum_rows in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode mats32_sum_rows( const sint32_mat *mat, void *intrp_ptr,
  sint32_mat *result )
{
  TMPL_MAT_SUM_ROWS( mats32, mat, intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_sum_rows in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode mats16_sum_rows( const sint16_mat *mat, void *intrp_ptr,
  sint16_mat *result)
{
  TMPL_MAT_SUM_ROWS( mats16, mat, intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_sum_rows in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu32_sum_rows( const uint32_mat *mat, void *intrp_ptr,
  uint32_mat *result )
{
  TMPL_MAT_SUM_ROWS( matu32, mat, intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_sum_rows in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu16_sum_rows( const uint16_mat *mat, void *intrp_ptr,
  uint16_mat *result)
{
  TMPL_MAT_SUM_ROWS( matu16, mat, intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_sum_rows in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu8_sum_rows( const uint8_mat *mat, void *intrp_ptr,
  uint8_mat *result)
{
  TMPL_MAT_SUM_ROWS( matu8, mat, intrp_ptr, result, uint8 );
}


/* Function documented in mat_summ.h */
/* Written by Qin Cai */
/* Edited to generalize for other data types by Vikram Chalana */
mutil_errcode matuniv_sum_cols( const univ_mat *mat, void *intrp_ptr,
  univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat;
  boolean       alloc_cast;

  MUTIL_TRACE("Start matuniv_sum_cols()");

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do standard casting */

  MUTIL_DO_STANDARD_CASTING( mat, result->type, cast_mat,
    intrp_ptr, alloc_cast, errcode );
  if( errcode ) {
    return errcode;
  }


  /* call the appropriate summary function by
     switching on the type of the output data */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_sum_cols( &(cast_mat.mat.dblmat), intrp_ptr,
         &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_sum_cols( &(cast_mat.mat.fltmat), intrp_ptr,
         &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_sum_cols( &(cast_mat.mat.u8mat), intrp_ptr,
         &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_sum_cols( &(cast_mat.mat.u16mat), intrp_ptr,
         &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_sum_cols( &(cast_mat.mat.u32mat), intrp_ptr,
         &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_sum_cols( &(cast_mat.mat.s16mat), intrp_ptr,
         &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_sum_cols( &(cast_mat.mat.s32mat), intrp_ptr,
         &(result->mat.s32mat) );
       break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrix */
  MUTIL_FREE_STANDARD_CASTING( cast_mat, alloc_cast );

  if( errcode ) {
    return errcode;
  }

  MUTIL_TRACE("matuniv_sum_cols() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_sum_cols in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matdbl_sum_cols( const double_mat *mat, void *intrp_ptr,
  double_mat *result)
{
  TMPL_MAT_SUM_COLS( matdbl, mat, intrp_ptr, result, double );
}


/* This function is documented under matuniv_sum_cols in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matflt_sum_cols( const float_mat *mat, void *intrp_ptr,
  float_mat *result)
{
  TMPL_MAT_SUM_COLS( matflt, mat, intrp_ptr, result, float );
}


/* This function is documented under matuniv_sum_cols in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode mats32_sum_cols( const sint32_mat *mat, void *intrp_ptr,
  sint32_mat *result )
{
  TMPL_MAT_SUM_COLS( mats32, mat, intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_sum_cols in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode mats16_sum_cols( const sint16_mat *mat, void *intrp_ptr,
  sint16_mat *result)
{
  TMPL_MAT_SUM_COLS( mats16, mat, intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_sum_cols in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu32_sum_cols( const uint32_mat *mat, void *intrp_ptr,
  uint32_mat *result )
{
  TMPL_MAT_SUM_COLS( matu32, mat, intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_sum_cols in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu16_sum_cols( const uint16_mat *mat, void *intrp_ptr,
  uint16_mat *result)
{
  TMPL_MAT_SUM_COLS( matu16, mat, intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_sum_cols in mat_summ.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu8_sum_cols( const uint8_mat *mat, void *intrp_ptr,
  uint8_mat *result)
{
  TMPL_MAT_SUM_COLS( matu8, mat, intrp_ptr, result, uint8 );
}


/* function documented in mat_summ.h */
/* written by Jennifer Hodgdon */
mutil_errcode matuniv_range( const univ_mat *mat, void *intrp_ptr,
  univ_scalar *minval, univ_scalar *maxval )
{
  mutil_errcode errcode;

  MUTIL_TRACE("Start matuniv_range()");

  if( !mat || !minval || !maxval ) {
    MUTIL_ERROR( "NULL pointer for operand");
    return MUTIL_ERR_NULL_POINTER;
  }

  switch(mat->type) {
    case MUTIL_DOUBLE:
      errcode = matdbl_range(&(mat->mat.dblmat), intrp_ptr,
        &(minval->num.dbl), &(maxval->num.dbl) );
      if(errcode) return errcode;
      minval->type = MUTIL_DOUBLE;
      maxval->type = MUTIL_DOUBLE;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_range(&(mat->mat.fltmat), intrp_ptr,
        &(minval->num.flt), &(maxval->num.flt) );
      if(errcode) return errcode;
      minval->type = MUTIL_FLOAT;
      maxval->type = MUTIL_FLOAT;
      break;

    case MUTIL_UINT8:
      errcode = matu8_range(&(mat->mat.u8mat), intrp_ptr,
        &(minval->num.u8), &(maxval->num.u8) );
      if(errcode) return errcode;
      minval->type = MUTIL_UINT8;
      maxval->type = MUTIL_UINT8;
      break;

    case MUTIL_UINT16:
      errcode = matu16_range(&(mat->mat.u16mat), intrp_ptr,
        &(minval->num.u16), &(maxval->num.u16) );
      if(errcode) return errcode;
      minval->type = MUTIL_UINT16;
      maxval->type = MUTIL_UINT16;
      break;

    case MUTIL_UINT32:
      errcode = matu32_range(&(mat->mat.u32mat), intrp_ptr,
        &(minval->num.u32), &(maxval->num.u32) );
      if(errcode) return errcode;
      minval->type = MUTIL_UINT32;
      maxval->type = MUTIL_UINT32;
      break;

    case MUTIL_SINT16:
      errcode = mats16_range(&(mat->mat.s16mat), intrp_ptr,
        &(minval->num.s16), &(maxval->num.s16) );
      if(errcode) return errcode;
      minval->type = MUTIL_SINT16;
      maxval->type = MUTIL_SINT16;
      break;

    case MUTIL_SINT32:
      errcode = mats32_range(&(mat->mat.s32mat), intrp_ptr,
        &(minval->num.s32), &(maxval->num.s32) );
      if(errcode) return errcode;
      minval->type = MUTIL_SINT32;
      maxval->type = MUTIL_SINT32;
      break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE("matuniv_range() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_range in mat_summ.h */
/* Re-written by Jennifer Hodgdon as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_range( const double_mat *mat,
  void *intrp_ptr, double *minval, double *maxval )
{
  TMPL_MAT_RANGE( matdbl, mat, minval, maxval, double, intrp_ptr );
}


/* This function is documented under matuniv_range in mat_summ.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matflt_range( const float_mat *mat,
  void *intrp_ptr, float *minval, float *maxval )
{
  TMPL_MAT_RANGE( matflt, mat, minval, maxval, float, intrp_ptr );
}


/* This function is documented under matuniv_range in mat_summ.h */
/* Written by Qin Cai, but pretty much copied
   from matdbl_range in mat_summ.c
*/
mutil_errcode mats32_range( const sint32_mat *mat,
  void *intrp_ptr, sint32 *minval, sint32 *maxval )
{
  TMPL_MAT_RANGE( mats32, mat, minval, maxval, sint32, intrp_ptr );
}


/* This function is documented under matuniv_range in mat_summ.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode mats16_range( const sint16_mat *mat,
  void *intrp_ptr, sint16 *minval, sint16 *maxval )
{
  TMPL_MAT_RANGE( mats16, mat, minval, maxval, sint16, intrp_ptr );
}


/* This function is documented under matuniv_range in mat_summ.h */
/* Written by Qin Cai, but pretty much copied
   from matdbl_range in mat_summ.c
*/
mutil_errcode matu32_range( const uint32_mat *mat,
  void *intrp_ptr, uint32 *minval, uint32 *maxval )
{
  TMPL_MAT_RANGE( matu32, mat, minval, maxval, uint32, intrp_ptr );
}


/* This function is documented under matuniv_range in mat_summ.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu16_range( const uint16_mat *mat,
  void *intrp_ptr, uint16 *minval, uint16 *maxval )
{
  TMPL_MAT_RANGE( matu16, mat, minval, maxval, uint16, intrp_ptr );
}


/* This function is documented under matuniv_range in mat_summ.h */
/* Written by Jennifer Hodgdon, uses template. */
mutil_errcode matu8_range( const uint8_mat *mat,
  void *intrp_ptr, uint8 *minval, uint8 *maxval )
{
  TMPL_MAT_RANGE( matu8, mat, minval, maxval, uint8, intrp_ptr );
}


/* function documented in mat_summ.h */
/* written by Jennifer Hodgdon */
mutil_errcode matuniv_range_robust( const univ_mat *mat,
  double lowexc, double highexc, void *intrp_ptr,
  univ_scalar *minval, univ_scalar *maxval )
{
  mutil_errcode errcode;
  float         minflt;
  float         maxflt;
  uint32        minu32;
  uint32        maxu32;
  uint16        minu16;
  uint16        maxu16;
  uint8         minu8;
  uint8         maxu8;
  sint16        mins16;
  sint16        maxs16;
  sint32        mins32;
  sint32        maxs32;

  MUTIL_TRACE("Start matuniv_range_robust()");

  if( !mat ) {
    MUTIL_ERROR( "NULL pointer for operand");
    return MUTIL_ERR_NULL_POINTER;
  }

  switch(mat->type) {
    case MUTIL_DOUBLE:
      errcode = matdbl_range_robust(&(mat->mat.dblmat),
        lowexc, highexc, intrp_ptr, &(minval->num.dbl), &(maxval->num.dbl) );
      if(errcode) return errcode;
      minval->type = MUTIL_DOUBLE;
      maxval->type = MUTIL_DOUBLE;
      break;


    case MUTIL_FLOAT:
      if( lowexc <= -MUTIL_FLOAT_MAX ) {
        minflt = (float) -MUTIL_FLOAT_MAX;
      }
      else if( lowexc > MUTIL_FLOAT_MAX ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else {
        minflt = (float) lowexc;
      }

      if( highexc <= -MUTIL_FLOAT_MAX ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else if( highexc > MUTIL_FLOAT_MAX ) {
        maxflt = (float) MUTIL_FLOAT_MAX;
      }
      else {
        maxflt = (float) highexc;
      }

      errcode = matflt_range_robust(&(mat->mat.fltmat),
        minflt, maxflt, intrp_ptr, &(minval->num.flt), &(maxval->num.flt) );
      if(errcode) return errcode;
      minval->type = MUTIL_FLOAT;
      maxval->type = MUTIL_FLOAT;
      break;


    case MUTIL_UINT8:

      /* convert input values to uint8 */
      /* the lowexc and highexc are exclusion zones, excluded with
         strict less than or greater than */

      if( lowexc <= 0.0 ) {
        minu8 = 0;
      }
      else if( lowexc > MUTIL_UINT8_MAX ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else {
        minu8 = (uint8) ceil( lowexc );
      }

      if( highexc <= 0.0 ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else if( highexc > MUTIL_UINT8_MAX ) {
        maxu8 = MUTIL_UINT8_MAX;
      }
      else {
        maxu8 = (uint8) floor( highexc );
      }

      errcode = matu8_range_robust(&(mat->mat.u8mat),
        minu8, maxu8, intrp_ptr, &(minval->num.u8), &(maxval->num.u8) );
      if(errcode) return errcode;
      minval->type = MUTIL_UINT8;
      maxval->type = MUTIL_UINT8;
      break;

    case MUTIL_UINT16:

      /* convert input values to uint16 */
      /* the lowexc and highexc are exclusion zones, excluded with
         strict less than or greater than */

      if( lowexc <= 0.0 ) {
        minu16 = 0;
      }
      else if( lowexc > MUTIL_UINT16_MAX ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else {
        minu16 = (uint16) ceil( lowexc );
      }

      if( highexc <= 0.0 ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else if( highexc > MUTIL_UINT16_MAX ) {
        maxu16 = MUTIL_UINT16_MAX;
      }
      else {
        maxu16 = (uint16) floor( highexc );
      }

      errcode = matu16_range_robust(&(mat->mat.u16mat),
        minu16, maxu16, intrp_ptr, &(minval->num.u16), &(maxval->num.u16) );
      if(errcode) return errcode;
      minval->type = MUTIL_UINT16;
      maxval->type = MUTIL_UINT16;
      break;

    case MUTIL_UINT32:

      /* convert input values to uint32 */
      /* the lowexc and highexc are exclusion zones, excluded with
         strict less than or greater than */

      if( lowexc <= 0.0 ) {
        minu32 = 0;
      }
      else if( lowexc > MUTIL_UINT32_MAX ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else {
        minu32 = (uint32) ceil( lowexc );
      }

      if( highexc <= 0.0 ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else if( highexc > MUTIL_UINT32_MAX ) {
        maxu32 = MUTIL_UINT32_MAX;
      }
      else {
        maxu32 = (uint32) floor( highexc );
      }

      errcode = matu32_range_robust(&(mat->mat.u32mat),
        minu32, maxu32, intrp_ptr, &(minval->num.u32), &(maxval->num.u32) );
      if(errcode) return errcode;
      minval->type = MUTIL_UINT32;
      maxval->type = MUTIL_UINT32;
      break;

    case MUTIL_SINT16:

      /* convert input values to sint16 */
      /* the lowexc and highexc are exclusion zones, excluded with
         strict less than or greater than */

      if( lowexc <= MUTIL_SINT16_MIN ) {
        mins16 = MUTIL_SINT16_MIN;
      }
      else if( lowexc > MUTIL_SINT16_MAX ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else {
        mins16 = (sint16) ceil( lowexc );
      }

      if( highexc <= MUTIL_SINT16_MIN ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else if( highexc > MUTIL_SINT16_MAX ) {
        maxs16 = MUTIL_SINT16_MAX;
      }
      else {
        maxs16 = (sint16) floor( highexc );
      }

      errcode = mats16_range_robust(&(mat->mat.s16mat),
        mins16, maxs16, intrp_ptr, &(minval->num.s16), &(maxval->num.s16) );
      if(errcode) return errcode;
      minval->type = MUTIL_SINT16;
      maxval->type = MUTIL_SINT16;
      break;

    case MUTIL_SINT32:

      /* convert input values to sint32 */
      /* the lowexc and highexc are exclusion zones, excluded with
         strict less than or greater than */

      if( lowexc <= MUTIL_SINT32_MIN ) {
        mins32 = MUTIL_SINT32_MIN;
      }
      else if( lowexc > MUTIL_SINT32_MAX ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else {
        mins32 = (sint32) ceil( lowexc );
      }

      if( highexc <= MUTIL_SINT32_MIN ) {
        MUTIL_ERROR( "No values within range" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }
      else if( highexc > MUTIL_SINT32_MAX ) {
        maxs32 = MUTIL_SINT32_MAX;
      }
      else {
        maxs32 = (sint32) floor( highexc );
      }

      errcode = mats32_range_robust(&(mat->mat.s32mat),
        mins32, maxs32, intrp_ptr, &(minval->num.s32), &(maxval->num.s32) );
      if(errcode) return errcode;
      minval->type = MUTIL_SINT32;
      maxval->type = MUTIL_SINT32;
      break;

   /* nothing else available now */
    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE("matuniv_range_robust() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_range_robust in mat_summ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matdbl_range_robust( const double_mat *mat,
  double lowexc, double highexc, void *intrp_ptr,
  double *minval, double *maxval )
{
  TMPL_MAT_RANGE_ROBUST( matdbl, mat, lowexc, highexc, minval, maxval,
    double, intrp_ptr );
}


/* This function is documented under matuniv_range_robust in mat_summ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matflt_range_robust( const float_mat *mat,
  float lowexc, float highexc, void *intrp_ptr,
  float *minval, float *maxval )
{
  TMPL_MAT_RANGE_ROBUST( matflt, mat, lowexc, highexc, minval, maxval,
    float, intrp_ptr );
}


/* This function is documented under matuniv_range_robust in mat_summ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats32_range_robust( const sint32_mat *mat,
  sint32 lowexc, sint32 highexc, void *intrp_ptr,
  sint32 *minval, sint32 *maxval )
{
  TMPL_MAT_RANGE_ROBUST( mats32, mat, lowexc, highexc, minval, maxval,
    sint32, intrp_ptr );
}


/* This function is documented under matuniv_range_robust in mat_summ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode mats16_range_robust( const sint16_mat *mat,
  sint16 lowexc, sint16 highexc, void *intrp_ptr,
  sint16 *minval, sint16 *maxval )
{
  TMPL_MAT_RANGE_ROBUST( mats16, mat, lowexc, highexc, minval, maxval,
    sint16, intrp_ptr );
}


/* This function is documented under matuniv_range_robust in mat_summ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu32_range_robust( const uint32_mat *mat,
  uint32 lowexc, uint32 highexc, void *intrp_ptr,
  uint32 *minval, uint32 *maxval )
{
  TMPL_MAT_RANGE_ROBUST( matu32, mat, lowexc, highexc, minval, maxval,
    uint32, intrp_ptr );
}


/* This function is documented under matuniv_range_robust in mat_summ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu16_range_robust( const uint16_mat *mat,
  uint16 lowexc, uint16 highexc, void *intrp_ptr,
  uint16 *minval, uint16 *maxval )
{
  TMPL_MAT_RANGE_ROBUST( matu16, mat, lowexc, highexc, minval, maxval,
    uint16, intrp_ptr );
}


/* This function is documented under matuniv_range_robust in mat_summ.h */
/* Written by Jennifer Hodgdon */
mutil_errcode matu8_range_robust( const uint8_mat *mat,
  uint8 lowexc, uint8 highexc, void *intrp_ptr,
  uint8 *minval, uint8 *maxval )
{
  TMPL_MAT_RANGE_ROBUST( matu8, mat, lowexc, highexc, minval, maxval,
    uint8, intrp_ptr );
}


