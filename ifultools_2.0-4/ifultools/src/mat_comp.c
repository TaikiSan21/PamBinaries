
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_comp.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "mat_tmpl.h"
#include "mat_comp.h"

#include "mat_cast.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "mat_usca.h"

#include "ut_debug.h"
#include "ut_intrn.h"

/*
   This file contains the definitions for comparison functions
   universal matrices declared in mat_comp.h such as
   number_equal, element by element min and max, ==, <, >.
*/

/* define local macros */

/** Template macro for calling the appropriate matrix comparison to scalar
 * functions, wrapping the result into a universal matrix (for the matched
 * values matrix only), and freeing the allocated memory upon encountering
 * an error.
 * Macro expanded in the body of the matuniv\_compare\_scalar() function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type,
 *                    such as dbl.
 * @param SCALAR_TYPE Relation scalar type, such as double.
 * @param SCALAR      Relation scalar value.
 * @param INTRP_PTR   Pointer for interrupt handling (function argument).
 * @param MAT_INDEX_MATCH_PTR Pointer to a sint32 matrix containing the
 *                    indices for those values of the input matrix which
 *                    satisfy the comparison relation. If the pointer is
 *                    not NULL,the memory for this matrix is allocated here,
 *                    otherwise no memory is allocated and no values are
 *                    returned.
 * @param MAT_MATCH_PTR Pointer to a matrix of the same type as the input,
 *                    containing the values of the input
 *                    matrix which satisfy the comparison relation. If the
 *                    pointer is not NULL, the memory for this matrix is
 *                    allocated herein, otherwise no memory is allocated
 *                    and no values are returned.
 * @param MUTIL_ERR   MUTIL error code variable.
 * @usage In the body of the matuniv\_compare\_scalar function:
 *     #LOCALDEF_MAT_COMPARE_SCALAR_AND_WRAP( dbl, double, scalar, intrp_ptr, match_index, match_value, err );#
 * @private
 */
#define LOCALDEF_MAT_COMPARE_SCALAR_AND_WRAP( MAT_FN_PREFIX, SCALAR_TYPE, \
  SCALAR, INTRP_PTR, MAT_INDEX_MATCH_PTR, MAT_MATCH_PTR, MUTIL_ERR ) \
  \
  MUTIL_ERR = mat ## MAT_FN_PREFIX ## _compare_scalar( \
    &(mat->mat. MAT_FN_PREFIX ## mat), \
    relation, SCAUNIV_CAST( SCALAR, SCALAR_TYPE ), INTRP_PTR, \
    MAT_INDEX_MATCH_PTR, &( MAT_MATCH_PTR ->mat. MAT_FN_PREFIX ## mat) ); \
  if ( MUTIL_ERR ) return( MUTIL_ERR ); \
  \
  if ( MAT_MATCH_PTR ) { \
    MUTIL_ERR = matuniv_wrap_matrix( MAT_MATCH_PTR, \
      &( MAT_MATCH_PTR ->mat. MAT_FN_PREFIX ## mat ), mat->type ); \
     if ( MUTIL_ERR ) { \
       MUTIL_FREE_WARN( mats32, MAT_INDEX_MATCH_PTR ); \
       MUTIL_FREE_WARN( matuniv, MAT_MATCH_PTR ); \
       return( MUTIL_ERR ); \
     } \
   }


/* Function documented in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matuniv_min( const univ_mat *mat1, const univ_mat *mat2,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat1;
  univ_mat      cast_mat2;
  boolean       alloc_cast1;
  boolean       alloc_cast2;

  MUTIL_TRACE( "Start matuniv_min()" );

  /* avoid lint warning */
  (void) whatssi;

  if( !mat1 || !mat2 || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do standard casting */

  MUTIL_DO_STANDARD_CASTING( mat1, result->type, cast_mat1,
    intrp_ptr, alloc_cast1, errcode );
  if( errcode ) {
    return errcode;
  }

  MUTIL_DO_STANDARD_CASTING( mat2, result->type, cast_mat2,
    intrp_ptr, alloc_cast2, errcode );
  if( errcode ) {
    MUTIL_FREE_STANDARD_CASTING( cast_mat1, alloc_cast1 );
    return errcode;
  }

  /* call the appropriate comparison function by
     switching on the type of the output matrix */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_min( &(cast_mat1.mat.dblmat),
        &(cast_mat2.mat.dblmat), intrp_ptr, &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_min( &(cast_mat1.mat.fltmat),
        &(cast_mat2.mat.fltmat), intrp_ptr, &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_min( &(cast_mat1.mat.u8mat),
        &(cast_mat2.mat.u8mat), intrp_ptr, &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_min( &(cast_mat1.mat.u16mat),
        &(cast_mat2.mat.u16mat), intrp_ptr, &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_min( &(cast_mat1.mat.u32mat),
        &(cast_mat2.mat.u32mat), intrp_ptr, &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_min( &(cast_mat1.mat.s16mat),
        &(cast_mat2.mat.s16mat), intrp_ptr, &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_min( &(cast_mat1.mat.s32mat),
        &(cast_mat2.mat.s32mat), intrp_ptr, &(result->mat.s32mat) );
       break;

    /* all other data types unsupported */
    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrices */
  MUTIL_FREE_STANDARD_CASTING( cast_mat1, alloc_cast1 );
  MUTIL_FREE_STANDARD_CASTING( cast_mat2, alloc_cast2 );

  if( errcode ) {
    return errcode;
  }

  MUTIL_TRACE( "matuniv_min() done" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_min in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matdbl_min( const double_mat *mat1, const double_mat *mat2,
  void *intrp_ptr, double_mat *result )
{
  TMPL_OPERATION_MAT_MAT( matdbl, _min, MUTIL_MIN, mat1, mat2,
    intrp_ptr, result, double );
}


/* This function is documented under matuniv_min in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matflt_min( const float_mat *mat1, const float_mat *mat2,
  void *intrp_ptr, float_mat *result )
{
  TMPL_OPERATION_MAT_MAT( matflt, _min, MUTIL_MIN, mat1, mat2,
    intrp_ptr, result, float );
}


/* This function is documented under matuniv_min in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode mats32_min( const sint32_mat *mat1, const sint32_mat *mat2,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_OPERATION_MAT_MAT( mats32, _min, MUTIL_MIN, mat1, mat2,
    intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_min in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode mats16_min( const sint16_mat *mat1, const sint16_mat *mat2,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_OPERATION_MAT_MAT( mats16, _min, MUTIL_MIN, mat1, mat2,
    intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_min in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu32_min( const uint32_mat *mat1, const uint32_mat *mat2,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_OPERATION_MAT_MAT( matu32, _min, MUTIL_MIN, mat1, mat2,
    intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_min in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu16_min( const uint16_mat *mat1, const uint16_mat *mat2,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_OPERATION_MAT_MAT( matu16, _min, MUTIL_MIN, mat1, mat2,
    intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_min in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu8_min( const uint8_mat *mat1, const uint8_mat *mat2,
  void *intrp_ptr, uint8_mat *result )
{
  TMPL_OPERATION_MAT_MAT( matu8, _min, MUTIL_MIN, mat1, mat2,
    intrp_ptr, result, uint8 );
}


/* Function documented in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matuniv_max( const univ_mat *mat1, const univ_mat *mat2,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat1;
  univ_mat      cast_mat2;
  boolean       alloc_cast1;
  boolean       alloc_cast2;

  MUTIL_TRACE( "Start matuniv_max()" );

  if( !mat1 || !mat2 || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  MUTIL_DO_STANDARD_CASTING( mat1, result->type, cast_mat1,
    intrp_ptr, alloc_cast1, errcode );
  if( errcode ) {
    return errcode;
  }

  MUTIL_DO_STANDARD_CASTING( mat2, result->type, cast_mat2,
    intrp_ptr, alloc_cast2, errcode );
  if( errcode ) {
    MUTIL_FREE_STANDARD_CASTING( cast_mat1, alloc_cast1 );
    return errcode;
  }


  /* call the appropriate comparison function by
     switching on the type of the output matrix */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_max( &(cast_mat1.mat.dblmat),
        &(cast_mat2.mat.dblmat), intrp_ptr, &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_max( &(cast_mat1.mat.fltmat),
        &(cast_mat2.mat.fltmat), intrp_ptr, &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_max( &(cast_mat1.mat.u8mat),
        &(cast_mat2.mat.u8mat), intrp_ptr, &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_max( &(cast_mat1.mat.u16mat),
        &(cast_mat2.mat.u16mat), intrp_ptr, &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_max( &(cast_mat1.mat.u32mat),
        &(cast_mat2.mat.u32mat), intrp_ptr, &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_max( &(cast_mat1.mat.s16mat),
        &(cast_mat2.mat.s16mat), intrp_ptr, &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_max( &(cast_mat1.mat.s32mat),
        &(cast_mat2.mat.s32mat), intrp_ptr, &(result->mat.s32mat) );
       break;

    /* all other data types unsupported */
    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrices */
  MUTIL_FREE_STANDARD_CASTING( cast_mat1, alloc_cast1 );
  MUTIL_FREE_STANDARD_CASTING( cast_mat2, alloc_cast2 );

  if( errcode ) {
    return errcode;
  }

  MUTIL_TRACE( "matuniv_max() done" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matdbl_max( const double_mat *mat1, const double_mat *mat2,
  void *intrp_ptr, double_mat *result )
{
  TMPL_OPERATION_MAT_MAT( matdbl, _max, MUTIL_MAX, mat1, mat2,
    intrp_ptr, result, double );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matflt_max( const float_mat *mat1, const float_mat *mat2,
  void *intrp_ptr, float_mat *result )
{
  TMPL_OPERATION_MAT_MAT( matflt, _max, MUTIL_MAX, mat1, mat2,
    intrp_ptr, result, float );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode mats32_max( const sint32_mat *mat1, const sint32_mat *mat2,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_OPERATION_MAT_MAT( mats32, _max, MUTIL_MAX, mat1, mat2,
    intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode mats16_max( const sint16_mat *mat1, const sint16_mat *mat2,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_OPERATION_MAT_MAT( mats16, _max, MUTIL_MAX, mat1, mat2,
    intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu32_max( const uint32_mat *mat1, const uint32_mat *mat2,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_OPERATION_MAT_MAT( matu32, _max, MUTIL_MAX, mat1, mat2,
    intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu16_max( const uint16_mat *mat1, const uint16_mat *mat2,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_OPERATION_MAT_MAT( matu16, _max, MUTIL_MAX, mat1, mat2,
    intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu8_max( const uint8_mat *mat1, const uint8_mat *mat2,
  void *intrp_ptr, uint8_mat *result )
{
  TMPL_OPERATION_MAT_MAT( matu8, _max, MUTIL_MAX, mat1, mat2,
    intrp_ptr, result, uint8 );
}


/* Function documented in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matuniv_number_equal( const univ_mat *mat1, const univ_mat *mat2,
  void *intrp_ptr, sint32 *result )
{
  mutil_errcode   errcode;

  MUTIL_TRACE( "Start matuniv_number_equal()" );

  if( !mat1 || !mat2 || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* casting doesn't work for this case since we should be able
     to compare all numbers */

  /* call the appropriate comparison function */
  if ( mat1->type == mat2->type ) { /* do fast routines */

    switch( mat1->type ) {
       case MUTIL_DOUBLE:
         errcode = matdbl_number_equal( &(mat1->mat.dblmat),
          &(mat2->mat.dblmat), intrp_ptr, result );
         break;

       case MUTIL_FLOAT:
         errcode = matflt_number_equal( &(mat1->mat.fltmat),
          &(mat2->mat.fltmat), intrp_ptr, result );
         break;

       case MUTIL_UINT8:
         errcode = matu8_number_equal( &(mat1->mat.u8mat),
          &(mat2->mat.u8mat), intrp_ptr, result );
         break;

       case MUTIL_UINT16:
         errcode = matu16_number_equal( &(mat1->mat.u16mat),
          &(mat2->mat.u16mat), intrp_ptr, result );
         break;

       case MUTIL_UINT32:
         errcode = matu32_number_equal( &(mat1->mat.u32mat),
          &(mat2->mat.u32mat), intrp_ptr, result );
         break;

       case MUTIL_SINT16:
         errcode = mats16_number_equal( &(mat1->mat.s16mat),
          &(mat2->mat.s16mat), intrp_ptr, result );
         break;

       case MUTIL_SINT32:
         errcode = mats32_number_equal( &(mat1->mat.s32mat),
          &(mat2->mat.s32mat), intrp_ptr, result );
         break;

      default:
        MUTIL_ERROR( "This matrix type is currently unsupported" );
        errcode = MUTIL_ERR_ILLEGAL_TYPE;
    }

    if ( errcode ) return errcode;
  }

  /* can't handle complex values at this time */
  else {
    MUTIL_ERROR( "Mixed matrix types are currently unsupported" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "matuniv_number_equal() done" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matdbl_number_equal( const double_mat *mat1,
  const double_mat *mat2, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_MAT_TALLY( matdbl, _number_equal, mat1, ==, mat2,
    intrp_ptr, result, double );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matflt_number_equal( const float_mat *mat1,
  const float_mat *mat2, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_MAT_TALLY( matflt, _number_equal, mat1, ==, mat2,
    intrp_ptr, result, float );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode mats32_number_equal( const sint32_mat *mat1,
  const sint32_mat *mat2, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_MAT_TALLY( mats32, _number_equal, mat1, ==, mat2,
    intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode mats16_number_equal( const sint16_mat *mat1,
  const sint16_mat *mat2, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_MAT_TALLY( mats16, _number_equal, mat1, ==, mat2,
    intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu32_number_equal( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_MAT_TALLY( matu32, _number_equal, mat1, ==, mat2,
    intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu16_number_equal( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_MAT_TALLY( matu16, _number_equal, mat1, ==, mat2,
    intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu8_number_equal( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_MAT_TALLY( matu8, _number_equal, mat1, ==, mat2,
    intrp_ptr, result, uint8 );
}


/* Function documented in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matuniv_number_equal_scalar( const univ_mat *mat,
  const univ_scalar scalar, void *intrp_ptr, sint32 *result )
{
  mutil_errcode   errcode;
  double          tmp;

  MUTIL_TRACE( "Start matuniv_number_equal_scalar()" );

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* casting doesn't work for this case since we should be able
     to compare all numbers */
  MUTIL_CHECK_STANDARD_CASTING_SCALAR( scalar, mat->type, errcode, tmp );
  if( errcode ) {
    return errcode;
  }

  switch( mat->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_number_equal_scalar( &(mat->mat.dblmat),
        SCAUNIV_CAST( scalar, double ), intrp_ptr, result );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_number_equal_scalar( &(mat->mat.fltmat),
        SCAUNIV_CAST( scalar, float ), intrp_ptr, result );
       break;

     case MUTIL_UINT8:
       errcode = matu8_number_equal_scalar( &(mat->mat.u8mat),
        SCAUNIV_CAST( scalar, uint8 ), intrp_ptr, result );
       break;

     case MUTIL_UINT16:
       errcode = matu16_number_equal_scalar( &(mat->mat.u16mat),
        SCAUNIV_CAST( scalar, uint16 ), intrp_ptr, result );
       break;

     case MUTIL_UINT32:
       errcode = matu32_number_equal_scalar( &(mat->mat.u32mat),
        SCAUNIV_CAST( scalar, uint32 ), intrp_ptr, result );
       break;

     case MUTIL_SINT16:
       errcode = mats16_number_equal_scalar( &(mat->mat.s16mat),
        SCAUNIV_CAST( scalar, sint16 ), intrp_ptr, result );
       break;

     case MUTIL_SINT32:
       errcode = mats32_number_equal_scalar( &(mat->mat.s32mat),
        SCAUNIV_CAST( scalar, sint32 ), intrp_ptr, result );
       break;

    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( errcode ) return errcode;

  MUTIL_TRACE( "matuniv_number_equal_scalar() done" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matdbl_number_equal_scalar( const double_mat *mat,
  const double scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( matdbl, _number_equal_scalar,
    mat, ==, scalar, intrp_ptr, result, double );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matflt_number_equal_scalar( const float_mat *mat,
  const float scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( matflt, _number_equal_scalar,
    mat, ==, scalar, intrp_ptr, result, float );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode mats32_number_equal_scalar( const sint32_mat *mat,
  const sint32 scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( mats32, _number_equal_scalar,
    mat, ==, scalar, intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode mats16_number_equal_scalar( const sint16_mat *mat,
  const sint16 scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( mats16, _number_equal_scalar,
    mat, ==, scalar, intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu32_number_equal_scalar( const uint32_mat *mat,
  const uint32 scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( matu32, _number_equal_scalar,
    mat, ==, scalar, intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu16_number_equal_scalar( const uint16_mat *mat,
  const uint16 scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( matu16, _number_equal_scalar,
    mat, ==, scalar, intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu8_number_equal_scalar( const uint8_mat *mat,
  const uint8 scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( matu8, _number_equal_scalar,
    mat, ==, scalar, intrp_ptr, result, uint8 );
}


/* Function documented in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matuniv_number_less_than_scalar( const univ_mat *mat,
  const univ_scalar scalar, void *intrp_ptr, sint32 *result )
{
  mutil_errcode   errcode;
  double          tmp;

  MUTIL_TRACE( "Start matuniv_number_less_than_scalar()" );

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* casting doesn't work for this case since we should be able
     to compare all numbers */
  MUTIL_CHECK_STANDARD_CASTING_SCALAR( scalar, mat->type, errcode, tmp );
  if( errcode ) {
    return errcode;
  }

  switch( mat->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_number_less_than_scalar( &(mat->mat.dblmat),
        SCAUNIV_CAST( scalar, double ), intrp_ptr, result );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_number_less_than_scalar( &(mat->mat.fltmat),
        SCAUNIV_CAST( scalar, float ), intrp_ptr, result );
       break;

     case MUTIL_UINT8:
       errcode = matu8_number_less_than_scalar( &(mat->mat.u8mat),
        SCAUNIV_CAST( scalar, uint8 ), intrp_ptr, result );
       break;

     case MUTIL_UINT16:
       errcode = matu16_number_less_than_scalar( &(mat->mat.u16mat),
        SCAUNIV_CAST( scalar, uint16 ), intrp_ptr, result );
       break;

     case MUTIL_UINT32:
       errcode = matu32_number_less_than_scalar( &(mat->mat.u32mat),
        SCAUNIV_CAST( scalar, uint32 ), intrp_ptr, result );
       break;

     case MUTIL_SINT16:
       errcode = mats16_number_less_than_scalar( &(mat->mat.s16mat),
        SCAUNIV_CAST( scalar, sint16 ), intrp_ptr, result );
       break;

     case MUTIL_SINT32:
       errcode = mats32_number_less_than_scalar( &(mat->mat.s32mat),
        SCAUNIV_CAST( scalar, sint32 ), intrp_ptr, result );
       break;

    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( errcode ) return errcode;

  MUTIL_TRACE( "matuniv_number_less_than_scalar() done" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matdbl_number_less_than_scalar( const double_mat *mat,
  const double scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( matdbl, _number_less_than_scalar,
    mat, <, scalar, intrp_ptr, result, double );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matflt_number_less_than_scalar( const float_mat *mat,
  const float scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( matflt, _number_less_than_scalar,
    mat, <, scalar, intrp_ptr, result, float );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode mats32_number_less_than_scalar( const sint32_mat *mat,
  const sint32 scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( mats32, _number_less_than_scalar,
    mat, <, scalar, intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode mats16_number_less_than_scalar( const sint16_mat *mat,
  const sint16 scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( mats16, _number_less_than_scalar,
    mat, <, scalar, intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu32_number_less_than_scalar( const uint32_mat *mat,
  const uint32 scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( matu32, _number_less_than_scalar,
    mat, <, scalar, intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu16_number_less_than_scalar( const uint16_mat *mat,
  const uint16 scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( matu16, _number_less_than_scalar,
    mat, <, scalar, intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_max in mat_comp.h */
/* Written by Jill Goldschneider */
mutil_errcode matu8_number_less_than_scalar( const uint8_mat *mat,
  const uint8 scalar, void *intrp_ptr, sint32 *result )
{
  TMPL_MAT_COMPARE_SCALAR_TALLY( matu8, _number_less_than_scalar,
    mat, <, scalar, intrp_ptr, result, uint8 );
}


/* Function documented in mat_comp.h */
/* Written by William Constantine    */
mutil_errcode matuniv_compare_scalar(
  const univ_mat       *mat,
  const mutil_relation  relation,
  const univ_scalar     scalar,
  void                 *intrp_ptr,
  sint32_mat           *match_index,
  univ_mat             *match_value )
{
  mutil_errcode   err;
  double          tmp;

  MUTIL_TRACE( "Start matuniv_compare_scalar()" );

  if ( !mat ) {
    MUTIL_ERROR( "NULL pointer for input matrix" );
    return MUTIL_ERR_NULL_POINTER;
  }
  if ( !match_index && !match_value ) {
    MUTIL_ERROR( "NULL pointer not allowed for both output matrices" );
    return MUTIL_ERR_NULL_POINTER;
  }

  switch( relation ) {
    case MUTIL_RELATION_LESS_THAN:
    case MUTIL_RELATION_LESS_THAN_OR_EQUAL:
    case MUTIL_RELATION_EQUAL:
    case MUTIL_RELATION_NOT_EQUAL:
    case MUTIL_RELATION_GREATER_THAN:
    case MUTIL_RELATION_GREATER_THAN_OR_EQUAL:
      break;
    default:
      MUTIL_ERROR( "Relation operator is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* casting doesn't work for this case since we should be able
     to compare all numbers */
  MUTIL_CHECK_STANDARD_CASTING_SCALAR( scalar, mat->type, err, tmp );
  if ( err ) return err;

  switch( mat->type ) {
    case MUTIL_DOUBLE:

      LOCALDEF_MAT_COMPARE_SCALAR_AND_WRAP( dbl, double, scalar, intrp_ptr,
        match_index, match_value, err );
      break;

    case MUTIL_FLOAT:
      LOCALDEF_MAT_COMPARE_SCALAR_AND_WRAP( flt, float, scalar, intrp_ptr,
        match_index, match_value, err );
      break;

    case MUTIL_UINT8:
      LOCALDEF_MAT_COMPARE_SCALAR_AND_WRAP( u8, uint8, scalar, intrp_ptr,
        match_index, match_value, err );
      break;

    case MUTIL_UINT16:
      LOCALDEF_MAT_COMPARE_SCALAR_AND_WRAP( u16, uint16, scalar, intrp_ptr,
        match_index, match_value, err );
      break;

    case MUTIL_UINT32:
      LOCALDEF_MAT_COMPARE_SCALAR_AND_WRAP( u32, uint32, scalar, intrp_ptr,
        match_index, match_value, err );
      break;

    case MUTIL_SINT16:
      LOCALDEF_MAT_COMPARE_SCALAR_AND_WRAP( s16, sint16, scalar, intrp_ptr,
        match_index, match_value, err );
      break;

    case MUTIL_SINT32:
      LOCALDEF_MAT_COMPARE_SCALAR_AND_WRAP( s32, sint32, scalar, intrp_ptr,
        match_index, match_value, err );
      break;

    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "matuniv_compare_scalar() done" );
  return MUTIL_ERR_OK;
}

/* This function is documented under matuniv_compare_scalar in mat_comp.h */
/* Written by William Constantine */
mutil_errcode matdbl_compare_scalar( const double_mat *mat,
  const mutil_relation relation, const double scalar,
  void *intrp_ptr, sint32_mat *match_index, double_mat *match_value )
{
  TMPL_MAT_COMPARE_SCALAR( matdbl, mat, relation,
    scalar, intrp_ptr, match_index, match_value );
}


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
/* Written by William Constantine */
mutil_errcode matflt_compare_scalar( const float_mat *mat,
  const mutil_relation relation, const float scalar,
  void *intrp_ptr, sint32_mat *match_index, float_mat *match_value )
{
  TMPL_MAT_COMPARE_SCALAR( matflt, mat, relation,
    scalar, intrp_ptr, match_index, match_value );
}


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
/* Written by William Constantine */
mutil_errcode mats32_compare_scalar( const sint32_mat *mat,
  const mutil_relation relation, const sint32 scalar,
  void *intrp_ptr, sint32_mat *match_index, sint32_mat *match_value )
{
  TMPL_MAT_COMPARE_SCALAR( mats32, mat, relation,
    scalar, intrp_ptr, match_index, match_value );
}


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
/* Written by William Constantine */
mutil_errcode mats16_compare_scalar( const sint16_mat *mat,
  const mutil_relation relation, const sint16 scalar,
  void *intrp_ptr, sint32_mat *match_index, sint16_mat *match_value )
{
  TMPL_MAT_COMPARE_SCALAR( mats16, mat, relation,
    scalar, intrp_ptr, match_index, match_value );
}


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
/* Written by William Constantine */
mutil_errcode matu32_compare_scalar( const uint32_mat *mat,
  const mutil_relation relation, const uint32 scalar,
  void *intrp_ptr, sint32_mat *match_index, uint32_mat *match_value )
{
  TMPL_MAT_COMPARE_SCALAR( matu32, mat, relation,
    scalar, intrp_ptr, match_index, match_value );
}


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
/* Written by William Constantine */
mutil_errcode matu16_compare_scalar( const uint16_mat *mat,
  const mutil_relation relation, const uint16 scalar,
  void *intrp_ptr, sint32_mat *match_index, uint16_mat *match_value )
{
  TMPL_MAT_COMPARE_SCALAR( matu16, mat, relation,
    scalar, intrp_ptr, match_index, match_value );
}


/* This function is documented under matuniv_compare_scalar in mat_comp.h */
/* Written by William Constantine */
mutil_errcode matu8_compare_scalar( const uint8_mat *mat,
  const mutil_relation relation, const uint8 scalar,
  void *intrp_ptr, sint32_mat *match_index, uint8_mat *match_value )
{
  TMPL_MAT_COMPARE_SCALAR( matu8, mat, relation,
    scalar, intrp_ptr, match_index, match_value );
}
