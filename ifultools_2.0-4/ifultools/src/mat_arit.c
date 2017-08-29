
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_arit.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "mat_arit.h"
#include "mat_tmpl.h"

#include "mat_assn.h"
#include "mat_cast.h"
#include "mat_summ.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "mat_usca.h"

#include "ut_debug.h"
#include "ut_limit.h"
#include "ut_intrn.h"

/*
   This file contains the definitions for the functions for
   universal matrices declared in mat_arit.h.
   This file contains functions normally associated with
   the arithmetic component of an ALU.
*/


/* Function documented in mat_arit.h */
/* Written by Qin Cai */
/* Edited to generalize for other data types by Vikram Chalana 3/29/1999 */
mutil_errcode matuniv_add( const univ_mat *mat1, const univ_mat *mat2,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat1;
  univ_mat      cast_mat2;
  boolean       alloc_cast1;
  boolean       alloc_cast2;

  MUTIL_TRACE( "Start matuniv_add()" );

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

  /* call the appropriate arithmetic function by
     switching on the type of the output matrix */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_add( &(cast_mat1.mat.dblmat),
        &(cast_mat2.mat.dblmat), intrp_ptr, &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_add( &(cast_mat1.mat.fltmat),
        &(cast_mat2.mat.fltmat), intrp_ptr, &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_add( &(cast_mat1.mat.u8mat),
        &(cast_mat2.mat.u8mat), intrp_ptr, &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_add( &(cast_mat1.mat.u16mat),
        &(cast_mat2.mat.u16mat), intrp_ptr, &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_add( &(cast_mat1.mat.u32mat),
        &(cast_mat2.mat.u32mat), intrp_ptr, &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_add( &(cast_mat1.mat.s16mat),
        &(cast_mat2.mat.s16mat), intrp_ptr, &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_add( &(cast_mat1.mat.s32mat),
        &(cast_mat2.mat.s32mat), intrp_ptr, &(result->mat.s32mat) );
       break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrices */
  MUTIL_FREE_STANDARD_CASTING( cast_mat1, alloc_cast1 );
  MUTIL_FREE_STANDARD_CASTING( cast_mat2, alloc_cast2 );

  if( errcode ) {
    return errcode;
  }

  MUTIL_TRACE("matuniv_add() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_add in mat_arit.h */
/* Re-written by Vikram Chalana as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_add( const double_mat *mat1, const double_mat *mat2,
  void *intrp_ptr, double_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matdbl, _add, mat1, +, mat2,
    intrp_ptr, result, double );
}


/* This function is documented under matuniv_add in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matflt_add( const float_mat *mat1, const float_mat *mat2,
  void *intrp_ptr, float_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matflt, _add, mat1, +, mat2,
    intrp_ptr, result, float );
}


/* This function is documented under matuniv_add in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode mats32_add( const sint32_mat *mat1, const sint32_mat *mat2,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( mats32, _add, mat1, +, mat2,
    intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_add in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode mats16_add( const sint16_mat *mat1, const sint16_mat *mat2,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( mats16, _add, mat1, +, mat2,
    intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_add in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu32_add( const uint32_mat *mat1, const uint32_mat *mat2,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matu32, _add, mat1, +, mat2,
    intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_add in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu16_add( const uint16_mat *mat1, const uint16_mat *mat2,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matu16, _add, mat1, +, mat2,
    intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_add in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu8_add( const uint8_mat *mat1, const uint8_mat *mat2,
  void *intrp_ptr, uint8_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matu8, _add, mat1, +, mat2,
    intrp_ptr, result, uint8 );
}


/* Function documented in mat_arit.h */
/* Written by Qin Cai */
/* Edited to generalize for other data types by Vikram Chalana 3/31/1999 */
mutil_errcode matuniv_add_scalar( const univ_mat *mat, univ_scalar scalar,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat;
  boolean       alloc_cast;
  double        tmp;

  MUTIL_TRACE("Start matuniv_add_scalar()");

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do standard casting */

  MUTIL_CHECK_STANDARD_CASTING_SCALAR( scalar, result->type, errcode, tmp );
  if( errcode ) {
    return errcode;
  }

  MUTIL_DO_STANDARD_CASTING( mat, result->type, cast_mat,
    intrp_ptr, alloc_cast, errcode );
  if( errcode ) {
    return errcode;
  }


  /* call the appropriate arithmetic function by
     switching on the type of the output data */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_add_scalar( &(cast_mat.mat.dblmat),
        SCAUNIV_CAST( scalar, double ), intrp_ptr, &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_add_scalar( &(cast_mat.mat.fltmat),
        SCAUNIV_CAST( scalar, float ), intrp_ptr, &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_add_scalar( &(cast_mat.mat.u8mat),
        SCAUNIV_CAST( scalar, uint8 ), intrp_ptr, &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_add_scalar( &(cast_mat.mat.u16mat),
        SCAUNIV_CAST( scalar, uint16 ), intrp_ptr, &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_add_scalar( &(cast_mat.mat.u32mat),
        SCAUNIV_CAST( scalar, uint32 ), intrp_ptr, &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_add_scalar( &(cast_mat.mat.s16mat),
        SCAUNIV_CAST( scalar, sint16 ), intrp_ptr, &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_add_scalar( &(cast_mat.mat.s32mat),
        SCAUNIV_CAST( scalar, sint32 ), intrp_ptr, &(result->mat.s32mat) );
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

  MUTIL_TRACE("matuniv_add_scalar() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_add_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matdbl_add_scalar( const double_mat *mat,
  double scalar, void *intrp_ptr, double_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( matdbl, _add_scalar, mat, +, scalar,
    intrp_ptr, result, double );
}


/* This function is documented under matuniv_add_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matflt_add_scalar( const float_mat *mat,
  float scalar, void *intrp_ptr, float_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( matflt, _add_scalar, mat, +, scalar,
    intrp_ptr, result, float );
}


/* This function is documented under matuniv_add_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode mats32_add_scalar( const sint32_mat *mat,
  sint32 scalar, void *intrp_ptr, sint32_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( mats32, _add_scalar, mat, +, scalar,
    intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_add_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode mats16_add_scalar( const sint16_mat *mat,
  sint16 scalar, void *intrp_ptr, sint16_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( mats16, _add_scalar, mat, +, scalar,
    intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_add_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu32_add_scalar( const uint32_mat *mat,
  uint32 scalar, void *intrp_ptr, uint32_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( matu32, _add_scalar, mat, +, scalar,
    intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_add_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu16_add_scalar( const uint16_mat *mat,
  uint16 scalar, void *intrp_ptr, uint16_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( matu16, _add_scalar, mat, +, scalar,
    intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_add_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu8_add_scalar( const uint8_mat *mat,
  uint8 scalar, void *intrp_ptr, uint8_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( matu8, _add_scalar, mat, +, scalar,
    intrp_ptr, result, uint8 );
}


/* Function documented in mat_arit.h */
/* Written by Qin Cai */
/* Edited to generalize for other data types by Vikram Chalana 3/29/1999 */
mutil_errcode matuniv_subtract( const univ_mat *mat1, const univ_mat *mat2,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat1;
  univ_mat      cast_mat2;
  boolean       alloc_cast1;
  boolean       alloc_cast2;

  MUTIL_TRACE("Start matuniv_subtract()");

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

  /* call the appropriate arithmetic function by
     switching on the type of the output matrix */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_subtract( &(cast_mat1.mat.dblmat),
        &(cast_mat2.mat.dblmat), intrp_ptr, &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_subtract( &(cast_mat1.mat.fltmat),
        &(cast_mat2.mat.fltmat), intrp_ptr, &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_subtract( &(cast_mat1.mat.u8mat),
        &(cast_mat2.mat.u8mat), intrp_ptr, &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_subtract( &(cast_mat1.mat.u16mat),
        &(cast_mat2.mat.u16mat), intrp_ptr, &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_subtract( &(cast_mat1.mat.u32mat),
        &(cast_mat2.mat.u32mat), intrp_ptr, &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_subtract( &(cast_mat1.mat.s16mat),
        &(cast_mat2.mat.s16mat), intrp_ptr, &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_subtract( &(cast_mat1.mat.s32mat),
        &(cast_mat2.mat.s32mat), intrp_ptr, &(result->mat.s32mat) );
       break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }


  /* finally, free the cast matrices */
  MUTIL_FREE_STANDARD_CASTING( cast_mat1, alloc_cast1 );
  MUTIL_FREE_STANDARD_CASTING( cast_mat2, alloc_cast2 );

  if( errcode ) {
    return errcode;
  }

  MUTIL_TRACE("matuniv_subtract() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_subtract in mat_arit.h */
/* Re-written by Vikram Chalana as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_subtract( const double_mat *mat1, const double_mat *mat2,
  void *intrp_ptr, double_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matdbl, _subtract, mat1, -, mat2,
    intrp_ptr, result, double );
}


/* This function is documented under matuniv_subtract in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matflt_subtract( const float_mat *mat1, const float_mat *mat2,
  void *intrp_ptr, float_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matflt, _subtract, mat1, -, mat2,
    intrp_ptr, result, float );
}


/* This function is documented under matuniv_subtract in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode mats32_subtract( const sint32_mat *mat1, const sint32_mat *mat2,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( mats32, _subtract, mat1, -, mat2,
    intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_subtract in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode mats16_subtract( const sint16_mat *mat1, const sint16_mat *mat2,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( mats16, _subtract, mat1, -, mat2,
    intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_subtract in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu32_subtract( const uint32_mat *mat1, const uint32_mat *mat2,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matu32, _subtract, mat1, -, mat2,
    intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_subtract in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu16_subtract( const uint16_mat *mat1, const uint16_mat *mat2,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matu16, _subtract, mat1, -, mat2,
    intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_subtract in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu8_subtract( const uint8_mat *mat1, const uint8_mat *mat2,
  void *intrp_ptr, uint8_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matu8, _subtract, mat1, -, mat2,
    intrp_ptr, result, uint8 );
}


/* Function documented in mat_arit.h */
/* Written by Qin Cai */
/* Edited to generalize for other data types by Vikram Chalana 3/29/1999 */
mutil_errcode matuniv_multiply( const univ_mat *mat1, const univ_mat *mat2,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat1;
  univ_mat      cast_mat2;
  boolean       alloc_cast1;
  boolean       alloc_cast2;

  MUTIL_TRACE("Start matuniv_multiply()");

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

  /* call the appropriate arithmetic function by
     switching on the type of the output matrix */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_multiply( &(cast_mat1.mat.dblmat),
         &(cast_mat2.mat.dblmat), intrp_ptr, &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_multiply( &(cast_mat1.mat.fltmat),
         &(cast_mat2.mat.fltmat), intrp_ptr, &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_multiply( &(cast_mat1.mat.u8mat),
         &(cast_mat2.mat.u8mat), intrp_ptr, &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_multiply( &(cast_mat1.mat.u16mat),
         &(cast_mat2.mat.u16mat), intrp_ptr, &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_multiply( &(cast_mat1.mat.u32mat),
         &(cast_mat2.mat.u32mat), intrp_ptr, &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_multiply( &(cast_mat1.mat.s16mat),
         &(cast_mat2.mat.s16mat), intrp_ptr, &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_multiply( &(cast_mat1.mat.s32mat),
        &(cast_mat2.mat.s32mat), intrp_ptr, &(result->mat.s32mat) );
       break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrices */
  MUTIL_FREE_STANDARD_CASTING( cast_mat1, alloc_cast1 );
  MUTIL_FREE_STANDARD_CASTING( cast_mat2, alloc_cast2 );

  if( errcode ) {
    return errcode;
  }

  MUTIL_TRACE("matuniv_multiply() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_multiply in mat_arit.h */
/* Re-written by Vikram Chalana as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_multiply( const double_mat *mat1, const double_mat *mat2,
  void *intrp_ptr, double_mat *result )
{
  TMPL_MAT_MULT( matdbl, mat1, mat2, intrp_ptr, result, double );
}


/* This function is documented under matuniv_multiply in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matflt_multiply( const float_mat *mat1, const float_mat *mat2,
  void *intrp_ptr, float_mat *result )
{
  TMPL_MAT_MULT( matflt, mat1, mat2, intrp_ptr, result, float );
}


/* This function is documented under matuniv_multiply in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode mats32_multiply( const sint32_mat *mat1, const sint32_mat *mat2,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_MAT_MULT( mats32, mat1, mat2, intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_multiply in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode mats16_multiply( const sint16_mat *mat1, const sint16_mat *mat2,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_MAT_MULT( mats16, mat1, mat2, intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_multiply in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu32_multiply( const uint32_mat *mat1, const uint32_mat *mat2,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_MAT_MULT( matu32, mat1, mat2, intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_multiply in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu16_multiply( const uint16_mat *mat1, const uint16_mat *mat2,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_MAT_MULT( matu16, mat1, mat2, intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_multiply in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu8_multiply( const uint8_mat *mat1, const uint8_mat *mat2,
  void *intrp_ptr, uint8_mat *result )
{
  TMPL_MAT_MULT( matu8, mat1, mat2, intrp_ptr, result, uint8 );
}


/* Function documented in mat_arit.h */
/* Written by Qin Cai */
/* Edited to generalize for other data types by Vikram Chalana 3/29/1999 */
mutil_errcode matuniv_multiply_elem( const univ_mat *mat1,
  const univ_mat *mat2, void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat1;
  univ_mat      cast_mat2;
  boolean       alloc_cast1;
  boolean       alloc_cast2;

  MUTIL_TRACE("Start matuniv_multiply_elem()");

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

  /* call the appropriate arithmetic function by
     switching on the type of the output matrix */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_multiply_elem( &(cast_mat1.mat.dblmat),
        &(cast_mat2.mat.dblmat), intrp_ptr, &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_multiply_elem( &(cast_mat1.mat.fltmat),
        &(cast_mat2.mat.fltmat), intrp_ptr, &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_multiply_elem( &(cast_mat1.mat.u8mat),
        &(cast_mat2.mat.u8mat), intrp_ptr, &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_multiply_elem( &(cast_mat1.mat.u16mat),
        &(cast_mat2.mat.u16mat), intrp_ptr, &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_multiply_elem( &(cast_mat1.mat.u32mat),
        &(cast_mat2.mat.u32mat), intrp_ptr, &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_multiply_elem( &(cast_mat1.mat.s16mat),
        &(cast_mat2.mat.s16mat), intrp_ptr, &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_multiply_elem( &(cast_mat1.mat.s32mat),
        &(cast_mat2.mat.s32mat), intrp_ptr, &(result->mat.s32mat) );
       break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrices */
  MUTIL_FREE_STANDARD_CASTING( cast_mat1, alloc_cast1 );
  MUTIL_FREE_STANDARD_CASTING( cast_mat2, alloc_cast2 );

  if( errcode ) {
    return errcode;
  }

  MUTIL_TRACE("matuniv_multiply_elem() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_multiply_elem in mat_arit.h */
/* Re-written by Vikram Chalana as macro derived from original
   function by Qin Cai. */
mutil_errcode matdbl_multiply_elem( const double_mat *mat1, const double_mat *mat2,
  void *intrp_ptr, double_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matdbl, _multiply_elem, mat1, *, mat2,
    intrp_ptr, result, double );
}


/* This function is documented under matuniv_multiply_elem in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matflt_multiply_elem( const float_mat *mat1, const float_mat *mat2,
  void *intrp_ptr, float_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matflt, _multiply_elem, mat1, *, mat2,
    intrp_ptr, result, float );
}


/* This function is documented under matuniv_multiply_elem in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode mats32_multiply_elem( const sint32_mat *mat1, const sint32_mat *mat2,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( mats32, _multiply_elem, mat1, *, mat2,
    intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_multiply_elem in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode mats16_multiply_elem( const sint16_mat *mat1, const sint16_mat *mat2,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( mats16, _multiply_elem, mat1, *, mat2,
    intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_multiply_elem in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu32_multiply_elem( const uint32_mat *mat1, const uint32_mat *mat2,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matu32, _multiply_elem, mat1, *, mat2,
    intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_multiply_elem in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu16_multiply_elem( const uint16_mat *mat1, const uint16_mat *mat2,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matu16, _multiply_elem, mat1, *, mat2,
    intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_multiply_elem in mat_arit.h */
/* Written by Vikram Chalana */
mutil_errcode matu8_multiply_elem( const uint8_mat *mat1, const uint8_mat *mat2,
  void *intrp_ptr, uint8_mat *result )
{
  TMPL_ALU_MAT_OPERATION_MAT( matu8, _multiply_elem, mat1, *, mat2,
    intrp_ptr, result, uint8 );
}


/* Function documented in mat_arit.h */
/* Written by Qin Cai */
/* Edited to generalize for other data types by Vikram Chalana 3/31/1999 */
mutil_errcode matuniv_multiply_scalar( const univ_mat *mat, univ_scalar scalar,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat;
  boolean       alloc_cast;
  double        tmp;

  MUTIL_TRACE("Start matuniv_multiply_scalar()");

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do standard casting */

  MUTIL_CHECK_STANDARD_CASTING_SCALAR( scalar, result->type, errcode, tmp );
  if( errcode ) {
    return errcode;
  }

  MUTIL_DO_STANDARD_CASTING( mat, result->type, cast_mat,
    intrp_ptr, alloc_cast, errcode );
  if( errcode ) {
    return errcode;
  }


  /* call the appropriate arithmetic function by
     switching on the type of the output data */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_multiply_scalar( &(cast_mat.mat.dblmat),
        SCAUNIV_CAST( scalar, double ), intrp_ptr, &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_multiply_scalar( &(cast_mat.mat.fltmat),
        SCAUNIV_CAST( scalar, float ), intrp_ptr, &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_multiply_scalar( &(cast_mat.mat.u8mat),
        SCAUNIV_CAST( scalar, uint8 ), intrp_ptr, &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_multiply_scalar( &(cast_mat.mat.u16mat),
        SCAUNIV_CAST( scalar, uint16 ), intrp_ptr, &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_multiply_scalar( &(cast_mat.mat.u32mat),
        SCAUNIV_CAST( scalar, uint32 ), intrp_ptr, &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_multiply_scalar( &(cast_mat.mat.s16mat),
        SCAUNIV_CAST( scalar, sint16 ), intrp_ptr, &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_multiply_scalar( &(cast_mat.mat.s32mat),
        SCAUNIV_CAST( scalar, sint32 ), intrp_ptr, &(result->mat.s32mat) );
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

  MUTIL_TRACE("matuniv_multiply_scalar() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_multiply_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matdbl_multiply_scalar( const double_mat *mat,
  double scalar, void *intrp_ptr, double_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( matdbl, _multiply_scalar, mat, *, scalar,
    intrp_ptr, result, double );
}


/* This function is documented under matuniv_multiply_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matflt_multiply_scalar( const float_mat *mat,
  float scalar, void *intrp_ptr, float_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( matflt, _multiply_scalar, mat, *, scalar,
    intrp_ptr, result, float );
}


/* This function is documented under matuniv_multiply_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode mats32_multiply_scalar( const sint32_mat *mat,
  sint32 scalar, void *intrp_ptr, sint32_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( mats32, _multiply_scalar, mat, *, scalar,
    intrp_ptr, result, sint32 );
}


/* This function is documented under matuniv_multiply_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode mats16_multiply_scalar( const sint16_mat *mat,
  sint16 scalar, void *intrp_ptr, sint16_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( mats16, _multiply_scalar, mat, *, scalar,
    intrp_ptr, result, sint16 );
}


/* This function is documented under matuniv_multiply_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu32_multiply_scalar( const uint32_mat *mat,
  uint32 scalar, void *intrp_ptr, uint32_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( matu32, _multiply_scalar, mat, *, scalar,
    intrp_ptr, result, uint32 );
}


/* This function is documented under matuniv_multiply_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu16_multiply_scalar( const uint16_mat *mat,
  uint16 scalar, void *intrp_ptr, uint16_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( matu16, _multiply_scalar, mat, *, scalar,
    intrp_ptr, result, uint16 );
}


/* This function is documented under matuniv_multiply_scalar in mat_arit.h */
/* Written by Qin Cai; Modified by Vikram Chalana to use templates */
mutil_errcode matu8_multiply_scalar( const uint8_mat *mat,
  uint8 scalar, void *intrp_ptr, uint8_mat *result )
{
  TMPL_ALU_MAT_OPERATION_SCALAR( matu8, _multiply_scalar, mat, *, scalar,
    intrp_ptr, result, uint8 );
}


/* Function documented in mat_arit.h */
/* written by Alan Gibbs, Jill Goldschneider */
mutil_errcode matuniv_divide_scalar( const univ_mat *mat, univ_scalar scalar,
  boolean mat_numerator, void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat;
  boolean       alloc_cast;
  double        tmp;

  MUTIL_TRACE("Start matuniv_divide_scalar()");

  if( !mat || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do standard casting */

  MUTIL_CHECK_STANDARD_CASTING_SCALAR( scalar, result->type, errcode, tmp );
  if( errcode ) return errcode;

  MUTIL_DO_STANDARD_CASTING( mat, result->type, cast_mat,
    intrp_ptr, alloc_cast, errcode );
  if( errcode ) return errcode;


  /* call the appropriate arithmetic function by
     switching on the type of the output data */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_divide_scalar( &(cast_mat.mat.dblmat),
        SCAUNIV_CAST( scalar, double ), mat_numerator, intrp_ptr,
        &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_divide_scalar( &(cast_mat.mat.fltmat),
        SCAUNIV_CAST( scalar, float ), mat_numerator, intrp_ptr,
        &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_divide_scalar( &(cast_mat.mat.u8mat),
        SCAUNIV_CAST( scalar, uint8 ), mat_numerator, intrp_ptr,
        &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_divide_scalar( &(cast_mat.mat.u16mat),
        SCAUNIV_CAST( scalar, uint16 ), mat_numerator, intrp_ptr,
        &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_divide_scalar( &(cast_mat.mat.u32mat),
        SCAUNIV_CAST( scalar, uint32 ), mat_numerator, intrp_ptr,
        &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_divide_scalar( &(cast_mat.mat.s16mat),
        SCAUNIV_CAST( scalar, sint16 ), mat_numerator, intrp_ptr,
        &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_divide_scalar( &(cast_mat.mat.s32mat),
        SCAUNIV_CAST( scalar, sint32 ), mat_numerator, intrp_ptr,
        &(result->mat.s32mat) );
       break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrix */
  MUTIL_FREE_STANDARD_CASTING( cast_mat, alloc_cast );

  if( errcode ) return errcode;

  MUTIL_TRACE("matuniv_divide_scalar() done");
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_divide_scalar in mat_arit.h */
/* written by Alan Gibbs, Jill Goldschneider */
mutil_errcode matdbl_divide_scalar( const double_mat *mat,
  double scalar, boolean mat_numerator, void *intrp_ptr, double_mat *result )
{
  TMPL_MAT_SCALAR_DIVIDE( matdbl, _divide_scalar, scalar, mat,
    mat_numerator, intrp_ptr, result, double )
}


/* This function is documented under matuniv_divide_scalar in mat_arit.h */
/* written by Alan Gibbs, Jill Goldschneider */
mutil_errcode matflt_divide_scalar( const float_mat *mat,
  float scalar, boolean mat_numerator, void *intrp_ptr, float_mat *result )
{
  TMPL_MAT_SCALAR_DIVIDE( matflt, _divide_scalar, scalar, mat,
    mat_numerator, intrp_ptr, result, float )
}


/* This function is documented under matuniv_divide_scalar in mat_arit.h */
/* written by Alan Gibbs, Jill Goldschneider */
mutil_errcode mats32_divide_scalar( const sint32_mat *mat,
  sint32 scalar, boolean mat_numerator, void *intrp_ptr, sint32_mat *result )
{
  TMPL_MAT_SCALAR_DIVIDE( mats32, _divide_scalar, scalar, mat,
    mat_numerator, intrp_ptr, result, sint32 )
}


/* This function is documented under matuniv_divide_scalar in mat_arit.h */
/* written by Alan Gibbs, Jill Goldschneider */
mutil_errcode mats16_divide_scalar( const sint16_mat *mat,
  sint16 scalar, boolean mat_numerator, void *intrp_ptr, sint16_mat *result )
{
  TMPL_MAT_SCALAR_DIVIDE( mats16, _divide_scalar, scalar, mat,
    mat_numerator, intrp_ptr, result, sint16 )
}


/* This function is documented under matuniv_divide_scalar in mat_arit.h */
/* written by Alan Gibbs, Jill Goldschneider */
mutil_errcode matu32_divide_scalar( const uint32_mat *mat,
  uint32 scalar, boolean mat_numerator, void *intrp_ptr, uint32_mat *result )
{
  TMPL_MAT_SCALAR_DIVIDE( matu32, _divide_scalar, scalar, mat,
    mat_numerator, intrp_ptr, result, uint32 )
}


/* This function is documented under matuniv_divide_scalar in mat_arit.h */
/* written by Alan Gibbs, Jill Goldschneider */
mutil_errcode matu16_divide_scalar( const uint16_mat *mat,
  uint16 scalar, boolean mat_numerator, void *intrp_ptr, uint16_mat *result )
{
  TMPL_MAT_SCALAR_DIVIDE( matu16, _divide_scalar, scalar, mat,
    mat_numerator, intrp_ptr, result, uint16 )
}


/* This function is documented under matuniv_divide_scalar in mat_arit.h */
/* written by Alan Gibbs, Jill Goldschneider */
mutil_errcode matu8_divide_scalar( const uint8_mat *mat,
  uint8 scalar, boolean mat_numerator, void *intrp_ptr, uint8_mat *result )
{
  TMPL_MAT_SCALAR_DIVIDE( matu8, _divide_scalar, scalar, mat,
    mat_numerator, intrp_ptr, result, uint8 )
}


/* Function documented in mat_arit.h */
/* written by Jill Goldschneider */
mutil_errcode matuniv_divide_elem( const univ_mat *mat1, const univ_mat *mat2,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode errcode;
  univ_mat      cast_mat1;
  univ_mat      cast_mat2;
  boolean       alloc_cast1;
  boolean       alloc_cast2;

  MUTIL_TRACE( "Start matuniv_divide_elem()" );

  if( !mat1 || !mat2 || !result ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do standard casting */

  MUTIL_DO_STANDARD_CASTING( mat1, result->type, cast_mat1,
    intrp_ptr, alloc_cast1, errcode );
  if( errcode ) return errcode;

  MUTIL_DO_STANDARD_CASTING( mat2, result->type, cast_mat2,
    intrp_ptr, alloc_cast2, errcode );
  if( errcode ) {
    MUTIL_FREE_STANDARD_CASTING( cast_mat1, alloc_cast1 );
    return errcode;
  }

  /* call the appropriate arithmetic function by
     switching on the type of the output data */
  switch( result->type ) {
     case MUTIL_DOUBLE:
       errcode = matdbl_divide_elem( &(cast_mat1.mat.dblmat),
        &(cast_mat2.mat.dblmat), intrp_ptr,
        &(result->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_divide_elem( &(cast_mat1.mat.fltmat),
        &(cast_mat2.mat.fltmat), intrp_ptr,
        &(result->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       errcode = matu8_divide_elem( &(cast_mat1.mat.u8mat),
        &(cast_mat2.mat.u8mat), intrp_ptr,
        &(result->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       errcode = matu16_divide_elem( &(cast_mat1.mat.u16mat),
        &(cast_mat2.mat.u16mat), intrp_ptr,
        &(result->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       errcode = matu32_divide_elem( &(cast_mat1.mat.u32mat),
        &(cast_mat2.mat.u32mat), intrp_ptr,
        &(result->mat.u32mat) );
       break;

     case MUTIL_SINT16:
       errcode = mats16_divide_elem( &(cast_mat1.mat.s16mat),
        &(cast_mat2.mat.s16mat), intrp_ptr,
        &(result->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_divide_elem( &(cast_mat1.mat.s32mat),
        &(cast_mat2.mat.s32mat), intrp_ptr,
        &(result->mat.s32mat) );
       break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrix */
  MUTIL_FREE_STANDARD_CASTING( cast_mat1, alloc_cast1 );
  MUTIL_FREE_STANDARD_CASTING( cast_mat2, alloc_cast2 );

  if( errcode ) return errcode;

  MUTIL_TRACE( "Done with matuniv_divide_elem()" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_divide_elem in mat_arit.h */
/* written by Jill Goldschneider */
mutil_errcode matdbl_divide_elem( const double_mat *mat1,
  const double_mat *mat2, void *intrp_ptr,
  double_mat *result )
{
  TMPL_MAT_DIVIDE_MAT_ELEM( matdbl, _divide_elem, mat1, mat2,
    intrp_ptr, result, double )
}


/* This function is documented under matuniv_divide_elem in mat_arit.h */
/* written by Jill Goldschneider */
mutil_errcode matflt_divide_elem( const float_mat *mat1,
  const float_mat *mat2, void *intrp_ptr,
  float_mat *result )
{
  TMPL_MAT_DIVIDE_MAT_ELEM( matflt, _divide_elem, mat1, mat2,
    intrp_ptr, result, float )
}


/* This function is documented under matuniv_divide_elem in mat_arit.h */
/* written by Jill Goldschneider */
mutil_errcode mats32_divide_elem( const sint32_mat *mat1,
  const sint32_mat *mat2, void *intrp_ptr,
  sint32_mat *result )
{
  TMPL_MAT_DIVIDE_MAT_ELEM( mats32, _divide_elem, mat1, mat2,
    intrp_ptr, result, sint32 )
}


/* This function is documented under matuniv_divide_elem in mat_arit.h */
/* written by Jill Goldschneider */
mutil_errcode mats16_divide_elem( const sint16_mat *mat1,
  const sint16_mat *mat2, void *intrp_ptr,
  sint16_mat *result )
{
  TMPL_MAT_DIVIDE_MAT_ELEM( mats16, _divide_elem, mat1, mat2,
    intrp_ptr, result, sint16 )
}


/* This function is documented under matuniv_divide_elem in mat_arit.h */
/* written by Jill Goldschneider */
mutil_errcode matu32_divide_elem( const uint32_mat *mat1,
  const uint32_mat *mat2, void *intrp_ptr,
  uint32_mat *result )
{
  TMPL_MAT_DIVIDE_MAT_ELEM( matu32, _divide_elem, mat1, mat2,
    intrp_ptr, result, uint32 )
}


/* This function is documented under matuniv_divide_elem in mat_arit.h */
/* written by Jill Goldschneider */
mutil_errcode matu16_divide_elem( const uint16_mat *mat1,
  const uint16_mat *mat2, void *intrp_ptr,
  uint16_mat *result )
{
  TMPL_MAT_DIVIDE_MAT_ELEM( matu16, _divide_elem, mat1, mat2,
    intrp_ptr, result, uint16 )
}


/* This function is documented under matuniv_divide_elem in mat_arit.h */
/* written by Jill Goldschneider */
mutil_errcode matu8_divide_elem( const uint8_mat *mat1,
  const uint8_mat *mat2, void *intrp_ptr,
  uint8_mat *result )
{
  TMPL_MAT_DIVIDE_MAT_ELEM( matu8, _divide_elem, mat1, mat2,
    intrp_ptr, result, uint8 )
}


/* Function documented in mat_arit.h */
/* Written by Jill Goldschneider */
mutil_errcode matuniv_abs( const univ_mat *mat_in,
  void *intrp_ptr, univ_mat *mat_out )
{
  mutil_errcode errcode;

  MUTIL_TRACE( "Start matuniv_abs()" );

  if( !mat_in || !mat_out ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* matrices must have same data types */
  if ( !MATUNIV_CHECK_TYPE( mat_in, mat_out ) ) {
    MUTIL_ERROR( "Data types of operand and result must be the same" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* call the appropriate logical function by
     switching on the type of the output data */
  switch( mat_out->type ) {
     case MUTIL_UINT8:
     case MUTIL_UINT16:
     case MUTIL_UINT32:
       MUTIL_WARN( "unsigned integer data type results in a copy" );
       errcode = matuniv_assign( mat_in, intrp_ptr, mat_out );
       break;

     case MUTIL_SINT16:
       errcode = mats16_abs( &(mat_in->mat.s16mat), intrp_ptr,
         &(mat_out->mat.s16mat) );
       break;

     case MUTIL_SINT32:
       errcode = mats32_abs( &(mat_in->mat.s32mat), intrp_ptr,
         &(mat_out->mat.s32mat) );
       break;

     case MUTIL_FLOAT:
       errcode = matflt_abs( &(mat_in->mat.fltmat), intrp_ptr,
         &(mat_out->mat.fltmat) );
       break;

     case MUTIL_DOUBLE:
       errcode = matdbl_abs( &(mat_in->mat.dblmat), intrp_ptr,
         &(mat_out->mat.dblmat) );
       break;

     default:
       MUTIL_ERROR( "This matrix type is currently unsupported" );
       errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( errcode ) return errcode;

  MUTIL_TRACE( "Done with matuniv_abs()" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_abs in mat_arit.h */
/* Written by Jill Goldschneider */
mutil_errcode matdbl_abs( const double_mat *mat_in,
 void *intrp_ptr, double_mat *mat_out )
{
  TMPL_ALU_OPERATION_MAT( matdbl, _abs, MUTIL_ABS, mat_in,
    intrp_ptr, mat_out, double );
}


/* This function is documented under matuniv_abs in mat_arit.h */
/* Written by Jill Goldschneider */
mutil_errcode matflt_abs( const float_mat *mat_in,
 void *intrp_ptr, float_mat *mat_out )
{
  TMPL_ALU_OPERATION_MAT( matflt, _abs, MUTIL_ABS, mat_in,
    intrp_ptr, mat_out, float );
}


/* This function is documented under matuniv_abs in mat_arit.h */
/* Written by Jill Goldschneider */
mutil_errcode mats32_abs( const sint32_mat *mat_in,
 void *intrp_ptr, sint32_mat *mat_out )
{
  TMPL_ALU_OPERATION_MAT( mats32, _abs, MUTIL_ABS, mat_in,
    intrp_ptr, mat_out, sint32 );
}


/* This function is documented under matuniv_abs in mat_arit.h */
/* Written by Jill Goldschneider */
mutil_errcode mats16_abs( const sint16_mat *mat_in,
  void *intrp_ptr, sint16_mat *mat_out )
{
  TMPL_ALU_OPERATION_MAT( mats16, _abs, MUTIL_ABS, mat_in,
    intrp_ptr, mat_out, sint16 );
}


/* function documented in mat_arit.h */
/* written by Krzysztof Koperski */
mutil_errcode matuniv_rescale( const univ_mat *mat_in,
  const univ_scalar min_val, const univ_scalar max_val,
  void *intrp_ptr, univ_mat *mat_out )
{

  mutil_errcode trouble;

  MUTIL_TRACE( "Start matuniv_rescale()" );

  if ( !mat_in || !mat_out ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
     return MUTIL_ERR_NULL_POINTER;
  }

  /* matrices must have same data types */
  if ( !MATUNIV_CHECK_TYPE( mat_in, mat_out ) ||
       !MATUNIV_CHECK_TYPE( &max_val, mat_out ) ||
       !MATUNIV_CHECK_TYPE( &min_val, mat_out ) ) {
    MUTIL_ERROR( "Data types of operands and result must be the same" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch ( mat_out->type ) {

    case MUTIL_DOUBLE:
      trouble = matdbl_rescale( &(mat_in->mat.dblmat),
        SCAUNIV_CAST( min_val, double ), SCAUNIV_CAST( max_val, double ),
        intrp_ptr, &(mat_out->mat.dblmat) );
       break;

    case MUTIL_FLOAT:
      trouble = matflt_rescale( &(mat_in->mat.fltmat),
        SCAUNIV_CAST( min_val, float ), SCAUNIV_CAST( max_val, float ),
        intrp_ptr, &(mat_out->mat.fltmat) );
       break;

    case MUTIL_UINT8:
      trouble = matu8_rescale( &(mat_in->mat.u8mat),
        SCAUNIV_CAST( min_val, uint8 ), SCAUNIV_CAST( max_val, uint8 ),
        intrp_ptr, &(mat_out->mat.u8mat) );
      break;

    case MUTIL_UINT16:
      trouble = matu16_rescale( &(mat_in->mat.u16mat),
        SCAUNIV_CAST( min_val, uint16 ), SCAUNIV_CAST( max_val, uint16 ),
        intrp_ptr, &(mat_out->mat.u16mat) );
      break;

    case MUTIL_UINT32:
      trouble = matu32_rescale( &(mat_in->mat.u32mat),
        SCAUNIV_CAST( min_val, uint32 ), SCAUNIV_CAST( max_val, uint32 ),
        intrp_ptr, &(mat_out->mat.u32mat) );
      break;

    case MUTIL_SINT16:
      trouble = mats16_rescale( &(mat_in->mat.s16mat),
        SCAUNIV_CAST( min_val, sint16 ), SCAUNIV_CAST( max_val, sint16 ),
        intrp_ptr, &(mat_out->mat.s16mat) );
      break;

    case MUTIL_SINT32:
      trouble = mats32_rescale( &(mat_in->mat.s32mat),
        SCAUNIV_CAST( min_val, sint32 ), SCAUNIV_CAST( max_val, sint32 ),
        intrp_ptr, &(mat_out->mat.s32mat) );
      break;

    default:
      MUTIL_ERROR("This matrix type is currently unsupported");
      trouble = MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( trouble ) return trouble;

  MUTIL_TRACE("matuniv_rescale() done");

  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_rescale in mat_arit.h */
/* Written by Krzysztof Koperski. */
mutil_errcode matdbl_rescale( const double_mat *mat_in,
  const double min_val, const double max_val,
  void *intrp_ptr, double_mat *mat_out )
{
  TMPL_MAT_RESCALE( matdbl, mat_in, min_val, max_val,
    intrp_ptr, mat_out, double );
}


/* This function is documented under matuniv_rescale in mat_arit.h */
/* Written by Krzysztof Koperski. */
mutil_errcode matflt_rescale( const float_mat *mat_in,
  const float min_val, const float max_val,
  void *intrp_ptr, float_mat *mat_out )
{
  TMPL_MAT_RESCALE( matflt, mat_in, min_val, max_val,
    intrp_ptr, mat_out, float );
}


/* This function is documented under matuniv_rescale in mat_arit.h */
/* Written by Krzysztof Koperski. */
mutil_errcode matu8_rescale( const uint8_mat *mat_in,
  const uint8 min_val, const uint8 max_val,
  void *intrp_ptr, uint8_mat *mat_out )
{
  TMPL_MAT_RESCALE( matu8, mat_in, min_val, max_val,
    intrp_ptr, mat_out, uint8 );
}


/* This function is documented under matuniv_rescale in mat_arit.h */
/* Written by Krzysztof Koperski. */
mutil_errcode matu16_rescale( const uint16_mat *mat_in,
  const uint16 min_val, const uint16 max_val,
  void *intrp_ptr, uint16_mat *mat_out )
{
  TMPL_MAT_RESCALE( matu16, mat_in, min_val, max_val,
    intrp_ptr, mat_out, uint16 );
}


/* This function is documented under matuniv_rescale in mat_arit.h */
/* Written by Krzysztof Koperski. */
mutil_errcode matu32_rescale( const uint32_mat *mat_in,
  const uint32 min_val, const uint32 max_val,
  void *intrp_ptr, uint32_mat *mat_out )
{
  TMPL_MAT_RESCALE( matu32, mat_in, min_val, max_val,
    intrp_ptr, mat_out, uint32 );
}


/* This function is documented under matuniv_rescale in mat_arit.h */
/* Written by Krzysztof Koperski. */
mutil_errcode mats16_rescale( const sint16_mat *mat_in,
  const sint16 min_val, const sint16 max_val,
  void *intrp_ptr, sint16_mat *mat_out )
{
  TMPL_MAT_RESCALE( mats16, mat_in, min_val, max_val,
    intrp_ptr, mat_out, sint16 );
}


/* This function is documented under matuniv_rescale in mat_arit.h */
/* Written by Krzysztof Koperski. */
mutil_errcode mats32_rescale( const sint32_mat *mat_in,
  const sint32 min_val, const sint32 max_val,
  void *intrp_ptr, sint32_mat *mat_out )
{
  TMPL_MAT_RESCALE( mats32, mat_in, min_val, max_val,
    intrp_ptr, mat_out, sint32 );
}
