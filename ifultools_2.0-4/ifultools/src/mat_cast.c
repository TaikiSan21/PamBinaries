
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_cast.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "mat_cast.h"

#include "mat_assn.h"

#include "mat_any.h"
#include "mat_univ.h"

#include "ut_limit.h"
#include "ut_debug.h"
#include "ut_math.h"
#include "ut_intrp.h"


/* This file contains implementations of functions in mat_cast.h
   which are functions for casting one matrix type to another.
*/


/* function documented in mat_cast.h */
/* written by Jennifer Hodgdon */
mutil_errcode matuniv_cast( const univ_mat *mat, boolean clip,
  void *intrp_ptr, univ_mat *result )
{
  mutil_errcode trouble;

  MUTIL_TRACE( "Start matuniv_cast()" );

  /* avoid lint warning */
  (void) whatssi;

  if(!mat || !result ) {
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if( mat->type == result->type ) {
    trouble = matuniv_assign( mat, intrp_ptr, result );
    if( trouble ) return trouble;
  }
  else { /* types are different */
    switch( mat->type ) {
      case MUTIL_DOUBLE:
        switch( result->type ) {

          case MUTIL_FLOAT:
            trouble = matdbl_cast_to_flt( &(mat->mat.dblmat),
              clip, intrp_ptr, &(result->mat.fltmat) );
            break;

          case MUTIL_DCOMPLEX:
            trouble = matdbl_cast_to_cpx( &(mat->mat.dblmat),
              intrp_ptr, &(result->mat.cpxmat) );
            break;

          case MUTIL_UINT8:
            trouble = matdbl_cast_to_u8( &(mat->mat.dblmat),
              clip, intrp_ptr, &(result->mat.u8mat) );
            break;

          case MUTIL_UINT16:
            trouble = matdbl_cast_to_u16( &(mat->mat.dblmat),
              clip, intrp_ptr, &(result->mat.u16mat) );
            break;

          case MUTIL_UINT32:
            trouble = matdbl_cast_to_u32( &(mat->mat.dblmat),
              clip, intrp_ptr, &(result->mat.u32mat) );
            break;

          case MUTIL_SINT16:
            trouble = matdbl_cast_to_s16( &(mat->mat.dblmat),
              clip, intrp_ptr, &(result->mat.s16mat) );
            break;

          case MUTIL_SINT32:
            trouble = matdbl_cast_to_s32( &(mat->mat.dblmat),
              clip, intrp_ptr, &(result->mat.s32mat) );
            break;

          default:
            MUTIL_ERROR( "This combination of matrix types is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
        }

        break; /* end of mat->type is double */

      case MUTIL_FLOAT:
        switch( result->type ) {

          case MUTIL_DOUBLE:
            trouble = matflt_cast_to_dbl( &(mat->mat.fltmat),
              intrp_ptr, &(result->mat.dblmat) );
            break;

          case MUTIL_DCOMPLEX:
            trouble = matflt_cast_to_cpx( &(mat->mat.fltmat),
              intrp_ptr, &(result->mat.cpxmat) );
            break;

          case MUTIL_UINT8:
            trouble = matflt_cast_to_u8( &(mat->mat.fltmat),
              clip, intrp_ptr, &(result->mat.u8mat) );
            break;

          case MUTIL_UINT16:
            trouble = matflt_cast_to_u16( &(mat->mat.fltmat),
              clip, intrp_ptr, &(result->mat.u16mat) );
            break;

          case MUTIL_UINT32:
            trouble = matflt_cast_to_u32( &(mat->mat.fltmat),
              clip, intrp_ptr, &(result->mat.u32mat) );
            break;

          case MUTIL_SINT16:
            trouble = matflt_cast_to_s16( &(mat->mat.fltmat),
              clip, intrp_ptr, &(result->mat.s16mat) );
            break;

          case MUTIL_SINT32:
            trouble = matflt_cast_to_s32( &(mat->mat.fltmat),
              clip, intrp_ptr, &(result->mat.s32mat) );
            break;

          default:
            MUTIL_ERROR( "This combination of matrix types is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
        }

        break; /* end of mat->type is float */

      case MUTIL_DCOMPLEX:

        switch( result->type ) {

          case MUTIL_DOUBLE:
            trouble = matcpx_cast_to_dbl( &(mat->mat.cpxmat),
              clip, intrp_ptr, &(result->mat.dblmat) );
            break;

          default:
            MUTIL_ERROR( "This combination of matrix types is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
        }

        break; /* end of mat->type is complex */

      case MUTIL_UINT8:

        switch( result->type ) {

          case MUTIL_DOUBLE:
            trouble = matu8_cast_to_dbl( &(mat->mat.u8mat),
              intrp_ptr, &(result->mat.dblmat) );
            break;

          case MUTIL_FLOAT:
            trouble = matu8_cast_to_flt( &(mat->mat.u8mat),
              intrp_ptr, &(result->mat.fltmat) );
            break;

          case MUTIL_UINT16:
            trouble = matu8_cast_to_u16( &(mat->mat.u8mat),
              intrp_ptr, &(result->mat.u16mat) );
            break;

          case MUTIL_UINT32:
            trouble = matu8_cast_to_u32( &(mat->mat.u8mat),
              intrp_ptr, &(result->mat.u32mat) );
            break;

          case MUTIL_SINT16:
            trouble = matu8_cast_to_s16( &(mat->mat.u8mat),
              intrp_ptr, &(result->mat.s16mat) );
            break;

          case MUTIL_SINT32:
            trouble = matu8_cast_to_s32( &(mat->mat.u8mat),
              intrp_ptr, &(result->mat.s32mat) );
            break;

          default:
            MUTIL_ERROR( "This combination of matrix types is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
        }


        break; /* end of mat->type is uint8 */

      case MUTIL_UINT16:

        switch( result->type ) {

          case MUTIL_DOUBLE:
            trouble = matu16_cast_to_dbl( &(mat->mat.u16mat),
              intrp_ptr, &(result->mat.dblmat) );
            break;

          case MUTIL_FLOAT:
            trouble = matu16_cast_to_flt( &(mat->mat.u16mat),
              intrp_ptr, &(result->mat.fltmat) );
            break;

          case MUTIL_UINT8:
            trouble = matu16_cast_to_u8( &(mat->mat.u16mat), clip,
              intrp_ptr, &(result->mat.u8mat) );
            break;

          case MUTIL_UINT32:
            trouble = matu16_cast_to_u32( &(mat->mat.u16mat),
              intrp_ptr, &(result->mat.u32mat) );
            break;

          case MUTIL_SINT16:
            trouble = matu16_cast_to_s16( &(mat->mat.u16mat), clip,
              intrp_ptr, &(result->mat.s16mat) );
            break;

          case MUTIL_SINT32:
            trouble = matu16_cast_to_s32( &(mat->mat.u16mat),
              intrp_ptr, &(result->mat.s32mat) );
            break;

          default:
            MUTIL_ERROR( "This combination of matrix types is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
        }

        break; /* end of mat->type is uint16 */

      case MUTIL_UINT32:

        switch( result->type ) {

          case MUTIL_DOUBLE:
            trouble = matu32_cast_to_dbl( &(mat->mat.u32mat),
              intrp_ptr, &(result->mat.dblmat) );
            break;

          case MUTIL_FLOAT:
            trouble = matu32_cast_to_flt( &(mat->mat.u32mat),
              intrp_ptr, &(result->mat.fltmat) );
            break;

          case MUTIL_UINT8:
            trouble = matu32_cast_to_u8( &(mat->mat.u32mat), clip,
              intrp_ptr, &(result->mat.u8mat) );
            break;

          case MUTIL_UINT16:
            trouble = matu32_cast_to_u16( &(mat->mat.u32mat), clip,
              intrp_ptr, &(result->mat.u16mat) );
            break;

          case MUTIL_SINT16:
            trouble = matu32_cast_to_s16( &(mat->mat.u32mat), clip,
              intrp_ptr, &(result->mat.s16mat) );
            break;

          case MUTIL_SINT32:
            trouble = matu32_cast_to_s32( &(mat->mat.u32mat), clip,
              intrp_ptr, &(result->mat.s32mat) );
            break;

          default:
            MUTIL_ERROR( "This combination of matrix types is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
        }

        break; /* end of mat->type is uint32 */

      case MUTIL_SINT16:

        switch( result->type ) {

          case MUTIL_DOUBLE:
            trouble = mats16_cast_to_dbl( &(mat->mat.s16mat),
              intrp_ptr, &(result->mat.dblmat) );
            break;

          case MUTIL_FLOAT:
            trouble = mats16_cast_to_flt( &(mat->mat.s16mat),
              intrp_ptr, &(result->mat.fltmat) );
            break;

          case MUTIL_UINT8:
            trouble = mats16_cast_to_u8( &(mat->mat.s16mat), clip,
              intrp_ptr, &(result->mat.u8mat) );
            break;

          case MUTIL_UINT16:
            trouble = mats16_cast_to_u16( &(mat->mat.s16mat), clip,
              intrp_ptr, &(result->mat.u16mat) );
            break;

          case MUTIL_UINT32:
            trouble = mats16_cast_to_u32( &(mat->mat.s16mat), clip,
              intrp_ptr, &(result->mat.u32mat) );
            break;

          case MUTIL_SINT32:
            trouble = mats16_cast_to_s32( &(mat->mat.s16mat),
              intrp_ptr, &(result->mat.s32mat) );
            break;

          default:
            MUTIL_ERROR( "This combination of matrix types is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
        }

        break; /* end of mat->type is sint32 */

      case MUTIL_SINT32:

        switch( result->type ) {

          case MUTIL_DOUBLE:
            trouble = mats32_cast_to_dbl( &(mat->mat.s32mat),
              intrp_ptr, &(result->mat.dblmat) );
            break;

          case MUTIL_FLOAT:
            trouble = mats32_cast_to_flt( &(mat->mat.s32mat),
              intrp_ptr, &(result->mat.fltmat) );
            break;

          case MUTIL_UINT8:
            trouble = mats32_cast_to_u8( &(mat->mat.s32mat), clip,
              intrp_ptr, &(result->mat.u8mat) );
            break;

          case MUTIL_UINT16:
            trouble = mats32_cast_to_u16( &(mat->mat.s32mat), clip,
              intrp_ptr, &(result->mat.u16mat) );
            break;

          case MUTIL_UINT32:
            trouble = mats32_cast_to_u32( &(mat->mat.s32mat), clip,
              intrp_ptr, &(result->mat.u32mat) );
            break;

          case MUTIL_SINT16:
            trouble = mats32_cast_to_s16( &(mat->mat.s32mat), clip,
              intrp_ptr, &(result->mat.s16mat) );
            break;

          default:
            MUTIL_ERROR( "This combination of matrix types is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
        }

        break; /* end of mat->type is sint32 */

      default:
        MUTIL_ERROR( "This matrix type is currently unsupported" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    } /* switch on matrix types */
  } /* if/else on mat and result type being same */

  if( trouble ) return trouble;

  MUTIL_TRACE( "matuniv_cast() done" );
  return MUTIL_ERR_OK;
}


/** Template macro for floating point type to integer cast function.
 * Macro that expands to the body of a function that casts
 * floating point type (double or float) matrices to integer types,
 * such as matdbl\_cast\_to\_u32.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_cast.c
 * @library matrix
 * @param MAT        Matrix pointer to cast (function argument).
 * @param INCODE     3-character code for input type.
 * @param RESULT     Matrix pointer to return (function argument).
 * @param CLIP       Flag saying to clip or return error for out of range
 *           (function argument).
 * @param OUTCODE    3-character code for output type.
 * @param OUTTYPE    Data type for output matrix.
 * @param MINVAL     Minimum legal value for this type.
 * @param MAXVAL     Maximum legal value for this type.
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_cast\_to\_u32 function:
 *     #TMPL_FPT_CAST_TO_INT( mat, flt, result, clip, u32, uint32, 0, MUTIL_UINT32_MAX, intrp_ptr );#
 * @private
 */
#define TMPL_FPT_CAST_TO_INT( MAT, INCODE, RESULT, CLIP, OUTCODE, OUTTYPE, \
  MINVAL, MAXVAL, INTRP_PTR ) \
  mutil_errcode trouble; \
  sint32        i; \
  double        tmp; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Starting mat" #INCODE "_cast_to_" #OUTCODE "()" ); \
  \
  trouble = mat ## INCODE ## _validate( MAT ); \
  if( trouble ) return trouble; \
  trouble = mat ## OUTCODE ## _validate( RESULT ); \
  if( trouble ) return trouble; \
  \
  if( !MATANY_EQUAL_DIM( MAT, RESULT ) ) { \
    MUTIL_ERROR( "Input and output matrices must have equal dimensions" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  for( i = 0; i < (MAT)->nelem; i++ ) { \
    tmp = MUTIL_ROUND( (MAT)->data[i] ); \
    if( tmp < (MINVAL) ) { \
      if( CLIP ) { \
        /* replaced the line below while reviewing */ \
        /* -- a possible cause of errors -- VC 4/27/99 */ \
        /*(RESULT)->data[i] = 0; \ */ \
        (RESULT)->data[i] = (MINVAL); \
      } \
      else { \
        MUTIL_ERROR( "Cannot cast -- value out of range" ); \
        return MUTIL_ERR_OVERFLOW; \
      } \
    } \
    else if( tmp > (MAXVAL) ) { \
      if( CLIP ) { \
        (RESULT)->data[i] = (MAXVAL); \
      } \
      else { \
        MUTIL_ERROR( "Cannot cast -- value out of range" ); \
        return MUTIL_ERR_OVERFLOW; \
      } \
    } \
    else { \
      (RESULT)->data[i] = ( OUTTYPE ) tmp; \
    } \
  } \
  if( MUTIL_INTERRUPT( (MAT)->nelem * 5.0, INTRP_PTR ) ) { \
    MUTIL_ERROR( "User interrupt" ); \
    return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "mat" #INCODE "_cast_to_" #OUTCODE "() done" ); \
  return MUTIL_ERR_OK


/** Template macro for unsigned integer to integer cast function.
 * Macro that expands to the body of a function that casts
 * unsigned integer matrices to other integer types, with range checking
 * only on the upper limit, such as matu32\_cast\_to\_u16.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_cast.c
 * @library matrix
 * @param MAT        Matrix pointer to cast (function argument).
 * @param RESULT     Matrix pointer to return (function argument).
 * @param CLIP       Flag saying to clip or return error for out of range
 *           (function argument).
 * @param INCODE     Three-character code for input type.
 * @param OUTCODE    Three-character code for output type.
 * @param OUTTYPE    Data type of output matrix.
 * @param MAXVAL     Maximum legal output value.
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the mats32\_cast\_to\_u32 function:
 *     #TMPL_UINT_CAST_TO_INT( mat, result, clip, u32, s32, sint32, MUTIL_SINT32_MAX, intrp_ptr );#
 * @private
 */
#define TMPL_UINT_CAST_TO_INT( MAT, RESULT, CLIP, INCODE, OUTCODE, \
  OUTTYPE, MAXVAL, INTRP_PTR ) \
  mutil_errcode trouble; \
  sint32        i; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Starting mat" #INCODE "_cast_to_" #OUTCODE "()" ); \
  \
  trouble = mat ## INCODE ## _validate( MAT ); \
  if( trouble ) return trouble; \
  trouble = mat ## OUTCODE ## _validate( RESULT ); \
  if( trouble ) return trouble; \
  \
  if( !MATANY_EQUAL_DIM( MAT, RESULT ) ) { \
    MUTIL_ERROR( "Input and output matrices must have equal dimensions" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  for( i = 0; i < (MAT)->nelem; i++ ) { \
    if( (MAT)->data[i] > ( MAXVAL ) ) { \
      if( clip ) { \
        ( RESULT )->data[i] = ( MAXVAL ); \
      } \
      else { \
        MUTIL_ERROR( "Cannot cast -- value out of range" ); \
        return MUTIL_ERR_OVERFLOW; \
      } \
    } \
    else { \
      ( RESULT )->data[i] = ( OUTTYPE ) mat->data[i]; \
    } \
  } \
  if( MUTIL_INTERRUPT( (MAT)->nelem * 5.0, INTRP_PTR ) ) { \
    MUTIL_ERROR( "User interrupt" ); \
    return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "mat" #INCODE "_cast_to_" #OUTCODE "() done" ); \
  return MUTIL_ERR_OK


/** Template macro for signed integer to unsigned integer cast function.
 * Macro that expands to the body of a function that casts
 * signed integer matrices to other unsigned integer types, with range checking
 * only on the lower limit, such as mats16\_cast\_to\_u16.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_cast.c
 * @library matrix
 * @param MAT        Matrix pointer to cast (function argument).
 * @param RESULT     Matrix pointer to return (function argument).
 * @param CLIP       Flag saying to clip or return error for out of range
 *           (function argument).
 * @param INCODE     Three-character code for input type.
 * @param OUTCODE    Three-character code for output type.
 * @param OUTTYPE    Data type of output matrix.
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the mats32\_cast\_to\_u32 function:
 *     #TMPL_INT_CAST_TO_UINT( mat, result, clip, u32, s32, sint32, intrp_ptr );#
 * @private
 */
#define TMPL_INT_CAST_TO_UINT( MAT, RESULT, CLIP, INCODE, OUTCODE, \
  OUTTYPE, INTRP_PTR ) \
  mutil_errcode trouble; \
  sint32        i; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Starting mat" #INCODE "_cast_to_" #OUTCODE "()" ); \
  \
  trouble = mat ## INCODE ## _validate( MAT ); \
  if( trouble ) return trouble; \
  trouble = mat ## OUTCODE ## _validate( RESULT ); \
  if( trouble ) return trouble; \
  \
  if( !MATANY_EQUAL_DIM( MAT, RESULT ) ) { \
    MUTIL_ERROR( "Input and output matrices must have equal dimensions" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  for( i = 0; i < (MAT)->nelem; i++ ) { \
    if( (MAT)->data[i] < 0 ) { \
      if( clip ) { \
        ( RESULT )->data[i] = 0; \
      } \
      else { \
        MUTIL_ERROR( "Cannot cast -- value out of range" ); \
        return MUTIL_ERR_OVERFLOW; \
      } \
    } \
    else { \
      ( RESULT )->data[i] = ( OUTTYPE ) mat->data[i]; \
    } \
  } \
  if( MUTIL_INTERRUPT( (MAT)->nelem * 5.0, INTRP_PTR ) ) { \
    MUTIL_ERROR( "User interrupt" ); \
    return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "mat" #INCODE "_cast_to_" #OUTCODE "() done" ); \
  return MUTIL_ERR_OK


/** Template macro for integer to integer cast function.
 * Macro that expands to the body of a function that casts
 * integer matrices to other integer types, with range checking,
 * such as mats32\_cast\_to\_u32.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_cast.c
 * @library matrix
 * @param MAT        Matrix pointer to cast (function argument).
 * @param RESULT     Matrix pointer to return (function argument).
 * @param CLIP       Flag saying to clip or return error for out of range
 *           (function argument).
 * @param INCODE     Three-character code for input type.
 * @param OUTCODE    Three-character code for output type.
 * @param OUTTYPE    Data type of output matrix.
 * @param MINVAL     Minimum legal output value.
 * @param MAXVAL     Maximum legal output value.
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the mats32\_cast\_to\_u32 function:
 *     #TMPL_INT_CAST_TO_INT( mat, result, clip, s32, u32, uint32, 0, MUTIL_UINT32_MAX, intrp_ptr );#
 * @private
 */
#define TMPL_INT_CAST_TO_INT( MAT, RESULT, CLIP, INCODE, OUTCODE, \
  OUTTYPE, MINVAL, MAXVAL, INTRP_PTR ) \
  mutil_errcode trouble; \
  sint32        i; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Starting mat" #INCODE "_cast_to_" #OUTCODE "()" ); \
  \
  trouble = mat ## INCODE ## _validate( MAT ); \
  if( trouble ) return trouble; \
  trouble = mat ## OUTCODE ## _validate( RESULT ); \
  if( trouble ) return trouble; \
  \
  if( !MATANY_EQUAL_DIM( MAT, RESULT ) ) { \
    MUTIL_ERROR( "Input and output matrices must have equal dimensions" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  for( i = 0; i < (MAT)->nelem; i++ ) { \
    if( (MAT)->data[i] < (MINVAL) ) { \
      if( CLIP ) { \
        (RESULT)->data[i] = (MINVAL); \
      } \
      else { \
        MUTIL_ERROR( "Cannot cast -- value out of range" ); \
        return MUTIL_ERR_OVERFLOW; \
      } \
    } \
    else if( (MAT)->data[i] > (MAXVAL) ) { \
      if( CLIP ) { \
        (RESULT)->data[i] = (MAXVAL); \
      } \
      else { \
        MUTIL_ERROR( "Cannot cast -- value out of range" ); \
        return MUTIL_ERR_OVERFLOW; \
      } \
    } \
    else { \
      (RESULT)->data[i] = ( OUTTYPE ) (MAT)->data[i]; \
    } \
  } \
  if( MUTIL_INTERRUPT( (MAT)->nelem * 5.0, INTRP_PTR ) ) { \
    MUTIL_ERROR( "User interrupt" ); \
    return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "mat" #INCODE "_cast_to_" #OUTCODE "() done" ); \
  return MUTIL_ERR_OK


/** Template macro for generic no-test cast function.
 * Macro that expands to the body of a function that casts
 * where there is no chance that overflow can occur.  For example,
 * signed or unsigned integer types to double, such as matu16\_cast\_to\_dbl.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_cast.c
 * @library matrix
 * @param MAT        Matrix pointer to cast (function argument).
 * @param RESULT     Matrix pointer to return (function argument).
 * @param OUTTYPE    Data type for output matrix.
 * @param INCODE     3-character code for input matrix type.
 * @param OUTCODE    3-character code for output matrix type.
 * @param INTRP_PTR  Pointer for interrupt handling (function argument).
 * @usage Body of the matu16\_cast\_to\_dbl function:
 *     #TMPL_NOCHECK_CAST(mat, result, double, u16, dbl, intrp_ptr);#
 * @private
 */
#define TMPL_NOCHECK_CAST( MAT, RESULT, OUTTYPE, INCODE, OUTCODE, INTRP_PTR ) \
  mutil_errcode trouble; \
  sint32        i; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Starting mat" #INCODE "_cast_to_" #OUTCODE "()" ); \
  \
  trouble = mat ## INCODE ## _validate( MAT ); \
  if( trouble ) return trouble; \
  trouble = mat ## OUTCODE ## _validate( RESULT ); \
  if( trouble ) return trouble; \
  \
  if( !MATANY_EQUAL_DIM( MAT, RESULT ) ) { \
    MUTIL_ERROR( "Input and output matrices must have equal dimensions" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  for( i = 0; i < (MAT)->nelem; i++ ) { \
    (RESULT)->data[i] = (OUTTYPE) (MAT)->data[i]; \
  } \
  if( MUTIL_INTERRUPT( (MAT)->nelem * 2.0, INTRP_PTR ) ) { \
    MUTIL_ERROR( "User interrupt" ); \
    return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "mat" #INCODE "_cast_to_" #OUTCODE "() done" ); \
  return MUTIL_ERR_OK


/* function documented in mat_cast.h */
/* written by Jennifer Hodgdon */
mutil_errcode matdbl_cast_to_u8( const double_mat *mat,
  boolean clip, void *intrp_ptr, uint8_mat *result )
{
  TMPL_FPT_CAST_TO_INT( mat, dbl, result, clip, u8, uint8,
    0, MUTIL_UINT8_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jennifer Hodgdon */
mutil_errcode matdbl_cast_to_u16( const double_mat *mat,
  boolean clip, void *intrp_ptr, uint16_mat *result )
{
  TMPL_FPT_CAST_TO_INT( mat, dbl, result, clip, u16, uint16,
    0, MUTIL_UINT16_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jennifer Hodgdon */
mutil_errcode matdbl_cast_to_u32( const double_mat *mat,
  boolean clip, void *intrp_ptr, uint32_mat *result )
{
  TMPL_FPT_CAST_TO_INT( mat, dbl, result, clip, u32, uint32,
    0, MUTIL_UINT32_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matdbl_cast_to_s16( const double_mat *mat,
  boolean clip, void *intrp_ptr, sint16_mat *result )
{
  TMPL_FPT_CAST_TO_INT( mat, dbl, result, clip, s16, sint16,
    MUTIL_SINT16_MIN, MUTIL_SINT16_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jennifer Hodgdon */
mutil_errcode matdbl_cast_to_s32( const double_mat *mat,
  boolean clip, void *intrp_ptr, sint32_mat *result )
{
  TMPL_FPT_CAST_TO_INT( mat, dbl, result, clip, s32, sint32,
    MUTIL_SINT32_MIN, MUTIL_SINT32_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Luca Cazzanti */
mutil_errcode matdbl_cast_to_cpx( const double_mat *mat,
  void *intrp_ptr, dcomplex_mat *result )
{
  mutil_errcode trouble;
  sint32 i;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Starting matdbl_cast_to_cpx()" );

  trouble = matdbl_validate( mat );
  if( trouble ) return trouble;

  trouble = matcpx_validate( result );
  if( trouble ) return trouble;

  if( !MATANY_EQUAL_DIM( mat, result ) ) {
    MUTIL_ERROR( "Input and output matrices must have equal dimensions" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  for( i = 0; i < mat->nelem; i++ ) {
      result->data[i].re = mat->data[i];
      result->data[i].im = 0;
  }

  if( MUTIL_INTERRUPT( mat->nelem * 4.0, intrp_ptr ) ) {
    MUTIL_ERROR( "User interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE( "matdbl_cast_to_cpx() done" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matdbl_cast_to_flt( const double_mat *mat,
  boolean clip, void *intrp_ptr, float_mat *result )
{
  mutil_errcode trouble;
  sint32        i;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Starting matdbl_cast_to_flt()" );

  trouble = matdbl_validate( mat );
  if( trouble ) return trouble;

  trouble = matflt_validate( result );
  if( trouble ) return trouble;

  if( !MATANY_EQUAL_DIM( mat, result ) ) {
    MUTIL_ERROR( "Input and output matrices must have equal dimensions" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  for( i = 0; i < mat->nelem; i++ ) {

    if( mat->data[i] < (-MUTIL_FLOAT_MAX) ) {
      if( clip ) {
        result->data[i] = (float) (-MUTIL_FLOAT_MAX);
      }
      else {
        MUTIL_ERROR( "Cannot cast -- value out of range" );
        return MUTIL_ERR_OVERFLOW;
      }
    }
    else if( mat->data[i] > MUTIL_FLOAT_MAX ) {
      if( clip ) {
        result->data[i] = (float) MUTIL_FLOAT_MAX;
      }
      else {
        MUTIL_ERROR( "Cannot cast -- value out of range" );
        return MUTIL_ERR_OVERFLOW;
      }
    }
    else {
      result->data[i] = ( float ) mat->data[i];
    }
  }

  if( MUTIL_INTERRUPT( mat->nelem * 4.0, intrp_ptr ) ) {
    MUTIL_ERROR( "User interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE( "matdbl_cast_to_flt() done" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_cast.h */
/* written by Luca Cazzanti */
mutil_errcode matcpx_cast_to_dbl( const dcomplex_mat *mat,
  boolean clip, void *intrp_ptr, double_mat *result )
{
  mutil_errcode trouble;
  sint32 i;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Starting matcpx_cast_to_dbl()" );

  trouble = matcpx_validate( mat );
  if( trouble ) return trouble;

  trouble = matdbl_validate( result );
  if( trouble ) return trouble;

  if( !MATANY_EQUAL_DIM( mat, result ) ) {
    MUTIL_ERROR( "Input and output matrices must have equal dimensions" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  for( i = 0; i < mat->nelem; i++ ) {
    if( !clip && mat->data[i].im ) {
      MUTIL_ERROR( "Cannot cast to double -- imaginary part non-zero" );
      return MUTIL_ERR_OVERFLOW;
    }
    result->data[i] = mat->data[i].re;
  }

  if( MUTIL_INTERRUPT( mat->nelem * 4.0, intrp_ptr ) ) {
    MUTIL_ERROR( "User interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE( "matcpx_cast_to_dbl() done" );
  return MUTIL_ERR_OK;
}


mutil_errcode matflt_cast_to_dbl( const float_mat *mat,
  void *intrp_ptr, double_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, double, flt, dbl, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matflt_cast_to_cpx( const float_mat *mat,
  void *intrp_ptr, dcomplex_mat *result )
{
  mutil_errcode trouble;
  sint32 i;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Starting matflt_cast_to_cpx()" );

  trouble = matflt_validate( mat );
  if( trouble ) return trouble;

  trouble = matcpx_validate( result );
  if( trouble ) return trouble;

  if( !MATANY_EQUAL_DIM( mat, result ) ) {
    MUTIL_ERROR( "Input and output matrices must have equal dimensions" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  for( i = 0; i < mat->nelem; i++ ) {
      result->data[i].re = (double) mat->data[i];
      result->data[i].im = 0.0;
  }

  if( MUTIL_INTERRUPT( mat->nelem * 4.0, intrp_ptr ) ) {
    MUTIL_ERROR( "User interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE( "matflt_cast_to_cpx() done" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matflt_cast_to_u8( const float_mat *mat,
  boolean clip, void *intrp_ptr, uint8_mat *result )
{
  TMPL_FPT_CAST_TO_INT( mat, flt, result, clip, u8, uint8,
    0, MUTIL_UINT8_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matflt_cast_to_u16( const float_mat *mat,
  boolean clip, void *intrp_ptr, uint16_mat *result )
{
  TMPL_FPT_CAST_TO_INT( mat, flt, result, clip, u16, uint16,
    0, MUTIL_UINT16_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matflt_cast_to_u32( const float_mat *mat,
  boolean clip, void *intrp_ptr, uint32_mat *result )
{
  TMPL_FPT_CAST_TO_INT( mat, flt, result, clip, u32, uint32,
    0, MUTIL_UINT32_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matflt_cast_to_s16( const float_mat *mat,
  boolean clip, void *intrp_ptr, sint16_mat *result )
{
  TMPL_FPT_CAST_TO_INT( mat, flt, result, clip, s16, sint16,
    MUTIL_SINT16_MIN, MUTIL_SINT16_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matflt_cast_to_s32( const float_mat *mat,
  boolean clip, void *intrp_ptr, sint32_mat *result )
{
  TMPL_FPT_CAST_TO_INT( mat, flt, result, clip, s32, sint32,
    MUTIL_SINT32_MIN, MUTIL_SINT32_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jennifer Hodgdon */
mutil_errcode matu8_cast_to_dbl( const uint8_mat *mat,
  void *intrp_ptr, double_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, double, u8, dbl, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matu8_cast_to_flt( const uint8_mat *mat,
  void *intrp_ptr, float_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, float, u8, flt, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Vikram Chalana */
mutil_errcode matu8_cast_to_u16( const uint8_mat *mat,
  void *intrp_ptr, uint16_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, uint16, u8, u16, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Vikram Chalana */
mutil_errcode matu8_cast_to_u32( const uint8_mat *mat,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, uint32, u8, u32, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matu8_cast_to_s16( const uint8_mat *mat,
  void *intrp_ptr, sint16_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, sint16, u8, s16, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jennifer Hodgdon */
mutil_errcode matu8_cast_to_s32( const uint8_mat *mat,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, sint32, u8, s32, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jennifer Hodgdon */
mutil_errcode matu16_cast_to_dbl( const uint16_mat *mat,
  void *intrp_ptr, double_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, double, u16, dbl, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matu16_cast_to_flt( const uint16_mat *mat,
  void *intrp_ptr, float_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, float, u16, flt, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Vikram Chalana */
mutil_errcode matu16_cast_to_u8( const uint16_mat *mat,
  boolean clip, void *intrp_ptr, uint8_mat *result )
{
  /* LINTED: Meant to convert to 8-bit value here */
  TMPL_UINT_CAST_TO_INT( mat, result, clip, u16, u8, uint8,
    MUTIL_UINT8_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Vikram Chalana */
mutil_errcode matu16_cast_to_u32( const uint16_mat *mat,
  void *intrp_ptr, uint32_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, uint32, u16, u32, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matu16_cast_to_s16( const uint16_mat *mat,
  boolean clip, void *intrp_ptr, sint16_mat *result )
{

  TMPL_UINT_CAST_TO_INT( mat, result, clip, u16, s16, sint16,
    MUTIL_SINT16_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jennifer Hodgdon */
mutil_errcode matu16_cast_to_s32( const uint16_mat *mat,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, sint32, u16, s32, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Vikram Chalana */
mutil_errcode matu32_cast_to_dbl( const uint32_mat *mat,
  void *intrp_ptr, double_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, double, u32, dbl, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matu32_cast_to_flt( const uint32_mat *mat,
  void *intrp_ptr, float_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, float, u32, flt, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Vikram Chalana */
mutil_errcode matu32_cast_to_u8( const uint32_mat *mat,
  boolean clip, void *intrp_ptr, uint8_mat *result )
{
  /* LINTED: Meant to convert to 8-bit value here */
  TMPL_UINT_CAST_TO_INT( mat, result, clip, u32, u8, uint8,
    MUTIL_UINT8_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Vikram Chalana */
mutil_errcode matu32_cast_to_u16( const uint32_mat *mat,
  boolean clip, void *intrp_ptr, uint16_mat *result )
{
  /* LINTED: Meant to convert to 16-bit value here */
  TMPL_UINT_CAST_TO_INT( mat, result, clip, u32, u16, uint16,
    MUTIL_UINT16_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matu32_cast_to_s16( const uint32_mat *mat,
  boolean clip, void *intrp_ptr, sint16_mat *result )
{
  /* LINTED: Meant to convert to 16-bit value here */
  TMPL_UINT_CAST_TO_INT( mat, result, clip, u32, s16, sint16,
    MUTIL_SINT16_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Vikram Chalana */
mutil_errcode matu32_cast_to_s32( const uint32_mat *mat,
  boolean clip, void *intrp_ptr, sint32_mat *result )
{
  TMPL_UINT_CAST_TO_INT( mat, result, clip, u32, s32, sint32,
    MUTIL_SINT32_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode mats16_cast_to_dbl( const sint16_mat *mat,
  void *intrp_ptr, double_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, double, s16, dbl, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode mats16_cast_to_flt( const sint16_mat *mat,
  void *intrp_ptr, float_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, float, s16, flt, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode mats16_cast_to_u8( const sint16_mat *mat,
  boolean clip, void *intrp_ptr, uint8_mat *result )
{
  /* LINTED: Meant to convert to 8-bit value here */
  TMPL_INT_CAST_TO_INT( mat, result, clip, s16, u8, uint8,
    0, MUTIL_UINT8_MAX, intrp_ptr );
}

/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode mats16_cast_to_u16( const sint16_mat *mat,
  boolean clip, void *intrp_ptr, uint16_mat *result )
{
  TMPL_INT_CAST_TO_UINT( mat, result, clip, s16, u16, uint16,
    intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode mats16_cast_to_u32( const sint16_mat *mat,
  boolean clip, void *intrp_ptr, uint32_mat *result )
{
  TMPL_INT_CAST_TO_UINT( mat, result, clip, s16, u32, uint32,
    intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode mats16_cast_to_s32( const sint16_mat *mat,
  void *intrp_ptr, sint32_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, sint32, s16, s32, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jennifer Hodgdon */
mutil_errcode mats32_cast_to_dbl( const sint32_mat *mat,
  void *intrp_ptr, double_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, double, s32, dbl, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode mats32_cast_to_flt( const sint32_mat *mat,
  void *intrp_ptr, float_mat *result )
{
  TMPL_NOCHECK_CAST(mat, result, float, s32, flt, intrp_ptr);
}


/* function documented in mat_cast.h */
/* written by Vikram Chalana */
mutil_errcode mats32_cast_to_u8( const sint32_mat *mat,
  boolean clip, void *intrp_ptr, uint8_mat *result )
{
  /* LINTED: Meant to convert to 8-bit value here */
  TMPL_INT_CAST_TO_INT( mat, result, clip, s32, u8, uint8,
    0, MUTIL_UINT8_MAX, intrp_ptr );
}

/* function documented in mat_cast.h */
/* written by Vikram Chalana */
mutil_errcode mats32_cast_to_u16( const sint32_mat *mat,
  boolean clip, void *intrp_ptr, uint16_mat *result )
{
  /* LINTED: Meant to convert to 16-bit value here */
  TMPL_INT_CAST_TO_INT( mat, result, clip, s32, u16, uint16,
    0, MUTIL_UINT16_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jennifer Hodgdon */
mutil_errcode mats32_cast_to_u32( const sint32_mat *mat,
  boolean clip, void *intrp_ptr, uint32_mat *result )
{
  TMPL_INT_CAST_TO_UINT( mat, result, clip, s32, u32, uint32,
    intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode mats32_cast_to_s16( const sint32_mat *mat,
  boolean clip, void *intrp_ptr, sint16_mat *result )
{
  /* LINTED: Meant to convert to 16-bit value here */
  TMPL_INT_CAST_TO_INT( mat, result, clip, s32, s16, sint16,
    MUTIL_SINT16_MIN, MUTIL_SINT16_MAX, intrp_ptr );
}


/* function documented in mat_cast.h */
/* written by Luca Cazzanti */
mutil_errcode matdbl_to_complex( const double_mat *real_part,
    const double_mat *imag_part, void *intrp_ptr, dcomplex_mat *result )
{
  mutil_errcode errcode;
  sint32        i;

  double        num_ops = 0;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start matdbl_to_complex()" );

  errcode = matcpx_validate( result );
  if( errcode ) return errcode;

  if( !real_part && !imag_part ) {
    MUTIL_ERROR( "Real and imaginary parts are both NULL pointers" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if( ( real_part && imag_part ) &&
    !MATANY_EQUAL_DIM( real_part, imag_part ) ) {
    MUTIL_ERROR( "Input matrices size mismatch" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if( real_part ) {
    if( !MATANY_EQUAL_DIM( real_part, result ) ) {
      MUTIL_ERROR( "Real part matrix and result matrix size mismatch" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    for( i = 0; i < result->nelem; i++ )
      result->data[ i ].re = real_part->data[ i ];
  }
  else {
    for( i = 0; i < result->nelem; i++ )
      result->data[ i ].re = 0;
  }

  num_ops += result->nelem;
  if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
    MUTIL_ERROR( "User Interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  if( imag_part ) {
    if( !MATANY_EQUAL_DIM( imag_part, result ) ) {
      MUTIL_ERROR( "Imaginary part matrix and result matrix size mismatch" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    for( i = 0; i < result->nelem; i++ )
      result->data[ i ].im = imag_part->data[ i ];
  }
  else {
    for( i = 0; i < result->nelem; i++ )
      result->data[ i ].im = 0;
  }
  num_ops += result->nelem;
  if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
    MUTIL_ERROR( "User Interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE( "matdbl_to_complex() done" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_cast.h */
/* written by Luca Cazzanti */
mutil_errcode matcpx_to_double( const dcomplex_mat *in_mat, void *intrp_ptr,
    double_mat *real_part, double_mat *imag_part )
{
  mutil_errcode errcode;
  sint32        i;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start matcpx_to_double()" );

  errcode = matcpx_validate( in_mat );
  if( errcode ) return errcode;

  if( !real_part && ! imag_part ) {
    MUTIL_ERROR( "NULL pointers for both real and imaginary parts matrices" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do real part */

  if ( real_part ) {
    if( !MATANY_EQUAL_DIM( in_mat, real_part ) ) {
      MUTIL_ERROR( "Complex matrix size does not match size of real part" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    for( i = 0; i < in_mat->nelem; i++ ) {
      real_part->data[i] = in_mat->data[i].re;
    }

    if( MUTIL_INTERRUPT( 2.0 * in_mat->nelem, intrp_ptr ) ) {
      MUTIL_ERROR( "User interrupt" );
      return MUTIL_ERR_INTERRUPT;
    }

  }

  /* do imaginary part */

  if( imag_part ) {
    if( !MATANY_EQUAL_DIM( in_mat, imag_part ) ) {
      MUTIL_ERROR( "Complex matrix size does not match size of "
        "imaginary part" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    for( i = 0; i < in_mat->nelem; i++ ) {
      imag_part->data[i] = in_mat->data[i].im;
    }

    if( MUTIL_INTERRUPT( 2.0 * in_mat->nelem, intrp_ptr ) ) {
      MUTIL_ERROR( "User interrupt" );
      return MUTIL_ERR_INTERRUPT;
    }

  }

  MUTIL_TRACE( "matcpx_to_double() done" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_cast.h */
/* written by Jill Goldschneider */
mutil_errcode matcpx_to_polar( const dcomplex_mat *in_mat,
  void *intrp_ptr, double_mat *magnitude, double_mat *phase )
{
  mutil_errcode trouble;
  double        radius;
  double        angle;
  double        xdif;
  double        ydif;
  sint32        i;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Starting matcpx_to_polar()" );

  trouble = matcpx_validate( in_mat );
  if( trouble ) return trouble;

  if( !magnitude && !phase ) {
    MUTIL_ERROR( "NULL pointers for both magnitude and phase parts matrices" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do magnitude */

  if ( magnitude ) {

    trouble = matdbl_validate( magnitude );
    if( trouble ) return trouble;

    if( !MATANY_EQUAL_DIM( in_mat, magnitude ) ) {
      MUTIL_ERROR( "Input and output matrices must have equal dimensions" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    for( i = 0; i < in_mat->nelem; i++ ) {
      xdif = in_mat->data[i].re;
      ydif = in_mat->data[i].im;
      radius = sqrt( xdif * xdif + ydif * ydif );
      magnitude->data[i] = radius;
    }

    if( MUTIL_INTERRUPT( in_mat->nelem * 10.0, intrp_ptr ) ) {
      MUTIL_ERROR( "User interrupt" );
      return MUTIL_ERR_INTERRUPT;
    }

  }

  /* do phase - this is identical to pgndbl_rectangular_to_polar() */

  if ( phase ) {

    trouble = matdbl_validate( phase );
    if( trouble ) return trouble;

    if( !MATANY_EQUAL_DIM( in_mat, phase ) ) {
      MUTIL_ERROR( "Input and output matrices must have equal dimensions" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    for ( i = 0; i < in_mat->nelem; i++ ) {

      xdif = in_mat->data[i].re;
      ydif = in_mat->data[i].im;

      /* check for case of 90 or -90 degrees */
      if ( xdif == 0.0 ) {
        if (  ydif >= 0.0 ) {
          angle = 0.5 * MUTIL_PI;
        }
        else {
          angle = 1.5 * MUTIL_PI;
        }
      }

      else {
        angle = atan( ydif / xdif );

        /* arctan goes from -pi/2 to pi/2 according to K&R */
        /* get the right quadrant */
        if ( xdif < 0.0 ) {
          angle += MUTIL_PI;
        }

        /* make it between 0 and 2 pi */
        if ( angle < 0.0 ) {
          angle += 2.0 * MUTIL_PI;
        }

      }

      phase->data[i] = angle;
    }

    if( MUTIL_INTERRUPT( in_mat->nelem * 10.0, intrp_ptr ) ) {
      MUTIL_ERROR( "User interrupt" );
      return MUTIL_ERR_INTERRUPT;
    }

  }
  MUTIL_TRACE( "Done with matcpx_to_polar()" );
  return MUTIL_ERR_OK;
}
