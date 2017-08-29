
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_usca.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_USCA_H
#define IN_MAT_USCA_H

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "mat_any.h"

/* This file contains functions and macros for universal matrix
   allocation and initialization operations.
 */

#ifdef __cplusplus
extern "C" {
#endif


/*
********************************
Macros for universal scalars
********************************
*/


/** Initialize a non-complex universal scalar.
 * Put a number and a type into a universal scalar structure.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_usca.h
 * @source mat\_usca.h
 * @library
 * @usage #SCAUNIV_INIT(uscal, MUTIL_DOUBLE, 3.5);#
 * @param   SCAL   Universal scalar to initialize.
 * @param   TYPE   Data type to use for universal scalar -- any type except
 *    MUTIL\_DCOMPLEX from \Ref{_mutil_data_type} will work, but the validity
 *    of the value is not checked.
 * @param   NUM  Number to put into universal scalar -- will be cast to
 *    correct type, but the number is not checked to be in range.
 * @see SCAUNIV_CAST
 * @see SCAUNIV_EQUAL
 * @see _mutil_data_type
 * @see _univ_scalar
 */
#define SCAUNIV_INIT( SCAL, TYPE, NUM ) \
  switch( TYPE ) { \
    case MUTIL_UINT8: \
      (SCAL).type   = (TYPE); \
      (SCAL).num.u8 = (uint8) (NUM); \
      break; \
    case MUTIL_SINT8: \
      (SCAL).type   = (TYPE); \
      (SCAL).num.s8 = (sint8) (NUM); \
      break; \
    case MUTIL_UINT16: \
      (SCAL).type    = (TYPE); \
      (SCAL).num.u16 = (uint16) (NUM); \
      break; \
    case MUTIL_SINT16: \
      (SCAL).type    = (TYPE); \
      (SCAL).num.s16 = (sint16) (NUM); \
      break; \
    case MUTIL_UINT32: \
      (SCAL).type    = (TYPE); \
      (SCAL).num.u32 = (uint32) (NUM); \
      break; \
    case MUTIL_SINT32: \
      (SCAL).type    = (TYPE); \
      (SCAL).num.s32 = (sint32) (NUM); \
      break; \
    case MUTIL_FLOAT: \
      (SCAL).type    = (TYPE); \
      (SCAL).num.flt = (float) (NUM); \
      break; \
    case MUTIL_DOUBLE: \
    default: \
      (SCAL).type    = (TYPE); \
      (SCAL).num.dbl = (double) (NUM); \
      break; \
  }


/** Extract a number from a non-complex universal scalar.
 * Get the number from a non-complex universal scalar and cast it to a
 * particular data type, without checking the number to see if it is in
 * range for the target data type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_usca.h
 * @source mat\_usca.h
 * @library
 * @usage #SCAUNIV_CAST(uscal, sint32);#
 * @param   SCAL   Universal scalar to extract number from -- assumed not
 *    to be of type dcomplex, but this is not checked.
 * @param   TYPE   Data type to cast to.
 * @see SCAUNIV_INIT
 * @see SCAUNIV_EQUAL
 * @see Scalar Data Types
 * @see _univ_scalar
 */
#define SCAUNIV_CAST( SCAL, TYPE ) \
   (TYPE) ( ( (SCAL).type == MUTIL_UINT8 )  ? \
     (uint8)( (SCAL).num.u8 ) : \
   ( (SCAL).type == MUTIL_SINT8 )  ? \
     (sint8)( (SCAL).num.s8 ) : \
   ( (SCAL).type == MUTIL_UINT16 ) ? \
     (uint16)( (SCAL).num.u16 ) : \
   ( (SCAL).type == MUTIL_SINT16 ) ? \
     (sint16)( (SCAL).num.s16 ) : \
   ( (SCAL).type == MUTIL_UINT32 ) ? \
     (uint32)( (SCAL).num.u32 ) : \
   ( (SCAL).type == MUTIL_SINT32 ) ? \
     (sint32)( (SCAL).num.s32 ) : \
   ( (SCAL).type == MUTIL_DOUBLE ) ? \
     (double)( (SCAL).num.dbl ) : \
     (float)( (SCAL).num.flt ))


/** Compares two non-complex universal scalars.
 * The macro returns TRUE if the two scalars have the same type and the same
 * value. Otherwise it returns FALSE.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_usca.h
 * @source mat\_usca.h
 * @library
 * @usage #SCAUNIV_EQUAL(uscal1, uscal2);#
 * @param   SCA1   First universal scalar for comparison -- assumed not
 *    to be of type dcomplex, but this is not checked.
 * @param   SCA2   Second universal scalar for comparison -- assumed not
 *    to be of type dcomplex, but this is not checked.
 * @see SCAUNIV_INIT
 * @see SCAUNIV_CAST
 * @see Scalar Data Types
 * @see _univ_scalar
 */
#define SCAUNIV_EQUAL( SCA1, SCA2 ) \
   ( ( (SCA1).type == (SCA2).type )  && \
     ( ((SCA1).type == MUTIL_UINT8  && (SCA1).num.u8  == (SCA2).num.u8 )  || \
     ((SCA1).type == MUTIL_UINT16 && (SCA1).num.u16 == (SCA2).num.u16 ) || \
     ((SCA1).type == MUTIL_UINT32 && (SCA1).num.u32 == (SCA2).num.u32 ) || \
     ((SCA1).type == MUTIL_SINT8  && (SCA1).num.s8  == (SCA2).num.s8 )  || \
     ((SCA1).type == MUTIL_SINT16 && (SCA1).num.s16 == (SCA2).num.s16 ) || \
     ((SCA1).type == MUTIL_SINT32 && (SCA1).num.s32 == (SCA2).num.s32 ) || \
     ((SCA1).type == MUTIL_FLOAT  && (SCA1).num.flt == (SCA2).num.flt ) || \
     ((SCA1).type == MUTIL_DOUBLE && (SCA1).num.dbl == (SCA2).num.dbl ) \
     ) )


#ifdef __cplusplus
}
#endif

#endif /* ifdef IN_MAT_USCA_H */
