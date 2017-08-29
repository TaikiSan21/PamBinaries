
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_tmpl.h $: $Revision: #1 $, $Date: 2008/03/21 $ */

/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_TMPL_H_
#define IN_MAT_TMPL_H_

/*
   This file contains template macros for various standard matrix
   functions, that are to be used for the function bodies.
*/

#include "mat_any.h"
#include "mat_comp.h"
#include "mat_type.h"

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_limit.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_math.h"

#include <string.h>


/** Template macro for matrix validate function.
 * Macro that expands to the body of a non-universal matrix
 * validate function, such as matdbl\_validate.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_PTR        Matrix pointer to validate (function argument).
 * @usage Body of the matdbl\_validate function:
 *     #TMPL_MAT_VALIDATE(matdbl, mat);#
 * @private
 */
#define TMPL_MAT_VALIDATE( MAT_FN_PREFIX, MAT_PTR ) \
  MUTIL_TRACE("Start " #MAT_FN_PREFIX "_validate()"); \
  \
  /* Sanity checks */ \
  \
  if( !(MAT_PTR) ) { \
    MUTIL_ERROR( "NULL matrix pointer" ); \
    return MUTIL_ERR_NULL_POINTER; \
  } \
  \
  if( (MAT_PTR)->nrow <= 0 || (MAT_PTR)->ncol <= 0 ) { \
    MUTIL_ERROR("Illegal matrix dimension value"); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  if( !(MAT_PTR)->data ) { \
    MUTIL_ERROR("NULL matrix data"); \
    return MUTIL_ERR_NULL_POINTER; \
  } \
  \
  if( (MAT_PTR)->nelem != ( (MAT_PTR)->nrow * (MAT_PTR)->ncol) ) { \
    MUTIL_ERROR("Number of elements not equal to rows * columns"); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  }\
  \
  if( (MAT_PTR)->nelem > MUTIL_SINT32_MAX ) { \
    MUTIL_ERROR("Number of elements is too large"); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  }\
  MUTIL_TRACE(#MAT_FN_PREFIX "_validate() done"); \
  return MUTIL_ERR_OK


/** Template macro for matrix wrap function.
 * Macro that expands to the body of a non-universal matrix
 * wrap function, such as matdbl\_wrap.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_PTR        Wrapper matrix pointer (function argument).
 * @param MAT_NROW       Number of rows (function argument).
 * @param MAT_NCOL       Number of columns (function argument).
 * @param MAT_DATA       Data pointer to wrap (function argument).
 * @usage Body of the matdbl\_wrap function:
 *     #TMPL_MAT_WRAP(matdbl, matrix, nrow, ncol, data);#
 * @private
 */
#define TMPL_MAT_WRAP( MAT_FN_PREFIX, MAT_PTR, MAT_NROW, MAT_NCOL, MAT_DATA ) \
  MUTIL_TRACE("Start " #MAT_FN_PREFIX "_wrap()"); \
  \
  /* avoid lint warning */ \
  (void) whatssi; \
  \
  /* Check matrix pointer; other checks done in matdbl_validate. */ \
  \
  if( !(MAT_PTR) )  { \
    MUTIL_ERROR("NULL matrix pointer"); \
    return MUTIL_ERR_NULL_POINTER; \
  }\
  \
  /* wrap by putting params into structure */ \
  \
  (MAT_PTR)->nrow  = (MAT_NROW); \
  (MAT_PTR)->ncol  = (MAT_NCOL); \
  (MAT_PTR)->nelem = (MAT_NROW) * (MAT_NCOL); \
  (MAT_PTR)->data  = (MAT_DATA); \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_wrap() done except validate" ); \
  return MAT_FN_PREFIX ## _validate(MAT_PTR)


/** Template macro for matrix malloc function.
 * Macro that expands to the body of a non-universal matrix
 * malloc function, such as matdbl\_malloc.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_PTR        Matrix pointer to allocate (function argument).
 * @param MAT_NROW       Number of rows (function argument).
 * @param MAT_NCOL       Number of columns (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matdbl\_malloc function:
 *     #TMPL_MAT_MALLOC(matdbl, matrix, nrow, ncol, double);#
 * @private
 */
#define TMPL_MAT_MALLOC( MAT_FN_PREFIX, MAT_PTR, MAT_NROW, MAT_NCOL, \
  DATA_TYPE ) \
  mutil_errcode errcode; \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_malloc()" ); \
  \
  /* Sanity checks */ \
  \
  if( !(MAT_PTR ) ) { \
    MUTIL_ERROR("NULL pointer for returning matrix"); \
    return MUTIL_ERR_NULL_POINTER; \
  } \
  if( (MAT_NROW) <= 0 || (MAT_NCOL) <= 0 ) { \
    MUTIL_ERROR("Number of rows or columns not positive"); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  }\
  \
  /* Fill matrix structure */ \
  \
  (MAT_PTR)->nrow  = (MAT_NROW); \
  (MAT_PTR)->ncol  = (MAT_NCOL); \
  (MAT_PTR)->nelem = (MAT_NROW) * (MAT_NCOL); \
  MUTIL_ASSERT( (MAT_PTR)->nelem == (MAT_PTR)->nrow * (MAT_PTR)->ncol ); \
  \
  /* Allocate space */ \
  \
  errcode = mutil_malloc( (MAT_PTR)->nelem * sizeof( DATA_TYPE ), \
    (void**) &(MAT_PTR)->data ); \
  if( errcode ) return errcode; \
  \
  MUTIL_ASSERT( (MAT_PTR)->data ); \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_malloc() done"); \
  return MUTIL_ERR_OK


/** Template macro for matrix realloc function.
 * Macro that expands to the body of a non-universal matrix
 * realloc function, such as matdbl\_realloc.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_PTR        Matrix pointer to reallocate (function argument).
 * @param MAT_NROW       Number of rows (function argument).
 * @param MAT_NCOL       Number of columns (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matdbl\_realloc function:
 *     #TMPL_MAT_REALLOC(matdbl, matrix, nrow, ncol, double);#
 * @private
 */
#define TMPL_MAT_REALLOC( MAT_FN_PREFIX, MAT_PTR, MAT_NROW, MAT_NCOL, \
  DATA_TYPE ) \
   mutil_errcode errcode; \
   sint32        old_nelem; \
   \
   MUTIL_TRACE("Start " #MAT_FN_PREFIX "_realloc()"); \
   \
   /* Sanity checks */ \
   errcode = MAT_FN_PREFIX ## _validate( MAT_PTR ); \
   if(errcode) return errcode; \
   \
   if( (MAT_NROW) <= 0 || (MAT_NCOL) <= 0 ) { \
     MUTIL_ERROR("New number of rows or columns must be positive"); \
     return MUTIL_ERR_ILLEGAL_SIZE; \
   } \
   \
   /* Fill matrix structure */ \
   \
   old_nelem        = (MAT_PTR)->nelem; \
   (MAT_PTR)->nrow  = (MAT_NROW); \
   (MAT_PTR)->ncol  = (MAT_NCOL); \
   (MAT_PTR)->nelem = (MAT_NROW) * (MAT_NCOL); \
   MUTIL_ASSERT((MAT_PTR)->nelem == (MAT_PTR)->nrow * (MAT_PTR)->ncol ); \
   \
   /* Reallocate data */ \
   \
   errcode = mutil_realloc( (void**) &((MAT_PTR)->data), \
     ((MAT_PTR)->nelem) * sizeof(DATA_TYPE), \
     old_nelem * sizeof(DATA_TYPE)); \
   if(errcode ) return errcode; \
   MUTIL_ASSERT( (MAT_PTR)->data ); \
   \
   MUTIL_TRACE( #MAT_FN_PREFIX "_realloc() done"); \
   return MUTIL_ERR_OK


/** Template macro for matrix free function.
 * Macro that expands to the body of a non-universal matrix
 * free function, such as matdbl\_free.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_PTR        Matrix pointer to free (function argument).
 * @param DATA_SIZE      Size of data type for this matrix type.
 * @usage Body of the matdbl\_free function:
 *     #TMPL_MAT_FREE(matdbl, mat, sizeof(double));#
 * @private
 */
#define TMPL_MAT_FREE( MAT_FN_PREFIX, MAT_PTR, DATA_SIZE ) \
   mutil_errcode errcode; \
   \
   MUTIL_TRACE("Start " #MAT_FN_PREFIX "_free()"); \
   \
   /* sanity checks */ \
   if( !(MAT_PTR) || !(MAT_PTR)->data ) { \
     MUTIL_WARN("Attempt to free NULL matrix"); \
     return MUTIL_ERR_OK; \
   } \
   \
   /* free data */ \
   \
   errcode = mutil_free((MAT_PTR)->data, (MAT_PTR)->nelem * DATA_SIZE ); \
   if(errcode) return errcode;  \
   (MAT_PTR)->data = NULL; \
   \
   MUTIL_TRACE( #MAT_FN_PREFIX "_free() done"); \
   return MUTIL_ERR_OK


/** Template macro for matrix assign function.
 * Macro that expands to the body of a non-universal matrix
 * assign function, such as matdbl\_assign.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to assign from (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to assign to (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_assign function:
 *     #TMPL_MAT_ASSIGN(matdbl, mat, result, double, intrp_ptr);#
 * @private
 */
#define TMPL_MAT_ASSIGN( MAT_FN_PREFIX, MAT_IN_PTR, MAT_OUT_PTR, DATA_TYPE,\
  INTRP_PTR ) \
   mutil_errcode errcode; \
   void *        outptr; \
   size_t        size_to_move; \
   sint32        tmpsize; \
   \
   MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
   \
   MUTIL_TRACE("Start " #MAT_FN_PREFIX "_assign()"); \
   \
   /* check that both matrices are valid */ \
   \
   errcode = MAT_FN_PREFIX ## _validate( MAT_IN_PTR); \
   if(errcode) return errcode; \
   \
   errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR); \
   if(errcode) return errcode; \
   \
   /* check the dimensions are the same */ \
   \
   if( !MATANY_EQUAL_DIM( MAT_IN_PTR, MAT_OUT_PTR ) ) { \
     MUTIL_ERROR("Dimensions inconsistent between result and operand"); \
     return MUTIL_ERR_ILLEGAL_SIZE; \
   }\
   \
   if( (MAT_OUT_PTR)->data == (MAT_IN_PTR)->data ) { \
     MUTIL_WARN("Copy requested to same location"); \
     return MUTIL_ERR_OK; \
   } \
   \
   tmpsize = ((MAT_IN_PTR)->nelem) * sizeof(DATA_TYPE); \
   if( tmpsize > MUTIL_SIZE_T_MAX ) { \
     MUTIL_ERROR( "Data size is too big to move" ); \
     return MUTIL_ERR_ILLEGAL_SIZE; \
   } \
   size_to_move = (size_t) tmpsize; \
   /* use memmove to do the work, so data can overlap */ \
   outptr = memmove( (void *) (MAT_OUT_PTR)->data, \
                     (void *) (MAT_IN_PTR)->data, size_to_move ); \
   \
   if( !outptr || !(MAT_OUT_PTR)->data ) { \
      MUTIL_ERROR("Error in memmove"); \
      return MUTIL_ERR_MEM_ALLOC; \
   } \
   if( MUTIL_INTERRUPT( (double) (MAT_IN_PTR)->nelem, INTRP_PTR )) { \
      MUTIL_ERROR("User interrupt"); \
      return MUTIL_ERR_INTERRUPT; \
   } \
   \
   MUTIL_TRACE( #MAT_FN_PREFIX "_assign() done"); \
   return MUTIL_ERR_OK


/** Template macro for matrix assign scalar function.
 * Macro that expands to the body of a non-universal matrix
 * scalar assign function, such as matdbl\_assign\_scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param SCALAR         Scalar value to assign (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to assign to (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_assign\_scalar function:
 *     #TMPL_MAT_ASSIGN_SCALAR(matdbl, scalar, mat, double, intrp_ptr);#
 * @private
 */
#define TMPL_MAT_ASSIGN_SCALAR( MAT_FN_PREFIX, SCALAR, MAT_OUT_PTR, DATA_TYPE,\
  INTRP_PTR ) \
  DATA_TYPE *     ret_data; \
  sint32          index; \
  sint32          nelem; \
  mutil_errcode   errcode; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE("Start " #MAT_FN_PREFIX "_assign_scalar()"); \
  \
  /* check that matrix is valid */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR); \
  if(errcode) return errcode; \
  \
  /* do copy on whole matrix as if flat */ \
  \
  ret_data = (MAT_OUT_PTR)->data; \
  nelem    = (MAT_OUT_PTR)->nelem; \
  \
  for( index = 0; index < nelem; index++ ) { \
    ret_data[index] = SCALAR; \
  } \
  \
  if( MUTIL_INTERRUPT( 2.0 * nelem, INTRP_PTR )) { \
      MUTIL_ERROR("User interrupt"); \
      return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_assign_scalar() done"); \
  return MUTIL_ERR_OK


/** Template macro for matrix assign-pad function.
 * Macro that expands to the body of a non-universal matrix
 * assign-pad function, such as matdbl\_assign\_zeropad.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to assign from (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to assign to (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_assign\_zeropad function:
 *     #TMPL_MAT_ASSIGN_PAD(matdbl, smallmat, result, intrp_ptr);#
 * @private
 */
#define TMPL_MAT_ASSIGN_PAD( MAT_FN_PREFIX, MAT_IN_PTR, MAT_OUT_PTR, \
  INTRP_PTR ) \
  mutil_errcode errcode; \
  sint32        index; \
  sint32        rowindex; \
  sint32        nrow_in; \
  sint32        ncol_in; \
  sint32        nrow_out; \
  sint32        ncol_out; \
  sint32        row; \
  sint32        col; \
  \
  double        num_ops = 0.0; \
  double        num_ops_add; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_assign_zeropad()" ); \
  \
  /* use assign_submat for assignment and validation */ \
  \
  errcode = MAT_FN_PREFIX ## _assign_submat( MAT_IN_PTR, 0, 0, \
    INTRP_PTR, MAT_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  ncol_in = (MAT_IN_PTR)->ncol; \
  ncol_out = (MAT_OUT_PTR)->ncol; \
  nrow_in = (MAT_IN_PTR)->nrow; \
  nrow_out = (MAT_OUT_PTR)->nrow; \
  \
  /* pad data */ \
  \
  num_ops_add =  4.0 * (ncol_out - ncol_in); \
  for( row = 0, rowindex = 0; \
       row < nrow_in;  \
       row++, rowindex += ncol_out ) { \
    index = rowindex + ncol_in; \
    for( col = ncol_in; col < ncol_out; col++ ) { \
      (MAT_OUT_PTR)->data[ index ] = 0; \
      index++; \
    } \
    num_ops += num_ops_add; \
    if( MUTIL_INTERRUPT( num_ops, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
    } \
  } \
  \
  num_ops_add = 4.0 * ncol_out; \
  for( row = nrow_in, rowindex = nrow_in * ncol_out; \
       row < nrow_out; \
       row++, rowindex += ncol_out ) { \
    index = rowindex; \
    for( col = 0; col < ncol_out; col++ ) { \
      (MAT_OUT_PTR)->data[ index ] = 0; \
      index++; \
    } \
    num_ops += num_ops_add; \
    if( MUTIL_INTERRUPT( num_ops, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
    } \
  } \
  MUTIL_TRACE( #MAT_FN_PREFIX "_assign_zeropad() done" ); \
  return MUTIL_ERR_OK


/** Template macro for matrix assign submatrix function.
 * Macro that expands to the body of a non-universal matrix
 * assign submatrix function, such as matdbl\_assign\_submat.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to assign from (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to assign to (function argument).
 * @param START_ROW      Row to start assigning at (function argument).
 * @param START_COL      Column to start assigning at (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_assign\_submat function:
 *     #TMPL_MAT_ASSIGN_SUBMAT(matdbl, smallmat, result, row, col, intrp_ptr);#
 * @private
 */
#define TMPL_MAT_ASSIGN_SUBMAT( MAT_FN_PREFIX, MAT_IN_PTR, MAT_OUT_PTR, \
  START_ROW, START_COL, INTRP_PTR ) \
  mutil_errcode errcode; \
  sint32        index_in; \
  sint32        rowindex_out; \
  sint32        nrow_in; \
  sint32        ncol_in; \
  sint32        nrow_out; \
  sint32        ncol_out; \
  sint32        row_indx; \
  sint32        col_indx; \
  double        num_ops = 0.0; \
  double        num_ops_add; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_assign_submat()" ); \
  \
  /* sanity checks */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_IN_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  ncol_in = (MAT_IN_PTR)->ncol; \
  ncol_out = (MAT_OUT_PTR)->ncol; \
  nrow_in = (MAT_IN_PTR)->nrow; \
  nrow_out = (MAT_OUT_PTR)->nrow; \
  \
  if((START_ROW) < 0 || \
     nrow_out < ( (START_ROW) + nrow_in ) || \
     (START_COL) < 0 || \
     ncol_out < ( (START_COL) + ncol_in )) { \
    MUTIL_ERROR( "Result, start position, and input incompatible" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  if( (MAT_IN_PTR)->data == (MAT_OUT_PTR)->data ) { \
    MUTIL_ERROR( "Input and output matrices must not share data arrays" ); \
    return MUTIL_ERR_ILLEGAL_ADDRESS; \
  } \
  \
  /* assign data */ \
  \
  num_ops_add = 8.0 * ncol_in; \
  index_in = 0; \
  rowindex_out = (START_ROW) * ncol_out + (START_COL); \
  for( row_indx = 0; row_indx < nrow_in; row_indx++ ) { \
    for( col_indx = 0; col_indx < ncol_in; col_indx++ ) { \
      (MAT_OUT_PTR)->data[ rowindex_out + col_indx ] = \
        (MAT_IN_PTR)->data[ index_in ]; \
      index_in++; \
    } \
    rowindex_out += ncol_out; \
    \
    num_ops += num_ops_add; \
    if( MUTIL_INTERRUPT( num_ops, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
    } \
  } \
  MUTIL_TRACE( #MAT_FN_PREFIX "_assign_submat() done" ); \
  return MUTIL_ERR_OK


/** Template macro for matrix extract submatrix function.
 * Macro that expands to the body of a non-universal matrix
 * extract submatrix function, such as matdbl\_extract.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to extract from (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to extract to (function argument).
 * @param START_ROW      Row to start extracting at (function argument).
 * @param START_COL      Column to start extracting at (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_extract function:
 *     #TMPL_MAT_EXTRACT(matdbl, mat, result, start_row, start_col, intrp_ptr);#
 * @private
 */
#define TMPL_MAT_EXTRACT( MAT_FN_PREFIX, MAT_IN_PTR, MAT_OUT_PTR, \
  START_ROW, START_COL, INTRP_PTR ) \
  mutil_errcode errcode; \
  sint32        index_in; \
  sint32        index_out; \
  sint32        rowindex_in; \
  sint32        ncol_in; \
  sint32        in_row; \
  sint32        in_col; \
  sint32        end_row; \
  sint32        end_col; \
  double        num_ops = 0.0; \
  double        num_ops_add; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_extract()" ); \
  \
  /* sanity checks */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_IN_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  if( (MAT_IN_PTR)->data == (MAT_OUT_PTR)->data ) { \
    MUTIL_ERROR( "Data pointers of result and source must differ" ); \
    return MUTIL_ERR_ILLEGAL_ADDRESS; \
  } \
  ncol_in = (MAT_IN_PTR)->ncol; \
  end_row = (START_ROW) + (MAT_OUT_PTR)->nrow - 1; \
  end_col = (START_COL) + (MAT_OUT_PTR)->ncol - 1; \
  \
  if( (START_ROW) < 0 || end_row >= (MAT_IN_PTR)->nrow || \
      (START_COL) < 0 || end_col >= (MAT_IN_PTR)->ncol ) { \
    MUTIL_ERROR( "Result, start position, and input incompatible" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  /* do the extraction */ \
  \
  num_ops_add = 8.0 * (end_row - (START_ROW)); \
  \
  index_out = 0; \
  rowindex_in = (START_ROW) * ncol_in; \
  for( in_row = (START_ROW); in_row <= end_row; in_row++ ) { \
    index_in = rowindex_in + START_COL; \
    for( in_col = (START_COL); in_col <= end_col; in_col++ ) { \
      (MAT_OUT_PTR)->data[ index_out ] =  (MAT_IN_PTR)->data[ index_in ]; \
      index_in++; \
      index_out++; \
    } \
    rowindex_in += ncol_in; \
    \
    num_ops += num_ops_add; \
    if( MUTIL_INTERRUPT( num_ops, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
    } \
  } \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_extract() done" ); \
  return MUTIL_ERR_OK


/** Template macro for matrix transpose function.
 * Macro that expands to the body of a non-universal matrix
 * transpose function, such as matdbl\_transpose.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to transpose (function argument).
 * @param MAT_OUT_PTR    Matrix pointer for result (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_transpose function:
 *     #TMPL_MAT_TRANSPOSE(matdbl, mat, result, intrp_ptr);#
 * @private
 */
#define TMPL_MAT_TRANSPOSE( MAT_FN_PREFIX, MAT_IN_PTR, MAT_OUT_PTR, \
  INTRP_PTR ) \
  mutil_errcode errcode; \
  sint32        tmp; \
  sint32        row; \
  sint32        col; \
  sint32        nrow; \
  sint32        ncol; \
  sint32        index_out; \
  sint32        row_index_in; \
  double        num_ops = 0.0; \
  double        num_ops_add; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_transpose()" ); \
  \
  /* sanity checks */ \
  \
  errcode = MAT_FN_PREFIX ## _validate(MAT_IN_PTR); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate(MAT_OUT_PTR); \
  if(errcode) return errcode; \
  \
  /* allow same input/output only for vectors */ \
  if( (MAT_IN_PTR)->data == (MAT_OUT_PTR)->data ) { \
    if( (MAT_IN_PTR)->nrow != 1 && (MAT_IN_PTR)->ncol != 1 ) { \
      MUTIL_ERROR( "In-place transpose only allowed for 1D matrices" ); \
      return MUTIL_ERR_ILLEGAL_SIZE; \
    } \
    \
    /* switch row and column number to transpose vector */ \
    \
    tmp                 = (MAT_OUT_PTR)->nrow; \
    (MAT_OUT_PTR)->nrow = (MAT_OUT_PTR)->ncol; \
    (MAT_OUT_PTR)->ncol = tmp; \
  } \
  else { /* result and input not same address */ \
    if((MAT_OUT_PTR)->nrow != (MAT_IN_PTR)->ncol || \
       (MAT_OUT_PTR)->ncol != (MAT_IN_PTR)->nrow ) { \
      MUTIL_ERROR( "Incompatible matrix dimensions" ); \
      return MUTIL_ERR_ILLEGAL_SIZE; \
    } \
    \
    nrow = (MAT_OUT_PTR)->nrow; \
    ncol = (MAT_OUT_PTR)->ncol; \
    \
    num_ops_add = 8.0 * ncol; \
    index_out = 0; \
    for( row = 0; row < nrow; row++ ) { \
      row_index_in = 0; \
      for( col = 0 ; col < ncol; col++ ) { \
        (MAT_OUT_PTR)->data[ index_out ] = \
          (MAT_IN_PTR)->data[ row_index_in + row ]; \
        index_out++; \
        row_index_in += nrow; \
      } \
      num_ops += num_ops_add; \
      if( MUTIL_INTERRUPT( num_ops, INTRP_PTR )) { \
        MUTIL_ERROR( "User interrupt" ); \
        return MUTIL_ERR_INTERRUPT; \
      } \
    } \
  } /* if same or different addresses */ \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_transpose() done" ); \
  return MUTIL_ERR_OK


/** Template macro for matrix range function.
 * Macro that expands to the body of a non-universal matrix
 * range function, such as matdbl\_range.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_PTR        Matrix pointer to take range of (function argument).
 * @param MIN_PTR        Pointer to put min value in (function argument).
 * @param MAX_PTR        Pointer to put max value in (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_range function:
 *     #TMPL_MAT_RANGE(matdbl, mat, minval, maxval, double, intrp_ptr);#
 * @private
 */
#define TMPL_MAT_RANGE( MAT_FN_PREFIX, MAT_PTR, MIN_PTR, MAX_PTR, \
  DATA_TYPE, INTRP_PTR ) \
  mutil_errcode      errcode; \
  const DATA_TYPE *  in_dat; \
  sint32             index; \
  sint32             nelem; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_range()" ); \
  \
  /* sanity checks */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_PTR ); \
  if(errcode) return errcode; \
  \
  if( !(MAX_PTR) || !(MIN_PTR) ) { \
    MUTIL_ERROR("NULL result pointer"); \
    return MUTIL_ERR_NULL_POINTER; \
  } \
  \
  /* do the range calculation */ \
  \
  in_dat = mat->data; \
  nelem = mat->nelem; \
  \
  *(MAX_PTR) = *(MIN_PTR) = in_dat[0]; \
  \
  for( index = 1; index < nelem; index++ ) { \
    if( in_dat[index] > *(MAX_PTR) ) *(MAX_PTR) = in_dat[index]; \
    if( in_dat[index] < *(MIN_PTR) ) *(MIN_PTR) = in_dat[index]; \
  } \
  \
  if( MUTIL_INTERRUPT( nelem * 5.0, INTRP_PTR )) { \
    MUTIL_ERROR( "User interrupt" ); \
    return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_range() done" ); \
  return MUTIL_ERR_OK



/** Template macro for matrix robust range function.
 * Macro that expands to the body of a non-universal matrix
 * robust range function, such as matdbl\_range\_robust.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_PTR        Matrix pointer to take range of (function argument).
 * @param LOWEXC         Exclude numbers below this value (function argument).
 * @param HIEXC          Exclude numbers above this value (function argument).
 * @param MIN_PTR        Pointer to put min value in (function argument).
 * @param MAX_PTR        Pointer to put max value in (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_range\_robust function:
 *     #TMPL_MAT_RANGE_ROBUST(matdbl, mat, lowexc, highexc, minval, maxval, double, intrp_ptr);#
 * @private
 */
#define TMPL_MAT_RANGE_ROBUST( MAT_FN_PREFIX, MAT_PTR, LOWEXC, HIEXC, \
  MIN_PTR, MAX_PTR, DATA_TYPE, INTRP_PTR ) \
  mutil_errcode     errcode; \
  const DATA_TYPE * in_dat; \
  sint32            index; \
  sint32            nelem; \
  boolean           valfound = FALSE; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_range_robust()" ); \
  \
  /* matrix validation */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_PTR ); \
  if(errcode) return errcode; \
  \
  if( !(MAX_PTR) || !(MIN_PTR) ) { \
    MUTIL_ERROR("NULL result pointer"); \
    return MUTIL_ERR_NULL_POINTER; \
  } \
  \
  if( HIEXC < LOWEXC ) { \
    MUTIL_ERROR( "Excluding all data, cannot compute range"); \
    return MUTIL_ERR_ILLEGAL_VALUE; \
  } \
  \
  in_dat = (MAT_PTR)->data; \
  nelem  = (MAT_PTR)->nelem; \
  \
  for( index = 0; index < nelem; index++ ) { \
    if( in_dat[index] < LOWEXC || in_dat[index] > HIEXC ) { \
      continue; \
    } \
    if( !valfound ) { \
      *(MAX_PTR) = *(MIN_PTR) = in_dat[index]; \
      valfound = TRUE; \
    } \
    else { \
      if( in_dat[index] > *(MAX_PTR) ) *(MAX_PTR) = in_dat[index]; \
      if( in_dat[index] < *(MIN_PTR) ) *(MIN_PTR) = in_dat[index]; \
    } \
  } /* for loop over pixels */ \
  \
  if( MUTIL_INTERRUPT( nelem * 8.0, INTRP_PTR )) { \
    MUTIL_ERROR( "User interrupt" ); \
    return MUTIL_ERR_INTERRUPT; \
  } \
  \
  if( !valfound ) { \
    MUTIL_ERROR( "No values found inside allowed range" ); \
    return MUTIL_ERR_ILLEGAL_VALUE; \
  } \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_range_robust() done" ); \
  return MUTIL_ERR_OK


/** Template macro for matrix translation function.
 * Macro that expands to the body of a non-universal matrix
 * translation function, such as matdbl\_translate.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Pointer to input matrix (function argument).
 * @param MAT_OUT_PTR    Pointer to translated matrix (function argument).
 * @param ROW_SHIFT      Number of rows to translate (function argument).
 * @param COL_SHIFT      Number of columns to translate (function argument).
 * @param PAD_VALUE      Value for padding the translated matrix
 *            (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_translate function:
 *     #TMPL_MAT_TRANSLATE(matdbl, mat, result, row_shift, col_shift, pad_value, intrp_ptr);#
 * @private
 */
#define TMPL_MAT_TRANSLATE( MAT_FN_PREFIX, MAT_IN_PTR, MAT_OUT_PTR, \
  ROW_SHIFT, COL_SHIFT, PAD_VALUE, INTRP_PTR ) \
  mutil_errcode errcode; \
  sint32        start_row_mat; \
  sint32        end_row_mat; \
  sint32        start_row_result; \
  sint32        start_col_mat; \
  sint32        end_col_mat; \
  sint32        start_col_result; \
  sint32        row_idx; \
  sint32        col_idx; \
  sint32        row_index_in; \
  sint32        row_index_out; \
  sint32        index_in; \
  sint32        index_out; \
  sint32        ncol; \
  double        num_ops = 0.0; \
  double        num_ops_add; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE("Start " #MAT_FN_PREFIX "_translate()"); \
  \
  /* case of no shift */ \
  \
  if( (ROW_SHIFT) == 0 && (COL_SHIFT) == 0) { \
     errcode = MAT_FN_PREFIX ## _assign(MAT_IN_PTR, INTRP_PTR, MAT_OUT_PTR); \
     if(errcode) return(errcode); \
     \
     MUTIL_TRACE(#MAT_FN_PREFIX "_translate() done"); \
     return MUTIL_ERR_OK; \
  } \
  \
  /* Matrix validation */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_IN_PTR); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate(MAT_OUT_PTR); \
  if(errcode) return errcode;  \
  \
  /* don't allow same input/output */ \
  if( (MAT_IN_PTR)->data == (MAT_OUT_PTR)->data ) { \
     MUTIL_ERROR( "In-place translation not allowed" ); \
     return(MUTIL_ERR_ILLEGAL_SIZE); \
  }\
  \
  /* Check for equal dimensions of matrices */ \
  if( !MATANY_EQUAL_DIM( MAT_IN_PTR, MAT_OUT_PTR ) ) { \
     MUTIL_ERROR("Dimensions inconsistent between result and operand"); \
     return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  ncol = (MAT_IN_PTR)->ncol; \
  \
  /* translate matrix */ \
  \
  /* default result is matrix of all pad_values */ \
  \
  errcode = MAT_FN_PREFIX ## _assign_scalar(PAD_VALUE, INTRP_PTR, \
    MAT_OUT_PTR); \
  if(errcode) return errcode; \
  \
  /* shift higher than dimension */ \
  if( (ROW_SHIFT) >= (MAT_IN_PTR)->nrow || \
    (COL_SHIFT) >= (MAT_IN_PTR)->ncol || \
    -(ROW_SHIFT) >= (MAT_IN_PTR)->nrow || \
    -(COL_SHIFT) >= (MAT_IN_PTR)->ncol   ) { \
    MUTIL_TRACE(#MAT_FN_PREFIX "_translate() done"); \
    return MUTIL_ERR_OK; \
  } \
  \
  /* determine starting and ending points for assigment between matrices */ \
  if( (ROW_SHIFT) >= 0) { \
    start_row_mat    = 0; \
    end_row_mat      = (MAT_IN_PTR)->nrow - (ROW_SHIFT); \
    start_row_result = (ROW_SHIFT); \
  } \
  else { \
    start_row_mat    = -(ROW_SHIFT); \
    end_row_mat      = (MAT_IN_PTR)->nrow; \
    start_row_result = 0; \
  } \
  \
  if( (COL_SHIFT) >= 0) { \
     start_col_mat    = 0; \
     end_col_mat      = (MAT_IN_PTR)->ncol - (COL_SHIFT); \
     start_col_result = (COL_SHIFT); \
  } \
  else { \
     start_col_mat    = -(COL_SHIFT); \
     end_col_mat      = (MAT_IN_PTR)->ncol; \
     start_col_result = 0; \
  } \
  \
  /* do translation finally */ \
  \
  num_ops_add = 8.0 * ( end_col_mat - start_col_mat + 1 ); \
  \
  /* row_idx_result = start_row_result;*/ \
  row_index_in = start_row_mat * ncol; \
  row_index_out = start_row_result * ncol; \
  for( row_idx = start_row_mat; \
       row_idx < end_row_mat; \
       row_idx++, row_index_in += ncol, row_index_out += ncol ) { \
    index_in = row_index_in + start_col_mat; \
    index_out = row_index_out + start_col_result; \
    for( col_idx = start_col_mat; \
         col_idx < end_col_mat; \
         col_idx++, index_in++, index_out++ ) { \
      \
      (MAT_OUT_PTR)->data[ index_out ] = (MAT_IN_PTR)->data[ index_in ]; \
    } \
    num_ops += num_ops_add; \
    if( MUTIL_INTERRUPT( num_ops, intrp_ptr )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
    } \
  } \
  \
  MUTIL_TRACE(#MAT_FN_PREFIX "_translate() done"); \
  return MUTIL_ERR_OK


/** Template macro for matrix flip left/right function.
 * Macro that expands to the body of a non-universal matrix
 * flip left/right function, such as matdbl\_flip\_left\_right.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to flip (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer for result (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matdbl\_flip\_left\_right function:
 *     #TMPL_MAT_FLIP_LEFT_RIGHT(matdbl, mat, intrp_ptr, result, double);#
 * @private
 */
#define TMPL_MAT_FLIP_LEFT_RIGHT( MAT_FN_PREFIX, MAT_IN_PTR, \
  INTRP_PTR, MAT_OUT_PTR, DATA_TYPE) \
  mutil_errcode  errcode; \
  DATA_TYPE     *in_data; \
  DATA_TYPE     *ret_data; \
  DATA_TYPE      tmp_num; \
  sint32         row; \
  sint32         col; \
  sint32         nrow; \
  sint32         ncol; \
  sint32         index; \
  sint32         index1; \
  sint32         index2; \
  double         num_ops = 0.0; \
  double         num_ops_add; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_flip_left_right()" ); \
  \
  /* sanity checks */ \
  \
  errcode = MAT_FN_PREFIX ## _validate(MAT_IN_PTR); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate(MAT_OUT_PTR); \
  if(errcode) return errcode; \
  \
  if( !MATANY_EQUAL_DIM( MAT_IN_PTR, MAT_OUT_PTR ) ) { \
    MUTIL_ERROR("Dimensions inconsistent between result and operand"); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  }\
  \
  nrow = (MAT_OUT_PTR)->nrow; \
  ncol = (MAT_OUT_PTR)->ncol; \
  in_data  = (MAT_IN_PTR)->data; \
  ret_data = (MAT_OUT_PTR)->data; \
  \
  index = 0; \
  num_ops_add = 8.0 * ncol; \
  for( row = 0; row < nrow; row++, index += ncol ) { \
    for( col = 0; col < ceil((double) ncol / 2.0); col++ ) { \
      index1 = index + col; \
      index2 = index + (ncol - 1 - col); \
      tmp_num = in_data[index1]; \
      ret_data[index1] = in_data[index2]; \
      ret_data[index2] = tmp_num; \
    } \
    \
    num_ops += num_ops_add; \
    if( MUTIL_INTERRUPT( num_ops, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
    } \
  } \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_flip_left_right() done" ); \
  return MUTIL_ERR_OK


/** Template macro for matrix flip up/down function.
 * Macro that expands to the body of a non-universal matrix
 * flip up/down function, such as matdbl\_flip\_up\_down.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to flip (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer for result (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matdbl\_flip\_up\_down function:
 *     #TMPL_MAT_FLIP_UP_DOWN(matdbl, mat, intrp_ptr, result, double);#
 * @private
 */
#define TMPL_MAT_FLIP_UP_DOWN( MAT_FN_PREFIX, MAT_IN_PTR,               \
  INTRP_PTR, MAT_OUT_PTR, DATA_TYPE)                                    \
  mutil_errcode  errcode;                                               \
  DATA_TYPE     *in_data;                                               \
  DATA_TYPE     *ret_data;                                              \
  DATA_TYPE      tmp_num;                                               \
  sint32         row;                                                   \
  sint32         col;                                                   \
  sint32         nrow;                                                  \
  sint32         ncol;                                                  \
  sint32         index1;                                                \
  sint32         index2;                                                \
  double         num_ops = 0.0;                                         \
  double         num_ops_add;                                           \
                                                                        \
  MUTIL_INTERRUPT_INIT( INTRP_PTR );                                    \
                                                                        \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_flip_up_down()" );             \
                                                                        \
  /* sanity checks */                                                   \
                                                                        \
  errcode = MAT_FN_PREFIX ## _validate(MAT_IN_PTR);                     \
  if(errcode) return errcode;                                           \
                                                                        \
  errcode = MAT_FN_PREFIX ## _validate(MAT_OUT_PTR);                    \
  if(errcode) return errcode;                                           \
                                                                        \
  if( !MATANY_EQUAL_DIM( MAT_IN_PTR, MAT_OUT_PTR ) ) {                  \
    MUTIL_ERROR("Dimensions inconsistent between result and operand");  \
    return MUTIL_ERR_ILLEGAL_SIZE;                                      \
  }                                                                     \
                                                                        \
  nrow = (MAT_OUT_PTR)->nrow;                                           \
  ncol = (MAT_OUT_PTR)->ncol;                                           \
  in_data  = (MAT_IN_PTR)->data;                                        \
  ret_data = (MAT_OUT_PTR)->data;                                       \
                                                                        \
  index1 = 0;                                                           \
  index2 = (nrow - 1) * ncol;                                           \
  num_ops_add = 8.0 * ncol;                                             \
  for( row = 0; row < ceil((double) nrow / 2.0);                        \
      row++, index2-=ncol ) {                                           \
    for( col = 0; col < ncol; col++, index1++  ) {                      \
      tmp_num = in_data[index1];                                        \
      ret_data[index1] = in_data[index2 + col];                         \
      ret_data[index2 + col] = tmp_num;                                 \
    }                                                                   \
                                                                        \
    num_ops += num_ops_add;                                             \
    if( MUTIL_INTERRUPT( num_ops, INTRP_PTR )) {                        \
      MUTIL_ERROR( "User interrupt" );                                  \
      return MUTIL_ERR_INTERRUPT;                                       \
    }                                                                   \
  }                                                                     \
                                                                        \
  MUTIL_TRACE( #MAT_FN_PREFIX "_flip_up_down() done" );                 \
  return MUTIL_ERR_OK


/** Template macro for comparing elements of two matrices and returning
 * the number of elements that satisfy the comparison.
 * Macro that expands to the body of a non-universal matrix
 * count of equal elements function, such as matdbl\_number\_equal.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_FN_NAME    With MAT\_FN\_PREFIX, the name of the function
 *     for this matrix type.
 * @param MAT1_PTR       Matrix pointer to first matrix (function argument).
 * @param COMPARATOR     Comparator symbol, such as ==.
 * @param MAT2_PTR       Matrix pointer to second matrix (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param OUT_VAL        Pointer to integer value to hold number of
 *     equal elements (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matdbl\_number\_equal function:
 *     #TMPL_MAT_COMPARE_MAT_TALLY(matdbl, _number_equal, mat1, ==, mat2, intrp_ptr, result, double);#
 * @private
 */
#define TMPL_MAT_COMPARE_MAT_TALLY( MAT_FN_PREFIX, MAT_FN_NAME, \
  MAT1_PTR, COMPARATOR, MAT2_PTR, \
  INTRP_PTR, OUT_VAL, DATA_TYPE ) \
  DATA_TYPE *     in_data1; \
  DATA_TYPE *     in_data2; \
  sint32          i; \
  sint32          nelem; \
  uint32          tally; \
  mutil_errcode   errcode; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  \
  /* check that input/output are valid */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT1_PTR ); \
  if ( errcode ) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT2_PTR ); \
  if ( errcode ) return errcode; \
  \
  if ( !OUT_VAL ) { \
    MUTIL_ERROR( "NULL pointer for operand or result" ); \
    return MUTIL_ERR_NULL_POINTER; \
  } \
  \
  /* check that dimensions are the same */ \
  \
  if ( !MATANY_EQUAL_DIM( MAT1_PTR, MAT2_PTR ) ) { \
    MUTIL_ERROR("Dimensions inconsistent between operands"); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  /* do comparison on whole matrix as if flat */ \
  \
  in_data1 = (MAT1_PTR)->data; \
  in_data2 = (MAT2_PTR)->data; \
  nelem    = (MAT1_PTR)->nelem; \
  \
  tally = 0; \
  \
  for( i = 0; i < nelem; i++ ) { \
    if ( ( in_data1[i] ) COMPARATOR ( in_data2[i]) ) { \
      tally++; \
    } \
  } \
  *OUT_VAL = tally; \
  \
  if ( MUTIL_INTERRUPT( 3.0 * nelem, INTRP_PTR )) { \
      MUTIL_ERROR("User interrupt"); \
      return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "Done with " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  return MUTIL_ERR_OK


/** Template macro for comparing elements of a matrix to a scalar and returning
 * the number of elements that satisfy the comparison.
 * Macro that expands to the body of a non-universal matrix
 * count of equal elements function, such as matdbl\_number\_equal\_scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_FN_NAME    With MAT\_FN\_PREFIX, the name of the function
 *     for this matrix type.
 * @param MAT_PTR        Matrix pointer to first matrix (function argument).
 * @param COMPARATOR     Comparator symbol, such as ==.
 * @param SCALAR         Scalar value (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param OUT_VAL        Pointer to integer value to hold number of
 *     equal elements (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matdbl\_compare\_scalar function:
 *     #TMPL_MAT_COMPARE_SCALAR_TALLY(matdbl, mat, scalar, intrp_ptr, result, double);#
 * @private
 */
#define TMPL_MAT_COMPARE_SCALAR_TALLY( MAT_FN_PREFIX, MAT_FN_NAME, \
   MAT_PTR, COMPARATOR, SCALAR, INTRP_PTR, \
   OUT_VAL, DATA_TYPE ) \
  DATA_TYPE *     in_data; \
  sint32          i; \
  sint32          nelem; \
  uint32          tally; \
  mutil_errcode   errcode; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  \
  /* check that input/output are valid */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_PTR ); \
  if ( errcode ) return errcode; \
  \
  if ( !OUT_VAL ) { \
    MUTIL_ERROR( "NULL pointer for operand or result" ); \
    return MUTIL_ERR_NULL_POINTER; \
  } \
  \
  /* do comparison on whole matrix as if flat */ \
  \
  in_data = (MAT_PTR)->data; \
  nelem   = (MAT_PTR)->nelem; \
  \
  tally = 0; \
  \
  for( i = 0; i < nelem; i++ ) { \
    if ( ( in_data[i] ) COMPARATOR ( SCALAR ) ) { \
      tally++; \
    } \
  } \
  *OUT_VAL = tally; \
  \
  if ( MUTIL_INTERRUPT( 3.0 * nelem, INTRP_PTR )) { \
      MUTIL_ERROR("User interrupt"); \
      return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "Done with " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  return MUTIL_ERR_OK


/** Template macro for matrix multiply function.
 * Macro that expands to the body of a non-universal matrix
 * multiply function, such as matdbl\_multiply.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT1_PTR       Matrix pointer to first matrix (function argument).
 * @param MAT2_PTR       Matrix pointer to second matrix (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to result matrix (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matdbl\_multiply function:
 *     #TMPL_MAT_MULT(matdbl, mat1, mat2, intrp_ptr, result, double);#
 * @see matuniv_multiply
 * @private
 */
#define TMPL_MAT_MULT( MAT_FN_PREFIX, MAT1_PTR, MAT2_PTR, INTRP_PTR, \
   MAT_OUT_PTR, DATA_TYPE ) \
  DATA_TYPE     tmp; \
  DATA_TYPE    *data1; \
  DATA_TYPE    *data2; \
  DATA_TYPE    *data_out; \
  sint32        row; \
  sint32        ncol1; \
  sint32        ncol2; \
  sint32        nrow_out; \
  sint32        ncol_out; \
  sint32        row_index1; \
  sint32        row_index2; \
  sint32        row_index_out; \
  sint32        col; \
  sint32        cc;  \
  mutil_errcode errcode; \
  double        num_ops = 0.0; \
  sint32        num_ops_add; \
  \
  MUTIL_INTERRUPT_INIT( intrp_ptr ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_multiply()" ); \
  \
  /* matrix validation */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR ); \
  if( errcode ) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT1_PTR ); \
  if( errcode ) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT2_PTR ); \
  if( errcode ) return errcode; \
  \
  if( ( MAT_OUT_PTR )->nrow != ( MAT1_PTR )->nrow || \
      ( MAT_OUT_PTR )->ncol != ( MAT2_PTR )->ncol || \
      ( MAT1_PTR )->ncol != ( MAT2_PTR )->nrow ) { \
    MUTIL_ERROR( "Inconsistent matrix dimensions" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  if( ( MAT1_PTR )->data == ( MAT_OUT_PTR )->data || \
      ( MAT2_PTR )->data == ( MAT_OUT_PTR )->data ) { \
    MUTIL_ERROR( "Input and output matrices must not share data arrays" ); \
    return MUTIL_ERR_ILLEGAL_ADDRESS; \
  } \
  \
  num_ops_add = ( MAT1_PTR )->ncol * ( MAT_OUT_PTR )->ncol * 8; \
  \
  ncol1 =    ( MAT1_PTR )->ncol; \
  ncol2    = ( MAT2_PTR )->ncol; \
  nrow_out = ( MAT_OUT_PTR )->nrow; \
  ncol_out = ( MAT_OUT_PTR )->ncol; \
  \
  data1 = ( MAT1_PTR )->data; \
  data2 = ( MAT2_PTR )->data; \
  data_out = ( MAT_OUT_PTR )->data; \
  \
  /* Simple, direct multiplication of matrices.  nothing fancy here. */ \
  \
  row_index1 = 0; \
  row_index_out = 0; \
  \
  for( row = 0; row < nrow_out; row++ ) { \
    for( col = 0; col < ncol_out; col++ ) { \
      tmp = 0; \
      row_index2 = 0; \
      for( cc = 0; cc < ncol1; cc++ ) { \
      \
        tmp += ( DATA_TYPE ) ( data1[ row_index1 + cc ] ) * \
               ( DATA_TYPE ) ( data2[ row_index2 + col ] ); \
        row_index2 += ncol2; \
      } \
      data_out[ row_index_out + col ] = tmp; \
      \
    } \
    num_ops += num_ops_add; \
    if( MUTIL_INTERRUPT( num_ops, intrp_ptr )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
    } \
    row_index1 += ncol1; \
    row_index_out += ncol_out; \
  } \
  \
  MUTIL_TRACE( "MAT_FN_PREFIX ## _multiply() done" ); \
  return MUTIL_ERR_OK


/** Template macro for matrix row summation function.
 * Macro that expands to the body of a non-universal matrix
 * row summation function, such as matdbl\_sum\_rows.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to input matrix (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to result matrix (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matdbl\_sum\_rows function:
 *     #TMPL_MAT_SUM_ROWS(matdbl, mat, intrp_ptr, result, double);#
 * @see matuniv_sum_rows
 * @private
 */
#define TMPL_MAT_SUM_ROWS( MAT_FN_PREFIX, MAT_IN_PTR, INTRP_PTR, \
   MAT_OUT_PTR, DATA_TYPE ) \
  DATA_TYPE       sum; \
  sint32          row; \
  sint32          col; \
  sint32          nrow; \
  sint32          ncol; \
  sint32          index; \
  mutil_errcode   errcode; \
  double          num_ops = 0.0; \
  double          num_ops_add; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_sum_rows()" ); \
  \
  /* check that matrices are valid */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_IN_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  /* check the dimensions are the correct */ \
  \
  if( (MAT_OUT_PTR)->nrow != (MAT_IN_PTR)->nrow || \
      (MAT_OUT_PTR)->ncol != 1 ) { \
    MUTIL_ERROR( "Dimensions inconsistent between result and operand" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  /* compute the row sums */ \
  \
  nrow = (MAT_IN_PTR)->nrow; \
  ncol = (MAT_IN_PTR)->ncol; \
  \
  num_ops_add = 4.0 * ncol; \
  index = 0; \
  for( row = 0; row < nrow; row++ ) { \
    sum=0; \
    for( col=0; col < ncol; col++ ) { \
      sum += (MAT_IN_PTR)->data[ index ]; \
      index++; \
    } \
    (MAT_OUT_PTR)->data[ row ] = sum; \
    \
    num_ops += num_ops_add; \
    if( MUTIL_INTERRUPT( num_ops, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
    } \
  } /* for loop over rows */ \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_sum_rows() done"); \
  return MUTIL_ERR_OK


/** Template macro for matrix col summation function.
 * Macro that expands to the body of a non-universal matrix
 * row summation function, such as matdbl\_sum\_cols.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to input matrix (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to result matrix (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matdbl\_sum\_rows function:
 *     #TMPL_MAT_SUM_COLS(matdbl, mat, intrp_ptr, result, double);#
 * @see matuniv_sum_cols
 * @private
 */
#define TMPL_MAT_SUM_COLS( MAT_FN_PREFIX, MAT_IN_PTR, INTRP_PTR, \
   MAT_OUT_PTR, DATA_TYPE ) \
  DATA_TYPE       sum; \
  sint32          row; \
  sint32          col; \
  sint32          nrow; \
  sint32          ncol; \
  sint32          row_index; \
  mutil_errcode   errcode; \
  double          num_ops = 0.0; \
  double          num_ops_add; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_sum_cols()" ); \
  \
  /* check that matrices are valid */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_IN_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  /* check the dimensions are the correct */ \
  \
  if( (MAT_OUT_PTR)->ncol != (MAT_IN_PTR)->ncol || \
      (MAT_OUT_PTR)->nrow != 1 ) { \
    MUTIL_ERROR( "Dimensions inconsistent between result and operand" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  /* compute the row sums */ \
  \
  nrow = (MAT_IN_PTR)->nrow; \
  ncol = (MAT_IN_PTR)->ncol; \
  \
  num_ops_add = 4.0 * ncol; \
  for( col = 0; col < ncol; col++ ) { \
    sum = 0; \
    for( row = 0, row_index = 0; row < nrow; row++, row_index += ncol ) { \
      sum += (MAT_IN_PTR)->data[ row_index + col ]; \
    } \
    (MAT_OUT_PTR)->data[ col ] = sum; \
    \
    num_ops += num_ops_add; \
    if( MUTIL_INTERRUPT( num_ops, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
    } \
  } /* for loop over rows */ \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_sum_cols() done" ); \
  return MUTIL_ERR_OK


/** Template macro for matrix summation function.
 * Macro that expands to the body of a non-universal matrix
 * summation function, such as matdbl\_sum.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX   Prefix for functions for this matrix type.
 * @param MAT_IN_PTR      Matrix pointer to input matrix (function argument).
 * @param INTRP_PTR       Pointer for interrupt handling (function argument).
 * @param DATA_OUT_PTR    Pointer to resulting sum (function argument).
 * @param DATA_TYPE       Data type of matrix data.
 * @usage Body of the matdbl\_sum\_rows function:
 *     #TMPL_MAT_SUM(matdbl, mat, intrp_ptr, result, double);#
 * @see matuniv_sum
 * @private
 */
#define TMPL_MAT_SUM( MAT_FN_PREFIX, MAT_IN_PTR, INTRP_PTR, \
   DATA_OUT_PTR, DATA_TYPE ) \
  DATA_TYPE       *in_data; \
  sint32          nelem; \
  sint32          i; \
  mutil_errcode   errcode; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_sum()" ); \
  \
  /* check that input matrix is valid */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_IN_PTR ); \
  if(errcode) return errcode; \
  \
  /* check output pointer is not null */ \
  if( !(DATA_OUT_PTR) ) { \
    MUTIL_ERROR( "Null result pointer" ); \
    return MUTIL_ERR_NULL_POINTER; \
  } \
  \
  in_data  = (MAT_IN_PTR)->data; \
  nelem    = (MAT_IN_PTR)->nelem; \
  *(DATA_OUT_PTR) = 0; \
  \
  for( i = 0; i < nelem; i++ ) { \
    *(DATA_OUT_PTR) += (DATA_TYPE) in_data[i]; \
  } \
  \
  if( MUTIL_INTERRUPT( 3.0 * nelem, INTRP_PTR )) { \
    MUTIL_ERROR( "User interrupt" ); \
    return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_sum() done"); \
  return MUTIL_ERR_OK


/** Template macro for cumulative sum function.
 * Macro that expands to the body of a non-universal matrix
 * validate function, such as matdbl\_validate.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Input matrix pointer (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Output matrix pointer (function argument).
 * @usage Body of the matdbl\_cumulative\_sum function:
 *     #TMPL_MAT_CUM_SUM(matdbl, mat);#
 * @private
 */
#define TMPL_MAT_CUM_SUM( MAT_FN_PREFIX, MAT_IN_PTR, INTRP_PTR, MAT_OUT_PTR ) \
  \
  mutil_errcode trouble; \
  double        num_ops = 0.0; \
  sint32        index; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE("Start " #MAT_FN_PREFIX "_cumulative_sum()"); \
  \
  /* Sanity checks done in mat*_assign call */ \
  \
  trouble = MAT_FN_PREFIX ## _assign( MAT_IN_PTR, INTRP_PTR, MAT_OUT_PTR ); \
  if ( trouble ) return trouble; \
  \
  for( index = 1; index < ( MAT_OUT_PTR )->nelem; index++ ) { \
    ( MAT_OUT_PTR )->data[index] += ( MAT_OUT_PTR )->data[index - 1]; \
  } \
  \
  num_ops += ( MAT_OUT_PTR )->nelem; \
  if ( MUTIL_INTERRUPT( num_ops, INTRP_PTR ) ) { \
    MUTIL_ERROR( "User interrupt" ); \
    return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE("Done with " #MAT_FN_PREFIX "_cumulative_sum()"); \
  return MUTIL_ERR_OK


/** Template macro for rescaling of matrices.
 * Macro expands to the body of a function that performs rescaling.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @usage #TMPL_RESCALE(double, matdbl, in_mat, min_val, max_val, intrp_ptr, result);#
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_IN_PTR     Input matrix pointer (function argument).
 * @param MIN_VAL        The minimum value in the resulting matrix.
 * @param MAX_VAL        The maximum value in the resulting matrix.
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Output matrix pointer (function argument).
 * @param DATA_TYPE      Type of the resulting matrix.
 * @private
 */
#define TMPL_MAT_RESCALE( \
  MAT_FN_PREFIX, \
  MAT_IN_PTR, \
  MIN_VAL, \
  MAX_VAL, \
  INTRP_PTR, \
  MAT_OUT_PTR, \
  DATA_TYPE ) \
  \
  DATA_TYPE     data_min; \
  DATA_TYPE     data_max; \
  double        scale; \
  mutil_errcode trouble; \
  sint32        i; \
  \
  /* Sanity checks */ \
  \
  trouble = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR ); \
  if ( trouble ) return trouble; \
  \
  /* range calls validate for input mat */ \
  trouble = MAT_FN_PREFIX ## _range( MAT_IN_PTR, INTRP_PTR, &data_min, &data_max ); \
  if ( trouble ) return trouble; \
  \
  /* Check for equal dimensions of matrices */ \
  if( !MATANY_EQUAL_DIM( MAT_IN_PTR, MAT_OUT_PTR ) ) { \
     MUTIL_ERROR("Matrices of inconsistent dimensions"); \
     return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  if ( MIN_VAL > MAX_VAL ) { \
    MUTIL_ERROR( "The maximum value of the new range must be larger than " \
      "or equal to the minimum value" ); \
    return MUTIL_ERR_ILLEGAL_VALUE; \
  } \
  \
  if ( data_min != data_max ) { \
    scale = ((double) MAX_VAL - (double) MIN_VAL) / \
      ((double) data_max - (double) data_min); \
  } \
  else { \
    scale = 0.0; \
  } \
  for ( i = 0; i < MAT_IN_PTR->nelem; i++ ) { \
    MAT_OUT_PTR->data[i] = \
      (DATA_TYPE) (((double) MAT_IN_PTR->data[i] - (double) data_min) * scale \
      + (double) MIN_VAL); \
  } \
  \
  return MUTIL_ERR_OK


/** Template macro for ALU operations between matrix and scalar.
 * Macro that expands to the body of a non-universal matrix-scalar
 * operation function, such as matu32\_bit\_shift\_left.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_FN_NAME    With MAT\_FN\_PREFIX, the name of the function for
 *     this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to input matrix (function argument).
 * @param OPERATION      The ALU operation.
 * @param SCALAR         The scalar data.
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to result matrix (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matu32\_bit\_and\_scalar function:
 *     #TMPL_ALU_MAT_OPERATION_SCALAR(matu32, _bit_and_scalar, mat, &, scalar, intrp_ptr, result, uint32);#
 * @private
 */
#define TMPL_ALU_MAT_OPERATION_SCALAR( MAT_FN_PREFIX,  MAT_FN_NAME, \
  MAT_IN_PTR, OPERATION, SCALAR, \
  INTRP_PTR, MAT_OUT_PTR, DATA_TYPE ) \
  DATA_TYPE      *ret_data; \
  DATA_TYPE      *in_data; \
  sint32          i; \
  sint32          nelem; \
  mutil_errcode   errcode; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  \
  /* check that matrices are valid */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_IN_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  /* check the dimensions are the same */ \
  \
  if( !MATANY_EQUAL_DIM( MAT_OUT_PTR, MAT_IN_PTR ) ) { \
    MUTIL_ERROR( "Dimensions inconsistent between result and operand" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  /* do operation on whole matrix as if flat */ \
  \
  ret_data = (MAT_OUT_PTR)->data; \
  in_data  = (MAT_IN_PTR)->data; \
  nelem    = (MAT_OUT_PTR)->nelem; \
  \
  for( i = 0; i < nelem; i++ ) { \
    ret_data[i] = ( in_data[i] ) OPERATION ( SCALAR ); \
  } \
  \
  if( MUTIL_INTERRUPT( 3.0 * nelem, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "Done with " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  return MUTIL_ERR_OK


/** Template macro for ALU operations between scalar and matrix.
 * Macro that expands to the body of a non-universal matrix-scalar
 * operation function, such as matu32\_divide\_scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_FN_NAME    With MAT\_FN\_PREFIX, the name of the function
 *     for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to input matrix (function argument).
 * @param OPERATION      The ALU operation.
 * @param SCALAR         The scalar data.
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to result matrix (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matu32\_divide\_scalar function:
 *     #TMPL_ALU_SCALAR_OPERATION_MAT(matu32, _divide_scalar, scalar, /, mat, intrp_ptr, result, uint32);#
 * @private
 */
#define TMPL_ALU_SCALAR_OPERATION_MAT( MAT_FN_PREFIX, MAT_FN_NAME, \
  SCALAR, OPERATION, MAT_IN_PTR, \
  INTRP_PTR, MAT_OUT_PTR, DATA_TYPE ) \
  DATA_TYPE      *ret_data; \
  DATA_TYPE      *in_data; \
  sint32          i; \
  sint32          nelem; \
  mutil_errcode   errcode; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  \
  /* check that matrices are valid */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_IN_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  /* check the dimensions are the same */ \
  \
  if( !MATANY_EQUAL_DIM( MAT_OUT_PTR, MAT_IN_PTR ) ) { \
    MUTIL_ERROR( "Dimensions inconsistent between result and operand" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  /* do operation on whole matrix as if flat */ \
  \
  ret_data = (MAT_OUT_PTR)->data; \
  in_data  = (MAT_IN_PTR)->data; \
  nelem    = (MAT_OUT_PTR)->nelem; \
  \
  for( i = 0; i < nelem; i++ ) { \
    ret_data[i] = ( SCALAR ) OPERATION ( in_data[i] ); \
  } \
  \
  if( MUTIL_INTERRUPT( 3.0 * nelem, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "Done with " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  return MUTIL_ERR_OK


/** Template macro for ALU operations between two matrices.
 * Macro that expands to the body of a non-universal matrix-matrix
 * operation function, such as matu32\_bit\_and.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_FN_NAME    With MAT\_FN\_PREFIX, the name of the function
 *     for this matrix type.
 * @param MAT1_PTR       Matrix pointer to first matrix (function argument).
 * @param OPERATION      The ALU operation.
 * @param MAT2_PTR       Matrix pointer to second matrix (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to result matrix (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matu32\_bit\_and  function:
 *     #TMPL_ALU_MAT_OPERATION_MAT(matu32, _bit_and, mat1, &, mat2, intrp_ptr, result, uint32);#
 * @private
 */
#define TMPL_ALU_MAT_OPERATION_MAT( MAT_FN_PREFIX,  MAT_FN_NAME, \
   MAT1_PTR, OPERATION, MAT2_PTR, INTRP_PTR, \
   MAT_OUT_PTR, DATA_TYPE ) \
  DATA_TYPE      *ret_data; \
  DATA_TYPE      *in_data1; \
  DATA_TYPE      *in_data2; \
  sint32          i; \
  sint32          nelem; \
  mutil_errcode   errcode; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  \
  /* check that matrices are valid */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT1_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT2_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  /* check that dimensions are the same */ \
  \
  if( !MATANY_EQUAL_DIM( MAT_OUT_PTR, MAT1_PTR ) ) { \
    MUTIL_ERROR( "Dimensions inconsistent between result and operand" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  if( !MATANY_EQUAL_DIM( MAT1_PTR, MAT2_PTR ) ) { \
    MUTIL_ERROR( "Dimensions inconsistent between operands" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  /* perform operation on whole matrix as if flat */ \
  \
  ret_data = (MAT_OUT_PTR)->data; \
  in_data1 = (MAT1_PTR)->data; \
  in_data2 = (MAT2_PTR)->data; \
  nelem    = (MAT_OUT_PTR)->nelem; \
  \
  for( i = 0; i < nelem; i++ ) { \
    ret_data[i] = ( in_data1[i] ) OPERATION ( in_data2[i] ); \
  } \
  \
  if( MUTIL_INTERRUPT( 3.0 * nelem, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "Done with " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  return MUTIL_ERR_OK


/** Template macro for ALU operations on one matrix.
 * Macro that expands to the body of a non-universal matrix
 * operation function, such as matu32\_bit\_not and  mats32\_abs..
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_FN_NAME    With MAT\_FN\_PREFIX, the name of the function
 *     for this matrix type.
 * @param OPERATION      The ALU operation.
 * @param MAT_PTR        Matrix pointer to first matrix (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to result matrix (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matu32\_bit\_and\_scalar and matu32\_abs functions:
 *     #TMPL_ALU_OPERATION_MAT(matu32, _bit_not, ~, mat, intrp_ptr, result, uint32);
 *     TMPL_ALU_OPERATION_MAT(mats32, _abs, MUTIL_ABS, mat, intrp_ptr, result, sint32);#
 * @private
 */
#define TMPL_ALU_OPERATION_MAT( MAT_FN_PREFIX, MAT_FN_NAME, \
   OPERATION, MAT_PTR, INTRP_PTR, \
   MAT_OUT_PTR, DATA_TYPE ) \
  DATA_TYPE      *ret_data; \
  DATA_TYPE      *in_data; \
  sint32          i; \
  sint32          nelem; \
  mutil_errcode   errcode; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  \
  /* check that matrices are valid */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  /* check that dimensions are the same */ \
  \
  if( !MATANY_EQUAL_DIM( MAT_OUT_PTR, MAT_PTR ) ) { \
    MUTIL_ERROR( "Dimensions inconsistent between result and operand" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  /* perform operation on whole matrix as if flat */ \
  \
  ret_data = (MAT_OUT_PTR)->data; \
  in_data  = (MAT_PTR)->data; \
  nelem    = (MAT_OUT_PTR)->nelem; \
  \
  for( i = 0; i < nelem; i++ ) { \
    ret_data[i] = OPERATION ( in_data[i] ); \
  } \
  \
  if( MUTIL_INTERRUPT( 3.0 * nelem, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "Done with " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  return MUTIL_ERR_OK


/** Template macro for various operations between two matrices.
 * Macro that expands to the body of a non-universal matrix-matrix
 * operation function, such as matu32\_max and matdbl\_min.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_FN_NAME    With MAT\_FN\_PREFIX, the name of the function
 *     for this matrix type.
 * @param OPERATION      The two valued operation.
 * @param MAT1_PTR       Matrix pointer to first matrix (function argument).
 * @param MAT2_PTR       Matrix pointer to second matrix (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to result matrix (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matu32\_min function:
 *     #TMPL_OPERATION_MAT_MAT(matu32, _min, MUTIL_MIN, mat1, mat2, intrp_ptr, result, uint32);#
 * @private
 */
#define TMPL_OPERATION_MAT_MAT( MAT_FN_PREFIX,  MAT_FN_NAME, \
   OPERATION, MAT1_PTR, MAT2_PTR, INTRP_PTR, \
   MAT_OUT_PTR, DATA_TYPE ) \
  DATA_TYPE      *ret_data; \
  DATA_TYPE      *in_data1; \
  DATA_TYPE      *in_data2; \
  sint32          i; \
  sint32          nelem; \
  mutil_errcode   errcode; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  \
  /* check that matrices are valid */ \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT1_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT2_PTR ); \
  if(errcode) return errcode; \
  \
  errcode = MAT_FN_PREFIX ## _validate( MAT_OUT_PTR ); \
  if(errcode) return errcode; \
  \
  /* check that dimensions are the same */ \
  \
  if( !MATANY_EQUAL_DIM( MAT_OUT_PTR, MAT1_PTR ) ) { \
    MUTIL_ERROR( "Dimensions inconsistent between result and operand" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  if( !MATANY_EQUAL_DIM( MAT1_PTR, MAT2_PTR ) ) { \
    MUTIL_ERROR( "Dimensions inconsistent between operands" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  /* perform operation on whole matrix as if flat */ \
  \
  ret_data = (MAT_OUT_PTR)->data; \
  in_data1 = (MAT1_PTR)->data; \
  in_data2 = (MAT2_PTR)->data; \
  nelem    = (MAT_OUT_PTR)->nelem; \
  \
  for( i = 0; i < nelem; i++ ) { \
    ret_data[i] = OPERATION ( in_data1[i], in_data2[i] ); \
  } \
  \
  if( MUTIL_INTERRUPT( 3.0 * nelem, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "Done with " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  return MUTIL_ERR_OK


/** Template macro to divide a scalar by each element of a matrix, or
 * each element of a matrix by a scalar. Macro expands to the body of a
 * non-universal matrix-scalar operation function,
 * such as matu32\_divide\_scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_FN_NAME    With MAT\_FN\_PREFIX, the name of the function
 *     for this matrix type.
 * @param MAT_IN_PTR     Matrix pointer to input matrix (function argument).
 * @param SCALAR         The scalar data.
 * @param MAT_NUM        Boolean flag indicating whether or not the matrix is
 *     the numerator
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to result matrix (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matu32\_divide\_scalar function:
 *     #TMPL_MAT_SCALAR_DIVIDE(matu32, _divide_scalar, scalar, mat, mat_numerator, intrp_ptr, result, uint32);#
 * @see matuniv_divide_elem
 * @private
 */
#define TMPL_MAT_SCALAR_DIVIDE( MAT_FN_PREFIX, MAT_FN_NAME, \
  SCALAR, MAT_IN_PTR, MAT_NUM, \
  INTRP_PTR, MAT_OUT_PTR, DATA_TYPE ) \
  sint32          nzero; \
  mutil_errcode   trouble; \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  \
  /* check that matrices are valid is done in the two calls below */ \
  \
  if ( !MAT_NUM ) { \
    trouble = MAT_FN_PREFIX ## _number_equal_scalar( MAT_IN_PTR, \
      (DATA_TYPE) 0.0, intrp_ptr, &nzero ); \
    \
    if ( trouble ) return trouble; \
    \
    if ( nzero >= 1 ) { \
      MUTIL_ERROR( "Denominator has one or more zero elements" ); \
      return MUTIL_ERR_DIVIDE_BY_ZERO; \
    } \
    \
    { \
      TMPL_ALU_SCALAR_OPERATION_MAT( MAT_FN_PREFIX, MAT_FN_NAME, \
        SCALAR, /, MAT_IN_PTR, \
        INTRP_PTR, MAT_OUT_PTR, DATA_TYPE ); \
    } \
  } \
  else { \
    if ( scalar == 0 ) { \
      MUTIL_ERROR( "Denominator has one or more zero elements" ); \
      return MUTIL_ERR_DIVIDE_BY_ZERO; \
    } \
    \
    { \
      TMPL_ALU_MAT_OPERATION_SCALAR( MAT_FN_PREFIX, MAT_FN_NAME, \
        MAT_IN_PTR, /, SCALAR, \
        INTRP_PTR, MAT_OUT_PTR, DATA_TYPE ); \
    } \
  }


/** Template macro to perform element by element matrix division.
 *  Macro expands to the body of a non-universal matrix-scalar
 * operation function, such as matu32\_divide\_elem.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_FN_NAME    With MAT\_FN\_PREFIX, the name of the function
 *     for this matrix type.
 * @param MAT1_IN_PTR    Matrix pointer to first input matrix (function argument).
 * @param MAT2_IN_PTR    Matrix pointer to second input matrix (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @param MAT_OUT_PTR    Matrix pointer to result matrix (function argument).
 * @param DATA_TYPE      Data type of matrix data.
 * @usage Body of the matu32\_divide\_scalar function:
 *     #TMPL_MAT_DIVIDE_MAT_ELEM(matu32, _divide_elem, mat1, mat2, intrp_ptr, result, uint32);#
 * @see matuniv_divide_scalar
 * @private
 */
#define TMPL_MAT_DIVIDE_MAT_ELEM( MAT_FN_PREFIX, MAT_FN_NAME, \
  MAT1_IN_PTR, MAT2_IN_PTR, INTRP_PTR, MAT_OUT_PTR, DATA_TYPE ) \
  sint32          nzero; \
  mutil_errcode   trouble; \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX #MAT_FN_NAME "()" ); \
  \
  /* check that matrices are valid is done in the two calls below */ \
  \
  trouble = MAT_FN_PREFIX ## _number_equal_scalar( MAT2_IN_PTR, \
    (DATA_TYPE) 0.0, intrp_ptr, &nzero ); \
  \
  if ( trouble ) return trouble; \
  \
  if ( nzero >= 1 ) { \
    MUTIL_ERROR( "Denominator has one or more zero elements" ); \
    return MUTIL_ERR_DIVIDE_BY_ZERO; \
  } \
  \
  { \
    TMPL_ALU_MAT_OPERATION_MAT( MAT_FN_PREFIX, MAT_FN_NAME, \
      MAT1_IN_PTR, /, MAT2_IN_PTR, \
      INTRP_PTR, MAT_OUT_PTR, DATA_TYPE ); \
  }


/** Template macro for comparing elements of a matrix to a scalar and returning
 * two row vectors containing those values which satisfy the comparison
 * and the (flattened) indices of those values.
 * Macro that expands to the body of a non-universal matrix
 * comparison function, such as matdbl\_compare\_scalar. The pointers
 * to the two output matrices can individually be NULL but the pointers
 * cannot both NULL. A NULL pointer for an output indicates that that
 * entity is not to be calculated.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @param MAT_PREFIX  Prefix for input matrix type like mats32.
 * @param MAT_IN_PTR  Matrix pointer to first matrix (function argument).
 * @param RELATION    MUTIL relational operator enum such as
 *                    MUTIL\_RELATION\_LESS\_THAN.
 * @param SCALAR      Scalar value (function argument).
 * @param INTRP_PTR   Pointer for interrupt handling (function argument).
 * @param MAT_INDEX_MATCH_PTR Pointer to a sint32 amtrix containing the indices for those
 *                    values of the input matrix which satify the comparison
 *                    relation. If the pointer is not NULL,the memory for
 *                    this matrix is allocated herein, otherwise no
 *                    memory is allocated and no values are returned.
 * @param MAT_MATCH_PTR Pointer to a matrix of the same type as the input,
 *                    containing the values of the input
 *                    matrix which satify the comparison relation. If the
 *                    pointer is not NULL, the memory for this matrix is
 *                    allocated herein, otherwise no memory is allocated
 *                    and no values are returned.
 * @limits Either pointer to the output matrices can be a NULL
 * @usage Body of the matdbl\_number\_equal\_scalar function:
 *     #TMPL_MAT_COMPARE_SCALAR(matdbl, mat, operator, scalar, intrp_ptr, result, double);#
 * @private
 */
#define TMPL_MAT_COMPARE_SCALAR( MAT_PREFIX, MAT_IN_PTR, RELATION, \
SCALAR, INTRP_PTR, MAT_INDEX_MATCH_PTR, MAT_MATCH_PTR ) \
  sint32          i; \
  uint32          count; \
  mutil_errcode   err; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  MUTIL_TRACE( "Start " #MAT_PREFIX  "_compare_scalar()" ); \
  \
  /* check that input/output are valid */ \
  \
  err = MAT_PREFIX ## _validate( MAT_IN_PTR ); \
  if ( err ) return err; \
  \
  if ( !MAT_INDEX_MATCH_PTR && !MAT_MATCH_PTR ) { \
    MUTIL_ERROR( "NULL pointer for both output matrices is not allowed" ); \
    return MUTIL_ERR_NULL_POINTER; \
  } \
  \
  /* allocate memory for the output */ \
  if ( MAT_INDEX_MATCH_PTR ) { \
    err = mats32_malloc( MAT_INDEX_MATCH_PTR, 1, 1 ); \
    if ( err ) return err; \
  } \
  \
  if ( MAT_MATCH_PTR ) { \
    err = MAT_PREFIX ## _malloc( MAT_MATCH_PTR, 1, 1 ); \
    if ( err ) { \
      if ( MAT_INDEX_MATCH_PTR ) MUTIL_FREE_WARN( mats32, \
        MAT_INDEX_MATCH_PTR ); \
      return err; \
    } \
  } \
  \
  /* do comparison on whole matrix as if flat */ \
  \
  count = 0; \
  \
  switch( RELATION ){ \
    case MUTIL_RELATION_LESS_THAN: \
      PERFORM_MAT_COMPARE_SCALAR( MAT_PREFIX, MAT_IN_PTR, <, SCALAR, \
        MAT_INDEX_MATCH_PTR, MAT_MATCH_PTR, count ) \
      break; \
    case MUTIL_RELATION_LESS_THAN_OR_EQUAL: \
      PERFORM_MAT_COMPARE_SCALAR( MAT_PREFIX, MAT_IN_PTR, <=, SCALAR, \
       MAT_INDEX_MATCH_PTR, MAT_MATCH_PTR, count ) \
      break; \
    case MUTIL_RELATION_EQUAL: \
      PERFORM_MAT_COMPARE_SCALAR( MAT_PREFIX, MAT_IN_PTR, ==, SCALAR, \
       MAT_INDEX_MATCH_PTR, MAT_MATCH_PTR, count ) \
      break; \
    case MUTIL_RELATION_NOT_EQUAL: \
      PERFORM_MAT_COMPARE_SCALAR( MAT_PREFIX, MAT_IN_PTR, !=, SCALAR, \
       MAT_INDEX_MATCH_PTR, MAT_MATCH_PTR, count ) \
      break; \
    case MUTIL_RELATION_GREATER_THAN: \
      PERFORM_MAT_COMPARE_SCALAR( MAT_PREFIX, MAT_IN_PTR, >, SCALAR, \
       MAT_INDEX_MATCH_PTR, MAT_MATCH_PTR, count ) \
      break; \
    case MUTIL_RELATION_GREATER_THAN_OR_EQUAL: \
      PERFORM_MAT_COMPARE_SCALAR( MAT_PREFIX, MAT_IN_PTR, >=, SCALAR, \
       MAT_INDEX_MATCH_PTR, MAT_MATCH_PTR, count ) \
      break; \
    default: \
      MUTIL_ERROR( "Relation operator is unsupported" ); \
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED; \
  } \
  \
  if ( count == 0 ){ \
    if ( MAT_INDEX_MATCH_PTR ){ \
       MUTIL_FREE_WARN( mats32, MAT_INDEX_MATCH_PTR ); \
       MAT_INDEX_MATCH_PTR->nelem = (sint32) 0; \
       MAT_INDEX_MATCH_PTR->nrow = (sint32) 0; \
       MAT_INDEX_MATCH_PTR->ncol = (sint32) 0; \
    } \
    if ( MAT_MATCH_PTR ){ \
       MUTIL_FREE_WARN( MAT_PREFIX, MAT_MATCH_PTR ); \
       MAT_MATCH_PTR->nelem = (sint32) 0; \
       MAT_MATCH_PTR->nrow = (sint32) 0; \
       MAT_MATCH_PTR->ncol = (sint32) 0; \
    } \
  }\
  \
  if ( MUTIL_INTERRUPT( 3.0 * count, INTRP_PTR )) { \
      MUTIL_ERROR("User interrupt"); \
      return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( "Done with " #MAT_PREFIX "_compare_scalar()" ); \
  return MUTIL_ERR_OK

/** Macro used exclusively in the TMPL\_MAT\_COMPARE\_SCALAR
 * template macro.
 * Macro which performs a relational comparison of each element
 * of a matrix to a scalar value and stores either (or both)
 * the value which satisfied the relationship and the
 * the flattened index of that value in the input matrix.
 *
 * @param MAT_PREFIX  Prefix for input matrix type like mats32.
 * @param MAT_IN_PTR  Matrix pointer to first matrix (function argument).
 * @param RELATION    Relational operator such as >= or !=.
 * @param SCALAR      Scalar value (function argument).
 * @param MAT_INDEX_MATCH_PTR Pointer to a sint32 amtrix containing the indices for those
 *                    values of the input matrix which satify the comparison
 *                    relation. If the pointer is not NULL,the memory for
 *                    this matrix is allocated herein, otherwise no
 *                    memory is allocated and no values are returned.
 * @param MAT_MATCH_PTR Pointer to a matrix of the same type as the input,
 *                    containing the values of the input
 *                    matrix which satify the comparison relation. If the
 *                    pointer is not NULL, the memory for this matrix is
 *                    allocated herein, otherwise no memory is allocated
 *                    and no values are returned.
 * @param COUNT       Variable to hold the current count of scalar matches.
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tmpl.h
 * @source mat\_tmpl.h
 * @library matrix
 * @usage Body of the TMPL\_MAT\_COMPARE\_SCALAR macro:
 *     #PERFORM_MAT_COMPARE_SCALAR( MAT_PREFIX, MAT_IN_PTR, <, SCALAR, MAT_INDEX_MATCH_PTR, MAT_MATCH_PTR, count )#
 * @private
 */
#define PERFORM_MAT_COMPARE_SCALAR( MAT_PREFIX, MAT_IN_PTR, RELATION, \
 SCALAR, MAT_INDEX_MATCH_PTR, MAT_MATCH_PTR, COUNT ) \
 for( i = 0; i < (MAT_IN_PTR)->nelem; i++ ) { \
 \
  if ( ( (MAT_IN_PTR)->data[i] ) RELATION ( SCALAR ) ) { \
      if ( MAT_MATCH_PTR ) { \
          \
          /* grow storage to accommodate new match */ \
          \
          err = MAT_PREFIX ## _realloc( MAT_MATCH_PTR, 1, COUNT + 1 ); \
          if ( err ) { \
            if ( MAT_INDEX_MATCH_PTR ) MUTIL_FREE_WARN( mats32, \
              MAT_INDEX_MATCH_PTR ); \
            if ( MAT_MATCH_PTR ) MUTIL_FREE_WARN( MAT_PREFIX, \
              MAT_MATCH_PTR ); \
            return err; \
          } \
          \
          (MAT_MATCH_PTR)->data[ COUNT ] = (MAT_IN_PTR)->data[ i ]; \
      } \
      \
      if ( MAT_INDEX_MATCH_PTR ) { \
          \
          /* grow storage to accommodate new match */ \
          \
          err = mats32_realloc( MAT_INDEX_MATCH_PTR, 1, COUNT + 1); \
          if ( err ) { \
            if ( MAT_INDEX_MATCH_PTR ) MUTIL_FREE_WARN( mats32, \
              MAT_INDEX_MATCH_PTR ); \
            if ( MAT_MATCH_PTR ) MUTIL_FREE_WARN( MAT_PREFIX, \
              MAT_MATCH_PTR ); \
            return err; \
          } \
          \
          (MAT_INDEX_MATCH_PTR)->data[ COUNT ] = (sint32) i; \
      } \
      \
      COUNT++; \
  } \
 }

#endif /* IN_MAT_TMPL_H_*/
