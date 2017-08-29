
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_iter.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_ITER_H_
#define IN_MAT_ITER_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "ut_err.h"
#include "mat_type.h"


/* This file contains struct and function declarations for iterating
 * over matrices and matrix sets  */

#ifdef __cplusplus
extern "C" {
#endif


/** Enum that describes the dimension of a mutil\_iterator.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_iter.h
 * @library matrix
 *
 * @same #typedef enum _mutil_iterator_dim mutil_iterator_dim;#
 *
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
enum _mutil_iterator_dim{

    /** column dimension */
    ITER_COL,

    /** row dimension */
    ITER_ROW,

    /** row dimension */
    ITER_SLICE,

    /** other dimension in 2 dimensions */
    ITER_OTHER,

    /** other dimension across slices */
    ITER_OTHER_SLICE
};

/* See above for documentation on this enum. */
typedef enum _mutil_iterator_dim mutil_iterator_dim;


/** Struct for iterating through mutils matrices and matrix sets.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_iter.h
 * @library matrix
 *
 * @same #typedef struct _mutil_iterator mutil_iterator;#
 *
 * @see _mutil_iterator_dim
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
struct _mutil_iterator{

  /** offsets.mat.s32mat.data contains the offsets to iterate through;
   *  offsets.mat.s32mat.nrow is equal to the max\_length of the iterator */
  univ_mat     offsets;

  /** vector univ\_mat used as temporary workspace */
  univ_mat    workspace;

  /** pointer to data to iterate */
  void *  base_ptr;

  /** pointer to address of the data pointer of offsets */
  sint32 * offsets_ptr;

  /** length over which to iterate */
  sint32 length;

  /** maximum length of the iterator which equal to
   * the number of data elements in the univ\_mat offsets */
  sint32 max_length;

  /** number of columns of matrix which initialized the iterator object */
  sint32 ncol;

  /** number of rows of matrix which initialized the iterator object */
  sint32 nrow;

  /** number of slices of mat\_set which initialized the iterator object */
  sint32 nslice;

  /** mutil\_data\_type of matrix or matrix set that initialized
   * the iterator object */
  mutil_data_type type;

  /** dimension over which the iterator iterates */
  mutil_iterator_dim dimension;

};

/* See above for documentation on this struct. */
typedef struct _mutil_iterator mutil_iterator;


/* ********************* */
/* malloc mutil_iterator */
/* ********************* */


/** Allocate an iterator for a univ\_mat.
 * This function constructs an iterator for stepping through
 * a mutils matrix in the column, row or a dimension of the
 * user's design (for example, a diagonal of the matrix).  The function
 * allocates
 * space of size (length * sint32) for the offsets univ\_mat and
 * sets the base\_pointer to the address of the data pointer of the
 * input univ\_mat to the function.
 *
 * To iterate over rows (i.e. from 0 to nrow by one row) and
 * start at data[j] where j < ncol of the matrix, you would allocate a
 * mutil\_iterator of type ITER\_ROW.
 * The offset data's ith entry is i * ncol.
 * Before beginning the iteration over rows of univ\_mat at data[j],
 * set the base pointer using
 * mutil\_iterator\_assign\_base with an offset of j.
 *
 * To iterate over columns (i.e. from 0 to ncol by one col) and
 * start at row j where j < nrow of the matrix, allocate a
 * mutil\_iterator of type ITER\_COL.  The offset data's ith entry is i.
 * Before beginning the iteration over columns of univ\_mat at row j,
 * set the base pointer using
 * mutil\_iterator\_assign\_base with an offset of j*ncol.
 *
 * If a dimension of ITER\_COL or ITER\_ROW is used
 * the offsets will be automatically computed.
 * If a dimension of ITER\_OTHER is used, the user must
 * set the offsets individually using MUTIL\_ITERATOR\_ASSIGN\_OFFSET.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @source mat\_iter.c
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = mutil_iterator_malloc(&iterator, ITER_COL, &umat, length);#
 *
 * @param  iterator     Pointer to iterator that is to be initialized
 *
 * @param  dim          Dimension to iterate over, which can be
 *                      ITER\_COL, ITER\_ROW or ITER\_OTHER.
 *
 * @param  umat         Pointer to matrix to iterate through.
 *
 * @param  length       Number of integer elements to iterate over.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
MUTIL_LIBEXPORT mutil_errcode mutil_iterator_malloc(
  mutil_iterator     *iterator,
  mutil_iterator_dim  dim,
  univ_mat           *umat,
  sint32              length );


/** Allocate an iterator for a mat\_set.
 * This function constructs an iterator for stepping through a mat\_set
 * in the column, row, slice dimension or a dimension of the
 * user's design (for example, a diagonal through the matrix set).
 * The function allocates
 * space for the offsets univ\_mat and sets the base\_pointer to the
 * the address of the data pointer of mats[0] of the input mat\_set.
 * All matrices of the mat\_set must be of the same size and type.
 *
 * To iterate over slices (i.e. from 0 to nslice by one slice) and start
 * at data[j] of mats[0] where j < nelem of mats[0] of the input mat\_set,
 * allocate a mutil\_iterator of type ITER\_SLICE.
 * The offset data's ith entry is the difference between the address
 * of the data pointer of mats[i] and the address of the data pointer of
 * mats[0].  Before beginning the
 * iteration over slices of univ\_mat at data[j] of mats[0],
 * set the base pointer using mutil\_iterator\_assign\_base\_for\_iter\_slice
 * with an offset of j.
 *
 * See the documentation for \Ref{mutil_iterator_malloc} for examples of
 * mutil\_iterators of type ITER\_COL and ITER\_ROW.
 *
 * If a dimension of ITER\_COL, ITER\_ROW, or ITER\_SLICE is used
 * the offsets will be automatically computed.
 * If a dimension of ITER\_OTHER or ITER\_OTHER\_SLICE is used, the user must
 * set the offsets individually using MUTIL\_ITERATOR\_ASSIGN\_OFFSET.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @source mat\_iter.c
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = mutil_iterator_malloc_with_matset(&iterator, ITER_COL, &matset, length);#
 *
 * @param  iterator     Pointer to iterator that is to be initialized
 *
 * @param  dim          Dimension to iterate over.
 *
 * @param  matset       Pointer to matrix set to iterate through.
 *                      All matrices of matset must be the same size and type.
 *
 * @param  length       Number of integer elements to iterate over.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
MUTIL_LIBEXPORT mutil_errcode mutil_iterator_malloc_with_matset(
  mutil_iterator     *iterator,
  mutil_iterator_dim  dim,
  const mat_set      *matset,
  sint32              length );


/** Allocate a temporary iterator.
 * The function allocates workspace field of the iterator with
 * univ\_mat of dimensions 1 row and n columns
 * where n is the value of the input parameter length.  The type of
 * the workspace univ\_mat is the input parameter mutil\_type.
 *
 * The function also allocates
 * space of size length * sint32 for the offsets univ\_mat and
 * sets the iterator field base\_pointer to
 * the address of the data pointer of the workspace univ\_mat.
 * The iterator is set to type ITER\_COL hence
 * the ith data entry in offsets is i.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @source mat\_iter.c
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = mutil_iterator_construct_temp(&iterator, length, MUTIL_UINT8);#
 *
 * @param  iterator     Pointer to iterator that is to be initialized
 *
 * @param  length       Number of integer elements to iterate over.
 *
 * @param  mutil_type   Type of the temporary iterator such as MUTIL\_SINT32.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
MUTIL_LIBEXPORT mutil_errcode mutil_iterator_malloc_temp(
  mutil_iterator *iterator, sint32 length, mutil_data_type mutil_type );


/* ******************* */
/* free mutil_iterator */
/* ******************* */

/** Free memory used by a matuniv\_iterator.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @source mat\_iter.c
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = mutil_iterator_free(&iterator);#
 *
 * @param  iterator  Pointer to iterator from which to free memory.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
MUTIL_LIBEXPORT mutil_errcode mutil_iterator_free( mutil_iterator *iterator );


/* ********************************************************* */
/* functions and macros that assign fields to mutil_iterator */
/* ********************************************************* */


/** Set the iterator base\_ptr.  The iterator base\_ptr is
 * set to the sum of the address of the data field of the univ\_mat
 * parameter mat\_base and the value of the argument offset.  The umat
 * argument must be the same size and type of the univ\_mat that
 * initialized the iterator.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @source mat\_iter.c
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = mutil_iterator_assign_base(&iterator, &umat, offset);#
 *
 * @param  iterator    Pointer to iterator whose
 *                     field base\_ptr is to be changed
 *
 * @param  umat        Pointer to univ\_mat that determines the
 *                     iterator base\_ptr.  Must be the same size and type
 *                     of univ\_mat that initialized the iterator.
 *
 * @param  offset      Offset to add to the data field of mat\_base to
 *                     assign to the iterator base\_ptr.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
MUTIL_LIBEXPORT mutil_errcode mutil_iterator_assign_base(
  mutil_iterator *iterator,
  const univ_mat *umat,
  sint32          offset );


/** Set the iterator base\_ptr for iterators of
 * dimension ITER\_SLICE.  The iterator base\_ptr is
 * set to the sum of the address of the data field of
 * matset\_base->mats[0] and the value of the argument offset.
 * It is better to use this function rather than
 * mutil\_iterator\_assign\_base since error checking is tailored
 * for iterators of type ITER\_SLICE.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @source mat\_iter.c
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = mutil_iterator_assign_base_for_iter_slice(&iterator, &matset, offset);#
 *
 * @param  iterator    Pointer to iterator whose
 *                     field base\_ptr is to be changed
 *
 * @param  matset      Pointer to mat\_set that determines the
 *                     iterator base\_ptr.  Must have same dimensions as
 *                     mat\_set that initialized the iterator.
 *
 * @param  offset      Offset to add to the data field of matset\_base->mats[0]
 *                     to assign to the iterator base\_ptr.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
MUTIL_LIBEXPORT mutil_errcode mutil_iterator_assign_base_for_iter_slice(
  mutil_iterator *iterator, const mat_set *matset, sint32 offset );


/** Set the length of the iterator.  The length must
 * always be less than or equal to the max\_length of the iterator.
 *
 * If the iterator is of type ITER\_COL then max\_length is the
 * ncol field of the matrix that initialized the iterator.
 * If the iterator is of type ITER\_ROW then max\_length is the
 * nrow field of the matrix that initialized the iterator.
 * If the iterator is of type ITER\_SLICE then max\_length is the
 * nelem field of the mat\_set that initialized the iterator.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = mutil_iterator_set_length(&iterator, length);#
 *
 * @param  iterator     Pointer to iterator.
 * @param  length       Value to assign to the iterator length.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
MUTIL_LIBEXPORT mutil_errcode mutil_iterator_set_length(
  mutil_iterator *iterator, sint32 length );


/** Reassign the base of an iterator given the index.
 * No type checking nor array bounds checking is done.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @library matrix
 * @usage #MUTIL_ITERATOR_MACRO_ASSIGN_BASE(&iterator, &umat, offset, err);#
 *
 * @param  ITERATOR     Pointer to iterator.
 * @param  MAT_BASE     Pointer to a univ\_mat whose data field address will
 *                      be used as the base.
 * @param  OFFSET       Offset from the address of the data field of the
 *                      MAT\_BASE for which to set the base of the iterator.
 * @param  ERR          Variable of type mutil\_errcode to hold possible error.
 *                      ERR is initially assigned MUTIL\_ERR\_OK.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
#define MUTIL_ITERATOR_MACRO_ASSIGN_BASE( ITERATOR, MAT_BASE, OFFSET, ERR ) \
  { \
  ERR = MUTIL_ERR_OK; \
  switch ( ( ITERATOR )->type ) { \
    case MUTIL_UINT8: \
      ( ITERATOR )->base_ptr =  \
        ( void * ) ( ( uint8 * ) ( ( MAT_BASE )->mat.u8mat.data + \
          OFFSET ) ); \
      break; \
      \
     case MUTIL_SINT8: \
      ( ITERATOR )->base_ptr =  \
        ( void * ) ( ( sint8 * ) ( ( MAT_BASE )->mat.s8mat.data + \
          OFFSET ) ); \
      break; \
      \
    case MUTIL_UINT16: \
      ( ITERATOR )->base_ptr =  \
        ( void * ) ( ( uint16 * ) ( ( MAT_BASE )->mat.u16mat.data + \
          OFFSET ) ); \
      break; \
      \
     case MUTIL_SINT16: \
      ( ITERATOR )->base_ptr =  \
        ( void * ) ( ( sint16 * ) ( ( MAT_BASE )->mat.s16mat.data + \
          OFFSET ) ); \
      break; \
      \
     case MUTIL_UINT32: \
      ( ITERATOR )->base_ptr =  \
        ( void * ) ( ( uint32 * ) ( ( MAT_BASE )->mat.u32mat.data + \
          OFFSET ) ); \
      break; \
      \
     case MUTIL_SINT32: \
      ( ITERATOR )->base_ptr =  \
        ( void * ) ( ( sint32 * ) ( ( MAT_BASE )->mat.s32mat.data + \
          OFFSET ) ); \
      break; \
      \
     case MUTIL_FLOAT: \
      ( ITERATOR )->base_ptr =  \
        ( void * ) ( ( float * ) ( ( MAT_BASE )->mat.fltmat.data + \
          OFFSET ) ); \
      break; \
      \
     case MUTIL_DOUBLE: \
      ( ITERATOR )->base_ptr =  \
        ( void * ) ( ( double * ) ( ( MAT_BASE )->mat.dblmat.data + \
          OFFSET ) ); \
      break; \
      \
    /* other data types (complex) are not supported right now */ \
    default: \
      MUTIL_ERROR( "This matrix type is currently unsupported" ); \
      ERR = MUTIL_ERR_ILLEGAL_TYPE; \
  } \
}


/** Reassign the base of an iterator given the index.
 * No type checking nor array bounds checking is done.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @library matrix
 * @usage #MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST(double, dbl, &iterator, &umat, offset);#
 *
 * @param  ITERATOR     Pointer to iterator.
 * @param  MAT_BASE     Pointer to a univ\_mat whose data field address will
 *                      be used as the base.
 * @param  OFFSET       Offset from the address of the data field of the
 *                      MAT\_BASE for which to set the base of the iterator.
 * @param  ERR          Variable of type mutil\_errcode to hold possible error.
 *                      ERR is initially assigned MUTIL\_ERR\_OK.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
#define MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST( TYPE, PREFIX, ITERATOR, \
  MAT_BASE, OFFSET ) \
  ( ITERATOR )->base_ptr =  \
    ( void * ) ( ( TYPE * ) ( ( MAT_BASE )->mat.## PREFIX ##mat.data + \
    OFFSET ) );


/** Assign to an element in the iterator given the index.
 * No type checking nor array bounds checking is done.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @library matrix
 * @usage #ITERATOR_ASSIGN_DATA(float, &iterator, index, value);#
 *
 * @param  TYPE         Type of the data pointer to by the iterator, such as
 *                      sint32 or double.
 * @param  ITERATOR     Pointer to iterator.
 * @param  INDEX        Index of the data to assign.
 * @param  VALUE        Value assigned.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
#define MUTIL_ITERATOR_ASSIGN_DATA( TYPE, ITERATOR, INDEX, VALUE )    \
  *( ( TYPE * ) ( ITERATOR )->base_ptr +                              \
          MUTIL_ITERATOR_GET_OFFSET( ITERATOR, INDEX ) ) = VALUE


/** Assign to the data of the offset field in the iterator, given an index.
 * No type checking nor array bounds checking is done.  The value assigned
 * should be of type sint32.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @library matrix
 * @usage #MUTIL_ITERATOR_ASSIGN_OFFSET(&iterator, index, value);#
 *
 * @param  ITERATOR     Pointer to iterator.
 * @param  INDEX        Index of the offset data to assign.
 * @param  VALUE        Value assigned which should be of type sint32.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
#define MUTIL_ITERATOR_ASSIGN_OFFSET( ITERATOR, INDEX, VALUE )   \
  ( ITERATOR )->offsets.mat.s32mat.data[ INDEX ] = VALUE;


/** Return the length of the iterator.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @library matrix
 * @returns The length (sint32) of the iterator.
 * @usage #length = MUTIL_ITERATOR_GET_LENGTH(&iterator);#
 *
 * @param  ITERATOR    Pointer to iterator
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
#define MUTIL_ITERATOR_GET_LENGTH( ITERATOR ) ( ITERATOR )->length


/** Return an offset of the iterator.
 * If the argument INDEX is set to i, the ith value of the offsets data
 * field is returned.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @library matrix
 * @returns The data value (sint32) of the offset of the iterator.
 * @usage offset = #MUTIL_ITERATOR_GET_OFFSET(&iterator, index);#
 *
 * @param  ITERATOR    Pointer to iterator.
 * @param  INDEX       Index of the offset to retrieve.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
#define MUTIL_ITERATOR_GET_OFFSET( ITERATOR, INDEX )     \
  ( ITERATOR )->offsets_ptr[ ( INDEX ) ]


/** Return the base address of an iterator.
 * The return is cast to a pointer of type TYPE.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @library matrix
 * @returns The base\_ptr of the iterator, cast to (TYPE *).
 * @usage #base = MUTIL_ITERATOR_GET_BASE_ADDRESS(sint16, &iterator);#
 *
 * @param  TYPE        Type of the base data, such as uint32 or double.
 * @param  ITERATOR    Pointer to iterator.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 * @see MUTIL_ITERATOR_GET_TYPE
 */
#define MUTIL_ITERATOR_GET_BASE_ADDRESS( TYPE, ITERATOR )     \
  ( TYPE * ) ( ( ITERATOR )->base_ptr )


/** Extract a data entry of an iterator given an index.
 *
 * Given the arguments ITERATOR, INDEX, and TYPE, the following is returned:
 *
 *  *(TYPE *) (ITERATOR->base\_ptr + iterator->offsets.mat.s32mat.data[INDEX])
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @library matrix
 * @returns A data entry of the iterator, cast to (TYPE).
 * @usage #data = MUTIL_ITERATOR_EXTRACT_DATA(sint8, &iterator, index);#
 *
 * @param  TYPE             Type of the base data, such as uint32 or double.
 * @param  ITERATOR         Pointer to iterator.
 * @param  INDEX            Index of the offset data to retrieve.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_GET_TYPE
 */
#define MUTIL_ITERATOR_EXTRACT_DATA( TYPE, ITERATOR, INDEX )     \
  * ( TYPE * ) ( ( TYPE * ) ( ITERATOR )->base_ptr +             \
          MUTIL_ITERATOR_GET_OFFSET( ITERATOR, INDEX ) )


/** Return the mutil\_data\_type of an iterator.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_iter.h
 * @library matrix
 * @returns The mutil\_data\_type the iterator.
 * @usage #data_type = MUTIL_ITERATOR_GET_TYPE(&iterator);#
 *
 * @param  ITERATOR    Pointer to iterator.
 *
 * @see _mutil_iterator_dim
 * @see _mutil_iterator
 * @see mutil_iterator_malloc
 * @see mutil_iterator_malloc_with_matset
 * @see mutil_iterator_malloc_temp
 * @see mutil_iterator_free
 * @see mutil_iterator_assign_base
 * @see mutil_iterator_assign_base_for_iter_slice
 * @see mutil_iterator_set_length
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE
 * @see MUTIL_ITERATOR_MACRO_ASSIGN_BASE_FAST
 * @see MUTIL_ITERATOR_ASSIGN_DATA
 * @see MUTIL_ITERATOR_ASSIGN_OFFSET
 * @see MUTIL_ITERATOR_GET_LENGTH
 * @see MUTIL_ITERATOR_GET_OFFSET
 * @see MUTIL_ITERATOR_GET_BASE_ADDRESS
 * @see MUTIL_ITERATOR_EXTRACT_DATA
 */
#define MUTIL_ITERATOR_GET_TYPE( ITERATOR ) ( ITERATOR )->type


#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_ITER_H_*/
