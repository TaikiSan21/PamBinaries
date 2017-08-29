
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_intr.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_INTR_H_
#define IN_MAT_INTR_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "ut_err.h"
#include "mat_iter.h"
#include "mat_type.h"


/*
   Internal (private) matrix functions (mat_intr), internal structs
   and functions for iterators of matrix sets and matrices.
   This file contains function declarations for private functions that
   are internal to the mutils library, but that are used by more than
   one c file, and therefore cannot be static functions.
*/

#ifdef __cplusplus
extern "C" {
#endif


/** Set the iterator offset\_ptr.
 * The iterator offset\_ptr is set to the sum of the address
 * of the data field of the iterator univ\_mat offsets and the offset argument.
 * Length will reset the iterator length.  If the length argument is
 * zero, length will be set to (iterator->max\_length - offset).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intr.h
 * @source mat\_intr.c
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = mutil_iterator_assign_offset_base(&iterator, offset, 0);#
 *
 * @param  iterator    Pointer to iterator that is to be initialized
 *
 * @param  offset      Offset to add to the data field of mat\_base to
 *                     assign to the iterator base\_ptr.
 *
 * @param  length      Number of integer elements to iterate over.
 *                     If (length==0), length will be set to
 *                     (iterator->mat\_length - offset).
 *
 * @see localfn_set_univ_mat_nelem_zero
 * @see localfn_iterator_set_workspace
 * @see mutil_iterator_zero_base_data
 * @see mutil_iterator_copy_base_data
 *
 * @private
 */
extern mutil_errcode mutil_iterator_assign_offset_base(
  mutil_iterator * iterator, sint32 offset, sint32 length );


/** Assign the value of zero to all the data elements of the iterator.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intr.h
 * @source mat\_intr.c
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = mutil_iterator_zero_base_data(iterator);#
 *
 * @param  iterator    Pointer to iterator
 *
 * @see localfn_set_univ_mat_nelem_zero
 * @see localfn_iterator_set_workspace
 * @see mutil_iterator_assign_offset_base
 * @see mutil_iterator_copy_base_data
 *
 * @private
 */
extern mutil_errcode mutil_iterator_zero_base_data( mutil_iterator *iterator );


/** Copy the data from one iterator to another.
 * All the data elements pointed to by iter\_in, are copied to iter\_out.
 * The iterators must have identical dimensions and mutil\_data\_types.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intr.h
 * @source mat\_intr.c
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = mutil_iterator_copy_base_data(&iter_in, &iter_out);#
 *
 * @param  iter_in    Pointer to input iterator
 *
 * @param  iter_out   Pointer to output iterator.  Data pointed to by
 *                    iter\_in is copied to iter\_out.
 *
 * @see localfn_set_univ_mat_nelem_zero
 * @see localfn_iterator_set_workspace
 * @see mutil_iterator_assign_offset_base
 * @see mutil_iterator_zero_base_data
 *
 * @private
 */
extern mutil_errcode mutil_iterator_copy_base_data( mutil_iterator *iter_in,
  mutil_iterator *iter_out );


/* ****************************** */
/* definition of mutil_vector     */
/* ****************************** */


/** Struct that stores information about a vector.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intr.h
 * @source mat\_intr.h
 * @library matrix
 *
 * @same #typedef struct _mutil_vector mutil_vector;#
 *
 * @see vector_construct
 * @see vector_destruct
 * @see VECTOR_ASSIGN_DATA
 * @see VECTOR_EXTRACT_DATA
 * @see VECTOR_GET_BASE_TYPE
 * @see VECTOR_GET_LENGTH
 * @see VECTOR_GET_BASE_ADDRESS
 *
 * @private
 */
struct _mutil_vector{

  /** a vector matrix with nelem = nrow = length */
  univ_mat base;

  /** length of the vector matrix */
  sint32   length;
};


/* See above for documentation on this struct. */
typedef struct _mutil_vector mutil_vector;


/* **************************** */
/* constructor for mutil_vector */
/* **************************** */


/** Constructor for a mutils vector.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intr.h
 * @source mat\_intr.c
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = vector_construct(&vector, MUTIL_DOUBLE);#
 *
 * @param  iterator     Pointer to the vector to be initialized.
 *
 * @param  length       Length of the vector.
 *
 * @param  mutil_type   Type of the vector such as MUTIL\_SINT32.
 *
 * @see _mutil_vector
 * @see vector_destruct
 * @see VECTOR_ASSIGN_DATA
 * @see VECTOR_EXTRACT_DATA
 * @see VECTOR_GET_BASE_TYPE
 * @see VECTOR_GET_LENGTH
 * @see VECTOR_GET_BASE_ADDRESS
 *
 * @private
 */
extern mutil_errcode vector_construct(
  mutil_vector *vector, sint32 length, mutil_data_type mutil_type );


/* *************************** */
/* destructor for mutil_vector */
/* *************************** */

/** Destructor for a mutils vector.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intr.h
 * @source mat\_intr.c
 * @library matrix
 * @return Standard mutils error/OK code.
 * @usage #err = vector_destruct(&vector);#
 *
 * @param  vector     Pointer to the vector to be deallocated.
 *
 * @see _mutil_vector
 * @see vector_construct
 * @see VECTOR_GET_BASE_TYPE
 * @see VECTOR_GET_LENGTH
 * @see VECTOR_GET_BASE_ADDRESS
 * @private
 */
extern mutil_errcode vector_destruct( mutil_vector *vector );


/* ********************************************************* */
/* macros that assign fields of mutil_vector */
/* ********************************************************* */

/** Assign data to a mutil\_vector. No type checking is done.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intr.h
 * @source mat\_intr.c
 * @library matrix
 * @usage #VECTOR_ASSIGN_DATA(&vector, index, data);#
 *
 * @param  VECTOR     Pointer to a vector.
 * @param  INDEX      Index of the data to assign.
 * @param  VALUE      Value to assign.
 *
 * @see _mutil_vector
 * @see vector_construct
 * @see vector_destruct
 * @see VECTOR_EXTRACT_DATA
 * @see VECTOR_GET_LENGTH
 * @see VECTOR_GET_BASE_ADDRESS
 *
 * @private
 */
#define VECTOR_ASSIGN_DATA( VECTOR, INDEX, VALUE ) \
  MATUNIV_ELEM_ASSIGN_BY_INDEX( &( (VECTOR)->base ), INDEX, VALUE );


/* ************************** */
/* accessors for mutil_vector */
/* ************************** */


/** Extract data from a vector given the index.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intr.h
 * @source mat\_intr.c
 * @library matrix
 * @returns A data entry of the vector,  No casting is done.
 * @usage #data = VECTOR_EXTRACT_DATA(&vector, index);#
 *
 * @param  VECTOR     Pointer to a vector.
 * @param  INDEX      Index of the data to retrieve.
 *
 * @see _mutil_vector
 * @see vector_construct
 * @see vector_destruct
 * @see VECTOR_ASSIGN_DATA
 * @see VECTOR_GET_LENGTH
 * @see VECTOR_GET_BASE_ADDRESS
 *
 * @private
 */
#define VECTOR_EXTRACT_DATA( VECTOR, INDEX ) \
  MATUNIV_ELEM_BY_INDEX( &( (VECTOR)->base ), INDEX );


/** Return the type of the vector.
 * Type returned is of type mutil\_data\_type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intr.h
 * @source mat\_intr.c
 * @library matrix
 * @return Returns a mutil\_data\_type.
 * @usage #type = VECTOR_GET_BASE_TYPE(&vector);#
 *
 * @param  vector     Pointer to a vector.
 *
 * @see _mutil_vector
 * @see vector_construct
 * @see vector_destruct
 * @see VECTOR_ASSIGN_DATA
 * @see VECTOR_EXTRACT_DATA
 * @see VECTOR_GET_LENGTH
 * @see VECTOR_GET_BASE_ADDRESS
 *
 * @private
 */
#define VECTOR_GET_BASE_TYPE( VECTOR ) ( VECTOR )->base.type


/** Return the length of the vector.
 * Type returned is of sint32.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intr.h
 * @source mat\_intr.c
 * @library matrix
 * @return The length (sint32) of the vector.
 * @usage #length = VECTOR_GET_LENGTH(&vector);#
 *
 * @param  VECTOR     Pointer to a vector.
 *
 * @see _mutil_vector
 * @see vector_construct
 * @see vector_destruct
 * @see VECTOR_ASSIGN_DATA
 * @see VECTOR_EXTRACT_DATA
 * @see VECTOR_GET_BASE_TYPE
 * @see VECTOR_GET_BASE_ADDRESS
 *
 * @private
 */
#define VECTOR_GET_LENGTH( VECTOR ) ( VECTOR )->length


/** Macro for getting the address of the data field of the vector.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intr.h
 * @source mat\_intr.c
 * @library matrix
 * @return The address of the data field of the vector cast to (TYPE *).
 * @usage #err = VECTOR_GET_BASE_ADDRESS(uint8, &vector);#
 *
 * @param  TYPE       Type of the vector such as sint32 or double.
 * @param  VECTOR     Pointer to a vector.
 *
 * @see _mutil_vector
 * @see vector_construct
 * @see vector_destruct
 * @see VECTOR_ASSIGN_DATA
 * @see VECTOR_EXTRACT_DATA
 * @see VECTOR_GET_BASE_TYPE
 * @see VECTOR_GET_LENGTH
 *
 * @private
 */
#define VECTOR_GET_BASE_ADDRESS( TYPE, VECTOR ) \
  ( TYPE * ) MATUNIV_DATA( &( ( VECTOR )->base ) )


#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_INTR_H_*/
