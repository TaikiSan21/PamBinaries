
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_alloc.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_ALLOC_H
#define IN_UT_ALLOC_H

#include "ut_plat.h"
#include "ut_err.h"
#include "ut_type.h"

#include <stdlib.h>

/*
   This file contains function declarations for dynamic memory
   management.  The functions are user-definable to let different
   applications have different dynamic memory management, and all
   dynamic memory allocation in the mutils library goes through these
   functions.
*/

#ifdef __cplusplus
extern "C" {
#endif


/** Allocate dynamic memory.
 * This function works like the ANSI C malloc function, to
 * allocate dynamic memory with no particular initialized values.
 * It is an error to request a zero or negative size.
 *
 * All dynamic memory allocation and deallocation in the mutils
 * library goes through the mutil\_malloc, mutil\_calloc,
 * mutil\_realloc, and mutil\_free functions.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_alloc.h
 * @source ut\_alloc.c
 * @library wrap
 * @usage  #err_code = mutil_malloc(3 * sizeof(int), &my_int_ptr);#
 * @return     Standard mutils error/OK code.
 * @param  size    Number of bytes of memory to allocate.
 * @param  data    Pointer to the allocated memory.
 * @see mutil_calloc
 * @see mutil_realloc
 * @see mutil_free
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_malloc( sint32 size, void **data );


/** Allocate dynamic memory, initialized to zero.
 * This function works like the ANSI C calloc function, to
 * allocate dynamic memory and initialize it to zero.
 * It is an error to request a zero or negative size.
 *
 * All dynamic memory allocation and deallocation in the mutils
 * library goes through the mutil\_malloc, mutil\_calloc,
 * mutil\_realloc, and mutil\_free functions.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_alloc.h
 * @source ut\_alloc.c
 * @library wrap
 * @usage  #err_code = mutil_calloc(3, sizeof(int), &my_int_ptr);#
 * @return    Standard mutils error/OK code.
 * @param  nobj    Number of objects to allocate.
 * @param  size    Size of each object to allocate.
 * @param  data    Pointer to the allocated memory.
 * @see mutil_malloc
 * @see mutil_realloc
 * @see mutil_free
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_calloc( sint32 nobj, size_t size,
  void **data );


/** Change the size of dynamic memory.
 * This function works like the ANSI C realloc function, to
 * change the size of the previously allocated memory. The contents
 * are left unchanged up to the minimum of the old and new sizes, and
 * if the new size is larger, the new space is uninitialized.  The new
 * space is not guaranteed to be at the same location as the old space.
 * It is an error to request a zero or negative size.
 *
 * All dynamic memory allocation and deallocation in the mutils
 * library goes through the mutil\_malloc, mutil\_calloc,
 * mutil\_realloc, and mutil\_free functions.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_alloc.h
 * @source ut\_alloc.c
 * @library wrap
 * @usage  #err_code = mutil_realloc(&my_int_ptr, 10 * sizeof(int), 3 * sizeof(int));#
 * @return    Standard mutils error/OK code.
 * @param   data       Pointer to memory to be reallocated, changed on
 *    return to new location.
 * @param   new_size   New size, in bytes.
 * @param   old_size   Old size, in bytes.
 * @see mutil_malloc
 * @see mutil_calloc
 * @see mutil_free
*/
MUTIL_WRAPEXPORT mutil_errcode mutil_realloc( void **data, sint32 new_size,
  sint32 old_size );


/** Deallocate previously allocated dynamic memory.
 * This function works like the ANSI C free function, to
 * deallocate previously allocated dynamic memory.
 *
 * All dynamic memory allocation and deallocation in the mutils
 * library goes through the mutil\_malloc, mutil\_calloc,
 * mutil\_realloc, and mutil\_free functions.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_alloc.h
 * @source ut\_alloc.c
 * @library wrap
 * @usage  #err_code = mutil_free(my_int_ptr, size);#
 * @return    Standard mutils error/OK code.
 * @param   data       Pointer to memory to be deallocated.
 * @param   old_size   Data size, in bytes.
 * @see mutil_malloc
 * @see mutil_calloc
 * @see mutil_realloc
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_free( void *data, sint32 old_size );


#ifdef __cplusplus
}
#endif

#endif /* ifndef IN_UT_ALLOC_H*/
