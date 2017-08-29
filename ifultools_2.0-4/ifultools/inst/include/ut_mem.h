
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_mem.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_MEM_H_
#define IN_UT_MEM_H_

#include "str_type.h"
#include <stdlib.h>
#include <stdio.h>

/* This file contains functions and macros for
   linearly linked memory lists               */

#ifdef __cplusplus
extern "C" {
#endif


/** Macro for freeing memory list on encountering an error.
 * This macro frees the allocated data pointed to by the
 * registered members of a memory list and frees the list itself
 * upon encountering an error.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.h
 * @library matrix
 * @usage #MEMLIST_FREE_ON_ERROR(err, list);#
 * @param ERR_CODE   Mutil error code.
 * @param MEM_LIST   Memory list.
 * @see memlist_free
 * @see _memlist
 * @private
 */
#define MEMLIST_FREE_ON_ERROR( ERR_CODE, MEM_LIST ) \
  if ( ERR_CODE ) {                                 \
    MUTIL_FREE_WARN( memlist, MEM_LIST );           \
    return ( ERR_CODE );                            \
  }


/** Initialize memory list.
 * This macro is needed to initialize the memory list.
 * It sets the pointer of the list root to NULL, signifying an empty list.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.h
 * @library matrix
 * @usage #MEMLIST_INIT(list);#
 * @param list Memory list.
 * @see _memlist
 */
#define MEMLIST_INIT( MEM_LIST ) {               \
  ( MEM_LIST ).root   = ( memlist_node * ) NULL; \
  ( MEM_LIST ).length = ( sint32 ) 0;            \
}


/** Universal matrix memory allocation with memory management registration.
 * Allocate uninitialized storage in dynamic memory for a
 * matrix data array of the desired size and type, and put it and the
 * size parameter into the supplied matrix structure.
 * A copy of the pointer to the allocated memory along
 * with the memory data type (specified with \Ref{_memlist_type})
 * is then registered with a memory manager.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = matuniv_malloc_register(&matrix, nrow, ncol, MUTIL_DOUBLE, &list);#
 * @return     Standard mutils error/OK code.
 * @param matrix   Pointer to matrix to initialize.
 * @param nrow     Number of rows to create.
 * @param ncol     Number of columns to create.
 * @param type     Data type to create (specific type versions of this
 *                 function omit this argument).
 * @param list     Pointer to memory list.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_malloc_register(double_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_malloc_register(float_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_malloc_register(uint8_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_malloc_register(uint16_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_malloc_register(uint32_mat *matrix, sint32 nrow,  sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_malloc_register(sint16_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_malloc_register(sint32_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matcpx_malloc_register(dcomplex_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \end{itemize}
 * @see memlist_free
 * @see _memlist_node
 * @see _memlist_type
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_malloc_register(
  univ_mat        *matrix,
  sint32           nrow,
  sint32           ncol,
  mutil_data_type  type,
  memlist         *list );


/* This function is documented under matuniv_malloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_malloc_register(
  double_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list );


/* This function is documented under matuniv_malloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matflt_malloc_register(
  float_mat *matrix,
  sint32 nrow,
  sint32 ncol,
  memlist *list );


/* This function is documented under matuniv_malloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matu8_malloc_register(
  uint8_mat *matrix,
  sint32 nrow,
  sint32 ncol,
  memlist *list );


/* This function is documented under matuniv_malloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matu16_malloc_register(
  uint16_mat *matrix,
  sint32 nrow,
  sint32 ncol,
  memlist *list );


/* This function is documented under matuniv_malloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matu32_malloc_register(
  uint32_mat *matrix,
  sint32 nrow,
  sint32 ncol,
  memlist *list );


/* This function is documented under matuniv_malloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode mats16_malloc_register(
  sint16_mat *matrix,
  sint32 nrow,
  sint32 ncol,
  memlist *list );


/* This function is documented under matuniv_malloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode mats32_malloc_register(
  sint32_mat *matrix,
  sint32 nrow,
  sint32 ncol,
  memlist *list );


/* This function is documented under matuniv_malloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_malloc_register(
  dcomplex_mat *matrix,
  sint32 nrow,
  sint32 ncol,
  memlist *list );


/** Universal matrix memory reallocation with memory management registration.
 * Reallocate storage in dynamic memory for a preallocate matrix. The data
 * member of the corresponding member in the (memory manager) list is updated.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = matuniv_realloc_register(&matrix, nrow, ncol, &list);#
 * @return     Standard mutils error/OK code.
 * @param matrix   Pointer to matrix to initialize.
 * @param nrow     Number of rows to create.
 * @param ncol     Number of columns to create.
 * @param list     Pointer to memory list.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_realloc_register(double_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_realloc_register(float_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_realloc_register(uint8_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_realloc_register(uint16_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_realloc_register(uint32_mat *matrix, sint32 nrow,  sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_realloc_register(sint16_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_realloc_register(sint32_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matcpx_realloc_register(dcomplex_mat *matrix, sint32 nrow, sint32 ncol, memlist *list);#
 *   \end{itemize}
 * @see matuniv_realloc
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_realloc_register(
  univ_mat *matrix,
  sint32    nrow,
  sint32    ncol,
  memlist  *list );


/* This function is documented under matuniv_realloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_realloc_register(
  double_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list );


/* This function is documented under matuniv_realloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matflt_realloc_register(
  float_mat *matrix,
  sint32     nrow,
  sint32     ncol,
  memlist   *list );


/* This function is documented under matuniv_realloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matu8_realloc_register(
  uint8_mat *matrix,
  sint32     nrow,
  sint32     ncol,
  memlist   *list );


/* This function is documented under matuniv_realloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matu16_realloc_register(
  uint16_mat *matrix,
  sint32     nrow,
  sint32     ncol,
  memlist   *list );


/* This function is documented under matuniv_realloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matu32_realloc_register(
  uint32_mat *matrix,
  sint32     nrow,
  sint32     ncol,
  memlist   *list );


/* This function is documented under matuniv_realloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode mats16_realloc_register(
  sint16_mat *matrix,
  sint32     nrow,
  sint32     ncol,
  memlist   *list );


/* This function is documented under matuniv_realloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode mats32_realloc_register(
  sint32_mat *matrix,
  sint32     nrow,
  sint32     ncol,
  memlist   *list );


/* This function is documented under matuniv_realloc_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_realloc_register(
  dcomplex_mat *matrix,
  sint32        nrow,
  sint32        ncol,
  memlist      *list );


/** Matrix set memory allocation with memory management registration.
 * Allocate uninitialized storage in dynamic memory for a
 * matrix set of the specified dimensions.  Memory is allocated
 * to store the dimensions and the universal matrix headers in
 * the matrix array, but no storage is created for the universal
 * matrix data (use \Ref{matset_malloc_matrices} or
 * \Ref{matset_malloc_matrices_arbitrary_size} for that).
 * All memory associated with the matrix set can be freed
 * with the \Ref{memlist_member_free} function, including memory for the
 * internal universal matrices.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = matset_malloc_register(&mset, ndim, dims);#
 * @return         Standard mutils error/OK code.
 * @param mset     Pointer to matrix set to initialize.
 * @param ndim     Number of dimensions in the matrix set.
 * @param dims     Dimensions of the matrix set.
 * @param list     Pointer to memory list.
 * @see memlist_free
 * @see memlist_member_free
 * @see _memlist_node
 * @see _memlist_type
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_malloc_register(
  mat_set      *mset,
  const sint32  ndim,
  const sint32 *dims,
  memlist      *list );


/** kd-tree structure memory allocation with memory management registration.
 * Allocate internal storage in dynamic memory for a kd-tree structure and
 * register the structure with the memory manager. A pointer to an existing
 * \Ref{_mutil_kdtree} structure is given and memory for the data fields is
 * allocated with a call to \Ref{mutil_kdtree_malloc}. A pointer to the
 * structure pointer is stored with the memory manager. When the memory
 * management list is freed, or the corresponding member of the kd-tree
 * is deleted via \Ref{memlist_member_free}, all structure memory will be
 * freed with a call to \Ref{mutil_kdtree_free}.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = mutil_kdtree_malloc_register(&kdtree,&data,bsz,&mlist);#
 * @return                Standard mutils error/OK code.
 * @param kdt             Pointer to a \Ref{_mutil_kdtree} structure for which
 *                        memory is to be allocated for the internal fields.
 * @param points          Pointer to a universal matrix
 *                        of type MUTIL\_DOUBLE containing the data points
 *                        to be represented.
 * @param bucket_size     A positive integer specifying the bucket size of
 *                        the kd-tree. This gives the maximum number of
 *                        points in any leaf of the tree. See
 *                        \Ref{_mutil_kdtree}.
 * @see _mutil_kdtree
 * @see memlist_free
 * @see memlist_member_free
 */
MUTIL_LIBEXPORT mutil_errcode mutil_kdtree_malloc_register(
  mutil_kdtree      *kdt,
  const univ_mat    *points,
  const sint32       bucket_size,
  memlist           *list );


/** Free entire memory list and corresponding dynamically allocated data.
 * Free the data pointed to by the registered members of a memory manager
 * list as well as the list itself.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err = memlist_free(&list);#
 * @return         Standard mutils error/OK code.
 * @param list     Pointer to memory list.
 * @see MEMLIST_FREE_ON_ERROR
 * @see memlist_member_free
 * @see _memlist_node
 */
MUTIL_LIBEXPORT mutil_errcode memlist_free( memlist *list );


/** Unregister a member from a memory manager list and free any
 * previously allocated memory associated with that member.
 * Unregisters a member and frees any corresponding dynamically
 * allocated memory associated with that member.
 * The function verifies that the memory pointed to by data
 * was previously registered with the the memory list prior to
 * any attempt to free the memory.
 *
 * For matrix sets, both the universal matrices within the matrix set and the matrix
 * set header will be freed. Attempts to free memory associated with any
 * member of a memory list should only be done through this function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = memlist_member_free(&data, &list);#
 * @return         Standard mutils error/OK code.
 * @param data   Pointer to data to be freed.
 * @param list   Pointer to memory list.
 * @see _memlist_node
 */
MUTIL_LIBEXPORT mutil_errcode memlist_member_free(
  void    *data,
  memlist *list );


/** Unregister a member from a memory manager list.
 * Unregisters a member from a memory manager list but does
 * not free any (dynamically allocated) memory associated with that member.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = memlist_member_unregister(&data, &list);#
 * @return         Standard mutils error/OK code.
 * @param data   Pointer to pre-registered dynamically allocated memory.
 * @param list   Pointer to memory list.
 * @see _memlist_node
 */
MUTIL_LIBEXPORT mutil_errcode memlist_member_unregister(
  void    *data,
  memlist *list );


/** Add a new member to a memory manager list.
 * A pointer to previously allocated data and corresponding
 * memory data type are used to register a new member in the
 * the memory manager list.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = memlist_member_register(&list, &data, type);#
 * @return         Standard mutils error/OK code.
 * @param list     Pointer to memory list.
 * @param data     Pointer to previously allocated data.
 * @param type     Data memory type.
 * @see _memlist_node
 * @see _memlist_type
 * @private
 */
MUTIL_LIBEXPORT mutil_errcode memlist_member_register(
  memlist      *list,
  void         *data,
  memlist_type  type );


/** Initialize a matrix with given size, type, and data
 * and register it with the memory manager.
 * This function takes a data pointer, size, and type, and
 * and puts this information into a matrix structure.
 * The given data array is incorporated directly,
 * and the data is not copied or reallocated. The allocated
 * universal matrix is then registered with the memory manager.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = matuniv_wrap_data_register(&my_umat, &my_dblarr, nrow, ncol, MUTIL_DOUBLE, &list);#
 * @return    Standard mutils error/OK code.
 * @param   matrix  Pointer to the matrix to initialize.
 * @param   data    Pointer to the data to be incorporated.
 * @param   nrow    Number of rows in the matrix.
 * @param   ncol    Number of columns in the matrix.
 * @param   type    Data type of matrix (specific type versions of this
 *   function omit this argument).
 * @param list   Pointer to memory list.
 * @same \begin{itemize}
 *    \item #MUTIL_LIBEXPORT mutil_errcode matdbl_wrap_data_register(double_mat *matrix, double *data, sint32 nrow, sint32 ncol, memlist *list);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode matflt_wrap_data_register(float_mat *matrix, float *data, sint32 nrow, sint32 ncol, memlist *list);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode matu8_wrap_data_register(uint8_mat *matrix, uint8 *data, sint32 nrow, sint32 ncol, memlist *list);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode matu16_wrap_data_register(uint16_mat *matrix, uint16 *data, sint32 nrow, sint32 ncol, memlist *list);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode matu32_wrap_data_register(uint32_mat *matrix, uint32 *data, sint32 nrow, sint32 ncol, memlist *list);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode mats16_wrap_data_register(sint16_mat *matrix, sint16 *data, sint32 nrow, sint32 ncol, memlist *list);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode mats32_wrap_data_register(sint32_mat *matrix, sint32 *data, sint32 nrow, sint32 ncol, memlist *list);#
 * \end{itemize}
 * @see matuniv_wrap_matrix
 * @see matuniv_wrap_matrix_register
 * @see matuniv_malloc
 * @see _mutil_data_type
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_wrap_data_register(
  univ_mat        *matrix,
  void            *data,
  sint32           nrow,
  sint32           ncol,
  mutil_data_type  type,
  memlist         *list );


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_wrap_data_register(
  double_mat *matrix,
  double     *data,
  sint32      nrow,
  sint32      ncol,
  memlist    *list );


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matflt_wrap_data_register(
  float_mat *matrix,
  float     *data,
  sint32     nrow,
  sint32     ncol,
  memlist   *list );


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matu8_wrap_data_register(
  uint8_mat *matrix,
  uint8     *data,
  sint32     nrow,
  sint32     ncol,
  memlist   *list );


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matu16_wrap_data_register(
  uint16_mat *matrix,
  uint16     *data,
  sint32      nrow,
  sint32      ncol,
  memlist    *list );


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode matu32_wrap_data_register(
  uint32_mat *matrix,
  uint32     *data,
  sint32      nrow,
  sint32      ncol,
  memlist    *list );


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode mats16_wrap_data_register(
  sint16_mat *matrix,
  sint16     *data,
  sint32      nrow,
  sint32      ncol,
  memlist    *list );


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
MUTIL_LIBEXPORT mutil_errcode mats32_wrap_data_register(
  sint32_mat *matrix,
  sint32     *data,
  sint32      nrow,
  sint32      ncol,
  memlist    *list );


/** Initialize a universal matrix from a non-universal matrix and
 * register it with the memory manager.
 * This function takes an existing (non-universal) matrix data structure,
 * and puts its information into a universal matrix header.
 * The matrix's underlying flat data array is incorporated directly,
 * and the data is not copied or reallocated. The universal matrix is then
 * registered with the memory manager as follows:
 * The address of the non-universal matrix is searched for in the memory list
 * to see if it has been registered. If so, the corresponding
 * node is updated (by replacing the data pointer)
 * with the new universal matrix pointer and memory type.
 * Otherwise, a new node for the universal matrix is created.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = matuniv_wrap_matrix_register(&my_umat, &my_dblmat, MUTIL_DOUBLE, &list);#
 * @return              Standard mutils error/OK code.
 * @param   umat        Pointer to the universal matrix to initialize.
 * @param   matrix_ptr  Pointer to the matrix to be incorporated.
 * @param   type        Data type of the matrix.
 * @param list   Pointer to memory list.
 * @see matuniv_malloc
 * @see matuniv_wrap_data
 * @see _mutil_data_type
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_wrap_matrix_register(
  univ_mat        *umat,
  void            *matrix_ptr,
  mutil_data_type  type,
  memlist         *list );


/** Initialize a universal matrix from a universal matrix and register
 * it with the memory manager.
 * This function takes an existing universal matrix data structure,
 * and puts its information into a universal matrix header.
 * The matrices underlying flat data array is incorporated directly,
 * and the data is not copied or reallocated. The universal matrix is then
 * registered with the memory manager.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #errcode = matuniv_wrap_univ_matrix(&umat1, &umat2, &list);#
 * @return              Standard mutils error/OK code.
 * @param   mat1        Pointer to the universal matrix to initialize.
 * @param   mat2        Pointer to the universal matrix to be incorporated.
 * @param list   Pointer to memory list.
 * @see matuniv_malloc
 * @see matuniv_wrap_data
 * @see _mutil_data_type
 * @see matuniv_wrap_matrix
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_wrap_univ_matrix_register(
  univ_mat *mat1,
  univ_mat *mat2,
  memlist  *list );


/** Allocate dynamic memory and register it with the memory manager.
 * This function works like the ANSI C malloc function, to
 * allocate dynamic memory with no particular initialized values.
 * It is an error to request a zero or negative size.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = mutil_malloc_register(3 * sizeof(int), &my_int_ptr, &list);#
 * @return     Standard mutils error/OK code.
 * @param  size    Number of bytes of memory to allocate.
 * @param  data    Pointer to the allocated memory.
 * @param list   Pointer to memory list.
 * @see mutil_malloc
 * @see mutil_calloc_register
 * @see mutil_realloc_register
 * @see memlist_free
 */
MUTIL_LIBEXPORT mutil_errcode mutil_malloc_register(
  sint32    size,
  void    **data,
  memlist  *list );


/** Allocate dynamic memory, initialized to zero and register
 * it with the memory manager.
 * This function works like the ANSI C calloc function, to
 * allocate dynamic memory and initialize it to zero.
 * It is an error to request a zero or negative size.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = mutil_calloc(3, sizeof(int), &my_int_ptr, &list);#
 * @return    Standard mutils error/OK code.
 * @param  nobj    Number of objects to allocate.
 * @param  size    Size of each object to allocate.
 * @param  data    Pointer to the allocated memory.
 * @param list   Pointer to memory list.
 * @see mutil_calloc
 * @see mutil_malloc_register
 * @see mutil_realloc_register
 * @see memlist_free
 */
MUTIL_LIBEXPORT mutil_errcode mutil_calloc_register(
  sint32    nobj,
  size_t    size,
  void    **data,
  memlist  *list );


/** Change the size of dynamic memory and register the memory
 * with the memory manager.
 * This function works like the ANSI C realloc function, to
 * change the size of the previously allocated memory. The contents
 * are left unchanged up to the minimum of the old and new sizes, and
 * if the new size is larger, the new space is uninitialized.  The new
 * space is not guaranteed to be at the same location as the old space.
 * It is an error to request a zero or negative size.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = mutil_realloc(&my_int_ptr, 10 * sizeof(int), 3 * sizeof(int), &list);#
 * @return    Standard mutils error/OK code.
 * @param  data       Pointer to memory to be reallocated, changed on
 *    return to new location.
 * @param  new_size   New size, in bytes.
 * @param  old_size   Old size, in bytes.
 * @param list   Pointer to memory list.
 * @see mutil_realloc
 * @see mutil_malloc_register
 * @see mutil_calloc_register
 * @see memlist_free
*/
MUTIL_LIBEXPORT mutil_errcode mutil_realloc_register(
  void    **data,
  sint32    new_size,
  sint32    old_size,
  memlist  *list );


/** Check the existence of a member in the memory manager.
 * This function is used to determine whether or not a specified address
 * has already been assigned to any member of a memory manger list.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #found = memlist_member_exist(&data, &list);#
 * @return        Boolean. If TRUE, the data address was found in the list.
 * @param  data   Pointer to data. The address of this memory is compared to
 * @param  list   Pointer to memory list.
 * @see _memlist_node
*/
MUTIL_LIBEXPORT boolean memlist_member_exist(
  void    *data,
  memlist *list );


/** Validate a memory list.
 * This function validates a memory list by checking that the length of the
 * list is valid, i.e. that the number valid list members is equal to the
 * length member of the \Ref{_memlist} structure. In addition, the memory
 * corresponding to each member of the list is checked by calling the
 * appropriate MUTILS validation function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err = memlist_validate(&list);#
 * @return  Standard mutils error/OK code.
 * @param  list   Pointer to memory list.
 * @see _memlist
 * @see matset_validate
 * @see matuniv_validate
*/
MUTIL_LIBEXPORT mutil_errcode memlist_validate( memlist *list );


/** Print a report of memory list contents.
 * This function prints a report of the contents of a memory list.
 * The member index, the registered member address, and the member type
 * are all printed. If the member is a matrix, the address of the
 * corresponding data block is also displayed. This function will typically
 * be used to help track down memory leaks.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err = memlist_report( out_file, &list );#
 * @return  Standard mutils error/OK code.
 * @param out_file Pointer to a file for output (use stdout for display to screen).
 * @param  list   Pointer to memory list.
 * @see _memlist
 * @see memlist_validate
*/
MUTIL_LIBEXPORT mutil_errcode memlist_print( FILE *out_file, memlist *list );


#ifdef __cplusplus
   }
#endif

#endif /* ifdef IN_UT_MEM_H_ */
