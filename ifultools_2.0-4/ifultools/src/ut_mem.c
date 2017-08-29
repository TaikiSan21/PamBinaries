
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/ut_mem.c $: $Revision: #1 $, $Date: 2008/03/21 $";

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_limit.h"
#include "ut_type.h"

#include "mat_type.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "mat_set.h"

#include "ut_mem.h"
#include "ut_kdtr.h"

/**********************************/
/* STATIC FUNCTIONS: DECLARATIONS */
/**********************************/

/* These functions documented below */


static boolean localfn_memlist_member_exist_node(
  void*,
  memlist*,
  memlist_node** );


static mutil_errcode localfn_memlist_member_validate(
  memlist_node *node );


static void localfn_memlist_member_data_free(
  memlist_node *node );


/**********************************/
/* STATIC MACRO DEFINITIONS       */
/**********************************/


/** Template macro for matrix wrap functions with memory
 * management registration.
 * Macro that expands to the body of a non-universal matrix
 * wrap\_data\_list function, such as matdbl\_wrap\_data\_list.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.h
 * @library matrix
 * @param MAT_FN_PREFIX    Prefix for functions for this matrix type.
 * @param MAT_PTR          Matrix pointer to allocate (function argument).
 * @param DATA_PTR         Data pointer to wrap (function argument).
 * @param MAT_NROW         Number of rows (function argument).
 * @param MAT_NCOL         Number of columns (function argument).
 * @param MEM_TYPE         Memory type (\Ref{_memlist_type}).
 * @param MEM_LIST_PTR     Pointer to memory list.
 * @usage Body of the matdbl\_wrap\_data\_list function:
 *     #TMPL_MAT_WRAP_DATA_LIST(matdbl, matrix, data, nrow, ncol, MEMTYPE_MATDBL, &list);#
 * @see _memlist_type
 * @see _memlist_node
 * @private
 */
#define TMPL_MAT_WRAP_DATA_LIST( MAT_FN_PREFIX, MAT_PTR, DATA_PTR,           \
  MAT_NROW, MAT_NCOL, MEM_TYPE, MEM_LIST_PTR )                               \
                                                                             \
  mutil_errcode err;                                                         \
  memlist_node *node;                                                        \
                                                                             \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_wrap_data_register()" );                \
                                                                             \
  /* Make sure the address of the list is not NULL */                        \
                                                                             \
  if ( MEM_LIST_PTR == ( memlist* ) NULL ) {                                 \
    MUTIL_ERROR( "Pointer to memory list is NULL" );                         \
    return MUTIL_ERR_NULL_POINTER;                                           \
  }                                                                          \
                                                                             \
  /* Validate the memory list */                                             \
  err = memlist_validate( MEM_LIST_PTR );                                    \
  if ( err ) {                                                               \
    MUTIL_ERROR( "Invalid memory list" );                                    \
    return err;                                                              \
  }                                                                          \
                                                                             \
  /* wrap data in matrix */                                                  \
                                                                             \
  err = MAT_FN_PREFIX ## _wrap_data( MAT_PTR,                                \
    DATA_PTR, MAT_NROW, MAT_NCOL );                                          \
  if ( err ) return err;                                                     \
                                                                             \
  /* register the new memory if it doesn't already exist */                  \
                                                                             \
  REPLACE_OR_ADD_MEMLIST_MEMBER( MAT_PTR, MEM_TYPE,                          \
    DATA_PTR, MEM_LIST_PTR, node, err );                                     \
                                                                             \
  MUTIL_TRACE( "Done with " #MAT_FN_PREFIX "_wrap_data_register()" );            \
  return MUTIL_ERR_OK


/** Template macro for matrix malloc list functions with memory
 * management registration.
 * Macro that expands to the body of a non-universal matrix
 * malloc\_list function, such as matdbl\_malloc\_list.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @param MAT_FN_PREFIX    Prefix for functions for this matrix type.
 * @param MAT_PTR          Matrix pointer to allocate (function argument).
 * @param MAT_NROW         Number of rows (function argument).
 * @param MAT_NCOL         Number of columns (function argument).
 * @param MEM_TYPE         Memory type (\Ref{_memlist_type}).
 * @param MEM_LIST_PTR     Pointer to memory list.
 * @usage Body of the matdbl\_malloc\_list function:
 *     #TMPL_MAT_MALLOC_LIST(matdbl, matrix, nrow, ncol, MEMTYPE_MATDBL, list);#
 * @see _memlist_type
 * @see _memlist_node
 * @private
 */
#define TMPL_MAT_MALLOC_LIST( MAT_FN_PREFIX, MAT_PTR, MAT_NROW, MAT_NCOL,    \
  MEM_TYPE, MEM_LIST_PTR )                                                   \
  mutil_errcode err;                                                         \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_malloc_register()" );                   \
                                                                             \
  /* Make sure the address of the list is not NULL */                        \
                                                                             \
  if ( MEM_LIST_PTR == ( memlist* ) NULL ) {                                 \
    MUTIL_ERROR( "Pointer to memory list is NULL" );                         \
    return MUTIL_ERR_NULL_POINTER;                                           \
  }                                                                          \
                                                                             \
  /* Validate the memory list */                                             \
  err = memlist_validate( MEM_LIST_PTR );                                    \
  if ( err ) {                                                               \
    MUTIL_ERROR( "Invalid memory list" );                                    \
    return err;                                                              \
  }                                                                          \
                                                                             \
  /* ensure that memory is not already registered */                         \
  /* with the memory manager */                                              \
                                                                             \
  if ( memlist_member_exist( MAT_PTR, MEM_LIST_PTR ) ) {                     \
    MUTIL_ERROR( "Memory already registered with the memory manager" );      \
    return MUTIL_ERR_ILLEGAL_ADDRESS;                                        \
  }                                                                          \
                                                                             \
  /* malloc space for matrix */                                              \
                                                                             \
  err = MAT_FN_PREFIX ## _malloc( MAT_PTR, MAT_NROW, MAT_NCOL );             \
  if ( err ) return err;                                                     \
                                                                             \
  err = memlist_member_register( MEM_LIST_PTR, MAT_PTR, MEM_TYPE );               \
  if ( err ) {                                                               \
    MUTIL_FREE_WARN( MAT_FN_PREFIX, MAT_PTR );                               \
    return err;                                                              \
  }                                                                          \
                                                                             \
  MUTIL_TRACE( "Done with " #MAT_FN_PREFIX "_malloc_register()" );               \
  return MUTIL_ERR_OK


/** Template macro for matrix realloc list functions with memory
 * management registration.
 * Macro that expands to the body of a non-universal matrix
 * realloc\_list function, such as matdbl\_realloc\_list.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @param MAT_FN_PREFIX    Prefix for functions for this matrix type.
 * @param MAT_PTR          Matrix pointer to allocate (function argument).
 * @param MAT_NROW         Number of rows (function argument).
 * @param MAT_NCOL         Number of columns (function argument).
 * @param MEM_LIST_PTR     Pointer to memory list.
 * @usage Body of the matdbl\_realloc\_list function:
 *     #TMPL_MAT_REALLOC_LIST(matdbl, matrix, nrow, ncol, list);#
 * @see _memlist_type
 * @see _memlist_node
 * @private
 */
#define TMPL_MAT_REALLOC_LIST( MAT_FN_PREFIX, MAT_PTR, MAT_NROW,             \
  MAT_NCOL, MEM_LIST_PTR )                                                   \
  mutil_errcode    err;                                                      \
  memlist_node    *node;                                                     \
                                                                             \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_realloc_register()" );                  \
                                                                             \
  /* Make sure the address of the list is not NULL */                        \
                                                                             \
  if ( MEM_LIST_PTR == ( memlist* ) NULL ) {                                 \
    MUTIL_ERROR( "Pointer to memory list is NULL" );                         \
    return MUTIL_ERR_NULL_POINTER;                                           \
  }                                                                          \
                                                                             \
  /* Validate the memory list */                                             \
  err = memlist_validate( MEM_LIST_PTR );                                    \
  if ( err ) {                                                               \
    MUTIL_ERROR( "Invalid memory list" );                                    \
    return err;                                                              \
  }                                                                          \
                                                                             \
  /* locate the appropriate node in the memory list */                       \
                                                                             \
  if ( localfn_memlist_member_exist_node( MAT_PTR,                           \
                      MEM_LIST_PTR,                                          \
                      &node ) ) {                                            \
                                                                             \
    /* realloc memory */                                                     \
                                                                             \
    err = MAT_FN_PREFIX ## _realloc( MAT_PTR, MAT_NROW, MAT_NCOL );          \
    if ( err ) return err;                                                   \
                                                                             \
    /* replace pointer in data member with new pointer */                    \
                                                                             \
    ( node )->data = MAT_PTR;                                                \
  }                                                                          \
  else{                                                                      \
                                                                             \
    MUTIL_ERROR( "The data to be replaced does not exist in memory list" );  \
    return MUTIL_ERR_ILLEGAL_ADDRESS;                                        \
  }                                                                          \
                                                                             \
  MUTIL_TRACE( "Done with " #MAT_FN_PREFIX "_realloc_register()" );              \
  return MUTIL_ERR_OK


/** Macro for updating the data in a memory list node.
 *
 * Previously allocated data is searched for in the memory list
 * to see if it has been registered. If so, the corresponding
 * node is updated with the new data pointer and memory type,
 * otherwise a new node is created.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @param NEW_DATA_PTR     Pointer to the new data.
 * @param NEW_MEMTYPE      Memory type of new data.
 * @param OLD_DATA_PTR     Pointer to the old data that may or may not
 *                         have already been registered with the memory
 *                         manager.
 * @param MEMLIST          Memory list.
 * @param MEMNODE          Memory list node.
 * @param MEMLIST          Memory list.
 * @param ERROR            MUTILS error code.

 * @usage #REPLACE_OR_ADD_MEMLIST_MEMBER( matrix, MEMTYPE_MATUNIV, data, list, node, err );#
 * @see _memlist_type
 * @see _memlist_node
 * @private
 */
#define REPLACE_OR_ADD_MEMLIST_MEMBER( NEW_DATA_PTR,         \
  NEW_MEMTYPE,                                               \
  OLD_DATA_PTR,                                              \
  MEMLIST,                                                   \
  MEMNODE,                                                   \
  ERROR )                                                    \
                                                             \
  if ( localfn_memlist_member_exist_node( OLD_DATA_PTR,      \
    MEMLIST, &( MEMNODE ) ) ) {                              \
    ( MEMNODE )->data = ( NEW_DATA_PTR );                    \
    ( MEMNODE )->type = ( memlist_type ) NEW_MEMTYPE;        \
  }                                                          \
  else{                                                      \
                                                             \
    /* Register allocated memory with memory list */         \
                                                             \
    ERROR = memlist_member_register( MEMLIST,                     \
      NEW_DATA_PTR, NEW_MEMTYPE );                           \
    if ( ERROR ) return ERROR;                               \
  }


/**********************************/
/* LIBRARY FUNCTIONS: DEFINITIONS */
/**********************************/


/* This function is documented under matuniv_malloc_register in ut_mem.h */
mutil_errcode matdbl_malloc_register(
  double_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_MALLOC_LIST( matdbl, matrix, nrow, ncol, MEMTYPE_MATDBL, list );
}


/* This function is documented under matuniv_malloc_register in ut_mem.h */
mutil_errcode matflt_malloc_register(
  float_mat *matrix,
  sint32     nrow,
  sint32     ncol,
  memlist   *list )
{
  TMPL_MAT_MALLOC_LIST( matflt, matrix, nrow, ncol, MEMTYPE_MATFLT, list );
}


/* This function is documented under matuniv_malloc_register in ut_mem.h */
mutil_errcode matu8_malloc_register(
  uint8_mat *matrix,
  sint32     nrow,
  sint32     ncol,
  memlist   *list )
{
  TMPL_MAT_MALLOC_LIST( matu8, matrix, nrow, ncol, MEMTYPE_MATU8, list );
}


/* This function is documented under matuniv_malloc_register in ut_mem.h */
mutil_errcode matu16_malloc_register(
  uint16_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_MALLOC_LIST( matu16, matrix, nrow, ncol, MEMTYPE_MATU16, list );
}


/* This function is documented under matuniv_malloc_register in ut_mem.h */
mutil_errcode matu32_malloc_register(
  uint32_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_MALLOC_LIST( matu32, matrix, nrow, ncol, MEMTYPE_MATU32, list );
}


/* This function is documented under matuniv_malloc_register in ut_mem.h */
mutil_errcode mats16_malloc_register(
  sint16_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_MALLOC_LIST( mats16, matrix, nrow, ncol, MEMTYPE_MATS16, list );
}


/* This function is documented under matuniv_malloc_register in ut_mem.h */
mutil_errcode mats32_malloc_register(
  sint32_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_MALLOC_LIST( mats32, matrix, nrow, ncol, MEMTYPE_MATS32, list );
}


/* This function is documented under matuniv_malloc_register in ut_mem.h */
mutil_errcode matcpx_malloc_register(
  dcomplex_mat *matrix,
  sint32        nrow,
  sint32        ncol,
  memlist      *list )
{
  TMPL_MAT_MALLOC_LIST( matcpx, matrix, nrow, ncol, MEMTYPE_MATCPX, list );
}


/* This function is documented under matuniv_realloc_register in ut_mem.h */
mutil_errcode matdbl_realloc_register(
  double_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_REALLOC_LIST( matdbl, matrix, nrow, ncol, list );
}


/* This function is documented under matuniv_realloc_register in ut_mem.h */
mutil_errcode matflt_realloc_register(
  float_mat *matrix,
  sint32     nrow,
  sint32     ncol,
  memlist   *list )
{
  TMPL_MAT_REALLOC_LIST( matflt, matrix, nrow, ncol, list );
}


/* This function is documented under matuniv_realloc_register in ut_mem.h */
mutil_errcode matu8_realloc_register(
  uint8_mat *matrix,
  sint32     nrow,
  sint32     ncol,
  memlist   *list )
{
  TMPL_MAT_REALLOC_LIST( matu8, matrix, nrow, ncol, list );
}


/* This function is documented under matuniv_realloc_register in ut_mem.h */
mutil_errcode matu16_realloc_register(
  uint16_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_REALLOC_LIST( matu16, matrix, nrow, ncol, list );
}


/* This function is documented under matuniv_realloc_register in ut_mem.h */
mutil_errcode matu32_realloc_register(
  uint32_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_REALLOC_LIST( matu32, matrix, nrow, ncol, list );
}


/* This function is documented under matuniv_realloc_register in ut_mem.h */
mutil_errcode mats16_realloc_register(
  sint16_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_REALLOC_LIST( mats16, matrix, nrow, ncol, list );
}


/* This function is documented under matuniv_realloc_register in ut_mem.h */
mutil_errcode mats32_realloc_register(
  sint32_mat *matrix,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_REALLOC_LIST( mats32, matrix, nrow, ncol, list );
}


/* This function is documented under matuniv_realloc_register in ut_mem.h */
mutil_errcode matcpx_realloc_register(
  dcomplex_mat *matrix,
  sint32        nrow,
  sint32        ncol,
  memlist      *list )
{
  TMPL_MAT_REALLOC_LIST( matcpx, matrix, nrow, ncol, list );
}


/* Matrix set memory allocation with */
/* memory management registration.   */
/* Documented in ut_mem.h            */
/* Written by Keith L. Davidson      */
mutil_errcode matset_malloc_register(
  mat_set      *mset,
  const sint32  ndim,
  const sint32 *dims,
  memlist      *list )
{
  mutil_errcode err;

  MUTIL_TRACE( "Start matset_malloc_register()" );

  /* avoid lint warning */
  (void) whatssi;

  /* Make sure the address of the list is not NULL */

  if ( list == ( memlist* ) NULL ) {
    MUTIL_ERROR( "Pointer to memory list is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* ensure that memory is not already registered
     with the memory manager */

  if ( memlist_member_exist( mset, list ) ) {
    MUTIL_ERROR( "Memory already registered with the memory manager" );
    return MUTIL_ERR_ILLEGAL_ADDRESS;
  }

  /* Allocate memory for matrix set */
  err = matset_malloc( mset, ndim, dims );
  if ( err ) return err;

  /* Register allocated memory with memory list */
  err = memlist_member_register( list, mset, MEMTYPE_MATSET );
  if ( err ) {
    MUTIL_FREE_WARN( matset, mset );
    return err;
  }

  MUTIL_TRACE( "Done with matset_malloc_register()" );
  return MUTIL_ERR_OK;
}


/* Matrix set memory allocation with */
/* memory management registration.   */
/* Documented in ut_mem.h            */
/* Written by Keith L. Davidson      */
MUTIL_LIBEXPORT mutil_errcode mutil_kdtree_malloc_register(
  mutil_kdtree      *kdt,
  const univ_mat    *points,
  const sint32       bucket_size,
  memlist           *list )
{
  mutil_errcode err;

  MUTIL_TRACE( "Start mutil_kdtree_malloc_register()" );

  /* Make sure the address of the list is not NULL */

  if ( list == (memlist*) NULL ) {
    MUTIL_ERROR( "Pointer to memory list is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* ensure that memory is not already registered
     with the memory manager */

  if ( memlist_member_exist( kdt, list ) ) {
    MUTIL_ERROR( "Memory already registered with the memory manager" );
    return MUTIL_ERR_ILLEGAL_ADDRESS;
  }

  /* Allocate memory for matrix set */
  err = mutil_kdtree_malloc( kdt, points, bucket_size );
  if ( err ) return err;

  /* Register allocated memory with memory list */
  err = memlist_member_register( list, kdt, MEMTYPE_KDTREE );
  if ( err ) {
    MUTIL_FREE_WARN( mutil_kdtree, kdt );
    return err;
  }

  MUTIL_TRACE( "Done with mutil_kdtree_malloc_register()" );
  return MUTIL_ERR_OK;
}



/* Universal matrix memory allocation with */
/* memory management registration.         */
/* Documented in ut_mem.h                  */
/* Written by William Constantine          */
mutil_errcode matuniv_malloc_register(
  univ_mat        *matrix,
  sint32           nrow,
  sint32           ncol,
  mutil_data_type  type,
  memlist         *list )
{
  mutil_errcode err;

  MUTIL_TRACE( "Start matuniv_malloc_register" );

  /* Make sure the address of the list is not NULL */

  if ( list == ( memlist* ) NULL ) {
    MUTIL_ERROR( "Pointer to memory list is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* ensure that memory is not already registered
     with the memory manager */

  if ( memlist_member_exist( matrix, list ) ) {
    MUTIL_ERROR( "Memory already registered with the memory manager" );
    return MUTIL_ERR_ILLEGAL_ADDRESS;
  }

  /* malloc space for universal matrix */

  err = matuniv_malloc( matrix, nrow, ncol, type );
  if ( err ) return err;

  /* Register allocated memory with memory list */

  err = memlist_member_register( list, matrix, MEMTYPE_MATUNIV );
  if ( err ) {
    MUTIL_FREE_WARN( matuniv, matrix );
    return err;
  }

  MUTIL_TRACE( "Done with matuniv_malloc_register()" );

  return MUTIL_ERR_OK;
}


/* Free memory list.              */
/* Documented in ut_mem.h         */
/* Written by William Constantine */
mutil_errcode memlist_free( memlist *list )
{
  mutil_errcode  err;
  memlist_node  *temp;

  MUTIL_TRACE( "Start memlist_free()" );

  /* The memory list is terminated with the
     final node pointing to a NULL. */

  while ( list->root != ( memlist_node* ) NULL ) {

    /* validate the current node */

    err = localfn_memlist_member_validate( list->root );

    /* if the node is invalid, issue warning and skip
       over it onto the next node. Otherwise, free the
       memory pointed to by the data member */

    if ( err ) {
      MUTIL_WARN( "Memory node is invalid and the memory pointed "
          "to by the data member cannot be freed" );
    }
    else localfn_memlist_member_data_free( list->root );

    /* make a copy of the pointer to the next node */

    temp = ( list->root )->next;

    /* free the current node */

    MUTIL_FREE_BUFFER_WARN( list->root, sizeof( memlist_node ) );

    /* reset root to be the next node */

    list->root = temp;

    /* decrement list length */

    ( list->length )--;
  }

  MUTIL_TRACE( "Done with memlist_free()" );
  return MUTIL_ERR_OK;
}


/* Delete a node in a memory list,       */
/* and free the associated memory.       */
/* Documented in ut_mem.h                */
/* Written by Keith L. Davidson          */
/* Body rewritten by William Constantine */
mutil_errcode memlist_member_free(
  void    *data,
  memlist *list )
{
  mutil_errcode  err;
  memlist_node  *node;
  memlist_node  *prev;

  MUTIL_TRACE( "Start memlist_member_free()" );

  /* Make sure the address of the list is not NULL */

  if ( list == ( memlist* ) NULL ) {
    MUTIL_ERROR( "Pointer to memory list is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* Make sure the address of the root node is not NULL */

  if ( list->root == ( memlist_node* ) NULL ) {
    MUTIL_ERROR( "Pointer to root node pointer is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* Make sure there is data to be freed */

  if ( data == ( void* ) NULL ) {
    MUTIL_ERROR( "Pointer to data to be freed is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* Make sure root node is valid. No need to check the entire list */

  err = localfn_memlist_member_validate ( list->root );
  if ( err ) {
    MUTIL_ERROR( "Invalid root node of memory list -- no data freed" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( ( list->root )->data == data ) {

    node = list->root;
  }
  else{

    prev = list->root;
    node = prev->next;

    while ( node != NULL ) {

      if ( node->data == data ) {

        prev->next = node->next;
        break;
      }
      else{

        node = node->next;
        prev = prev->next;
      }
    }
  }

  if ( node == NULL ) {

    MUTIL_ERROR( "Data is not registered with "
      "memory manager or list is empty" );
    return MUTIL_ERR_ILLEGAL_ADDRESS;
  }
  else{

    if ( node == list->root ) {
      list->root = node->next;
    }

    /* Free memory associated with the deleted node */

    localfn_memlist_member_data_free( node );

    /* Free memory for the node itself */

    MUTIL_FREE_BUFFER_WARN( node, sizeof( memlist_node ) );

  }

  /* Decrement length of list */

  ( list->length )--;

  MUTIL_TRACE( "Done with memlist_member_free()" );
  return MUTIL_ERR_OK;
}


/* Delete a node in a memory list     */
/* but not the associated memory.     */
/* Documented in ut_mem.h             */
/* William Constantine                */
mutil_errcode memlist_member_unregister(
  void    *data,
  memlist *list )
{
  mutil_errcode  err;
  memlist_node  *node;
  memlist_node  *prev;

  MUTIL_TRACE( "Start memlist_member_unregister()" );

  /* Make sure the address of the list is not NULL */

  if ( list == ( memlist* ) NULL ) {
    MUTIL_ERROR( "Pointer to memory list is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* Make sure the address of the root node is not NULL */

  if ( list->root == ( memlist_node* ) NULL ) {
    MUTIL_ERROR( "Pointer to root node pointer is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* Make sure there is data to be freed */

  if ( data == ( void* ) NULL ) {
    MUTIL_ERROR( "Pointer to data to be freed is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* Make sure root node is valid. No need to check the entire list */

  err = localfn_memlist_member_validate ( list->root );
  if ( err ) {
    MUTIL_ERROR( "Invalid root node of memory list -- no data freed" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( ( list->root )->data == data ) {
    node = list->root;
  }

  else{

    prev = list->root;
    node = prev->next;

    while ( node != NULL ) {

      if ( node->data == data ) {

        prev->next = node->next;
        break;
      }
      else{

        node = node->next;
        prev = prev->next;
      }
    }
  }

  if ( node == NULL ) {

    MUTIL_ERROR( "Data is not registered with "
      "memory manager or list is empty" );
    return MUTIL_ERR_ILLEGAL_ADDRESS;
  }

  else{

    if ( node == list->root ) {
      list->root = node->next;
    }

    /* Free memory for the node itself */

    MUTIL_FREE_BUFFER_WARN( node, sizeof( memlist_node ) );
  }

  /* Decrement length of list */

  ( list->length )--;

  MUTIL_TRACE( "Done with memlist_member_unregister()" );
  return MUTIL_ERR_OK;
}


/* Prepend new member to a memory list. */
/* Documented in ut_mem.h               */
/* Written by William Constantine       */
mutil_errcode memlist_member_register(
  memlist      *list,
  void         *new_data,
  memlist_type  mem_type )
{
  mutil_errcode  err;
  memlist_node  *new_node = NULL;

  MUTIL_TRACE( "Start memlist_member_register()" );

  /* validate root node */

  err = localfn_memlist_member_validate( list->root );
  if ( err ) return err;

  /* check to see if memory is already registered */

  if ( localfn_memlist_member_exist_node( new_data, list, &new_node ) ) {
    MUTIL_ERROR( "The data is already registered with the memory manager" );
    return MUTIL_ERR_ILLEGAL_ADDRESS;
  }

  /* allocate memory for new node */

  err = mutil_malloc( sizeof( memlist_node ), ( void** ) &new_node );
  if ( err ) {
    MUTIL_ERROR( "Could not allocate memory for memlist node" );
    return MUTIL_ERR_MEM_ALLOC;
  }

  /* update information in new node */

  new_node->data = new_data;
  new_node->type = mem_type;

  if ( list->root == ( memlist_node* ) NULL ) {

    /* this is the first node in the list:
       (1) point root to the new node and
       (2) terminate the list by pointing
       the "next" member of the new node
       to NULL.                         */

    list->root     = new_node;
    new_node->next = ( memlist_node* ) NULL;
  }

  else {

    /* make the new node's "next" member point
       to the beginning of the current list,
       i.e. point to the location that the
       current root points. this "pushes" the
       new node into the frontmost position
       in the list.                        */

    new_node->next = list->root;
    list->root     = new_node;
  }

  /* increment length member of memory list */

  ( list->length )++;

  MUTIL_TRACE( "Done with memlist_member_register()" );
  return MUTIL_ERR_OK;
}


/* Reallocate universal matrix and  */
/* update memory list               */
/* Documented in ut_mem.h           */
/* Written by William Constantine   */
mutil_errcode matuniv_realloc_register(
  univ_mat *matrix,
  sint32    nrow,
  sint32    ncol,
  memlist  *list )
{
  mutil_errcode  err;
  memlist_node  *node;

  MUTIL_TRACE( "Start matuniv_realloc_register()" );

  if ( localfn_memlist_member_exist_node( matrix, list, &node ) ) {

    /* realloc memory */

    err = matuniv_realloc( matrix, nrow, ncol );
    if ( err ) return err;

    /* replace pointer in data member with new pointer */

    ( node )->data = matrix;
  }
  else {

    MUTIL_ERROR( "The data to be replaced does not exist in memory list" );
    return MUTIL_ERR_ILLEGAL_ADDRESS;
  }

  MUTIL_TRACE( "Done matuniv_realloc_register()" );
  return MUTIL_ERR_OK;
}


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
mutil_errcode matdbl_wrap_data_register(
  double_mat *matrix,
  double     *data,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_WRAP_DATA_LIST( matdbl, matrix, data, nrow, ncol, \
    MEMTYPE_MATDBL, list );
}


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
mutil_errcode matflt_wrap_data_register(
  float_mat *matrix,
  float     *data,
  sint32     nrow,
  sint32     ncol,
  memlist   *list )
{
  TMPL_MAT_WRAP_DATA_LIST( matflt, matrix, data, nrow, ncol, \
    MEMTYPE_MATFLT, list );
}


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
mutil_errcode matu8_wrap_data_register(
  uint8_mat *matrix,
  uint8     *data,
  sint32     nrow,
  sint32     ncol,
  memlist   *list )
{
  TMPL_MAT_WRAP_DATA_LIST( matu8, matrix, data, nrow, ncol, \
    MEMTYPE_MATU8, list );
}


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
mutil_errcode matu16_wrap_data_register(
  uint16_mat *matrix,
  uint16     *data,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_WRAP_DATA_LIST( matu16, matrix, data, nrow, ncol, \
    MEMTYPE_MATU16, list );
}


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
mutil_errcode matu32_wrap_data_register(
  uint32_mat *matrix,
  uint32     *data,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_WRAP_DATA_LIST( matu32, matrix, data, nrow, ncol, \
    MEMTYPE_MATU32, list );
}


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
mutil_errcode mats16_wrap_data_register(
  sint16_mat *matrix,
  sint16     *data,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_WRAP_DATA_LIST( mats16, matrix, data, nrow, ncol, \
    MEMTYPE_MATS16, list );
}


/* This function is documented under matuniv_wrap_data_register in ut_mem.h */
mutil_errcode mats32_wrap_data_register(
  sint32_mat *matrix,
  sint32     *data,
  sint32      nrow,
  sint32      ncol,
  memlist    *list )
{
  TMPL_MAT_WRAP_DATA_LIST( mats32, matrix, data, nrow, ncol, \
    MEMTYPE_MATS32, list );
}


/* Initialize a universal matrix with given size, type, */
/* and data and register it with the memory manager.    */
/* Documented in ut_mem.h                               */
/* Written by William Constantine                       */
mutil_errcode matuniv_wrap_data_register(
  univ_mat        *matrix,
  void            *data,
  sint32           nrow,
  sint32           ncol,
  mutil_data_type  type,
  memlist         *list )
{
  mutil_errcode err;
  memlist_node *node;

  MUTIL_TRACE( "Start matuniv_wrap_data_register()" );

  /* wrap data in universal matrix */

  err = matuniv_wrap_data( matrix, data, nrow, ncol, type );
  if ( err ) return err;

  /* register the new memory if it doesn't already exist */

  REPLACE_OR_ADD_MEMLIST_MEMBER( matrix, MEMTYPE_MATUNIV,
    data, list, node, err );

  MUTIL_TRACE( "Done with matuniv_wrap_data_register()" );

  return MUTIL_ERR_OK;
}


/* Initialize a universal matrix from a non-universal */
/* matrix and register it with the memory manager.    */
/* Documented in ut_mem.h                             */
/* Written by William Constantine                     */
mutil_errcode matuniv_wrap_matrix_register(
  univ_mat        *umat,
  void            *matrix,
  mutil_data_type  type,
  memlist         *list )
{
  mutil_errcode  err;
  memlist_node  *node;

  MUTIL_TRACE( "Start matuniv_wrap_matrix_register()" );

  /* wrap matrix in universal matrix */

  err = matuniv_wrap_matrix( umat, matrix, type );
  if ( err ) return err;

  /* register the new memory if it doesn't already exist */

  REPLACE_OR_ADD_MEMLIST_MEMBER( umat, MEMTYPE_MATUNIV,
                                 matrix, list, node, err );

  MUTIL_TRACE( "Done with matuniv_wrap_matrix_register()" );

  return MUTIL_ERR_OK;
}

/* Initialize a universal matrix from a universal     */
/* matrix and register it with the memory manager.    */
/* Documented in ut_mem.h                             */
/* Written by William Constantine                     */
mutil_errcode matuniv_wrap_univ_matrix_register(
  univ_mat *umat1,
  univ_mat *umat2,
  memlist  *list )
{
  mutil_errcode  err;
  memlist_node  *node;

  MUTIL_TRACE( "Start matuniv_wrap_univ_matrix_register()" );

  /* wrap universal matrix 1 with universal matrix 2 */

  err = matuniv_wrap_univ_matrix( umat1, umat2 );
  if ( err ) return err;

  /* register the new memory if it doesn't already exist */

  REPLACE_OR_ADD_MEMLIST_MEMBER( umat1, MEMTYPE_MATUNIV,
                                 umat2, list, node, err );

  MUTIL_TRACE( "Done with matuniv_wrap_univ_matrix_register()" );

  return MUTIL_ERR_OK;
}


/* Allocate dynamic memory and register with memory manager. */
/* Documented in ut_mem.h                                    */
/* Written by William Constantine                            */
mutil_errcode mutil_malloc_register(
  sint32    size,
  void    **data,
  memlist  *list )
{
  mutil_errcode err;

  MUTIL_TRACE( "Start mutil_malloc_register()" );

  /* Make sure the address of the list is not NULL */

  if ( list == ( memlist* ) NULL ) {
    MUTIL_ERROR( "Pointer to memory list is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* allocate memory */

  err = mutil_malloc( size, data );
  if ( err ) return err;

  /* Register allocated memory with memory list */

  err = memlist_member_register( list, *data, MEMTYPE_BUFFER );
  if ( err ) {
    MUTIL_FREE_BUFFER_WARN( data, size );
    return err;
  }

  /* assign buffer size in corresponding member of current node */

  ( list->root )->buffer_size = size;

  MUTIL_TRACE( "Done with mutil_malloc_register()" );

  return MUTIL_ERR_OK;
}


/* Allocate dynamic memory, initialized to zero and register */
/* the memory with the memory manager.                       */
/* Documented in ut_mem.h                                    */
/* Written by William Constantine                            */
mutil_errcode mutil_calloc_register(
  sint32    nobj,
  size_t    size,
  void    **data,
  memlist  *list )
{
  mutil_errcode err;

  MUTIL_TRACE( "Start mutil_calloc_register()" );

  /* Make sure the address of the list is not NULL */

  if ( list == ( memlist* ) NULL ) {
    MUTIL_ERROR( "Pointer to memory list is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* allocate memory */

  err = mutil_calloc( nobj, size, data );
  if ( err ) return err;

  /* Register allocated memory with memory list */

  err = memlist_member_register( list, *data, MEMTYPE_BUFFER );
  if ( err ) {
    MUTIL_FREE_BUFFER_WARN( data, ( size * nobj ) );
    return err;
  }

  /* assign buffer size in corresponding member of current node */

  ( list->root )->buffer_size = (sint32) ( size * nobj );

  MUTIL_TRACE( "Done with mutil_calloc_register()" );

  return MUTIL_ERR_OK;
}


/* Change the size of dynamic memory and register the memory */
/* with the memory manager.                                  */
/* Documented in ut_mem.h                                    */
/* Written by William Constantine                            */
mutil_errcode mutil_realloc_register(
  void    **data,
  sint32    new_size,
  sint32    old_size,
  memlist  *list )
{
  mutil_errcode  err;
  memlist_node  *node;

  MUTIL_TRACE( "Start mutil_realloc_register()" );

  /* Make sure the address of the list is not NULL */

  if ( list == ( memlist* ) NULL ) {
    MUTIL_ERROR( "Pointer to memory list is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* first locate the appropriate node in the memory list */

  if ( localfn_memlist_member_exist_node( *data, list, &node ) ) {

    /* realloc memory */

    err = mutil_realloc( data, new_size, old_size );
    if ( err ) return err;

    /* replace pointer in data member with new pointer */

    ( node )->data = *data;

    /* update buffer size */

    ( node )->buffer_size = new_size;
  }

  else{
    MUTIL_ERROR( "The data to be replaced does not exist in memory list" );
    return MUTIL_ERR_ILLEGAL_ADDRESS;
  }

  MUTIL_TRACE( "Done with mutil_realloc_register()" );

  return MUTIL_ERR_OK;
}


/* Check the existence of pre-registered memory in the */
/* memory manager.                                     */
/* Documented in ut_mem.h                              */
/* Written by William Constantine                      */
boolean memlist_member_exist(
  void    *data,
  memlist *list )
{
  memlist_node  *node;
  boolean       found = FALSE;

  MUTIL_TRACE( "Start memlist_member_exist()" );

  /* find the address in the list */

  for ( node = list->root; node != (memlist_node *) NULL; node = node->next ) {

    if ( ( node )->data == data ) {
      found = TRUE;
      break;
    }
  }

  MUTIL_TRACE( "Done with memlist_member_exist()" );
  return found;
}


/* Validate a memory list.      */
/* Documented in ut_mem.h       */
/* Written by Keith L. Davidson */
mutil_errcode memlist_validate( memlist *list )
{
  memlist_node*   curr_node;
  sint32          count;
  mutil_errcode   err;

  MUTIL_TRACE( "Start memlist_validate()" );

  /* Check for NULL pointer */

  if ( list == ( memlist* ) NULL ) {
    MUTIL_ERROR( "NULL memlist pointer passed to function" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* Check length of list */

  if ( list->length < 1 ) {

    /* make sure the list isn't empty (newly initialized) */

    if ( !( ( list->length == 0 ) &&
      ( list->root == ( memlist_node* ) NULL ) ) ) {
      MUTIL_ERROR( "Length member of list structure has value less "
        "than 1,\n or zero length member and non-NULL root node" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
  }

  curr_node = list->root;
  count = 0;
  while( ( curr_node != ( memlist_node* ) NULL ) &&
         ( count < MUTIL_SINT32_MAX ) ) {
    count++;
    if ( ( count == list->length ) &&
      ( curr_node->next != ( memlist_node* ) NULL ) ) {
      MUTIL_ERROR( "Number of list members larger than list length" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
    curr_node = curr_node->next;
  }

  if ( count != list->length ) {
    if ( count == MUTIL_SINT32_MAX ) {
      MUTIL_ERROR( "Non-terminating memory list" );
    }
    else {
      MUTIL_ERROR( "Number of list members smaller than list length" );
    }
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* Since the list length is valid, we now validate the
     memory of each list member */

  curr_node = list->root;
  for ( count = 0; count < list->length; count++ ) {
    err = localfn_memlist_member_validate( curr_node );
    if ( err ) {
      MUTIL_ERROR( "Memory list has invalid member" );
      return err;
    }
    curr_node = curr_node->next;
  }

  MUTIL_TRACE( "Done with memlist_validate()" );
  return MUTIL_ERR_OK;
}

/* Print a report of the contents */
/* of a memory list               */
/* Documented in ut_mem.h         */
/* Written by William Constantine */
mutil_errcode memlist_print( FILE *out_file, memlist *list )
{
  mutil_errcode  err;
  memlist_node  *temp;
  sint32         count = 0;
  sint32         type;
  void          *data;
  sint32         nrow;
  sint32         ncol;

  char *memtype[]={
    (char *) "double matrix",
    (char *) "float matrix",
    (char *) "uint8 matrix",
    (char *) "uint16 matrix",
    (char *) "uint32 matrix",
    (char *) "sint16 matrix",
    (char *) "sint32 matrix",
    (char *) "complex matrix",
    (char *) "universal matrix",
    (char *) "matrix set",
    (char *) "kd-tree structure",
    (char *) "contiguous buffer" };


  char *umtype[]={
    (char *) "uint8 matrix",
    (char *) "sint8 matrix",
    (char *) "uint16 matrix",
    (char *) "sint16 matrix",
    (char *) "uint32 matrix",
    (char *) "sint32 matrix",
    (char *) "float matrix",
    (char *) "double matrix",
    (char *) "complex matrix"
  };


  MUTIL_TRACE( "Start memlist_report()" );

  /* validate the root node */

  err = localfn_memlist_member_validate( list->root );
  if ( err ){
    MUTIL_ERROR( "Root node in memory list is invalid. A report cannot be generated." );
    return ( err );
  }

  temp = list->root;

  (void) fprintf( out_file, "\nNODE\tMEMBER ADDRESS\tMEMBER TYPE\n" );
  (void) fprintf( out_file, "----\t--------------\t-------------------------------\n" );

  /* The memory list is terminated with the
     final node pointing to a NULL. */

  while ( temp != ( memlist_node* ) NULL ) {

    type = (sint32) ( temp->type );

    switch( type ){

    case 0:
      data = (void *) ( (double_mat *) ( temp->data ) )->data;
      nrow = (sint32) ( (double_mat *) ( temp->data ) )->nrow;
      ncol = (sint32) ( (double_mat *) ( temp->data ) )->ncol;
      break;
    case 1:
      data = (void *) ( (float_mat *) ( temp->data ) )->data;
      nrow = (sint32) ( (float_mat *) ( temp->data ) )->nrow;
      ncol = (sint32) ( (float_mat *) ( temp->data ) )->ncol;
      break;
    case 2:
      data = (void *) ( (uint8_mat *) ( temp->data ) )->data;
      nrow = (sint32) ( (uint8_mat *) ( temp->data ) )->nrow;
      ncol = (sint32) ( (uint8_mat *) ( temp->data ) )->ncol;
      break;
    case 3:
      data = (void *) ( (uint16_mat *) ( temp->data ) )->data;
      nrow = (sint32) ( (uint16_mat *) ( temp->data ) )->nrow;
      ncol = (sint32) ( (uint16_mat *) ( temp->data ) )->ncol;
      break;
    case 4:
      data = (void *) ( (uint32_mat *) ( temp->data ) )->data;
      nrow = (sint32) ( (uint32_mat *) ( temp->data ) )->nrow;
      ncol = (sint32) ( (uint32_mat *) ( temp->data ) )->ncol;
      break;
    case 5:
      data = (void *) ( (sint16_mat *) ( temp->data ) )->data;
      nrow = (sint32) ( (sint16_mat *) ( temp->data ) )->nrow;
      ncol = (sint32) ( (sint16_mat *) ( temp->data ) )->ncol;
      break;
    case 6:
      data = (void *) ( (sint32_mat *) ( temp->data ) )->data;
      nrow = (sint32) ( (sint32_mat *) ( temp->data ) )->nrow;
      ncol = (sint32) ( (sint32_mat *) ( temp->data ) )->ncol;
      break;
    case 7:
      data = (void *) ( (dcomplex_mat *) ( temp->data ) )->data;
      nrow = (sint32) ( (dcomplex_mat *) ( temp->data ) )->nrow;
      ncol = (sint32) ( (dcomplex_mat *) ( temp->data ) )->ncol;
      break;
    case 8:
      data = MATUNIV_DATA( (univ_mat *) ( temp->data ) );
      nrow = MATUNIV_NROW( (univ_mat *) ( temp->data ) );
      ncol = MATUNIV_NCOL( (univ_mat *) ( temp->data ) );
      break;
    default:
      data = (void *) NULL;
    }

    if ( data ){

      (void) fprintf( out_file, "%ld\t0x%x\t%s (data: 0x%x, dim: %ld x %ld",
		      (long) count++, temp->data, memtype[ type ], data, (long) nrow, (long) ncol );

      if ( type == 8 ){
	(void) fprintf( out_file, ", type: %s )\n", umtype[ (sint32) ( (univ_mat *) ( temp->data ) )->type ] );
      }
      else{
	(void) fprintf( out_file, " )\n" );
      }

    }
    else{

      (void) fprintf( out_file, "%ld\t0x%x\t%s\n",
		      (long) count++, temp->data, memtype[ type ] );
    }

    /* increment to next node */

    temp = temp->next;
  }

  (void) fprintf( out_file, "\n" );

  MUTIL_TRACE( "Done with memlist_free()" );
  return MUTIL_ERR_OK;
}

/**********************************/
/* STATIC FUNCTIONS: DEFINITIONS  */
/**********************************/

/** Check the existence of pre-registered memory in the
 * memory manager and return corresponding node.
 * If found, a boolean (representing the search
 * result) is returned. In addition, the node corresponding to the
 * data found in the memory list is returned through the node
 * variable in the argument list.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #localfn_memlist_member_exist_node( &data, &list, &node );#
 * @return       Boolean. If TRUE, the data has already been registered
 *               in the memory list.
 * @param   data Pointer to data.
 * @param   list Pointer to memory list.
 * @param   node Pointer to pointer of a memory list node. If the
 *               return boolean is TRUE, this argument points to
 *               the node in the memory list which points to the data.
 * @see _memlist_node
 * @private
 */
static boolean localfn_memlist_member_exist_node(
  void          *data,
  memlist       *list,
  memlist_node **node )
{
  boolean found_node = FALSE;

  MUTIL_TRACE( "Start localfn_memlist_member_exist_node()" );

  /* find the address in the list */

  for ( *node = list->root; *node != ( memlist_node* ) NULL;
  ( *node ) = ( *node )->next ) {

    if ( ( *node )->data == data ) {
      found_node = TRUE;
      break;
    }
  }

  MUTIL_TRACE( "Done with localfn_memlist_member_exist_node()" );
  return found_node;
}


/** Memory list node validation.
 * Validates a node in a memory list. A node is considered valid if the
 * data member does not point to NULL and the memory type
 * member is supported. In addition, except in the case the node is of
 * type MEMTYPE\_BUFFER, the function calls the appropriate
 * library function to validate the memory corresponding to the
 * node.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #err_code = localfn_memlist_member_validate( node );#
 * @return         Standard mutils error/OK code.
 * @param node     Pointer to a node in memory list.
 * @see _memlist_node
 * @see _memlist_type
 * @see matuniv_validate
 * @see matset_validate
 * @private
 */
static mutil_errcode localfn_memlist_member_validate(
  memlist_node *node )
{
  mutil_errcode   err;

  MUTIL_TRACE( "Start localfn_memlist_member_validate()" );

  /* if node does not point to NULL then perform checks */

  if ( node != ( memlist_node* ) NULL ) {

    /* test for NULL pointer in data member */

    if ( !( node->data ) ) {

      MUTIL_ERROR( "NULL pointer in memory list node for data member" );
      return MUTIL_ERR_NULL_POINTER;
    }

#define TMPL_MEMLIST_MEMBER_VALIDATE_CALL( FUNCPRE, PCAST )    \
  err = FUNCPRE ## _validate( ( PCAST* ) node->data );         \
  if ( err ) return err

    /* test for supported memory types */

    switch ( node->type ) {
    case MEMTYPE_MATCPX:
      TMPL_MEMLIST_MEMBER_VALIDATE_CALL( matcpx, dcomplex_mat );
      break;

    case MEMTYPE_MATDBL:
      TMPL_MEMLIST_MEMBER_VALIDATE_CALL( matdbl, double_mat );
      break;

    case MEMTYPE_MATFLT:
      TMPL_MEMLIST_MEMBER_VALIDATE_CALL( matflt, float_mat );
      break;

    case MEMTYPE_MATS16:
      TMPL_MEMLIST_MEMBER_VALIDATE_CALL( mats16, sint16_mat );
      break;

    case MEMTYPE_MATS32:
      TMPL_MEMLIST_MEMBER_VALIDATE_CALL( mats32, sint32_mat );
      break;

    case MEMTYPE_MATU16:
      TMPL_MEMLIST_MEMBER_VALIDATE_CALL( matu16, uint16_mat );
      break;

    case MEMTYPE_MATU32:
      TMPL_MEMLIST_MEMBER_VALIDATE_CALL( matu32, uint32_mat );
      break;

    case MEMTYPE_MATU8:
      TMPL_MEMLIST_MEMBER_VALIDATE_CALL( matu8, uint8_mat );
      break;

    case MEMTYPE_MATUNIV:
      TMPL_MEMLIST_MEMBER_VALIDATE_CALL( matuniv, univ_mat );
      break;

    case MEMTYPE_MATSET:
      TMPL_MEMLIST_MEMBER_VALIDATE_CALL( matset, mat_set );
      break;

    case MEMTYPE_KDTREE:
      TMPL_MEMLIST_MEMBER_VALIDATE_CALL( mutil_kdtree, mutil_kdtree );
      break;

    case MEMTYPE_BUFFER:
      break;

    default:
      MUTIL_ERROR( "Memory type is unsupported" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
  }

#undef TMPL_MEMLIST_MEMBER_VALIDATE_CALL

  MUTIL_TRACE( "Done with localfn_memlist_member_validate()" );
  return MUTIL_ERR_OK;
}


/** Free data pointed to by current node of memory list.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_mem.h
 * @source ut\_mem.c
 * @library matrix
 * @usage  #memlist_member_data_free( node );#
 * @return       Nothing (void).
 * @param root  Pointer to current node in memory list.
 * @param data  Pointer to data (any valid type). For matrix sets, memory
 *   will be freed for the universal matrices within the set (if necessary),
 *   in addition to the universal matrix headers.
 * @param type   Enumerated memory type (\Ref{_memlist_type}).
 * @see memlist_free
 * @see MEMLIST_FREE_ON_ERROR
 * @see _memlist_type
 * @see _memlist_node
 * @private
 */
static void localfn_memlist_member_data_free(
  memlist_node *node )
{
  MUTIL_TRACE( "Start localfn_memlist_member_data_free()" );

  switch( node->type ) {

  case MEMTYPE_MATDBL:
    MUTIL_FREE_WARN( matdbl, ( double_mat* ) node->data );
    break;

  case MEMTYPE_MATFLT:
    MUTIL_FREE_WARN( matflt, ( float_mat* ) node->data );
    break;

  case MEMTYPE_MATU8:
    MUTIL_FREE_WARN( matu8, ( uint8_mat* ) node->data );
    break;

  case MEMTYPE_MATU16:
    MUTIL_FREE_WARN( matu16, ( uint16_mat* ) node->data );
    break;

  case MEMTYPE_MATU32:
    MUTIL_FREE_WARN( matu32, ( uint32_mat* ) node->data );
    break;

  case MEMTYPE_MATS16:
    MUTIL_FREE_WARN( mats16, ( sint16_mat* ) node->data );
    break;

  case MEMTYPE_MATS32:
    MUTIL_FREE_WARN( mats32, ( sint32_mat* ) node->data );
    break;

  case MEMTYPE_MATCPX:
    MUTIL_FREE_WARN( matcpx, ( dcomplex_mat* ) node->data );
    break;

  case MEMTYPE_MATUNIV:
    MUTIL_FREE_WARN( matuniv, ( univ_mat* ) node->data );
    break;

  case MEMTYPE_MATSET:
    MUTIL_FREEALL_MATSET_WARN( ( mat_set* ) node->data );
    break;
  case MEMTYPE_KDTREE:
    MUTIL_FREE_WARN( mutil_kdtree, ( mutil_kdtree* ) node->data );
    break;
  case MEMTYPE_BUFFER:
    MUTIL_FREE_BUFFER_WARN( node->data, node->buffer_size );
    break;

  default:
    MUTIL_ERROR( "Memory type not supported" );
    break;
  }

  MUTIL_TRACE( "Done with localfn_memlist_member_data_free()" );
}
