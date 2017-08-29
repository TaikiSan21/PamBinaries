
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/str_type.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_STR_TYPE_H_
#define IN_STR_TYPE_H_

#include "mat_type.h"

/* This file contains typedefs and structs for the non-matrix
   algorithmic-based data-structure types and memory management
   for the mutils library.
*/

#ifdef __cplusplus
extern "C" {
#endif


/** Enum of memory allocation types.
 * The allocated memory types are used to identify the
 * type of data pointed to by the data member of the
 * \Ref{_memlist_node} structures.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include str\_type.h
 * @source str\_type.h
 * @library matrix
 * @same #typedef enum _memlist_type#
 */
enum _memlist_type
{
  /** Double matrix */
  MEMTYPE_MATDBL,

  /** Float matrix */
  MEMTYPE_MATFLT,

  /** Unsigned integer 8-bit matrix */
  MEMTYPE_MATU8,

  /** Unsigned integer 16-bit matrix */
  MEMTYPE_MATU16,

  /** Unsigned integer 32-bit matrix */
  MEMTYPE_MATU32,

  /** Signed integer 16-bit matrix */
  MEMTYPE_MATS16,

  /** Signed integer 32-bit matrix */
  MEMTYPE_MATS32,

  /** Complex matrix */
  MEMTYPE_MATCPX,

  /** Universal matrix */
  MEMTYPE_MATUNIV,

  /** Matrix set */
  MEMTYPE_MATSET,

  /** A kd-tree structure */
  MEMTYPE_KDTREE,

  /** Contiguous memory buffer */
  MEMTYPE_BUFFER
};

/* See above documentation for _memlist_type for explanation. */
typedef enum _memlist_type memlist_type;


/** Structure for memory node in linked list.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include str\_type.h
 * @source str\_type.h
 * @library matrix
 * @same
 *   #typedef struct _memlist_node memlist_node;#
 * @see _memlist
 * @see _memlist_type
 */
struct _memlist_node {

    /** pointer to data (any type) */
    void *data;

    /** memory type */
    memlist_type type;

    /** buffer size */
    sint32 buffer_size;

    /** pointer to the next node in the linked memory list */
    struct _memlist_node *next;
};


/* See above documentation for _memlist_node */
typedef struct _memlist_node memlist_node;


/** Structure for memory manager.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include str\_type.h
 * @source str\_type.h
 * @library matrix
 * @same
 *   #typedef struct _memlist memlist;#
 * @see _memlist_node
 */
struct _memlist {

    /** pointer to root node */
    memlist_node *root;

    /** number of nodes in memory list */
    sint32 length;
};


/* See above documentation for _memlist */
typedef struct _memlist memlist;


/** Structure holding a kd-tree which represents a set of data points in
 *  k-dimensional space.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_kdtr.h
 * @source ut\_kdtr.h
 * @library matrix
 * @same
 *   #typedef struct _mutil_kdtree mutil_kdtree;#
 *
 */
struct _mutil_kdtree {

  /** Pointer to the data represented by the kd-tree */

  double_mat     points;

  /** Number of points represented by the tree */

  sint32         npoints;

  /** Dimension of the space from which the represented points are drawn */

  sint32         dim;

  /** Bucket size */

  sint32         bucket_size;

  /** Split indices */

  sint32_mat     split_index;

  /** Point index */

  sint32_mat     point_index;

  /** Computed medians used to traverse the kd-tree */

  double_mat     medians;

  /** Memory management list */

  memlist        mlist;
};

/* See above documentation for _mutil_kdtree for explanation. */
typedef struct _mutil_kdtree mutil_kdtree;


#ifdef __cplusplus
}
#endif

#endif /* IN_STR_TYPE_H_ */
