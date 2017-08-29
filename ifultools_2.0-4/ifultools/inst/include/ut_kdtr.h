
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_kdtr.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_KDTR_H_
#define IN_UT_KDTR_H_


/* This file contains data types, function prototypes and macros which   */
/* provide support for a kd-tree. A kd-tree is a data structure that may */
/* be used for fast nearest-neighbor searches.                           */


#include "str_type.h"

#ifdef __cplusplus
extern "C" {
#endif


/** Create a kd-tree representing a set of points in k-dimensional space.
 *
 * This function creates a kd-tree from a given set of points in k-dimensional
 * Euclidian space. A kd-tree is a data structure containing a partition of a
 * set of data points that may be used for fast nearest-neighbor searches.
 * See \Ref{_mutil_kdtree} for further details.
 *
 * The resulting \Ref{_mutil_kdtree} structure will
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_kdtr.h
 * @source ut\_kdtr.c
 * @library matrix
 * @usage #err = mutil_kdtree_malloc(&pts,1,intrp_ptr,&kd_struct);#
 * @return Standard mutils error/OK code.
 * @param points          Pointer to a universal matrix
 *                        of type MUTIL\_DOUBLE containing the data
 *                        to be represented.  Each column of the
 *                        matrix respresents a single coordinate and each
 *                        row denotes the coordiantes of a single point.
 * @param bucket_size     A positive integer specifying the bucket size of
 *                        the kd-tree. This gives the maximum number of
 *                        points in any leaf of the tree.
 * @param kdt             Pointer to a structure of type \Ref{_mutil_kdtree}
 *                        which contains attributes that define the kd-tree.
 *                        All memory will be allocated by this function.
 *                        See \Ref{_mutil_kdtree} for more details on
 *                        kd-trees.
 *
 * @see _mutil_kdtree
 */
MUTIL_LIBEXPORT mutil_errcode mutil_kdtree_malloc(
  mutil_kdtree     *kdt,
  const univ_mat   *points,
  const sint32      bucket_size );


/** Free allocated memory within a kd-tree structure.
 *
 * This function frees the memory allocated for the fields of the
 * \Ref{_mutil_kdtree} structure.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_kdtr.h
 * @source ut_kdtr.c
 * @library matrix
 * @usage #err = mutil_kdtree_free(&kd_struct);#
 * @return Standard mutils error/OK code.
 * @param kdt       Pointer to a \Ref{_mutil_kdtree} structure for
 *                  which memory is to be freed.
 *
 * @see _mutil_kdtree
 * @see mutil_kdtree_malloc
 */
MUTIL_LIBEXPORT mutil_errcode mutil_kdtree_free( mutil_kdtree *kdt );


/** Test the validity of a kd-tree structure.
 *
 * This functions tests the validity of a kd-tree structure. A structure is
 * considered valid if the internal matrices are valid, the dimensions are
 * consistent, and all memory is registered with it's internal memory
 * management list.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_kdtr.h
 * @source ut\_kdtr.c
 * @library matrix
 * @usage #err = mutil_kdtree_validate(&kd_struct);#
 * @return Standard mutils error/OK code.
 * @param kdt       Pointer to a \Ref{_mutil_kdtree} structure. The
 *                  structure will be considered valid if its internal
 *                  matrices and memory list are successfully validated,
 *                  and other fields have consistent values.
 *
 * @see _mutil_kdtree
 */
MUTIL_LIBEXPORT mutil_errcode mutil_kdtree_validate(
  const mutil_kdtree *kdt );


#ifdef __cplusplus
 }
#endif

#endif /* ifdef IN_UT_KDTR_H_ */
