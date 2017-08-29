
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_neig.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_NEIG_H_
#define IN_FRA_NEIG_H_


#include "fra_type.h"
#include "mat_type.h"

#include "ut_plat.h"
#include "ut_type.h"
#include "ut_err.h"

#include "ut_kdtr.h"

/* This file contains function declarations for nearest neighbor searches   */
/* in the MUTILS fractal library. Code from the S+Spatial module was        */
/* adapted for use here in MUTILS.                                          */

#ifdef __cplusplus
extern "C" {
#endif

/** Find the nearest neighbors of a multidimensional embedding
 * of data in the phase space.
 *
 * Finds a user-specified number of nearest neighbors of a multivariate space
 * defined by the coordinates of the input matrix whose columns are the
 * coordinates. Alternatively, the user can specify a maximum radius over
 * which to search for neighbors. An efficient recursive algorithm is used to
 * find all nearest neighbors. First, a kd-tree created from the input is
 * traversed to find the leaf with medians nearest the point for which
 * neighbors are sought. Second, observations in this leaf are searched to
 * find nearest neighbors. Finally, if necessary, adjoining leaves are
 * searched for nearest neighbors.
 *
 * References:
 *
 * 1. Friedman, J., Bentley, J. L., and Finkel, R. A. (1977). An algorithm
 *    for finding best matches in logarithmic expected time. ACM
 *    Transactions on Mathematical Software 3, 209-226.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_neig.h
 * @source fra\_neig.c
 * @library fractal
 * @usage #err = frauniv_neighbor_find(&embed,nn,dmax,FRA_DISTANCE_L1,NULL,TRUE,orb_lag,intrp_ptr,&original_indices,&neighbor_indices,&neighbor_distances);#
 * @return Standard mutils error/OK code.
 * @param embedding       Pointer to a pre-allocated universal matrix
 *                        of type MUTIL\_DOUBLE containing the embedding
 *                        matrix. Each column of the matrix respresents
 *                        a single coordinate of the embedding and each
 *                        row denotes the coordiantes of a single point
 *                        in the embedding.
 * @param n_neighbor      The number of neighbors to find for each point in
 *                        the embedding. If you wish to use distance\_max
 *                        instead, you must set n\_neighbor <= 0.
 * @param distance_max    Used an alternative to n\_neighbor, use this
 *                        parameter to specify the maximum distance to
 *                        search relative to the current point in the
 *                        phase space. Set  n\_neighbor <= 0 to use this
 *                        option.
 * @param distance_metric The metric used to define the distance between
 *                        points in the embedding. The argument is of
 *                        type \Ref{_fra_distance_metric}.
 * @param pdmatrix        Pointer to a pre-allocated universal matrix
 *                        of type MUTIL\_DOUBLE containing a square positive
 *                        definite matrix. The number of columns (rows) must
 *                        match the column dimension of {\bf points}. This
 *                        input is ignored if {\bf distance_metric} is
 *                        different from FRA\_DISTANCE\_MAHALANOBIS. Otherwise,
 *                        the matrix will be used to implement the Mahalanobis
 *                        distance metric during neighbor searching. NOTE: The
 *                        matrix is not checked for positive definiteness.
 * @param sort_distances  If TRUE then the neighbor distances will be sorted
 *                        in ascending order, and the corresponding neighbor
 *                        indices are arranged in accordance with the sort.
 * @param orbital_lag     Time (in lags) progresses along the rows of the
 *                        input embedding matrix. This parameter is used to
 *                        find neighbers close in the phase space but not
 *                        close in time (lag). If greater than zero then
 *                        no returned neighbors of a particular point will
 *                        be less than obrbital_delay lags from that point.
 *                        E.g., setting this parameter to 1 prevents a point
 *                        from being its own neighbor.
 * @param intrp_ptr       Pointer for implementation of interrupt checking.
 * @param original_indices Pointer to a single-column matrix of type
 *                        MUTIL\_SINT32 which (upon return) will contain the
 *                        indices (rows) of the points in the input matrix for
 *                        which neighbors were found. The memory for this
 *                        matrix is allocated within the function.
 * @param neighbor_indices Pointer to a single-column matrix of type
 *                        MUTIL\_SINT32 which (upon return) will contain the
 *                        indices of the neighbors found for those listed in
 *                        the original\_indices vector. The memory for this
 *                        matrix is allocated within the function.
 * @param neighbor_distances Pointer to a single-column matrix of type
 *                        MUTIL\_DOUBLE which (upon return) will contain the
 *                        distances between the points denoted by the
 *                        original\_indices vectors and those denoted by the
 *                        correspoding neighbor\_indices vector.
 *                        The memory for this matrix is allocated within the
 *                        function. The distance is based on the metric
 *                        specified by distance\_metric.
 *
 * @see _fra_distance_metric
 * @see frauniv_embed_neighbors
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_neighbor_find(
  const univ_mat            *embedding,
  const sint32               n_neighbor,
  const double               distance_max,
  const fra_distance_metric  distance_metric,
  const univ_mat            *pdmatrix,
  const boolean              sort_distances,
  const sint32               orbital_lag,
  void                      *intrp_ptr,
  univ_mat                  *original_indices,
  univ_mat                  *neighbor_indices,
  univ_mat                  *neighbor_distances );


/** Find the nearest neighbors of an arbitrary set of points, from an embedding
 *in the phase space.
 *
 * For a given set of points in the phase space this function finds nearest
 * neighbors for each point, drawing neighbors from an embedding in the phase
 * space. This function is more general than
 * \Ref{frauniv_neighbor_find}, where neighbors of each point in an
 * embedding (or any set of points in the phase space)
 * are drawn from the embedding itself. Here the points for which neighbors
 * are found may be different than the set of points the neighbors are
 * drawn from. This neighbor set is represented by a kd-tree, which is
 * typically the result of a call to \Ref{mutil_kdtree_malloc}.
 *
 * The function return three vectors, all with the same number of elements.
 * The first, {\bf original_indices}, contains indices of
 * points in {\bf points} for which neighbors were found. The number of times
 * an index is repeated is equal to the number of neighbors found for the
 * corresponding point in {\bf points}. The indices are grouped contiguously
 * and the groups are ordered according to {\bf points}. Thus, if the index
 * the first row of points in no the first element of {\bf original_indices}
 * then no neighbors were found for that point. The second output vector,
 * {\bf neigbor_indices}, contains indices of the neighbors found in
 * the kd-tree. NOTE: these indices are the row numbers of the embedding
 * matrix that was passed to \Ref{mutil_kdtree_malloc} to create the kd-tree.
 * Zero-based indexing is used. The third output vector,
 * {\bf neighbor_distances}, is the distance between points and corresponding
 * neighbors. As an example, using array notation for convenience, suppose
 * {\bf original_indices[i]}=3 and {\bf neighbor_indices[i]}=7. Then the
 * point given by the fourth row of
 * {\bf points} has the eighth point of {\bf kdtree} as a neighbor and the
 * distance between these two points is neighbor_distances[i]. The metric
 * used is determined by {\bf distance_metric}. This example assumes the
 * "natural" indexing of the rows of {bf points}, i.e., 0, 1, 2, ...
 *
 * In the special case where {\bf points} is an arbitrary subset of the
 * embedding used to create {\bf kdtree}, it may be desirable to use
 * the original indices instead of assuming the "natural" indexing scheme.
 * The input {\bf timestamp} may then be used to pass these original indices.
 * {\bf timestamp} is a row or colum vector with length equal to the number
 * of rows in {\bf points}. This vector contains the (zero-based) original
 * indices of {\bf points}.
 *
 * Reference:
 *
 *    Friedman, J., Bentley, J. L., and Finkel, R. A. (1977). An algorithm
 *    for finding best matches in logarithmic expected time. ACM
 *    Transactions on Mathematical Software 3, 209-226.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_neig.h
 * @source fra\_neig.c
 * @library fractal
 * @usage #err = frauniv_neighbor_find_arbitrary(&embed,&kdt,n_neighbor,distance_max,FRA_DISTANCE_L2,NULL,false,0,NULL,intrp_ptr,&original_indices,&neighbor_indices,&neighbor_distances);#
 * @return Standard mutils error/OK code.
 * @param points          Pointer to a pre-allocated universal matrix
 *                        of type MUTIL\_DOUBLE containing points in the
 *                        phase space. Each column of the matrix
 *                        respresents a single coordinate and each
 *                        row denotes the coordiantes of a single point.
 *                        The matrix may have a single row, i.e., represent
 *                        a single point in the phase space.
 * @param kdtree          A \Ref{_mutil_kdtree}, typically resulting from a
 *                        call to \Ref{mutil_kdtree_malloc}, representing an
 *                        embedding in the phase space. The
 *                        dimension of the phase space must be identical to
 *                        the number of columns in {\bf points}.
 * @param n_neighbor      The number of neighbors to find for each point.
 *                        If you wish to use {\bf distance_max}
 *                        instead, to search for neighbors within a given
 *                        distance, you must set {\bf n_neighbor} to zero.
 * @param distance_max    For each point, tells the function to find all
 *                        neighbors within this distance of the point. Each
 *                        point in {\bf points} will be the center of a "ball"
 *                        in the phase space with radius {\bf distance_max}.
 *                        The neighbors of that point will be all points in
 *                        {\bf kdtree} which are contained in this ball.
 *                        This parameter is ignored if {\bf n_neighbor} is
 *                        non-zero.
 * @param distance_metric The metric used to define the distance between
 *                        points in the phase space. See
 *                        \Ref{_fra_distance_metric} for supported metrics.
 * @param pdmatrix        Pointer to a pre-allocated universal matrix
 *                        of type MUTIL\_DOUBLE containing a square positive
 *                        definite matrix. The number of columns (rows) must
 *                        match the column dimension of {\bf points}. This
 *                        input is ignored if {\bf distance_metric} is
 *                        different from FRA\_DISTANCE\_MAHALANOBIS. Otherwise,
 *                        the matrix will be used to implement the Mahalanobis
 *                        distance metric during neighbor searching. NOTE: The
 *                        matrix is not checked for positive definiteness.
 * @param sort_distances  If TRUE, the neighbor distances will be sorted
 *                        in ascending order, and the corresponding neighbor
 *                        indices permutated accordingly.
 * @param orbital_lag     NOTE: This parameter is only meaningful if
 *                        {\bf points} is a subset of the points used to
 *                        create {\bf kdtree}. Time (in lags) progresses
 *                        along the rows of {\bf points}. This parameter is
 *                        used to find neighbers spatially close in the phase
 *                        space but not close in time (lag). No returned
 *                        neighbors of a particular point will be less than
 *                        {\bf orbital_lag} lags from that point. Setting
 *                        this parameter to 1 prevents a point from being its
 *                        own neighbor; setting it to 0 allows for neighbors
 *                        close in space {\it and} time.
 * @param timestamp       NOTE: This parameter is only meaningful if
 *                        {\bf points} is a subset of the points used to
 *                        create {\bf kdtree}. This parameter is a column
 *                        or row vector with length equal to the number of
 *                        rows {\bf points}. The elements of timestamp
 *                        are the (zero-based) indices of the points in
 *                        {\bf points}, with respect to the points represented
 *                        by {\bf kdtree}.
 * @param intrp_ptr       Pointer for implementation of interrupt checking.
 * @param original_indices Pointer to a single-column matrix of type
 *                        MUTIL\_SINT32 which will contain the indices
 *                        of the points in {\bf points} for which neighbors
 *                        were found. See detailed description above. Memory
 *                        for this matrix is allocated within the function.
 * @param neighbor_indices Pointer to a single-column matrix of type
 *                        MUTIL\_SINT32 which will contain the neighbor
 *                        indices for the indices in {\bf origianl_indices}.
 *                        See detailed description above. Memory
 *                        for this matrix is allocated within the
 *                        function.
 * @param neighbor_distances Pointer to a single-column matrix of type
 *                        MUTIL\_DOUBLE which will contain the
 *                        distances between the points in {\bf points} and
 *                        correspoding neighbors. The distance is based on
 *                        the metric specified by {\bf distance_metric}.See
 *                        description above. Memory for this matrix is
 *                        allocated within the function.
 *
 * @see frauniv_embed_neighbors
 * @see _fra_distance_metric
 * @see _mutil_kdtree
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_neighbor_find_arbitrary(
  const univ_mat            *points,
  const mutil_kdtree        *kdtree,
  const sint32               n_neighbor,
  const double               distance_max,
  const fra_distance_metric  distance_metric,
  const univ_mat            *pdmatrix,
  const boolean              sort_distances,
  const sint32               orbital_lag,
  const univ_mat            *timestamp,
  void                      *intrp_ptr,
  univ_mat                  *original_indices,
  univ_mat                  *neighbor_indices,
  univ_mat                  *neighbor_distances );



#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_NEIG_H_ */
