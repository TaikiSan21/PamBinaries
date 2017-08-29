
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_util.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_UTIL_H_
#define IN_FRA_UTIL_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

/* This file contains function declarations the mutils fractal library. */

#ifdef __cplusplus
extern "C" {
#endif

/** Check delay coordinate embedding arguments.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_util.h
 * @source fra\_util.c
 * @library fractal
 * @usage #err = frauniv_check_embedding_inputs( &data, dim, time_lag, orbital_lag, intrp_ptr, &is_delay_embedding, &n_embed );#
 * @return Standard mutils error/OK code.
 * @param  data    Pointer to either (1) a pre-allocated single-row or
 *                 single-column universal matrix of type
 *                 MUTIL\_DOUBLE containing the time series, or (2)
 *                 an embedding matrix whose columns contain the data
 *                 for individual coordinates.
 * @param  dim         An integer denoting the embedding dimension.
 * @param  time_lag  An integer denoting the integer delay to test
 *                     for a delay embedding. This is only tested if
 *                     the input matrix (data) is a vector.
 * @param orbital_lag  An integer denoting the lag between points
 *                     along a trajectory/orbit.
 * @param intrp_ptr    Pointer for implementation of interrupt checking.
 * @param is_delay_embedding Pointer to a boolean denoting whether the
 *                     embedding is a delay embedding or a regular embedding.
 * @param n_embed      Pointer to an integer denoting the number of points
 *                     in the embedding.
 * @see frauniv_neighbor_partition
 * @see frauniv_embed_neighbors
 * @private
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_check_embedding_inputs(
  const univ_mat     *data,
  sint32              dim,
  sint32              time_lag,
  sint32              orbital_lag,
  void               *intrp_ptr,
  boolean            *is_delay_embedding,
  sint32             *n_embed );

/** Check 2-D partition histogram argument.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_util.h
 * @source fra\_util.c
 * @library fractal
 * @usage #err = frauniv_check_partition( &partition );#
 * @return Standard mutils error/OK code.
 * @param  partition      Pointer to a pre-allocated matrix set
 *                        containing two single-row or
 *                        single-column universal matrices of type
 *                        MUTIL\_SINT32. The matrices contain, respectively,
 *                        the flattened 2-D phase space partition
 *                        and a vector of indices corresponding to
 *                        the points in the delayed coordinate phase space which fall
 *                        into the 2-D partitioned space. The partitioned
 *                        space `matrix' contains integers which
 *                        point to the location within the index vector
 *                        where the collection of points corresponding
 *                        the current box are located.
 * @param  scale          The scale of the partition.
 * @see frauniv_neighbor_partition
 * @see frauniv_embed_neighbors
 * @private
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_check_partition(
  const mat_set *partition,
  double         scale );

#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_UTIL_H_*/









