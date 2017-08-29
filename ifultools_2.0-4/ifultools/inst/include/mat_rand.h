
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_rand.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_RAND_H_
#define IN_MAT_RAND_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "fra_type.h"

/* This file contains function declarations
   for the creation of randomized matrices
   and belongs to the mutils matrix library. */

#ifdef __cplusplus
extern "C" {
#endif

/** Random index generation with optional replacement.
 * A set of integers drawn from a unifomrly distributed
 * random distribution is created. The user can specify
 * the draw to be with replacement in which case a given
 * index is allowed to be repeated. The range of returned
 * indices is on [0, N-1] where N is the number of elements
 * in the output matrix (indirectly specified by the user).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_rand.c
 * @library matrix
 * @usage #err = matuniv_random_uniform_indices( nrow, ncol, replacement, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param  nrow The number of rows in the output matrix.
 * @param  ncol The number of columns in the output matrix.
 * @param  replacement A logical flag. If TRUE, the
 *         draw from the uniform random distribution is with
 *         replacement, i.e. indices may be repeated.
 * @param  seed Input value for seeding random number generator.
 *              If the input for the seed is a NULL pointer,
 *              this function will set the seed to a pseudo-randomized number,
 *              based on the date and time of day or some other somewhat
 *              random input, depending on the implementation.
 * @param  intrp_ptr Pointer for implementation of interrupt checking.
 * @param  result Pointer to an nrow by ncol universal matrix of
 *         type MUTIL\_SINT32
 *         which, upon return, will contain the random indices.
 *         This memory for this matrix is automatically
 *         allocated within the function.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_random_uniform_indices( const sint32 nrow, const sint32 ncol, const boolean replacement, void *intrp_ptr, sint32_mat *result );#
 *  \end{itemize}
 * @see mutil_rand_uniform
 * @see matuniv_permute
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_random_uniform_indices(
  const sint32    nrow,
  const sint32    ncol,
  const boolean   replacement,
  void           *rand_ptr,
  void           *intrp_ptr,
  univ_mat       *result );

/* This function is documented under matuniv_random_uniform_indices in mat_rand.h */
MUTIL_LIBEXPORT mutil_errcode mats32_random_uniform_indices(
  const sint32    nrow,
  const sint32    ncol,
  const boolean   replacement,
  void           *rand_ptr,
  void           *intrp_ptr,
  sint32_mat     *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_RAND_H_ */
