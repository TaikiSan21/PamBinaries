
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_pca.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_PCA_H_
#define IN_MAT_PCA_H_

#include <math.h>

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

/* This file contains function declarations for matrix
   principal component analysis functions.
*/

#ifdef __cplusplus
extern "C" {
#endif


/* bill - would you go over this and explain?
   why is it htat at least two matrices must be input?
   why won't this work on one matrix?  it pca done
   across one multiband pixel? of on each matrix independently */

/* bill: also what is all of this scaling business? shouldn't we
   return doubles, and then let the user decide what to do? */

/** Principal component analysis (PCA)
 * Function to transform a matrix set to the
 * set of principal components values. A fast PCA technique is used.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_pca.h
 * @source mat\_pca.c
 * @library matrix
 * @usage  #err_code = matset_principal_components_analysis( &input, &intrp_ptr, &output );#
 * @return             Standard mutils error/OK code.
 * @param in_matrices  Pointer to multiband matrix set.
 *                     Each band is represented by one matrix.
 * @param intrp_ptr    Pointer for implementation of interrupt checking.
 * @param pca_matrices Pointer to the matrix set that contains
 *                     principal components.  The data structure should be
 *                     created outside this function.
 *                     The type of pca\_matrices should be double.
 *                     The number of matrices should be the same as the number of bands.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_principal_components_analysis(const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_principal_components_analysis(const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_principal_components_analysis(const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_principal_components_analysis(const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_principal_components_analysis(const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_principal_components_analysis(const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_principal_components_analysis(const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices);#
 * \end{itemize}
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matset_principal_components_analysis(
  const mat_set *in_matrices,  void *intrp_ptr,  mat_set *pca_matrices );


/* function documented with matrix set above */
MUTIL_LIBEXPORT mutil_errcode matu8_principal_components_analysis(
  const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices );


/* function documented with matrix set above */
MUTIL_LIBEXPORT mutil_errcode matu16_principal_components_analysis(
  const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices );


/* function documented with matrix set above */
MUTIL_LIBEXPORT mutil_errcode matu32_principal_components_analysis(
  const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices );


/* function documented with matrix set above */
MUTIL_LIBEXPORT mutil_errcode mats16_principal_components_analysis(
  const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices );


/* function documented with matrix set above */
MUTIL_LIBEXPORT mutil_errcode mats32_principal_components_analysis(
  const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices );


/* function documented with matrix set above */
MUTIL_LIBEXPORT mutil_errcode matflt_principal_components_analysis(
  const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices );


/* function documented with matrix set above */
MUTIL_LIBEXPORT mutil_errcode matdbl_principal_components_analysis(
  const mat_set *in_matrices, void *intrp_ptr, mat_set *pca_matrices );


#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_PCA_H_*/
