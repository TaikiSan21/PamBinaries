
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_kde.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_KDE_H_
#define IN_FRA_KDE_H_


#include "mat_type.h"


/* This file contains function declarations for computing */
/* kernel density estimates */

#ifdef __cplusplus
extern "C" {
#endif


/** Non-parametric kernel density estimation.
 * Non-parametric estimation of a multi-dimensional probability density
 * function at a given set of points. This function performs kernel density
 * estimation using the Epanechnikov kernel to compute the density
 * estimate.
 *
 * Input {\bf data} is an $N \times D$ universal matrix of (real) points in
 * $D$-dimensional Euclidian space. These points are used to compute kernel
 * estimates. Input {\bf points} is a $M \times D$ universal matrix of points
 * at which the kernel will be estimated. {\bf points} may be NULL, in which
 * case the kernel will be estimated at the points in {\bf data}.
 *
 * The kernel bandwidth is determined by first computing the minimum variance
 * of all dimensions of {\bf data} (over all columns). This minimum variance
 * is then used in Scott's Rule to compute the final bandwidth.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_kde.h
 * @source fra\_kde.c
 * @library fractal
 * @usage #err = frauniv_kernel_density_estimate(&data,&points,intrp_ptr,&kde);#
 * @return Standard mutils error/OK code.
 *
 * @param data        Pointer to a universal matrix of type MUTIL\_DOUBLE
 *                    which holds the points from which the estimate will
 *                    be computed. The dimension of the space is determined
 *                    by the number of columns.
 *
 * @param points      Pointer to a universal matrix of type MUTIL\_DOUBLE
 *                    which holds the points at which the density estimate
 *                    is to be evaluated. The number of columns in this
 *                    matrix must be identical to that of {\bf data}. This
 *                    parameter may also be NULL, in which case the density
 *                    will be evaluated at the points given in {\bf data}.
 *
 * @param intrp_ptr    Parameter used for user interrupts.
 *
 * @param kde         The density estimates at the points given in
 *                    {\bf points}. This length of this matrix will be equal
 *                    to the number of rows in {\bf points}, or {\bf data} if
 *                    {\bf points} is NULL.
 * @see frauniv_embed
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_kernel_density_estimate(
  const univ_mat         *data,
  const univ_mat         *points,
  void                   *intr_ptr,
  univ_mat               *kde );

#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_KDE_H_ */
