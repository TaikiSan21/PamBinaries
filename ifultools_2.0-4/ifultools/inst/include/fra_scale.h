
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_scale.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_SCALE_H_
#define IN_FRA_SCALE_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

/* This file contains function declarations the mutils fractal library. */

#ifdef __cplusplus
extern "C" {
#endif

/** Piecewise linear segmentation of a time series.
 * Locates the change-points of time series based on a piecewise linear
 * segmentation algorithm. Given a window size (n\_fit) and an angle tolerance
 * (angle\_tolerance), the segmentation algorithm starts by finding the
 * slope of the first n\_fit points of the series via least squares
 * regression. The window is slid over one point to the right, the points within
 * the new window are regressed, and the new slope is compared to the old slope.
 * If the change in slope exceeds the specified angle\_tolerance, a change-point
 * is recorded as the rightmost point of the previous iteration's window. The
 * routine then picks up again starting at the point just to the right of the
 * change-point. If the change in slope does not exceed the specified
 * angle\_tolerance, then the old slope is updated (in a running average sense),
 * and the algorithm continues as usual. The window is slid along the
 * until its rightmost point reaches the edge of the time series.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_scale.h
 * @source fra\_scale.c
 * @library fractal
 * @usage #err = frauniv_piecwise_linear_segmentation( &xdata, &ydata, n_fit, angle_tolerance, intrp_ptr, &result );#
 * @return Standard mutils error/OK code.
 * @param xdata Pointer to a pre-allocated single-row or
 *              single-column universal matrix of type MUTIL\_DOUBLE
 *              containing the independent variable of the time series.
 * @param ydata Pointer to a pre-allocated single-row or
 *              single-column universal matrix of type MUTIL\_DOUBLE
 *              containing the dependent variables of the time series.
 * @param n_fit An integer denoting the window size, not to exceed the
 *              number of samples in the time series.
 * @param angle_tolerance The maximum angle in degrees that the running
 *              average of the slopes in the current set of points
 *              can change relative to the slope of the data calculated
 *              in the most current (rightmost) window before a
 *              change-point is recorded.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param result Pointer to a universal matrix of type MUTIL\_SINT32
 *               which (upon return) will contain a single-column
 *               matrix containing the change-points.
 *               If no change points are found via the segmentation,
 *               then N - 1 is returned as the only change-point
 *               where N is the number of samples in the time
 *               series. The memory for this matrix is allocated
 *               within the function.
 * @see frauniv_embed
 */
MUTIL_LIBEXPORT mutil_errcode frauniv_piecwise_linear_segmentation(
 const univ_mat *xdata,
 const univ_mat *ydata,
 const sint32    n_fit,
 const double    angle_tolerance,
 void           *intrp_ptr,
 univ_mat       *result );

#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_SCALE_H_ */
