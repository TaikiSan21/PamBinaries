
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_tdmi.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_FRA_TDMI_H_
#define IN_FRA_TDMI_H_


#include "mat_type.h"


/* This file contains function declarations for computing mutual */
/* information between a time series and its delayed versions    */

#ifdef __cplusplus
extern "C" {
#endif

/** Time delayed mutual information.
 *
 * This function finds the mutual information between a time series and its
 * delayed versions. The time delayed mutual information (TDMI) generalizes
 * the correlation function. For a time series, $x_n,$ the TDMI, denoted by
 * $I(k),$ is a measure of the information we have about $x_{n+k}$ if we
 * know $x_n.$
 *
 * Input {\bf time\_series} is a row or column vector, of any real type,
 * containing a scalar time series. Input {\bf lags} is a vector of
 * non-negative integers containing lags $(k)$ at which the TDMI, $I(k),$
 * will be computed.
 *
 * This function calls the kernel density estimation routine
 * \Ref{frauniv_kernel_density_estimate} for estimating the joint
 * probability density function of a two dimensional and single
 * dimensional delayed embedding of a single variable time series.
 * These density estimates are used
 * in the formulation of the mutual information.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_tdmi.h
 * @source fra\_tdmi.c
 * @library fractal
 * @usage #err = frauniv_time_delayed_mutual_information(&ts,&tau,int_ptr,&mi);#
 * @return Standard mutils error/OK code.
 * @param time_series Pointer to a universal matrix of any real type,
 *                    which contains a row or column vector representing a
 *                    time series.
 * @param lags      Pointer to a universal matrix of type MUTIL\_SINT32,
 *                    which specifies the integer lags for which
 *                    {\bf time\_series} will be delayed. The mutual
 *                    information will be computed between
 *                    {\bf time\_series} and each of its delayed
 *                    versions. This must be a row or column vector.
 * @param intrp_ptr    Parameter used for user interrupts.
 * @param tdmi        Pointer to an unallocated universal matrix. Upon
 *                    return, a matrix of type MUTIL\_DOUBLE which holds
 *                    the computed mutual information values. The dimensions
 *                    are identical to those of {\bf lags}.
 * @see frauniv_kernel_density_estimate
 */
 MUTIL_LIBEXPORT mutil_errcode frauniv_time_delayed_mutual_information(
   const univ_mat   *time_series,
   const univ_mat   *lags,
   void             *intrp_ptr,
   univ_mat         *tdmi );


#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_TDMI_H_ */
