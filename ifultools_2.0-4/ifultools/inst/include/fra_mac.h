
/* $File: //depot/Research/ifultools/pkg/ifultools/inst/include/fra_mac.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file */

#ifndef IN_FRA_MACRO_H_
#define IN_FRA_MACRO_H_

/* This file contains fractal related macro definitions */

/* Headers to define types used in this file */
#include "ut_type.h"

#ifdef __cplusplus
extern "C" {
#endif


#define NBOX 256
#define LOG2_NBOX 8
#define II 100000000
#define INDEX(a,b) (NBOX*(((sint32)(a)+II)%NBOX) + ((sint32)(b)+II)%NBOX)
#define MAXSAMPLE 50000
#define MAXDIM 20
#define MAX(a,b) (a)>(b) ? (a) : (b)
#define MIN(a,b) (a)<(b) ? (a) : (b)
#define LOG(a) (((a)>1e-20)? log((a)): log(1e-20))


/** Phase space point separation distance.
 * This macro is used to calculate the point separation of two
 * points in the phase space. The separation is not a
 * Euclidean distance. Rather it is simply the difference
 * in a single coordinate of the embedded points, based on
 * a delay coordinate embedding. The user supplies a pointer
 * to the time series, the indices of the two points (index base zero),
 * and the time delay corresponding to the desired coordinate.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include fra\_mac.h
 * @source fra\_mac.h
 * @library fractal
 * @usage #POINT_SEPARATION( &time_series, 20, 30, 10 );#
 * @param SERIES_PTR Pointer to the time series.
 * @param TIMEA Index of the first point.
 * @param TIMEB Index of the second point.
 * @param TIME_LAG The time delay corresponding to the coordinate
 *                 over which the distance is calculated.
 * @see frauniv_embed_neighbors
 */

#undef POINT_SEPARATION
#define POINT_SEPARATION( SERIES_PTR, TIMEA, TIMEB, TIME_LAG )  \
(double) MUTIL_ABS( (double) ( SERIES_PTR[ TIMEA + TIME_LAG ] - \
	  SERIES_PTR[ TIMEB + TIME_LAG ] ) )


#ifdef __cplusplus
}
#endif

#endif /* IN_FRA_MACRO_H_ */
