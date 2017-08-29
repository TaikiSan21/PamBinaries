/* This file contains wrapper functions, callable from R, for the LAPACK library. */

#ifndef IN_SAPA_LAPACK_H_
#define IN_SAPA_LAPACK_H_

#include "R.h"
#ifdef USING_R
#include "Rinternals.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Discrete prolate spheroidal sequences for multitaper SDF estimation.
 * @source sapa\_lapack.c
 * @author William Constantine
 * @usage #.Call( "R_sapa_dpss", nsample, ntaper, bandwidth)#
 * @return          An R ... containing ...
 * @param nsample   Number of points in each taper (integer)
 * @param ntaper    Number of tapers (integer)
 * @param bandwidth Resolution bandwidth (double)
*/
SEXP R_sapa_dpss(SEXP nsample, SEXP ntaper, SEXP bandwidth);

#ifdef __cplusplus
}
#endif
#endif
#endif /* IN_SAPA_LAPACK_H_ */
