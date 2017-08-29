
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_sig_win.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS signal library.

   Functions wrapped:

   siguniv_taper()
*/

#include "sig_win.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"


/** Oracle function for obtaining a particular taper/window.
 * @source RS\_sig\_win.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_signal_taper", taper, nrow, ncol, param, normalize))#
 * @return           An R ... containing ...
 * @param taper      Pointer to an R object containing ... taper
 * @param nrow       Pointer to an R object containing ... nrow
 * @param ncol       Pointer to an R object containing ... ncol
 * @param param      Pointer to an R object containing ... param
 * @param normalize  Pointer to an R object containing ... normalize
 * @see siguniv_window_rectangle
 * @see siguniv_window_triangle
 * @see siguniv_window_raised_cosine
 * @see siguniv_window_hanning
 * @see siguniv_window_hamming
 * @see siguniv_window_blackman
 * @see siguniv_window_nuttall
 * @see siguniv_window_gaussian
 * @see siguniv_window_kaiser
 * @see siguniv_window_chebyshev
 * @see siguniv_window_born_jordan
 * @see siguniv_window_sinusoidal
 * @see siguniv_window_parzen
 * @see siguniv_window_papoulis
 * @see siguniv_window_daniell
*/
EXTERN_R SEXP RS_signal_taper(
 SEXP pr_taper,
 SEXP pr_nrow,
 SEXP pr_ncol,
 SEXP pr_param,
 SEXP pr_normalize )
{
  SEXP             pr_ret_result;   
  boolean          normalize;       
  double           param;           
  mutil_errcode    err;             
  sig_taper_type   taper;           
  sint32           ncol;            
  sint32           nrow;            
  univ_mat         result;          
  void             *VPNULL = NULL;  

  /* Avoid lint warning */

  (void) whatssi;

  /* Conversion of input data ... */

  /* ... pr_taper to taper */
  err = sig_taper_type_from_R( pr_taper, &taper );
  if ( err ){
    PROBLEM "Unable to convert sig_taper_type type argument pr_taper to taper" ERROR;
  }

  /* ... pr_nrow to nrow */
  err = sint32_from_R( pr_nrow, &nrow );
  if ( err ){
    PROBLEM "Unable to convert sint32 type argument pr_nrow to nrow" ERROR;
  }

  /* ... pr_ncol to ncol */
  err = sint32_from_R( pr_ncol, &ncol );
  if ( err ){
    PROBLEM "Unable to convert sint32 type argument pr_ncol to ncol" ERROR;
  }

  /* ... pr_param to param */
  err = double_from_R( pr_param, &param );
  if ( err ){
    PROBLEM "Unable to convert double type argument pr_param to param" ERROR;
  }

  /* ... pr_normalize to normalize */
  err = boolean_from_R( pr_normalize, &normalize );
  if ( err ){
    PROBLEM "Unable to convert boolean type argument pr_normalize to normalize" ERROR;
  }

  /* Call the function */
  err = siguniv_taper(
    taper,
    nrow,
    ncol,
    param,
    normalize,
    VPNULL,
    &result );
  if ( err ){
    PROBLEM "Problem calling siguniv_taper() function" ERROR;
  }

  /* create the output R object */
  err = matuniv_to_R( &result, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_result );
  MUTIL_FREE_WARN( matuniv, &result );
  if ( err ) {
      PROBLEM "Unable to convert output data to R format" ERROR;
  }

  return pr_ret_result;
}

