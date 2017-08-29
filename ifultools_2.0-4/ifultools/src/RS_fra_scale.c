
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_scale.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS fractal library.

   Functions wrapped:

   frauniv_piecwise_linear_segmentation()
*/

#include "fra_scale.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Piecewise linear segmentation of a time series.
 * @source RS\_fra\_scale.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_piecwise_linear_segmentation", xdata, ydata, n.fit, angle.tolerance))#
 * @return                 An R ... containing ...
 * @param xdata            Pointer to an R object containing ... xdata
 * @param ydata            Pointer to an R object containing ... ydata
 * @param n.fit            Pointer to an R object containing ... n.fit
 * @param angle.tolerance  Pointer to an R object containing ... angle.tolerance
 * @see frauniv_embed
*/
EXTERN_R SEXP RS_fractal_piecwise_linear_segmentation(
 SEXP pr_xdata,
 SEXP pr_ydata,
 SEXP pr_n_fit,
 SEXP pr_angle_tolerance )
{
  SEXP              pr_ret_result;    
  double            angle_tolerance;  
  mutil_data_type   type;             
  mutil_errcode     err;              
  sint32            n_fit;            
  univ_mat          result;           
  univ_mat          xdata;            
  univ_mat          ydata;            
  void              *VPNULL = NULL;   
  memlist           list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_xdata to xdata */
  READ_MATRIX_REGISTER( pr_xdata, &xdata );

  /* ... pr_ydata to ydata */
  READ_MATRIX_REGISTER( pr_ydata, &ydata );

  /* ... pr_n_fit to n_fit */
  SINT32_FROM_R( pr_n_fit, &n_fit );

  /* ... pr_angle_tolerance to angle_tolerance */
  err = double_from_R( pr_angle_tolerance, &angle_tolerance );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert double type argument pr_angle_tolerance to angle_tolerance" );

  /* Call the function */
  err = frauniv_piecwise_linear_segmentation(
    &xdata,
    &ydata,
    n_fit,
    angle_tolerance,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_piecwise_linear_segmentation, &result, &pr_ret_result );
}

