
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_kde.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS fractal library.

   Functions wrapped:

   frauniv_kernel_density_estimate()
*/

#include "fra_kde.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Non-parametric kernel density estimation.
 * @source RS\_fra\_kde.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_kernel_density_estimate", data, points))#
 * @return        An R ... containing ...
 * @param data    Pointer to an R object containing ... data
 * @param points  Pointer to an R object containing ... points
 * @see frauniv_embed
*/
EXTERN_R SEXP RS_fractal_kernel_density_estimate(
 SEXP pr_data,
 SEXP pr_points )
{
  SEXP             pr_ret_result;      
  mutil_data_type  type;            
  mutil_errcode    err;             
  univ_mat         data;            
  univ_mat         result;             
  univ_mat         points;          
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_data to data */
  READ_MATRIX_REGISTER( pr_data, &data );

  /* ... pr_points to points */
  READ_MATRIX_REGISTER( pr_points, &points );

  /* Call the function */
  err = frauniv_kernel_density_estimate(
    &data,
    &points,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_kernel_density_estimate, &result, &pr_ret_result );
}

