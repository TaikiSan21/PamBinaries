
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_sig_conv.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from RS
   for the MUTILS signal library.

   Functions wrapped:

   siguniv_correlate()
*/

#include "sig_conv.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** One- or two-dimensional signal correlation, with arbitrary step and phase.
 * @source RS\_sig\_conv.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_signal_correlate", in.sig, kernel, row.step, col.step, row.overlap, col.overlap, boundary
 * @return             An RS ... containing ...
 * @param in.sig       Pointer to an RS object containing ... in.sig
 * @param kernel       Pointer to an RS object containing ... kernel
 * @param row.step     Pointer to an RS object containing ... row.step
 * @param col.step     Pointer to an RS object containing ... col.step
 * @param row.overlap  Pointer to an RS object containing ... row.overlap
 * @param col.overlap  Pointer to an RS object containing ... col.overlap
 * @param boundary     Pointer to an RS object containing ... boundary
 * @see siguniv_convolve
 * @see _mutil_boundary_type
 * @see Matrix Data Types
 * @see Interrupt Handling
*/
EXTERN_R SEXP RS_signal_correlate(
 SEXP pr_series,
 SEXP pr_kernel,
 SEXP pr_boundary)
{
  SEXP                 pr_ret_result;     
  mutil_boundary_type  boundary;            
  mutil_data_type      type;                
  mutil_errcode        err;                 
  univ_mat             series;              
  univ_mat             kernel;              
  univ_mat             result;             
  void                 *VPNULL = NULL;      
  memlist              list;
  sint32               n_series;
  sint32               n_kernel;
  sint32               n_out;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_series to series */
  READ_MATRIX_REGISTER( pr_series, &series );

  /* ... pr_kernel to kernel */
  READ_MATRIX_REGISTER( pr_kernel, &kernel );

  /* ... pr_boundary to boundary */
  err = mutil_boundary_type_from_R( pr_boundary, &boundary );
  CHECK_CONVERSION( mutil_boundary_type, pr_boundary, &boundary );

  /* initialize parameters */
  n_series = MATUNIV_NELEM( &series );
  n_kernel = MATUNIV_NELEM( &kernel );
  n_out    = n_series + n_kernel - 1;

  /* force the input vectors to be single-rows */
  series.mat.dblmat.ncol  = n_series;
  series.mat.dblmat.nrow  = 1;
  kernel.mat.dblmat.ncol  = n_kernel;
  kernel.mat.dblmat.nrow  = 1;

  /* allocate space for the result */
  err = matuniv_malloc_register( &result, 1, n_out, MUTIL_DOUBLE, &list );
  if ( err ){
    MUTIL_FREE_WARN( memlist,  &list );
    PROBLEM "Problem allocating and registering memory" ERROR;
  }

  /* Call the function */
  err = siguniv_convolve( &series, &kernel, 1, 1, 1, 1, boundary, VPNULL, &result );
  if ( err ){
    MUTIL_FREE_WARN( memlist, &list );
    PROBLEM "Problem calling siguniv_convolve() function" ERROR;
  }

  /* create the output RS object */
  err = matuniv_to_R( &result, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_result );
  MUTIL_FREE_WARN( memlist, &list );
  if ( err ) {
      PROBLEM "Unable to convert output data to RS format" ERROR;
  }

  return pr_ret_result;
}
