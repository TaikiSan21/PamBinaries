
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_wav_mrd.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS wavelets library.

   Functions wrapped:

   wavuniv_transform_packet_detail
 */

#include "wav_modw.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Calculate the detail sequence for a specified crystal of a 1-D wavelet transform.
 * @source RS\_wav\_mrd.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_wavelets_transform_packet_detail", transform, filters, level, node, xformtype))#
 * @return           An R ... containing ...
 * @param transform  Pointer to an R object containing ... transform
 * @param filters    Pointer to an R object containing ... filters
 * @param level      Pointer to an R object containing ... level
 * @param node       Pointer to an R object containing ... node
 * @param xformtype  Pointer to an R object containing ... xformtype
 * @see wavuniv_transform_discrete_wavelet_convolution
 * @see wavuniv_transform_packet
 * @see wavuniv_transform_maximum_overlap
 * @see wavuniv_transform_maximum_overlap_packet
*/
EXTERN_R SEXP RS_wavelets_transform_packet_detail(
 SEXP pr_transform,
 SEXP pr_filters,
 SEXP pr_level,
 SEXP pr_node,
 SEXP pr_xformtype )
{
  SEXP             pr_ret_result;   
  mat_set          filters;         
  mat_set          transform;       
  mutil_errcode    err;             
  sint32           level;           
  sint32           node;            
  univ_mat         result;          
  void             *VPNULL = NULL;  
  wav_transform    xformtype;       
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_transform to transform */
  READ_MATSET_REGISTER( pr_transform, MUTIL_DOUBLE, &transform );

  /* ... pr_filters to filters */
  READ_MATSET_REGISTER( pr_filters, MUTIL_DOUBLE, &filters );

  /* ... pr_level to level */
  err = sint32_from_R( pr_level, &level );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert sint32 type argument pr_level to level" );

  /* ... pr_node to node */
  err = sint32_from_R( pr_node, &node );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert sint32 type argument pr_node to node" );

  /* ... pr_xformtype to xformtype */
  WAV_TRANSFORM_FROM_R( pr_xformtype, &xformtype );

  /* Call the function */
  err = wavuniv_transform_packet_detail(
    &transform,
    &filters,
    level,
    node,
    xformtype,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( wavuniv_transform_packet_detail, &result, &pr_ret_result );
}
