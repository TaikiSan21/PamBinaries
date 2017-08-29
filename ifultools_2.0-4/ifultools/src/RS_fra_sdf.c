
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_sdf.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS fractal library.

   Functions wrapped:

   frauniv_spectral_density_function_direct()
   frauniv_spectral_density_function_lag_window()
   frauniv_spectral_density_function_wosa()
   frauniv_spectral_density_function_multitaper()
*/

#include "fra_sdf.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Nonparametric cross-spectral density function estimation.
 * @source RS\_fra\_sdf.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_spectral_density_function_direct", time.series, taper, center, recenter, single.sided, npad))#
 * @return              An R ... containing ...
 * @param time.series   Pointer to an R object containing ... time.series
 * @param taper         Pointer to an R object containing ... taper
 * @param center        Pointer to an R object containing ... center
 * @param recenter      Pointer to an R object containing ... recenter
 * @param single.sided  Pointer to an R object containing ... single.sided
 * @param npad          Pointer to an R object containing ... npad
 * @see sig_taper
 * @see _fra_sdf_taper
 * @see fra_spectral_density_function_lag_window
 * @see fra_spectral_density_function_wosa
 * @see fra_spectral_density_function_multitaper
*/
EXTERN_R SEXP RS_fractal_spectral_density_function_direct(
 SEXP pr_time_series,
 SEXP pr_taper,
 SEXP pr_center,
 SEXP pr_recenter,
 SEXP pr_single_sided,
 SEXP pr_npad )
{
  SEXP             pr_ret_result;   
  boolean          center;          
  boolean          recenter;        
  boolean          single_sided;    
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           npad;            
  univ_mat         result;          
  univ_mat         taper;           
  univ_mat         time_series;     
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_taper to taper */
  READ_MATRIX_REGISTER( pr_taper, &taper );

  /* ... pr_center to center */
  BOOLEAN_FROM_R( pr_center, &center );

  /* ... pr_recenter to recenter */
  BOOLEAN_FROM_R( pr_recenter, &recenter );

  /* ... pr_single_sided to single_sided */
  BOOLEAN_FROM_R( pr_single_sided, &single_sided );

  /* ... pr_npad to npad */
  SINT32_FROM_R( pr_npad, &npad );

  /* Call the function */
  err = frauniv_spectral_density_function_direct(
    &time_series,
    &taper,
    center,
    recenter,
    single_sided,
    npad,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_spectral_density_function_direct, &result, &pr_ret_result );
}

/** Nonparametric cross-spectral density function estimation.
 * @source RS\_fra\_sdf.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_spectral_density_function_lag_window", time.series, lag.window, taper, center, recenter, single.sided, npad))#
 * @return              An R ... containing ...
 * @param time.series   Pointer to an R object containing ... time.series
 * @param lag.window    Pointer to an R object containing ... lag.window
 * @param taper         Pointer to an R object containing ... taper
 * @param center        Pointer to an R object containing ... center
 * @param recenter      Pointer to an R object containing ... recenter
 * @param single.sided  Pointer to an R object containing ... single.sided
 * @param npad          Pointer to an R object containing ... npad
 * @see sig_taper
 * @see _fra_sdf_taper
 * @see fra_spectral_density_function_direct
 * @see fra_spectral_density_function_wosa
 * @see fra_spectral_density_function_multitaper
*/
EXTERN_R SEXP RS_fractal_spectral_density_function_lag_window(
 SEXP pr_time_series,
 SEXP pr_lag_window,
 SEXP pr_taper,
 SEXP pr_center,
 SEXP pr_recenter,
 SEXP pr_single_sided,
 SEXP pr_npad )
{
  SEXP             pr_ret_result;   
  boolean          center;          
  boolean          recenter;        
  boolean          single_sided;    
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           npad;            
  univ_mat         lag_window;      
  univ_mat         result;          
  univ_mat         taper;           
  univ_mat         time_series;     
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_lag_window to lag_window */
  READ_MATRIX_REGISTER( pr_lag_window, &lag_window );

  /* ... pr_taper to taper */
  READ_MATRIX_REGISTER( pr_taper, &taper );

  /* ... pr_center to center */
  BOOLEAN_FROM_R( pr_center, &center );

  /* ... pr_recenter to recenter */
  BOOLEAN_FROM_R( pr_recenter, &recenter );

  /* ... pr_single_sided to single_sided */
  BOOLEAN_FROM_R( pr_single_sided, &single_sided );

  /* ... pr_npad to npad */
  SINT32_FROM_R( pr_npad, &npad );

  /* Call the function */
  err = frauniv_spectral_density_function_lag_window(
    &time_series,
    &lag_window,
    &taper,
    center,
    recenter,
    single_sided,
    npad,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_spectral_density_function_lag_window, &result, &pr_ret_result );
}

/** Nonparametric cross-spectral density function estimation.
 * @source RS\_fra\_sdf.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_spectral_density_function_wosa", time.series, taper, overlap, center, recenter, single.sided, npad))#
 * @return              An R ... containing ...
 * @param time.series   Pointer to an R object containing ... time.series
 * @param taper         Pointer to an R object containing ... taper
 * @param overlap       Pointer to an R object containing ... overlap
 * @param center        Pointer to an R object containing ... center
 * @param recenter      Pointer to an R object containing ... recenter
 * @param single.sided  Pointer to an R object containing ... single.sided
 * @param npad          Pointer to an R object containing ... npad
 * @see sig_taper
 * @see _fra_sdf_taper
 * @see fra_spectral_density_function_direct
 * @see fra_spectral_density_function_lag_window
 * @see fra_spectral_density_function_multitaper
*/
EXTERN_R SEXP RS_fractal_spectral_density_function_wosa(
 SEXP pr_time_series,
 SEXP pr_taper,
 SEXP pr_overlap,
 SEXP pr_center,
 SEXP pr_recenter,
 SEXP pr_single_sided,
 SEXP pr_npad )
{
  SEXP             pr_ret_result;   
  boolean          center;          
  boolean          recenter;        
  boolean          single_sided;    
  double           overlap;         
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           npad;            
  univ_mat         result;          
  univ_mat         taper;           
  univ_mat         time_series;     
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_taper to taper */
  READ_MATRIX_REGISTER( pr_taper, &taper );

  /* ... pr_overlap to overlap */
  DOUBLE_FROM_R( pr_overlap, &overlap );

  /* ... pr_center to center */
  BOOLEAN_FROM_R( pr_center, &center );

  /* ... pr_recenter to recenter */
  BOOLEAN_FROM_R( pr_recenter, &recenter );

  /* ... pr_single_sided to single_sided */
  BOOLEAN_FROM_R( pr_single_sided, &single_sided );

  /* ... pr_npad to npad */
  SINT32_FROM_R( pr_npad, &npad );

  /* Call the function */
  err = frauniv_spectral_density_function_wosa(
    &time_series,
    &taper,
    overlap,
    center,
    recenter,
    single_sided,
    npad,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_spectral_density_function_wosa, &result, &pr_ret_result );
}

/** Nonparametric cross-spectral density function estimation.
 * @source RS\_fra\_sdf.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_spectral_density_function_multitaper", time.series, taper, center, recenter, single.sided, npad, npad))#
 * @return              An R ... containing ...
 * @param time.series   Pointer to an R object containing ... time.series
 * @param taper         Pointer to an R object containing ... taper
 * @param center        Pointer to an R object containing ... center
 * @param recenter      Pointer to an R object containing ... recenter
 * @param single.sided  Pointer to an R object containing ... single.sided
 * @param npad          Pointer to an R object containing ... npad
 * @see sig_taper
 * @see _fra_sdf_taper
 * @see fra_spectral_density_function_direct
 * @see fra_spectral_density_function_lag_window
 * @see fra_spectral_density_function_multitaper
*/
EXTERN_R SEXP RS_fractal_spectral_density_function_multitaper(
 SEXP pr_time_series,
 SEXP pr_taper,
 SEXP pr_center,
 SEXP pr_recenter,
 SEXP pr_single_sided,
 SEXP pr_npad )
{
  SEXP             pr_ret_result;   
  boolean          center;          
  boolean          recenter;        
  boolean          single_sided;    
  mutil_data_type  type;            
  mutil_errcode    err;             
  sint32           npad;            
  univ_mat         result;          
  univ_mat         taper;           
  univ_mat         time_series;     
  void             *VPNULL = NULL;  
  memlist          list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_time_series, &time_series );

  /* ... pr_taper to taper */
  READ_MATRIX_REGISTER( pr_taper, &taper );

  /* ... pr_center to center */
  BOOLEAN_FROM_R( pr_center, &center );

  /* ... pr_recenter to recenter */
  BOOLEAN_FROM_R( pr_recenter, &recenter );

  /* ... pr_single_sided to single_sided */
  BOOLEAN_FROM_R( pr_single_sided, &single_sided );

  /* ... pr_npad to npad */
  SINT32_FROM_R( pr_npad, &npad );

  /* Call the function */
  err = frauniv_spectral_density_function_multitaper(
    &time_series,
    &taper,
    center,
    recenter,
    single_sided,
    npad,
    VPNULL,
    &result );
  CONVERT_MATRIX_AND_RETURN( frauniv_spectral_density_function_multitaper, &result, &pr_ret_result );
}

