
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mth_var.c $: $Revision: #1 $, $Date: 2008/03/21 $";

#include "mth_var.h"
#include "mat_stat.h"
#include "mat_type.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "sig_tran.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_math.h"
#include "ut_mem.h"
#include <math.h>

#undef LOCALDEF_CHECK_NULL_MATRIX_POINTER
#define LOCALDEF_CHECK_NULL_MATRIX_POINTER( DATA_PTR, DATA_TYPE,  \
                                     TYPE_PREFIX )                \
   err = TYPE_PREFIX ## _validate( DATA_PTR );                    \
   if ( err ) return err;                                         \
   if ( DATA_PTR == ( DATA_TYPE * ) NULL ){                       \
     MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" );    \
     return MUTIL_ERR_NULL_POINTER;                               \
   }

/* One-sided autocovariance sequence using DFTs */
/* Documented in mth_var.h                      */
/* Written by William Constantine               */

mutil_errcode mthuniv_acvs(
  const univ_mat   *time_series,
  const boolean     biased,
  const boolean     recenter,
  void             *intrp_ptr,
  univ_mat         *result )
{
  dcomplex      *pdft;            /* pointer to the DFT data                                          */
  dcomplex      *ppadded;         /* pointer to the data in the zero padded series                    */
  double         mean = 0.0;      /* the mean of the input series                                     */
  double         variance;        /* the variance of the input series                                 */
  double        *presult;         /* pointer to the data in the result matrix                         */
  double        *pseries;         /* pointer to the data in the series matrix                         */
  mutil_errcode  err;             /* MUTIL error code                                                 */
  sint32         n_sample;        /* number of points in series                                       */
  sint32         i;               /* counting index                                                   */
  sint32         zeropad_length;  /* total number of points to feed to DFT routine                    */
  univ_mat       dft;             /* complex vector containing the DFT data                           */
  univ_mat       idft;            /* complex vector containing the inverse DFT data                   */
  univ_mat       padded_series;   /* vector containing the zero padded series sent to the DFT routine */
  memlist        list;            /* dynamically allocated memory manager list                        */

  MUTIL_TRACE( "Start mthuniv_acvs()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* begin I/O checking */
  /* validate matrices  */

  LOCALDEF_CHECK_NULL_MATRIX_POINTER( time_series, univ_mat, matuniv  );

  /* ... if the type is double */

  if ( time_series->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Time series matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... if it is a vector */

  if ( !MATANY_IS_VEC( &( time_series->mat.dblmat ) ) ){
    MUTIL_ERROR( "Time series matrix must be a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* obtain and set sizes */

  n_sample       = MATUNIV_NELEM( time_series );
  zeropad_length = 2 * n_sample;

  /* recenter the data if requested */

  if (recenter){
    err = matuniv_mean_variance( time_series, FALSE, intrp_ptr, &mean, &variance );
    if (err) return(err);
  }
  else mean = 0.0;

  /* malloc space for matrices */

  err = matuniv_malloc_register( &padded_series, zeropad_length, 1, MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &dft, zeropad_length, 1, MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &idft, zeropad_length, 1, MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( result, MATUNIV_NROW( time_series ),
    MATUNIV_NCOL( time_series ), MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* set pointers */

  pseries = (double *) MATUNIV_DATA( time_series );
  presult = (double *) MATUNIV_DATA( result );
  ppadded = (dcomplex *) MATUNIV_DATA( &padded_series );
  pdft    = (dcomplex *) MATUNIV_DATA( &dft );

  /* explicitly fill the zero padded dft vector.
     this operation zeros out the imaginary portion, recenters the data,
     zero pads the real
     portion beyond the series length, and ensures a column vector
     migration of the input series (which may be a row vector) */

  for (i = 0; i < zeropad_length; i++){

    if (i < n_sample) ppadded[ i ].re = pseries[ i ] - mean;
    else ppadded[ i ].re = 0;

    ppadded[ i ].im = 0;
  }

  /* perform DFT */

  err = siguniv_transform_discrete_fourier( &padded_series, FALSE, intrp_ptr, &dft );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* replace elements of dft by the product of its conjugate and itself */

  for (i = 0; i < zeropad_length; i++){
    pdft[ i ].re = (pdft[ i ].re * pdft[ i ].re  + pdft[ i ].im * pdft[ i ].im);
    pdft[ i ].im = 0;
  }

  /* perform the inverse dft */

  err = siguniv_transform_discrete_fourier( &dft, TRUE, intrp_ptr, &idft );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* reassign pointers */

  pdft = idft.mat.cpxmat.data;

  /* fill up result array and normalize...

     A NOTE ON NORMALIZATION:

     Without any DFT normalization in either direction, the normalization factor
     would be 2 * N. However, the MUTILS inverse DFT does normalize by the number of Fourier
     coefficients which (in our case with the padded series) is 2 * N. So, by coincidence
     only, we do not have to normalize the inverse DFT output. however, the standard normalization
     for the acvs is used and its value is N for the biased case (for all lags) and
     N - |t| for the unbiased case where t = 0, ..., N - 1                             */

  for ( i = 0; i < n_sample; i++ ){

    if ( biased ) presult[ i ] = pdft[ i ].re / (double) n_sample;
    else presult[ i ] = pdft[ i ].re / (double) ( n_sample - i );
  }

  /* free nodes corresponding to registered
     modwpt matrix memory, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free malloced space and corresponding linked list structures */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with mthuniv_acvs()" );

  return MUTIL_ERR_OK;
}
