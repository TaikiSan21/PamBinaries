
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/wav_math.c $: $Revision: #1 $, $Date: 2008/03/21 $";

#include "wav_math.h"
#include "mat_stat.h"
#include "mat_type.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrp.h"
#include "ut_math.h"
#include <math.h>

/* Table lookup using linear interpolation for 1D matrices */
/* Function documented in wav_math.h                       */
/* Written by William Constantine                          */
mutil_errcode wavuniv_statistic_interpolation_linear(
  const univ_mat *y,
  const univ_mat *x,
  const univ_mat *xi,
  void            *intrp_ptr,
  univ_mat        *yi )
{

  double        dx;               /* sampling interval in x                  */
  double        num_ops = 0.0;    /* variable used for interrupt handling    */
  double       *px;               /* pointer to x data                       */
  double       *pxi;              /* pointer to xi data                      */
  double       *py;               /* pointer to y data                       */
  double       *pyi;              /* pointer to yi                           */
  mutil_errcode err;              /* MUTILS error code                       */
  sint32        N;                /* number of coefficients in x or y        */
  sint32        Ninterp;          /* number of coefficients in xi or yi      */
  sint32        i;                /* counting variable                       */
  sint32        ilow;             /* lowest index in fit interval            */

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_statistic_interpolation_linear()" );

  /* avoid lint warning */

  (void) whatssi;

  /* begin I/O checking */
  /* validate matrices  */

  err = matuniv_validate( x );
  if ( err ) return err;

  err = matuniv_validate( y );
  if ( err ) return err;

  err = matuniv_validate( xi );
  if ( err ) return err;

  err = matuniv_validate( yi );
  if ( err ) return err;

  if ( (x->type  != MUTIL_DOUBLE) ||
       (y->type  != MUTIL_DOUBLE) ||
       (xi->type != MUTIL_DOUBLE) ||
       (yi->type != MUTIL_DOUBLE) ){
    MUTIL_ERROR( "All matrices must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( (MATUNIV_NROW( x ) != 1) &  (MATUNIV_NCOL( x ) != 1) ){
    MUTIL_ERROR( "Input matrix x must be a single row or column vector." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( (MATUNIV_NROW( xi ) != 1) &  (MATUNIV_NCOL( xi ) != 1) ){
    MUTIL_ERROR( "Input matrix xi must be a single row or column vector." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( MATUNIV_NELEM( y ) != MATUNIV_NELEM( x ) ){
    MUTIL_ERROR( "Input matrix y must be the same length as input matrix x." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( MATUNIV_NELEM( yi ) != MATUNIV_NELEM( xi ) ){
    MUTIL_ERROR( "Output matrix yi must be the same length as input matrix xi." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( MATUNIV_NELEM( x ) < 2 ){
    MUTIL_ERROR( "Number of elements in input matrix x must be greater than 1." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* obtain size information */

  N       = MATUNIV_NELEM( x  );
  Ninterp = MATUNIV_NELEM( xi );

  /* set pointers */

  px  = &(x->mat.dblmat.data[0]);
  py  = &(y->mat.dblmat.data[0]);
  pxi = &(xi->mat.dblmat.data[0]);
  pyi = &(yi->mat.dblmat.data[0]);

  /* initialize variables */

  dx = px[1] - px[0];

  for (i = 0; i < Ninterp; i++){

    if ( pxi[i] <= px[0] ) pyi[i] = py[0];
    else if ( pxi[i] >= px[N - 1] ) pyi[i] = py[N - 1];
    else{

      ilow   = (sint32) floor( (pxi[i] - px[0]) / dx );
      pyi[i] = ( pxi[i] - px[ilow] ) * ( py[ilow + 1] - py[ilow] ) / dx +
        py[ilow];
    }

    /* Check for interrupts */
    num_ops += 10.0;

    if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      return MUTIL_ERR_INTERRUPT;
    }
  }

  MUTIL_TRACE( "Done with wavuniv_statistic_interpolation_linear()" );

  return MUTIL_ERR_OK;
}



