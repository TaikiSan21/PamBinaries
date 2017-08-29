
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/sig_win.c $: $Revision: #1 $, $Date: 2008/03/21 $";

/* This is a self-documenting doc++ file */

/*
  This file contains implementations of functions for creating
  window tapers frequently used for windowing signals
  in signal processing.
*/

#include "sig_win.h"
#include "sig_type.h"

#include "mat_arit.h"
#include "mat_any.h"
#include "mat_assn.h"
#include "mat_summ.h"
#include "mat_univ.h"

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrp.h"
#include "ut_limit.h"
#include "ut_math.h"
#include "ut_intrn.h"

/* inverse hyperbolic cosine */
#define COSHIN( arg ) log( (arg) + sqrt( (arg) * (arg) - 1.0 ) )

/* macro to assert that type is legal */
#define ASSERT_LEGAL_TYPE(TYPE)  MUTIL_ASSERT((TYPE) >= MUTIL_UINT8 && \
                                             (TYPE) <= MUTIL_DCOMPLEX )

/* allow for in-place complex multiplication */

#define LOCALDEF_COMPLEX_MULTIPLY(a,a1,a2) {\
	dcomplex Ctemp;\
	Ctemp.re = (a1)->re * (a2)->re - (a1)->im * (a2)->im;\
	Ctemp.im = (a1)->im * (a2)->re + (a1)->re * (a2)->im;\
	*(a) = Ctemp;}

/* Static functions declared here and defined at end of file */

static mutil_errcode localfn_acos_complex( dcomplex *z );
static mutil_errcode localfn_cosh_complex( dcomplex *z );
static mutil_errcode localfn_acosh_complex( dcomplex *z );
static mutil_errcode localfn_taper_normalize( double_mat *taper );
static mutil_errcode localfn_lag_window_normalize( double_mat *taper );
static mutil_errcode localfn_acosh( double x, double *result );
static double localfn_hypot( double x, double y );

/*
****************
STATIC FUNCTIONS
****************
*/

/**  Zeroth-order modified Bessel function of the first kind.
 * Calculates the zeroth-order modified Bessel function
 * of the first kind by polynomial approximation. Takes a pointer
 * to an array of doubles and the array length, and returns the
 * result in-place.
 *
 * @usage #localfn_sigdbl_besselI0(array_ptr, array_len, Intrp_ptr);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param x Pointer to array of doubles.
 * @param len Length of array.
 * @param intrp_ptr Pointer to implement user interrrupt.
 * @private
*/
static mutil_errcode localfn_sigdbl_besselI0( double *x,
    sint32 len, void *intrp_ptr )
{
    sint32 i;
    double abs_x;
    double y;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start localfn_sigdbl_besselI0()" );

    for (i = 0 ; i < len ; i++) {
        abs_x = fabs( x[i] );
        if (abs_x < 3.75) {
            y = abs_x * abs_x / 14.0625 ;
            x[i] = 1.0 + y * ( 3.5156229 + y *
                               ( 3.0899424 + y *
                                 ( 1.2067492 + y *
                                   ( 0.2659732 + y *
                                     ( 0.0360768 + y * 0.0045813 ) ) ) ) );
        }
        else {
            y = 3.75 / abs_x ;
            x[i] = exp(abs_x) *
                ( 0.39894228 + y *
                  ( 0.01328592 + y *
                    ( 0.00225319 + y *
                      ( -0.00157565 + y *
                       ( 0.00916281 + y *
                         ( -0.02057706 + y *
                           ( 0.02635537 + y *
                             ( -0.01647633 + y * 0.00392377 )
                               ) ) ) ) ) ) ) / sqrt(abs_x);
        }
        num_ops += 500;
    }

    if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
        MUTIL_ERROR( "User Interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    MUTIL_TRACE( "localfn_sigdbl_besselI0() done" );

    return MUTIL_ERR_OK;

}

/*
***********************
Double Matrix Functions
***********************
*/


/* Function documented in sig_win.h */
/* Written by James Pitton. */
/* Adapted for mutils library by Luca Cazzanti */
mutil_errcode sigdbl_window_rectangle( void *intrp_ptr, double_mat *win )
{

    mutil_errcode errcode;
    sint32        i;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start sigdbl_window_rectangle()" );

    /* avoid lint warning */
    (void) whatssi;

    /* validate matrix */

    errcode = matdbl_validate(win);
    if(errcode) return(errcode);

    /* matrix must be a vector */

    if( !MATANY_IS_VEC(win) ) {
        MUTIL_ERROR( "Input matrix must be a 1D vector" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* build window */
    /* rectangular window is all ones for length of window */

    for( i = 0; i < win->nelem; i++ ) {
        win->data[i] = 1.0;
    }

    num_ops += 5 * win->nelem;
    if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
        MUTIL_ERROR( "User Interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    MUTIL_TRACE( "sigdbl_window_rectangle() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by James Pitton */
/* Adapted for mutils library by Luca Cazzanti */
mutil_errcode sigdbl_window_triangle( void *intrp_ptr, double_mat *win )
{
    mutil_errcode errcode;

    double tmp;
    sint32 i;
    sint32 mid_pnt;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start sigdbl_window_triangle()" );

    /* check matrix */

    errcode = matdbl_validate(win);
    if(errcode) return(errcode);

    /* matrix must be a vector */

    if( !MATANY_IS_VEC(win) ) {
        MUTIL_ERROR( "Input matrix must be a vector" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* build window */
    /* triangle window is symmetric, so calculate only half */
    /* the window and assign rest by symmetry */

    tmp = 2.0 / ( win->nelem + 1 );
    mid_pnt = win->nelem / 2;

    for( i = 0; i < mid_pnt ; i++ ) {
        win->data[i] = tmp * ( i + 1 );
        win->data[ win->nelem - i - 1 ] = win->data[i];
    }

    if( win->nelem % 2 ) {
        win->data[mid_pnt] = 1;
    }

    num_ops += 5 * win->nelem;
    if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
        MUTIL_ERROR( "User Interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    MUTIL_TRACE( "sigdbl_window_triangle() done" );

    return MUTIL_ERR_OK;

}

/* Function documented in sig_win.h */
/* Written by James Pitton */
/* Adapted for mutils library by Luca Cazzanti */
mutil_errcode sigdbl_window_raised_cosine( double percent, void *intrp_ptr,
    double_mat *win )
{
    mutil_errcode errcode;

    double scale;
    double low;
    double angle;
    sint32 i;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start sigdbl_window_raised_cosine()" );

    /* check arguments*/

    if( percent < 0 || percent > 1 ) {
        MUTIL_ERROR( "Fraction of raised cosine part must be between [0,1]" );
        return MUTIL_ERR_ILLEGAL_VALUE;
    }

    errcode = matdbl_validate(win);
    if(errcode) return(errcode);

    /* matrix must be a vector */

    if( !MATANY_IS_VEC(win) ) {
        MUTIL_ERROR( "Input matrix must be a vector" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* build window */

    /* central part of window is 1 */
    /* set all values to one for now */

    for( i = 0; i < win->nelem; i++ ) {
        win->data[i] = 1.0;
    }

    num_ops += 5 * win->nelem;
    if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
        MUTIL_ERROR( "User Interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    if( percent == 0 ) {
        MUTIL_TRACE( "sigdbl_window_raised_cosine() done" );
        return MUTIL_ERR_OK;
    }

    /* rest of the window is raised cosine */

    scale = floor( percent * win->nelem );
    low   = floor( scale / 2.0 );
    angle = 2 * MUTIL_PI / ( scale + 1 );

    for( i = 0; i < low; i++ ) {
        win->data[i] = 0.5 * ( 1.0 - cos( angle * ( i + 1 ) ) );
        win->data[ win->nelem - i - 1 ] = win->data[i];
    }

    num_ops += 10 * win->nelem;
    if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
        MUTIL_ERROR( "User Interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    MUTIL_TRACE( "sigdbl_window_raised_cosine() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by James Pitton */
/* Adapted to mutils library by Luca Cazzanti */
mutil_errcode sigdbl_window_hanning( void *intrp_ptr, double_mat *win )
{
    mutil_errcode errcode;

    double angle;
    sint32 i;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start sigdbl_window_hanning()" );

    /* check matrix */

    errcode = matdbl_validate(win);
    if(errcode) return(errcode);

    /* matrix must be a vector */

    if( !MATANY_IS_VEC(win) ) {
        MUTIL_ERROR( "Input matrix must be a vector" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* build window */

    angle = 2 * MUTIL_PI / ( win->nelem + 1 );
    for( i = 0; i < win->nelem; i++ ) {
        win->data[i] = 0.5 * ( 1.0 - cos( angle * ( i + 1 ) ) );
    }

    num_ops += 50 * win->nelem;
    if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
        MUTIL_ERROR( "User Interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    MUTIL_TRACE( "sigdbl_window_hanning() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by James Pitton */
/* Adapted for mutils library by Luca Cazzanti */
mutil_errcode sigdbl_window_hamming( void *intrp_ptr, double_mat *win )
{
    mutil_errcode errcode;

    double angle;
    sint32 i;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start sigdbl_window_hamming()" );

    /* check matrix */

    errcode = matdbl_validate(win);
    if(errcode) return(errcode);

    /* matrix must be a vector */

    if( !MATANY_IS_VEC(win) ) {
        MUTIL_ERROR( "Input matrix must be a vector" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* build window */

    angle = 2 * MUTIL_PI / ( win->nelem - 1 );
    for( i = 0; i < win->nelem; i++ ) {
        win->data[i] = 0.54 - 0.46 * cos( angle * i );
    }

    num_ops += 50 * win->nelem;
    if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
        MUTIL_ERROR( "User Interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    MUTIL_TRACE( "sigdbl_window_hamming() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Writen by James Pitton */
/* Adapted to mutils library by Luca Cazzanti */
mutil_errcode sigdbl_window_blackman( void *intrp_ptr, double_mat *win )
{
    mutil_errcode errcode;

    double angle;
    sint32 i;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start sigdbl_window_hamming()" );

    /* check matrix */

    errcode = matdbl_validate(win);
    if(errcode) return(errcode);

    /* matrix must be a vector */

    if( !MATANY_IS_VEC(win) ) {
        MUTIL_ERROR( "Input matrix must be a vector" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* build window */

    angle = 2 * MUTIL_PI / ( win->nelem + 1 );
    for(i = 0; i < win->nelem; i++) {
        win->data[i] = 0.42 - 0.5 * cos( angle * ( i + 1 ) ) +
            0.08 * cos( 2 * angle * ( i + 1 ) );
    }

    num_ops += 100 * win->nelem;
    if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
        MUTIL_ERROR( "User Interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    MUTIL_TRACE( "sigdbl_window_blackman() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by James Pitton */
/* Adapted for mutils library by Luca Cazzanti */
mutil_errcode sigdbl_window_nuttall( void *intrp_ptr,double_mat *win )
{
    mutil_errcode errcode;

    sint32 i;
    sint32 len1; /* length of window minus one */
    double arg;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start sigdbl_window_nuttall()" );

    /* check matrix */

    errcode = matdbl_validate(win);
    if(errcode) return(errcode);

    /* matrix must be a vector */

    if( !MATANY_IS_VEC(win) ) {
        MUTIL_ERROR( "Input matrix must be a vector" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* build window */

    len1              = win->nelem - 1;
    win->data[0]      = 0.0;
    win->data[ len1 ] = 0.0;
    for ( i = 1; i < len1; i++ ) {
        arg = 2.0 * MUTIL_PI * (double) i / (double)( len1 );
        win->data[i] =   0.3635819
            - 0.4891775 * cos( arg )
            + 0.1365995 * cos( arg + arg )
            - 0.0106411 * cos( arg + arg + arg );
    }

    num_ops += 200 * win->nelem;
    if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
        MUTIL_ERROR( "User Interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    MUTIL_TRACE( "sigdbl_window_nuttall() done" );

    return MUTIL_ERR_OK;
}

/* Function decumented in sig_win.h */
/* Written by James Pitton */
/* Adapted for mutils library by Luca Cazzanti */
mutil_errcode sigdbl_window_gaussian ( double alpha, void *intrp_ptr,
    double_mat *win )
{
    mutil_errcode errcode;

    sint32 i;
    sint32 half_len;
    sint32 win_len;
    double arg;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start sigdbl_window_gaussian()" );

    /* check matrix */

    errcode = matdbl_validate(win);
    if(errcode) return(errcode);

    /* matrix must be a vector */

    if( !MATANY_IS_VEC(win) ) {
        MUTIL_ERROR( "Input matrix must be a vector" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* non-positive alphas are not allowed */

    if( alpha <= 0.0) {
        MUTIL_ERROR( "Alpha parameter must be positive" );
        return MUTIL_ERR_ILLEGAL_VALUE;
    }

    /* build window */

    win_len  = win->nelem;
    half_len = win_len / 2;
    for (i = 0; i < half_len; i++) {
        arg = 2.0 * alpha * ( 0.5 - (double) i / (double) ( win_len - 1 ) );
        arg                          = -arg * arg / 2.0;
        win->data[i]                 = exp( arg );
        win->data[ win_len - i - 1 ] = win->data[i];
    }

    if ( half_len + half_len != win_len) {
        win->data[ half_len ] = 1.0;
    }

    num_ops += 100 * win_len;
    if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
        MUTIL_ERROR( "User Interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    MUTIL_TRACE( "sigdbl_window_gaussian() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by William Constantine   */
mutil_errcode sigdbl_window_kaiser( double beta, void *intrp_ptr, double_mat *win )
{
  double         alpha;
  double         alphanorm;
  double         norm = beta;
  double        *pd_taper;
  mutil_errcode  err;
  sint32         i;
  sint32         nelem;

  double num_ops = 0;

  MUTIL_INTERRUPT_INIT(intrp_ptr);

  MUTIL_TRACE( "Start sigdbl_window_kaiser()" );

  /* check matrix */

  err = matdbl_validate( win );
  if( err ) return( err );

  /* matrix must be a vector */

  if( !MATANY_IS_VEC(win) ) {
    MUTIL_ERROR( "Input matrix must be a vector" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* resolution muast be positive */

  if( beta < 0.0 ) {
    MUTIL_ERROR( "Kaiser shape parameter (beta) must be positive" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* build window */

  nelem = win->nelem;
  alphanorm = (double) (nelem * nelem);

  pd_taper = win->data;

  for ( i = 0; i < nelem; i++ ){

    alpha = (double) (2 * i + 1 - nelem);

    *pd_taper = beta * sqrt( 1.0 - alpha * alpha / alphanorm );

    pd_taper++;
  }

  num_ops += 200 * win->nelem;
  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
    MUTIL_ERROR( "User Interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  err = localfn_sigdbl_besselI0( win->data, nelem, intrp_ptr );
  if ( err ) return err;

  /* normalize */

  /* normalizing factor */

  err = localfn_sigdbl_besselI0( &norm, 1, intrp_ptr );
  if( err ) return err;

  pd_taper = win->data;

  for ( i = 0; i < nelem; i++ ){

    (*pd_taper) /= norm;

    pd_taper++;
  }

  num_ops += 10 * nelem;

  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
    MUTIL_ERROR( "User Interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE( "sigdbl_window_kaiser() done" );

  return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* William Constantine */
mutil_errcode sigdbl_window_chebyshev( double sidelobe, void *intrp_ptr,
  double_mat *win )
{
  dcomplex       z;
  dcomplex      *pd_z;
  dcomplex_mat   B;
  double         M;
  double         R;
  double         fac;
  double         num_ops = 0;
  double         x0;
  double        *pd_taper;
  double        *pd_w;
  double        *pd_x;
  double_mat     w;
  double_mat     x;
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         imax;
  sint32         j;
  sint32         jmax;
  sint32         nelem;
  sint32         rem;
  double         absw;
  double         sum;
  double         norm;
  double         acoshR;

  MUTIL_INTERRUPT_INIT(intrp_ptr);

  MUTIL_TRACE( "Start sigdbl_window_chebyshev()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check matrix */

  err = matdbl_validate(win);
  if(err) return(err);

  /* matrix must be a vector */

  if( !MATANY_IS_VEC(win) ) {
    MUTIL_ERROR( "Input matrix must be a vector" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( sidelobe < 0.0 ){
    MUTIL_ERROR( "Sidelobe level must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  nelem = win->nelem;

  M = (double) nelem - 1;

  rem  = nelem % 2;
  imax = (sint32) floor( (double) ( nelem + 1 ) / 2.0 );
  jmax = (sint32) ( nelem - rem ) / 2;

  /* convert sidelobe level from dB to linear scale */

  R = MUTIL_POW( 10.0, sidelobe / 20.0 );

  err = localfn_acosh( R, &acoshR );
  if ( err ) return err;

  x0 = cosh( acoshR / M );

  /* allocate memory */

  err = matdbl_malloc_register( &x, 1, nelem, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matcpx_malloc_register( &B, 1, nelem, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &w, 1, imax, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  pd_x = x.data;

  for ( i = 0; i < nelem; i++ ){

    *pd_x = x0 * cos( MUTIL_PI * (double) ( i + 1 ) / (double) nelem );

    pd_x++;
  }

  /* develop filter transfer function */

  pd_x = x.data;
  pd_z = B.data;

  for ( i = 0; i < nelem; i++ ){

    if ( MUTIL_ABS( *pd_x ) > 1.0 ){

      z.re = *pd_x;
      z.im = 0.0;

      err = localfn_acosh_complex( &z );
      MEMLIST_FREE_ON_ERROR( err, &list );

      (z.re) *= M;
      (z.im) *= M;

      err = localfn_cosh_complex( &z );
      MEMLIST_FREE_ON_ERROR( err, &list );

      pd_z->re = z.re /= R;
      pd_z->im = z.im /= R;
    }
    else{

      pd_z->re = cos( M * acos( *pd_x ) ) / R;
      pd_z->im = 0.0;
    }

    /* increment pointers */

    pd_x++;
    pd_z++;
  }

  /* get coefficients via DFT */

  fac = MUTIL_PI / (double) nelem;

  pd_w = w.data;

  norm = 0.0;

  for ( i = 1; i <= imax; i++ ){

    *pd_w = 1.0;

    pd_z = B.data;

    for ( j = 1; j <= jmax; j++ ){

      (*pd_w) += pd_z->re * 2.0 *
	cos( fac * ( 2.0 * (double) ( i - 1 ) +  (double) ( 1.0 - rem ) ) * (double) j );

      pd_z++;
    }

    absw = MUTIL_ABS( (*pd_w) );

    if ( absw > norm ) norm = absw;

    pd_w++;
  }

  if ( norm > 0.0 ){

    err = matdbl_divide_scalar( &w, norm, TRUE, intrp_ptr, &w );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* set pointer to last element of w vector */

  pd_w     = w.data + imax - 1;
  pd_taper = win->data;

  for ( i = 0; i < imax - rem; i++ ){

    *pd_taper = *pd_w;

    pd_taper++;
    pd_w--;
  }

  /* set pointer to first element of w vector */

  pd_w = w.data;

  for ( i = 0; i < imax; i++ ){

    *pd_taper = *pd_w;

    pd_taper++;
    pd_w++;
  }

  /* normalize */

  err = matdbl_sum( win, intrp_ptr, &sum );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( sum > 0.0 ){

    err = matdbl_divide_scalar( win, sum, TRUE, intrp_ptr, win );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  num_ops += 5 * nelem;
  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
    MUTIL_ERROR( "User Interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  /* free other malloced space and
     corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "sigdbl_window_chebyshev() done" );

  return MUTIL_ERR_OK;
}


/* Function documented in sig_win.h */
/* Written Luca Cazzanti from old code by James Pitton */
mutil_errcode sigdbl_window_born_jordan( void *intrp_ptr, double_mat *win )
{
    mutil_errcode errcode;

    sint32 i;
    sint32 mid_pnt;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start sigdbl_window_born_jordan()" );

    /* check matrix */

    errcode = matdbl_validate(win);
    if(errcode) return(errcode);

    /* matrix must be a vector */

    if( !MATANY_IS_VEC(win) ) {
        MUTIL_ERROR( "Input matrix must be a vector" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* build window */
    /* Born-Jordan window is symmetric, so calculate only half */
    /* the window and assign rest by symmetry */

    mid_pnt = ( win->nelem - 1 ) / 2;

    for( i = 0; i <= mid_pnt ; i++ ) {
        win->data[i] = ( 1.0 / ( mid_pnt - i + 1.0 ) );
        win->data[ win->nelem - i - 1 ] = win->data[i];
    }

    num_ops += 5 * win->nelem;
    if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
        MUTIL_ERROR( "User Interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    MUTIL_TRACE( "sigdbl_window_born_jordan() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by William Constantine   */
mutil_errcode sigdbl_window_sinusoidal( void *intrp_ptr, double_mat *win )
{
  double        *pd_taper;
  mutil_errcode  err;
  sint32         t;
  sint32         k;
  sint32         ntaper;
  sint32         nsample;
  double         fac;
  double         amp;
  double         num_ops = 0;

  MUTIL_INTERRUPT_INIT(intrp_ptr);

  MUTIL_TRACE( "Start sigdbl_window_sinusoidal()" );

  /* check matrix */

  err = matdbl_validate(win);
  if ( err ) return( err );

  /* create the taper */

  nsample  = win->ncol;
  ntaper   = win->nrow;
  pd_taper = win->data;

  if ( ntaper > nsample ){
    MUTIL_ERROR( "Number of tapers is limited to the number of samples in each taper" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  amp = sqrt( 2.0 / (double) ( nsample + 1 ) );
  fac = MUTIL_PI / (double) ( nsample + 1 );

  for ( k = 1; k <= ntaper; k++ ){

    for ( t = 1; t <= nsample; t++ ){

      *pd_taper = amp * sin( fac * (double) ( t * k ) );

      pd_taper++;
    }
  }

  num_ops += 5 * win->nelem;
  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
    MUTIL_ERROR( "User Interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE( "sigdbl_window_sinusoidal() done" );

  return MUTIL_ERR_OK;
}


/* Function documented in sig_win.h */
/* Written by William Constantine   */
mutil_errcode sigdbl_window_parzen( sint32 cutoff, void *intrp_ptr, double_mat *win )
{
  double         fac;
  double        *pd_taper;
  mutil_errcode  err;
  sint32         i;
  double         num_ops = 0;

  MUTIL_INTERRUPT_INIT(intrp_ptr);

  MUTIL_TRACE( "Start sigdbl_window_parzen()" );

  /* check matrix */

  err = matdbl_validate(win);
  if(err) return(err);

  /* matrix must be a vector */

  if( !MATANY_IS_VEC(win) ) {
    MUTIL_ERROR( "Input matrix must be a vector" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( cutoff < 1 || cutoff > win->nelem ){
    MUTIL_ERROR( "Parzen cutoff must be positive and less than or equal to the sample size." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* create the taper */

  pd_taper = win->data;

  for ( i = 0; i < win->nelem; i++ ){

    fac = (double) i / (double) cutoff;

    if ( i <= cutoff / 2 ){

      *pd_taper = 1.0 - 6.0 * fac * fac * ( 1.0 - fac );
    }
    else if ( i > cutoff / 2 && i <= cutoff ){

      *pd_taper = 2.0 * MUTIL_POW( 1.0 - fac, 3.0 );
    }
    else *pd_taper = 0.0;

    pd_taper++;
  }

  num_ops += 5 * win->nelem;
  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
    MUTIL_ERROR( "User Interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE( "sigdbl_window_parzen() done" );

  return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by William Constantine   */
mutil_errcode sigdbl_window_papoulis( sint32 cutoff, void *intrp_ptr, double_mat *win )
{
  double         fac;
  double         fac2;
  double         amp;
  double        *pd_taper;
  mutil_errcode  err;
  sint32         i;
  double         num_ops = 0;

  MUTIL_INTERRUPT_INIT(intrp_ptr);

  MUTIL_TRACE( "Start sigdbl_window_papoulis()" );

  /* check matrix */

  err = matdbl_validate(win);
  if(err) return(err);

  /* matrix must be a vector */

  if( !MATANY_IS_VEC(win) ) {
    MUTIL_ERROR( "Input matrix must be a vector" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( cutoff < 1 || cutoff > win->nelem ){
    MUTIL_ERROR( "Papoulis cutoff must be positive and less than or equal to the sample size." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* create the taper */

  pd_taper = win->data;

  amp = 1.0 / MUTIL_PI;

  for ( i = 0; i < win->nelem; i++ ){

    fac  = (double) i / (double) cutoff;
    fac2 = MUTIL_PI * fac;

    if ( i < cutoff ){

      *pd_taper = amp * MUTIL_ABS( sin( fac2 ) ) + ( 1.0 - fac ) * cos( fac2 );
    }
    else *pd_taper = 0.0;

    pd_taper++;
  }

  num_ops += 5 * win->nelem;
  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
    MUTIL_ERROR( "User Interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE( "sigdbl_window_papoulis() done" );

  return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by William Constantine   */
mutil_errcode sigdbl_window_daniell( double roughness, void *intrp_ptr, double_mat *win )
{
  double         fac;
  double         fac2;
  double        *pd_taper;
  mutil_errcode  err;
  sint32         i;
  double         num_ops = 0;

  MUTIL_INTERRUPT_INIT(intrp_ptr);

  MUTIL_TRACE( "Start sigdbl_window_daniell()" );

  /* check matrix */

  err = matdbl_validate(win);
  if(err) return(err);

  /* matrix must be a vector */

  if( !MATANY_IS_VEC(win) ) {
    MUTIL_ERROR( "Input matrix must be a vector" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( roughness <= 0.0 ){
    MUTIL_ERROR( "Daniell roughness factor must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* create the taper */

  pd_taper = win->data;

  fac = MUTIL_PI / roughness;

  *pd_taper = 1.0;

  pd_taper++;

  for ( i = 1; i < win->nelem; i++ ){

    fac2 = fac * (double) i;

    *pd_taper = sin( fac2 ) / fac2;

    pd_taper++;
  }

  num_ops += 5 * win->nelem;
  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
    MUTIL_ERROR( "User Interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }

  MUTIL_TRACE( "sigdbl_window_daniell() done" );

  return MUTIL_ERR_OK;
}



/*
**************************
Universal Matrix Functions
**************************
*/

/* Function documented in sig_win.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_window_rectangle( void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_rectangle()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_rectangle(intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_rectangle() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_window_triangle( void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_triangle()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_triangle(intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_triangle() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_window_raised_cosine( double percent, void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_raised_cosine()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_raised_cosine( percent, intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_raised_cosine() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_window_hanning( void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_hanning()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_hanning( intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_hanning() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_window_hamming( void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_hamming()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_hamming(intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_hamming() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_window_blackman( void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_blackman()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_blackman(intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_blackman() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_window_nuttall( void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_nuttall()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_nuttall( intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
    }

    MUTIL_TRACE( "siguniv_window_nuttall() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_window_gaussian( double alpha, void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_gaussian()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_gaussian( alpha, intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_gaussian() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_window_kaiser( double res, void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_kaiser()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_kaiser( res, intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_kaiser() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_window_chebyshev( double ripple, void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_chebyshev()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_chebyshev( ripple, intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_chebyshev() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_window_born_jordan( void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_born_jordan()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_born_jordan(intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_born_jordan() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by William Constantine   */
mutil_errcode siguniv_window_sinusoidal( void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_sinusoidal()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_sinusoidal(intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_sinusoidal() done" );

    return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by William Constantine   */
mutil_errcode siguniv_window_parzen( sint32 cutoff, void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_parzen()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_parzen( cutoff, intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_parzen() done" );

    return MUTIL_ERR_OK;
}


/* Function documented in sig_win.h */
/* Written by William Constantine   */
mutil_errcode siguniv_window_papoulis( sint32 cutoff, void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_papoulis()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_papoulis( cutoff, intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_papoulis() done" );

    return MUTIL_ERR_OK;
}


/* Function documented in sig_win.h */
/* Written by William Constantine   */
mutil_errcode siguniv_window_daniell( double roughness, void *intrp_ptr, univ_mat *win )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_window_daniell()" );

    if( !win ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    ASSERT_LEGAL_TYPE(win->type);

    switch(win->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_window_daniell( roughness, intrp_ptr, &(win->mat.dblmat) );
            if(errcode) return errcode;
            break;

        default:
            MUTIL_ERROR( "This matrix type not yet supported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_window_daniell() done" );

    return MUTIL_ERR_OK;
}


/* Obtain taper and possibly normalize. */
/*                                      */
/* Documented in sig_win.h              */
/* Written by William Constantine       */

mutil_errcode sigdbl_taper(
  const sig_taper_type  taper,
  const sint32          nrow,
  const sint32          ncol,
  const double          param,
  const boolean         normalize,
  void                 *intrp_ptr,
  double_mat           *result )
{
  memlist        list;
  mutil_errcode  err;

  MUTIL_TRACE( "Start sigdbl_taper()" );

  if( !result ) {
    MUTIL_ERROR( "NULL pointer for input matrix" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  switch( taper ){
    case SIG_TAPER_RECTANGULAR:
    case SIG_TAPER_TRIANGLE:
    case SIG_TAPER_RAISED_COSINE:
    case SIG_TAPER_HANNING:
    case SIG_TAPER_HAMMING:
    case SIG_TAPER_BLACKMAN:
    case SIG_TAPER_NUTTALL:
    case SIG_TAPER_GAUSSIAN:
    case SIG_TAPER_KAISER:
    case SIG_TAPER_CHEBYSHEV:
    case SIG_TAPER_BORN_JORDAN:
    case SIG_TAPER_SINUSOIDAL:
    case SIG_TAPER_PARZEN:
    case SIG_TAPER_PAPOULIS:
    case SIG_TAPER_DANIELL:
      break;
    default:
      MUTIL_ERROR( "Taper type is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* with the exception of multitapers,
     insist that at least one dimension be unity */

  if ( taper != SIG_TAPER_SINUSOIDAL ){

    if ( nrow != 1 && ncol != 1 ){
      MUTIL_ERROR( "One of the taper dimensions must be unity." );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
    }
  }

  /* allocate memory */

  err = matdbl_malloc_register( result, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  switch( taper ){
    case SIG_TAPER_RECTANGULAR:
      err = sigdbl_window_rectangle( intrp_ptr, result );
      break;
    case SIG_TAPER_TRIANGLE:
      err = sigdbl_window_triangle( intrp_ptr, result );
      break;
    case SIG_TAPER_RAISED_COSINE:
      err = sigdbl_window_raised_cosine( param, intrp_ptr, result );
      break;
    case SIG_TAPER_HANNING:
      err = sigdbl_window_hanning( intrp_ptr, result );
      break;
    case SIG_TAPER_HAMMING:
      err = sigdbl_window_hamming( intrp_ptr, result );
      break;
    case SIG_TAPER_BLACKMAN:
      err = sigdbl_window_blackman( intrp_ptr, result );
      break;
    case SIG_TAPER_NUTTALL:
      err = sigdbl_window_nuttall( intrp_ptr, result );
      break;
    case SIG_TAPER_GAUSSIAN:
      err = sigdbl_window_gaussian( param, intrp_ptr, result );
      break;
    case SIG_TAPER_KAISER:
      err = sigdbl_window_kaiser( param, intrp_ptr, result );
      break;
    case SIG_TAPER_CHEBYSHEV:
      err = sigdbl_window_chebyshev( param, intrp_ptr, result );
      break;
    case SIG_TAPER_BORN_JORDAN:
      err = sigdbl_window_born_jordan( intrp_ptr, result );
      break;
    case SIG_TAPER_SINUSOIDAL:
      err = sigdbl_window_sinusoidal( intrp_ptr, result );
      break;
    case SIG_TAPER_PARZEN:
      err = sigdbl_window_parzen( (sint32) param, intrp_ptr, result );
      break;
    case SIG_TAPER_PAPOULIS:
      err = sigdbl_window_papoulis( (sint32) param, intrp_ptr, result );
      break;
    case SIG_TAPER_DANIELL:
      err = sigdbl_window_daniell( param, intrp_ptr, result );
      break;
    default:
      MUTIL_ERROR( "Taper type is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  MEMLIST_FREE_ON_ERROR( err, &list );

  if ( normalize ){

    switch( taper ){
      case SIG_TAPER_PARZEN:
      case SIG_TAPER_PAPOULIS:
      case SIG_TAPER_DANIELL:
	err = localfn_lag_window_normalize( result );
	break;

      default:
	err = localfn_taper_normalize( result );
	break;
    }

    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* free nodes corresponding to registered
     memory for the result, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
     corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with sigdbl_taper()" );

  return MUTIL_ERR_OK;
}

/* Function documented in sig_win.h */
/* Written by William Constantine   */
mutil_errcode siguniv_taper(
  const sig_taper_type  taper,
  const sint32          nrow,
  const sint32          ncol,
  const double          param,
  const boolean         normalize,
  void                 *intrp_ptr,
  univ_mat             *result )
{
  mutil_errcode err;

  MUTIL_TRACE( "Start frauniv_taper()" );

  if( !result ) {
    MUTIL_ERROR( "NULL pointer for input matrix" );
    return MUTIL_ERR_NULL_POINTER;
  }

  result->type = MUTIL_DOUBLE;

  err = sigdbl_taper( taper, nrow, ncol, param, normalize,
    intrp_ptr, &(result->mat.dblmat) );
  if( err ) return err;

  MUTIL_TRACE( "Done with frauniv_taper()" );

  return MUTIL_ERR_OK;
}


/* STATIC FUNCTION DEFINITIONS */

/**  Normalizes a taper to have unit energy
 * Normalizes the taper so that the sum of the squares of the
 * taper coefficients is unity.
 *
 * @usage #localfn_taper_normalize( &taper );#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param taper Pointer to array of doubles.
 * @private
*/
static mutil_errcode localfn_taper_normalize( double_mat *taper )
{
  sint32         i;
  sint32         j;
  double         norm = 0.0;
  double        *pd_taper;
  double_mat     temp;
  mutil_errcode  err;

  MUTIL_TRACE( "Start localfn_taper_normalize()" );

  /* if multitaper case, normalize each row (each taper)
     of the input */

  if ( taper->nrow > 1 && taper->ncol > 1 ){

    temp.data  = taper->data;
    temp.nelem = taper->ncol;
    temp.ncol  = taper->ncol;
    temp.nrow  = 1;

    /* recursively call the normalization function,
       one call for each row */

    for ( j = 0; j < taper->nrow; j++ ){

      err = localfn_taper_normalize( &temp );
      if ( err ) return err;

      /* update data pointer */

      temp.data += taper->ncol;
    }
  }
  else{

    pd_taper = taper->data;

    for ( i = 0; i < taper->nelem; i++ ){

      norm += (*pd_taper) * (*pd_taper);

      pd_taper++;
    }

    norm = 1.0 / sqrt( norm );

    /* normalize the taper to have unit energy */

    pd_taper = taper->data;

    for ( i = 0; i < taper->nelem; i++ ){

      (*pd_taper) *= norm;

      pd_taper++;
    }
  }
  MUTIL_TRACE( "Done with localfn_taper_normalize()" );

  return MUTIL_ERR_OK;
}



/**  Normalizes a lag window so that the zeroth lag is unity.
 *
 * @usage #localfn_lag_window_normalize( &taper );#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source sig\_win.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param taper Pointer to array of doubles.
 * @private
*/
static mutil_errcode localfn_lag_window_normalize( double_mat *taper )
{
  sint32         i;
  sint32         j;
  double        *pd_taper;
  double         norm;
  double_mat     temp;
  mutil_errcode  err;

  MUTIL_TRACE( "Start localfn_lag_window_normalize()" );

  /* if multitaper case, normalize each row (each taper)
     of the input */

  if ( taper->nrow > 1 && taper->ncol > 1 ){

    temp.data  = taper->data;
    temp.nelem = taper->ncol;
    temp.ncol  = taper->ncol;
    temp.nrow  = 1;

    /* recursively call the normalization function,
       one call for each row */

    for ( j = 0; j < taper->nrow; j++ ){

      err = localfn_taper_normalize( &temp );
      if ( err ) return err;

      /* update data pointer */

      temp.data += taper->ncol;
    }
  }

  /* normalize the taper to have unit energy at lag zero */

  pd_taper = taper->data;

  norm = 1.0 / *pd_taper;

  for ( i = 0; i < taper->nelem; i++ ){

    (*pd_taper) *= norm;

    pd_taper++;
  }

  MUTIL_TRACE( "Done with localfn_lag_window_normalize()" );

  return MUTIL_ERR_OK;
}


/* code set for acosh(z) where z is a complex vector */

static mutil_errcode localfn_acosh_complex( dcomplex *z )
{
  boolean        lower_half;
  dcomplex       I;
  mutil_errcode  err;

  I.re = 0.0;
  I.im = 1.0;

  lower_half = (boolean) ( z->im < 0.0 );

  err = localfn_acos_complex( z );
  if ( err ) return err;

  LOCALDEF_COMPLEX_MULTIPLY( z, z, &I );

  if( lower_half ){
    z->re = -z->re;
    z->im = -z->im;
  }

  return MUTIL_ERR_OK;
}

/* code set for cosh(z) where z is a complex number */

static mutil_errcode localfn_cosh_complex( dcomplex *z )
{
  double         dExpLimit = log( MUTIL_DOUBLE_MAX ) * ( 1 + MUTIL_DOUBLE_EPSILON );
  double         dTrigLimit = 0.1 / MUTIL_DOUBLE_EPSILON;
  double         im = z->im;
  double         re = z->re;

  if( im < -dTrigLimit || im > dTrigLimit ){
    MUTIL_ERROR("Uncalculable complex cosh() value.");
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if( MUTIL_ABS( re ) > dExpLimit ) {
    MUTIL_ERROR("Infinite complex cosh() value.");
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  z->re = cosh( re ) * cos( im );

  z->im = sinh( re ) * sin( im );

  return MUTIL_ERR_OK;
}


static mutil_errcode localfn_acos_complex( dcomplex *z )
{
  double re = z->re;
  double im = z->im;
  double x;
  double y;
  double alpha;
  double beta;

  x = localfn_hypot(re+1, im);
  y = localfn_hypot(re-1, im);

  alpha = ( x + y ) / 2.0;
  beta  = ( x - y ) / 2.0;

  /* beta should be <= 1 in absolute value: force for acos() */

  if( beta > 1.0 ){
    beta = 1.0;
  }
  else if( beta < -1 ){
    beta = -1;
  }

  z->re = acos( beta );
  z->im = -log( alpha + sqrt( alpha * alpha - 1 ) );

  if( im < 0 || ( im == 0 && re >= 1 ) ){
    z->im = -z->im;
  }

  return MUTIL_ERR_OK;
}

/* Windows does not have acosh in any standard library
   so , here it is
*/

static mutil_errcode localfn_acosh( double x, double *result )
{
  double z;

  if( x < 1.0 ){

    MUTIL_ERROR( "acosh argument must be < 1.0" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if( x > 1500.0 ){
    *result = log( x ) + 0.69314718055994529;
  }

  z = x - 1.0;

  if ( z < 0.5 ){
    z = (((( 0.0017596881071 * z
	  - 0.0075272886713) * z
	  + 0.026454905019) * z
	  - 0.11784741703) * z
	  + 1.4142135263) * sqrt( z );
  }
  else{
    z = sqrt( z * ( x + 1.0 ) );
    z = log( x  +  z );
  }

  *result = z;

  return MUTIL_ERR_OK;
}

/* hypotenuse function */

static double localfn_hypot( double x, double y )
{
  return( (double) sqrt( x * x + y * y ) );
} 



