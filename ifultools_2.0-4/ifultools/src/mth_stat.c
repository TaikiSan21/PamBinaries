
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mth_stat.c $: $Revision: #1 $, $Date: 2008/03/21 $";

#include "mth_stat.h"
#include "mat_stat.h"
#include "mat_type.h"
#include "ut_debug.h"
#include "ut_limit.h"
#include "ut_math.h"
#include <math.h>

/* define macros */

#undef LOCALDEF_CHECK_NULL_MATRIX_POINTER
#define LOCALDEF_CHECK_NULL_MATRIX_POINTER( DATA_PTR, DATA_TYPE,  \
                                     TYPE_PREFIX )                \
   err = TYPE_PREFIX ## _validate( DATA_PTR );                    \
   if ( err ) return err;                                         \
   if ( DATA_PTR == ( DATA_TYPE * ) NULL ){                       \
     MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" );    \
     return MUTIL_ERR_NULL_POINTER;                               \
   }

#define SQRTPI 1.7724538509055160272981674833411451827976

#define DPOLYD( y, p, q ) \
 for ( n = d = 0, i = sizeof( p ) / sizeof( p[ 0 ] ); --i >= 0; ) \
   { n = n * y + p[i]; d = d * y + q[ i ]; }

/* define static functions */


/* define static arrays used in calculating the erf() and erfc()
   functions */

static double p1[] = {
  2.4266795523053175e2,
  2.1979261618294152e1,
  6.9963834886191355e0,
  -3.5609843701815385e-2
};

static double q1[]  = {
  2.1505887586986120e2,
  9.1164905404514901e1,
  1.5082797630407787e1,
  1.0
};

static double p2[]  = {
  3.004592610201616005e2,
  4.519189537118729422e2,
  3.393208167343436870e2,
  1.529892850469404039e2,
  4.316222722205673530e1,
  7.211758250883093659e0,
  5.641955174789739711e-1,
  -1.368648573827167067e-7
};

static double q2[]  = {
  3.004592609569832933e2,
  7.909509253278980272e2,
  9.313540948506096211e2,
  6.389802644656311665e2,
  2.775854447439876434e2,
  7.700015293522947295e1,
  1.278272731962942351e1,
  1.0
};

static double p3[] = {
  -2.99610707703542174e-3,
  -4.94730910623250734e-2,
  -2.26956593539686930e-1,
  -2.78661308609647788e-1,
  -2.23192459734184686e-2
};

static double q3[] = {
  1.06209230528467918e-2,
  1.91308926107829841e-1,
  1.05167510706793207e0,
  1.98733201817135256e0,
  1.0
};

/* The gamma function.               */
/* Function documented in mth_stat.h */
/* Written by William Constantine    */

double mth_gamma( double x )
{
  double w;
  double y;
  sint32 k;
  sint32 n;

  /* avoid lint warning */

  (void) whatssi;

  MUTIL_TRACE( "Start mth_gamma()" );

  /* check for positive value of input */

  if ( x < 0.0 ){
    MUTIL_ERROR("Input must be positive");
    return ((double) 0.0);
  }

  n = x < 1.5 ? -((sint32) floor(2.5 - x)) : (sint32) floor(x - 1.5);
  w = x - (n + 2);
  y = ((((((((((((-1.99542863674e-7 * w + 1.337767384067e-6) * w -
    2.591225267689e-6) * w - 1.7545539395205e-5) * w +
    1.45596568617526e-4) * w - 3.60837876648255e-4) * w -
    8.04329819255744e-4) * w + 0.008023273027855346) * w -
    0.017645244547851414) * w - 0.024552490005641278) * w +
    0.19109110138763841) * w - 0.233093736421782878) * w -
    0.422784335098466784) * w + 0.99999999999999999;

  if ( n > 0 ){
    w = x - 1;

    for ( k = 2; k <= n; k++ ){
      w *= x - k;
    }
  }
  else{
    w = 1;

    for ( k = 0; k > n; k-- ){
      y *= x - k;
    }
  }

  MUTIL_TRACE( "Done with mth_gamma()" );

  return ( w / y );
}

/* The digamma function.             */
/* Function documented in mth_stat.h */
/* Written by William Constantine    */

double mth_digamma( double x )
{
  /* set constants, SN = Nth Stirling coefficient, D1 = DIGAMMA( 1.0 ) */

  double result;
  double S  = 1E-05;
  double C  = 8.5;
  double S3 = 1.0 / 12.0;
  double S4 = 1.0 / 120.0;
  double S5 = 3.968253968E-03;
  double D1 = -0.5772156649;
  double R;

  /* avoid lint warning */

  (void) whatssi;

  MUTIL_TRACE( "Start mth_digamma()" );

  /* initialize result */

  result = 0.0;

  /* check for positive value of input */

  if (x < 0.0){
    MUTIL_ERROR("Input must be positive");
    return result;
  }

  /* use approximation if argument <= S */

  if ( x < S ){
    result = D1 - 1.0 / x;
    return result;
  }

  /* reduce to result( X + N ) where ( X + N ) >= C */

  while (x < C){
    result -= 1.0 / x;
    x += 1.0;
  }

  /* use de Moivre's expansion */

  R       = 1.0 / x;
  result += log(x) - 0.5 * R;
  R       = R * R;
  result -= R * ( S3 - R * ( S4 - R * S5 ) );

  MUTIL_TRACE( "Done with mth_digamma()" );

  return( result );
}

/* The trigamma function.             */
/* Function documented in mth_stat.h */
/* Written by William Constantine    */

double mth_trigamma( double x )
{
  double result = 0.0;
  double a  = 1e-04;
  double b  = 5.0;
  double b2 = 1.0 / 6.0;
  double b4 = -1.0 / 30.0;
  double b6 = 1.0 / 42.0;
  double b8 = b4;
  double y;

  /* avoid lint warning */

  (void) whatssi;

  MUTIL_TRACE( "Start mth_trigamma()" );

  /* check for positive value of input */

  if ( x < 0.0 ){
    MUTIL_ERROR("Input must be positive");
    return result;
  }

  /* use small value approximation if x <= a */

  if ( x <= a ){
    result = 1.0 / ( x * x );
  }
  else{

    while ( x < b ){

      result += 1.0 / ( x * x );
      x += 1.0;
    }

    /* apply asymptotic formula */

    y       = 1.0 / ( x * x );
    result += 0.5 * y + ( 1.0 + y * ( b2 + y * ( b4 + y * ( b6 + y * b8 ) ) ) ) / x;
  }

  MUTIL_TRACE( "Done with mth_trigamma()" );

  return result;
}

/* The error function.               */
/* Function documented in mth_stat.h */
/* Written by William Constantine    */

double mth_erf( double x )
{
  double d;
  double n;
  double xsq;
  double result;
  register int i;

  /* avoid lint warning */

  (void) whatssi;

  MUTIL_TRACE( "Start mth_erf()" );

  if ( x < 0.0 ){

    result = - mth_erf( - x );
  }
  else if ( x == 0.0 ){

    result = 0.0;
  }
  else if ( x > 0.0 && x < .46875 ){

    xsq = x * x;
    DPOLYD( xsq, p1, q1 );
    result =  x * n / d;
  }
  else if ( x >= .46875 && x < 4.0 ){

    result = 1.0 - mth_erfc( x );
  }
  else {

    if ( x > sqrt( log( MUTIL_DOUBLE_MAX ) ) ){

      result =  1.0;
    }
    else{
      result = 1.0 - mth_erfc( x );
    }
  }

  MUTIL_TRACE( "Done with mth_erf()" );

  return( result );
}

/* The complementary error function   */
/* Function documented in mth_stat.h */
/* Written by William Constantine    */

double mth_erfc(double x)
{
  double d;
  double n;
  double ooxsq;
  double result;
  double xsq;
  register int i;

  /* avoid lint warning */

  (void) whatssi;

  MUTIL_TRACE( "Start mth_erfc()" );

  if ( x < 0.0 ){

    result = 2.0 - mth_erfc( -x );
  }
  else if ( x == 0.0 ){

    result = 1.0;
  }
  else if ( x > 0.0 && x < .46875 ){

    result = 1.0 - mth_erf( x );
  }
  else if ( x >= .46875 && x < 4.0 ){

    DPOLYD(x, p2, q2);
    result = exp( - x * x ) * n / d;
  }
  else{

    if ( x > sqrt( log( MUTIL_DOUBLE_MAX ) ) ){
      result = 0.0;
    }
    else{
      xsq   = x * x;
      ooxsq = 1.0 / xsq;
      DPOLYD( ooxsq, p3, q3 );
      result = ( exp( - xsq ) / x ) * ( ( 1.0 / SQRTPI ) + ( 1.0 / xsq ) * ( n / d ) );
    }
  }

  MUTIL_TRACE( "Done with mth_erfc()" );

  return( result );
}




