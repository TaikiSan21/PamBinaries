
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mth_dist.c $: $Revision: #1 $, $Date: 2008/03/21 $";

#include "mth_dist.h"
#include "mth_stat.h"
#include "mat_stat.h"
#include "mat_type.h"
#include "ut_debug.h"
#include "ut_limit.h"
#include "ut_math.h"
#include <math.h>

/* define static variables */

static double gamma_normal_approx_min_df = 1.0e8;
static double goobie = 0.9189385332046727417803297;

static double a_1[] = {
	 1.264616e-2,
	-1.425296e-2,
	 1.400483e-2,
	-5.886090e-3,
	-1.091214e-2,
	-2.304527e-2,
	 3.135411e-3,
	-2.728484e-4,
	-9.699681e-3,
	 1.316872e-2,
	 2.618914e-2,
	-2.222222e-1,
	 5.406674e-5,
	 3.483789e-5,
	-7.274761e-4,
	 3.292181e-3,
	-8.729713e-3,
	 4.714045e-1,
	 1.000000
};
static double p_1[] = {
	0.83333333333333101837e-1,
	-.277777777735865004e-2,
	0.793650576493454e-3,
	-.5951896861197e-3,
	0.83645878922e-3,
	-.1633436431e-2,
};
static double p_2[] = {
	-.42353689509744089647e5,
	-.20886861789269887364e5,
	-.87627102978521489560e4,
	-.20085274013072791214e4,
	-.43933044406002567613e3,
	-.50108693752970953015e2,
	-.67449507245925289918e1,
	0.0,
};
static double q_2[] = {
	-.42353689509744090010e5,
	-.29803853309256649932e4,
	0.99403074150827709015e4,
	-.15286072737795220248e4,
	-.49902852662143904834e3,
	0.18949823415702801641e3,
	-.23081551524580124562e2,
	0.10000000000000000000e1,
};


static double num[] = {
	2.515517,
	0.802853,
	0.010328
};

static double den[] = {
	1.000000,
	1.432788,
	0.189269,
	0.001308
};
/* define macros */

#define ROOT_2	     1.4142135623730950488016887242096980785697 /*2^1/2*/
#define ROOT_2PI     2.5066282746310005024157652848110452530070 /*2pi^1/2*/
#define MVAR	     (sizeof(p_1)/sizeof(double))
#define NVAR	     (sizeof(p_2)/sizeof(double))
#define DOUBLE_EPS   2.2204460492503131E-16
#define DOUBLE_XMIN  2.2250738585072014E-308
#define DOUBLE_XMAX  1.7976931348623157E+308
#define	IGAMMA_LARGE 1.0e30

/* declare static functions */

static double localfn_fmax(double a, double b);
static double localfn_gamma_asym(double x);
static double localfn_gamma_neg(double x);
static double localfn_gamma_pos(double x);
static double localfn_igamma(double x, double df);
static double localfn_igamma_normal(double x, double df);
static double localfn_igauss(double x);
static double localfn_invigamma(double p, double df);
static double localfn_invigauss(double p);
static double localfn_invigauss_quick(double p);
static double localfn_S_log_gamma(double x);
static double localfn_log_gamma(double x);
static double localfn_S_pnorm(double q);


/* The cumulative distribution function */
/* for a random normal variable.        */
/* Function documented in mth_dist.h    */
/* Written by William Constantine       */

double mth_pnorm( double x )
{
  double result;

  /* avoid lint warning */

  (void) whatssi;

  MUTIL_TRACE( "Start mth_pnorm()" );

  if( x == 0 ){

    result = 0.5;
  }
  else if( x > 0 ){

    result = ( 1 + mth_erf( x / ROOT_2 ) ) / 2;
  }
  else{

    result = mth_erfc( - x / ROOT_2 ) / 2;
  }

  MUTIL_TRACE( "Done with mth_pnorm()" );

  return( result );
}


/* Quantiles for a chi-square distribution */
/* with a specified number of degrees of   */
/* freedom.                                */
/* Function documented in mth_dist.h       */
/* Written by William Constantine          */

double mth_qchisq( double probability, double dof )
{
  return( 2.0 * localfn_invigamma( probability, dof / 2 ) );
}


/* Quantiles for a normal distribution     */
/* with a specified number of degrees of   */
/* freedom.                                */
/* Function documented in mth_dist.h       */
/* Written by William Constantine          */
double mth_qnorm(double probability)
{
	return( localfn_invigauss( probability ) );
}


/* define static functions */

static double localfn_invigamma(double p, double df)
{
  int i, sign;
  double q, p1, f1, f2, pdif, pabs, del, dl, pl;

  del = 0.0;
  dl  = 0.0;
  pl  = 0.0;

  if(p == 0)
    return(0.0);
  if(df < 0.5) {
    /* following is lower bound for q, very close for small p */
    q = exp( (log(p*df) + localfn_log_gamma(df))/df ) ;
  } else if(df == 0.5) {
    q = localfn_invigauss((1-p)/2);
    return(q*q/2);
  } else if(df == 1)
    return(-log(1-p));
  else {
    f1 = 0.5 / df;
    f2 = sqrt(f1) * localfn_invigauss(p);
    q = (((a_1[0]+a_1[1]*f2)*f1+(((a_1[2]+a_1[3]*f2)*f2+a_1[4])*
	    f2+a_1[5]))*f1+(((((a_1[6]+a_1[7]*f2)*f2+a_1[8])*f2+
			     a_1[9])*f2+a_1[10])*f2+a_1[11]))*f1+(((((a_1[12]*f2+
			     a_1[13])*f2+a_1[14])*f2+a_1[15])*f2+a_1[16])*f2*f2+
			       a_1[17])*f2+a_1[18];
    q = q * q * q * df;
    if(df >= 30 ||
      (df >= 15 && p > 0.49 && p <= 0.99) ||
      (p >= 0.01 && p <= 0.99 && df > 25.2-20.8*p))
      return(q);
  }
  for(i = 0; i < 50; i++) {
    p1 = localfn_igamma(q, df);
    pdif = p - p1;
    pabs = fabs(pdif);
    if(pabs <= DOUBLE_EPS)
      break;
    if(i && pabs >= pl)
      del *= -0.5;
    else {	/* del is a newton-raphson correction */
      del = log(pabs) - (df-1)*log(q) + q + localfn_log_gamma(df);
      sign = pdif >= 0 ? 1 : -1;
      del = sign * exp(del);
      if(-del > q)
	del = -0.9 * q;
      if(sign * del < 0)
	del = sign * dl / 2;
    }
    dl = fabs(del);
    pl = pabs;
    q = localfn_fmax(q+del, DOUBLE_XMIN);
    if(dl <= q*DOUBLE_EPS)
      break;
  }
  if (i==50 && pabs > 32 * DOUBLE_EPS) {
    /* pabs>32*DOUBLE_EPS is to deal with NR iterations that
     * bounce around right value but never move to center.
     * pabs>DOUBLE_EPS about 8/1000 times with random p,df.
     */

    MUTIL_ERROR( "Nonconvergence in localfn_invigamma " );
    q = - 1.0;

    /*     PROBLEM "Nonconvergence in localfn_invigamma (p=%g,df=%g)", p, df WARNING(NULL_ENTRY) ; */
    /*     na_set3(&q, DOUBLE,To_NaN); */
  }
  return(q);
}

/*
 * Inverse Gaussian distribution function.
 * Use rational approximation to get close,
 * then bisection and secant to refine.
 * Assumption: 0 < p < 1.
 */

static double localfn_invigauss(double p)
{
	int i;
	double ql, qr, qm, qdiff;
	double pl, pr, pm, pdiff;

	int dum = 1; /* used in the while(dum == 1) code below to avoid lint warnings */

	qm = 0.0;

	if(p == 0.5)
		return(0.0);

	/* initialize [ql,qr] containing the root */
	ql = qr = localfn_invigauss_quick(p);
	pl = pr = localfn_igauss(ql);
	if(pl == p)
		return(ql);
	if(pl < p)

		while(dum == 1) {
			qr += 0.001;
			pr = localfn_igauss(qr);
			if(pr == p)
				return(pr);
			if(pr > p)
				break;
		}
	else
		while(dum == 1) {
			ql -= 0.001;
			pl = localfn_igauss(ql);
			if(pl == p)
				return(pl);
			if(pl < p)
				break;
		}

	/* a few steps of bisection */
	for(i = 0; i < 5; i++) {
		qm = (ql + qr) / 2;
		pm = localfn_igauss(qm);
		qdiff = qr - ql;
		pdiff = pm - p;
		if(fabs(qdiff) < DOUBLE_EPS*qm || fabs(pdiff) < DOUBLE_EPS)
			return(qm);
		if(pdiff < 0) {
			ql = qm;
			pl = pm;
		} else {
			qr = qm;
			pr = pm;
		}
	}

	/* a few steps of secant */
	for(i = 0; i < 20; i++) {
		qm = ql + (p-pl)*(qr-ql)/(pr-pl);
		pm = localfn_igauss(qm);
		qdiff = qr - ql;
		pdiff = pm - p;
		if(fabs(qdiff) < DOUBLE_EPS*qm || fabs(pdiff) < DOUBLE_EPS)
			return(qm);
		if(pdiff < 0) {
			ql = qm;
			pl = pm;
		} else {
			qr = qm;
			pr = pm;
		}
	}

	/* no convergence */
	return(qm);
}

/*
 * Rational approximation to inverse Gaussian distribution.
 * Absolute error is bounded by 4.5e-4.
 * Reference: Abramowitz and Stegun, page 933.
 * Assumption: 0 < p < 1.
 */

static double localfn_invigauss_quick(double p)
{
	int lower;
	double t, n, d, q;

	if(p == 0.5)
		return(0.0);
	lower = p < 0.5;
	p = lower ? p : 1 - p;
	t = sqrt(-2 * log(p));
	n = (num[2]*t + num[1])*t + num[0];
	d = ((den[3]*t + den[2])*t + den[1])*t + den[0];
	q = lower ? n/d - t : t - n/d;
	return(q);
}
/*
 * Gaussian distribution function.
 */
static double localfn_igauss(double x)
{
	if(x == 0)
		return(0.5);
	else if(x > 0)
		return((1 + mth_erf(x/ROOT_2))/2);
	else
		return(mth_erfc(-x/ROOT_2)/2);
}




static double localfn_igamma_normal(double x, double df)
{
	/* normal approximation with adjustment for skewness.
	 * See Abramowitz&Stegun (1970) 26.2.48 and 26.1.32.
	 * The latter shows skewness is 2/sqrt(df), the former
	 * is general formula for the adjustment, -skewness/6*z2.
	 */
	double val ;
	double z2 ;
	x = (x - df) / sqrt(df) ;
	val = localfn_S_pnorm(x) ;
	z2 = exp(- x*x/2.0) ;
	if (z2>0) /* avoid 0*Inf */
		z2 = z2 * (x*x - 1.0) / ROOT_2PI;
	/* z2 is now 2nd derivative of normal density at x */
	val = val - (2.0/sqrt(df))/6.0 * z2 ;
	if (val < 0.0) /* happens for sufficiently small x, less so for larger df */
		return 0.0;
	return val ;
}

static double localfn_igamma(double x, double df)
{
	double factor, term, gintegral, pn[6], rn, ak, bk;
	double increment, df1;
	int i, count, k;

	if (x <= 0.0) return(0.0);

	if ( df >= MUTIL_DOUBLE_MAX ) return 0.0;

	if(df>=gamma_normal_approx_min_df)
	  return localfn_igamma_normal(x, df) ;

	if (df < 1.0) {
		increment = exp(df*log(x) - x - localfn_log_gamma(df + 1.0));
		df1 = df + 1.0;
	} else {
		increment = 0.0;
		df1 = df;
	}

	factor = exp(df1*log(x) - x - localfn_log_gamma(df1));

	if (x > 1.0 && x >= df1) {
		pn[0] = 0.0;
		pn[2] = pn[1] = 1.0;
		pn[3] = x;
		count = 1;
		rn = 1.0 / x;
		do {
			count++;
			k = count / 2;
			gintegral = rn;
			if (count%2 == 0) {
				bk = 1.0;
				ak = (double)k - df1;
			} else {
				bk = x;
				ak = (double)k;
			}
			pn[4] = bk*pn[2] + ak*pn[0];
			pn[5] = bk*pn[3] + ak*pn[1];
			rn = pn[4] / pn[5];
			for (i=0; i<4; i++)
				pn[i] = pn[i+2];
			if (pn[4] > IGAMMA_LARGE)
				for (i=0; i<4; i++)
					pn[i] /= IGAMMA_LARGE;
		} while (fabs(gintegral-rn) > DOUBLE_EPS*rn);
		gintegral = 1.0 - factor*rn;
	} else {
		gintegral = term = 1.0;
		rn = df1;
		do {
			rn += 1.0;
			term *= x/rn;
			gintegral += term;
		} while (term > DOUBLE_EPS*gintegral);
		gintegral *= factor/df1;
	}
	return(increment + gintegral);
}


/*
 * Log gamma for real argument.
 * The coefficients for the expansion around zero
 * are #5243 from Hart & Cheney; for the expansion
 * around infinity they are #5404.
 */

static double localfn_log_gamma(double x)
{
	return localfn_S_log_gamma(x) ;
}

static double localfn_S_log_gamma(double x)
{
/* 	SignGam = 1; */
	if(x <= 0)
		return(localfn_gamma_neg(x));
	if(x > 8)
		return(localfn_gamma_asym(x));
	return(log(localfn_gamma_pos(x)));
}

static double localfn_gamma_asym(double x)
{
	double n, xx;
	int i;

	xx = 1/(x*x);
	for(n=0, i=MVAR-1; i>=0; i--)
		n = n*xx + p_1[i];
	return((x-0.5)*log(x) - x + goobie + n/x);
}

static double localfn_gamma_neg(double x)
{
	double temp;

	x = -x;
	temp = sin(MUTIL_PI*x);
	if(temp == 0)
		return(DOUBLE_XMAX);
	if(temp < 0)
		temp = -temp;
/* 	else */
/* 		SignGam = -1; */
	return(-log(x*localfn_gamma_pos(x)*temp/MUTIL_PI));
}

static double localfn_gamma_pos(double x)
{
	double n, d, s;
	register int i;

	if(x < 2)
		return(localfn_gamma_pos(x+1)/x);
	if(x > 3)
		return((x-1)*localfn_gamma_pos(x-1));
	s = x - 2;
	for(n=0, d=0, i=NVAR-1; i>=0; i--) {
		n = n*s + p_2[i];
		d = d*s + q_2[i];
	}
	return(n/d);
}

static double localfn_fmax(double a, double b)
{
	return(a > b ? a : b);
}

static double localfn_S_pnorm(double q)
{
	return(localfn_igauss(q));
}
