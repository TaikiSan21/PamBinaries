
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/sig_tran.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";

/* This is a self-documenting doc++ file */

#include "sig_tran.h"

#include "mat_any.h"
#include "mat_assn.h"
#include "mat_cast.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_math.h"

/* This file contains implementations of the functions in the
   sig_tran.h header file, which are functions for signal processing
   involving transformations such as discrete cosine transform
   and Fourier transform.
*/

/* square root of 2 */

#define SQRT_2 1.414213562373095

/* helps with code writing for interrupt checking */

#define CLEAN_ALL \
   MUTIL_FREE_BUFFER_WARN( at_dummy, vec_size ); \
   MUTIL_FREE_BUFFER_WARN( ck_dummy, vec_size ); \
   MUTIL_FREE_BUFFER_WARN( bt_dummy, vec_size ); \
   MUTIL_FREE_BUFFER_WARN( sk_dummy, vec_size ); \
   MUTIL_FREE_BUFFER_WARN( np_dummy, nvec * sizeof( sint32) )


/* Static functions are declared here. The definitions and
 * documentation are at the bottom of the file.
 */

static void localfn_factorize_siglen( sint32 n, sint32 *m, sint32 *i,
    sint32 *k, sint32 *k1, sint32 *kt, sint32 *nfac, sint32 *j, sint32 *jj,
    sint32 *maxf, double *num_ops );

static mutil_errcode localfn_cpxdft( dcomplex *z, sint32 ntot, sint32 n,
    sint32 nspan, sint32 isn, void *intrp_ptr );

static mutil_errcode localfn_dctII( double_mat *sig, void *intrp_ptr);

static mutil_errcode localfn_idctII( double_mat *sig, void *intrp_ptr);

/*
***********************
Double Matrix Functions
***********************
*/


/* Function documented in sig_tran.h */
/* written by Luca Cazzanti */
mutil_errcode sigdbl_transform_discrete_cosine_II( const double_mat *sig,
    void *intrp_ptr, double_mat *result )
{
    mutil_errcode errcode;

    boolean trans_flag = FALSE; /* tells if matrix was transposed */
                                /* used for 1D vector inputs      */

    MUTIL_TRACE( "Start sigdbl_transform_discrete_cosine_II()" );

    /* avoid lint warning */
    (void) whatssi;

    /* validate input */

    errcode = matdbl_validate(sig);
    if(errcode) return(errcode);

    errcode = matdbl_validate(result);
    if(errcode) return(errcode);

    /* matrices must have same dimensions */

    if( !MATANY_EQUAL_DIM(sig, result) ) {
        MUTIL_ERROR( "Matrices must have same dimensions" );
        return(errcode);
    }

    /* copy input into result and then use in-place routines */

    errcode = matdbl_assign(sig, intrp_ptr, result);
    if(errcode) return(errcode);

    /* trivial case : signal length is 1 */

    if( sig->nelem == 1 ) {
	return MUTIL_ERR_OK;
    }

    /* 1D vectors should be in column vector form */
    /* OK to use same matrix pointer for 1D vectors in transpose */

    if(sig->nrow == 1) {
        errcode = matdbl_transpose(result, intrp_ptr, result);
        if(errcode) {
            return(errcode);
        }
        trans_flag = TRUE;
    }

    /* compute dct */

    errcode = localfn_dctII(result, intrp_ptr);
    if(errcode) return(errcode);

    /* must restore original result matrix size if it was transposed */
    /* OK to use same matrix for 1D transpose */

    if(trans_flag == TRUE) {
        errcode = matdbl_transpose(result, intrp_ptr, result);
        if(errcode) {
            return(errcode);
        }
    }

    MUTIL_TRACE( "sigdbl_transform_discrete_cosine_II() done" );
    return MUTIL_ERR_OK;

}


/* Function documented in sig_tran.h */
/* written by Luca Cazzanti */
mutil_errcode sigdbl_transform_discrete_cosine_II_inverse(
    const double_mat *sig, void *intrp_ptr, double_mat *result )
{
    mutil_errcode errcode;

    boolean trans_flag = FALSE; /* tells if matrix was transposed */
                                /* used for 1D vector inputs      */

    MUTIL_TRACE( "Start sigdbl_transform_discrete_cosine_II_inverse()" );

    /* validate input */

    errcode = matdbl_validate(sig);
    if(errcode) return(errcode);

    errcode = matdbl_validate(result);
    if(errcode) return(errcode);

    /* must have same dimensions */

    if( !MATANY_EQUAL_DIM(sig, result) ) {
        MUTIL_ERROR( "Matrices must have same dimensions" );
        return(errcode);
    }

    /* copy signal to result */

    errcode = matdbl_assign(sig, intrp_ptr, result);
    if(errcode) return(errcode);

    /* trivial case : signal length is 1 */

    if( sig->nelem == 1 ) {
	return MUTIL_ERR_OK;
    }

    /* 1D vectors should be in column vector form */
    /* OK to use same matrix pointer for 1D vectors in transpose */

    if(sig->nrow == 1) {
        errcode = matdbl_transpose(result, intrp_ptr, result);
        if(errcode) {
            return(errcode);
        }
        trans_flag = TRUE;
    }

    /* compute inverse dct */

    errcode = localfn_idctII(result, intrp_ptr);
    if(errcode) return(errcode);

    /* must restore original result matrix size if it was transposed */
    /* OK to use same matrix for 1D transpose */

    if(trans_flag == TRUE) {
        errcode = matdbl_transpose(result, intrp_ptr, result);
        if(errcode) {
            return(errcode);
        }
    }

    MUTIL_TRACE( "sigdbl_transform_discrete_cosine_II_inverse() done" );
    return(MUTIL_ERR_OK);
}


/* Function documented in sig_tran.h */
/* written by Luca Cazzanti */
mutil_errcode sigdbl_transform_discrete_cosine_II_2d( const double_mat *sig,
    boolean inverse_flag, void *intrp_ptr, double_mat *result)
{
    mutil_errcode errcode;
    double_mat    tmp_mat;

    MUTIL_TRACE( "Start sigdbl_transform_discrete_cosine_II_2d()" );

    /* validate arguments */

    errcode = matdbl_validate(sig);
    if(errcode) {
        return(errcode);
    }

    errcode = matdbl_validate(result);
    if(errcode) {
        return errcode;
    }

    /* matrices must have same dimensions */

    if( !MATANY_EQUAL_DIM(sig, result) ) {
        MUTIL_ERROR( "Matrices must have the same size" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* trivial case: signal length is 1 */

    if( sig->nelem == 1 ) {
	result->data[ 0 ] = sig->data[ 0 ];
	return MUTIL_ERR_OK;
    }

    /* make sure input is 2D */

    if( MATANY_IS_VEC(sig) ) {
        MUTIL_ERROR( "Input must be a 2D matrix" );
        return(MUTIL_ERR_ILLEGAL_SIZE);
    }

    /* need this matrix for temporary storage of transposed matrices */

    errcode = matdbl_malloc(&tmp_mat, sig->ncol, sig->nrow);
    if(errcode) return errcode;

    /* if computing the direct transform, first transform
     * the rows of the original matrix and then the columns
     * of the row-transformed matrix */

    if( inverse_flag == FALSE ) {
        errcode = matdbl_transpose(sig, intrp_ptr, &tmp_mat);
        if(errcode) {
            MUTIL_FREE_WARN(matdbl, &tmp_mat);
            return(errcode);
        }

        errcode = localfn_dctII(&tmp_mat, intrp_ptr);
        if(errcode) {
            MUTIL_FREE_WARN(matdbl, &tmp_mat);
            return(errcode);
        }

        errcode = matdbl_transpose(&tmp_mat, intrp_ptr, result);
        if(errcode) {
            MUTIL_FREE_WARN(matdbl, &tmp_mat);
            return(errcode);
        }

        errcode = localfn_dctII(result, intrp_ptr);
        if(errcode) {
            MUTIL_FREE_WARN(matdbl, &tmp_mat);
            return(errcode);
        }

    } else {

        /* inverse transform: inverse-transform the columns of the
         * original matrix first, and then the rows of the
         * columnwise inverse-transformed matrix */

        errcode = matdbl_assign(sig, intrp_ptr, result);
        if(errcode) {
            MUTIL_FREE_WARN(matdbl, &tmp_mat);
            return(errcode);
        }

        errcode = localfn_idctII(result, intrp_ptr);
        if(errcode) {
            MUTIL_FREE_WARN(matdbl, &tmp_mat);
            return(errcode);
        }

        errcode = matdbl_transpose(result, intrp_ptr, &tmp_mat);
        if(errcode) {
            MUTIL_FREE_WARN(matdbl, &tmp_mat);
            return(errcode);
        }

        errcode = localfn_idctII(&tmp_mat, intrp_ptr);
        if(errcode) {
            MUTIL_FREE_WARN(matdbl, &tmp_mat);
            return(errcode);
        }

        errcode = matdbl_transpose(&tmp_mat, intrp_ptr, result);
        if(errcode) {
            MUTIL_FREE_WARN(matdbl, &tmp_mat);
            return(errcode);
        }
    } /* if-else for inverse-forward transform */

    /* clean up */

    MUTIL_FREE_WARN(matdbl, &tmp_mat);

    MUTIL_TRACE( "sigdbl_transform_discrete_cosine_II_2d() done" );
    return(MUTIL_ERR_OK);
}


/* Function documented in sig_tran.h */
/* Written by Luca Cazzanti */
mutil_errcode sigdbl_transform_discrete_fourier( const double_mat *sig,
    boolean inverse_flag, void *intrp_ptr, dcomplex_mat *result )
{

    mutil_errcode errcode;
    dcomplex_mat  sigcpx;

    MUTIL_TRACE( "Start sigdbl_transform_discrete_fourier()" );

    /* check arguments */

    errcode = matdbl_validate(sig);
    if(errcode) return errcode;

    errcode = matcpx_validate(result);
    if(errcode) return errcode;

    /* matrices must have same dimensions */

    if( !MATANY_EQUAL_DIM(sig, result) ) {
        MUTIL_ERROR( "Input and result matrices must have the same size" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* for now we use a dft for complex signals,
     * so we cast input to complex matrix
     */

    errcode = matcpx_malloc(&sigcpx, sig->nrow, sig->ncol);
    if(errcode) {
        return errcode;
    }

    errcode = matdbl_cast_to_cpx(sig, intrp_ptr, &sigcpx);
    if(errcode) {
        MUTIL_FREE_WARN(matcpx, &sigcpx);
        return errcode;
    }

    /* calculate using function for complex matrices.
     * in the future, we should write a function specific to
     * double matrices.
     */

    errcode = sigcpx_transform_discrete_fourier(&sigcpx, inverse_flag,
	intrp_ptr, result);
    MUTIL_FREE_WARN(matcpx, &sigcpx);
    if(errcode) {
        return errcode;
    }

    /* clean up */

    MUTIL_TRACE( "sigdbl_transform_discrete_fourier() done" );
    return MUTIL_ERR_OK;

}


/*
************************
Complex Matrix Functions
************************
*/

/* Function documented in sig_tran.h */
/* written by Luca Cazzanti */
mutil_errcode sigcpx_transform_discrete_fourier( const dcomplex_mat *sig,
    boolean inverse_flag, void * intrp_ptr, dcomplex_mat *result )
{

    mutil_errcode errcode;

    dcomplex_mat tmp_mat;

    sint32 isn; /* -1=forward transform; 1=inverse transform */
    sint32 row; /* indices */
    sint32 col;
    sint32 idx;

    MUTIL_TRACE( "Start sigcpx_transform_discrete_fourier()" );

    /* check arguments */

    errcode = matcpx_validate(sig);
    if(errcode) {
        return(errcode);
    }

    errcode = matcpx_validate(result);
    if(errcode) {
        return errcode;
    }

    /* dimensions must be equal */

    if( !MATANY_EQUAL_DIM(sig, result) ) {
        MUTIL_ERROR( "Input and result matrix must have same dimensions" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }


    /* accomplishes two things: copy the data from input to temporary
     * matrix for in-place localfn functions and puts data in
     * column-major format
     */

    errcode = matcpx_malloc(&tmp_mat, sig->nrow, sig->ncol);
    if(errcode) return(errcode);

    idx = 0;
    for( col = 0; col < tmp_mat.ncol; col++ ) {
        for( row=0; row < tmp_mat.nrow; row++ ) {
            tmp_mat.data[idx].re = sig->data[MATANY_INDEX( sig, row, col )].re;
            tmp_mat.data[idx].im = sig->data[MATANY_INDEX( sig, row, col )].im;
            idx++;
        }
    }

    /* set inverse/forward transform flag */

    if( inverse_flag == TRUE ) {
        isn =  1;
    }
    else {
        isn = -1;
    }

    errcode = localfn_cpxdft( tmp_mat.data, tmp_mat.nelem, tmp_mat.nrow,
	tmp_mat.nrow, isn, intrp_ptr );
    if(errcode) {
         MUTIL_FREE_WARN(matcpx, &tmp_mat);
        return errcode;
    }


     /* restore original data format */

    idx = 0;
    for( col = 0; col < result->ncol; col++ ) {
        for( row = 0; row < result->nrow; row++ ) {
            result->data[ MATANY_INDEX( result, row, col ) ].re =
                tmp_mat.data[ idx ].re;
            result->data[ MATANY_INDEX( result, row, col ) ].im =
                tmp_mat.data[ idx ].im;
            idx++;
        }
    }

    MUTIL_FREE_WARN(matcpx, &tmp_mat);

    MUTIL_TRACE( "sigcpx_transform_discrete_fourier() done" );

    return(MUTIL_ERR_OK);
}

/*
*****************
Static Functions
*****************
*/

/** Factorize length of signal for mixed-radix FFT.
 * Factorizes the length of a signal into shorter lengths, which can be
 * used to calculate the fast Fourier transform with the mixed-radix algorithm.
 * See \Ref{localfn_cpxdft} for references to the algorithm.
 *
 * This function is meant to be called only by the function localfn\_cpxdft,
 * and it exists only to improve its readability.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source sig\_tran.c
 * @library signal
 * @return void
 * @see localfn_cpxdft
 * @private
 */
static void localfn_factorize_siglen(sint32 n, sint32 *m, sint32 *i, sint32 *k,
    sint32 *k1, sint32 *kt, sint32 *nfac, sint32 *j, sint32 *jj, sint32 *maxf,
    double *num_ops )
{

    MUTIL_TRACE( "Start localfn_factorize_siglen()" );

    *m = 0;
    for ( *k = n; ( *k % 16 ) == 0; *k /= 16 ) {
	(*m)++;
	(*k1)++;
	nfac[*m] = 4;

	*num_ops += 50;
    }

    *j  = 2;
    *jj = 4;

    *i = *k / *jj;
    do{
	while ( ( *i = *k / *jj ) * *jj == *k ) {
	    (*m)++;
	    nfac[*m] = *j;
	    *k = *i;
	}

	*num_ops += 100;

	if ( *j == 2 )
	    *j = 3;
	else
	    *j += 2;

	*jj = *j * *j;
    } while ( *jj < *k );

  *kt = *m;

  for ( *j = 2; *j <= *k; ) {
      while ( ( *i = *k / *j ) * *j == *k ) {
	  (*m)++;
	  nfac[*m] = *j;
	  *k = *i;
      }

      *num_ops += 50;

      if ( *j == 2 )
	  *j =  3;
      else
	  *j += 2;
  }

  *i = *m - *kt;

  /* sets maxf to the largest factor + 1 */

  if ( nfac[*kt] > nfac[*m] )
      *maxf = 1 + ( ( nfac[*kt] > nfac[1] ) ? nfac[*kt] : nfac[1] );
  else
      *maxf = 1 + ( ( nfac[*m] > nfac[1] ) ? nfac[*m] : nfac[1] );

  /* get remaining parts of squared factors */

  for ( *j = *kt; *j >= 1; (*j)-- ) {
      (*m)++;
      nfac[*m] = nfac[*j];
  }

  *num_ops += 100;

  MUTIL_TRACE( "localfn_factorize_siglen() done" );

  return;
}


 /** Discrete Fourier transform and its inverse.
 * Calculates the DFT (and inverse DFT) of a signal for data of type complex.
 * It can calculate the DFT of a multichannel signal, by specifying a row
 * length less than the number of total points (see below).
 * It may be used to calculate the DFT of an N-dimensional signal, by
 * calling it repeatedly, each time with the appropriate values for the
 * arguments (see below).
 *
 * The inverse is calculated if the
 * paramenter isn is set to 1; the forward transform is computed if isn
 * is set to -1. The calculation is in-place, and it is assumed that the
 * data in array z is indexed in column-major format. The forward transform
 * does not scale the data; the inverse transform scales the data by the
 * inverse of parameter n.
 *
 * This function implements the mixed-radix algorithm by R. C. Singleton:
 * "An Algorithm for Computing the mixed radix fast Fourier transform",
 * IEEE Trans. Audio Electroacoust. AU-17(2), 93-103 (1969). Its current form
 * was taken from the S+SigPro signal processing module
 * for Splus. Modified by Jeff Silverman for Splus and
 * later by Luca Cazzanti for the mutils library.
 *
 * Example: univariate, single-channel signal:
 *
 * #err = localfn_cpxdft( &data, num_pts, num_pts, num_pts, isn, &intrp_ptr);#
 *
 * Example: multi-channel signal, with one channel per column:
 *
 * #err = localfn_cpxdft( &data, num_rows*num_cols, num_rows, num_rows, isn, &intrp_ptr );#
 *
 * Example: trivariate signal, with dimensions n1, n2, and n3:
 *
 * #err = localfn_cpxdft( &data, n1*n2*n3, n1, n1, isn, &intrp_ptr);#
 * #err = localfn_cpxdft( &data, n1*n2*n3, n2, n1*n2, isn, &intrp_ptr);#
 * #err = localfn_cpxdft( &data, n1*n2*n3, n3, n1*n2*n3, isn, &intrp_ptr);#
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source sig\_tran.c
 * @library signal
 * @usage #err = localfn_cpxdft(&data, tot_pts, num_dim, col_span, inverse_flag, &intrp_ptr);#
 * @return Standard mutils error/OK code.
 * @param z Pointer to complex array containing the original data.
 *          Will point to transformed data upon exit.
 * @param ntot Total number of data points.
 * @param n Dimension of current variable being transformed.
 * @param nspan Number of data points of current dimension.
 *              Note: the ratio nspan/n is the spacing between consecutive
 *              data points for the current dimension.
 * @param isn Set to -1 for forward transform, and 1 for inverse transform.
 * @param intrp_ptr Pointer to implement user interrupt.
 * @private
 */
static mutil_errcode localfn_cpxdft( dcomplex *z, sint32 ntot, sint32 n,
    sint32 nspan, sint32 isn, void *intrp_ptr )
{
    mutil_errcode errcode;

    /* angles and scaling coefficients (tweedle factors in FFT lingo)*/
    double rad;
    double s72;
    double c72;
    double s120;
    double radf;
    double sd;
    double cd;
    double c1;
    double c2;
    double c3;
    double s1;
    double s2;
    double s3;
    double akp;
    double akm;
    double ajp;
    double ajm;
    double bkp;
    double bkm;
    double bjp;
    double bjm;
    double aa;
    double bb;
    double ak;
    double bk;
    double aj;
    double bj;

    /* work arrays and arrays of coefficients for scaling */
    double *at;
    double *ck;
    double *bt;
    double *sk;
    sint32 *np;
    double *at_dummy = NULL; /* dummy variables to avoid lint warnings */
    double *ck_dummy = NULL;
    double *bt_dummy = NULL;
    double *sk_dummy = NULL;
    sint32 *np_dummy = NULL;
    sint32 nfac[20]; /* array of factors; size is arbitrary;
		        20 is historical value */
    sint32 vec_size; /* stores size of most vectors */

    /* counters, indices, helpers to keep track of array position, etc. */
    sint32 inc;
    sint32 nt;
    sint32 ks;
    sint32 kspan;
    sint32 nn;
    sint32 jc;
    sint32 i;
    sint32 jf;
    sint32 j;
    sint32 jj;
    sint32 m;
    sint32 k;
    sint32 kt;
    sint32 kk;
    sint32 k1;
    sint32 k2;
    sint32 k3;
    sint32 k4;
    sint32 kspnn;
    sint32 maxf;
    sint32 nvec;

    double normalizer; /* scaling factor for inverse fft */

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start localfn_cpxdft()" );



    /* do nothing in this case */

    if (n < 2)
        return(MUTIL_ERR_OK);

    /* set variables */

    c2 = 0;
    c3 = 0;
    s2 = 0;
    s3 = 0;
    k3 = 0;

    rad  = 6.28318530717958647692528e00;               /* 2 * pi */
    s120 = .8660254037844386467e00;                   /* sqrt(3)/2 */
    s72  = rad / 5.0;
    c72  = cos(s72);
    s72  = sin(s72);
    inc  = isn;

    if( isn < 0 ) {
	s72  = -s72;
	s120 = -s120;
	rad  = -rad;
	inc  = -inc;
    }

    nt      = inc * ntot;
    ks      = inc * nspan;
    kspan   = ks;
    nn      = nt - inc;
    jc      = ks / n;
    radf    = rad * jc * 0.5;
    jf      = 0;
    nfac[0] = 1;
    k1      = 0;

    /* factorize length into smaller factors */

    localfn_factorize_siglen( n, &m, &i, &k, &k1, &kt, nfac, &j, &jj, &maxf,
	&num_ops );

    /* after the above factorization:
     * kt = number of squared factors
     * m = number of factors
     * *nfac = vector of factors
     * maxf = largest factor
     * i = number of non-squared factors
     */

    /* setup work arrays */

    nvec = maxf;
    vec_size = (sint32) ( nvec * sizeof(double) );

    errcode = mutil_malloc ( vec_size, (void**) &at_dummy );
    if(errcode) {
	return(errcode);
    }
    at = at_dummy;

    errcode = mutil_malloc ( vec_size, (void**) &ck_dummy );
    if(errcode) {
	MUTIL_FREE_BUFFER_WARN( at_dummy, vec_size );
	return(errcode);
    }
    ck = ck_dummy;

    errcode = mutil_malloc ( vec_size, (void**) &bt_dummy );
    if(errcode) {
	MUTIL_FREE_BUFFER_WARN( at_dummy, vec_size );
	MUTIL_FREE_BUFFER_WARN( ck_dummy, vec_size );
	return(errcode);
    }
    bt = bt_dummy;

    errcode = mutil_malloc ( vec_size, (void**) &sk_dummy );
    if(errcode) {
	MUTIL_FREE_BUFFER_WARN( at_dummy, vec_size );
	MUTIL_FREE_BUFFER_WARN( ck_dummy, vec_size );
	MUTIL_FREE_BUFFER_WARN( bt_dummy, vec_size );
	return(errcode);
    }
    sk = sk_dummy;

    jj = 1;
     if ( i > 1 ) {
	 for ( j = 1; j <= i; j++ )
	     jj *= nfac[ kt + j ];
     }
     nvec = ( jj> m + 4 * k1 ) ? ( jj + 1 ):( m + 4 * k1 + 2 );

     errcode = mutil_malloc ( nvec * sizeof(sint32), (void**) &np_dummy );
     if(errcode) {
	 MUTIL_FREE_BUFFER_WARN( at_dummy, vec_size );
	 MUTIL_FREE_BUFFER_WARN( ck_dummy, vec_size );
	 MUTIL_FREE_BUFFER_WARN( bt_dummy, vec_size );
	 MUTIL_FREE_BUFFER_WARN( sk_dummy, vec_size );
         return(errcode);
     }
     np = np_dummy;


     /* error check
      * error check was moved just above, after each
      * allocation of memory
      */

     num_ops += 50;

     /* compute Fourier transform */

     i = 0;
  l1:
     if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }


     sd = radf / kspan;
     cd = 2.0 * sin(sd) * sin(sd);
     sd = sin( sd * 2.0 );
     kk = 0;
     i++;

     num_ops += 10;

     /* transform for factor of 2 */

     if( nfac[i] != 2 )
          goto l6;
     kspan /= 2;
     k1     = kspan + 1;

  l2: if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }


     k2        = kk + kspan;
     ak        = z[k2].re;
     bk        = z[k2].im;
     z[k2].re  = z[kk].re - ak;
     z[k2].im  = z[kk].im - bk;
     z[kk].re += ak;
     z[kk].im += bk;
     kk        = k2 + kspan;

     num_ops  += 15;

     if ( kk < nn)
          goto l2;
     kk -= nn;
     if (kk < jc)
          goto l2;
     if (kk >= kspan)
         goto l20;

  l3:  c1 = 1.0 - cd;
       s1 = sd;

l4:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k2        = kk + kspan;
     ak        = z[kk].re - z[k2].re;
     bk        = z[kk].im - z[k2].im;
     z[kk].re += z[k2].re;
     z[kk].im += z[k2].im;
     z[k2].re  = c1*ak - s1*bk;
     z[k2].im  = s1*ak + c1*bk;
     kk        = k2 + kspan;

     num_ops += 20;

     if (kk < nt-1)
          goto l4;
     k2 = kk - nt;
     c1 = -c1;
     kk = k1 - k2 - 1;

     num_ops += 10;

     if ( kk > k2)
          goto l4;
     ak  = c1 - (cd*c1 + sd*s1);
     s1 += (sd*c1-cd*s1);
     c1  = 2.0 - (ak*ak+s1*s1);
     s1  = c1*s1;
     c1  = c1*ak;
     kk += jc;

     num_ops += 25;

     if ( kk < k2)
          goto l4;
     k1 += 2*inc;
     kk  = (k1 - kspan)/2 + jc;
     if (kk < 2*jc)
          goto l3;
     goto l1;

     /* transform for factor of 3 (optional code) */

l5:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k1       = kk + kspan;
     k2       = k1 + kspan;
     ak       = z[kk].re;
     bk       = z[kk].im;
     aj       = z[k1].re + z[k2].re;
     bj       = z[k1].im + z[k2].im;
     z[kk].re = ak + aj;
     z[kk].im = bk + bj;
     ak      -= 0.5 * aj;
     bk      -= 0.5 * bj;
     aj       = (z[k1].re - z[k2].re) * s120;
     bj       = (z[k1].im - z[k2].im) * s120;
     z[k1].re = ak - bj;
     z[k1].im = bk + aj;
     z[k2].re = ak + bj;
     z[k2].im = bk - aj;
     kk =      k2 + kspan;

     num_ops += 50;

     if ( kk < nn-1)
          goto l5;
     kk -= nn;
     if (kk < kspan)
          goto l5;
     goto l16;

     /* transform for factor of 4 */

l6:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
      if (nfac[i] != 4)
          goto l11;
     kspnn = kspan;

#ifdef lint
        /* apparently the value just assigned to kspnn is never used! */
        kspan = kspnn;
#endif

     kspan /= 4;
l7:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     c1 = 1.0;
     s1 = 0.0;
l8:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k1       = kk + kspan;
     k2       = k1 + kspan;
     k3       = k2 + kspan;
     akp      = z[kk].re + z[k2].re;
     akm      = z[kk].re - z[k2].re;
     ajp      = z[k1].re + z[k3].re;
     ajm      = z[k1].re - z[k3].re;
     z[kk].re = akp + ajp;
     ajp      = akp - ajp;
     bkp      = z[kk].im + z[k2].im;
     bkm      = z[kk].im - z[k2].im;
     bjp      = z[k1].im + z[k3].im;
     bjm      = z[k1].im - z[k3].im;
     z[kk].im = bkp + bjp;
     bjp      = bkp - bjp;
     if (isn < 0) {
          akp  = akm + bjm;
          akm -= bjm;
          bkp  = bkm - ajm;
          bkm += ajm;
     }
     else {
          akp  = akm - bjm;
          akm += bjm;
          bkp  = bkm + ajm;
          bkm -= ajm;
     }
     if (s1 == 0.0) {             /*set up integer flag for zero */
          z[k1].re = akp;
          z[k1].im = bkp;
          z[k2].re = ajp;
          z[k2].im = bjp;
          z[k3].re = akm;
          z[k3].im = bkm;
     }
     else {
          z[k1].re = akp*c1 - bkp*s1;
          z[k1].im = akp*s1 + bkp*c1;
          z[k2].re = ajp*c2 - bjp*s2;
          z[k2].im = ajp*s2 + bjp*c2;
          z[k3].re = akm*c3 - bkm*s3;
          z[k3].im = akm*s3 + bkm*c3;
     }
     kk = k3 + kspan;

     num_ops += 100;

     if ( kk < nt)
          goto l8;
     c2  = c1 - (cd*c1 + sd*s1);
     s1 += (sd*c1 - cd*s1);
     c1  = 2.0 - (c2*c2 + s1*s1);
     s1  = c1*s1;
     c1  = c1*c2;
     c2  = c1*c1 - s1*s1;
     s2  = 2.0*c1*s1;
     c3  = c2*c1 - s2*s1;
     s3  = c2*s1 + s2*c1;
     kk -= (nt -jc);

     num_ops += 50;

     if (kk < kspan)
          goto l8;
     kk -= (kspan-inc);
     if (kk < jc)
          goto l7;
     if (kspan == jc)
          goto l20;
     goto l1;

          /* transform for factor of 5 (optional code) */

l9:   if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
    c2 = c72 * c72 - s72 * s72;
    s2 = 2.0 * c72 * s72;

l10:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k1 = kk + kspan;
     k2 = k1 + kspan;
     k3 = k2 + kspan;
     k4 = k3 + kspan;
     akp = z[k1].re + z[k4].re;
     akm = z[k1].re - z[k4].re;
     bkp = z[k1].im + z[k4].im;
     bkm = z[k1].im - z[k4].im;
     ajp = z[k2].re + z[k3].re;
     ajm = z[k2].re - z[k3].re;
     bjp = z[k2].im + z[k3].im;
     bjm = z[k2].im - z[k3].im;
     aa = z[kk].re;
     bb = z[kk].im;
     z[kk].re = aa + akp + ajp;
     z[kk].im = bb + bkp + bjp;
     ak = akp*c72 + ajp*c2 + aa;
     bk = bkp*c72 + bjp*c2 + bb;
     aj = akm*s72 + ajm*s2;
     bj = bkm*s72 + bjm*s2;
     z[k1].re = ak - bj;
     z[k4].re = ak + bj;
     z[k1].im = bk + aj;
     z[k4].im = bk - aj;
     ak = akp*c2 + ajp*c72 + aa;
     bk = bkp*c2 + bjp*c72 + bb;
     aj = akm*s2 - ajm*s72;
     bj = bkm*s2 - bjm*s72;
     z[k2].re = ak - bj;
     z[k3].re = ak + bj;
     z[k2].im = bk + aj;
     z[k3].im = bk - aj;
     kk = k4 + kspan;

     num_ops += 100;

     if ( kk < nn-1)
          goto l10;
     kk -= nn;
     if (kk < kspan)
          goto l10;
     goto l16;

          /* transform for odd factors */

l11:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k = nfac[i];
     kspnn = kspan;
     kspan /= k;
     if (k == 3)
          goto l5;
     if (k == 5)
          goto l9;
     if (k == jf)
          goto l12;
     jf = k;
     s1 = rad / k;
     c1 = cos(s1);
     s1 = sin(s1);

     num_ops += 20;

          /* error check */

     if (jf >= maxf) {
         MUTIL_ERROR( "Error in size of ck, sk, at, bt" );
         CLEAN_ALL;
         return MUTIL_ERR_ILLEGAL_SIZE;
     }

     ck[jf] = 1.0;
     sk[jf] = 0.0;
     j = 1;
     do {
          ck[j] = ck[k]*c1 + sk[k]*s1;
          sk[j] = ck[k]*s1 - sk[k]*c1;
          k--;
          ck[k] = ck[j];
          sk[k] = -sk[j];
          j++;

          num_ops += 1;

     } while (j < k);

l12:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k1 = kk;
     k2 = kk + kspnn;
     aa = z[kk].re;
     bb = z[kk].im;
     ak = aa;
     bk = bb;
     j = 1;
     k1 += kspan;

     num_ops += 10;

l13: if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k2 -= kspan;
     j++;
     at[j] = z[k1].re + z[k2].re;
     ak += at[j];
     bt[j] = z[k1].im + z[k2].im;
     bk += bt[j];
     j++;
     at[j] = z[k1].re - z[k2].re;
     bt[j] = z[k1].im - z[k2].im;
     k1 += kspan;

     num_ops += 25;

     if (k1 < k2)
          goto l13;
     z[kk].re = ak;
     z[kk].im = bk;
     k1 = kk;
     k2 = kk + kspnn;
     j = 1;

     num_ops += 10;

l14: if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
      k1 += kspan;
     k2 -= kspan;
     jj = j;
     ak = aa;
     bk = bb;
     aj = 0.0;
     bj = 0.0;
     k = 1;

     num_ops += 10;

l15:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k++;
     ak += at[k]*ck[jj];
     bk += bt[k]*ck[jj];
     k++;
     aj += at[k]*sk[jj];
     bj += bt[k]*sk[jj];
     jj += j;
     if (jj > jf)
         jj -= jf;

     num_ops += 20;

     if (k < jf)
          goto l15;
     k = jf - j;
     z[k1].re = ak - bj;
     z[k1].im = bk + aj;
     z[k2].re = ak + bj;
     z[k2].im = bk - aj;
     j++;

     num_ops += 10;

     if (j < k)
          goto l14;
     kk += kspnn;
     if (kk < nn)
          goto l12;
     kk -= nn;
     if (kk < kspan)
          goto l12;

     /* multiplying by rotation factor (except for factors of 2 and 4) */

l16:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
    if (i == m)
          goto l20;
     kk = jc;
l17:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     c2 = 1.0 - cd;
     s1 = sd;
l18:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     c1 = c2;
     s2 = s1;
     kk += kspan;
l19:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     ak = z[kk].re;
     z[kk].re = c2*ak - s2*z[kk].im;
     z[kk].im = s2*ak + c2*z[kk].im;
     kk += kspnn;

     num_ops += 20;

     if (kk < nt)
          goto l19;

     ak = s1*s2;
     s2 = s1*c2 + c1*s2;
     c2 = c1 * c2 - ak;
     kk -= (nt - kspan);

     num_ops += 10;

     if (kk < kspnn)
          goto l19;
     c2 = c1 - (cd*c1 + sd*s1);
     s1 += (sd*c1 - cd*s1);
     c1 = 2.0 - (c2*c2 + s1*s1);
     s1 = c1*s1;
     c2 = c1*c2;
     kk -= (kspnn - jc);

     num_ops += 20;

     if (kk < kspan)
         goto l18;

     kk -= (kspan - jc - inc);
     if (kk < 2*jc)
          goto l17;
     goto l1;

     /* permute results to normal order---done in two stages */

     /*permutation for square factors of n */

l20:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     np[1] = ks;
     if (kt == 0)
          goto l28;
     k = 2*kt + 1;
     if ( m < k)
          k--;
     j = 1;
     np[k+1] = jc;
     do {
          np[j+1] = np[j]/nfac[j];
          np[k] = np[k+1]*nfac[j];
          j++;
          k--;

          num_ops += 10;

     } while (j < k);
     k3 = np[k+1];
     kspan = np[2];
     kk = jc;
     k2 = kspan;
     j = 1;
     if (n != ntot)
          goto l24;

     /* permutation for single-variate transform (optional code) */

l21:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     do {
          ak = z[kk].re;
          z[kk].re = z[k2].re;
          z[k2].re = ak;
          bk = z[kk].im;
          z[kk].im = z[k2].im;
          z[k2].im = bk;
          kk += inc;
          k2 += kspan;

          num_ops += 10;

     } while (k2 < ks-1);
l22:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k2 -= np[j];
     j++;
     k2 += np[j+1];

     num_ops += 10;

     if ( k2 >= np[j])
          goto l22;
     j=1;
l23:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     if (kk < k2)
          goto l21;
     kk += inc;
     k2 += kspan;
     if (k2 < ks-1)
          goto l23;
     if (kk < ks-1)
          goto l22;
     jc = k3;
     goto l28;

     /* permutation for multivariate transform */

l24: k = kk + jc;
l25:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     ak = z[kk].re;
     z[kk].re = z[k2].re;
     z[k2].re = ak;
     bk = z[kk].im;
     z[kk].im = z[k2].im;
     z[k2].im = bk;
     kk += inc;
     k2 += inc;

     num_ops += 10;

     if ( kk < k)
          goto l25;
     kk += ks - jc;
     k2 += ks - jc;
     if (kk < nt-1)
          goto l24;
     k2 -= nt - kspan;
     kk -= nt - jc;
     if (k2 < ks-1)
          goto l24;
l26:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k2 -= np[j];
     j++;
     k2 += np[j+1];
     if (k2 >= np[j])
          goto l26;
     j = 1;
l27:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     if (kk < k2)
          goto l24;
     kk += jc;
     k2 += kspan;
     if (k2 < ks-1)
          goto l27;
     if (kk < ks-1)
          goto l26;
     jc = k3;
l28:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     if (2*kt+1 >= m)
          goto l42;
     kspnn = np[kt+1];

     /* permutation for square free factors of n */

     j = m - kt;
     nfac[j+1] = 1;
l29:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     nfac[j] = nfac[j] * nfac[j+1];
     j--;
     if (j != kt)
          goto l29;
     kt++;
     nn = nfac[kt] - 1;

     /* error check */

     if ( nn >= nvec) {
         MUTIL_ERROR( "Error in size of np" );
         CLEAN_ALL;
         return MUTIL_ERR_ILLEGAL_SIZE;
     }

     jj = j = 0;
     goto l32;
l30:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     jj -= k2;
     k2 = kk;
     k++;
     kk = nfac[k];
l31:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     jj += kk;
     if (jj >= k2)
          goto l30;
     np[j] = jj;
l32: k2 = nfac[kt];
     k = kt + 1;
     kk = nfac[k];
     j++;
     if (j <= nn)
          goto l31;

     /* determine the permutation cycles of length greater than one */

     j = 0;
     goto l34;
l33:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k = kk;
     kk = np[k];
     np[k] = -kk;

     num_ops += 10;

     if (kk != j)
          goto l33;
     k3 = kk;
l34:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     j++;
     kk = np[j];
     if (kk < 0)
          goto l34;
     if (kk != j)
          goto l33;
     np[j] = -j;

     num_ops += 10;

     if (j != nn)
          goto l34;
     maxf = inc * maxf;

     /* reordering a and b according to the permutation cycles */

     goto l41;
l35:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     j--;
     if (np[j] < 0)
          goto l35;
     jj = jc;
l36:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     kspan = jj;
     if ( jj > maxf)
          kspan = maxf;
     jj -= kspan;
     k = np[j];
     kk = jc * k + i + jj - 1;
     k1 = kk + kspan;
     k2 = -1;

     num_ops += 10;

l37:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k2++;
     at[k2] = z[k1].re;
     bt[k2] = z[k1].im;
     k1 -= inc;

     num_ops += 10;

     if (k1 != kk)
          goto l37;
l38:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k1 = kk + kspan;
     k2 = k1 - jc * (k + np[k]);
     k = -np[k];

     num_ops += 10;

l39:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     z[k1].re = z[k2].re;
     z[k1].im = z[k2].im;
     k1 -= inc;
     k2 -= inc;
     if (k1 != kk)
          goto l39;
     kk = k2;
     if (k != j)
          goto l38;
     k1 = kk + kspan;
     k2 = -1;
l40:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     k2++;
     z[k1].re = at[k2];
     z[k1].im = bt[k2];
     k1 -= inc;

     if (k1 != kk)
          goto l40;
     if (jj != 0)
          goto l36;
     if (j != 1)
          goto l35;
l41:  if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
         MUTIL_ERROR( "User interrupt" );
         CLEAN_ALL;
         return(errcode);
     }
     j = k3 + 1;
     nt -= kspnn;
     i = nt - inc + 1;
     if (nt >= 0)
          goto l35;

  l42:
     CLEAN_ALL;



/* normalize the inverse transform */
     if ( isn > 0 ) {
         normalizer = 1/(double) n;
         for ( j=0; j < ntot; j++ ) {
             z[j].re *= normalizer;
             z[j].im *= normalizer;
         }
         num_ops += 3 * ntot;
         if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
	     MUTIL_ERROR( "User interrupt" );
	     return(errcode);
	 }

     }

     MUTIL_TRACE( "localfn_cpxdft() done" );

     return MUTIL_ERR_OK;
}

/** Discrete cosine transform of type II.
 * Computes the discrete cosine transform of type II for
 * a signal of any length. For 1D signals, the data should be
 * arranged in a column vector. If the input matrix is 2D, it
 * computes the transform of each column. This is useful for
 * multi-channel signals, with one channel per column.
 *
 * Uses a discrete Fourier transform algorithm to compute the transform,
 * exploiting the close relationship between the DFT and the DCT-II.
 * The algorithm calculates
 * the DFT of a sequence which is twice as long as the original one,
 * and zero padded at the end. Then to go from the DFT to the DCT, the
 * following formula is used for the first N elements of the DFT, where
 * N is the original signal length:
 *
 * $G(n) = \sqrt{\frac{1}{N}} (\cos(\alpha) \Re(G_{dft}(n)) +
 *          \sin(\alpha) \Im(G_{dft}(n)))$
 *
 * where $\alpha$ is a sequence of angles with values
 * $\frac{n \pi}{2 N}$
 *
 * Based on original S code by Andrew Bruce, appearing in the
 * S+Wavelets module. Adapted by Luca Cazzanti to conform to
 * the mutils standards.
 *
 * @usage #err = localfn_dctII(&signal, &intrp_ptr);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source sig\_tran.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param sig Pointer to double matrix containing data to be transformed.
 *            Will point to transformed data upon return.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @see siguniv_transform_discrete_cosine_II
 * @see localfn_idctII
 * @see localfn_dctII_pow2
 * @see localfn_idctII_pow2
 * @private
 */
static mutil_errcode localfn_dctII(double_mat *sig, void *intrp_ptr)
{
    mutil_errcode errcode;

    dcomplex_mat tmp_mat;
    dcomplex_mat tmp_mat2;
    double_mat angles_mat;
    double scale_factor;
    sint32 row_idx;
    sint32 col_idx;
    sint32 tmp;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start localfn_dctII()" );

    /* no argument checking because called only by
     * functions in this file */

    /* allocate temporary matrix needed for row-col-major change */

    errcode = matcpx_malloc(&tmp_mat2, 2 * sig->nrow, sig->ncol);
    if(errcode) {
        MUTIL_ERROR( "Could not allocate temporary matrix" );
        return(errcode);
    }

    /* must zero-pad the rows */

    for(row_idx = sig->nrow; row_idx < tmp_mat2.nrow; row_idx++) {
        for(col_idx = 0; col_idx < tmp_mat2.ncol; col_idx++) {
            tmp_mat2.data[ MATANY_INDEX(&tmp_mat2, row_idx, col_idx) ].re = 0;
            tmp_mat2.data[ MATANY_INDEX(&tmp_mat2, row_idx, col_idx) ].im = 0;
        }
        num_ops += 3 * tmp_mat2.ncol;
        if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
            MUTIL_ERROR( "User interrupt" );
            MUTIL_FREE_WARN(matcpx, &tmp_mat2);
            return MUTIL_ERR_INTERRUPT;
        }
    }

    /* assign real data to complex matrix (imaginary part is zero) */

    for(col_idx = 0; col_idx < sig->ncol; col_idx++) {
        for(row_idx = 0; row_idx < sig->nrow; row_idx++) {
            tmp_mat2.data[ MATANY_INDEX( &tmp_mat2, row_idx, col_idx ) ].re=
                sig->data[ MATANY_INDEX( sig, row_idx, col_idx ) ];

            tmp_mat2.data[ MATANY_INDEX( &tmp_mat2, row_idx, col_idx ) ].im= 0;
        }
        num_ops += 2 * sig->nrow;
        if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
            MUTIL_ERROR( "User interrupt" );
            MUTIL_FREE_WARN(matcpx, &tmp_mat2);
            return MUTIL_ERR_INTERRUPT;
        }
    }

    /* accomplishes two things: copies data to tmp_mat and
     * puts data in column-major format for dft routine
     */

    errcode = matcpx_malloc(&tmp_mat, 2 * sig->nrow, sig->ncol);
    if(errcode) {
        MUTIL_FREE_WARN(matcpx, &tmp_mat2);
        return(errcode);
    }

    tmp = 0;
    for(col_idx = 0; col_idx < tmp_mat.ncol; col_idx++) {
        for(row_idx = 0; row_idx < tmp_mat.nrow; row_idx++) {
            tmp_mat.data[ tmp ].re =
                tmp_mat2.data[ MATANY_INDEX(&tmp_mat2, row_idx, col_idx) ].re;
            tmp_mat.data[ tmp ].im =
                tmp_mat2.data[ MATANY_INDEX(&tmp_mat2, row_idx, col_idx) ].im;
            tmp++;
        }
        num_ops += 3 * tmp_mat.ncol;
        if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
            MUTIL_ERROR( "User interrupt" );
            MUTIL_FREE_WARN(matcpx, &tmp_mat2);
            return MUTIL_ERR_INTERRUPT;
        }
    }

    /* calculate DFT */

    errcode = localfn_cpxdft( tmp_mat.data, tmp_mat.nelem,
                              tmp_mat.nrow, tmp_mat.nrow,
                              -1, intrp_ptr );
    if(errcode) {
        MUTIL_ERROR( "Error computing Fourier transform" );
        MUTIL_FREE_WARN(matcpx, &tmp_mat);
        MUTIL_FREE_WARN(matcpx, &tmp_mat2);
        return(errcode);
    }

    /* restore row-major format and free some memory
     * that is not needed any more
     */

    tmp = 0;
    for(col_idx = 0; col_idx < tmp_mat.ncol; col_idx++) {
        for(row_idx = 0; row_idx < tmp_mat.nrow; row_idx++) {
            tmp_mat2.data[ MATANY_INDEX(&tmp_mat2, row_idx, col_idx) ].re=
                tmp_mat.data[ tmp ].re;
            tmp_mat2.data[ MATANY_INDEX( &tmp_mat2, row_idx, col_idx) ].im=
                tmp_mat.data[ tmp ].im;
            tmp++;
        }
        num_ops += 3 * tmp_mat2.ncol;
        if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
            MUTIL_ERROR( "User interrupt" );
            MUTIL_FREE_WARN(matcpx, &tmp_mat2);
            return MUTIL_ERR_INTERRUPT;
        }
    }

    MUTIL_FREE_WARN(matcpx, &tmp_mat);

    /* allocate matrix of angles for conversion from DFT to DCT */

    errcode = matdbl_malloc(&angles_mat, sig->nrow, 1);
    if(errcode) {
        MUTIL_ERROR( "Could not allocate matrix for conversion angles" );
        MUTIL_FREE_WARN( matcpx, &tmp_mat2 );
        return(errcode);
    }

    for(row_idx = 0; row_idx < sig->nrow; row_idx++) {
        angles_mat.data[row_idx] = MUTIL_PI *
            (double) row_idx / (double) ( 2 * sig->nrow );
    }
    num_ops += 3 * sig->nrow;

    /* go from DFT to DCT and put result in original matrix */

    /* the first elements of each column get a special scaling factor */

    scale_factor = 1.0 / sqrt( (double) sig->nrow);
    for(col_idx = 0; col_idx < sig->ncol; col_idx++) {
        sig->data[ MATANY_INDEX( sig, 0, col_idx ) ] = scale_factor *
                ( cos( angles_mat.data[ 0 ] ) *
                  tmp_mat2.data[ MATANY_INDEX( &tmp_mat2, 0, col_idx ) ].re
                  +
                  sin( angles_mat.data[ 0 ] ) *
                  tmp_mat2.data[ MATANY_INDEX( &tmp_mat2, 0, col_idx ) ].im
                    );
    }
    num_ops += 5 * sig->nrow;

    /* the rest of the elements of each column get a
     * different caling factor */

    scale_factor *= SQRT_2;
    for(col_idx = 0; col_idx < sig->ncol; col_idx++) {
        for(row_idx = 1; row_idx < sig->nrow; row_idx++) {
            sig->data[ MATANY_INDEX( sig, row_idx, col_idx ) ] = scale_factor *
                ( cos( angles_mat.data[ row_idx ] ) *
                  tmp_mat2.data[ MATANY_INDEX( &tmp_mat2,
                                               row_idx, col_idx ) ].re
                  +
                  sin( angles_mat.data[ row_idx ] ) *
                  tmp_mat2.data[MATANY_INDEX( &tmp_mat2,
                                              row_idx, col_idx ) ].im
                 );
        }
        num_ops += 5 * sig->nrow;
        if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
            MUTIL_ERROR( "User interrupt" );
            MUTIL_FREE_WARN(matcpx, &tmp_mat2);
            MUTIL_FREE_WARN(matdbl, &angles_mat);
            return MUTIL_ERR_INTERRUPT;
        }
    }


    /* clean up */

    MUTIL_FREE_WARN(matdbl, &angles_mat);
    MUTIL_FREE_WARN(matcpx, &tmp_mat2);

    MUTIL_TRACE( "localfn_dctII() done" );
    return(MUTIL_ERR_OK);
}

/** Inverse discrete cosine transform of type II.
 * Computes the inverse discrete cosine transform of type II for
 * a signal of any length. For 1D signals, the data should be
 * arranged in a column vector. If the input matrix is 2D, it
 * computes the transform of each column. This is useful for
 * multi-channel signals, with one channel per column.
 *
 * Uses a discrete Fourier transform algorithm to compute the transform,
 * exploiting the close relationship between the DFT and the DCT-II.
 * The algorithm calculates
 * the DFT of a sequence which is four times as long as the original one,
 * and zero padded at the end. Then, to go from the DFT to the original
 * signal, the following formula is used for the first N elements of
 * the DFT, where N is the original signal length:
 *
 * $g(n) = \sqrt{\frac{1}{N}} \Re(G_{dft}(2n+1))$
 *
 * Based on original S code by Andrew Bruce appearing in the
 * S+Wavelets module. Adapted by Luca Cazzanti to conform to
 * the mutils standards.
 *
 * @usage #err = localfn_idctII(&signal, &intrp_ptr);#
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source sig\_tran.c
 * @library signal
 * @return Standard mutils error/OK code.
 * @param sig Pointer to double matrix containing data to be transformed.
 *     Will point to transformed data upon return.
 * @param intrp_ptr Pointer for implementing user interrupt.
 * @see siguniv_transform_discrete_cosine_II_inverse
 * @see localfn_dctII
 * @see localfn_idctII_pow2
 * @see localfn_dctII_pow2
 * @private
 */
static mutil_errcode localfn_idctII(double_mat *sig, void *intrp_ptr)
{
    mutil_errcode errcode;

    dcomplex_mat tmp_mat;
    dcomplex_mat tmp_mat2;
    double scale_factor;
    sint32 row_idx;
    sint32 col_idx;
    sint32 tmp;

    double num_ops = 0;

    MUTIL_INTERRUPT_INIT(intrp_ptr);

    MUTIL_TRACE( "Start localfn_idctII()" );

    /* allocate temporary complex matrix */

    errcode = matcpx_malloc(&tmp_mat2, 4*(sig->nrow), sig->ncol);
    if(errcode) {
        MUTIL_ERROR( "Could not allocate temporary matrix" );
        return(errcode);
    }

    /* must zero-pad temporary matrix */

    for(col_idx = 0; col_idx < tmp_mat2.ncol; col_idx++) {
        for(row_idx = sig->nrow; row_idx < tmp_mat2.nrow; row_idx++) {
            tmp_mat2.data[ MATANY_INDEX(&tmp_mat2, row_idx, col_idx) ].re = 0;
            tmp_mat2.data[ MATANY_INDEX(&tmp_mat2, row_idx, col_idx) ].im = 0;
        }
        num_ops += 3 * tmp_mat2.nelem;
        if( MUTIL_INTERRUPT(num_ops, intrp_ptr) ) {
            MUTIL_ERROR( "User interrupt" );
            MUTIL_FREE_WARN(matcpx, &tmp_mat2);
            return MUTIL_ERR_INTERRUPT;
        }
    }

    /* copy real data into temporary matrix (imaginary part set to zero) */

    for( col_idx = 0; col_idx < sig->ncol; col_idx++ ) {
        for( row_idx = 0; row_idx < sig->nrow; row_idx++ ) {
            tmp_mat2.data[ MATANY_INDEX( &tmp_mat2, row_idx, col_idx ) ].re =
                sig->data[ MATANY_INDEX( sig, row_idx, col_idx ) ];

            tmp_mat2.data[ MATANY_INDEX( &tmp_mat2, row_idx, col_idx ) ].im =0;

        }
        num_ops += 2 * sig->nrow;
        if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
            MUTIL_ERROR( "User interrupt" );
            MUTIL_FREE_WARN(matcpx, &tmp_mat2);
            return MUTIL_ERR_INTERRUPT;
        }
    }

    /* Scaling required by DCT-through-DFT approach. */

    for(col_idx = 0; col_idx < tmp_mat2.ncol; col_idx++) {
        tmp_mat2.data[ MATANY_INDEX( &tmp_mat2, 0, col_idx ) ].re /= 2.0;
    }
    for( col_idx = 0; col_idx < sig->ncol; col_idx++ ) {
        for( row_idx = 1; row_idx < sig->nrow; row_idx++ ) {
            tmp_mat2.data[ MATANY_INDEX( &tmp_mat2, row_idx, col_idx ) ].re
                /= SQRT_2;
        }
        num_ops += 2 * sig->nrow;
        if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
            MUTIL_ERROR( "User interrupt" );
            MUTIL_FREE_WARN(matcpx, &tmp_mat2);
            return MUTIL_ERR_INTERRUPT;
        }
    }

    /* copy data in temporary matrix which is in column-major format */

    errcode = matcpx_malloc(&tmp_mat, tmp_mat2.nrow, tmp_mat2.ncol);
    if(errcode) {
        MUTIL_FREE_WARN(matcpx, &tmp_mat2);
        return errcode;
    }

    /* the dft function assumes column major order */

    tmp = 0;
    for(col_idx = 0; col_idx < tmp_mat.ncol; col_idx++) {
        for(row_idx = 0; row_idx < tmp_mat.nrow; row_idx++) {
            tmp_mat.data[ tmp ].re =
                tmp_mat2.data[ MATANY_INDEX(&tmp_mat2, row_idx, col_idx) ].re;
            tmp_mat.data[ tmp ].im =
                tmp_mat2.data[ MATANY_INDEX(&tmp_mat2, row_idx, col_idx) ].im;
            tmp++;
        }
        num_ops += 3 * tmp_mat2.ncol;
        if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
            MUTIL_ERROR( "User interrupt" );
            MUTIL_FREE_WARN(matcpx, &tmp_mat2);
            return MUTIL_ERR_INTERRUPT;
        }
    }

    /* calculate DFT */

    errcode = localfn_cpxdft( tmp_mat.data, tmp_mat.nelem,
                              tmp_mat.nrow, tmp_mat.nrow,
                              -1, intrp_ptr );
    if(errcode) {
        MUTIL_ERROR( "Error computing Fourier transform" );
        MUTIL_FREE_WARN( matcpx, &tmp_mat );
        return(errcode);
    }

    /* restore row-major format */

    tmp = 0;
    for(col_idx = 0; col_idx < tmp_mat.ncol; col_idx++) {
        for(row_idx = 0; row_idx < tmp_mat.nrow; row_idx++) {
            tmp_mat2.data[ MATANY_INDEX(&tmp_mat2, row_idx, col_idx) ].re=
                tmp_mat.data[ tmp ].re;
             tmp_mat2.data[ MATANY_INDEX(&tmp_mat2, row_idx, col_idx) ].im=
                tmp_mat.data[ tmp ].im;
            tmp++;
        }
        num_ops += 3 * tmp_mat2.ncol;
        if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
            MUTIL_ERROR( "User interrupt" );
            MUTIL_FREE_WARN(matcpx, &tmp_mat2);
            return MUTIL_ERR_INTERRUPT;
        }
    }

    MUTIL_FREE_WARN(matcpx, &tmp_mat);

    /* go from DFT to DCT */

    scale_factor = 2.0 / sqrt( (double) sig->nrow );
    for( col_idx = 0; col_idx < sig->ncol; col_idx++ ) {
        for( row_idx = 0; row_idx< sig->nrow; row_idx++) {
            sig->data[ MATANY_INDEX( sig, row_idx, col_idx ) ] =
                scale_factor *
                tmp_mat2.data[ MATANY_INDEX( &tmp_mat2,
                                             (2 * row_idx + 1), col_idx )].re;
        }
        num_ops += 5 * sig->nrow;
        if( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
            MUTIL_ERROR( "User interrupt" );
            MUTIL_FREE_WARN(matcpx, &tmp_mat2);
            return MUTIL_ERR_INTERRUPT;
        }
    }

    /* clean up */

    MUTIL_FREE_WARN(matcpx, &tmp_mat2);

    MUTIL_TRACE( "localfn_idctII() done" );
    return(MUTIL_ERR_OK);
}

/*
**************************
Universal Matrix Functions
**************************
*/

/* Function documented in sig_tran.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_transform_discrete_cosine_II(const univ_mat *sig,
    void *intrp_ptr, univ_mat *result)
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_transform_discrete_cosine_II()" );

    /* check arguments */

    if( !sig || !result ) {
        MUTIL_ERROR( "NULL pointer for input or operand" );
        return MUTIL_ERR_NULL_POINTER;
    }

    if( !MATUNIV_CHECK_TYPE( sig, result ) ) {
        MUTIL_ERROR( "Type mismatch between input and output matrices" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }

    switch(sig->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_transform_discrete_cosine_II( &(sig->mat.dblmat), intrp_ptr,
                &(result->mat.dblmat) );
            if(errcode) {
                return errcode;
            }
            break;

        default:

            /* only some matrix types supported for now */

            MUTIL_ERROR( "This matrix type is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_transform_discrete_cosine_II() done" );
    return MUTIL_ERR_OK;

}

/* Function documented in sig_tran.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_transform_discrete_cosine_II_inverse(
    const univ_mat *sig, void *intrp_ptr, univ_mat *result )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_transform_discrete_cosine_II_inverse()" );

    /* check arguments */

    if( !sig || !result ) {
        MUTIL_ERROR( "NULL pointer for input or operand" );
        return MUTIL_ERR_NULL_POINTER;
    }

    if( !MATUNIV_CHECK_TYPE( sig, result ) ) {
        MUTIL_ERROR( "Type mismatch between input and output matrices" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }

    switch(sig->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_transform_discrete_cosine_II_inverse( &(sig->mat.dblmat), intrp_ptr,
                &(result->mat.dblmat) );
            if(errcode) {
                return errcode;
            }
            break;

        default:

            /* only some matrix types supported for now */

            MUTIL_ERROR( "This matrix type is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_transform_discrete_cosine_II_inverse() done" );
    return MUTIL_ERR_OK;

}

/* Function documented in sig_tran.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_transform_discrete_cosine_II_2d( const univ_mat *sig,
    boolean inverse_flag, void *intrp_ptr, univ_mat *result )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_transform_discrete_cosine_II_2d()" );

    /* check arguments */

    if( !sig || !result ) {
        MUTIL_ERROR( "NULL pointer for input or operand" );
        return MUTIL_ERR_NULL_POINTER;
    }

    if( !MATUNIV_CHECK_TYPE( sig, result ) ) {
        MUTIL_ERROR( "Type mismatch between input and output matrices" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }

    switch(sig->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_transform_discrete_cosine_II_2d(
		&(sig->mat.dblmat), inverse_flag, intrp_ptr,
		&(result->mat.dblmat) );
            if(errcode) {
                return errcode;
            }
            break;

        default:

            /* only some matrix types supported for now */

            MUTIL_ERROR( "This matrix type is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_transform_discrete_cosine_II_2d() done" );
    return MUTIL_ERR_OK;

}

/* Function documented in sig_tran.h */
/* Written by Luca Cazzanti */
mutil_errcode siguniv_transform_discrete_fourier( const univ_mat *sig,
    boolean inverse_flag, void *intrp_ptr, univ_mat *result )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start siguniv_transform_discrete_fourier()" );

    /* check arguments */

    if( !sig || !result ) {
        MUTIL_ERROR( "NULL pointer for input or operand" );
        return MUTIL_ERR_NULL_POINTER;
    }

    /* the result must be a double-precision complex matrix */

    if( result->type != MUTIL_DCOMPLEX ) {
        MUTIL_ERROR( "Output data must be complex" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }

    switch(sig->type) {
        case MUTIL_DOUBLE:
            errcode = sigdbl_transform_discrete_fourier( &(sig->mat.dblmat),
		inverse_flag, intrp_ptr, &(result->mat.cpxmat) );
            if(errcode) {
                return errcode;
            }
            break;

        case MUTIL_DCOMPLEX:
            errcode = sigcpx_transform_discrete_fourier( &(sig->mat.cpxmat),
		inverse_flag, intrp_ptr, &(result->mat.cpxmat) );
            if(errcode) {
                return errcode;
            }
            break;

        default:

            /* only some matrix types supported for now */

            MUTIL_ERROR( "This matrix type is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "siguniv_transform_discrete_fourier() done" );
    return MUTIL_ERR_OK;
}

