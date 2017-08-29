
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/fra_fdp.c $: $Revision: #1 $, $Date: 2008/03/21 $";

#include "fra_fdp.h"

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_limit.h"
#include "ut_math.h"
#include "ut_mem.h"
#include "ut_rand.h"

#include "mat_any.h"
#include "mat_arit.h"
#include "mat_assn.h"
#include "mat_set.h"
#include "mat_sort.h"
#include "mat_summ.h"
#include "mat_type.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "mth_mac.h"

/***********************************************/
/* STATIC SIMULATION FUNCTION DECLARATIONS     */
/***********************************************/

static mutil_errcode localfn_cumsum_lower_triangle(
  double_mat *matrix);

static mutil_errcode localfn_upper_triangular_back_substitution(
  double_mat *left,
  sint32      index,
  double_mat *result);

static mutil_errcode localfn_fdp_Cholesky_matrix(
  double      delta,
  double_mat *result );

static mutil_errcode localfn_mse_prediction(
  double      delta,
  double      innovation_variance,
  double_mat *L,
  double_mat *result );

static mutil_errcode localfn_fdp_simulate_input_check(
  const univ_mat *delta,
  const univ_mat *innovation_variance);

static mutil_errcode localfn_noise_gaussian(
  double_mat *result );

/***********************************************/
/* STATIC SIMULATION MACRO DEFINITIONS         */
/***********************************************/

#define INDEX( row, col, N ) ( (row) * (N) ) + (col)

#define LOCALDEF_CHECK_NULL_POINTER_FDP( DATA_PTR, DATA_TYPE,  \
                                     TYPE_PREFIX )             \
   err = TYPE_PREFIX ## _validate( DATA_PTR );                 \
   if ( err ) return err;                                      \
   if ( DATA_PTR == ( DATA_TYPE * ) NULL ){                    \
     MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
     return MUTIL_ERR_NULL_POINTER;                            \
   }

/***********************************************/
/* SIMULATION FUNCTION DEFINITIONS             */
/***********************************************/

/* Function to simulate an FD process */
/* whose model parameters vary as a   */
/* function of time.                  */
/* Function documented in fra_fdp.h   */
/* Written by William Constantine     */

mutil_errcode frauniv_fdp_simulate(
  const univ_mat *delta,
  const univ_mat *innovation_variance,
  void           *intrp_ptr,
  univ_mat       *result)
{
  double        cumsum;
  double_mat    noise;
  memlist       list;
  mutil_errcode err;
  sint32        N;
  sint32        k;
  sint32        row;
  sint32        t;
  univ_mat      weights;

  /* initialize interrupt pointer */

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_fdp_simulate() ::" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize the memory list */

  MEMLIST_INIT( list );

  /* validate the argument list inputs */

  err = localfn_fdp_simulate_input_check(
    delta, innovation_variance );
  if ( err ) return err;

  /* obtain the realization length */

  N = MATUNIV_NELEM( delta );

  /* allocate space for the output */

  err = matuniv_malloc_register( result, N, 1, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* allocate space for common matrices */

  err = matuniv_malloc_register( &weights, N, N,
                             MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &noise, N, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate the FD weights matrix */

  err = frauniv_fdp_simulate_weights(
    delta,
    innovation_variance,
    intrp_ptr,
    &weights);
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* generate a vector of independent and
     identically distributed random variables
     with zero mean and unit variance */

  err = localfn_noise_gaussian( &noise );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* form a cumulative summation over the
     product of the weights with the
     Gaussian distributed pseudo-random
     noise */

  for ( t = 0; t < N; t++ ){

    row = t * N;
    cumsum = 0.0;

    for ( k = 0; k <= t; k++ ){
      cumsum += weights.mat.dblmat.data[ row + k ] * noise.data[ k ];
    }

    result->mat.dblmat.data[ t ] = cumsum;

    /* Check for interrupts */

    if ( MUTIL_INTERRUPT( 3.0 * N, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );

      return MUTIL_ERR_INTERRUPT;
    }
  }

  /* free the node corresponding to the
     output in the memory list but not the
     memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free the remaining registered memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( ":: Finished frauniv_fdp_simulate()" );

  return MUTIL_ERR_OK;
}

/* Function to generate the FD weights */
/* for a time varying FD simulation.   */
/* Function documented in fra_fdp.h    */
/* Written by William Constantine      */

mutil_errcode frauniv_fdp_simulate_weights(
  const univ_mat *delta,
  const univ_mat *innovation_variance,
  void           *intrp_ptr,
  univ_mat       *result )
{
  double        delta_t;
  double_mat    L;
  double_mat    mse;
  double_mat    weight_row;
  double_mat    weights;
  memlist       list;
  mutil_errcode err;
  sint32        N;
  sint32        k;
  sint32        kstart;
  sint32        nsums;
  sint32        t;

  /* initialize interrupt pointer */

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start frauniv_fdp_simulate_weights() ::" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize the memory list */

  MEMLIST_INIT( list );

  /* validate the argument list */

  err = localfn_fdp_simulate_input_check(
    delta, innovation_variance);
  if ( err ) return err;

  /* obtain the realization length */

  N = MATUNIV_NELEM( delta );

  /* allocate space for the output */

  err = matuniv_malloc( result, N, N, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate the weights for -1/2 <= delta < 1/2 */

  for ( t = 0; t < N; t++ ){

    /* register output matrix with the memory manager.*/

    err = memlist_member_register( &list, result, MEMTYPE_MATUNIV );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* allocate space for local variables */

    err = matdbl_malloc_register( &weights, t + 1, t + 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( &L, t + 1, t + 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( &mse, 1, t + 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = matdbl_malloc_register( &weight_row, 1, t + 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* if delta is in the nonstationary region,
    calculate the number of cumulative
    summations that we will have to perform
    over the corresponding stationary delta
    to properly generate the nonstationary
    weights */

    delta_t = delta->mat.dblmat.data[ t ];
    nsums   = 0;

    while ( delta_t >= 0.5 ){
      delta_t -= 1.0;
      nsums++;
    }

    /* develop partial autocorrelation
    sequence for FD process */

    err = localfn_fdp_Cholesky_matrix(
      delta_t, &L );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* calculate the mean square prediction
    error for an FD process */

    err = localfn_mse_prediction(
      delta_t,
      innovation_variance->mat.dblmat.data[ t ],
      &L,
      &mse );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* take the sqrt of the mean square prediction
    error to form mean prediction error */

    for ( k = 0; k < ( t + 1 ); k++ ){

      mse.data[ k ] = sqrt( mse.data[ k ] );
    }

    /* create the weight matrix for the current time ... */

    /* ... in the case where we have to perform cumulative
    summations to generate weights for a nonstationary
    FD parameter (delta), we need to create all rows of the
    weight matrix using the correpsonding stationary
    delta and then do the cumulative summations down each
    column. if, however, delta is already in the stationary
    range then we need only generate the final row of the
    weight matrix since no further processing is needed. The
    line below tests nsums, avariable which will have a
    positive value if delta is nonstationary. */

    if ( nsums > 0 ) kstart = 0;
    else kstart = t;

    for ( k = kstart; k <= t; k++ ){

    /* Solve for the current row in the
    weight matrix by solving a
    linear system of equations
    defined by the matrix relation:

      L * r = i

        where L is an t x t Levinson-Durbin
        matrix for an FD(delta[t], innovation_variance[t])
        process, r is an Nx1 vector of weights,
        and i is an tx1 vector whose p-th element
        is unity while all other elements are zero.
        Since the L matrix is upper triangular by
        design, we will use a simple yet efficient back
        substitution solver.

          Note: the r vector is not normalized and
          needs to be multiplied by the square root
          of the mle vector to calculate the true
          weights.
      */

      err = localfn_upper_triangular_back_substitution(
        &L, k, &weight_row );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* normalize the weight vector */

      err = matdbl_multiply_elem( &weight_row, &mse, intrp_ptr, &weight_row );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* copy the solution of the current row into the corresponding row
      of the weight matrix */

      err = matdbl_assign_submat( &weight_row, k, 0, intrp_ptr, &weights );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    /* form cumulative summations of stationary
    weights for nonstationary delta weights */

    for ( k = 0; k < nsums; k++ ){

      err = localfn_cumsum_lower_triangle( &weights );
      if ( err ) return err;
    }

    /* assign the last row of the local weights
    matrix to the global weights matrix */

    for ( k = 0; k < ( t + 1 ); k++ ){

      result->mat.dblmat.data[ t * N + k ] =
        weights.data[ t * ( t + 1 ) + k ];
    }

    /* free memory list (this is done for each loop
    and is faster than reallocating the memory).
    be sure to first free the node (and NOT the
    associated data) which corresponds to the
    result matrix.  */

    err = memlist_member_unregister( result, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    MUTIL_FREE_WARN( memlist, &list );
    MEMLIST_INIT( list );

    /* Check for interrupts */

    if ( MUTIL_INTERRUPT( 3.0 * N, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );

      return MUTIL_ERR_INTERRUPT;
    }
  }

  /* explicitly zero out the upper triangle
  of the weights matrix */

  for ( k = 0; k < N; k++ ){

    for ( t = k + 1; t < N; t++ ){

      result->mat.dblmat.data[ k * N + t ] = 0.0;
    }
  }

  /* free the registered memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( ":: Finished frauniv_fdp_simulate_weights()" );

  return MUTIL_ERR_OK;
}

/***********************************************/
/* STATIC SIMULATION FUNCTION DEFINITIONS      */
/***********************************************/

/** Input argument list check for the frauniv\_fdp\_simulate()
 * and frauniv\_fdp\_simulate\_weights() functions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_fdp_simulate_input_check( &delta, &innovation_variance );#
 * @return Standard mutils error/OK code.
 * @param  delta        Pointer to a pre-allocated single-column or
 *                      single-row universal matrix of type MUTIL\_DOUBLE
 *                      representing the FD parameter.
 * @param  innovation_variance Pointer to a pre-allocated single-column or
 *                      single-row universal matrix of type MUTIL\_DOUBLE.
 *                      representing the FD innovation variance.
 *
 * @see frauniv_fdp_simulate
 * @see frauniv_fdp_simulate_weights
 * @private
 */
static mutil_errcode localfn_fdp_simulate_input_check(
  const univ_mat *delta,
  const univ_mat *innovation_variance)
{
  mutil_errcode err;

  MUTIL_TRACE( "Start localfn_fdp_simulate_input_check()" );

  /*** check delta series ... ***/

  /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( delta, univ_mat, matuniv );

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( &(delta->mat.dblmat) ) ){
    MUTIL_ERROR( "Delta series matrix must be a "
                 "single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if number of elements is positive */

  if ( MATUNIV_NELEM( delta ) < 1 ){
    MUTIL_ERROR( "Number of elements in delta series "
                 "must be positive." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if the type is double */

  if ( delta->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Delta series matrix must be of type "
                 "MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /*** check innovation_variance series ... ***/

   /* ... for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( innovation_variance,
                                   univ_mat, matuniv );

  /* ... to see if it is a vector */

  if ( !MATANY_IS_VEC( &(innovation_variance->mat.dblmat) ) ){
    MUTIL_ERROR( "Innovation variance series matrix must be "
                 "a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if number of elements is positive */

  if ( MATUNIV_NELEM( innovation_variance ) < 1 ){
    MUTIL_ERROR( "Number of elements in innovation_variance "
                 "series must be positive." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if the type is double */

  if ( innovation_variance->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Innovation variance series matrix must be "
                 "of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* perform a cross length check on the two FD series */

  if ( MATUNIV_NELEM( innovation_variance ) !=
       MATUNIV_NELEM( delta ) ){
    MUTIL_ERROR( "The innovation variance series matrix "
                 "must have the same number of elements "
                 " as does the delta series matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  MUTIL_TRACE( "Done with localfn_fdp_simulate_input_check()" );

  return MUTIL_ERR_OK;
}

/** Transpose of the lower triangular Choleksy matrix
 * built with the partial autocorrelation sequence for an FD process.
 *
 * Creates an upper triangular matrix whose diagonal elements
 * are set to unity and whose upper triangle is filled
 * with Levinson-Durbin recursions of the partial autocorrelation
 * sequence (PACS) for a fractionally differenced (FD) process.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_fdp_Cholesky_matrix( delta, &result );#
 * @return Standard mutils error/OK code.
 * @param  delta        The FD parameter.
 * @param  result       Pointer to a pre-allocated square double matrix.
 *
 * @see frauniv_fdp_simulate_weights
 * @private
 */
static mutil_errcode localfn_fdp_Cholesky_matrix(
  double delta,
  double_mat *result )
{
  double        den;
  mutil_errcode err;
  sint32        N;
  sint32        k;
  sint32        row;
  sint32        t;

  /* check sizes */

  err = matdbl_validate( result );
  if ( err ) return err;

  N = result->nrow;

  if ( N != result->ncol ){
    MUTIL_ERROR( "The result matrix must be square." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* form the partial autocorrelation sequence
     and put into the first row of the result
     matrix starting from the second column */

  for ( t = 0; t < N - 1; t++ ){

    den = (double) t + 1.0 - delta;

    if ( den != 0.0 ){
      result->data[ t + 1 ] = delta / den;
    }
    else return MUTIL_ERR_DIVIDE_BY_ZERO;
  }

  /* form the upper triangle of the result
     matrix, sans the negative sign */

  for ( t = 1; t < N - 1; t++ ){

    for( k = 0; k < t; k++ ){

      result->data[ INDEX( t - k, t + 1, N ) ] =
      result->data[ INDEX( t - 1 - k, t, N ) ] -
      result->data[ t + 1 ] *
      result->data[ INDEX( k, t, N ) ];
    }
  }

  /* fill in the lower triangle with
     zeros, the diagonal with ones,
     and negate the entire upper
     triangle */

  for ( t = 0; t < N; t++ ){

    row = t * N;

    for( k = 0; k < N; k++ ){

      if ( k == t ){
        result->data[ row + k ] = 1.0;
      }
      else if ( k < t ){
        result->data[ row + k ] = 0.0;
      }
      else{
        result->data[ row + k ] = - result->data[ row + k ];
      }
    }
  }

  return MUTIL_ERR_OK;
}

/** Mean square prediction error for an FD process.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_mse_prediction( delta, &L, &result );#
 * @return Standard mutils error/OK code.
 * @param  delta        The FD parameter.
 * @param  innovation_variance The FD innovation variance.
 * @param  L            Pointer to a pre-allocated square double matrix
 *                      containing the Levinson-Durbin matrix.
 * @param  result       Pointer to a pre-allocated single-row or
 *                      single-column double matrix. The number of
 *                      elements in the matrix must be equal to the
 *                      number of rows (or columns) in the L matrix.
 *
 * @see frauniv_fdp_simulate_weights
 * @see localfn_fdp_Cholesky_matrix
 * @private
 */
static mutil_errcode localfn_mse_prediction(
  double      delta,
  double      innovation_variance,
  double_mat *L,
  double_mat *result )
{
  double        cumprod = 1.0;
  sint32        t;
  sint32        N;
  mutil_errcode err;

  /* validate matrices */

  err = matdbl_validate( L );
  if ( err ) return err;

  err = matdbl_validate( result );
  if ( err ) return err;

  /* check sizes */

  if ( !MATANY_IS_VEC( result ) ){
    MUTIL_ERROR( "Result matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  N = result->nelem;

  if ( ( L->nrow != N ) || ( L->ncol != N ) ){
    MUTIL_ERROR( "The L matrix must be square and have the same "
                 "number of rows or columns as the solution vector." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* calculate the mean square prediction error */

  for ( t = 0; t < N; t++ ){

    if ( t == 0 ){

      /* calculate the process variance and store
         as first element in the result matrix */

      result->data[ t ] = innovation_variance *
      MUTIL_GAMMA( 1.0 - 2 * delta ) /
      MUTIL_POW( MUTIL_GAMMA( 1.0 - delta ), 2.0 );

      cumprod = 1.0;
    }
    else{

      cumprod *= ( 1.0 - MUTIL_POW( L->data[ t ], 2.0 ) );

      result->data[ t ] = result->data[ 0 ] * cumprod;
    }
  }

  return MUTIL_ERR_OK;
}

/** Solution of a system of linear algebraic equations through back
 * substitution.
 * Solves for the vector x in the equation A * x = b, where A is an
 * upper triangular square matrix and b is an N element solution
 * vector. The solution vector is defined to have all zeros
 * with the exception of one element which is unity. The only nonzero
 * element is specified by the user via an input argument.
 * There are no upper-triangular checks made on the A matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_upper_triangular_back_substitution( &A, index, &x );#
 * @return Standard mutils error/OK code.
 * @param  A            Pointer to a pre-allocated NxN double matrix.
 * @param  index        The index of the only nonzero element of the solution
 *                      vector (index base is zero).
 * @param  result       Pointer to a pre-allocated single column or row double
 *                      matrix with N elements.
 *
 * @see frauniv_fdp_simulation_weights
 * @private
 */
static mutil_errcode localfn_upper_triangular_back_substitution(
  double_mat *A,
  sint32      index,
  double_mat *result)
{
  double        backsum;
  mutil_errcode err;
  sint32        N;
  sint32        i;
  sint32        j;
  sint32        row;

  /* validate matrices */

  err = matdbl_validate( A );
  if ( err ) return err;

  err = matdbl_validate( result );
  if ( err ) return err;

  /* check sizes */

  if ( !MATANY_IS_VEC( result ) ){
    MUTIL_ERROR( "Result matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  N = result->nelem;

  if ( ( A->nrow != N ) || ( A->ncol != N ) ){
    MUTIL_ERROR( "The A matrix must be square and have the same "
                 "number of rows or columns as the solution vector." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( index > N - 1 ){
    MUTIL_ERROR( "The nonzero index of the solution matrix must be "
                 "less than the number of rows (or columns) in "
                 "the left matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* solve the linear system via back-substitution */

  for ( i = N - 1; i >= 0; i-- ){

    if ( i < index ){

      backsum = 0.0;
      row     = i * N;

      for ( j = i + 1; j < N; j++ ){
        backsum += A->data[ row + j ] * result->data[ j ];
      }

      if ( A->data[ row + i ] != 0.0 ){
        result->data[ i ]  = - backsum / A->data[ row + i ];
      }
      else return MUTIL_ERR_DIVIDE_BY_ZERO;
    }
    else if( i > index ) result->data[ i ] = 0.0;
    else result->data[ i ] = 1.0 / A->data[ i * ( N + 1 ) ];
  }

  return MUTIL_ERR_OK;
}

/** Generate a series of normally distributed random deviates with
 * zero mean and unit variance.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_noise_gaussian( &result );#
 * @return Standard mutils error/OK code.
 * @param  result       Pointer to pre-allocated double matrix
 *                      of any size. The number of elements
 *                      of this matrix will be used to
 *                      determine the number of deviates to produce.
 *
 * @see frauniv_fdp_simulation
 * @private
 */
static mutil_errcode localfn_noise_gaussian(
  double_mat *result )
{
  mutil_errcode  err;
  sint32         t;
  void          *rand;

  /* check result for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_POINTER_FDP( result, double_mat, matdbl );

  /* .. initiate random number generator */

  err = mutil_rand_begin( &rand );
  if ( err ) return err;

  /* ... generate a Gaussian random deviate,
     one for each point in the realization */

  for ( t = 0; t < result->nelem; t++ ){
    err = mutil_rand_normal( rand, result->data + t );
    if ( err ) return err;
  }

  /* ... stop the random number generator */

  err = mutil_rand_end( rand );
  if ( err ) return err;

  return MUTIL_ERR_OK;
}

/** In-place cumumlative summation for each column
 * of a lower triangular matrix.
 *
 * For each column of a lower triangular matrix, a
 * cumulative sum is calculated over the rows. The results
 * replace the original matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_fdp.c
 * @library wavelets
 * @usage #err = localfn_cumsum_lower_triangle( &matrix );#
 * @return Standard mutils error/OK code.
 * @param  matrix            Pointer to a pre-allocated
 *                           double matrix. The original
 *                           matrix will be replaced with the result.
 *
 * @see frauniv_fdp_simulation_weights
 * @private
 */
static mutil_errcode localfn_cumsum_lower_triangle(
  double_mat *matrix)
{
  double         sum;
  mutil_errcode  err;
  sint32         col;
  sint32         row;
  sint32         index;

  err = matdbl_validate( matrix );
  if ( err ) return err;

  for ( col = 0; col < matrix->ncol; col++ ){

    sum = 0.0;

    for ( row = col; row < matrix->nrow; row++ ){

      index = row * matrix->ncol + col;

      sum += matrix->data[ index ];

      matrix->data[ index ] = sum;
    }
  }

  return MUTIL_ERR_OK;
}
