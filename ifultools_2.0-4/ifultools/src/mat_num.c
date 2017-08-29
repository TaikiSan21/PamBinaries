
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_num.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "mat_num.h"

#include "mat_any.h"
#include "mat_assn.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "ut_debug.h"
#include "ut_math.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_limit.h"

#define TINY 1.0e-20

/* This file contains implementations of matrix numerical functions
    in mat_num.h.
    */

/* function documented in mat_num.h */
/* written by Qin Cai */
mutil_errcode matuniv_lu_decomposition( const univ_mat *in_mat,
  void *intrp_ptr, univ_mat *indx, univ_mat *out_mat )
{
    mutil_errcode err;

    MUTIL_TRACE( "Start matuniv_lu_decomposition()" );

    /* avoid lint warning */
    (void) whatssi;

    /* sanity checks */
    if( !in_mat ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }
    if( !indx ) {
        MUTIL_ERROR( "NULL pointer for index vector" );
        return MUTIL_ERR_NULL_POINTER;
    }
    if( !out_mat ) {
        MUTIL_ERROR( "NULL pointer for output matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }

    if( !MATUNIV_CHECK_TYPE( in_mat, out_mat ) ) {
        MUTIL_ERROR( "Data types of input and output matrices are different" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }

    if( indx->type != MUTIL_SINT32) {
        MUTIL_ERROR( "Data type of index vector must be sint32" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }

    switch( in_mat->type ) {
        case MUTIL_DOUBLE:
            err = matdbl_lu_decomposition( &(in_mat->mat.dblmat),
              intrp_ptr, &(indx->mat.s32mat), &(out_mat->mat.dblmat) );
            if( err ) {
                return err;
            }
            break;

            /* nothing but doubles available now */
        default:
            MUTIL_ERROR( "This matrix type is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "Done with matuniv_lu_decomposition()" );
    return MUTIL_ERR_OK;
}


/* function documented in mat_num.h */
/* written by Qin Cai */
mutil_errcode matdbl_lu_decomposition( const double_mat *in_mat,
  void *intrp_ptr, sint32_mat *indx, double_mat *out_mat )
{
    mutil_errcode err;

    sint32        nrow;
    sint32        i;
    sint32        j;
    sint32        k;
    sint32        imax;

    double_mat    vv;
    double        big;
    double *      out;
    double        sum;
    double        tmp;

    double        num_ops = 0;

    MUTIL_INTERRUPT_INIT( intrp_ptr );

    MUTIL_TRACE( "Start matdbl_lu_decomposition()" );

    err = matdbl_validate( in_mat );
    if( err ) {
        return err;
    }

    err = mats32_validate( indx );
    if( err ) {
        return err;
    }

    err = matdbl_validate( out_mat );
    if( err ) {
        return err;
    }

    nrow = in_mat->nrow;
    if( nrow != in_mat->ncol) {
        MUTIL_ERROR( "The input matrix must be a square matrix" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    if( indx->nrow != nrow || indx->ncol != 1 ) {
        MUTIL_ERROR( "The index must be a column vector with same num rows" );
        MUTIL_ERROR( " as the input matix" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    if( !MATANY_EQUAL_DIM( in_mat, out_mat )) {
        MUTIL_ERROR( "Input and output matrices must be same size" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* initialize out_mat */
    if( in_mat->data != out_mat->data ) {
      err = matdbl_assign( in_mat, intrp_ptr, out_mat );
      if( err ) {
        return err;
      }
    }
    out = out_mat->data;

    err = matdbl_malloc( &vv, nrow, 1 );
    if( err ) {
        return err;
    }

    /* Loop over rows to get the implicit scaling */
    for( i = 0; i < nrow; i++) {
        big = 0.0;
        for( j = 0; j < nrow; j++) {
            if( (tmp = fabs( out[ MATANY_INDEX( out_mat, i, j )])) > big ) {
              big = tmp;
            }
        }
        if( !big ) {
            MUTIL_ERROR( "Singular matrix in routine matdbl_lu_decomposition()" );
            MUTIL_FREE_WARN( matdbl, &vv );
            return MUTIL_ERR_SINGULAR_MATRIX;
        }

        /* No nonzero largest element */
        vv.data[i] = 1.0 / big; /* save the scaling */
    }

    for( j = 0; j < nrow; j++ ) {   /* Loop over columns of Crout's method */
        for( i = 0; i < j; i++ ) {  /* i<j portion */
            sum = out[ MATANY_INDEX( out_mat, i, j )];
            for( k = 0; k < i; k++ ) {
                sum -= out[ MATANY_INDEX( out_mat, i, k ) ] *
                  out[ MATANY_INDEX( out_mat, k, j )];
            }
            out[ MATANY_INDEX( out_mat, i, j )] = sum;
        }
        big = 0.0;
        imax = j;
        for( i = j; i < nrow; i++) { /* i>=j portion */
            sum = out[ MATANY_INDEX( out_mat, i, j )];
            for( k = 0; k < j; k++ ) {
              sum -= out[ MATANY_INDEX( out_mat, i, k ) ] *
                out[ MATANY_INDEX( out_mat, k, j )];
            }
            out[ MATANY_INDEX( out_mat, i, j ) ] = sum;
            tmp = vv.data[i] * fabs( sum ) ;
            if( tmp >= big ) {
                big = tmp;
                imax = i;
            }
        }
        /* do we need to interchange rows */
        if( j != imax ) {
            /* yes */
            for( k = 0; k < nrow; k++ ) {
                tmp = out[ MATANY_INDEX( out_mat, imax, k ) ];
                out[ MATANY_INDEX( out_mat, imax, k )] =
                  out[ MATANY_INDEX( out_mat, j, k ) ];
                out[ MATANY_INDEX( out_mat, j, k )] = tmp;
            }
            vv.data[ imax ] = vv.data[ j ];
        }
        indx->data[j] = imax;
        if( !out[ MATANY_INDEX( out_mat, j, j )] ) {
            out[ MATANY_INDEX( out_mat, j, j )] = TINY;
        }

        /* if the pivot element is zero, the matrix is singular (at least
           to the precision of the algorithm).
           For some application, it is
           desirable to subsitute TINY for zero */

        if( j != nrow - 1) {
            tmp = 1.0 / out[ MATANY_INDEX( out_mat, j, j )];
            for( i = j + 1; i < nrow; i++ ) {
                out[ MATANY_INDEX( out_mat, i, j ) ] *= tmp;
            }
        }

        num_ops += nrow * nrow * 20.0;
        if( MUTIL_INTERRUPT( num_ops, intrp_ptr )) {
          MUTIL_ERROR( "User interrupt" );
          return MUTIL_ERR_INTERRUPT;
        }
    } /* Go back for the next column in the reduction */

    MUTIL_FREE_WARN(matdbl, &vv);

    MUTIL_TRACE( "Done with matdbl_lu_decomposition()" );
    return MUTIL_ERR_OK;
}


/* function documented in mat_num.h */
/* written by Qin Cai */
mutil_errcode matuniv_lu_solve( const univ_mat *lu_mat,
  const univ_mat *in_vec, const univ_mat *indx, void *intrp_ptr,
  univ_mat *out_vec )
{
    mutil_errcode err;

    MUTIL_TRACE( "Start matuniv_lu_solve()" );

    /* sanity checks */
    if( !lu_mat ) {
        MUTIL_ERROR( "NULL pointer for input matrix" );
        return MUTIL_ERR_NULL_POINTER;
    }
    if( !in_vec ) {
        MUTIL_ERROR( "NULL pointer for input vector" );
        return MUTIL_ERR_NULL_POINTER;
    }
    if( !indx ) {
        MUTIL_ERROR( "NULL pointer for index vector" );
        return MUTIL_ERR_NULL_POINTER;
    }
    if( !out_vec ) {
        MUTIL_ERROR( "NULL pointer for output parameter vector" );
        return MUTIL_ERR_NULL_POINTER;
    }

    if( !MATUNIV_CHECK_TYPE( lu_mat, in_vec ) ) {
        MUTIL_ERROR( "Data types of input matrix and vector are different" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }

    if( indx->type != MUTIL_SINT32) {
        MUTIL_ERROR( "Data type of index vector must be sint32" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }
    if( !MATUNIV_CHECK_TYPE( lu_mat, out_vec ) ) {
        MUTIL_ERROR( "Data types of input matrix and output parameter vector " );
        MUTIL_ERROR( "are different" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }

    switch( lu_mat->type ) {
        case MUTIL_DOUBLE:
            err = matdbl_lu_solve( &(lu_mat->mat.dblmat),
              &(in_vec->mat.dblmat), &(indx->mat.s32mat), intrp_ptr,
              &(out_vec->mat.dblmat) );
            if( err ) {
              return err;
            }
            break;

        default:
            MUTIL_ERROR( "This matrix type is currently unsupported" );
            return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "Done with matuniv_lu_solve()" );

    return MUTIL_ERR_OK;
}


/* function documented in mat_num.h */
/* written by Qin Cai */
mutil_errcode matdbl_lu_solve( const double_mat *lu_mat,
  const double_mat *in_vec, const sint32_mat *indx,
  void *intrp_ptr, double_mat *out_vec )
{
    sint32        i;
    sint32        ii;
    sint32        j;
    sint32        nrow;

    sint32        perm_indx;

    double        sum;
    double *      outdat;
    double *      ludat;

    mutil_errcode err;
    double        num_ops = 0;

    MUTIL_INTERRUPT_INIT( intrp_ptr );

    MUTIL_TRACE( "Start matdbl_lu_solve()" );

    err = matdbl_validate( lu_mat );
    if( err ) {
        return err;
    }

    err = matdbl_validate( in_vec );
    if( err ) {
        return err;
    }

    err = mats32_validate( indx );
    if( err ) {
        return err;
    }

    err = matdbl_validate( out_vec );
    if( err ) {
        return err;
    }

    if( lu_mat->nrow != lu_mat->ncol ) {
        MUTIL_ERROR( "The input matrix must be a square matrix" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }
    nrow = lu_mat->nrow;

    if( in_vec->nrow != nrow || in_vec->ncol != 1 ) {
        MUTIL_ERROR( "The input vector must be a row vector whose length " );
        MUTIL_ERROR( "is equal to the number of rows of the input matix" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }
    if( !MATANY_EQUAL_DIM( in_vec, indx ) ) {
        MUTIL_ERROR( "The index and input vectors must have the same dimension" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }
    if( !MATANY_EQUAL_DIM( in_vec, out_vec ) ) {
        MUTIL_ERROR( "The input and output vectors must have the same dimension" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* initialize the output vector */
    if( in_vec->data != out_vec->data ) {
      err = matdbl_assign( in_vec, intrp_ptr, out_vec );
      if( err ) {
        return err;
      }
    }
    outdat = out_vec->data;
    ludat  = lu_mat->data;

    /* when ii is set to a positive value, it will become the index
       of the first nonvanishing element of in_vec. We now do the forward
       substitution. The only new wrinkle is to un-scramble the permutation
       as we go.
    */
    ii = -1;
    for( i = 0; i < nrow; i++ ) {
        perm_indx      = indx->data[ i ];
        sum            = outdat[ perm_indx ];
        outdat[perm_indx] = outdat[ i ];

        if( ii >= 0 ) {
            for( j = ii; j <= i - 1; j++ ) {
                sum -= ludat[ MATANY_INDEX( lu_mat, i, j )] * outdat[j];
            }
        }
        else if( fabs( sum ) > TINY ) {
            ii = i;
        }
        outdat[i] = sum;
    }

    num_ops += nrow * ( nrow * 3.0 + 5.0 );
    if( MUTIL_INTERRUPT( num_ops, intrp_ptr )) {
        MUTIL_ERROR( "User interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    for( i = nrow - 1; i >= 0; i-- ) {
        sum = outdat[i];
        for( j = i+1; j < nrow; j++ ) {
            sum -= ludat[ MATANY_INDEX( lu_mat, i, j ) ] * outdat[j];
        }
        if( !ludat[ MATANY_INDEX( lu_mat, i, i )] ) {
          MUTIL_ERROR( "Illegal LU decomposition, cannot solve" );
          return MUTIL_ERR_ILLEGAL_VALUE;
        }
        outdat[i] = sum / ludat[ MATANY_INDEX( lu_mat, i, i )];
    }

    num_ops += nrow * nrow * 3.0;
    if( MUTIL_INTERRUPT( num_ops, intrp_ptr )) {
        MUTIL_ERROR( "User interrupt" );
        return MUTIL_ERR_INTERRUPT;
    }

    MUTIL_TRACE( "Done with matdbl_lu_solve()" );

    return MUTIL_ERR_OK;
}


/* function documented in mat_num.h */
/* written by Qin Cai */
mutil_errcode matuniv_inverse( const univ_mat *in_mat, void *intrp_ptr,
  univ_mat *out_mat )
{
    mutil_errcode err;

    MUTIL_TRACE( "Start matuniv_inverse()" );

    /* sanity checks */

    if( !in_mat ) {
      MUTIL_ERROR( "NULL pointer for input matrix" );
      return MUTIL_ERR_NULL_POINTER;
    }

    if( !out_mat ) {
      MUTIL_ERROR( "NULL pointer for output matrix" );
      return MUTIL_ERR_NULL_POINTER;
    }

    if( !MATUNIV_CHECK_TYPE( in_mat, out_mat ) ) {
      MUTIL_ERROR( "Data types of input and output matrices are different" );
      return MUTIL_ERR_ILLEGAL_TYPE;
    }

    switch( in_mat->type ) {
      case MUTIL_DOUBLE:
        err = matdbl_inverse( &(in_mat->mat.dblmat), intrp_ptr,
          &(out_mat->mat.dblmat) );
        if( err ) {
          return err;
        }
        break;

      default:
        MUTIL_ERROR( "This matrix type is currently unsupported" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }

    MUTIL_TRACE( "Done with matuniv_inverse()" );
    return MUTIL_ERR_OK;
}


/* function documented in mat_num.h */
/* written by Qin Cai */
mutil_errcode matdbl_inverse( const double_mat *in_mat,
  void *intrp_ptr, double_mat *out_mat )
{
    mutil_errcode err;

    sint32_mat    indx;
    double_mat    tmp_vec;
    double_mat    lu_mat;

    sint32        nrow;
    sint32        col;

    MUTIL_TRACE( "Start matdbl_inverse()" );

    err = matdbl_validate( in_mat );
    if( err ) {
        return err;
    }

    err = matdbl_validate( out_mat );
    if( err ) {
        return err;
    }

    if( in_mat->nrow != in_mat->ncol ) {
        MUTIL_ERROR( "The input matrix must be a square matrix" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }
    nrow = in_mat->nrow;

    if( !MATANY_EQUAL_DIM( in_mat, out_mat ) ) {
        MUTIL_ERROR( "The input and output must have the same dimension" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    /* allocate the space for the lu decomposition and index  */

    err = matdbl_malloc( &lu_mat, nrow, nrow );
    if( err ) {
        return err;
    }

    err = mats32_malloc( &indx, nrow, 1 );
    if( err ) {
        MUTIL_FREE_WARN( matdbl, &lu_mat );
        return err;
    }

    /* perform LU decomposition */

    err = matdbl_lu_decomposition( in_mat, intrp_ptr, &indx, &lu_mat );
    if( err ) {
        MUTIL_FREE_WARN( matdbl, &lu_mat );
        MUTIL_FREE_WARN( mats32, &indx );
        return err;
    }

    /* allocate the space for a column vector */

    err = matdbl_malloc( &tmp_vec, nrow, 1 );
    if( err ) {
        MUTIL_FREE_WARN( matdbl, &lu_mat );
        MUTIL_FREE_WARN( mats32, &indx );
        return err;
    }

#define CLEANUP \
   MUTIL_FREE_WARN(matdbl, &lu_mat);\
   MUTIL_FREE_WARN(mats32, &indx); \
   MUTIL_FREE_WARN(matdbl, &tmp_vec)

    /* call matdbl_lu_solve column by column */

    for( col = 0; col < nrow; col++ ) {

        /* set up vector with 1 non-zero element */

        err = matdbl_assign_scalar( 0.0, intrp_ptr, &tmp_vec );
        if( err ) {
            CLEANUP;
            return err;
        }
        tmp_vec.data[col] = 1.0;

        /* solve LU * vec = out */

        err = matdbl_lu_solve( &lu_mat, &tmp_vec, &indx, intrp_ptr, &tmp_vec );
        if( err ) {
            CLEANUP;
            return err;
        }

        /* put this vector of values into the output matrix */

        err = matdbl_assign_submat( &tmp_vec, 0, col, intrp_ptr, out_mat );
        if( err ) {
            CLEANUP;
            return err;
        }

    } /* end of loop over columns */

    CLEANUP;

    /* note that each function in loop had interrupts so we do
       not have to check again */

    MUTIL_TRACE( "Done with matdbl_inverse()" );
    return MUTIL_ERR_OK;
}


#define EIGEN_ROTATE( mymat, size, i, j, k, l ) \
    gval = (mymat)->mat.dblmat.data[ i * size + j ]; \
    hval = (mymat)->mat.dblmat.data[ k * size + l ]; \
    (mymat)->mat.dblmat.data[ i * size + j ] = gval - \
      sval * ( hval + gval * tau ); \
    (mymat)->mat.dblmat.data[ k * size + l ] = hval + \
      sval * ( gval - hval * tau )

/* function documented in mat_num.h */
/* written by Vikram Chalana */
mutil_errcode matuniv_eigen_jacobi( const univ_mat *in_mat,
  void *intrp_ptr, univ_mat *eigenvalues, univ_mat *eigenvectors )
{
  mutil_errcode err;
  sint32        row;
  sint32        col;
  sint32        idx;
  sint32        size;

  univ_mat      amat;
  univ_mat      bvec;
  univ_mat      zvec;

  double        val;
  sint32        nrot;
  sint32        iter;

  double        sum;
  double        thresh;
  double        gval;
  double        theta;
  double        tau;
  double        tval;
  double        sval;
  double        hval;
  double        cval;

  double        num_ops = 0;
  
  MUTIL_INTERRUPT_INIT( intrp_ptr );



  /* Sanity checks */

  err = matuniv_validate( in_mat );
  if( err ) return err;

  err = matuniv_validate( eigenvalues );
  if( err ) return err;

  err = matuniv_validate( eigenvectors );
  if( err ) return err;

  /* Only matrices of type doubles supported now */
  if ( ( in_mat->type != MUTIL_DOUBLE ) ||
       ( eigenvalues->type != MUTIL_DOUBLE ) ||
       ( eigenvectors->type != MUTIL_DOUBLE ) ) {
    MUTIL_ERROR( "Input and output matrix type must be double" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* Check sizes */
  if( !MATANY_EQUAL_DIM( &in_mat->mat.dblmat, &eigenvectors->mat.dblmat ) ) {
    MUTIL_ERROR( "The input and output matrices must have the "
      "same dimension" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if( MATUNIV_NROW( in_mat ) != MATUNIV_NCOL( in_mat ) ) {
    MUTIL_ERROR( "Only square input matrices allowed" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  size = MATUNIV_NROW( in_mat );
  if ( MATUNIV_NELEM( eigenvalues ) != size ) {
    MUTIL_ERROR( "Number of elements of eigenvector should be same as "
      "number of rows in input matrix" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }


  /* Initialize the output matrix to identity */
  for ( row = 0; row < size ; row++ ) {
    for ( col = 0; col < size; col++ ) {
      eigenvectors->mat.dblmat.data[ row * size + col ] = 0;
    }
    eigenvectors->mat.dblmat.data[ row * size + row ] = 1;
  }

  /* Allocate space for a copy of in_mat and make a copy */

  err = matuniv_malloc( &amat, size, size, MUTIL_DOUBLE );
  if ( err ) {
    return err;
  }

  /* amat = copy of in_mat */
  err = matuniv_assign( in_mat, intrp_ptr, &amat );
  if ( err ) {
    MUTIL_FREE_WARN( matuniv, &amat);
    return err;
  }

  /* Allocate space for two temporary storage matrices */
  err = matuniv_malloc( &bvec, size, 1, MUTIL_DOUBLE );
  if ( err ) {
    MUTIL_FREE_WARN( matuniv, &amat);
    return err;
  }
  err = matuniv_malloc( &zvec, size, 1, MUTIL_DOUBLE );
  if ( err ) {
    MUTIL_FREE_WARN( matuniv, &amat);
    MUTIL_FREE_WARN( matuniv, &bvec);
    return err;
  }

  /* Initialize the various vectors */
  for ( idx = 0; idx < size; idx++ ) {

    /* Get the diagonal entries */
    val = amat.mat.dblmat.data[ idx * size +  idx ];
    eigenvalues->mat.dblmat.data[ idx ] = val;
    bvec.mat.dblmat.data[ idx ]         = val;
    zvec.mat.dblmat.data[ idx ]         = 0.0;
  }

  num_ops += ( size * 3 );
  if( MUTIL_INTERRUPT( num_ops, intrp_ptr )) {
    MUTIL_ERROR( "User interrupt" );
    return MUTIL_ERR_INTERRUPT;
  }


  nrot = 0;
  for ( iter = 0; iter < 50; iter ++ ) {

    /* Compute the sum of off-diagonal entries */
    sum = 0.0;
    for ( row = 0; row < size - 1; row++ ) {
      for ( col = row + 1; col < size; col++ ) {
        sum += fabs( amat.mat.dblmat.data[ row * size + col ] );
      }
    }

    /* The normal return relies of quadratic convergence to machine
    underflow */
    if ( sum < MUTIL_DOUBLE_EPSILON ) {
      MUTIL_FREE_WARN( matuniv, &amat);
      MUTIL_FREE_WARN( matuniv, &bvec);
      MUTIL_FREE_WARN( matuniv, &zvec);
      return MUTIL_ERR_OK;
    }

    if ( iter < 3 )
      thresh = ( 0.2 * sum ) / (size * size);
    else
      thresh=0.0;

    for ( row = 0; row < size - 1; row++ ) {
      for ( col = row + 1; col < size; col++ ) {
        gval = 100.0 * fabs( amat.mat.dblmat.data[ row * size + col ] );

        /* After four sweeps,
           skip the rotation if the off-diagonal element is small */
        if ( ( iter > 3 ) &&
          MUTIL_EQUAL_TO( ( fabs( eigenvalues->mat.dblmat.data[ row ] ) + gval ),
            fabs( eigenvalues->mat.dblmat.data[ row ] ), 4 ) &&
          MUTIL_EQUAL_TO ( ( fabs( eigenvalues->mat.dblmat.data[ col ] ) + gval ),
            fabs( eigenvalues->mat.dblmat.data[ col ] ), 4 ) ) {

          amat.mat.dblmat.data[ row * size + col ] = 0.0;
        }
        else if ( fabs( amat.mat.dblmat.data[ row * size + col ] ) > thresh) {

          hval = eigenvalues->mat.dblmat.data[ col ] -
            eigenvalues->mat.dblmat.data[ row ];

          if ( MUTIL_EQUAL_TO( fabs( hval ) + gval, fabs( hval ), 4 ) ) {
            tval = amat.mat.dblmat.data[ row * size + col ] / hval;
          }
          else {
            theta = 0.5 * hval / amat.mat.dblmat.data[ row * size + col ];
            tval  = 1.0 / ( fabs( theta ) + sqrt( 1.0 + theta * theta) );
            if ( theta < 0.0 )
              tval = -tval;
          }

          cval = 1.0 / sqrt( 1 + tval * tval );
          sval = tval * cval;
          tau  = sval / ( 1.0 + cval );
          hval = tval * amat.mat.dblmat.data[ row * size + col ];

          zvec.mat.dblmat.data[ row ] -= hval;
          zvec.mat.dblmat.data[ col ] += hval;

          eigenvalues->mat.dblmat.data[ row ] -= hval;
          eigenvalues->mat.dblmat.data[ col ] += hval;

          amat.mat.dblmat.data[ row * size + col ] = 0.0;

          /* Case 0 <= idx < row */
          for ( idx = 0; idx <= row - 1; idx++ ) {
            EIGEN_ROTATE( &amat, size, idx, row, idx, col );
          }

          /* Case row < idx < col */
          for ( idx = row + 1; idx <= col-1; idx++ ) {
            EIGEN_ROTATE( &amat, size, row, idx, idx, col );
          }

          /* Case col < idx < size */
          for ( idx = col + 1; idx < size; idx++ ) {
            EIGEN_ROTATE( &amat, size, row, idx, col, idx );
          }

          for ( idx = 0; idx < size; idx++ ) {
            EIGEN_ROTATE( eigenvectors, size, idx, row, idx, col );
          }
          ++nrot;
        } /* if ( fabs( amat.mat.dblmat.data[ row * size + col ] ) > thresh) */
      } /* for ( col = row + 1; col < size; col++ )  */
    } /* for ( row = 0; row < size - 1; row++ ) */

    for ( idx = 0; idx < size; idx++ ) {
      bvec.mat.dblmat.data[ idx ]        += zvec.mat.dblmat.data[ idx ];
      eigenvalues->mat.dblmat.data[ idx ] = bvec.mat.dblmat.data[ idx ];
      zvec.mat.dblmat.data[ idx ]         = 0.0;
    }

    num_ops += ( size * size * 100 );
    if( MUTIL_INTERRUPT( num_ops, intrp_ptr )) {
      MUTIL_ERROR( "User interrupt" );
      return MUTIL_ERR_INTERRUPT;
    }
  } /* for ( iter = 0; iter < 50; iter ++ ) */

  MUTIL_FREE_WARN( matuniv, &amat);
  MUTIL_FREE_WARN( matuniv, &bvec);
  MUTIL_FREE_WARN( matuniv, &zvec);
  MUTIL_ERROR( "Too many iterations in routine jacobi" );
  return MUTIL_ERR_NOT_CONVERGING;
}


/* function documented in mat_num.h */
/* written by Vikram Chalana */
mutil_errcode matuniv_eigen_sort( const univ_mat *in_eigenvalues,
  const univ_mat *in_eigenvectors, void *intrp_ptr,
  univ_mat *out_eigenvalues, univ_mat *out_eigenvectors )
{
  mutil_errcode err;

  sint32        idx;
  sint32        jdx;
  sint32        kdx;
  sint32        fastindex;
  sint32        size;
  double        pval;

  double        *edata;

  double        num_ops = 0;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  /* Sanity checks */

  err = matuniv_validate( in_eigenvalues );
  if( err ) return err;
  err = matuniv_validate( in_eigenvectors );
  if( err ) return err;
  err = matuniv_validate( out_eigenvalues );
  if( err ) return err;
  err = matuniv_validate( out_eigenvectors );
  if( err ) return err;

  if ( ( in_eigenvalues->type   != MUTIL_DOUBLE ) ||
       ( in_eigenvectors->type  != MUTIL_DOUBLE ) ||
       ( out_eigenvalues->type  != MUTIL_DOUBLE ) ||
       ( out_eigenvectors->type != MUTIL_DOUBLE ) ) {
    MUTIL_ERROR( "Input and output matrix type must be double" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( !MATANY_EQUAL_DIM( &in_eigenvectors->mat.dblmat,
    &out_eigenvectors->mat.dblmat ) ) {
    MUTIL_ERROR( "The input and output eigenvector matrices must have "
      "the same dimension" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if( !MATANY_EQUAL_DIM( &in_eigenvalues->mat.dblmat,
    &out_eigenvalues->mat.dblmat ) ) {
    MUTIL_ERROR( "The input and output eigenvalue matrices must have "
      "the same dimension" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if( MATUNIV_NROW( in_eigenvectors ) != MATUNIV_NCOL( in_eigenvectors ) ) {
    MUTIL_ERROR( "Only square eigenvector matrices allowed" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  size = MATUNIV_NROW( in_eigenvectors );
  if ( MATUNIV_NELEM( in_eigenvalues ) != size ) {
    MUTIL_ERROR( "Number of elements of eigenvalues should be same as "
      "number of rows in input matrix" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* Copy inputs to outputs, and operate in-place after that */

  err = matuniv_assign( in_eigenvalues, intrp_ptr, out_eigenvalues );
  if ( err ) return err;

  err = matuniv_assign( in_eigenvectors, intrp_ptr, out_eigenvectors );
  if ( err ) return err;

  edata = out_eigenvalues->mat.dblmat.data;

  for ( idx = 0; idx < size; idx ++ ) {

    /* find largest eigenvalue */
    pval = edata[ idx ];
    kdx = idx;
    for ( jdx = idx + 1; jdx < size; jdx ++ ) {
      if ( edata[ jdx ] > pval ) {
        pval = edata[ jdx ];
        kdx = jdx;
      }
    }

    num_ops += 4.0 * ( size - idx );

    /* rearrange if necessary */
    if ( kdx != idx ) {

      edata[ kdx ] = edata[ idx ];
      edata[ idx ] = pval;

      for ( jdx = 0; jdx < size; jdx++ ) {
        fastindex = jdx * size;
        pval = out_eigenvectors->mat.dblmat.data[ fastindex + idx ];
        out_eigenvectors->mat.dblmat.data[ fastindex + idx ] =
          out_eigenvectors->mat.dblmat.data[ fastindex + kdx ];
        out_eigenvectors->mat.dblmat.data[ fastindex + kdx ] = pval;
      }

      num_ops += 12.0 * size;

    } /* if ( kdx != idx ) */

    if( MUTIL_INTERRUPT( num_ops, intrp_ptr )) {
      MUTIL_ERROR( "User interrupt" );
      return MUTIL_ERR_INTERRUPT;
    }

  } /* for ( idx ... ) */

  return MUTIL_ERR_OK;
}
