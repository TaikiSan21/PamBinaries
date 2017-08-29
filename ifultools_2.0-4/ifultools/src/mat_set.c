
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_set.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "mat_set.h"

#include "mat_assn.h"
#include "mat_cast.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"

/*
   This file contains the definitions for the functions for
   sets of matrices declared in mat_set.h.
*/


/* function documented in mat_set.h */
/* written by Jennifer Hodgdon */
mutil_errcode matset_malloc( mat_set *matset, sint32 ndim,
  const sint32 *dims )
{
  mutil_errcode trouble;
  sint32        dim_idx;

  MUTIL_TRACE( "Start matset_malloc()" );

  /* avoid lint warning */
  (void) whatssi;

  if(!matset) {
    MUTIL_ERROR( "NULL matrix set pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if(!dims) {
    MUTIL_ERROR( "NULL dimensions pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if( ndim <= 0 ) {
    MUTIL_ERROR( "Number of dimensions must be positive" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* allocate space for dims */

  trouble = mutil_malloc( ndim * sizeof( sint32 ), (void **) &(matset->dims));
  if( trouble ) return trouble;
  MUTIL_ASSERT( matset->dims );

  /* copy in dims and calculate total elems */

  matset->nelem = 1;
  matset->ndim  = ndim;

  for( dim_idx = 0; dim_idx < ndim; dim_idx++ ) {
    if( dims[dim_idx] <= 0 ) {
      MUTIL_FREE_BUFFER_WARN( matset->dims, ndim * sizeof( sint32 ));
      MUTIL_ERROR( "All dimensions must be positive" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    matset->dims[dim_idx] = dims[dim_idx];
    matset->nelem  *= dims[dim_idx];
  }

  MUTIL_ASSERT( matset->nelem > 0 );

  trouble = mutil_malloc( matset->nelem * sizeof( univ_mat ),
    (void **) &(matset->mats));
  if( trouble ) {
    MUTIL_FREE_BUFFER_WARN( matset->dims, ndim * sizeof( sint32 ));
    return trouble;
  }

  MUTIL_ASSERT( matset->mats );

  /* set contiguous flag to false */
  matset->contiguous = FALSE;


  MUTIL_TRACE( "Done with matset_malloc()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Jennifer Hodgdon */
mutil_errcode matset_malloc_matrices( mat_set *matset,
  sint32 nrow, sint32 ncol, mutil_data_type type )
{
  mutil_errcode trouble;
  sint32        elem_idx;
  sint32        inner_idx;
  univ_mat *    mat;

  MUTIL_TRACE( "Start matset_malloc_matrices()" );

  trouble = matset_validate( matset );
  if( trouble ) return trouble;

  for( elem_idx = 0; elem_idx < matset->nelem; elem_idx++ ) {
    mat = &(matset->mats[elem_idx]);
    trouble = matuniv_malloc( mat, nrow, ncol, type );
    if( trouble ) {

      /* try to recover memory */
      for( inner_idx = 0; inner_idx < elem_idx; inner_idx++ ) {
        mat = &(matset->mats[inner_idx]);
        MUTIL_FREE_WARN( matuniv, mat );
      }
      return trouble;
    } /* if trouble */
  } /* for loop to allocate all matrices */

  /* set contiguous flag to false */
  matset->contiguous = FALSE;

  MUTIL_TRACE( "Done with matset_malloc_matrices()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Lydia Ng */
mutil_errcode matset_malloc_matrices_contiguous( mat_set *matset,
  sint32 nrow, sint32 ncol, mutil_data_type type )
{
  mutil_errcode trouble;
  void *        dataptr;
  sint32        typesize;
  sint32        elem_idx;
  univ_mat *    mat;
  sint32        nmatelem;

  MUTIL_TRACE( "Start matset_malloc_matrices_contiguous()" );

  if( nrow <= 0 || ncol <= 0 ) {
    MUTIL_ERROR( "Number of rows or columns not positive" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  trouble = matset_validate( matset );
  if( trouble ) return trouble;

  switch ( type ) {
    case MUTIL_DOUBLE:
      typesize = sizeof( double );
      break;

    case MUTIL_FLOAT:
      typesize = sizeof( float );
      break;

    case MUTIL_UINT8:
      typesize = sizeof( uint8 );
      break;

    case MUTIL_UINT16:
      typesize = sizeof( uint16 );
      break;

    case MUTIL_UINT32:
      typesize = sizeof( uint32 );
      break;

    case MUTIL_SINT16:
      typesize = sizeof( sint16 );
      break;

    case MUTIL_SINT32:
      typesize = sizeof( sint32 );
      break;

    case MUTIL_DCOMPLEX:
      typesize = sizeof( dcomplex );
      break;

    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* allocated contiguous block of memory */
  nmatelem = nrow * ncol;
  trouble = mutil_malloc( nmatelem * matset->nelem * typesize,
    (void **) &dataptr );
  if( trouble ) return trouble;

  MUTIL_ASSERT( dataptr );

  /* initialize the matrices with the allocated memory */
  for( elem_idx = 0; elem_idx < matset->nelem; elem_idx++ ) {

    mat = &(matset->mats[elem_idx]);
    trouble = matuniv_wrap_data( mat, dataptr, nrow, ncol, type );

    /* try to recover memory */
    if( trouble ) {
      MUTIL_FREE_BUFFER_WARN( dataptr,
        nrow * ncol * matset->nelem * typesize );
      return trouble;
    }

    /* advance pointer to the next matrix */
    switch( type ) {
      case MUTIL_DOUBLE:
        dataptr = (void *) ( ( (double *) dataptr ) + nmatelem );
        break;

      case MUTIL_FLOAT:
        dataptr = (void *) ( ( (float *) dataptr ) + nmatelem );
        break;

      case MUTIL_UINT8:
        dataptr = (void *) ( ( (uint8 *) dataptr ) + nmatelem );
        break;

      case MUTIL_UINT16:
        dataptr = (void *) ( ( (uint16 *) dataptr ) + nmatelem );
        break;

      case MUTIL_UINT32:
        dataptr = (void *) ( ( (uint32 *) dataptr ) + nmatelem );
        break;

      case MUTIL_SINT16:
        dataptr = (void *) ( ( (sint16 *) dataptr ) + nmatelem );
        break;

      case MUTIL_SINT32:
        dataptr = (void *) ( ( (sint32 *) dataptr ) + nmatelem );
        break;

      case MUTIL_DCOMPLEX:
        dataptr = (void *) ( ( (dcomplex *) dataptr ) + nmatelem );
        break;

      default:
        MUTIL_ERROR( "This matrix type is currently unsupported" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }

  } /* for loop to wrap all matrices */

  /* set contiguous flag to true */
  matset->contiguous = TRUE;

  MUTIL_TRACE( "Done with matset_malloc_matrices_contiguous()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Luca Cazzanti */
mutil_errcode matset_malloc_matrices_arbitrary_size( mat_set *matset,
    const sint32_mat *nrow_mat, const sint32_mat *ncol_mat,
    mutil_data_type type )
{
  mutil_errcode trouble;
  sint32        elem_idx;
  sint32        inner_idx;
  univ_mat *    mat;

  MUTIL_TRACE( "Start matset_malloc_matrices_arbitrary_size()" );

  trouble = matset_validate( matset );
  if( trouble ) return trouble;

  trouble = mats32_validate( nrow_mat );
  if( trouble ) return trouble;

  trouble = mats32_validate( ncol_mat );
  if( trouble) return trouble;

  if( nrow_mat->nelem != ncol_mat->nelem ) {
      MUTIL_ERROR( "Matrices containing dimensions must have the same size" );
      return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if( matset->nelem != nrow_mat->nelem ) {
      MUTIL_ERROR( "Number of matrices does not match number of dimensions" );
      return MUTIL_ERR_ILLEGAL_SIZE;
  }

  for( elem_idx = 0; elem_idx < matset->nelem; elem_idx++ ) {
    mat = &(matset->mats[elem_idx]);
    trouble = matuniv_malloc( mat, nrow_mat->data[elem_idx],
        ncol_mat->data[elem_idx], type );
    if( trouble ) {

      /* try to recover memory */
      for( inner_idx = 0; inner_idx < elem_idx; inner_idx++ ) {
        mat = &(matset->mats[inner_idx]);
        MUTIL_FREE_WARN( matuniv, mat );
      }
      return trouble;
    } /* if trouble */
  } /* for loop to allocate all matrices */

  /* set contiguous flag to false */
  matset->contiguous = FALSE;

  MUTIL_TRACE( "Done with matset_malloc_matrices_arbitrary_size()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Jill Goldschneider */
mutil_errcode matset_malloc_flat(
  mat_set         *matset,
  sint32           nelem,
  sint32           nrow,
  sint32           ncol,
  mutil_data_type  type,
  boolean          contiguous )
{
  mutil_errcode err;

  MUTIL_TRACE( "Start matset_malloc_flat()" );

  err = matset_malloc( matset, 1, &nelem );
  if ( err ) return err;

  if ( contiguous ) {
    err = matset_malloc_matrices_contiguous( matset, nrow, ncol, type );
  }
  else {
    err = matset_malloc_matrices( matset, nrow, ncol, type );
  }
  if ( err ) {
    MUTIL_FREE_WARN( matset, matset );
    return err;
  }

  MUTIL_TRACE( "Done with matset_malloc_flat()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Jennifer Hodgdon */
mutil_errcode matset_free( mat_set *matset )
{
  mutil_errcode trouble;

  MUTIL_TRACE( "Start matset_free()" );

  if( !matset ) {
    MUTIL_ERROR( "NULL matrix set pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }

  trouble = mutil_free( (void *) matset->mats,
    matset->nelem * sizeof( univ_mat ));
  if( trouble ) return trouble;
  matset->mats = NULL;
  matset->nelem = 0;

  trouble = mutil_free( matset->dims, matset->ndim * sizeof( sint32 ));
  if( trouble ) return trouble;

  matset->dims = NULL;
  matset->ndim = 0;

  MUTIL_TRACE( "Done with matset_free()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Jennifer Hodgdon, Lydia Ng */
mutil_errcode matset_matrices_free( mat_set *matset )
{
  mutil_errcode trouble;
  sint32        elem_idx;
  sint32        typesize;
  sint32        totalelem;
  univ_mat *    mat;

  MUTIL_TRACE( "Start matset_matrices_free()" );

  trouble = matset_validate( matset );
  if ( trouble ) return trouble;

  if ( matset->contiguous == FALSE ) {

    MUTIL_TRACE( "Freeing non-contiguous matset memory" );

    for( elem_idx = 0; elem_idx < matset->nelem; elem_idx++ ) {
      mat = &(matset->mats[elem_idx]);
      trouble = matuniv_free( mat );
      if ( trouble ) return trouble;
    }
  }

  else {

    MUTIL_TRACE( "Freeing contiguous matset memory" );

    trouble = matset_verify_allsame( matset );
    if ( trouble ) return trouble;

    switch( matset->mats[0].type ) {
      case MUTIL_DOUBLE:
        typesize = sizeof( double );
        break;

      case MUTIL_FLOAT:
        typesize = sizeof( float );
        break;

      case MUTIL_UINT8:
        typesize = sizeof( uint8 );
        break;

      case MUTIL_UINT16:
        typesize = sizeof( uint16 );
        break;

      case MUTIL_UINT32:
        typesize = sizeof( uint32 );
        break;

      case MUTIL_SINT16:
        typesize = sizeof( sint16 );
        break;

      case MUTIL_SINT32:
        typesize = sizeof( sint32 );
        break;

      case MUTIL_DCOMPLEX:
        typesize = sizeof( dcomplex );
        break;

      default:
        MUTIL_ERROR( "This matrix type is currently unsupported" );
        return MUTIL_ERR_ILLEGAL_TYPE;
    }

    /* calculate total memory used */
    totalelem = matset->nelem * MATUNIV_NROW( &(matset->mats[0]) ) *
       MATUNIV_NCOL( &(matset->mats[0]) );

    /* free the memory */
    trouble = mutil_free( MATUNIV_DATA( &(matset->mats[0]) ),
      totalelem * typesize );
    if( trouble ) return trouble;

    /* set the contiguous flag to false */
    matset->contiguous = FALSE;
  }

  MUTIL_TRACE( "Done with matset_matrices_free()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Luca Cazzanti */
mutil_errcode matset_resize( sint32 ndim, sint32 *dims, mat_set *matset)
{
  mutil_errcode trouble;
  sint32        new_nelem;
  sint32        *new_dims;
  sint32        idx;

  MUTIL_TRACE( "Start matset_resize()" );

  trouble = matset_validate( matset );
  if( trouble ) return trouble;

  if( !dims ) {
    MUTIL_ERROR( "NULL pointer for dims array" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if( ndim < 1 ) {
    MUTIL_ERROR( "Number of dimensions must be at least 1" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  new_nelem = 1;
  for( idx = 0; idx < ndim; idx++ ) {
    new_nelem *= dims[ idx ];
  }

  if( new_nelem != matset->nelem ) {
    MUTIL_ERROR( "New dimensions do not match old number of matrices" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  trouble = mutil_malloc( (sint32) (ndim * sizeof(sint32) ),
    (void **) &new_dims );
  if( trouble ) {
    return trouble;
  }

  for( idx = 0; idx < ndim; idx++ ) {
    new_dims[ idx ] = dims[ idx ];
  }

  MUTIL_FREE_BUFFER_WARN( matset->dims, ndim * sizeof(sint32) );
  matset->dims = new_dims;
  matset->ndim = ndim;

  MUTIL_TRACE( "Done with matset_resize()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Jennifer Hodgdon */
mutil_errcode matset_validate( const mat_set *matset )
{
  sint32 dim_idx;
  sint32 nelem;

  MUTIL_TRACE( "Start matset_validate()" );

  if( !matset ) {
    MUTIL_ERROR( "NULL matrix set pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if( !matset->dims || !matset->mats ) {
    MUTIL_ERROR( "NULL matrix set member pointer" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if( matset->ndim <= 0 || matset->nelem <= 0 ) {
    MUTIL_ERROR( "Number of dimensions and elements must be positive" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  nelem = 1;
  for( dim_idx = 0; dim_idx < matset->ndim; dim_idx++ ) {
    if( matset->dims[dim_idx] <= 0 ) {
      MUTIL_ERROR( "All dimensions must be positive" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    nelem *= matset->dims[dim_idx];
  }

  if( nelem != matset->nelem ) {
    MUTIL_ERROR( "Number of elements disagrees with dimensions" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  MUTIL_TRACE( "Done with matset_validate()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Jill Goldschneider */
mutil_errcode matset_validate_matrices( const mat_set *matset )
{
  mutil_errcode   trouble;
  sint32          elem_idx;
  univ_mat       *mat;

  MUTIL_TRACE( "Start matset_validate_matrices()" );

  trouble = matset_validate( matset );
  if( trouble ) return trouble;

  for( elem_idx = 0; elem_idx < matset->nelem; elem_idx++ ) {
    mat = &(matset->mats[elem_idx]);

    trouble = matuniv_validate( mat );
    if( trouble ) return trouble;

  }

  MUTIL_TRACE( "Done with matset_validate_matrices()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Jennifer Hodgdon */
mutil_errcode matset_verify_allsame( const mat_set *matset )
{
  sint32          elem_idx;
  sint32          nrow;
  sint32          ncol;
  mutil_errcode   trouble;
  mutil_data_type type;
  univ_mat *      mat;

  MUTIL_TRACE( "Start matset_verify_allsame()" );

  trouble = matset_validate_matrices( matset );
  if( trouble ) return trouble;

  mat  = &(matset->mats[0]);
  type = mat->type;
  nrow = MATUNIV_NROW( mat );
  ncol = MATUNIV_NCOL( mat );

  for( elem_idx = 1; elem_idx < matset->nelem; elem_idx++ ) {
    mat = &(matset->mats[elem_idx]);
    if( type != mat->type ) {
      MUTIL_ERROR( "Not all matrices are same type" );
      return MUTIL_ERR_ILLEGAL_TYPE;
    }

    if( nrow != MATUNIV_NROW( mat ) ||
      ncol != MATUNIV_NCOL( mat )) {
      MUTIL_ERROR( "Not all matrices are same size" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

  }

  MUTIL_TRACE( "Done with matset_verify_allsame()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Jill Goldschneider */
mutil_errcode matset_verify_aresame( const mat_set *matset1,
    const mat_set *matset2 )
{
  sint32           dim_idx;
  sint32           elem_idx;
  sint32           nrow1;
  sint32           ncol1;
  sint32           nrow2;
  sint32           ncol2;
  mutil_errcode    trouble;
  mutil_data_type  type1;
  mutil_data_type  type2;
  univ_mat         *mat1;
  univ_mat         *mat2;

  MUTIL_TRACE( "Start matset_verify_aresame()" );

  trouble = matset_validate_matrices( matset1 );
  if ( trouble ) {
    return trouble;
  }

  trouble = matset_validate_matrices( matset2 );
  if ( trouble ) {
    return trouble;
  }

  if ( matset1->ndim != matset2->ndim ) {
    MUTIL_ERROR( "Matrix sets have different dimensions" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  for ( dim_idx = 0; dim_idx < matset1->ndim; dim_idx++ ) {
    if ( matset1->dims[dim_idx] != matset2->dims[dim_idx] ) {
        MUTIL_ERROR( "Matrix sets have different dimensions" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }
  }

  if ( matset1->nelem != matset2->nelem ) {
    MUTIL_ERROR( "Matrix sets have different number of matrices" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  for ( elem_idx = 0; elem_idx < matset1->nelem; elem_idx++ ) {

    mat1  = &(matset1->mats[elem_idx]);
    type1 = mat1->type;
    nrow1 = MATUNIV_NROW( mat1 );
    ncol1 = MATUNIV_NCOL( mat1 );

    mat2  = &(matset2->mats[elem_idx]);
    type2 = mat2->type;
    nrow2 = MATUNIV_NROW( mat2 );
    ncol2 = MATUNIV_NCOL( mat2 );

    if ( type1 != type2 ) {
      MUTIL_ERROR( "Not all matrices are same type" );
      return MUTIL_ERR_ILLEGAL_TYPE;
    }

    if ( nrow1 != nrow2 ||
         ncol1 != ncol2 ) {
      MUTIL_ERROR( "Not all matrices are same size" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

  } /* end of loop elem_idx over all matrices in set */

  MUTIL_TRACE( "Done with matset_verify_aresame()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Jill Goldschneider */
mutil_errcode matset_verify_same_dims( const mat_set *matset1,
  const mat_set *matset2 )
{
  mutil_errcode    trouble;
  sint32           dim_idx;
  sint32           nrow1;
  sint32           ncol1;
  sint32           nrow2;
  sint32           ncol2;

  MUTIL_TRACE( "Start matset_verify_same_dims()" );

  trouble = matset_verify_allsame( matset1 );
  if ( trouble ) {
    return trouble;
  }

  trouble = matset_verify_allsame( matset2 );
  if ( trouble ) {
    return trouble;
  }

  if ( matset1->ndim != matset2->ndim ) {
    MUTIL_ERROR( "Matrix sets have different dimensions" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  for ( dim_idx = 0; dim_idx < matset1->ndim; dim_idx++ ) {
    if ( matset1->dims[dim_idx] != matset2->dims[dim_idx] ) {
      MUTIL_ERROR( "Matrix sets have different dimensions" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
  }

  nrow1 = MATUNIV_NROW( &(matset1->mats[0]) );
  ncol1 = MATUNIV_NCOL( &(matset1->mats[0]) );

  nrow2 = MATUNIV_NROW( &(matset2->mats[0]) );
  ncol2 = MATUNIV_NCOL( &(matset2->mats[0]) );

  if ( nrow1 != nrow2 || ncol1 != ncol2 ) {
    MUTIL_ERROR( "Not all matrices are same size" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  MUTIL_TRACE( "Done with matset_verify_same_dims()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Paul Reinholdtsen */
mutil_errcode matset_cast(
  const mat_set  *matset,
  boolean         clip,
  void           *intrp_ptr,
  mat_set        *result )
{
  mutil_errcode  errcode;
  sint32         index;
  sint32         dim_idx;


  MUTIL_TRACE( "Start matset_cast()" );

  errcode = matset_validate( matset );
  if ( errcode ) return errcode;

  errcode = matset_validate( result );
  if ( errcode ) return errcode;

  if ( matset->ndim != result->ndim ) {
    MUTIL_ERROR( "Matrix sets have different dimensions" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  for ( dim_idx = 0; dim_idx < matset->ndim; dim_idx++ ) {
    if ( matset->dims[dim_idx] != result->dims[dim_idx] ) {
      MUTIL_ERROR( "Matrix sets have different dimensions" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
  }

  if ( matset == result ) return MUTIL_ERR_OK;

  for ( index = 0; index < matset->nelem; index++ ) {
    errcode = matuniv_cast( &( matset->mats[ index ] ), clip,
      intrp_ptr, &( result->mats[ index ] ) );
    if ( errcode ) return errcode;
  }

  MUTIL_TRACE( "Done with matset_cast()" );

  return errcode;
}


/* function documented in mat_set.h */
/* written by Jill Goldschneider */
mutil_errcode matset_assign(
  const mat_set  *matset,
  void           *intrp_ptr,
  mat_set        *result )
{
  mutil_errcode  errcode;
  sint32         index;
  sint32         dim_idx;


  MUTIL_TRACE( "Start matset_assign()" );

  errcode = matset_validate( matset );
  if ( errcode ) return errcode;

  errcode = matset_validate( result );
  if ( errcode ) return errcode;

  if ( matset == result ) return MUTIL_ERR_OK;

  if ( matset->ndim != result->ndim ) {
    MUTIL_ERROR( "Matrix sets have different dimensions" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  for ( dim_idx = 0; dim_idx < matset->ndim; dim_idx++ ) {
    if ( matset->dims[dim_idx] != result->dims[dim_idx] ) {
      MUTIL_ERROR( "Matrix sets have different dimensions" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
  }

  for ( index = 0; index < matset->nelem; index++ ) {
    errcode = matuniv_assign( &( matset->mats[ index ] ),
      intrp_ptr, &( result->mats[ index ] ) );
    if ( errcode ) return errcode;
  }

  MUTIL_TRACE( "Done with matset_assign()" );

  return errcode;
}


/* function documented in mat_set.h */
/* written by Jennifer Hodgdon */
mutil_errcode matset_tile( const mat_set *matset,
  void *intrp_ptr, univ_mat *umat )
{
  mutil_errcode trouble;
  univ_mat *    pmat;
  sint32        tile_i;
  sint32        tile_j;
  sint32        set_i;
  sint32        set_j;
  sint32        setrows;
  sint32        setcols;
  sint32        umatrows;
  sint32        umatcols;
  sint32        subrows;
  sint32        subcols;

  double        num_ops = 0;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start of matset_tile()" );

  /* sanity checks */

  trouble = matset_verify_allsame( matset );
  if( trouble ) return trouble;

  trouble = matuniv_validate( umat );
  if( trouble ) return trouble;

  pmat = &(matset->mats[0]);
  if( pmat->type != umat->type ) {
    MUTIL_ERROR( "Matrix set and output matrix must be same type" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( matset->ndim > 2 ) {
    MUTIL_ERROR( "Matrix set can be no more than two-dimensional" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* if 1D, use as a single row */

  if( matset->ndim == 1 ) {
    setcols = matset->dims[0];
    setrows = 1;
  }
  else {
    setrows = matset->dims[0];
    setcols = matset->dims[1];
  }

  umatrows = MATUNIV_NROW( umat );
  umatcols = MATUNIV_NCOL( umat );
  subrows  = MATUNIV_NROW( pmat );
  subcols  = MATUNIV_NCOL( pmat );

  if(( subrows * setrows != umatrows )
    || ( subcols * setcols != umatcols )) {
    MUTIL_ERROR( "Matrix set and output matrix have incompatible sizes" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* OK, everything checks out, do the interlacing */

#define DO_INTERLACE( MAT_TYPE ) \
    for( tile_i = 0; tile_i < subrows; tile_i++ ) { \
      for( tile_j = 0; tile_j < subcols; tile_j++ ) { \
        for( set_i = 0; set_i < setrows; set_i++ ) { \
          for( set_j = 0; set_j < setcols; set_j++ ) { \
            umat->mat. MAT_TYPE .data[ \
                 (tile_i * setrows + set_i) * umatcols + \
                 (tile_j * setcols + set_j )] = \
              matset->mats[set_i * setcols + set_j \
                 ].mat. MAT_TYPE .data[ tile_i * subcols + tile_j ]; \
          } \
          num_ops += 11 * setcols; \
          if( MUTIL_INTERRUPT( num_ops, intrp_ptr )) { \
            MUTIL_ERROR( "User interrupt" ); \
            return MUTIL_ERR_INTERRUPT; \
          } \
        } \
      } \
    }

  switch( umat->type ) {
    case MUTIL_UINT8:
      DO_INTERLACE( u8mat );
      break;
    case MUTIL_UINT16:
      DO_INTERLACE( u16mat );
      break;
    case MUTIL_UINT32:
      DO_INTERLACE( u32mat );
      break;
    case MUTIL_SINT8:
      DO_INTERLACE( s8mat );
      break;
    case MUTIL_SINT16:
      DO_INTERLACE( s16mat );
      break;
    case MUTIL_SINT32:
      DO_INTERLACE( s32mat );
      break;
    case MUTIL_DOUBLE:
      DO_INTERLACE( dblmat );
      break;
    case MUTIL_FLOAT:
      DO_INTERLACE( fltmat );
      break;
    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "Done with matset_tile()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_set.h */
/* written by Jennifer Hodgdon */
mutil_errcode matset_untile( const univ_mat *umat,
  void *intrp_ptr, mat_set *matset )
{
  mutil_errcode trouble;
  univ_mat *    pmat;
  sint32        tile_i;
  sint32        tile_j;
  sint32        set_i;
  sint32        set_j;
  sint32        setrows;
  sint32        setcols;
  sint32        umatrows;
  sint32        umatcols;
  sint32        subrows;
  sint32        subcols;

  double        num_ops = 0;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start of matset_untile()" );

  /* sanity checks */

  trouble = matset_verify_allsame( matset );
  if( trouble ) return trouble;

  trouble = matuniv_validate( umat );
  if( trouble ) return trouble;

  pmat = &(matset->mats[0]);

  if( pmat->type != umat->type ) {
    MUTIL_ERROR( "Matrix set and input matrix must be same type" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( matset->ndim > 2 ) {
    MUTIL_ERROR( "Matrix set can be no more than two-dimensional" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* if 1D, use as a single row */
  if( matset->ndim == 1 ) {
    setcols = matset->dims[0];
    setrows = 1;
  }
  else {
    setrows = matset->dims[0];
    setcols = matset->dims[1];
  }

  umatrows = MATUNIV_NROW( umat );
  umatcols = MATUNIV_NCOL( umat );
  subrows  = MATUNIV_NROW( pmat );
  subcols  = MATUNIV_NCOL( pmat );

  if(( subrows * setrows != umatrows )
    || ( subcols * setcols != umatcols )) {
    MUTIL_ERROR( "Matrix set and output matrix have incompatible sizes" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* OK, everything checks out, do the unlacing */

#define DO_UNLACE( MAT_TYPE ) \
    for( tile_i = 0; tile_i < subrows; tile_i++ ) { \
      for( tile_j = 0; tile_j < subcols; tile_j++ ) { \
        for( set_i = 0; set_i < setrows; set_i++ ) { \
          for( set_j = 0; set_j < setcols; set_j++ ) { \
            matset->mats[set_i * setcols + set_j \
                  ].mat. MAT_TYPE .data[ tile_i * subcols + tile_j ] = \
               umat->mat. MAT_TYPE .data[ \
                  (tile_i * setrows + set_i) * umatcols + \
                  (tile_j * setcols + set_j )]; \
          } \
          num_ops += 11 * setcols; \
          if( MUTIL_INTERRUPT( num_ops, intrp_ptr )) { \
            MUTIL_ERROR( "User interrupt" ); \
            return MUTIL_ERR_INTERRUPT; \
          } \
        } \
      } \
    }


  switch( umat->type ) {
    case MUTIL_UINT8:
      DO_UNLACE( u8mat );
      break;
    case MUTIL_UINT16:
      DO_UNLACE( u16mat );
      break;
    case MUTIL_UINT32:
      DO_UNLACE( u32mat );
      break;
    case MUTIL_SINT8:
      DO_UNLACE( s8mat );
      break;
    case MUTIL_SINT16:
      DO_UNLACE( s16mat );
      break;
    case MUTIL_SINT32:
      DO_UNLACE( s32mat );
      break;
    case MUTIL_DOUBLE:
      DO_UNLACE( dblmat );
      break;
    case MUTIL_FLOAT:
      DO_UNLACE( fltmat );
      break;
    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "Done with matset_untile()" );
  return MUTIL_ERR_OK;
}
