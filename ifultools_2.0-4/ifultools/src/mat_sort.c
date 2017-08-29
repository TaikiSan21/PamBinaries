
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_sort.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

/* This file contains implementations for functions in mat_sort.h,
   which are functions for sorting matrices */

#include "mat_sort.h"

#include "mat_any.h"
#include "mat_assn.h"
#include "mat_summ.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_limit.h"
#include "ut_math.h"


/* suggested value for when to do insertion sorts from Knuth */
#define INSERTION_SORT_LIMIT 9

/* headers for local functions, defined and documented at bottom of file */

static void localfn_dodbl_inssort( const double *arr, sint32 istart,
  sint32 iend, sint32 *index );

static mutil_errcode localfn_dodbl_quicksort( const double *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr );

static void localfn_doflt_inssort( const float *arr, sint32 istart,
  sint32 iend, sint32 *index );

static mutil_errcode localfn_doflt_quicksort( const float *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr );

static void localfn_dou8_inssort( const uint8 *arr, sint32 istart,
  sint32 iend, sint32 *index );

static mutil_errcode localfn_dou8_quicksort( const uint8 *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr );

static void localfn_dou16_inssort( const uint16 *arr, sint32 istart,
  sint32 iend, sint32 *index );

static mutil_errcode localfn_dou16_quicksort( const uint16 *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr );

static void localfn_dos16_inssort( const sint16 *arr, sint32 istart,
  sint32 iend, sint32 *index );

static mutil_errcode localfn_dos16_quicksort( const sint16 *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr );

static void localfn_dou32_inssort( const uint32 *arr, sint32 istart,
  sint32 iend, sint32 *index );

static mutil_errcode localfn_dou32_quicksort( const uint32 *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr );

static void localfn_dos32_inssort( const sint32 *arr, sint32 istart,
  sint32 iend, sint32 *index );

static mutil_errcode localfn_dos32_quicksort( const sint32 *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr );


static mutil_errcode localfn_sort_in_place( univ_mat *mat, void *intrp_ptr );


/******************
 Universal matrix functions
 ******************/


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode matuniv_sort_index_partial( const univ_mat *mat,
  const univ_mat *needed, void *intrp_ptr, univ_mat *index )
{
  mutil_errcode       errcode;
  const sint32_mat *  neededmat;

  MUTIL_TRACE( "Start matuniv_sort_index_partial()" );

  /* avoid lint warning */
  (void) whatssi;

  if( !mat || !index ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if( index->type != MUTIL_SINT32 ) {
    MUTIL_ERROR( "Data type for returning index must be sint32" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( needed ) {
    if( needed->type != MUTIL_SINT32 ) {
      MUTIL_ERROR( "Data type for needed indices must be sint32" );
      return MUTIL_ERR_ILLEGAL_TYPE;
    }
    neededmat = &(needed->mat.s32mat);
  }
  else {
    neededmat = NULL;
  }

  switch( mat->type ) {
    case MUTIL_DOUBLE:
      errcode = matdbl_sort_index_partial( &(mat->mat.dblmat),
        neededmat, intrp_ptr, &(index->mat.s32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_sort_index_partial( &(mat->mat.fltmat),
        neededmat, intrp_ptr, &(index->mat.s32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_sort_index_partial( &(mat->mat.u8mat),
        neededmat, intrp_ptr, &(index->mat.s32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_sort_index_partial( &(mat->mat.u16mat),
        neededmat, intrp_ptr, &(index->mat.s32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_sort_index_partial( &(mat->mat.s16mat),
        neededmat, intrp_ptr, &(index->mat.s32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_sort_index_partial( &(mat->mat.u32mat),
        neededmat, intrp_ptr, &(index->mat.s32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_sort_index_partial( &(mat->mat.s32mat),
        neededmat, intrp_ptr, &(index->mat.s32mat) );
      if(errcode) return errcode;
      break;

      /* not all types available yet */
    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "matuniv_sort_index_partial() done" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode matuniv_permute( const univ_mat *mat,
  const univ_mat *index, void *intrp_ptr, univ_mat *out )
{
  mutil_errcode       errcode;

  MUTIL_TRACE( "Start matuniv_permute()" );

  if( !mat || !index || !out ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if( index->type != MUTIL_SINT32 ) {
    MUTIL_ERROR( "Data type of index must be sint32" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if( !MATUNIV_CHECK_TYPE( mat, out )) {
    MUTIL_ERROR( "Data type of input and output must be same" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch( mat->type ) {
    case MUTIL_DOUBLE:
      errcode = matdbl_permute( &(mat->mat.dblmat),
        &(index->mat.s32mat), intrp_ptr, &(out->mat.dblmat) );
      if(errcode) return errcode;
      break;

    case MUTIL_FLOAT:
      errcode = matflt_permute( &(mat->mat.fltmat),
        &(index->mat.s32mat), intrp_ptr, &(out->mat.fltmat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT8:
      errcode = matu8_permute( &(mat->mat.u8mat),
        &(index->mat.s32mat), intrp_ptr, &(out->mat.u8mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT16:
      errcode = matu16_permute( &(mat->mat.u16mat),
        &(index->mat.s32mat), intrp_ptr, &(out->mat.u16mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_SINT16:
      errcode = mats16_permute( &(mat->mat.s16mat),
        &(index->mat.s32mat), intrp_ptr, &(out->mat.s16mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_UINT32:
      errcode = matu32_permute( &(mat->mat.u32mat),
        &(index->mat.s32mat), intrp_ptr, &(out->mat.u32mat) );
      if(errcode) return errcode;
      break;

    case MUTIL_SINT32:
      errcode = mats32_permute( &(mat->mat.s32mat),
        &(index->mat.s32mat), intrp_ptr, &(out->mat.s32mat) );
      if(errcode) return errcode;
      break;

      /* not all types available yet */
    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "matuniv_permute() done" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode matuniv_sort( const univ_mat *mat, void *intrp_ptr,
  univ_mat *out )
{
  mutil_errcode trouble;
  univ_mat    index;

  MUTIL_TRACE( "Start matuniv_sort()" );

  /* sanity checks */

  if( !mat ) {
    MUTIL_ERROR( "NULL pointer for input" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* calculate the sorting index */

  trouble = matuniv_malloc( &index, MATUNIV_NROW( mat ),
    MATUNIV_NCOL( mat ), MUTIL_SINT32 );
  if( trouble ) return trouble;

  trouble = matuniv_sort_index_partial( mat, NULL, intrp_ptr, &index );
  if( trouble ) {
    MUTIL_FREE_WARN( matuniv, &index );
    return trouble;
  }

  trouble = matuniv_permute( mat, &index, intrp_ptr, out );
  MUTIL_FREE_WARN( matuniv, &index );
  if( trouble ) return trouble;

  MUTIL_TRACE( "matuniv_sort() done" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_sort.h */
/* written by Jill Goldschneider */
mutil_errcode matuniv_table_lookup(
  const univ_mat *table,
  const univ_mat *index,
  void           *intrp_ptr,
  univ_mat       *result )
{
  mutil_errcode trouble;

  MUTIL_TRACE( "Start matuniv_table_lookup()" );

  if ( !table || !index || !result ) {
    MUTIL_ERROR( "NULL pointer for input or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if ( !MATUNIV_CHECK_TYPE( table, result ) ) {
    MUTIL_ERROR( "Data types of operand and result are inconsistent" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  switch( result->type ) {

    case MUTIL_DOUBLE:
      trouble = matdbl_table_lookup( &( table->mat.dblmat ),
        index, intrp_ptr, &( result->mat.dblmat ) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_FLOAT:
      trouble = matflt_table_lookup( &( table->mat.fltmat ),
        index, intrp_ptr, &( result->mat.fltmat ) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_UINT8:
      trouble = matu8_table_lookup( &( table->mat.u8mat ),
        index, intrp_ptr, &( result->mat.u8mat ) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_UINT16:
      trouble = matu16_table_lookup( &( table->mat.u16mat ),
        index, intrp_ptr, &( result->mat.u16mat ) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_SINT16:
      trouble = mats16_table_lookup( &( table->mat.s16mat ),
        index, intrp_ptr, &( result->mat.s16mat ) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_UINT32:
      trouble = matu32_table_lookup( &( table->mat.u32mat ),
        index, intrp_ptr, &( result->mat.u32mat ) );
      if ( trouble ) return trouble;
      break;

    case MUTIL_SINT32:
      trouble = mats32_table_lookup( &( table->mat.s32mat ),
        index, intrp_ptr, &( result->mat.s32mat ) );
      if ( trouble ) return trouble;
      break;

      /* not all types available yet */
    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "Done with matuniv_table_lookup()" );
  return MUTIL_ERR_OK;
}


/* function documented in mat_sort.h */
/* written by Andrea Borning and Paul Reinholdtsen */
mutil_errcode matuniv_unique( const univ_mat *mat, boolean sort,
  void *intrp_ptr, univ_mat *unique_vector )
{
  mutil_errcode errcode;

  MUTIL_TRACE( "Start matuniv_unique()" );

  /* sanity checks */

  if ( !mat || !unique_vector ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* assign the correct data type to unique_vector->type */
  unique_vector->type = mat->type;

  switch( mat->type ) {

    case MUTIL_UINT8:
      errcode = matu8_unique( &( mat->mat.u8mat ), sort, intrp_ptr,
        &( unique_vector->mat.u8mat ));
      break;

    case MUTIL_UINT16:
      errcode = matu16_unique( &( mat->mat.u16mat ), sort, intrp_ptr,
        &( unique_vector->mat.u16mat ));
      break;

    case MUTIL_SINT16:
      errcode = mats16_unique( &(mat->mat.s16mat), sort, intrp_ptr,
        &( unique_vector->mat.s16mat ));
      break;

    case MUTIL_UINT32:
      errcode = matu32_unique( &(mat->mat.u32mat), sort, intrp_ptr,
        &( unique_vector->mat.u32mat ));
      break;

    case MUTIL_SINT32:
      errcode = mats32_unique( &(mat->mat.s32mat), sort, intrp_ptr,
        &( unique_vector->mat.s32mat ));
      break;

    case MUTIL_FLOAT:
      errcode = matflt_unique( &(mat->mat.fltmat), sort, intrp_ptr,
        &( unique_vector->mat.fltmat ));
      break;

    case MUTIL_DOUBLE:
      errcode = matdbl_unique( &(mat->mat.dblmat), sort, intrp_ptr,
        &( unique_vector->mat.dblmat ));
      break;

    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      errcode = MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( !errcode ) {
    MUTIL_TRACE( "Done with matuniv_unique()" );
  }

  return errcode;
}


/****************
 templates for specific matrix functions and specific matrix functions
 ***************/


/** Template macro for matrix permute function.
 * Macro that expands to the body of a non-universal matrix
 * permute function, such as matdbl\_permute.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_sort.c
 * @library matrix
 * @param MAT_FN_PREFIX  Prefix for functions for this matrix type.
 * @param MAT_PTR        Matrix pointer to permute (function argument).
 * @param IDX_PTR        Matrix pointer for indices (function argument).
 * @param OUT_PTR        Matrix pointer for output (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_permute function:
 *     #TMPL_MAT_PERMUTE(matdbl, mat, index, out, intrp_ptr);#
 * @private
 */
#define TMPL_MAT_PERMUTE( MAT_FN_PREFIX, MAT_PTR, IDX_PTR, OUT_PTR, \
  INTRP_PTR ) \
  mutil_errcode trouble; \
  sint32        i;  \
\
  MUTIL_INTERRUPT_INIT( intrp_ptr ); \
  \
  MUTIL_TRACE( "Start " #MAT_FN_PREFIX "_permute()" ); \
  \
  /* sanity checks */ \
  \
  trouble = MAT_FN_PREFIX ## _validate( MAT_PTR ); \
  if( trouble ) return trouble; \
  \
  trouble = MAT_FN_PREFIX ## _validate( OUT_PTR ); \
  if( trouble ) return trouble; \
  \
  trouble = mats32_validate( IDX_PTR ); \
  if( trouble ) return trouble; \
  \
  if( !MATANY_EQUAL_DIM( IDX_PTR, OUT_PTR )) { \
    MUTIL_ERROR( "Index and output matrices must have same dimensions" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  /* apply index to get return array */  \
  \
  for( i = 0; i < (OUT_PTR)->nelem; i++ ) { \
    if( (IDX_PTR)->data[i] < 0 || \
      (IDX_PTR)->data[i] >= (MAT_PTR)->nelem ) { \
      MUTIL_ERROR( "Illegal index value" ); \
      return MUTIL_ERR_ILLEGAL_VALUE; \
    } \
    (OUT_PTR)->data[i] = (MAT_PTR)->data[ (IDX_PTR)->data[i] ]; \
  } \
  \
  if( MUTIL_INTERRUPT( 2.0 * (OUT_PTR)->nelem, INTRP_PTR )) { \
    MUTIL_ERROR( "User interrupt" ); \
    return MUTIL_ERR_INTERRUPT; \
  } \
  \
  MUTIL_TRACE( #MAT_FN_PREFIX "_permute() done" ); \
  return MUTIL_ERR_OK


/******************
 specific matrix functions
 *****************/

/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode matdbl_permute( const double_mat *mat,
const sint32_mat *index, void *intrp_ptr, double_mat *out )
{
  TMPL_MAT_PERMUTE( matdbl, mat, index, out, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode matflt_permute( const float_mat *mat,
  const sint32_mat *index, void *intrp_ptr, float_mat *out )
{
  TMPL_MAT_PERMUTE( matflt, mat, index, out, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode matu8_permute( const uint8_mat *mat,
  const sint32_mat *index, void *intrp_ptr, uint8_mat *out )
{
  TMPL_MAT_PERMUTE( matu8, mat, index, out, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode matu16_permute( const uint16_mat *mat,
  const sint32_mat *index, void *intrp_ptr, uint16_mat *out )
{
  TMPL_MAT_PERMUTE( matu16, mat, index, out, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode mats16_permute( const sint16_mat *mat,
  const sint32_mat *index, void *intrp_ptr, sint16_mat *out )
{
  TMPL_MAT_PERMUTE( mats16, mat, index, out, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon, Jill Goldschneider */
mutil_errcode matu32_permute( const uint32_mat *mat,
  const sint32_mat *index, void *intrp_ptr, uint32_mat *out )
{
  TMPL_MAT_PERMUTE( matu32, mat, index, out, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon, Jill Goldschneider */
mutil_errcode mats32_permute( const sint32_mat *mat,
  const sint32_mat *index, void *intrp_ptr, sint32_mat *out )
{
  TMPL_MAT_PERMUTE( mats32, mat, index, out, intrp_ptr );
}


/** Template macro for matrix partial index sort function.
 * Macro that expands to the body of a non-universal matrix
 * partial index sort function, such as matdbl\_sort\_index\_partial.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_sort.c
 * @library matrix
 * @param MAT_ABBREV     Abbreviation for this matrix type.
 * @param MAT_PTR        Matrix pointer to sort (function argument).
 * @param NEEDED_PTR     Matrix pointer for needed indices (function argument).
 * @param IDX_PTR        Matrix pointer for output indices (function argument).
 * @param INTRP_PTR      Pointer for interrupt handling (function argument).
 * @usage Body of the matdbl\_sort\_index\_partial function:
 *     #TMPL_MAT_SORT_INDEX_PARTIAL(dbl, mat, needed, index, intrp_ptr);#
 * @private
 */
#define TMPL_MAT_SORT_INDEX_PARTIAL( MAT_ABBREV, MAT_PTR, NEEDED_PTR, \
  IDX_PTR, INTRP_PTR ) \
  mutil_errcode trouble; \
  sint32        nelem; \
  sint32        i; \
  \
  MUTIL_TRACE( "Start mat" #MAT_ABBREV "_sort_index_partial()" ); \
  \
  /* sanity checks */ \
  \
  trouble = mat ## MAT_ABBREV ## _validate( MAT_PTR ); \
  if( trouble ) return trouble; \
  \
  trouble = mats32_validate( IDX_PTR ); \
  if( trouble ) return trouble; \
  \
  if( NEEDED_PTR ) { /* needed is NULL if we need a full sort */ \
    trouble = mats32_validate( NEEDED_PTR ); \
    if( trouble ) return trouble; \
  } \
  \
  if( !MATANY_EQUAL_DIM( MAT_PTR, IDX_PTR )) { \
    MUTIL_ERROR( "Matrix input and index output must have same dimensions" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  /* needed indices vector must be ascending and valid indices */ \
  \
  nelem = (MAT_PTR)->nelem; \
  \
  if( NEEDED_PTR ) { \
    for( i = 0; i < (NEEDED_PTR)->nelem; i++ ) { \
      if( (NEEDED_PTR)->data[i] < 0 || (NEEDED_PTR)->data[i] >= nelem ) { \
        MUTIL_ERROR( "Invalid index in needed indices" ); \
        return MUTIL_ERR_ILLEGAL_VALUE; \
      } \
    } \
  } \
  \
  /* initialize the index to 1:n */ \
  for( i = 0; i < nelem; i++ ) { \
    (IDX_PTR)->data[i] = i; \
  } \
  \
  /* perform the sort -- interrupt or success should be result */  \
  \
  trouble = localfn_do ## MAT_ABBREV ## _quicksort( (MAT_PTR)->data, \
    NEEDED_PTR, 0, nelem - 1, (IDX_PTR)->data, INTRP_PTR ); \
  if( trouble == MUTIL_ERR_INTERRUPT ) { \
    MUTIL_ERROR( "User interrupt" ); \
    return trouble; \
  } \
  \
  MUTIL_ASSERT( !trouble ); \
  \
  MUTIL_TRACE( "mat" #MAT_ABBREV "_sort_index_partial() done" ); \
  return MUTIL_ERR_OK


/******************
 specific matrix functions
 *****************/

/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode matdbl_sort_index_partial( const double_mat *mat,
const sint32_mat *needed, void *intrp_ptr, sint32_mat *index )
{
  TMPL_MAT_SORT_INDEX_PARTIAL( dbl, mat, needed, index, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode matflt_sort_index_partial( const float_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index )
{
  TMPL_MAT_SORT_INDEX_PARTIAL( flt, mat, needed, index, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode matu8_sort_index_partial( const uint8_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index )
{
  TMPL_MAT_SORT_INDEX_PARTIAL( u8, mat, needed, index, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode matu16_sort_index_partial( const uint16_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index )
{
  TMPL_MAT_SORT_INDEX_PARTIAL( u16, mat, needed, index, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon */
mutil_errcode mats16_sort_index_partial( const sint16_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index )
{
  TMPL_MAT_SORT_INDEX_PARTIAL( s16, mat, needed, index, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon, Jill Goldschneider */
mutil_errcode matu32_sort_index_partial( const uint32_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index )
{
  TMPL_MAT_SORT_INDEX_PARTIAL( u32, mat, needed, index, intrp_ptr );
}


/* function documented in mat_sort.h */
/* written by Jennifer Hodgdon, Jill Goldschneider */
mutil_errcode mats32_sort_index_partial( const sint32_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index )
{
  TMPL_MAT_SORT_INDEX_PARTIAL( s32, mat, needed, index, intrp_ptr );
}


/** Template macro for table lookup function.
 * Macro that expands the body of a non-universal image
 * table lookup function, such as TMPL\_MAT\_TBL\_LKP\_TYPE or
 * matdbl\_table\_lookup.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_sort.c
 * @library matrix
 * @param MAT_FN_PREFIX    Prefix for functions for this matrix type.
 * @param INDEX_FN_PREFIX  Prefix for functions for this index matrix type.
 * @param INDEX_TYPE       Index matrix data type.
 * @param TABLE_PTR        Matrix pointer to lookup table (function argument).
 * @param INDEX_PTR        Matrix pointer to indices (function argument).
 * @param INTRP_PTR        Interrupt pointer (function argument).
 * @param RESULT_PTR       Matrix pointer to decoded table values
 * (function argument).
 * @usage Body of the TMPL\_MAT\_TBL\_LKP\_TYPE lookup function:
 *   #TMPL_MAT_TBL_LKP(matdbl, matu32, uint32, table, index, intrp_ptr, result);#
 * @private
 */
#define TMPL_MAT_TBL_LKP( MAT_FN_PREFIX, INDEX_FN_PREFIX, INDEX_TYPE, \
   TABLE_PTR, INDEX_PTR, INTRP_PTR, RESULT_PTR ) \
  \
  INDEX_TYPE        minval; \
  INDEX_TYPE        maxval; \
  sint32            ncol; \
  sint32            icounter; \
  sint32            r_row; \
  sint32            t_row; \
  sint32            jcounter; \
  double            num_ops = 0.0; \
  double            num_ops_row; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  /* sanity checks */ \
  trouble = MAT_FN_PREFIX ## _validate( TABLE_PTR ); \
  if ( trouble ) return trouble; \
  \
  trouble = MAT_FN_PREFIX ## _validate( RESULT_PTR ); \
  if ( trouble ) return trouble; \
  \
  /* The result matrix can have more rows then the number of indices */ \
  if ( ( RESULT_PTR )->nrow < ( INDEX_PTR )->nelem ) { \
    MUTIL_ERROR( "Result matrix not capable of holding all decoded values" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  ncol = ( TABLE_PTR )->ncol; \
  if ( ncol !=  ( RESULT_PTR )->ncol ) { \
    MUTIL_ERROR( "Number of columns in the table and result must be equal" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  if ( !MATANY_IS_VEC( INDEX_PTR) ) { \
    MUTIL_ERROR( "The index must be a vector" ); \
    return MUTIL_ERR_ILLEGAL_SIZE; \
  } \
  \
  /* check range of index values */ \
  trouble = INDEX_FN_PREFIX ## _range( INDEX_PTR, INTRP_PTR, &minval, &maxval ); \
  if ( trouble ) return trouble; \
  \
  /* maxval cannot exceed nrow  */ \
  if ( ( (maxval > (INDEX_TYPE) 0) && \
    (uint32) maxval >= (uint32) MUTIL_SINT32_MAX ) || \
    ( (sint32) maxval >= ( TABLE_PTR )->nrow ) ) { \
    MUTIL_ERROR( "Index range exceeds length of table" ); \
    return MUTIL_ERR_OUT_OF_BOUNDS; \
  } \
  \
  if ( ((sint32) -minval) > 0 ) { \
    MUTIL_ERROR( "Index values cannot be negative" ); \
    return MUTIL_ERR_OUT_OF_BOUNDS; \
  } \
  \
  num_ops_row= 4.0 * ncol; \
  r_row = t_row = 0; \
  /* decode (assign a lookup row from table to result) */ \
  for ( icounter = t_row; icounter < ( INDEX_PTR )->nelem; icounter++ ) { \
    t_row = ( INDEX_PTR )->data[icounter] * ncol; \
    for ( jcounter = 0; jcounter < ncol; jcounter++ ) { \
      (RESULT_PTR)->data[r_row + jcounter] = \
        (TABLE_PTR)->data[ t_row + jcounter ]; \
    } \
    r_row += ncol; \
    num_ops += num_ops_row; \
    if ( MUTIL_INTERRUPT( num_ops, INTRP_PTR ) ) { \
      MUTIL_ERROR( "User interrupt" ); \
      return MUTIL_ERR_INTERRUPT; \
    } \
  }


/** Template macro for table lookup specific type functions function.
 * Macro that expands the body of a non-universal image
 * table lookup function, such as matdbl\_table\_lookup.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_sort.c
 * @library matrix
 * @param MAT_TYPE         Data type of table and result.
 * @usage Body of the matdbl\_table lookup function:
 *     #TMPL_MAT_TBL_LKP_TYPE(matdbl);#
 * @private
 */
#define TMPL_MAT_TBL_LKP_TYPE( MAT_TYPE ) \
  \
  mutil_errcode trouble; \
  \
  MUTIL_TRACE( "Start " #MAT_TYPE "_table_lookup()" ); \
  \
  trouble = matuniv_validate( index ); \
  if ( trouble ) return trouble; \
  \
  if ( index->type == MUTIL_UINT32 ) { \
    TMPL_MAT_TBL_LKP( MAT_TYPE, matu32, uint32, table, \
      &( index->mat.u32mat ), intrp_ptr, result ); \
  } \
  \
  else if ( index->type == MUTIL_UINT16 ) { \
    TMPL_MAT_TBL_LKP( MAT_TYPE, matu16, uint16, table, \
      &( index->mat.u16mat ), intrp_ptr, result ); \
  } \
  \
  else if ( index->type == MUTIL_UINT8 ) { \
    TMPL_MAT_TBL_LKP( MAT_TYPE, matu8, uint8, table, \
      &( index->mat.u8mat ), intrp_ptr, result ); \
  } \
  \
  else if ( index->type == MUTIL_SINT32 ) { \
    TMPL_MAT_TBL_LKP( MAT_TYPE, mats32, sint32, table, \
      &( index->mat.s32mat ), intrp_ptr, result ); \
  } \
  \
  else if ( index->type == MUTIL_SINT16 ) { \
    TMPL_MAT_TBL_LKP( MAT_TYPE, mats16, sint16, table, \
      &( index->mat.s16mat ), intrp_ptr, result ); \
  } \
  \
  else { \
    MUTIL_ERROR( "Index type must be of integer type" ); \
    return MUTIL_ERR_ILLEGAL_TYPE; \
  } \
  \
  MUTIL_TRACE( "Done with " #MAT_TYPE "_table_lookup()" ); \
  return MUTIL_ERR_OK;


/******************
 specific matrix functions
 *****************/

/* function documented in mat_sort.h */
/* written by Jill Goldschneider */
mutil_errcode matdbl_table_lookup(
  const double_mat *table,
  const univ_mat   *index,
  void             *intrp_ptr,
  double_mat       *result )
{
  TMPL_MAT_TBL_LKP_TYPE( matdbl )
}


/* function documented in mat_sort.h */
/* written by Jill Goldschneider */
mutil_errcode matflt_table_lookup(
  const float_mat *table,
  const univ_mat  *index,
  void            *intrp_ptr,
  float_mat       *result )
{
  TMPL_MAT_TBL_LKP_TYPE( matflt )
}


/* function documented in mat_sort.h */
/* written by Jill Goldschneider */
mutil_errcode matu32_table_lookup(
  const uint32_mat *table,
  const univ_mat   *index,
  void             *intrp_ptr,
  uint32_mat       *result )
{
  TMPL_MAT_TBL_LKP_TYPE( matu32 )
}


/* function documented in mat_sort.h */
/* written by Jill Goldschneider */
mutil_errcode mats32_table_lookup(
  const sint32_mat *table,
  const univ_mat   *index,
  void             *intrp_ptr,
  sint32_mat       *result )
{
  TMPL_MAT_TBL_LKP_TYPE( mats32 )
}


/* function documented in mat_sort.h */
/* written by Jill Goldschneider */
mutil_errcode matu16_table_lookup(
  const uint16_mat *table,
  const univ_mat   *index,
  void             *intrp_ptr,
  uint16_mat       *result )
{
  TMPL_MAT_TBL_LKP_TYPE( matu16 )
}


/* function documented in mat_sort.h */
/* written by Jill Goldschneider */
mutil_errcode mats16_table_lookup(
  const sint16_mat *table,
  const univ_mat   *index,
  void             *intrp_ptr,
  sint16_mat       *result )
{
  TMPL_MAT_TBL_LKP_TYPE( mats16 )
}


/* function documented in mat_sort.h */
/* written by Jill Goldschneider */
mutil_errcode matu8_table_lookup(
  const uint8_mat *table,
  const univ_mat  *index,
  void            *intrp_ptr,
  uint8_mat       *result )
{
  TMPL_MAT_TBL_LKP_TYPE( matu8 )
}


/** Template macro for constructing the unique elements of a vector.
 * Determines the unique elements of a matrix by doing the equivalent
 * of sorting the elements into ascending order and removing duplicate
 * entries.  Any duplicate entries are contiguous in the sorted data,
 * making their removal quick.  The actual algorithm uses a level of
 * indirection by computing a sorting permutation vector and removing
 * indices that point to duplicate entries.  This indirection is
 * required to keep the unique values ordered according to their first
 * appearance in the input matrix, if the sort parameter is false.  If
 * the sort parameter is true, the original data is permuted using the
 * sorted permutation vector into the output vector. If the sort
 * parameter is false, the permutation indices are sorted into
 * ascending order.  If there are multiple permutation indices that
 * point to equal indexed values, the lowest index is kept. This
 * feature ensures that if the sort parameter is false, the unique
 * values are ordered according to their first appearance in the input
 * matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_uniq.c
 * @library matrix
 * @param FN_PREFIX      Prefix for functions for this matrix type.
 *                         (e.g. dbl, s16, u32, etc.)
 * @param TYPE           Type of each element in the data of the matrix.
 * @param MAT_IN_PTR     Matrix pointer to input matrix.
 * @param SORT           Boolean to indicate whether to sort the unique
 *                         elements.
 * @param INTRP_PTR      Interrupt pointer (function argument).
 * @param VECTOR_OUT_PTR Matrix pointer to a universal matrix that has not
 *                         been previously allocated which will receive the
 *                         list of unique elements of MAT\_IN\_PTR, either
 *                         sorted or unsorted depending on the value of SORT.
 * @usage Body of the matdbl\_unique function:
 *     #TMPL_MAT_UNIQUE(dbl, double, mat, sort, intrp_ptr, unique_vector);#
 * @private
 */
#define TMPL_MAT_UNIQUE( FN_PREFIX, TYPE, MAT_IN_PTR,                           \
  SORT, INTRP_PTR, VECTOR_OUT_PTR )                                             \
                                                                                \
  TYPE *        in_data;     /* input data,                                     \
                                used to reduce indirection. */                  \
                                                                                \
  univ_mat      ordering;    /* permutation vector of the data */               \
                                                                                \
  sint32 *      order;       /* raw permutation vector,                         \
                                used to reduce indirection. */                  \
                                                                                \
  sint32        index;       /* index into the permutation vector */            \
                                                                                \
  sint32        in_nelem;    /* number of elements in input matrix,             \
                                used to reduce indirection. */                  \
                                                                                \
  sint32        n_unique;    /* number of elements in unique data */            \
                                                                                \
  mutil_errcode errcode;                                                        \
                                                                                \
  MUTIL_TRACE( "Start mat" #FN_PREFIX "_unique()" );                            \
                                                                                \
  /* sanity checks */                                                           \
                                                                                \
  errcode = mat ## FN_PREFIX ## _validate( MAT_IN_PTR );                        \
  if ( errcode ) return errcode;                                                \
                                                                                \
  if( !VECTOR_OUT_PTR ) {                                                       \
    MUTIL_ERROR( "NULL pointer for result" );                                   \
    return MUTIL_ERR_NULL_POINTER;                                              \
  }                                                                             \
                                                                                \
  in_nelem = ( MAT_IN_PTR )->nelem;                                             \
  in_data = ( MAT_IN_PTR )->data;                                               \
                                                                                \
  errcode = matuniv_malloc( &ordering,                                          \
    MAT_IN_PTR->nrow, MAT_IN_PTR->ncol, MUTIL_SINT32 );                         \
  if ( errcode ) return errcode;                                                \
                                                                                \
  /* compute the sorted permutation vector */                                   \
  errcode = mat ## FN_PREFIX ## _sort_index_partial(                            \
    MAT_IN_PTR, NULL, INTRP_PTR, &ordering.mat.s32mat );                        \
  if ( errcode ) goto cleanup;                                                  \
                                                                                \
  /* remove the permutation indices pointing to duplicate items */              \
  order = ordering.mat.s32mat.data;                                             \
  n_unique = 0;                                                                 \
  for ( index = 1; index < in_nelem; index++ ) {                                \
                                                                                \
    /* if the next value is unique,                                             \
     * add the permutation index */                                             \
    if ( in_data[ order [ n_unique ] ] != in_data[ order [ index ] ] ) {        \
      n_unique++;                                                               \
      order[ n_unique ] = order[ index ];                                       \
    }                                                                           \
                                                                                \
    /* if the value is not unique,                                              \
     * choose the permutation index with the lower valued index */              \
    else {                                                                      \
      order[ n_unique ] = MUTIL_MIN( order[ n_unique ], order[ index ] );       \
    }                                                                           \
  }                                                                             \
  n_unique++;                                                                   \
                                                                                \
  /* resize the permutation vector to the number of unique elements */          \
  errcode = matuniv_realloc( &ordering, 1, n_unique );                          \
  if ( errcode ) {                                                              \
    goto cleanup;                                                               \
  }                                                                             \
                                                                                \
  /* sorting the permutation vector, unsorts the sorted elements */             \
  if ( !sort ) {                                                                \
    errcode = localfn_sort_in_place( &ordering, INTRP_PTR );                    \
    if ( errcode ) {                                                            \
      goto cleanup;                                                             \
    }                                                                           \
  }                                                                             \
                                                                                \
  /* copy the (un)sorted elements to the output using the permutations */       \
  errcode = mat ## FN_PREFIX ## _malloc( VECTOR_OUT_PTR, 1, n_unique );         \
  if ( !errcode ) {                                                             \
    errcode = mat ## FN_PREFIX ## _permute(                                     \
      MAT_IN_PTR, &ordering.mat.s32mat, INTRP_PTR, VECTOR_OUT_PTR );            \
                                                                                \
    /* if error, keep the output vector un-allocated */                         \
    if ( errcode ) {                                                            \
      MUTIL_FREE_WARN( mat ## FN_PREFIX, VECTOR_OUT_PTR );                      \
    }                                                                           \
  }                                                                             \
                                                                                \
  /* goto here on errors that occur after ordering has been malloced */         \
cleanup:                                                                        \
                                                                                \
  MUTIL_FREE_WARN( matuniv, &ordering );                                        \
                                                                                \
  if ( !errcode ) {                                                             \
    MUTIL_TRACE( "Done with mat" #FN_PREFIX "_unique()" );                      \
  }                                                                             \
                                                                                \
  return errcode


/******************
 specific matrix functions
 *****************/


/* function documented in mat_sort.h */
/* written by Andrea Borning and Paul Reinholdtsen */
mutil_errcode matu8_unique( const uint8_mat *mat, boolean sort,
  void *intrp_ptr, uint8_mat *unique_vector )
{
  TMPL_MAT_UNIQUE( u8, uint8, mat, sort, intrp_ptr, unique_vector );
}


/* function documented in mat_sort.h */
/* written by Andrea Borning and Paul Reinholdtsen */
mutil_errcode matu16_unique( const uint16_mat *mat, boolean sort,
  void *intrp_ptr, uint16_mat *unique_vector )
{
  TMPL_MAT_UNIQUE( u16, uint16, mat, sort, intrp_ptr, unique_vector );
}


/* function documented in mat_sort.h */
/* written by Andrea Borning and Paul Reinholdtsen */
mutil_errcode matu32_unique( const uint32_mat *mat, boolean sort,
  void *intrp_ptr, uint32_mat *unique_vector )
{
  TMPL_MAT_UNIQUE( u32, uint32, mat, sort, intrp_ptr, unique_vector );
}


/* function documented in mat_sort.h */
/* written by Andrea Borning and Paul Reinholdtsen */
mutil_errcode mats16_unique( const sint16_mat *mat, boolean sort,
  void *intrp_ptr, sint16_mat *unique_vector )
{
  TMPL_MAT_UNIQUE( s16, sint16, mat, sort, intrp_ptr, unique_vector );
}


/* function documented in mat_sort.h */
/* written by Andrea Borning and Paul Reinholdtsen */
mutil_errcode mats32_unique( const sint32_mat *mat, boolean sort,
  void *intrp_ptr, sint32_mat *unique_vector )
{
  TMPL_MAT_UNIQUE( s32, sint32, mat, sort, intrp_ptr, unique_vector );
}


/* function documented in mat_sort.h */
/* written by Andrea Borning and Paul Reinholdtsen */
mutil_errcode matflt_unique( const float_mat *mat, boolean sort,
  void *intrp_ptr, float_mat *unique_vector )
{
  TMPL_MAT_UNIQUE( flt, float, mat, sort, intrp_ptr, unique_vector );
}


/* function documented in mat_sort.h */
/* written by Andrea Borning and Paul Reinholdtsen */
mutil_errcode matdbl_unique( const double_mat *mat, boolean sort,
  void *intrp_ptr, double_mat *unique_vector )
{
  TMPL_MAT_UNIQUE( dbl, double, mat, sort, intrp_ptr, unique_vector );
}


/*******************
 Local functions and Macros to generate them
 *******************/


/** Template macro for local insertion sort function.
 * Macro that expands to the body of a non-universal matrix
 * insertion sort function, such as localfn\_dodbl\_inssort.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_sort.c
 * @library matrix
 * @param DATATYPE      Data type for this function.
 * @param ARR_IN        Array to sort (function argument).
 * @param ISTART        Begin index (function argument).
 * @param IEND          End index (function argument).
 * @param ARR_IDX       Index array (function argument).
 * @usage Body of the localfn\_dodbl\_inssort function:
 *     #TMPL_LOCAL_INSSORT(double, arr, istart, iend, index);#
 * @private
 */
#define TMPL_LOCAL_INSSORT( DATATYPE, ARR_IN, ISTART, IEND, ARR_IDX ) \
  sint32   elemind; \
  sint32   cmpind; \
  sint32   elemnum; \
  DATATYPE value; \
  \
  /* go through positions starting with second */ \
  /* in each loop iteration, everything from istart to elemind - 1 \
     is in order already */ \
  \
  for( elemind = istart + 1; elemind <= iend; elemind++ ) { \
    /* store value (and index) we are inserting in correct position */ \
    \
    elemnum = (ARR_IDX)[elemind]; \
    value   = (ARR_IN)[elemnum]; \
    \
    /* go through elements (already in order) up to this position, \
       in reverse order */ \
    /* this one belongs when we get to first smaller element */ \
    /* if we don't find one, make room to insert to left */ \
    \
    for( cmpind = elemind - 1; cmpind >= istart; cmpind-- ) { \
      if( (ARR_IN)[ (ARR_IDX)[ cmpind ]] <= value ) { \
        break; \
      } \
      (ARR_IDX)[ cmpind + 1 ] = (ARR_IDX)[cmpind]; \
    } \
    \
    /* wherever we got to is where elem belongs */ \
    (ARR_IDX)[ cmpind + 1 ] = elemnum; \
  }


/** Local macro to swap values in index array.
 * This is a simple swap used by quicksort.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_sort.c
 * @library matrix
 * @param A       Index of data in index array to be swapped.
 * @param B       Index of data in index array to be swapped.
 * @param IDX     Index array of type sint32.
 * @param TMP     Temporary variable needed for swapping.
 * @usage
 *     #INDEXSWAP((ISTART) + 1, iright, ARR_IDX, iswap);#
 * @private
 */
#define INDEXSWAP(A,B,IDX,TMP) (TMP) = (IDX)[A];\
  (IDX)[A] = (IDX)[B];\
  (IDX)[B]=(TMP)


/** Local macro used by quicksort to check if data can be skipped.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_sort.c
 * @library matrix
 * @param CANSKIP       Boolean to indicate if data can be skipped or not.
 * @param MAT_NEEDED    Pointer to matrix (of type sint32)
 *                      of needed indices for partial sorting.
 * @param ISTART        Begin index (function argument).
 * @param IEND          End index (function argument).
 * @usage Used in the TMPL\_LOCAL\_QUICKSORT macro:
 *     #TMPL_LOCAL_CHECK_MAT_NEEDED(canskip, MAT_NEEDED, ISTART, IEND);#
 * @private
 */
#define TMPL_LOCAL_CHECK_MAT_NEEDED( CANSKIP, MAT_NEEDED, ISTART, IEND ) \
{ \
  sint32 idx; \
 \
  CANSKIP = TRUE; \
  for( idx = 0; idx < (MAT_NEEDED)->nelem; idx++ ) {  \
    if( (MAT_NEEDED)->data[ idx ] >= (ISTART) &&  \
        (MAT_NEEDED)->data[ idx ] <= (IEND) ) { \
      CANSKIP = FALSE; \
      break; \
    } \
  } \
}


/** Template macro for local quicksort function.
 * Macro that expands to the body of a non-universal matrix
 * quicksort function, such as localfn\_dodbl\_quicksort.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_sort.c
 * @library matrix
 * @param DATATYPE     Data type for this function.
 * @param MAT_ABBREV   Matrix abbreviation for this matrix type.
 * @param ARR_IN       Array to sort (function argument).
 * @param ISTART       Begin index (function argument).
 * @param IEND         End index (function argument).
 * @param ARR_IDX      Index array (function argument).
 * @param INTRP_PTR    Pointer for interrupt handling (function argument).
 * @param MAT_NEEDED   Pointer to matrix of needed indices for partial sorting.
 * @usage Body of the localfn\_dodbl\_quicksort function:
 *     #TMPL_LOCAL_QUICKSORT(double, dbl, arr, istart, iend, index, intrp_ptr, needed);#
 * @private
 */
#define TMPL_LOCAL_QUICKSORT( DATATYPE, MAT_ABBREV, ARR_IN, ISTART, IEND, \
  ARR_IDX, INTRP_PTR, MAT_NEEDED ) \
  boolean       canskip; \
  boolean       while_bool_outer; \
  boolean       while_bool_inner; \
  DATATYPE      part_elem; \
  mutil_errcode trouble; \
  sint32        iswap; \
  sint32        midpt; \
  sint32        ileft; \
  sint32        iright; \
  sint32        stacksize; \
  sint32        stackptr; \
  sint32       *stack; \
  \
  MUTIL_INTERRUPT_INIT( INTRP_PTR ); \
  \
  /* if small, do as insertion sort */ \
  if( (IEND) - (ISTART) + 1 <= INSERTION_SORT_LIMIT ) { \
     \
    canskip = FALSE; \
     \
    if( MAT_NEEDED ) { \
       \
      /* verify that we need at least one index in range */ \
      TMPL_LOCAL_CHECK_MAT_NEEDED( canskip, MAT_NEEDED, ISTART, IEND ); \
    } \
    if ( !canskip ) { \
      localfn_do ## MAT_ABBREV ## _inssort( ARR_IN, ISTART, \
        IEND, ARR_IDX ); \
    } \
    return MUTIL_ERR_OK; \
  } \
  \
  MUTIL_ASSERT( (IEND) - (ISTART) + 1 > INSERTION_SORT_LIMIT ); \
  \
  /* The algorithm requires an auxiliary array of storage, \
     of length 2*log_{2}(array_size), which is uses as a push-down stack for \
     keeping track of the pending subarrays. \
     Numerical Recipes in C, Second Edition, pg. 333 */ \
  /* iend/istart are sint32 values */ \
  stacksize = (sint32) ( 2 * floor( MUTIL_LOG2( IEND - ISTART ) ) ); \
  trouble = mutil_malloc( stacksize * sizeof( sint32 ), ( void** ) &stack ); \
  if ( trouble ) return trouble; \
  \
  stackptr = 0; \
  \
  /* avoid lint error of constant in conditional context */ \
  while_bool_outer = TRUE; \
  \
  while( while_bool_outer ) { \
    while ( IEND > ISTART ) { \
       \
      if( MAT_NEEDED ) { \
        \
        /* verify that we need at least one index in range */ \
        TMPL_LOCAL_CHECK_MAT_NEEDED( canskip, MAT_NEEDED, ISTART, IEND ); \
        if( canskip ) { \
          \
          /* if ( canskip ) then we don't need anything in this range \
             so exit inner loop, while ( IEND > ISTART ) */ \
          break; \
        } \
      } \
      \
      /* if small, do as insertion sort */ \
      if( (IEND) - (ISTART) + 1 <= INSERTION_SORT_LIMIT ) { \
        localfn_do ## MAT_ABBREV ## _inssort( ARR_IN, ISTART, \
          IEND, ARR_IDX ); \
        \
        /* exit inner loop, while ( IEND > ISTART ) */ \
        break; \
      } \
      \
      /* do a quicksort */ \
      \
      /* To avoid trouble with already-ordered arrays, standard thing \
         is to pick partitioning element as median of start, end, middle. \
         Put these in positions (ISTART), (ISTART) + 1, (IEND) in order. */ \
      \
      midpt = ((ISTART) + (IEND))/2; \
      INDEXSWAP( (ISTART) + 1, midpt, ARR_IDX, iswap ); \
      if((ARR_IN)[ (ARR_IDX)[ ISTART ]] >= (ARR_IN)[ (ARR_IDX)[ IEND ]] ) { \
        INDEXSWAP( ISTART, IEND, ARR_IDX, iswap ); \
      } \
      if((ARR_IN)[ (ARR_IDX)[ ISTART ]] >=  \
        (ARR_IN)[ (ARR_IDX)[ (ISTART) + 1 ]] ) { \
        INDEXSWAP( ISTART, (ISTART) + 1, ARR_IDX, iswap ); \
      } \
      \
      if((ARR_IN)[ (ARR_IDX)[ (ISTART) + 1 ]] >= \
        (ARR_IN)[ (ARR_IDX)[ IEND ]] ) { \
        INDEXSWAP( (ISTART) + 1, IEND, ARR_IDX, iswap ); \
      } \
      \
      part_elem = (ARR_IN)[ (ARR_IDX)[ (ISTART) + 1 ]]; \
       \
      /* now we move all elements <= part_elem left, and >= part_elem right \
         by swapping */ \
      \
      ileft = (ISTART) + 1; \
      iright = (IEND); \
      \
      /* avoid lint error of constant in conditional context */ \
      while_bool_inner = TRUE; \
      while( while_bool_inner ) { \
        \
        /* find next element on left that belongs right and \
           note that this will definitely terminate at (IEND) */ \
        for( ileft++; (ARR_IN)[ (ARR_IDX)[ ileft ]] < part_elem; ) { \
          ileft++; \
        } \
        \
        /* find next element on right that belongs left and \
           note that this will definitely terminate at (ISTART) */ \
        for( iright--; (ARR_IN)[ (ARR_IDX)[ iright ]] > part_elem; ) { \
          iright--; \
        } \
         \
        if( ileft >= iright ) { \
          while_bool_inner = FALSE; \
        } \
        \
        /* exchange these two elements */ \
        else { \
          INDEXSWAP( ileft, iright, ARR_IDX, iswap ); \
        } \
      } \
      \
      /* now everything has been moved -- put partitioning elem in place */ \
      INDEXSWAP( (ISTART) + 1, iright, ARR_IDX, iswap ); \
      \
      /* and sort the two halves \
         sorting the smaller first to reduce stack size */ \
      if( iright - 1 - (ISTART) > (IEND) - ileft ) { \
        \
        if ( stackptr + 2 >= stacksize ) { \
          MUTIL_ERROR( "Stack too small to complete sort" ); \
          break; \
        } \
        stack[ stackptr++ ] = ISTART; \
        stack[ stackptr++ ] = iright - 1; \
        \
        ISTART = ileft; \
        IEND = IEND; \
      } \
      else { \
         \
        MUTIL_ASSERT( stackptr + 2 < stacksize ); \
        \
        stack[ stackptr++ ] = ileft; \
        stack[ stackptr++ ] = IEND; \
         \
        ISTART = ISTART; \
        IEND = iright - 1; \
      } \
      \
    } /* end: while ( IEND > ISTART ) */ \
    if ( stackptr >= 2 ) { \
      IEND = stack[ --stackptr ]; \
      ISTART = stack[ --stackptr ]; \
    } \
    else { \
      while_bool_outer = FALSE; \
    } \
    \
    /* just add by some constant here since we don't want to slow this \
       down with another log calculation */ \
    if( MUTIL_INTERRUPT( interrupt_next_check + 10.0, INTRP_PTR )) { \
      MUTIL_ERROR( "User interrupt" ); \
      MUTIL_FREE_BUFFER_WARN( stack, stacksize ); \
      return MUTIL_ERR_INTERRUPT; \
    } \
  } /* end: while( while_bool_outer ) */ \
  \
  /* stacksize should be larger than zero, since \
     (IEND) - (ISTART) + 1 > INSERTION_SORT_LIMIT and \
     stacksize = 2 * MUTIL_LOG2( IEND - ISTART ) */ \
  MUTIL_FREE_BUFFER_WARN( stack, stacksize ); \
  \
  return trouble


/** Indexed insertion sort for double array.
 * Local function to do an insertion sort on part of an array,
 * where the result is permuted indices and the array itself is
 * not changed.
 *
 * Meant only to be called from localfn\_dodbl\_quicksort, so no
 * error checking is done on the inputs, and it always returns
 * MUTIL\_ERR\_OK.
 *
 * See Numerical Recipes chapter 8 for algorithm.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_sort.c
 * @library matrix
 * @usage  #localfn_dodbl_inssort(arr, i, j, index);#
 * @return        Standard mutils error/OK code.
 * @param arr        Flat array to sort.
 * @param istart     Beginning of index range to be sorted.
 * @param iend       End of index range to be sorted.
 * @param index      Array of indices to change into sorted order.
 * @same \begin{itemize}
 *  \item #static void localfn_dou8_inssort( const uint8 *arr, sint32 istart, sint32 iend, sint32 *index );#
 *  \item #static void localfn_dou16_inssort( const uint16 *arr, sint32 istart, sint32 iend, sint32 *index );#
 *  \item #static void localfn_dos16_inssort( const sint16 *arr, sint32 istart, sint32 iend, sint32 *index );#
 *  \item #static void localfn_dou32_inssort( const uint32 *arr, sint32 istart, sint32 iend, sint32 *index );#
 *  \item #static void localfn_dos32_inssort( const sint32 *arr, sint32 istart, sint32 iend, sint32 *index );#
 *  \item #static void localfn_doflt_inssort( const float *arr, sint32 istart, sint32 iend, sint32 *index );#
 *  \item #static void localfn_dodbl_inssort( const double *arr, sint32 istart, sint32 iend, sint32 *index );#
 *  \end{itemize}
 * @see localfn_dodbl_quicksort
 * @private
 */
static void localfn_dodbl_inssort( const double *arr, sint32 istart,
  sint32 iend, sint32 *index )
{
  TMPL_LOCAL_INSSORT( double, arr, istart, iend, index );
}

static void localfn_doflt_inssort( const float *arr, sint32 istart,
  sint32 iend, sint32 *index )
{
  TMPL_LOCAL_INSSORT( float, arr, istart, iend, index );
}

/* Function documented with the double matrix version, above */
static void localfn_dou16_inssort( const uint16 *arr, sint32 istart,
  sint32 iend, sint32 *index )
{
  TMPL_LOCAL_INSSORT( uint16, arr, istart, iend, index );
}

/* Function documented with the double matrix version, above */
static void localfn_dos16_inssort( const sint16 *arr, sint32 istart,
  sint32 iend, sint32 *index )
{
  TMPL_LOCAL_INSSORT( sint16, arr, istart, iend, index );
}

/* Function documented with the double matrix version, above */
static void localfn_dou8_inssort( const uint8 *arr, sint32 istart,
  sint32 iend, sint32 *index )
{
  TMPL_LOCAL_INSSORT( uint8, arr, istart, iend, index );
}

/* Function documented with the double matrix version, above */
static void localfn_dos32_inssort( const sint32 *arr, sint32 istart,
  sint32 iend, sint32 *index )
{
  TMPL_LOCAL_INSSORT( sint32, arr, istart, iend, index );
}

/* Function documented with the double matrix version, above */
static void localfn_dou32_inssort( const uint32 *arr, sint32 istart,
  sint32 iend, sint32 *index )
{
  TMPL_LOCAL_INSSORT( uint32, arr, istart, iend, index );
}


/** Indexed iterative partial quicksort for double array.
 * Local function to do a partial quicksort on part of an array,
 * where the result is permuted indices and the array itself is
 * not changed.
 *
 * No error checking is done on the inputs, and it always returns
 * MUTIL\_ERR\_OK unless there is an interrupt.  It checks for
 * interrupt at the beginning of each recursive call.
 *
 * See Numerical Recipes chapter 8 for algorithm.  Numerical Recipes
 * does not discuss how to do a partial sort.  To do a
 * partial sort, you simply choose not to sort sub-arrays if you
 * don't require any indices in them.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_sort.c
 * @library matrix
 * @usage  #err_code = localfn_dodbl_quicksort(arr, needed, i, j, index);#
 * @return        Standard mutils error/OK code.
 * @param arr        Flat array to sort.
 * @param needed     List of needed indices for partial sorting.
 * @param istart     Beginning of index range to be sorted.
 * @param iend       End of index range to be sorted.
 * @param index      Array of indices to change into sorted order.
 * @param intrp_ptr  Pointer for interrupt handling.
 * @same \begin{itemize}
 *  \item #static mutil_errcode localfn_dou8_quicksort( const uint8 *arr, const sint32_mat *needed, sint32 istart, sint32 iend, sint32 *index, void *intrp_ptr );#
 *  \item #static mutil_errcode localfn_dou16_quicksort( const uint16 *arr, const sint32_mat *needed, sint32 istart, sint32 iend, sint32 *index, void *intrp_ptr );#
 *  \item #static mutil_errcode localfn_dos16_quicksort( const sint16 *arr, const sint32_mat *needed, sint32 istart, sint32 iend, sint32 *index, void *intrp_ptr );#
 *  \item #static mutil_errcode localfn_dou32_quicksort( const uint32 *arr, const sint32_mat *needed, sint32 istart, sint32 iend, sint32 *index, void *intrp_ptr );#
 *  \item #static mutil_errcode localfn_dos32_quicksort( const sint32 *arr, const sint32_mat *needed, sint32 istart, sint32 iend, sint32 *index, void *intrp_ptr );#
 *  \item #static mutil_errcode localfn_doflt_quicksort( const float *arr, const sint32_mat *needed, sint32 istart, sint32 iend, sint32 *index, void *intrp_ptr );#
 *  \item #static mutil_errcode localfn_dodbl_quicksort( const double *arr, const sint32_mat *needed, sint32 istart, sint32 iend, sint32 *index, void *intrp_ptr );#
 * \end{itemize}
 * @see localfn_dodbl_inssort
 * @private
 */
static mutil_errcode localfn_dodbl_quicksort( const double *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr )
{
  TMPL_LOCAL_QUICKSORT( double, dbl, arr, istart, iend, index,
    intrp_ptr, needed );
}

/* Function documented with the double matrix version, above */
static mutil_errcode localfn_doflt_quicksort( const float *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr )
{
  TMPL_LOCAL_QUICKSORT( float, flt, arr, istart, iend, index,
    intrp_ptr, needed );
}

/* Function documented with the double matrix version, above */
static mutil_errcode localfn_dou16_quicksort( const uint16 *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr )
{
  TMPL_LOCAL_QUICKSORT( uint16, u16, arr, istart, iend, index,
    intrp_ptr, needed );
}

/* Function documented with the double matrix version, above */
static mutil_errcode localfn_dos16_quicksort( const sint16 *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr )
{
  TMPL_LOCAL_QUICKSORT( sint16, s16, arr, istart, iend, index,
    intrp_ptr, needed );
}

/* Function documented with the double matrix version, above */
static mutil_errcode localfn_dou8_quicksort( const uint8 *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr )
{
  TMPL_LOCAL_QUICKSORT( uint8, u8, arr, istart, iend, index,
    intrp_ptr, needed );
}

/* Function documented with the double matrix version, above */
static mutil_errcode localfn_dos32_quicksort( const sint32 *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr )
{
  TMPL_LOCAL_QUICKSORT( sint32, s32, arr, istart, iend, index,
    intrp_ptr, needed );
}

/* Function documented with the double matrix version, above */
static mutil_errcode localfn_dou32_quicksort( const uint32 *arr,
  const sint32_mat *needed, sint32 istart, sint32 iend,
  sint32 *index, void *intrp_ptr )
{
  TMPL_LOCAL_QUICKSORT( uint32, u32, arr, istart, iend, index,
    intrp_ptr, needed );
}


/** Sorts a universal matrix in place.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source mat\_uniq.c
 * @library matrix
 * @usage  #err_code = localfn_sort_in_place(&mymat, intrp_ptr);#
 * @return          Standard mutils error/OK code.
 * @param mat       Pointer to the matrix to sort.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @see matuniv_sort
 *
 * @private
 */
static mutil_errcode localfn_sort_in_place( univ_mat *mat, void *intrp_ptr )
{
  mutil_errcode errcode;

  univ_mat sorted_mat;

  errcode = matuniv_malloc( &sorted_mat, 1, MATUNIV_NELEM( mat ), mat->type );
  if ( errcode ) return errcode;

  errcode = matuniv_sort( mat, intrp_ptr, &sorted_mat );
  if ( errcode ) {
    MUTIL_FREE_WARN( matuniv, &sorted_mat );
    return errcode;
  }

  errcode = matuniv_assign( &sorted_mat, intrp_ptr, mat );
  MUTIL_FREE_WARN( matuniv, &sorted_mat );
  if ( errcode ) return errcode;

  return MUTIL_ERR_OK;
}
