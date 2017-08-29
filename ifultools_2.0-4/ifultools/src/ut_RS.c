
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/ut_RS.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";

/* This is a self-documenting doc++ file. */

/* This file contains definitions of functions used to convert
   R objects to and from MUTILS data structures, so that
   R functions can access MUTILS libraries.  The functions are
   declared in ut_RS.h.
*/

#include "ut_RS.h"

#ifdef __cplusplus
extern "C"{
#endif

#include "R.h"
#include "Rinternals.h"

#ifdef __cplusplus
}
#endif

#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_limit.h"
#include "ut_plat.h"

#include "mat_assn.h"
#include "mat_usca.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "mat_set.h"
#include "mat_cast.h"
#include "wav_type.h"
#include "fra_type.h"
#include "mat_type.h"
#include "sig_win.h"

/********************************/
/*                              */
/* STATIC FUNCTION DECLARATIONS */
/*                              */
/********************************/

static boolean R_extract_complex( SEXP robj, sint32 *ret_length, dcomplex **pz_data );
static boolean R_extract_integer( SEXP robj, sint32 *ret_length, sint32 **ps_data );
static boolean R_extract_logical( SEXP robj, sint32 *ret_length, sint32 **pl_data );
static boolean R_extract_numeric( SEXP robj, sint32 *ret_length, double **pd_data );

static boolean matnum_get_pieces( SEXP robj, sint32 *nrow,
  sint32 *ncol, mutil_data_type *type, void **data_ptr );
static mutil_errcode vecnum_create( sint32 nelem,
  mutil_data_type type, SEXP *robj, void **data_ptr );
static int arrnum_get_pieces( SEXP robj, sint32 *ndim,
  mutil_data_type *type, void **data_ptr, sint32 **dims );
static void matset_index( sint32 R_arr_indx, sint32 nmats,
  sint32 ndim, sint32 *dims, sint32 *mat_set_indx );
static int matnum_create( sint32 nrow, sint32 ncol,
  mutil_data_type type, SEXP *robj, void **data_ptr );

/*******************************/
/*                             */
/*      MACRO DEFINITIONS      */
/*                             */
/*******************************/

#undef ISNAN
#define ISNAN(x) _isnan(x)
#undef R_FINITE
#define R_FINITE(x) _finite(x)

/*******************************/
/*                             */
/*    FUNCTION DEFINITIONS     */
/*                             */
/*******************************/


/* wav_white_test enum mapping      */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode wav_white_test_from_R( SEXP robj,
  wav_white_test *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the enum */

  switch( asInteger(robj) ){
    case 0:
      *type = WAV_WHITE_TEST_PORTMANTEAU_I;
      break;
    case 1:
      *type = WAV_WHITE_TEST_PORTMANTEAU_II;
      break;
    case 2:
      *type = WAV_WHITE_TEST_PORTMANTEAU_III;
      break;
    case 3:
      *type = WAV_WHITE_TEST_CUMULATIVE_PERIODOGRAM;
      break;
    default:
      PROBLEM "Unsupported white noise test type" ERROR;
      break;
  }

  return MUTIL_ERR_OK;
}

/* wav_transform_peak enum mapping      */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode wav_transform_peak_from_R( SEXP robj,
  wav_transform_peak *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the enum */

  switch( asInteger(robj) ){
    case 0:
      *type = WAV_TRANSFORM_PEAK_EXTREMA;
      break;
    case 1:
      *type = WAV_TRANSFORM_PEAK_MAXIMA;
      break;
    case 2:
      *type = WAV_TRANSFORM_PEAK_MINIMA;
      break;
    default:
      PROBLEM "Unsupported wavelet transform local peak type" ERROR;
      break;
  }

  return MUTIL_ERR_OK;
}

/* sig_taper_type enum mapping      */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode sig_taper_type_from_R( SEXP robj,
  sig_taper_type *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the enum */

  switch( asInteger(robj) ){
    case 0:
      *type = SIG_TAPER_RECTANGULAR;
      break;
    case 1:
      *type = SIG_TAPER_TRIANGLE;
      break;
    case 2:
      *type = SIG_TAPER_RAISED_COSINE;
      break;
    case 3:
      *type = SIG_TAPER_HANNING;
      break;
    case 4:
      *type = SIG_TAPER_HAMMING;
      break;
    case 5:
      *type = SIG_TAPER_BLACKMAN;
      break;
    case 6:
      *type = SIG_TAPER_NUTTALL;
      break;
    case 7:
      *type = SIG_TAPER_GAUSSIAN;
      break;
    case 8:
      *type = SIG_TAPER_KAISER;
      break;
    case 9:
      *type = SIG_TAPER_CHEBYSHEV;
      break;
    case 10:
      *type = SIG_TAPER_BORN_JORDAN;
      break;
    case 11:
      *type = SIG_TAPER_SINUSOIDAL;
      break;
    case 12:
      *type = SIG_TAPER_PARZEN;
      break;
    case 13:
      *type = SIG_TAPER_PAPOULIS;
      break;
    case 14:
      *type = SIG_TAPER_DANIELL;
      break;
    default:
      PROBLEM "Unsupported taper" ERROR;
      break;
  }

  return MUTIL_ERR_OK;
}


/* fra_distance_metric enum mapping */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode fra_distance_metric_from_R( SEXP robj,
  fra_distance_metric *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the enum */

  switch( asInteger(robj) ){
    case 0:
      *type = FRA_DISTANCE_L1;
      break;
    case 1:
      *type = FRA_DISTANCE_L2;
      break;
    case 2:
      *type = FRA_DISTANCE_LINFINITY;
      break;
    default:
      PROBLEM "Unsupported distance metric" ERROR;
      break;
  }

  return MUTIL_ERR_OK;
}

/* fra_extrema_type enum mapping    */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode fra_extrema_type_from_R( SEXP robj,
  fra_extrema_type *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the enum */

  switch( asInteger(robj) ){
    case 0:
      *type = FRA_EXTREMA_MINIMA;
      break;
    case 1:
      *type = FRA_EXTREMA_MAXIMA;
      break;
    case 2:
      *type = FRA_EXTREMA_ALL;
      break;
    default:
      PROBLEM "Unsupported extrema type" ERROR;
      break;
  }

  return MUTIL_ERR_OK;
}

/* wav_filter_type enum mapping     */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode wav_filter_type_from_R( SEXP robj,
  wav_filter_type *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the white noise test */

  switch( asInteger(robj) ){
    case 0:
      *type =  WAV_FILTER_EXTREMAL_PHASE;
      break;
    case 1:
      *type =  WAV_FILTER_LEAST_ASYMMETRIC;
      break;
    case 2:
      *type =  WAV_FILTER_BEST_LOCALIZED;
      break;
    case 3:
      *type =  WAV_FILTER_COIFLET;
      break;
    case 4:
      *type =  WAV_FILTER_GAUSSIAN_I;
      break;
    case 5:
      *type =  WAV_FILTER_GAUSSIAN_II;
      break;
    case 6:
      *type =  WAV_FILTER_MORLET;
      break;
    case 7:
      *type =  WAV_FILTER_HAAR;
      break;
    default:
      PROBLEM "Unsupported filter type" ERROR;
      break;
  }

  return MUTIL_ERR_OK;
}


/* fra_surrogate enum mapping       */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode fra_surrogate_from_R( SEXP robj,
  fra_surrogate *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the white noise test */

  switch( asInteger(robj) ){
  case 0:
    *type = FRA_SURROGATE_RANDOM_PHASE;
    break;
  case 1:
    *type = FRA_SURROGATE_AAFT;
    break;
  default:
    PROBLEM "Unsupported surrogate type" ERROR;
    break;
  }

  return MUTIL_ERR_OK;
}

/* mutil_boundary_type enum mapping */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode mutil_boundary_type_from_R( SEXP robj,
  mutil_boundary_type *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the white noise test */

  switch( asInteger(robj) ){
  case 0:
    *type = MUTIL_BOUNDARY_ZERO;
    break;
  case 1:
    *type = MUTIL_BOUNDARY_PERIODIC;
    break;
  case 2:
    *type = MUTIL_BOUNDARY_REFLECT;
    break;
  case 3:
    *type = MUTIL_BOUNDARY_CONTINUE;
    break;
  default:
    PROBLEM "Unsupported boundary type" ERROR;
    break;
  }

  return MUTIL_ERR_OK;
}

/* wav_transform enum mapping       */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode wav_transform_from_R( SEXP robj,
  wav_transform *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the white noise test */

  switch( asInteger(robj) ){
    case 0:
      *type = WAV_TRANSFORM_MODWT;
      break;
    case 1:
      *type = WAV_TRANSFORM_MODWPT;
      break;
    case 2:
      *type = WAV_TRANSFORM_DWT;
      break;
    case 3:
      *type = WAV_TRANSFORM_DWPT;
      break;
  default:
    PROBLEM "Unsupported transform type" ERROR;
    break;
  }

  return MUTIL_ERR_OK;
}

/* wav_shrink_threshold enum mapping       */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode wav_shrink_threshold_from_R( SEXP robj,
  wav_shrink_threshold *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the white noise test */

  switch( asInteger(robj) ){
    case 0:
      *type = WAV_SHRINK_THRESHOLD_UNIVERSAL;
      break;
    case 1:
      *type = WAV_SHRINK_THRESHOLD_MINIMAX;
      break;
    case 2:
      *type = WAV_SHRINK_THRESHOLD_ADAPTIVE;
      break;
    default:
      PROBLEM "Unsupported wavelet shrinkage threshold function" ERROR;
      break;
  }

  return MUTIL_ERR_OK;
}

/* wav_shrink_function enum mapping       */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode wav_shrink_function_from_R( SEXP robj,
  wav_shrink_function *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the white noise test */

  switch( asInteger(robj) ){
    case 0:
      *type = WAV_SHRINK_FUNCTION_HARD;
      break;
    case 1:
      *type = WAV_SHRINK_FUNCTION_SOFT;
      break;
    case 2:
      *type = WAV_SHRINK_FUNCTION_MID;
      break;
    default:
      PROBLEM "Unsupported wavelet shrinkage function" ERROR;
      break;
  }

  return MUTIL_ERR_OK;
}

/* wav_fdp_estimator enum mapping       */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode wav_fdp_estimator_from_R( SEXP robj,
  wav_fdp_estimator *type )
{
  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* map the white noise test */

  switch( asInteger(robj) ){
    case 0:
      *type = WAV_FDP_MAXIMUM_LIKELIHOOD;
      break;
    case 1:
      *type = WAV_FDP_LEAST_SQUARES;
      break;
    default:
      PROBLEM "Unsupported wavelet shrinkage function" ERROR;
      break;
  }

  return MUTIL_ERR_OK;
}

/* Verify and obtain supported R    */
/* object class                     */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode mutil_type_from_R_class( SEXP robj,
  mutil_R_class_type *class_type )
{
  if( !robj || !class_type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if ( isNewList(robj) ){
    *class_type = MUTIL_R_LIST;
  }
  else if ( isArray(robj) ){
    *class_type = MUTIL_R_ARRAY;
  }
  else if ( isMatrix(robj) ){
    *class_type = MUTIL_R_MATRIX;
  }
  else if ( isInteger(robj) ||
            isLogical(robj) ||
            isReal(robj) ||
            isNumeric(robj) ||
            isComplex(robj) ){
    *class_type = MUTIL_R_VECTOR;
  }
  else{
    MUTIL_ERROR( "R class type not supported" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  return MUTIL_ERR_OK;
}


/* Return equivalent MUTILS data type */
/*                                    */
/* Documented in ut_RS.h               */
/* Written by William Constantine     */
mutil_errcode mutil_R_type( const SEXP robj,
  mutil_data_type *type )
{
  MUTIL_TRACE( "Start mutil_R_type()" );

  if( !robj || !type ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* see if it is a list */

  if( isNewList(robj) ){
    if( !length(robj) ){
      MUTIL_ERROR( "List has zero length, unable to find type" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    MUTIL_TRACE( "Done with mutil_R_type() -- except recursive call" );
    return mutil_R_type( VECTOR_ELT( robj, 0), type );

  }

  /* see if it is a vector, matrix, or array */

  if ( isInteger(robj) || isLogical(robj) || isReal(robj) || isComplex(robj)
    || isMatrix(robj) || isArray(robj) ){

    if( isInteger(robj) || isLogical(robj) ){
      *type = MUTIL_SINT32;
      MUTIL_TRACE( "Done with mutil_R_type()" );
      return MUTIL_ERR_OK;
    }

    if( isReal(robj) ){
      *type = MUTIL_DOUBLE;
      MUTIL_TRACE( "Done with mutil_R_type()" );
      return MUTIL_ERR_OK;
    }

    if( isComplex(robj) ){
      *type = MUTIL_DCOMPLEX;
      MUTIL_TRACE( "Done with mutil_R_type()" );
      return MUTIL_ERR_OK;
    }

  }

  /* do not know what this is */

  MUTIL_ERROR( "Unable to find numeric data type" );
  return MUTIL_ERR_ILLEGAL_TYPE;
}


/* Convert an R matrix or vector    */
/* to a universal matrix.           */
/*                                  */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode matuniv_from_R( SEXP robj,
  mutil_data_type out_type, univ_mat *umat )
{
  void            *data_ptr;
  sint32           nrow;
  sint32           ncol;
  mutil_data_type  in_type;
  mutil_errcode    trouble;
  univ_mat         tmpmat;
  univ_mat         tmpmat2;
  void            *intrp_ptr = NULL;
  double          *dbldatptr;
  sint32          *lngdatptr;
  dcomplex        *cpxdatptr;
  boolean          need_free;

  MUTIL_TRACE( "Start matuniv_from_R()" );

  /* avoid lint warning */
  (void) whatssi;

  if( !umat || !robj ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* figure out what class the input object is and read into tmpmat2 */

  need_free = FALSE;

  if( isMatrix(robj) ){


    /* it's a matrix -- extract attributes */

    if( !matnum_get_pieces( robj, &nrow, &ncol, &in_type, &data_ptr ) ){
      MUTIL_ERROR( "Unable to read data in matrix" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

    /* we have to transpose the data -- R uses column order and
       univ mat uses row order.  So reverse nrow, ncol in wrap_data  */

    trouble = matuniv_wrap_data( &tmpmat, data_ptr, ncol, nrow, in_type );
    if( trouble ){
      return trouble;
    }

    trouble = matuniv_malloc( &tmpmat2, nrow, ncol, in_type );
    if( trouble ){
      return trouble;
    }

    trouble = matuniv_transpose( &tmpmat, intrp_ptr, &tmpmat2 );
    if( trouble ){
      MUTIL_FREE_WARN( matuniv, &tmpmat2 );
      return trouble;
    }

    need_free = TRUE;

  } /* end of if( it was a matrix ) */
  /* try to read it as a numeric double vector */
  else if( R_extract_numeric( robj, &nrow, &dbldatptr ) ){


    trouble = matuniv_wrap_data( &tmpmat2, dbldatptr, nrow, 1, MUTIL_DOUBLE );
    if( trouble ){
      return trouble;
    }

  }
  else if( R_extract_integer( robj, &nrow, &lngdatptr ) ||
           R_extract_logical( robj, &nrow, &lngdatptr ) ){


    trouble = matuniv_wrap_data( &tmpmat2, lngdatptr, nrow, 1, MUTIL_SINT32 );
    if( trouble ){
      return trouble;
    }


  } /* end of if( it was integer vector ) */
  /* try to read it as a complex vector */
  else if( R_extract_complex( robj, &nrow, &cpxdatptr ) ){
    trouble = matuniv_wrap_data( &tmpmat2, cpxdatptr, nrow, 1, MUTIL_DCOMPLEX );
    if( trouble ){
      return trouble;
    }
  } /* end of if( it was a vector of complex numbers ) */
  else{
    /* nothing else supported now */
    MUTIL_ERROR( "Object unable to be interpreted as a matrix" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }


  /* convert to correct output type */

  if( tmpmat2.type == out_type && need_free ){


    /* data is correct type and we have already allocated it,
       so we can just wrap it into umat */

    /* we do this by assigning the umat header */
    *umat = tmpmat2;
  }
  else{
    /* either data is not right type or we need to allocate it */

    trouble = matuniv_malloc( umat, MATUNIV_NROW( &tmpmat2 ),
      MATUNIV_NCOL( &tmpmat2 ), out_type );
    if( trouble ){
      if( need_free ){
        MUTIL_FREE_WARN( matuniv, &tmpmat2 );
      }
      return trouble;
    }


    trouble = matuniv_cast( &tmpmat2, TRUE, intrp_ptr, umat );
    if( need_free ){
      MUTIL_FREE_WARN( matuniv, &tmpmat2 );
    }

    if( trouble ){
      MUTIL_FREE_WARN( matuniv, umat );
      return trouble;
    }

  } /* end of if/else on cast or wrap data */

  MUTIL_TRACE( "Done with matuniv_from_R()" );
  return MUTIL_ERR_OK;
}

/* Convert an R matrix or vector  */
/* to a universal matrix.         */
/*                                */
/* Documented in ut_RS.h           */
/* Written by William Constantine */
mutil_errcode matuniv_to_R( univ_mat *umat,
  mutil_R_class_type class_type, SEXP *robj )
{
  mutil_errcode trouble;

  if ( class_type == MUTIL_R_VECTOR )
    trouble = matuniv_to_R_vector( umat, robj );
  else if ( class_type == MUTIL_R_MATRIX )
    trouble = matuniv_to_R_matrix( umat, robj );
  else
    trouble = MUTIL_ERR_ILLEGAL_TYPE;

  return trouble;
}

/* Convert matuniv to an R vector. */
/*                                 */
/* Documented in ut_RS.h            */
/* Written by William Constantine  */
mutil_errcode matuniv_to_R_vector( const univ_mat *umat, SEXP *robj )
{
  mutil_errcode trouble;
  void         *out_data_ptr;
  univ_mat      tmpmat;
  univ_mat      tmpmat2;
  boolean       did_alloc;
  void         *intrp_ptr = NULL;

  MUTIL_TRACE( "Start matuniv_to_R_vector()" );

  if( !robj ){
    MUTIL_ERROR( "NULL pointer for output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  trouble = matuniv_validate( umat );
  if( trouble ){
    return trouble;
  }
  if ( !( MATUNIV_NCOL( umat ) == 1 || MATUNIV_NROW( umat ) == 1 ) ){
    MUTIL_ERROR( "Object unable to be interpreted as a vector" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* Begin by making data into either double or sint32,
     or keep it as is if already double, sint32, or complex */

  did_alloc = FALSE;

  switch( umat->type ){
    case MUTIL_DOUBLE   :
    case MUTIL_FLOAT    :
    case MUTIL_SINT32   :
    case MUTIL_DCOMPLEX :
      /* use matrix header/type directly --
         do by assignment of matrix header */
      tmpmat = *umat;
      break;

    case MUTIL_SINT16 :
    case MUTIL_SINT8  :
    case MUTIL_UINT16 :
    case MUTIL_UINT8  :
    case MUTIL_UINT32 :

      /* cast to sint32 */
      trouble = matuniv_malloc( &tmpmat,
        MATUNIV_NROW( umat ), MATUNIV_NCOL( umat ), MUTIL_SINT32 );
      if( trouble ){
        return trouble;
      }
      did_alloc = TRUE;

      trouble = matuniv_cast( umat, TRUE, intrp_ptr, &tmpmat );
      if( trouble ){
        MUTIL_FREE_WARN( matuniv, &tmpmat );
        return trouble;
      }
      break;

    default:
      MUTIL_ERROR( "Unsupported type of vector -- unable to convert to R" );
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* create output R matrix */

  trouble = vecnum_create( MATUNIV_NELEM( umat ),
    tmpmat.type, robj, &out_data_ptr );
  if( trouble ){
    MUTIL_ERROR( "Unable to allocate data for output vector" );
    if( did_alloc ){
      MUTIL_FREE_WARN( matuniv, &tmpmat );
    }
    return MUTIL_ERR_MEM_ALLOC;
  }

  /* wrap it into a univmat and assign data to it */

  trouble = matuniv_wrap_data( &tmpmat2, out_data_ptr,
    MATUNIV_NROW( umat ), MATUNIV_NCOL( umat ), tmpmat.type );
  if( trouble ){
    if( did_alloc ){
      MUTIL_FREE_WARN( matuniv, &tmpmat );
    }
    return trouble;
  }

  trouble = matuniv_assign( &tmpmat, intrp_ptr, &tmpmat2 );
  if( did_alloc ){
    MUTIL_FREE_WARN( matuniv, &tmpmat );
  }
  if( trouble ){
    return trouble;
  }

  MUTIL_TRACE( "Done with matuniv_to_R_vector()" );
  return MUTIL_ERR_OK;
}

/* Convert matuniv to an R matrix. */
/*                                 */
/* Documented in ut_RS.h            */
/* Written by William Constantine  */
mutil_errcode matuniv_to_R_matrix( const univ_mat *umat,
  SEXP *robj )
{
  mutil_errcode trouble;
  void         *out_data_ptr;
  univ_mat      tmpmat;
  univ_mat      tmpmat2;
  boolean       did_alloc;
  void         *intrp_ptr = NULL;

  MUTIL_TRACE( "Start matuniv_to_R_matrix()" );

  if( !robj ){
    MUTIL_ERROR( "NULL pointer for output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  trouble = matuniv_validate( umat );
  if( trouble ){
    return trouble;
  }

  /* Begin by making data into either double or sint32,
     or keep it as is if already double, sint32, or complex */

  did_alloc = FALSE;

  switch( umat->type ){
    case MUTIL_DOUBLE:
    case MUTIL_FLOAT:
    case MUTIL_SINT32:
    case MUTIL_DCOMPLEX:
      /* use matrix header/type directly --
         do by assignment of matrix header */
      tmpmat = *umat;
      break;

    case MUTIL_SINT16:
    case MUTIL_SINT8:
    case MUTIL_UINT16:
    case MUTIL_UINT8:
    case MUTIL_UINT32:

      /* cast to sint32 */
      trouble = matuniv_malloc( &tmpmat, MATUNIV_NROW( umat ),
        MATUNIV_NCOL( umat ), MUTIL_SINT32 );
      if( trouble ){
        return trouble;
      }
      did_alloc = TRUE;

      trouble = matuniv_cast( umat, TRUE, intrp_ptr, &tmpmat );
      if( trouble ){
        MUTIL_FREE_WARN( matuniv, &tmpmat );
        return trouble;
      }
      break;

    default:
      MUTIL_ERROR( "Unsupported type of matrix -- unable to convert to R matrix" );
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* create output R matrix */

  if( !matnum_create( MATUNIV_NROW( umat ), MATUNIV_NCOL( umat ),
    tmpmat.type, robj, &out_data_ptr ) ){
    MUTIL_ERROR( "Unable to allocate data for output matrix" );
    if( did_alloc ){
      MUTIL_FREE_WARN( matuniv, &tmpmat );
    }
    return MUTIL_ERR_MEM_ALLOC;
  }

  /* wrap it into a univmat and assign data to it */
  /* Note that the R data is transposed. */

  trouble = matuniv_wrap_data( &tmpmat2, out_data_ptr,
    MATUNIV_NCOL( umat ), MATUNIV_NROW( umat ), tmpmat.type );
  if( trouble ){
    if( did_alloc ){
      MUTIL_FREE_WARN( matuniv, &tmpmat );
    }
    return trouble;
  }

  trouble = matuniv_transpose( &tmpmat, intrp_ptr, &tmpmat2 );
  if( did_alloc ){
    MUTIL_FREE_WARN( matuniv, &tmpmat );
  }
  if( trouble ){
    return trouble;
  }

  MUTIL_TRACE( "Done with matuniv_to_R_matrix()" );
  return MUTIL_ERR_OK;
}


/* Convert an R scalar to a universal */
/* scalar of a specified type.        */
/*                                    */
/* Documented in ut_RS.h               */
/* Written by William Constantine     */
mutil_errcode scauniv_from_R( SEXP robj,
  mutil_data_type out_type, univ_scalar *usca )
{
  double          *dbldat;
  sint32          *lngdat;
  sint32          num_elem;

  MUTIL_TRACE( "Start scauniv_from_R()" );

  if( !robj || !usca ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* figure out what class the input object is and read into tmpmat2 */

  if( R_extract_numeric( robj, &num_elem, &dbldat ) ){

    /* make sure it's a 1-element vector and put it into return object */

    if( num_elem != 1 ){
      MUTIL_ERROR( "R object for scalar must have exactly 1 element" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    MUTIL_ASSERT( dbldat );
    SCAUNIV_INIT( *usca, out_type, *dbldat );
  }
  else if( R_extract_integer( robj, &num_elem, &lngdat ) ||
           R_extract_logical( robj, &num_elem, &lngdat ) ){

    /* make sure it's a 1-element vector and put it into return object */

    if( num_elem != 1 ){
      MUTIL_ERROR( "R object for scalar must have exactly 1 element" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    MUTIL_ASSERT( lngdat );

    /* create univ_scalar and convert to correct output type */
    /*LINTED: cast OK, checked range */

    SCAUNIV_INIT( *usca, out_type, *lngdat );
  }
  else{

    /* no other types are supported */

    MUTIL_ERROR( "Unable to read R object as a scalar" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "Done with scauniv_from_R()" );
  return MUTIL_ERR_OK;
}

/* Convert an R scalar to a double */
/*                                 */
/* Documented in ut_RS.h            */
/* Written by William Constantine  */
mutil_errcode double_from_R( SEXP robj,
  double *num )
{
  double    *dbldat;
  sint32    *lngdat;
  sint32    num_elem;

  MUTIL_TRACE( "Start double_from_R()" );

  if( !robj || !num ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* try to read it as an S numeric, integer, or logical */

  if( R_extract_numeric( robj, &num_elem, &dbldat ) ){

    /* make sure it's a 1-element vector and put it into return object */

    if( num_elem != 1 ){
      MUTIL_ERROR( "R object for scalar must have exactly 1 element" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    MUTIL_ASSERT( dbldat );
    *num = *dbldat;
  }
  else if( R_extract_integer( robj, &num_elem, &lngdat ) ||
           R_extract_logical( robj, &num_elem, &lngdat ) ){

    /* make sure it's a 1-element vector and put it into return object */

    if( num_elem != 1 ){
      MUTIL_ERROR( "R object for scalar must have exactly 1 element" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    MUTIL_ASSERT( lngdat );
    *num = (double) *lngdat;
  }
  else{

    /* no other types are supported */

    MUTIL_ERROR( "Unable to read R object as a scalar" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "Done with double_from_R()" );
  return MUTIL_ERR_OK;
}

/* Convert an R scalar to a dcomplex */
/*                                   */
/* Documented in ut_RS.h              */
/* Written by William Constantine    */
mutil_errcode dcomplex_from_R( SEXP robj,
  dcomplex *num )
{
  sint32        num_elem;

  MUTIL_TRACE( "Start dcomplex_from_R()" );

  if( !robj || !num ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* try to read it as an R complex  */

  if( R_extract_complex( robj, &num_elem, &num ) ){
    /* make sure it's a 1-element vector and put it into return object */
    if( num_elem != 1 ){
      MUTIL_ERROR( "R object for scalar must have exactly 1 element" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    MUTIL_ASSERT( num );
  }
  else{
    /* no other types are supported */
    MUTIL_ERROR( "Unable to read R object as a scalar" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "Done with dcomplex_from_R()" );
  return MUTIL_ERR_OK;
}

/* Convert an R scalar to a boolean  */
/*                                   */
/* Documented in ut_RS.h              */
/* Written by William Constantine    */
mutil_errcode boolean_from_R( SEXP robj,
  boolean *num )
{
  double    *dbldat;
  sint32    *lngdat;
  sint32    num_elem;

  MUTIL_TRACE( "Start boolean_from_R()" );

  if( !robj || !num ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* try to read it as an R numeric, integer, or logical */

  if( R_extract_numeric( robj, &num_elem, &dbldat ) ){
    /* make sure it's a 1-element vector and put it into return object */
    if( num_elem != 1 ){
      MUTIL_ERROR( "R object for scalar must have exactly 1 element" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    MUTIL_ASSERT( dbldat );
    if( *dbldat ){
      *num = TRUE;
    }
    else{
      *num = FALSE;
    }
  }
  else if( R_extract_integer( robj, &num_elem, &lngdat ) ||
           R_extract_logical( robj, &num_elem, &lngdat ) ){
    /* make sure it's a 1-element vector and put it into return object */
    if( num_elem != 1 ){
      MUTIL_ERROR( "R object for scalar must have exactly 1 element" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    MUTIL_ASSERT( lngdat );
    if( *lngdat ){
      *num = TRUE;
    }
    else{
      *num = FALSE;
    }
  }
  else{
    /* no other types are supported */
    MUTIL_ERROR( "Unable to read R object as a scalar" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "Done with boolean_from_R()" );
  return MUTIL_ERR_OK;
}

/* Convert an R scalar to a sint32   */
/*                                   */
/* Documented in ut_RS.h              */
/* Written by William Constantine    */
mutil_errcode sint32_from_R( SEXP robj,
  sint32 *num )
{
  double    *dbldat;
  sint32    *lngdat;
  sint32    num_elem;

  MUTIL_TRACE( "Start sint32_from_R()" );

  if( !robj || !num ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* try to read it as an R numeric, integer, or logical */

  if( R_extract_numeric( robj, &num_elem, &dbldat ) ){
    /* make sure it's a 1-element vector and put it into return object */
    if( num_elem != 1 ){
      MUTIL_ERROR( "R object for scalar must have exactly 1 element" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    MUTIL_ASSERT( dbldat );
    *num = (sint32) *dbldat;
  }
  else if( R_extract_integer( robj, &num_elem, &lngdat ) ||
           R_extract_logical( robj, &num_elem, &lngdat ) ){
    /* make sure it's a 1-element vector and put it into return object */
    if( num_elem != 1 ){
      MUTIL_ERROR( "R object for scalar must have exactly 1 element" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }
    MUTIL_ASSERT( lngdat );
    *num = *lngdat;
  }
  else{
    /* no other types are supported */
    MUTIL_ERROR( "Unable to read R object as a scalar" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "Done with sint32_from_R()" );
  return MUTIL_ERR_OK;
}

/* Convert a universal scalar to an  */
/* R scalar                          */
/* Documented in ut_RS.h              */
/* Written by William Constantine    */
mutil_errcode scauniv_to_R( univ_scalar usca, SEXP *robj )
{
  double    *dbldat;
  sint32    *lngdat;

  MUTIL_TRACE( "Start scauniv_to_R()" );

  if( !robj ){
    MUTIL_ERROR( "NULL pointer for input" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* create an R numeric vector, length 1, and initialize it */

  switch( usca.type ){
    case MUTIL_DOUBLE:

      /* cast to double */

      if( vecnum_create( (sint32) 1, usca.type, robj, (void **) &dbldat ) ){
        MUTIL_ERROR( "Unable to allocate R numeric object" );
        return MUTIL_ERR_MEM_ALLOC;
      }

      MUTIL_ASSERT( dbldat );
      *dbldat = SCAUNIV_CAST( usca, double );
      break;

    case MUTIL_SINT32:
    case MUTIL_SINT16:
    case MUTIL_SINT8:
    case MUTIL_UINT32:
    case MUTIL_UINT16:
    case MUTIL_UINT8:

      /* cast to sint32 */
      if( vecnum_create( (sint32) 1, MUTIL_SINT32, robj, (void **) &lngdat ) ){
        MUTIL_ERROR( "Unable to allocate R numeric object" );
        return MUTIL_ERR_MEM_ALLOC;
      }

      MUTIL_ASSERT( lngdat );
      *lngdat = SCAUNIV_CAST( usca, sint32 );
      break;

    default:
      MUTIL_ERROR( "Unsupported type -- unable to convert to R" );
      return MUTIL_ERR_ILLEGAL_TYPE;
  }

  MUTIL_TRACE( "Done with scauniv_to_R()" );
  return MUTIL_ERR_OK;
}

/* Convert a double to an R scalar */
/* Documented in ut_RS.h            */
/* Written by William Constantine  */
mutil_errcode double_to_R( double num, SEXP *robj )
{
  double *dbldat;

  MUTIL_TRACE( "Start double_to_R()" );

  if( !robj ){
    MUTIL_ERROR( "NULL pointer for input" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* create an R numeric vector, length 1, and initialize it */

  if( vecnum_create( (sint32) 1, MUTIL_DOUBLE, robj, (void **) &dbldat ) ){
    MUTIL_ERROR( "Unable to allocate R numeric object" );
    return MUTIL_ERR_MEM_ALLOC;
  }

  MUTIL_ASSERT( dbldat );
  *dbldat = num;

  MUTIL_TRACE( "Done with double_to_R()" );
  return MUTIL_ERR_OK;
}

/* Convert a dcomplex to an R scalar */
/* Documented in ut_RS.h              */
/* Written by William Constantine    */
mutil_errcode dcomplex_to_R( dcomplex num, SEXP *robj )
{
  dcomplex *cpxdat;

  MUTIL_TRACE( "Start dcomplex_to_R()" );

  if( !robj ){
    MUTIL_ERROR( "NULL pointer for input" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* create an R numeric vector, length 1, and initialize it */

  if( vecnum_create( (sint32) 1, MUTIL_DCOMPLEX, robj, (void **) &cpxdat ) ){
    MUTIL_ERROR( "Unable to allocate R numeric object" );
    return MUTIL_ERR_MEM_ALLOC;
  }

  MUTIL_ASSERT( cpxdat );
  *cpxdat = num;

  MUTIL_TRACE( "Done with dcomplex_to_R()" );
  return MUTIL_ERR_OK;
}

/* Convert a sint32 to an R scalar   */
/* Documented in ut_RS.h              */
/* Written by William Constantine    */
mutil_errcode sint32_to_R( sint32 num, SEXP *robj )
{
  sint32    *lngdat;

  MUTIL_TRACE( "Start sint32_to_R()" );

  if( !robj ){
    MUTIL_ERROR( "NULL pointer for input" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* create an R integer vector, length 1, and initialize it */

  if( vecnum_create( (sint32) 1, MUTIL_SINT32, robj, (void **) &lngdat ) ){
    MUTIL_ERROR( "Unable to allocate R integer object" );
    return MUTIL_ERR_MEM_ALLOC;
  }

  MUTIL_ASSERT( lngdat );
  *lngdat = num;

  MUTIL_TRACE( "Done with sint32_to_R()" );
  return MUTIL_ERR_OK;
}

/* Convert a boolean to an R scalar  */
/* Documented in ut_RS.h              */
/* Written by William Constantine    */
mutil_errcode boolean_to_R( boolean num, SEXP *robj )
{
  sint32 *lngdat;

  MUTIL_TRACE( "Start boolean_to_R()" );

  if( !robj ){
    MUTIL_ERROR( "NULL pointer for input" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* create an R logical vector, length 1, and initialize it */

  if( vecnum_create( (sint32) 1, MUTIL_SINT32, robj, (void **) &lngdat ) ){
    MUTIL_ERROR( "Unable to allocate R logical object" );
    return MUTIL_ERR_MEM_ALLOC;
  }

  MUTIL_ASSERT( lngdat );
  *lngdat = num;

  MUTIL_TRACE( "Done with boolean_to_R()" );
  return MUTIL_ERR_OK;
}

/* Convert an R list to a matset  */
/* Documented in ut_RS.h           */
/* Written by William Constantine */
mutil_errcode matset_from_R( SEXP robj, mutil_data_type type,
  mat_set *outset )
{
  void            *intrp_ptr = NULL;
  mutil_errcode    trouble;
  sint32           nmats;
  sint32           imat;

  mutil_data_type  intype;
  sint32           ndim;
  sint32           nrow;
  sint32           ncol;
  sint32          *dims;
  void            *dataptr;

  univ_mat         tmpmat;
  univ_mat         tmpmat2;
  sint32           setindx;

  MUTIL_TRACE( "Start matset_from_R()" );

  if( !outset || !robj ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* try to read it as a list */

  if( isNewList(robj) ){

    /* figure out how big the list is, and create a set */

    nmats = length(robj);
    if( nmats < 1 ){
      MUTIL_ERROR( "Input list is empty" );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    trouble = matset_malloc( outset, 1, &nmats );
    if( trouble ){
      return trouble;
    }

    /* for each element of the list, copy it to a univmat in the set */

    for( imat = 0; imat < nmats; imat++ ){


      trouble = matuniv_from_R( VECTOR_ELT( robj, imat ), type,
        &(outset->mats[imat]) );


      if( trouble ){

        /* free allocated matrices and return */

        nmats = imat;

        for ( imat = 0; imat < nmats; imat++ ){
          MUTIL_FREE_WARN( matuniv, &(outset->mats[imat] ) );
        }

        MUTIL_FREE_WARN( matset, outset );
        return trouble;
      }
    }
  } /* end of if( it is a list ) */
  else if( arrnum_get_pieces( robj, &ndim, &intype, &dataptr, &dims )){


    /* it is a multidimensional array */

    MUTIL_ASSERT( ndim > 0 );

    /* allocate a matrix set */

    nmats = 1;
    if( ndim < 3 ){
      /* only 1 matrix in set */
      trouble = matset_malloc( outset, 1, &nmats );
      if( trouble ){
        return trouble;
      }
    }
    else{
      trouble = matset_malloc( outset, ndim - 2, &(dims[2]) );
      if( trouble ){
        return trouble;
      }
      nmats = outset->nelem;
    }

    /* allocate space for the matset matrices */

    nrow = dims[0];
    if( ndim > 1 ){
      ncol = dims[1];
    }
    else{
      ncol = 1;
    }

    trouble = matset_malloc_matrices( outset, nrow, ncol, type );
    if( trouble ){
      MUTIL_FREE_WARN( matset, outset );
      return trouble;
    }

    /* create temporary space for casting and transposing */
    if( intype != type ){
      trouble = matuniv_malloc( &tmpmat2, ncol, nrow, type );
      if( trouble ){
        MUTIL_FREEALL_MATSET_WARN( outset );
        return trouble;
      }
    }

    /* Copy in the data, changing the order.
       See documentation for static function matset_index
       for more information on the order
       */

    /* go through the R array in order */

    for( imat = 0; imat < nmats; imat++ ){


        /* figure out where this goes in the matset */

      matset_index( imat, nmats, ndim - 2, &(dims[2]), &setindx );

      /* temporarily wrap matrix data */

      switch( intype ){

        case MUTIL_SINT32:
          trouble = matuniv_wrap_data( &tmpmat,
            &(( (sint32 *) dataptr )[ nrow * ncol * imat ]),
            ncol, nrow, MUTIL_SINT32 );
          break;

        case MUTIL_DOUBLE:
          trouble = matuniv_wrap_data( &tmpmat,
            &(( (double *) dataptr )[ nrow * ncol * imat ]),
            ncol, nrow, MUTIL_DOUBLE );
          break;

        case MUTIL_FLOAT:
          trouble = matuniv_wrap_data( &tmpmat,
            &(( (float *) dataptr )[ nrow * ncol * imat ]),
            ncol, nrow, MUTIL_FLOAT );
          break;

        case MUTIL_DCOMPLEX:
          trouble = matuniv_wrap_data( &tmpmat,
            &(( (dcomplex *) dataptr )[ nrow * ncol * imat ]),
            ncol, nrow, MUTIL_DCOMPLEX );
          break;

        default:
          MUTIL_ERROR( "Unsupported type of matrix -- cannot read" );
          trouble = MUTIL_ERR_ILLEGAL_TYPE;
      } /* end of switch on intype */

      if( trouble ){
        MUTIL_FREEALL_MATSET_WARN( outset );
        if( intype != type ){
          MUTIL_FREE_WARN( matuniv, &tmpmat2 );
        }
        return trouble;
      }


      /* transpose it into its resting place in the matset */

      if( intype != type ){
        trouble = matuniv_cast( &tmpmat, TRUE, intrp_ptr, &tmpmat2 );
        if( trouble ){
          MUTIL_FREEALL_MATSET_WARN( outset );
          MUTIL_FREE_WARN( matuniv, &tmpmat2 );
          return trouble;
        }

        trouble = matuniv_transpose( &tmpmat2, intrp_ptr,
          &(outset->mats[setindx] ) );
        if( trouble ){
          MUTIL_FREEALL_MATSET_WARN( outset );
          MUTIL_FREE_WARN( matuniv, &tmpmat2 );
          return trouble;
        }
      }
      else{
        trouble = matuniv_transpose( &tmpmat, intrp_ptr,
          &(outset->mats[setindx] ) );
        if( trouble ){
          MUTIL_FREEALL_MATSET_WARN( outset );
          return trouble;
        }
      }


    } /* end of loop over data */


    if( intype != type ){
      MUTIL_FREE_WARN( matuniv, &tmpmat2 );
    }

  } /* end of if( it was an array ) */

  /* last resort: try to read it as a vector */
  else{

    nmats = 1;
    trouble = matset_malloc( outset, 1, &nmats );
    if( trouble ){
      return trouble;
    }

    trouble = matuniv_from_R( robj, type, &(outset->mats[0]) );
    if( trouble ){
      MUTIL_FREE_WARN( matset, outset );
      MUTIL_ERROR( "Unable to read data from R object" );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

  }


  MUTIL_TRACE( "Done with matset_from_R()" );
  return MUTIL_ERR_OK;
}

/* Convert a matset to an R vector, */
/* matrix, array, or list.          */
/* Documented in ut_RS.h             */
/* Written by William Constantine   */
mutil_errcode matset_to_R( mat_set *set,
  mutil_R_class_type class_type, SEXP *robj )
{
  mutil_errcode trouble;

  if ( class_type == MUTIL_R_VECTOR )
    trouble = matuniv_to_R_vector( &(set->mats[0]), robj );
  else if ( class_type == MUTIL_R_MATRIX )
    trouble = matuniv_to_R_matrix( &(set->mats[0]), robj );
//  else if ( class_type == MUTIL_R_ARRAY )
//    trouble = matset_to_R_array( set, robj );
  else if ( class_type == MUTIL_R_LIST )
    trouble = matset_to_R_list( set, robj );
  else
    trouble = MUTIL_ERR_ILLEGAL_TYPE;

  return trouble;
}

/* Convert a matset to an R list. */
/* Documented in ut_RS.h           */
/* Written by William Constantine */
mutil_errcode matset_to_R_list( const mat_set *inset,
  SEXP *robj )
{
  mutil_errcode trouble;
  sint32        nmats;
  sint32        imat;
  SEXP          tmpobj;

  MUTIL_TRACE( "Start matset_to_R_list()" );

  if( !robj ){
    MUTIL_ERROR( "NULL pointer for input or output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  trouble = matset_validate( inset );
  if( trouble ){
    return trouble;
  }

  nmats = inset->nelem;

  /* create a list R object */

  PROTECT(*robj = allocVector(VECSXP, nmats));

  MUTIL_ASSERT( robj && *robj );
  MUTIL_ASSERT( nmats < MUTIL_SINT32_MAX );

  /* create matrix R objects and put them into the list */

  for( imat = 0; imat < nmats; imat++ ){
    trouble = matuniv_to_R_matrix( &(inset->mats[imat]),
      &tmpobj );
    if( trouble ){
      return trouble;
    }
    /*LINTED: cast OK, checked range */
    SET_VECTOR_ELT( *robj, (int) imat, tmpobj );
  }

  UNPROTECT(1);

  MUTIL_TRACE( "Done with matset_to_R_list()" );

  return MUTIL_ERR_OK;
}

/* Extract boundary type from R object. */
/* Documented in ut_RS.h                 */
/* Written by William Constantine       */
mutil_errcode mutil_boundary_from_R( SEXP robj,
    mutil_boundary_type *boundary )
{
  mutil_errcode err;
  sint32        s32_boundary;

  /* error checking done by sint32_from_R */

  MUTIL_TRACE( "Start mutil_boundary_from_R()" );

  err = sint32_from_R( robj, &s32_boundary );
  if( err ){
      return err;
  }
  switch( s32_boundary ){
    case 0:
      *boundary = MUTIL_BOUNDARY_ZERO;
      break;
    case 1:
      *boundary = MUTIL_BOUNDARY_PERIODIC;
      break;
    case 2:
      *boundary = MUTIL_BOUNDARY_REFLECT;
      break;
    case 3:
      *boundary = MUTIL_BOUNDARY_CONTINUE;
      break;
    default:
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

  MUTIL_TRACE( "Done with mutil_boundary_from_R()" );
  return MUTIL_ERR_OK;

}

/* Determine if an R object contains NULL. */
/* Documented in ut_RS.h                   */
/* Written by William Constantine          */
mutil_errcode null_object_from_R( SEXP robj, boolean *is_null )
{

  MUTIL_TRACE( "Start null_object_from_R()" );

  if( !robj ){
    MUTIL_ERROR( "NULL pointer for input" );
    return MUTIL_ERR_NULL_POINTER;
  }

  *is_null = FALSE;
  if ( isNull(robj) ){
    *is_null = TRUE;
  }

  MUTIL_TRACE( "Done with null_object_from_R()" );
  return MUTIL_ERR_OK;
}

/*******************************/
/*                             */
/* STATIC FUNCTION DEFINITIONS */
/*                             */
/*******************************/

static boolean R_extract_complex( SEXP robj, sint32 *ret_length, dcomplex **pz_data )
{

  /* test for type dcomplex */

  if ( !isComplex(robj) ){
    return (boolean) FALSE;
  }

  /* obtain length and data pointer */

  *ret_length = (sint32) length(robj);
  *pz_data    = (dcomplex *) COMPLEX(robj);

  return (boolean) TRUE;
}

static boolean R_extract_integer( SEXP robj, sint32 *ret_length, sint32 **ps_data )
{

  /* test for type integer */

  if ( !isInteger(robj) ){
    return (boolean) FALSE;
  }

  /* obtain length and data pointer */

  *ret_length = (sint32) length(robj);
  *ps_data    = (sint32 *) INTEGER(robj);

  return (boolean) TRUE;
}

static boolean R_extract_logical( SEXP robj, sint32 *ret_length, sint32 **pl_data )
{

  /* test for type logical */

  if ( !isLogical(robj) ){
    return (boolean) FALSE;
  }

  /* obtain length and data pointer */

  *ret_length = (sint32) length(robj);
  *pl_data    = (sint32 *) LOGICAL(robj);

  return (boolean) TRUE;
}

static boolean R_extract_numeric( SEXP robj, sint32 *ret_length, double **pd_data )
{

  /* test for type double */

  if ( !isReal(robj) ){
    return (boolean) FALSE;
  }

  /* obtain length and data pointer */

  *ret_length = (sint32) length(robj);
  *pd_data    = (double *) REAL(robj);

  return (boolean) TRUE;
}

/** Extract data characteristics from an R matrix.
 * From an R matrix object, extract the number of rows
 * and columns, the data type, and a pointer to the flat array
 * containing the data.  Fails if the input is not a numeric
 * or integer or complex R matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source ut\_R.c
 * @library wrap
 * @usage  #ok = matnum_get_pieces(robj, &nrow, &ncol, &type, &data);#
 * @return   TRUE/FALSE for success/failure.
 * @param robj     The input R object containing a matrix.
 * @param nrow     Pointer for returning number of rows.
 * @param ncol     Pointer for returning number of columns.
 * @param type     Pointer for returning type of data
 *                 (numeric or integer or complex).
 * @param data_ptr Pointer to the data in the object.
 */
static boolean matnum_get_pieces( SEXP robj, sint32 *nrow,
  sint32 *ncol, mutil_data_type *type, void **data_ptr )
{
  sint32     lngth;
  SEXP       dim;

  if( !robj || !nrow || !ncol || !type || !data_ptr ){
    return FALSE;
  }

  /* make sure it's a matrix */

  if ( !isMatrix(robj) ){
    return FALSE;
  }

  /* get the dimensions */

  PROTECT( dim = getAttrib( robj, install("dim") ) );

  *nrow = INTEGER( dim )[0];
  *ncol = INTEGER( dim )[1];

  UNPROTECT(1);

  /* try to extract as numeric - double */

  if( R_extract_numeric( robj, &lngth, (double **) data_ptr ) ){
    *type = MUTIL_DOUBLE;
    MUTIL_ASSERT( *data_ptr );
    if( lngth != *nrow * *ncol ){
      return FALSE;
    }
    return TRUE;
  }

  /* try to extract as integer */

  if( R_extract_integer( robj, &lngth, (sint32 **) data_ptr ) ){
    *type = MUTIL_SINT32;
    MUTIL_ASSERT( *data_ptr );
    if( lngth != *nrow * *ncol ){
      return FALSE;
    }
    return TRUE;
  }

   /* try to extract as complex */

  if( R_extract_complex( robj, &lngth, (dcomplex **) data_ptr ) ){
    *type = MUTIL_DCOMPLEX;
    MUTIL_ASSERT( *data_ptr );
    if( lngth != *nrow * *ncol ){
      return FALSE;
    }
    return TRUE;
  }

  return FALSE;
}

/** Extract data from an R array.
 * From an R array (or matrix) object,
 * extract the number of dimensions,
 * the data type, and pointers to the flat array
 * containing the data and the flat array containing the
 * dimensions.  Fails if the input is not a numeric
 * or integer R array.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source ut\_R.c
 * @library wrap
 * @usage  #ok = arrnum_get_pieces(robj, &ndim, &type, &data, &dims);#
 * @return   TRUE/FALSE for success/failure.
 * @param robj     The input R object containing an array.
 * @param nrow     Pointer for returning number of dimensions.
 * @param type     Pointer for returning type of data (numeric or integer).
 * @param data_ptr Pointer to the data in the object.
 * @param dims     Pointer to the dimensions vector of the object.
 */
static int arrnum_get_pieces( SEXP robj, sint32 *ndim,
  mutil_data_type *type, void **data_ptr, sint32 **dims )
{
  sint32     elem;
  sint32     lngth;
  sint32     dimprod;
  SEXP       rdim;

  if( !robj || !ndim || !type || !data_ptr || !dims ){
    return FALSE;
  }

  /* make sure it's an array or matrix */

  if( !isArray(robj) && !isMatrix(robj)){
    return FALSE;
  }

  /* get the dimensions and check validity */

  /* get the dimensions */

  PROTECT( rdim = getAttrib( robj, install("dim") ) );

  if( !rdim ||
    !R_extract_integer( rdim, ndim, dims ) ||
    ( *ndim < 1 ) ){
    return FALSE;
  }
  MUTIL_ASSERT( *dims );

  dimprod = 1;
  for( elem = 0; elem < *ndim; elem++ ){
    if( (*dims)[elem] <= 0 ){
      return FALSE;
    }
    dimprod *= (*dims)[elem];
  }


  /* try to extract data as numeric - double */

  if( R_extract_numeric( robj, &lngth, (double **) data_ptr ) ){
    *type = MUTIL_DOUBLE;
    MUTIL_ASSERT( *data_ptr );
    if( lngth != dimprod ){
      return FALSE;
    }
    return TRUE;
  }

  /* try to extract data as integer */

  if( R_extract_integer( robj, &lngth, (sint32 **) data_ptr ) ){
    *type = MUTIL_SINT32;
    MUTIL_ASSERT( *data_ptr );
    if( lngth != dimprod ){
      return FALSE;
    }
    return TRUE;
  }

   /* try to extract data as complex */

  if( R_extract_complex( robj, &lngth, (dcomplex **) data_ptr ) ){
    *type = MUTIL_DCOMPLEX;
    MUTIL_ASSERT( *data_ptr );
    if( lngth != dimprod ){
      return FALSE;
    }
    return TRUE;
  }

  UNPROTECT(1);
  return FALSE;
}

/** Create an R vector.
 * Create an R vector object with the given length
 * and the given data type. Return the object and
 * the pointer to the flat array containing the data.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source ut\_R.c
 * @library wrap
 * @usage  #err = vecnum_create(nelem, type, &robj, &data);#
 * @return   TRUE/FALSE for success/failure.
 * @param nelem    Number of elements in vector.
 * @param type     Type of data (numeric or integer).
 * @param robj     Pointer for returning the created R matrix object.
 * @param data_ptr Pointer to the data in the object.
 */
static mutil_errcode vecnum_create( sint32 nelem,
  mutil_data_type type, SEXP *robj, void **data_ptr )
{

  /* create the data slot and extract its pointer */

  if( type == MUTIL_DOUBLE ){


    PROTECT( *robj = allocVector(REALSXP, nelem) );
    *data_ptr = REAL(*robj);

    MUTIL_ASSERT( *data_ptr );
  }

  else if( type == MUTIL_DCOMPLEX ){

    PROTECT( *robj = allocVector(CPLXSXP, nelem) );
    *data_ptr = (dcomplex *) (COMPLEX(*robj));
    MUTIL_ASSERT( *data_ptr );
  }

  else if( type == MUTIL_SINT32 ){

    PROTECT( *robj = allocVector(INTSXP, nelem) );
    *data_ptr = (sint32 *) (INTEGER(*robj));
    MUTIL_ASSERT( *data_ptr );
  }

  /* unsupported type */
  else{
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  UNPROTECT(1);
  return MUTIL_ERR_OK;
}

/** Calculate matrix set index corresponding to R array index.
 * For a given sequential index into an R array, figure out
 * where in the sequential index of the corresponding MUTILS
 * matrix set the matrix belongs.
 *
 * R arrays have a flat data segment, and matrix sets are stored
 * with each matrix in the set in its own data segment.  The meaning of
 * the order of the dimensions is the same for R and matrix sets --
 * e.g. when we talk of the [k,l,m,n] element of the matrix, we mean
 * the kth row, lth column, and the (m,n) member of the set.
 * In R arrays, the fastest-changing
 * dimension is the row, then the column, then the first set dimension,
 * the second set dimension, and so on.  In the matrix set, the
 * last dimension varies the fastest within the set, then the second
 * to last, and so on; and within each matrix, the column dimension
 * varies the fastest and then the row.  So, in both cases, the matrices
 * within the set have contiguous data, but the data are transposed.
 *
 * For example, a 2x3x4x5 R array means that it is a
 * 4x5 set of 2x3 matrices. The (i,j) matrix comes at position
 * (i + 4 * j) in the R array,
 * position (5 * i + j) in the mat set.  So, for i=1 and j=3,
 * this function would take input index 13 for the R array,
 * and calculate output index 8 in the matrix set.
 * For a 2x3x4x5x6 array, the (i,j,k) matrix is at position
 * i + 4 * j + 20 * k in R, and 6 * i + 30 * j + k in mat set.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source ut\_R.c
 * @library wrap
 * @usage  #matset_indx(spl_indx, nmats, ndim, &dims, &matset_indx);#
 * @param R_arr_indx     Sequential index of R array (starting at 0).
 * @param nmats            Total number of matrices in set.
 * @param ndim             Number of dimensions of matrix set.
 * @param dims             Dimensions of the set.
 * @param mat_set_indx     Pointer for returning sequential index into matrix
 *                         set (starting at 0).
 */
static void matset_index( sint32 R_arr_indx,
  sint32 nmats, sint32 ndim, sint32 *dims, sint32 *mat_set_indx )
{
  sint32  dimprod;
  sint32  remaining;
  sint32  idim;
  sint32  itmp;

  /* we will add up the matrix set index, dimension by dimension */

  *mat_set_indx  = 0;
  dimprod        = nmats;
  remaining      = R_arr_indx;

  for( idim = 0; idim < ndim; idim++ ){

    /* figure out the actual R index in this dimension */

    itmp      = remaining % dims[ idim ];
    remaining = ( remaining - itmp ) / dims[ idim ];

    /* add that into the matset dimension index */

    dimprod /= dims[ idim ];
    *mat_set_indx += dimprod * itmp;
  }
}

/** Create an R matrix.
 * Create an R matrix object with the given number of rows
 * and columns, and the given data type.  Return the object and
 * the pointer to the flat array containing the data.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source ut\_R.c
 * @library wrap
 * @usage  #ok = matnum_create(nrow, ncol, type, &robj, &data);#
 * @return TRUE/FALSE for success/failure.
 * @param nrow     Number of rows.
 * @param ncol     Number of columns.
 * @param type     Type of data (numeric or integer).
 * @param robj     Pointer for returning the created R matrix object.
 * @param data_ptr Pointer to the data in the object.
 */
static int matnum_create( sint32 nrow, sint32 ncol,
  mutil_data_type type, SEXP *robj, void **data_ptr )
{
  if( nrow <= 0 || ncol <= 0 || !robj || !data_ptr ){
    return FALSE;
  }

  /* create the data slot and extract its pointer */

  if( type == MUTIL_DOUBLE ){

    PROTECT( *robj = allocMatrix( REALSXP, (int) nrow, (int) ncol ) );
    *data_ptr = REAL(*robj);

  }
  else if( type == MUTIL_DCOMPLEX ){
    PROTECT( *robj = allocMatrix( CPLXSXP, (int) nrow, (int) ncol ) );
    *data_ptr = (dcomplex *) (COMPLEX(*robj));
  }
  else if( type == MUTIL_SINT32 ){
    PROTECT( *robj = allocMatrix( INTSXP, (int) nrow, (int) ncol ) );
    *data_ptr = (sint32 *) (INTEGER(*robj));
  }
  else {
    /* unsupported type */
    return FALSE;
  }

  MUTIL_ASSERT( *data_ptr );

  UNPROTECT( 1 );

  return TRUE;
}
