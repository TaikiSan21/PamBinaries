
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_RS.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_RS_H
#define IN_UT_RS_H

/* This file contains definitions of functions used to convert
   R and S-PLUS objects to and from MUTILS data structures, so that
   R and S-PLUS functions can access MUTILS libraries.  The functions are
   defined in ut_RS.c.
*/

#ifdef __cplusplus
extern "C" {
#endif

// COMMENTED OUT
//
// /* S.h has a typedef for boolean -- ours is different.  If you
//   want to use the S.h version, call it LOCAL_S5_boolean in
//   any functions that include this header file */
//
// /* Also, make sure S.h is included within the __cplusplus ifdef,
// and that it is the first include */
//
// #define boolean LOCAL_S5_boolean
// #define NO_NEWIO
// #define USE_NEWIO
// #include "S.h"
// #undef boolean

/* MESSAGE, PROBLEM, ERROR, RECOVER, WARNING, LOCAL_EVALUATOR, NULL_ENTRY, and WARN
are defined in R_ext/RS.h which is included in R.h. These are only accessible if
STRICT_R_HEADERS is not defined, so make sure that we undefine in the case where
it is defined */

#ifdef STRICT_R_HEADERS
#undef STRICT_R_HEADERS
#endif

/* the compiler won't let us cast
   to/from Rcomplex/dcomplex,
   although the two structures are identical.
   So we use this trick. */

//#define dcomplex tmp_dcomplex
//#include "ut_type.h"
//#undef dcomplex
//#define dcomplex Rcomplex

#include "R.h"
#include "Rinternals.h"
#include "mat_type.h"
#include "ut_err.h"
#include "wav_type.h"
#include "fra_type.h"
#include "mat_type.h"
#include "sig_type.h"

/** R PROBLEM/ERROR messaging for the err message returned
 * by a call to any of the MUTILS libraries.
 * This macro looks up the appropriate error string message
 * and passes the message to R macros PROBLEM and ERROR.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.h
 * @library wrap\_RS
 * @usage  #if(err) R_PROBLEM_ERROR_MSG(err);#
 * @see mutil_err_string
 * @see _mutil_errcode
 */
#define R_PROBLEM_ERROR_MSG( TROUBLE ) \
{ \
  const char    *msgstr; \
  mutil_errcode  othererr; \
  othererr = mutil_err_string(TROUBLE, &msgstr); \
  if( othererr || !msgstr) \
    PROBLEM "error messaging seems to have failed" ERROR; \
  else \
    PROBLEM msgstr ERROR; \
}

/** R PROBLEM/WARN messaging for the err message returned
 * by a call to any of the MUTILS libraries.
 * This macro looks up the appropriate error string message
 * and passes the message to R macros PROBLEM and WARN.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.h
 * @library wrap\_RS
 * @usage  #if(err) R_PROBLEM_WARN_MSG(err);#
 * @see mutil_err_string
 * @see _mutil_errcode
 */
#define R_PROBLEM_WARN_MSG( TROUBLE ) \
{ \
  const char    *msgstr; \
  mutil_errcode  othererr; \
  othererr = mutil_err_string(TROUBLE, &msgstr); \
  if( othererr || !msgstr) \
    PROBLEM "warning messaging seems to have failed" ERROR; \
  else \
    PROBLEM msgstr WARN; \
}

/** sig_taper_type enum mapping
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = sig_taper_type(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode sig_taper_type_from_R( SEXP robj,
  sig_taper_type *type );

/** wav_white_test enum mapping.
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = wav_white_test_from_R(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode wav_white_test_from_R( SEXP robj,
  wav_white_test *type );

/** wav_transform_peak enum mapping.
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = wav_transform_peak_from_R(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode wav_transform_peak_from_R( SEXP robj,
  wav_transform_peak *type );

/** fra_distance_metric enum mapping.
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = fra_distance_metric(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode fra_distance_metric_from_R( SEXP robj,
  fra_distance_metric *type );

/** fra_extrema_type enum mapping.
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = fra_extrema_type(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode fra_extrema_type_from_R( SEXP robj,
  fra_extrema_type *type );

/** wav_filter_type enum mapping.
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = wav_filter_type_from_R(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode wav_filter_type_from_R( SEXP robj,
  wav_filter_type *type );

/** fra_surrogate enum mapping.
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = fra_surrogate_from_R(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode fra_surrogate_from_R( SEXP robj,
  fra_surrogate *type );

/** mutil_boundary_type enum mapping.
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = mutil_boundary_type_from_R(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_boundary_type_from_R( SEXP robj,
  mutil_boundary_type *type );

/** wav_transform_from_R enum mapping.
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = wav_transform_from_R(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode wav_transform_from_R( SEXP robj,
  wav_transform *type );

/** wav_shrink_threshold enum mapping.
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = wav_shrink_threshold_from_R(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode wav_shrink_threshold_from_R( SEXP robj,
  wav_shrink_threshold *type );

/** wav_shrink_function enum mapping.
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = wav_shrink_function_from_R(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode wav_shrink_function_from_R( SEXP robj,
  wav_shrink_function *type );

/** wav_fdp_estimator enum mapping.
 * Maps an integer from R to the corresponding enum type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = wav_fdp_estimator_from_R(robj, &type);#
 * @return      Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The corresponding enum type.
 */
MUTIL_WRAPEXPORT mutil_errcode wav_fdp_estimator_from_R( SEXP robj,
  wav_fdp_estimator *type );

/** Find numeric type of R data.
 * Find the underlying type of an R object.  First check to
 * see if it is an integer, logical, or numeric vector, matrix, or array.
 * Then, if it is a list, look at the first element of the list
 * and make the same check, calling this function recursively.
 * The function will return an error if the data is not numbers.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = mutil_R_type(robj, &type);#
 * @return        Standard mutils error/OK code.
 * @param robj  The R object to check.
 * @param type  The type of the object (generally either sint32 or double).
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_R_type( const SEXP robj,
  mutil_data_type *type );

/** Enum of R types supported by MUTILS.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.h
 * @library wrap\_RS
 * @same #typedef enum _mutil_R_class_type mutil_R_class_type;#
 */
enum _mutil_R_class_type
{
  /** R vector */
  MUTIL_R_VECTOR,

  /** R matrix */
  MUTIL_R_MATRIX,

  /** R array */
  MUTIL_R_ARRAY,

  /** R list */
  MUTIL_R_LIST
};


/* See above for documentation on this enum. */
typedef enum _mutil_R_class_type mutil_R_class_type;


/** Collect the R class as supported by MUTILS.
 * Take an R object representing a list, array, matrix,
 * or a numeric vector
 * and return which \Ref{_mutil_R_class_type} it is.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = mutil_type_from_R_class(robj, &cl_type);#
 * @return     Standard MUTILS error/OK code.
 * @param robj Pointer to the input R object.
 * @param type Data class of robj.
 * @see matset_from_R
 * @see matuniv_from_R
 * @see scauniv_from_R
 * @see matuniv_to_R
 * @see _mutil_R_class_type
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_type_from_R_class( SEXP robj,
  mutil_R_class_type *class_type );

/** Convert an R matrix or vector to a universal matrix.
 * Take an R object representing a matrix, or a numeric
 * vector, and convert and store it in a universal matrix
 * which should not be allocated by the calling functions.
 * The universal matrix is allocated in this
 * function to the appropriate size for the input object,
 * and the data are cast to the requested type.
 * The calling function is responsible for freeing the matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = matuniv_from_R(robj, MUTIL_DOUBLE, &outmat);#
 * @return           Standard MUTILS error/OK code.
 * @param robj       Pointer to the input R object.
 * @param out_type   Data type to cast to.
 * @param umat       Pointer to a unallocated universal matrix
 *    for output.
 * @see matset_from_R
 * @see scauniv_from_R
 * @see matset_to_R
 * @see matset_to_R_array
 * @see matset_to_R_list
 * @see matuniv_to_R
 * @see matuniv_to_R_vector
 * @see matuniv_to_R_matrix
 * @see scauniv_to_R
 * @see _mutil_R_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matuniv_from_R( SEXP robj,
  mutil_data_type out_type, univ_mat *umat );

/** Convert an R universal matrix to an R object vector or
 * matrix class. This function returns an error if the class is unknown
 * or conversion doesn't work.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err = matuniv_to_R(&umat, class_type, &ret_obj);#
 * @return           Standard MUTILS error/OK code.
 * @param umat       Pointer to universal matrix.
 * @param class_type R class.
 * @param robj       Pointer to the input R object.
 * @see matset_from_R
 * @see matuniv_from_R
 * @see scauniv_from_R
 * @see matset_to_R
 * @see matset_to_R_array
 * @see matset_to_R_list
 * @see matuniv_to_R_vector
 * @see matuniv_to_R_matrix
 * @see scauniv_to_R
 * @see _mutil_R_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matuniv_to_R( univ_mat *umat,
  mutil_R_class_type class_type, SEXP *robj );

/** Convert a universal matrix to an R vector.
 * Take a universal matrix and convert it to
 * an R object representing a flat vector.
 * The R object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = matuniv_to_R_vector(&mymat, &robj);#
 * @return        Standard MUTILS error/OK code.
 * @param umat    The input universal matrix.
 * @param robj    Pointer for returning the newly allocated R object.
 * @see matset_from_R
 * @see matuniv_from_R
 * @see scauniv_from_R
 * @see matset_to_R
 * @see matset_to_R_array
 * @see matset_to_R_list
 * @see matuniv_to_R
 * @see matuniv_to_R_matrix
 * @see scauniv_to_R
 * @see _mutil_R_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matuniv_to_R_vector( const univ_mat *umat,
  SEXP *robj );

/** Convert a universal matrix to an R matrix.
 * Take a universal matrix and convert it to
 * an R object representing the same matrix.
 * The R object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = matuniv_to_R_matrix(&mymat, &robj);#
 * @return        Standard MUTILS error/OK code.
 * @param umat    The input universal matrix.
 * @param robj    Pointer for returning the newly allocated R object.
 * @see matuniv_to_R_vector
 * @see matuniv_from_R
 * @see scauniv_to_R
 */
MUTIL_WRAPEXPORT mutil_errcode matuniv_to_R_matrix( const univ_mat *umat,
  SEXP *robj );

/** Read an R scalar as a universal scalar.
 * Take an R object containing a scalar, read the number,
 * and cast it to the requested type.  Returns an error if the input is not
 * a scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = scauniv_from_R(robj, type, &usca);#
 * @return        Standard MUTILS error/OK code.
 * @param robj    Pointer to the input R object.
 * @param out_type   Data type to cast to.
 * @param usca     Pointer for returning the universal scalar.
 * @see scauniv_to_R
 * @see double_to_R
 * @see sint32_to_R
 * @see boolean_to_R
 * @see dcomplex_to_R
 * @see double_from_R
 * @see sint32_from_R
 * @see boolean_from_R
 * @see dcomplex_from_R
 */
MUTIL_WRAPEXPORT mutil_errcode scauniv_from_R( SEXP robj,
  mutil_data_type out_type, univ_scalar *usca );

/** Read an R scalar as a double.
 * Take an R object containing a scalar, read the number,
 * and cast it to double.  Returns an error if the input is not
 * a scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = double_from_R(robj, &num);#
 * @return        Standard MUTILS error/OK code.
 * @param robj    Pointer to the input R object.
 * @param num     Pointer for returning the number.
 * @see scauniv_to_R
 * @see double_to_R
 * @see sint32_to_R
 * @see boolean_to_R
 * @see dcomplex_to_R
 * @see scauniv_from_R
 * @see sint32_from_R
 * @see boolean_from_R
 * @see dcomplex_from_R
 */
MUTIL_WRAPEXPORT mutil_errcode double_from_R( SEXP robj,
  double *num );

/** Read an R scalar as a dcomplex.
 * Take an R object containing a scalar complex number,
 * read the number, and cast it to dcomplex.
 * Returns an error if the input is not a complex scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = dcomplex_from_R(robj, &num);#
 * @return        Standard MUTILS error/OK code.
 * @param robj    Pointer to the input R object.
 * @param num     Pointer for returning the number.
 * @see boolean_to_R
 * @see scauniv_to_R
 * @see double_to_R
 * @see sint32_to_R
 * @see boolean_from_R
 * @see scauniv_from_R
 * @see double_from_R
 * @see sint32_from_R
 * @see dcomplex_to_R
 */
MUTIL_WRAPEXPORT mutil_errcode dcomplex_from_R( SEXP robj,
    dcomplex *num );

/** Read an R scalar as an sint32.
 * Take an R object containing a scalar, read the number,
 * and cast it to sint32.  Returns an error if the input is not
 * a scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = sint32_from_R(robj, &num);#
 * @return        Standard MUTILS error/OK code.
 * @param robj    Pointer to the input R object.
 * @param num     Pointer for returning the number.
 * @see scauniv_to_R
 * @see double_to_R
 * @see sint32_to_R
 * @see boolean_to_R
 * @see dcomplex_to_R
 * @see scauniv_from_R
 * @see double_from_R
 * @see boolean_from_R
 * @see dcomplex_from_R
 */
MUTIL_WRAPEXPORT mutil_errcode sint32_from_R( SEXP robj,
  sint32 *num );

/** Read an R scalar as a boolean.
 * Take an R object containing a scalar, read the number,
 * and cast it to boolean.  Returns an error if the input is not
 * a scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = boolean_from_R(robj, &num);#
 * @return        Standard MUTILS error/OK code.
 * @param robj    Pointer to the input R object.
 * @param num     Pointer for returning the number.
 * @see scauniv_to_R
 * @see double_to_R
 * @see sint32_to_R
 * @see boolean_to_R
 * @see dcomplex_to_R
 * @see scauniv_from_R
 * @see double_from_R
 * @see sint32_from_R
 * @see dcomplex_from_R
 */
MUTIL_WRAPEXPORT mutil_errcode boolean_from_R( SEXP robj,
  boolean *num );

/** Convert a universal scalar to an R scalar.
 * Take a number and convert it to
 * an R object representing the same number.
 * The R object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = scauniv_to_R(mynum, &robj);#
 * @return        Standard MUTILS error/OK code.
 * @param usca     The universal scalar.
 * @param robj    Pointer for returning the newly allocated R object.
 * @see double_to_R
 * @see sint32_to_R
 * @see boolean_to_R
 * @see dcomplex_to_R
 * @see scauniv_from_R
 * @see double_from_R
 * @see sint32_from_R
 * @see boolean_from_R
 * @see dcomplex_from_R
 */
MUTIL_WRAPEXPORT mutil_errcode scauniv_to_R( univ_scalar usca,
  SEXP *robj );

/** Convert a double number to an R scalar.
 * Take a number and convert it to
 * an R object representing the same number.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = double_to_R(mynum, &robj);#
 * @return        Standard MUTILS error/OK code.
 * @param num     The number.
 * @param robj    Pointer for returning the newly allocated R object.
 * @see scauniv_to_R
 * @see sint32_to_R
 * @see boolean_to_R
 * @see dcomplex_to_R
 * @see scauniv_from_R
 * @see double_from_R
 * @see sint32_from_R
 * @see boolean_from_R
 * @see dcomplex_from_R
 */
MUTIL_WRAPEXPORT mutil_errcode double_to_R( double num,
  SEXP *robj );

/** Convert a complex number to an R scalar.
 * Take a number and convert it to
 * an R object representing the same number.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = dcomplex_to_R(mynum, &robj);#
 * @return        Standard MUTILS error/OK code.
 * @param num     The number.
 * @param robj    Pointer for returning the newly allocated R object.
 * @see scauniv_to_R
 * @see double_to_R
 * @see sint32_to_R
 * @see boolean_to_R
 * @see scauniv_from_R
 * @see double_from_R
 * @see sint32_from_R
 * @see boolean_from_R
 * @see dcomplex_from_R
 */
MUTIL_WRAPEXPORT mutil_errcode dcomplex_to_R( dcomplex num,
  SEXP *robj );

/** Convert an sint32 number to an R scalar.
 * Take a number and convert it to
 * an R object representing the same number.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = sint32_to_R(mynum, &robj);#
 * @return        Standard MUTILS error/OK code.
 * @param num     The number.
 * @param robj    Pointer for returning the newly allocated R object.
 * @see scauniv_to_R
 * @see double_to_R
 * @see boolean_to_R
 * @see dcomplex_to_R
 * @see scauniv_from_R
 * @see double_from_R
 * @see sint32_from_R
 * @see boolean_from_R
 * @see dcomplex_from_R
 */
MUTIL_WRAPEXPORT mutil_errcode sint32_to_R( sint32 num,
  SEXP *robj );

/** Convert a boolean number to an R scalar.
 * Take a number and convert it to
 * an R object representing the same number.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = boolean_to_R(mynum, &robj);#
 * @return        Standard MUTILS error/OK code.
 * @param num     The number.
 * @param robj    Pointer for returning the newly allocated R object.
 * @see scauniv_to_R
 * @see double_to_R
 * @see sint32_to_R
 * @see dcomplex_to_R
 * @see scauniv_from_R
 * @see double_from_R
 * @see sint32_from_R
 * @see boolean_from_R
 * @see dcomplex_from_R
 */
MUTIL_WRAPEXPORT mutil_errcode boolean_to_R( boolean num,
  SEXP *robj );

/** Convert an R list of matrices or multidimensional array
 * or matrix or vector to a matrix set.
 * Take an R object representing a list of matrices,
 * a multidimensional array, a matrix, or a vector
 * representing a multidimensional dataset,
 * and convert and store it in a matrix set that has not been previously
 * allocated.  The matrix set is allocated in this
 * function to the appropriate size for the input object,
 * and the calling function is responsible for freeing the set.
 * The data are cast to the requested type.
 *
 * If the input is a list of matrices, the output is a
 * flat matrix set of the same length containing those matrices.
 * Nested lists are not supported.
 *
 * If the input is a multidimensional array, then the first two
 * dimensions are the dimensions of the matrices in the set, and
 * the rest of the dimensions become the dimensions of the set itself
 * (e.g. an R array with dimensions (2,3,4,5) would transform into
 * a 4 by 5 set of 20 distinct 2 row by 3 column matrices).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = matset_from_R(robj, type, &outset);#
 * @return        Standard MUTILS error/OK code.
 * @param robj    Pointer to the input R object.
 * @param type    Data type to cast to.
 * @param outset  Pointer to a unallocated matrix set
 *    for output.
 * @see matuniv_from_R
 * @see scauniv_from_R
 * @see matset_to_R
 * @see matset_to_R_array
 * @see matset_to_R_list
 * @see matuniv_to_R
 * @see matuniv_to_R_vector
 * @see matuniv_to_R_matrix
 * @see scauniv_to_R
 * @see _mutil_R_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matset_from_R( SEXP robj,
  mutil_data_type type, mat_set *outset );

/** Convert a matrix set to an R object of vector, matrix,
 * array or list class.
 * This function returns an error if the class is unknown or conversion
 * does not work.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.h
 * @library wrap\_RS
 * @usage  #err = matset_to_R(&set_out, class_type, &ret_obj);#
 * @return        Standard MUTILS error/OK code.
 * @param set        Pointer to matrix set.
 * @param class_type R class.
 * @param robj       Pointer to the input R object.
 * @see matset_from_R
 * @see matuniv_from_R
 * @see scauniv_from_R
 * @see matset_to_R_array
 * @see matset_to_R_list
 * @see matuniv_to_R
 * @see matuniv_to_R_vector
 * @see matuniv_to_R_matrix
 * @see scauniv_to_R
 * @see _mutil_R_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matset_to_R( mat_set *set,
  mutil_R_class_type class_type, SEXP *robj );

/** Convert a matrix set to a list of R matrices.
 * Take a matrix set and convert it to
 * an R object representing the same list of matrices.
 * This function treats the matrix set
 * as if it were one-dimensional, regardless of its dimensions. Therefore,
 * nested lists are not supported. In particular, if you use
 * \Ref{matset_from_R} on the list object returned by this function, you
 * will not recover the original matrix set. Instead, you will have a matrix
 * set with ndims = 1 and *dims = nelem. You might want to use
 * \Ref{matset_resize} to restore the original matrix set size.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = matset_to_R_list(&myset, &robj);#
 * @return        Standard MUTILS error/OK code.
 * @param inmat   The input matrix set.
 * @param robj    Pointer for returning the newly allocated R object.
 * @see matset_from_R
 * @see matuniv_from_R
 * @see scauniv_from_R
 * @see matset_to_R
 * @see matuniv_to_R
 * @see matuniv_to_R_vector
 * @see matuniv_to_R_matrix
 * @see scauniv_to_R
 * @see _mutil_R_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matset_to_R_list( const mat_set *inset,
  SEXP *robj );

/** Extract boundary type from R.
 * Extract the \Ref{_mutil_boundary_type} from an R object.
 * This function calls \Ref{sint32_from_R} to extract an sint32
 * from the R objects, and then executes a switch statements on that
 * number.
 *
 * @limits Care should be taken in the S wrappers to match the enums in the
 * same order as in \_mutil\_boundary\_type. Functions that use only some
 * boundary conditions out of order will return either an error or will end
 * up using the wrong boundary.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = mutil_boundary_from_R(robj,&boundary);#
 * @return       Standard MUTILS error/OK code.
 * @param robj     The R object containing the boundary type
 *    (an integer in R).
 * @param boundary The MUTILS boundary type.
 * @see mutil_R_type
 * @see matset_from_R
 * @see _mutil_boundary_type
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_boundary_from_R( SEXP robj,
    mutil_boundary_type *boundary );

/** Determine if an R object contains NULL.
 * Returns an error if the robj is NULL.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.c
 * @library wrap\_RS
 * @usage  #err_code = null_object_from_R(robj, &bool_null);#
 * @return        Standard MUTILS error/OK code.
 * @param robj    Pointer to the input R object.
 * @param in_null Pointer for returning the boolean flag.
 */
MUTIL_WRAPEXPORT mutil_errcode null_object_from_R( SEXP robj,
  boolean *is_null );

#ifdef __cplusplus
}
#endif

#ifndef __cplusplus

/** Define for C function exported to R.
 * Use this define in place of the extern keyword
 * in the declaration of any function that is to be called from R.
 * It expands to extern for C and extern "C" for C++.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_RS.h
 * @source ut\_RS.h
 * @library wrap\_RS
 * @usage  #EXTERN_R s_object *myfunc( s_object arg1 )#
 */
#if defined(__unix)
#define EXTERN_R extern
#elif defined(__MACH__) && defined(__APPLE__)
#define EXTERN_R extern
#else
#define EXTERN_R extern __declspec(dllexport)
#endif /* ifdef/ifndef __unix */

#else /* ifndef __cplusplus */

#if defined(__unix)
#define EXTERN_R extern "C"
#elif defined(__MACH__) && defined(__APPLE__)
#define EXTERN_R extern "C"
#else
#define EXTERN_R extern "C" __declspec(dllexport)
#endif /* ifdef/ifndef __unix */

#endif /* ifdef/ifndef __cplusplus */


#endif /* ifndef IN_UT_RS_H*/
