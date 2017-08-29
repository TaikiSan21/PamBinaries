
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_spls5.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_SPLS5_H
#define IN_UT_SPLS5_H

/* This file contains declarations of functions used to convert
   Splus 5.x objects to and from mutils data structures, so that
   Splus functions can use the mutils library.  The functions are
   defined in ut_spls5.c.
*/

#ifdef __cplusplus
extern "C" {
#endif

/* S.h has a typedef for boolean -- ours is different.  If you
   want to use the S.h version, call it LOCAL_S5_boolean in
   any functions that include this header file */

/* Also, make sure S.h is included within the __cplusplus ifdef,
and that it is the first include */

#define boolean LOCAL_S5_boolean
#define NO_NEWIO
#define USE_NEWIO
#include "S.h"
#undef boolean

/* the compiler won't let us cast
   to/from S_complex/dcomplex,
   although the two structures are identical.
   So we use this trick. */

#define dcomplex tmp_dcomplex
#include "ut_type.h"
#undef dcomplex
#define dcomplex S_complex

#include "mat_type.h"
#include "ut_err.h"


/** Splus PROBLEM/ERROR messaging for the trouble message returned
 * by a call to any of the mutils libraries.
 * This macro looks up the appropriate error string message
 * and passes the message to Splus macros PROBLEM and ERROR.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.h
 * @library wrap\_spl
 * @usage  #if(trouble) SPLUS_PROBLEM_ERROR_MSG(trouble);#
 * @see mutil_err_string
 * @see _mutil_errcode
 */
#define SPLUS_PROBLEM_ERROR_MSG( TROUBLE ) \
  { \
    const char    *msgstr; \
    mutil_errcode  othertrouble; \
    othertrouble = mutil_err_string(TROUBLE, &msgstr); \
    if( othertrouble || !msgstr) \
      PROBLEM "error messaging seems to have failed" ERROR; \
    else \
      PROBLEM msgstr ERROR; \
  }


/** Splus PROBLEM/WARN messaging for the trouble message returned
 * by a call to any of the mutils libraries.
 * This macro looks up the appropriate error string message
 * and passes the message to Splus macros PROBLEM and WARN.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.h
 * @library wrap\_spl
 * @usage  #if(trouble) SPLUS_PROBLEM_WARN_MSG(trouble);#
 * @see mutil_err_string
 * @see _mutil_errcode
 */
#define SPLUS_PROBLEM_WARN_MSG( TROUBLE ) \
  { \
    const char    *msgstr; \
    mutil_errcode  othertrouble; \
    othertrouble = mutil_err_string(TROUBLE, &msgstr); \
    if( othertrouble || !msgstr) \
      PROBLEM "warning messaging seems to have failed" ERROR; \
    else \
      PROBLEM msgstr WARN; \
  }


/** Enum of Splus types supported by mutils.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.h
 * @library wrap\_spl
 * @same #typedef enum _mutil_splus_class_type mutil_splus_class_type;#
 */
enum _mutil_splus_class_type
{
  /** Splus vector */
  MUTIL_SPLUS_VECTOR,

  /** Splus matrix */
  MUTIL_SPLUS_MATRIX,

  /** Splus array */
  MUTIL_SPLUS_ARRAY,

  /** Splus list */
  MUTIL_SPLUS_LIST
};


/* See above for documentation on this enum. */
typedef enum _mutil_splus_class_type mutil_splus_class_type;


/** Collect the Splus class as supported by mutils.
 * Take an Splus 5.x object representing a list, array, matrix,
 * or a numeric vector
 * and return which \Ref{_mutil_splus_class_type} it is.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = mutil_type_from_splus_class(sobj, &cl_type);#
 * @return        Standard mutils error/OK code.
 * @param sobj       Pointer to the input Splus object.
 * @param type       Data class of sobj.
 * @see matset_from_splus
 * @see matuniv_from_splus
 * @see scauniv_from_splus
 * @see matuniv_to_splus
 * @see _mutil_splus_class_type
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_type_from_splus_class( s_object *sobj,
  mutil_splus_class_type *class_type );


/** Convert an Splus matrix or vector to a universal matrix.
 * Take an Splus 5.x object representing a matrix,
 * or a numeric vector,
 * and convert and store it in a universal matrix which should not be
 * allocated
 * by the calling functions.  The universal matrix is allocated in this
 * function to the appropriate size for the input object,
 * and the data are cast to the requested type.
 * The calling function is responsible for freeing the matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = matuniv_from_splus(sobj, MUTIL_DOUBLE, &outmat);#
 * @return        Standard mutils error/OK code.
 * @param sobj       Pointer to the input Splus object.
 * @param out_type   Data type to cast to.
 * @param umat       Pointer to a unallocated universal matrix
 *    for output.
 * @see matset_from_splus
 * @see scauniv_from_splus
 * @see matset_to_splus
 * @see matset_to_splus_array
 * @see matset_to_splus_list
 * @see matuniv_to_splus
 * @see matuniv_to_splus_vector
 * @see matuniv_to_splus_matrix
 * @see scauniv_to_splus
 * @see _mutil_splus_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matuniv_from_splus( s_object *sobj,
  mutil_data_type out_type, univ_mat *umat );


/** Convert an Splus universal matrix to an Splus object vector or
 * matrix class. This function returns an error if the class is unknown
 * or conversion doesn't work.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err = matuniv_to_splus(&umat, class_type, &ret_obj);#
 * @return        Standard mutils error/OK code.
 * @param umat       Pointer to universal matrix.
 * @param class_type Splus class.
 * @param sobj       Pointer to the input Splus object.
 * @see matset_from_splus
 * @see matuniv_from_splus
 * @see scauniv_from_splus
 * @see matset_to_splus
 * @see matset_to_splus_array
 * @see matset_to_splus_list
 * @see matuniv_to_splus_vector
 * @see matuniv_to_splus_matrix
 * @see scauniv_to_splus
 * @see _mutil_splus_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matuniv_to_splus( univ_mat *umat,
  mutil_splus_class_type class_type, s_object **sobj );


/** Convert a universal matrix to an Splus vector.
 * Take a universal matrix and convert it to
 * an Splus 5.x object representing a flat vector.
 * The Splus object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = matuniv_to_splus_vector(&mymat, &sobj);#
 * @return        Standard mutils error/OK code.
 * @param umat    The input universal matrix.
 * @param sobj    Pointer for returning the newly allocated Splus object.
 * @see matset_from_splus
 * @see matuniv_from_splus
 * @see scauniv_from_splus
 * @see matset_to_splus
 * @see matset_to_splus_array
 * @see matset_to_splus_list
 * @see matuniv_to_splus
 * @see matuniv_to_splus_matrix
 * @see scauniv_to_splus
 * @see _mutil_splus_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matuniv_to_splus_vector( const univ_mat *umat,
  s_object **sobj );


/** Convert a universal matrix to an Splus matrix.
 * Take a universal matrix and convert it to
 * an Splus 5.x object representing the same matrix.
 * The Splus object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = matuniv_to_splus_matrix(&mymat, &sobj);#
 * @return        Standard mutils error/OK code.
 * @param umat    The input universal matrix.
 * @param sobj    Pointer for returning the newly allocated Splus object.
 * @see matuniv_to_splus_vector
 * @see matuniv_from_splus
 * @see scauniv_to_splus
 */
MUTIL_WRAPEXPORT mutil_errcode matuniv_to_splus_matrix( const univ_mat *umat,
  s_object **sobj );


/** Read an Splus scalar as a universal scalar.
 * Take an Splus 5.x object containing a scalar, read the number,
 * and cast it to the requested type.  Returns an error if the input is not
 * a scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = scauniv_from_splus(sobj, type, &usca);#
 * @return        Standard mutils error/OK code.
 * @param sobj    Pointer to the input Splus object.
 * @param out_type   Data type to cast to.
 * @param usca     Pointer for returning the universal scalar.
 * @see scauniv_to_splus
 * @see double_to_splus
 * @see float_to_splus
 * @see sint32_to_splus
 * @see boolean_to_splus
 * @see dcomplex_to_splus
 * @see double_from_splus
 * @see float_from_splus
 * @see sint32_from_splus
 * @see boolean_from_splus
 * @see dcomplex_from_splus
 */
MUTIL_WRAPEXPORT mutil_errcode scauniv_from_splus( s_object *sobj,
  mutil_data_type out_type, univ_scalar *usca );


/** Read an Splus scalar as a double.
 * Take an Splus 5.x object containing a scalar, read the number,
 * and cast it to double.  Returns an error if the input is not
 * a scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = double_from_splus(sobj, &num);#
 * @return        Standard mutils error/OK code.
 * @param sobj    Pointer to the input Splus object.
 * @param num     Pointer for returning the number.
 * @see scauniv_to_splus
 * @see double_to_splus
 * @see float_to_splus
 * @see sint32_to_splus
 * @see boolean_to_splus
 * @see dcomplex_to_splus
 * @see scauniv_from_splus
 * @see float_from_splus
 * @see sint32_from_splus
 * @see boolean_from_splus
 * @see dcomplex_from_splus
 */
MUTIL_WRAPEXPORT mutil_errcode double_from_splus( s_object *sobj,
  double *num );


/** Read an Splus scalar as a float.
 * Take an Splus 5.x object containing a scalar, read the number,
 * and cast it to float.  Returns an error if the input is not
 * a scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = float_from_splus(sobj, &num);#
 * @return        Standard mutils error/OK code.
 * @param sobj    Pointer to the input Splus object.
 * @param num     Pointer for returning the number.
 * @see scauniv_to_splus
 * @see double_to_splus
 * @see float_to_splus
 * @see sint32_to_splus
 * @see boolean_to_splus
 * @see dcomplex_to_splus
 * @see scauniv_from_splus
 * @see sint32_from_splus
 * @see double_from_splus
 * @see boolean_from_splus
 * @see dcomplex_from_splus
 */
MUTIL_WRAPEXPORT mutil_errcode float_from_splus( s_object *sobj,
  float *num );


/** Read an Splus scalar as a dcomplex.
 * Take an Splus 5.x object containing a scalar complex number,
 * read the number, and cast it to dcomplex.
 * Returns an error if the input is not a complex scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = dcomplex_from_splus(sobj, &num);#
 * @return        Standard mutils error/OK code.
 * @param sobj    Pointer to the input Splus object.
 * @param num     Pointer for returning the number.
 * @see boolean_to_splus
 * @see scauniv_to_splus
 * @see double_to_splus
 * @see float_to_splus
 * @see sint32_to_splus
 * @see boolean_from_splus
 * @see scauniv_from_splus
 * @see double_from_splus
 * @see float_from_splus
 * @see sint32_from_splus
 * @see dcomplex_to_splus
 */
MUTIL_WRAPEXPORT mutil_errcode dcomplex_from_splus( s_object *sobj,
    dcomplex *num );


/** Read an Splus scalar as an sint32.
 * Take an Splus 5.x object containing a scalar, read the number,
 * and cast it to sint32.  Returns an error if the input is not
 * a scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = sint32_from_splus(sobj, &num);#
 * @return        Standard mutils error/OK code.
 * @param sobj    Pointer to the input Splus object.
 * @param num     Pointer for returning the number.
 * @see scauniv_to_splus
 * @see double_to_splus
 * @see float_to_splus
 * @see sint32_to_splus
 * @see boolean_to_splus
 * @see dcomplex_to_splus
 * @see scauniv_from_splus
 * @see double_from_splus
 * @see float_from_splus
 * @see boolean_from_splus
 * @see dcomplex_from_splus
 */
MUTIL_WRAPEXPORT mutil_errcode sint32_from_splus( s_object *sobj,
  sint32 *num );


/** Read an Splus scalar as a boolean.
 * Take an Splus 5.x object containing a scalar, read the number,
 * and cast it to boolean.  Returns an error if the input is not
 * a scalar.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = boolean_from_splus(sobj, &num);#
 * @return        Standard mutils error/OK code.
 * @param sobj    Pointer to the input Splus object.
 * @param num     Pointer for returning the number.
 * @see scauniv_to_splus
 * @see double_to_splus
 * @see float_to_splus
 * @see sint32_to_splus
 * @see boolean_to_splus
 * @see dcomplex_to_splus
 * @see scauniv_from_splus
 * @see double_from_splus
 * @see float_from_splus
 * @see sint32_from_splus
 * @see dcomplex_from_splus
 */
MUTIL_WRAPEXPORT mutil_errcode boolean_from_splus( s_object *sobj,
  boolean *num );


/** Determine if an Splus object contains NULL.
 * Returns an error if the sobj is NULL.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = null_object_from_splus(sobj, &bool_null);#
 * @return        Standard mutils error/OK code.
 * @param sobj    Pointer to the input Splus object.
 * @param in_null Pointer for returning the boolean flag.
 */
MUTIL_WRAPEXPORT mutil_errcode null_object_from_splus( s_object *sobj,
  boolean *is_null );

/** Convert a universal scalar to an Splus scalar.
 * Take a number and convert it to
 * an Splus 5.x object representing the same number.
 * The Splus object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = scauniv_to_splus(mynum, &sobj);#
 * @return        Standard mutils error/OK code.
 * @param usca     The universal scalar.
 * @param sobj    Pointer for returning the newly allocated Splus object.
 * @see double_to_splus
 * @see float_to_splus
 * @see sint32_to_splus
 * @see boolean_to_splus
 * @see dcomplex_to_splus
 * @see scauniv_from_splus
 * @see double_from_splus
 * @see float_from_splus
 * @see sint32_from_splus
 * @see boolean_from_splus
 * @see dcomplex_from_splus
 */
MUTIL_WRAPEXPORT mutil_errcode scauniv_to_splus( univ_scalar usca,
  s_object **sobj );


/** Convert a double number to an Splus scalar.
 * Take a number and convert it to
 * an Splus 5.x object representing the same number.
 * The Splus object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = double_to_splus(mynum, &sobj);#
 * @return        Standard mutils error/OK code.
 * @param num     The number.
 * @param sobj    Pointer for returning the newly allocated Splus object.
 * @see scauniv_to_splus
 * @see float_to_splus
 * @see sint32_to_splus
 * @see boolean_to_splus
 * @see dcomplex_to_splus
 * @see scauniv_from_splus
 * @see double_from_splus
 * @see float_from_splus
 * @see sint32_from_splus
 * @see boolean_from_splus
 * @see dcomplex_from_splus
 */
MUTIL_WRAPEXPORT mutil_errcode double_to_splus( double num,
  s_object **sobj );


/** Convert a float number to an Splus scalar.
 * Take a number and convert it to
 * an Splus 5.x object representing the same number.
 * The Splus object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = float_to_splus(mynum, &sobj);#
 * @return        Standard mutils error/OK code.
 * @param num     The number.
 * @param sobj    Pointer for returning the newly allocated Splus object.
 * @see scauniv_to_splus
 * @see double_to_splus
 * @see sint32_to_splus
 * @see boolean_to_splus
 * @see dcomplex_to_splus
 * @see scauniv_from_splus
 * @see double_from_splus
 * @see float_from_splus
 * @see sint32_from_splus
 * @see boolean_from_splus
 * @see dcomplex_from_splus
 */
MUTIL_WRAPEXPORT mutil_errcode float_to_splus( float num,
  s_object **sobj );


/** Convert a complex number to an Splus scalar.
 * Take a number and convert it to
 * an Splus 5.x object representing the same number.
 * The Splus object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = dcomplex_to_splus(mynum, &sobj);#
 * @return        Standard mutils error/OK code.
 * @param num     The number.
 * @param sobj    Pointer for returning the newly allocated Splus object.
 * @see scauniv_to_splus
 * @see double_to_splus
 * @see float_to_splus
 * @see sint32_to_splus
 * @see boolean_to_splus
 * @see scauniv_from_splus
 * @see double_from_splus
 * @see float_from_splus
 * @see sint32_from_splus
 * @see boolean_from_splus
 * @see dcomplex_from_splus
 */
MUTIL_WRAPEXPORT mutil_errcode dcomplex_to_splus( dcomplex num,
  s_object **sobj );


/** Convert an sint32 number to an Splus scalar.
 * Take a number and convert it to
 * an Splus 5.x object representing the same number.
 * The Splus object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = sint32_to_splus(mynum, &sobj);#
 * @return        Standard mutils error/OK code.
 * @param num     The number.
 * @param sobj    Pointer for returning the newly allocated Splus object.
 * @see scauniv_to_splus
 * @see double_to_splus
 * @see float_to_splus
 * @see boolean_to_splus
 * @see dcomplex_to_splus
 * @see scauniv_from_splus
 * @see double_from_splus
 * @see float_from_splus
 * @see sint32_from_splus
 * @see boolean_from_splus
 * @see dcomplex_from_splus
 */
MUTIL_WRAPEXPORT mutil_errcode sint32_to_splus( sint32 num,
  s_object **sobj );


/** Convert a boolean number to an Splus scalar.
 * Take a number and convert it to
 * an Splus 5.x object representing the same number.
 * The Splus object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = boolean_to_splus(mynum, &sobj);#
 * @return        Standard mutils error/OK code.
 * @param num     The number.
 * @param sobj    Pointer for returning the newly allocated Splus object.
 * @see scauniv_to_splus
 * @see double_to_splus
 * @see float_to_splus
 * @see sint32_to_splus
 * @see dcomplex_to_splus
 * @see scauniv_from_splus
 * @see double_from_splus
 * @see float_from_splus
 * @see sint32_from_splus
 * @see boolean_from_splus
 * @see dcomplex_from_splus
 */
MUTIL_WRAPEXPORT mutil_errcode boolean_to_splus( boolean num,
  s_object **sobj );


/** Convert an Splus list of matrices or multidimensional array
 *    or matrix or vector to a matrix set.
 * Take an Splus 5.x object representing a list of matrices,
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
 * (e.g. an Splus array with dimensions (2,3,4,5) would transform into
 * a 4 by 5 set of 20 distinct 2 row by 3 column matrices).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = matset_from_splus(sobj, type, &outset);#
 * @return        Standard mutils error/OK code.
 * @param sobj    Pointer to the input Splus object.
 * @param type    Data type to cast to.
 * @param outset  Pointer to a unallocated matrix set
 *    for output.
 * @see matuniv_from_splus
 * @see scauniv_from_splus
 * @see matset_to_splus
 * @see matset_to_splus_array
 * @see matset_to_splus_list
 * @see matuniv_to_splus
 * @see matuniv_to_splus_vector
 * @see matuniv_to_splus_matrix
 * @see scauniv_to_splus
 * @see _mutil_splus_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matset_from_splus( s_object *sobj,
  mutil_data_type type, mat_set *outset );


/** Convert an Splus matrix  set to an Splus object of vector, matrix,
 * array or list class.
 * This function returns an error if the class is unknown or conversion
 * doesn't work.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.h
 * @library wrap\_spl
 * @usage  #err = matset_to_splus(&set_out, class_type, &ret_obj);#
 * @return        Standard mutils error/OK code.
 * @param set        Pointer to matrix set.
 * @param class_type Splus class.
 * @param sobj       Pointer to the input Splus object.
 * @see matset_from_splus
 * @see matuniv_from_splus
 * @see scauniv_from_splus
 * @see matset_to_splus_array
 * @see matset_to_splus_list
 * @see matuniv_to_splus
 * @see matuniv_to_splus_vector
 * @see matuniv_to_splus_matrix
 * @see scauniv_to_splus
 * @see _mutil_splus_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matset_to_splus( mat_set *set,
  mutil_splus_class_type class_type, s_object **sobj );


/** Convert a matrix set to a list of Splus matrices.
 * Take a matrix set and convert it to
 * an Splus 5.x object representing the same list of matrices.
 * The Splus object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.  This function treats the matrix set
 * as if it were one-dimensional, regardless of its dimensions. Therefore,
 * nested lists are not supported. In particular, if you use
 * \Ref{matset_from_splus} on the list object returned by this function, you
 * will not recover the original matrix set. Instead, you will have a matrix
 * set with ndims = 1 and *dims = nelem. You might want to use
 * \Ref{matset_resize} to restore the original matrix set size.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = matset_to_splus_list(&myset, &sobj);#
 * @return        Standard mutils error/OK code.
 * @param inmat   The input matrix set.
 * @param sobj    Pointer for returning the newly allocated Splus object.
 * @see matset_from_splus
 * @see matuniv_from_splus
 * @see scauniv_from_splus
 * @see matset_to_splus
 * @see matset_to_splus_array
 * @see matuniv_to_splus
 * @see matuniv_to_splus_vector
 * @see matuniv_to_splus_matrix
 * @see scauniv_to_splus
 * @see _mutil_splus_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matset_to_splus_list( const mat_set *inset,
  s_object **sobj );


/** Convert a matrix set to a multidimensional Splus array.
 * Take a matrix set and convert it to
 * an Splus 5.x object representing the matrices as a multidimensional
 * array. The Splus object is allocated in this
 * function using the NEW macros from S.h, so that it will automatically
 * be freed when the frame is released.
 *
 * The matrices in the set must all have the same dimension and type
 * to be represented as a multidimensional array.  The dimensions
 * of the matrices will become the first two dimensions of the array,
 * and the dimensions of the set will be the remaining dimensions of
 * the array
 * (e.g. an Splus array with dimensions (2,3,4,5) would be created from
 * a 4 by 5 set of 20 distinct 2 row by 3 column matrices).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = matset_to_splus_array(&myset, &sobj);#
 * @return        Standard mutils error/OK code.
 * @param inset   The input matrix set.
 * @param sobj    Pointer for returning the newly allocated Splus object.
 * @see matset_from_splus
 * @see matuniv_from_splus
 * @see scauniv_from_splus
 * @see matset_to_splus
 * @see matset_to_splus_list
 * @see matuniv_to_splus
 * @see matuniv_to_splus_vector
 * @see matuniv_to_splus_matrix
 * @see scauniv_to_splus
 * @see _mutil_splus_class_type
 * @see _mutil_errcode
 */
MUTIL_WRAPEXPORT mutil_errcode matset_to_splus_array( const mat_set *inset,
  s_object **sobj );


/** Find numeric type of Splus data.
 * Find the underlying type of an Splus object.  First check to
 * see if it is an integer, logical, or numeric vector, matrix, or array.
 * Then, if it is a list, look at the first element of the list
 * and make the same check, calling this function recursively.
 * The function will return an error if the data is not numbers.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = mutil_splus_type(sobj, &type);#
 * @return        Standard mutils error/OK code.
 * @param sobj  The S object to check.
 * @param type  The type of the object (generally either sint32 or double).
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_splus_type( const s_object *sobj,
  mutil_data_type *type );


/** Extract boundary type from Splus.
 * Extract the \Ref{_mutil_boundary_type} from an Splus object.
 * This function calls \Ref{sint32_from_splus} to extract an sint32
 * from the Splus objects, and then executes a switch statements on that
 * number.
 *
 * @limits Care should be taken in the S wrappers to match the enums in the
 * same order as in \_mutil\_boundary\_type. Functions that use only some
 * boundary conditions out of order will return either an error or will end
 * up using the wrong boundary.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.c
 * @library wrap\_spl
 * @usage  #err_code = mutil_boundary_from_splus(sobj,&boundary);#
 * @return       Standard mutils error/OK code.
 * @param sobj     The S object containing the boundary type
 *    (an integer in Splus).
 * @param boundary The mutils boundary type.
 * @see mutil_splus_type
 * @see matset_from_splus
 * @see _mutil_boundary_type
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_boundary_from_splus( s_object *sobj,
    mutil_boundary_type *boundary );


#ifdef __cplusplus
}
#endif

#ifndef __cplusplus

/** Define for C function exported to Splus.
 * Use this define in place of the extern keyword
 * in the declaration of any function that is to be called from Splus.
 * It expands to extern for C and extern "C" for C++.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_spls5.h
 * @source ut\_spls5.h
 * @library wrap\_spl
 * @usage  #EXTERN_SPL s_object *myfunc( s_object arg1 )#
 */
#if defined(__unix)
#define EXTERN_SPL extern
#else
#define EXTERN_SPL extern __declspec(dllexport)
#endif /* ifdef/ifndef __unix */

#else /* ifndef __cplusplus */

#if defined(__unix)
#define EXTERN_SPL extern "C"
#else
#define EXTERN_SPL extern "C" __declspec(dllexport)
#endif /* ifdef/ifndef __unix */

#endif /* ifdef/ifndef __cplusplus */


#endif /* ifndef IN_UT_SPLS5_H*/
