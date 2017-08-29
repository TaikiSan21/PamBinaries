
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_univ.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_UNIV_MAT_H
#define IN_UNIV_MAT_H

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "mat_any.h"

/* This file contains functions and macros for universal matrix
   allocation and initialization operations.
 */

#ifdef __cplusplus
extern "C" {
#endif


/** Initialize a universal matrix from a non-universal matrix.
 * This function takes an existing (non-universal) matrix data structure,
 * and puts its information into a universal matrix header.
 * The matrix's underlying flat data array is incorporated directly,
 * and the data is not copied or reallocated.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_univ.h
 * @source mat\_univ.c
 * @library matrix
 * @usage  #err_code = matuniv_wrap_matrix(&my_umat, &my_dblmat, MUTIL_DOUBLE);#
 * @return        Standard mutils error/OK code.
 * @param   umat        Pointer to the universal matrix to initialize.
 * @param   matrix_ptr  Pointer to the matrix to be incorporated.
 * @param   type        Data type of the matrix.
 * @see matuniv_malloc
 * @see matuniv_wrap_data
 * @see _mutil_data_type
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_wrap_matrix( univ_mat *umat,
  void *matrix_ptr, mutil_data_type type );


/** Initialize a universal matrix from a universal matrix.
 * This function takes an existing universal matrix data structure,
 * and puts its information into a universal matrix header.
 * The matrices underlying flat data array is incorporated directly,
 * and the data is not copied or reallocated.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_univ.h
 * @source mat\_univ.c
 * @usage  #errcode = matuniv_wrap_univ_matrix(&umat1, &umat2);#
 * @return        Standard mutils error/OK code.
 * @param   mat1        Pointer to the universal matrix to initialize.
 * @param   mat2        Pointer to the universal matrix to be incorporated.
 * @see matuniv_malloc
 * @see matuniv_wrap_data
 * @see _mutil_data_type
 * @see matuniv_wrap_matrix
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_wrap_univ_matrix( univ_mat *mat1,
  univ_mat *mat2 );


/** Initialize a matrix with given size, type, and data.
 * This function takes a data pointer, size, and type, and
 * and puts this information into a matrix structure.
 * The given data array is incorporated directly,
 * and the data is not copied or reallocated.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_univ.h
 * @source mat\_univ.c
 * @library matrix
 * @usage  #err_code = matuniv_wrap_data(&my_umat, &my_dblarr, nrow, ncol, MUTIL_DOUBLE);#
 * @return    Standard mutils error/OK code.
 * @param   matrix  Pointer to the matrix to initialize.
 * @param   data    Pointer to the data to be incorporated.
 * @param   nrow    Number of rows in the matrix.
 * @param   ncol    Number of columns in the matrix.
 * @param   type    Data type of matrix (specific type versions of this
 *   function omit this argument).
 * @same \begin{itemize}
 *    \item #MUTIL_LIBEXPORT mutil_errcode matcpx_wrap_data(dcomplex_mat *matrix, dcomplex *data, sint32 nrow, sint32 ncol);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode matdbl_wrap_data(double_mat *matrix, double *data, sint32 nrow, sint32 ncol);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode matflt_wrap_data(float_mat *matrix, float *data, sint32 nrow, sint32 ncol);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode matu8_wrap_data(uint8_mat *matrix, uint8 *data, sint32 nrow, sint32 ncol);#
*     \item #MUTIL_LIBEXPORT mutil_errcode matu16_wrap_data(uint16_mat *matrix, uint16 *data, sint32 nrow, sint32 ncol);#
*     \item #MUTIL_LIBEXPORT mutil_errcode matu32_wrap_data(uint32_mat *matrix, uint32 *data, sint32 nrow, sint32 ncol);#
*     \item #MUTIL_LIBEXPORT mutil_errcode mats16_wrap_data(sint16_mat *matrix, sint16 *data, sint32 nrow, sint32 ncol);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode mats32_wrap_data(sint32_mat *matrix, sint32 *data, sint32 nrow, sint32 ncol);#
 * \end{itemize}
 * @see matuniv_wrap_matrix
 * @see matuniv_malloc
 * @see _mutil_data_type
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_wrap_data( univ_mat *matrix,
  void *data, sint32 nrow, sint32 ncol, mutil_data_type type );


/* This function is documented under matuniv_wrap_data in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_wrap_data( dcomplex_mat *matrix,
  dcomplex *data, sint32 nrow, sint32 ncol );


/* This function is documented under matuniv_wrap_data in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_wrap_data( double_mat *matrix,
  double *data, sint32 nrow, sint32 ncol );


/* This function is documented under matuniv_wrap_data in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matflt_wrap_data( float_mat *matrix,
  float *data, sint32 nrow, sint32 ncol );


/* This function is documented under matuniv_wrap_data in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode mats32_wrap_data( sint32_mat *matrix,
  sint32 *data, sint32 nrow, sint32 ncol );


/* This function is documented under matuniv_wrap_data in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode mats16_wrap_data( sint16_mat *matrix,
  sint16 *data, sint32 nrow, sint32 ncol );


/* This function is documented under matuniv_wrap_data in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu32_wrap_data( uint32_mat *matrix,
  uint32 *data, sint32 nrow, sint32 ncol );


/* This function is documented under matuniv_wrap_data in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu16_wrap_data( uint16_mat *matrix,
  uint16 *data, sint32 nrow, sint32 ncol );


/* This function is documented under matuniv_wrap_data in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu8_wrap_data( uint8_mat *matrix,
  uint8 *data, sint32 nrow, sint32 ncol );


/** Initialize a matrix by allocating memory.
 * Allocate uninitialized storage in dynamic memory for a
 * matrix data array of the desired size and type, and put it and the
 * size parameters into the supplied matrix structure.
 * The memory allocated by this function
 * can be freed using the \Ref{matuniv_free} function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_univ.h
 * @source mat\_univ.c
 * @library matrix
 * @usage  #err_code = matuniv_malloc(&matrix, nrow, ncol, MUTIL_DOUBLE);#
 * @return     Standard mutils error/OK code.
 * @param matrix   Pointer to matrix to initialize.
 * @param nrow     Number of rows to create.
 * @param ncol     Number of columns to create.
 * @param type     Data type to create (specific type versions of this
 *   function omit this argument).
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_malloc(double_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_malloc(float_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_malloc(uint8_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_malloc(uint16_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_malloc(uint32_mat *matrix, sint32 nrow,  sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_malloc(sint16_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_malloc(sint32_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matcpx_malloc(dcomplex_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \end{itemize}
 * @see matuniv_wrap_data
 * @see matuniv_realloc
 * @see matuniv_free
 * @see _mutil_data_type
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_malloc( univ_mat *matrix, sint32 nrow,
  sint32 ncol, mutil_data_type type );


/* This function is documented under matuniv_malloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_malloc( double_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_malloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matflt_malloc( float_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_malloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode mats32_malloc( sint32_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_malloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode mats16_malloc( sint16_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_malloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu32_malloc( uint32_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_malloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu16_malloc( uint16_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_malloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu8_malloc( uint8_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_malloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_malloc( dcomplex_mat *matrix, sint32 nrow,
  sint32 ncol );


/** Reallocate a previously allocated matrix.
 * Change the size of a matrix whose data array was
 * previously allocated in dynamic memory, such as one
 * initialized with the \Ref{matuniv_malloc} function.  The data
 * pointer inside the matrix may be moved.  If the total
 * number of matrix elements grows, the new ones will be
 * uninitialized; otherwise, the element values (in row-major
 * order) previously in the matrix will be preserved up to
 * the new size of the matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_univ.h
 * @source mat\_univ.c
 * @library matrix
 * @usage  #err_code = matuniv_realloc(&mymat, newrows, newcols);#
 * @return    Standard mutils error/OK code.
 * @param matrix  Pointer to matrix to be resized.
 * @param nrow    New number of rows for matrix.
 * @param ncol    New number of columns for matrix.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_realloc(double_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_realloc(float_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_realloc(uint8_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_realloc(uint16_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_realloc(uint32_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_realloc(sint16_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_realloc(sint32_mat *matrix, sint32 nrow, sint32 ncol);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matcpx_realloc(dcomplex_mat *matrix, sint32 nrow, sint32 ncol);#
 *  \end{itemize}
 * @see matuniv_malloc
 * @see matuniv_free
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_realloc( univ_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_realloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_realloc( double_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_realloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matflt_realloc( float_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_realloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode mats32_realloc( sint32_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_realloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode mats16_realloc( sint16_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_realloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu32_realloc( uint32_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_realloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu16_realloc( uint16_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_realloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu8_realloc( uint8_mat *matrix, sint32 nrow,
  sint32 ncol );


/* This function is documented under matuniv_realloc in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_realloc( dcomplex_mat *matrix, sint32 nrow,
  sint32 ncol );


/** Free the data in a previously allocated matrix.
 * This function frees previously allocated dynamic memory for
 * the data array in a matrix, typically allocated with
 * the \Ref{matuniv_malloc} function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_univ.h
 * @source mat\_univ.c
 * @library matrix
 * @usage  #err_code = matuniv_free(&mymat);#
 * @return  Standard mutils error/OK code.
 * @param mat  Pointer to matrix to be freed.
 * @same \begin{itemize}
 *     \item #MUTIL_LIBEXPORT mutil_errcode matdbl_free(double_mat *mat);#
 *     \item #MUTIL_LIBEXPORT mutil_errcode matflt_free(float_mat *mat);#
 *     \item #MUTIL_LIBEXPORT mutil_errcode matu8_free(uint8_mat *mat);#
 *     \item #MUTIL_LIBEXPORT mutil_errcode matu16_free(uint16_mat *mat);#
 *     \item #MUTIL_LIBEXPORT mutil_errcode matu32_free(uint32_mat *mat);#
 *     \item #MUTIL_LIBEXPORT mutil_errcode mats16_free(sint16_mat *mat);#
 *     \item #MUTIL_LIBEXPORT mutil_errcode mats32_free(sint32_mat *mat);#
 *     \item #MUTIL_LIBEXPORT mutil_errcode matcpx_free(dcomplex_mat *mat);#
 *  \end{itemize}
 * @see matuniv_malloc
 * @see matuniv_realloc
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_free( univ_mat *mat );


/* This function is documented under matuniv_free in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_free( double_mat *mat );


/* This function is documented under matuniv_free in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matflt_free( float_mat *mat );


/* This function is documented under matuniv_free in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode mats32_free( sint32_mat *mat );


/* This function is documented under matuniv_free in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode mats16_free( sint16_mat *mat );


/* This function is documented under matuniv_free in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu32_free( uint32_mat *mat );


/* This function is documented under matuniv_free in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu16_free( uint16_mat *mat );


/* This function is documented under matuniv_free in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu8_free( uint8_mat *mat );


/* This function is documented under matuniv_free in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_free( dcomplex_mat *mat );


/** Check if a matrix structure is valid.
 * This function checks to see whether a matrix is valid and
 * can be used freely in functions that access its members.  It is
 * valid if it is not null, has a valid type, has a positive number
 * of rows and columns, if the number of elements agrees with the
 * number of rows and columns, and if its data array is not NULL.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_univ.h
 * @source mat\_univ.c
 * @library matrix
 * @usage  #err_code = matuniv_validate(mat_ptr);#
 * @return  Standard mutils error/OK code.
 * @param mat  Pointer to the matrix to be checked.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_validate(const double_mat *mat);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_validate(const float_mat *mat);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_validate(const uint8_mat *mat);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_validate(const uint16_mat *mat);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_validate(const uint32_mat *mat);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_validate(const sint16_mat *mat);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_validate(const sint32_mat *mat);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matcpx_validate(const dcomplex_mat *mat);#
 *  \end{itemize}
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_validate( const univ_mat *mat );


/* This function is documented under matuniv_validate in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matdbl_validate( const double_mat *mat );


/* This function is documented under matuniv_validate in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matflt_validate( const float_mat *mat );


/* This function is documented under matuniv_validate in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode mats32_validate( const sint32_mat *mat );


/* This function is documented under matuniv_validate in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode mats16_validate( const sint16_mat *mat );


/* This function is documented under matuniv_validate in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu32_validate( const uint32_mat *mat );


/* This function is documented under matuniv_validate in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu16_validate( const uint16_mat *mat );


/* This function is documented under matuniv_validate in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matu8_validate( const uint8_mat *mat );


/* This function is documented under matuniv_validate in mat_univ.h */
MUTIL_LIBEXPORT mutil_errcode matcpx_validate( const dcomplex_mat *mat );


/** Verify two universal matrices are the same.
 * This function checks to see whether two universal matrices have the same
 * dimensions and type.
 * It also calls \Ref{matuniv_validate} for each matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_univ.h
 * @source mat\_univ.c
 * @library matrix
 * @usage  #err_code = matuniv_verify_aresame(&mat1, &mat2);#
 * @return        Standard mutils error/OK code.
 * @param mat1  Pointer to the first matrix to be checked.
 * @param mat2  Pointer to the second matrix to be checked.
 * @see matuniv_validate
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_verify_aresame( const univ_mat *mat1,
  const univ_mat *mat2 );


#ifdef __cplusplus
}
#endif

#endif /* ifdef IN_UNIV_MAT_H */
