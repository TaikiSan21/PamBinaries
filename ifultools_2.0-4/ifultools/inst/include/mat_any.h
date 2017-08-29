
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_any.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_ANY_MAT_H
#define IN_ANY_MAT_H

/* This file contains macros that can be used on any of the
   basic matrix types.
*/


#ifdef __cplusplus
extern "C" {
#endif


/** Check if a matrix row number is valid.
 * Check to see if a number to be used as a matrix row number
 * is within the bounds of the matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_any.h
 * @source mat\_any.h
 * @library
 * @return True (in bounds) or false (out of bounds).
 * @param matrix Pointer to a matrix.
 * @param row    The row number to check.
 * @usage #ok = MATANY_CHECK_ROW(&my_uint8_mat, 7);#
 * @see Matrix Data Types
 * @see MATANY_CHECK_COL
 */
#define MATANY_CHECK_ROW(matrix, row) ((row) >= 0 && (row) < (matrix)->nrow)


/** Check if a matrix column number is valid.
 * Check to see if a number to be used as a matrix column number
 * is within the bounds of the matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_any.h
 * @source mat\_any.h
 * @library
 * @return True (in bounds) or false (out of bounds).
 * @param matrix Pointer to a matrix.
 * @param col    The column number to check.
 * @usage #ok = MATANY_CHECK_COL(&my_uint8_mat, 7);#
 * @see Matrix Data Types
 * @see MATANY_CHECK_ROW
 */
#define MATANY_CHECK_COL(matrix, col) ((col) >= 0 && (col) < (matrix)->ncol)


/** Calculate matrix index from row and column.
 * Compose an index into a matrix's flat data array from a row and
 * column number. Matrices are always stored in flat arrays in
 * row-major (raster scan) order.  No checking is done on whether
 * or not the given row and column are within the array bounds.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_any.h
 * @source mat\_any.h
 * @library
 * @return The desired index.
 * @param matrix Pointer to a matrix.
 * @param row    The row number.
 * @param col    The column number.
 * @usage #x31 = my_dbl_mat.data[MATANY_INDEX(&my_dbl_mat, 3, 1)]#
 * @see Matrix Data Types
 * @see MATANY_CHECK_ROW
 * @see MATANY_CHECK_COL
 */
#define MATANY_INDEX(matrix, row, col) ( (row) * (matrix)->ncol + (col))

/** Check if two matrices have the same number of rows and columns.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_any.h
 * @source mat\_any.h
 * @library
 * @return True (same dimensions) or false (different dimensions)
 * @param mat1 Pointer to first matrix.
 * @param mat2 Pointer to second matrix.
 * @usage #test = MATANY_EQUAL_DIM(&mat1, &mat2);#
 * @see Matrix Data Types
 */
#define MATANY_EQUAL_DIM(mat1, mat2) ((mat1)->nrow == (mat2)->nrow && \
                                   (mat1)->ncol == (mat2)->ncol )


/** Check if a matrix is a row vector.
 * Takes a pointer to any matrix type and checks if the
 * matrix is a row (horizontal) vector.
 * One-element matrices are considered vectors.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_any.h
 * @source mat\_any.h
 * @library
 * @return TRUE (row vector) or FALSE (not a row vector).
 * @param mat Pointer to matrix.
 * @usage #test = MATANY_IS_VEC_ROW(&mat);#
 * @see Matrix Data Types
 * @see MATANY_IS_VEC
 * @see MATANY_IS_VEC_COL
 */
#define MATANY_IS_VEC_ROW(mat) ( (mat)->nrow == 1 && (mat)->ncol >= 1 )


/** Check if a matrix is a column vector.
 * Takes a pointer to any matrix type and checks if the
 * matrix is a column (vertical) vector.
 * One-element matrices are considered vectors.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_any.h
 * @source mat\_any.h
 * @library
 * @return TRUE (column vector) or FALSE (not a column vector).
 * @param mat Pointer to matrix.
 * @usage #test = MATANY_IS_VEC_COL(&mat);#
 * @see Matrix Data Types
 * @see MATANY_IS_VEC
 * @see MATANY_IS_VEC_ROW
 */
#define MATANY_IS_VEC_COL(mat) ( (mat)->nrow >= 1 && (mat)->ncol == 1 )


/** Check if a matrix is a one-dimensional vector.
 * Takes a pointer to any matrix type and checks if the
 * matrix is either a column or row vector.
 * One-element matrices are considered vectors.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_any.h
 * @source mat\_any.h
 * @library
 * @return True (row or column vector) or false (at least 2x2 dimensions).
 * @param mat Pointer to matrix.
 * @usage #test = MATANY_IS_VEC(&mat);#
 * @see Matrix Data Types
 * @see MATANY_IS_VEC_ROW
 * @see MATANY_IS_VEC_COL
 */
#define MATANY_IS_VEC(mat) ( MATANY_IS_VEC_ROW(mat) || MATANY_IS_VEC_COL(mat) )


#ifdef __cplusplus
}
#endif

#endif /* ifdef IN_ANY_MAT_H */
