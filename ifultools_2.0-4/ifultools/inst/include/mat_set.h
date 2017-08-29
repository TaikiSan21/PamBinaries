
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_set.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef MAT_SET_H_
#define MAT_SET_H_

#include "ut_plat.h"
#include "mat_type.h"

/* This file contains macros and function declarations for
   sets of matrices for the mutils library.
*/

#ifdef __cplusplus
extern "C" {
#endif

/*
*************************
Macros for matrix sets
*************************
*/

/** Calculate matrix set index from dimensional indices.
 * Compose an index into a matrix set's flat data array from a
 * set of dimensional indices.  Matrix sets are always stored in
 * flat arrays in an order where the last dimension varies the fastest.
 * This macro expands to several C statements that compute an
 * index that can be used to pick out a single universal matrix
 * from the set.  No checking is done to see whether or not the given
 * indices are legal.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.h
 * @library matrix
 * @param MATSET   Pointer to a set of matrices.
 * @param INDICES  Array of indices into the set, one per dimension.
 * @param OUTINDEX Variable (sint32) to store the calculated index in.
 * @param ITERATOR Variable (sint32) to use as a loop iterator in code block.
 * @param TMP      Variable (sint32) to use as temporary storage in code block.
 * @usage #MATSET_COMPUTE_INDEX( &myset, myindices, ind, i, tmp); mymatptr = &myset.mats[ind];#
 * @see _mat_set
 */
#define MATSET_COMPUTE_INDEX( MATSET, INDICES, OUTINDEX, ITERATOR, TMP ) \
    OUTINDEX = 0; TMP = 1; \
    for( ITERATOR = (MATSET)->ndim - 1; ITERATOR >= 0; ITERATOR-- ) \
    { \
        OUTINDEX += TMP * (INDICES)[ITERATOR]; \
        TMP *= (MATSET)->dims[ITERATOR]; \
    }


/*
*************************
Memory allocation functions for matrix sets
*************************
*/


/** Initialize a matrix set by allocating memory.
 * Allocate uninitialized storage in dynamic memory for a
 * matrix set of the desired dimensions.  Memory is allocated
 * to store the dimensions and the universal matrix headers in
 * the matrix array, but no storage is created for the universal
 * matrix data (use \Ref{matset_malloc_matrices} for that).
 * The memory can be freed using the \Ref{matset_free} function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_malloc(&matset, ndim, dims);#
 * @return         Standard mutils error/OK code.
 * @param matset   Pointer to matrix set to initialize.
 * @param ndim     Number of dimensions in matrix set.
 * @param dims     Dimensions of matrix set, array of length ndim.
 * @see matset_free
 * @see matset_malloc_matrices
 * @see matset_malloc_matrices_contiguous
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_malloc( mat_set *matset, sint32 ndim,
  const sint32 *dims );


/** Initialize matrices in a matrix set by allocating memory.
 * Allocate uninitialized storage in dynamic memory for the
 * data in all of the universal matrices in a previously initialized (using
 * \Ref{matset_malloc}) matrix set.
 * The memory can be freed using the \Ref{matset_matrices_free} function.
 *
 * Individual matrices in the matrix set can also be indexed
 * using the \Ref{MATSET_COMPUTE_INDEX} macro, and then initialized
 * using the \Ref{matuniv_wrap_matrix}, \Ref{matuniv_wrap_data},
 * or \Ref{matuniv_malloc} functions, for instance.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_malloc_matrices(&matset, nrow, ncol, MUTIL_DOUBLE);#
 * @return         Standard mutils error/OK code.
 * @param matset   Pointer to matrix set to initialize.
 * @param nrow     Number of rows to create.
 * @param ncol     Number of columns to create.
 * @param type     Data type to create.
 * @see matset_malloc
 * @see matset_malloc_matrices_contiguous
 * @see matset_matrices_free
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_malloc_matrices( mat_set *matset,
  sint32 nrow, sint32 ncol, mutil_data_type type );


/** Initialize matrices in a matrix set by allocating memory.
 * Allocate uninitialized storage in dynamic memory for the
 * data in all of the universal matrices in a previously initialized (using
 * \Ref{matset_malloc}) matrix set.
 *
 * Unlike \Ref{matset_malloc_matrices}, the storage for the matrices are
 * allocated together in one contiguous memory block.
 *
 * The memory must be freed using the \Ref{matset_matrices_free}
 * function. Memory should never be freed or realloc'ed per individual
 * matrices in the matrix set.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_malloc_matrices_contiguous(&matset, nrow, ncol, MUTIL_DOUBLE);#
 * @return         Standard mutils error/OK code.
 * @param matset   Pointer to matrix set to initialize.
 * @param nrow     Number of rows to create.
 * @param ncol     Number of columns to create.
 * @param type     Data type to create.
 * @see matset_malloc
 * @see matset_malloc_matrices
 * @see matset_matrices_free
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_malloc_matrices_contiguous(
  mat_set *matset, sint32 nrow, sint32 ncol, mutil_data_type type );


/** Initialize matrices in a matrix set by allocating memory.
 * Allocate uninitialized storage in dynamic memory for the
 * data in all of the universal matrices in a previously initialized (using
 * \Ref{matset_malloc}) matrix set. The matrices may have arbitrary
 * dimensions, so their sizes need not be all the same as in
 * \Ref{matset_malloc_matrices}.
 * The memory can be freed using the \Ref{matset_matrices_free} function.
 *
 * Individual matrices in the matrix set can also be indexed
 * using the \Ref{MATSET_COMPUTE_INDEX} macro, and then initialized
 * using the \Ref{matuniv_wrap_matrix}, \Ref{matuniv_wrap_data},
 * or \Ref{matuniv_malloc} functions, for instance.
 *
 * The number of rows and columns are passed to the function
 * in two matrices, which must have the same number of elements
 * as the number of matrices in the matrix set.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_malloc_matrices_arbitrary_size(&matset, &nrow, &ncol, MUTIL_DOUBLE);#
 * @return         Standard mutils error/OK code.
 * @param matset   Pointer to matrix set to initialize.
 * @param nrow_mat Pointer to an sint32 matrix containing the number
 *                 of rows of each matrix to create.
 * @param ncol_mat Pointer to an sint32 matrix containing the number
 *                 of columns of each matrix to create.
 * @param type     Data type to create.
 * @see matset_malloc_matrices
 * @see matset_malloc_matrices_contiguous
 * @see matset_matrices_free
 * @see matset_malloc
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_malloc_matrices_arbitrary_size(
    mat_set *matset, const sint32_mat *nrow_mat, const sint32_mat *ncol_mat,
    mutil_data_type type );


/** Initialize a flat (1d) matrix set by allocating memory for it and its
 * matrices. Allocate uninitialized storage in dynamic memory for a
 * one dimensional matrix set.  Memory is allocated
 * to store the dimensions and the universal matrix headers in
 * the matrix array using \Ref{matset_malloc}, and
 * storage is created for the universal matrix data using
 * \Ref{matset_malloc_matrices}.
 * The memory can be freed using the \Ref{matset_matrices_free}
 * and \Ref{matset_free} functions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_malloc_flat(&matset, nslice, nrow, ncol, MUTIL_DOUBLE, TRUE);#
 * @return         Standard mutils error/OK code.
 * @param matset   Pointer to matrix set to initialize.
 * @param nelem    Number of matrices in 1-d matrix set.
 * @param nrow     Number of rows for each matrix in the set.
 * @param ncol     Number of columns for each matrix in the set.
 * @param type     Data type to create.
 * @param contiguous   If TRUE, storage for the matrices are
 *     allocated in one contiguous memory block.
 *
 * @see matset_malloc
 * @see matset_malloc_matrices
 * @see matset_malloc_matrices_contiguous
 * @see matset_matrices_free
 * @see matset_free
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_malloc_flat(
  mat_set         *matset,
  sint32           nelem,
  sint32           nrow,
  sint32           ncol,
  mutil_data_type  type,
  boolean          contiguous );


/** Free a dynamically-allocated matrix set.
 * Free a matrix set previously allocated by \Ref{matset_malloc}.
 * The storage for the dimensions and the universal matrix headers in
 * the matrix array is freed, but not the storage for the matrix
 * data in the individual universal matrices (that is freed
 * using \Ref{matset_matrices_free}.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_free(&matset);#
 * @return         Standard mutils error/OK code.
 * @param matset   Pointer to matrix set to free.
 * @see matset_malloc
 * @see matset_matrices_free
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_free( mat_set *matset );


/** Free data storage for matrices in a matrix set.
 * Free dynamic memory previously allocated using \Ref{matset_malloc_matrices}
 * for the data in all of the universal matrices matrix set.  The
 * matrix set's memory itself can be freed using the \Ref{matset_free}
 * function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_matrices_free(&matset);#
 * @return         Standard mutils error/OK code.
 * @param matset   Pointer to matrix set to free matrices of.
 * @see matset_free
 * @see matset_malloc_matrices
 * @see matset_malloc_matrices_contiguous
 * @see matset_malloc_matrices_arbitrary_size
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_matrices_free( mat_set *matset );


/** Change the dimensions of a matrix set.
 * Change the dimensions of a matrix set, while keeping the
 * total number of matrices the same. This function reallocates the memory
 * in the dims array of the \Ref{_mat_set} data structure to accommodate
 * the new dimensions and changes the value of dims in the matrix set. It
 * checks that the total number of matrices does not change.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_resize(ndims, dims, &matset);#
 * @return        Standard mutils error/OK code.
 * @param ndims  Number of new dimensions.
 * @param dims   Pointer to array containing the new dimensions. Its length
 *               must be equal to ndims.
 * @param matset Pointer to the matrix set to be resized.
 * @see matset_malloc
 * @see matset_malloc_matrices
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_resize( sint32 ndims, sint32 *dims,
  mat_set *matset);


/*
*************************
Validation functions for matrix sets
*************************
*/


/** Check if the members of a matrix set structure are valid.
 * This function checks to see whether a matrix set is valid and
 * can be used freely in functions that access its members.  It is
 * valid if it is not null, has non-null matrix array and dimensions
 * pointers, its total number of elements is equal to the product
 * of the dimensions, and the dimensions are all positive.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_validate(matset_ptr);#
 * @return        Standard mutils error/OK code.
 * @param matset  Pointer to the matrix set to be checked.
 * @see matset_validate_matrices
 * @see matset_verify_allsame
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_validate( const mat_set *matset );


/** Check if the members of a matrix set structure are valid
 * and if all matrices in a matrix set are valid.
 * This function calls \Ref{matset_validate} to check the matrix
 * set structure, and \Ref{matuniv_validate} to check each
 * universal matrix in the matrix set.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_validate_matrices(matset_ptr);#
 * @return        Standard mutils error/OK code.
 * @param matset  Pointer to the matrix set to be checked.
 * @see matset_validate
 * @see matuniv_validate
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_validate_matrices( const mat_set *matset );


/** Verify that a matrix set is all the same type and size.
 * This function checks to see whether a matrix set contains
 * matrices that are all the same type and size.  This need not
 * be true in general in order to have a valid matrix set,
 * but it is required by some matrix set functions.  It calls
 * \Ref{matset_validate_matrices} also, so it is not necessary to call both
 * functions to validate a same-size set.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_verify_allsame(&myset);#
 * @return        Standard mutils error/OK code.
 * @param matset  Pointer to the matrix set to be checked.
 * @see matset_validate_matrices
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_verify_allsame( const mat_set *matset );


/** Verify two matrix sets have same dimensions and same matrices.
 * This function checks to see whether two matrix sets have the same
 * dimensions and contain matrices that are all the same type and size.
 * It calls \Ref{matset_validate_matrices} also, so it is not necessary to
 * call both functions to validate that two matrix sets are all the same.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_verify_aresame(&myset1, &myset2);#
 * @return        Standard mutils error/OK code.
 * @param matset1  Pointer to the first matrix set to be checked.
 * @param matset2  Pointer to the second matrix set to be checked.
 * @see matset_validate_matrices
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_verify_aresame( const mat_set *matset1,
  const mat_set *matset2 );


/** Verify two matrix sets have same dimensions and same matrix dimensions;
 * the two set types may differ.
 * This function checks to see whether each matrix set has matrices that
 * are all the same size and type, and whether the two matrix sets have the
 * same dimensions and matrix dimensions.  The type of matset1 may be
 * different from that of matset2.
 *
 * It calls \Ref{matset_verify_allsame} also, so it is not necessary to
 * call both functions to validate that two matrix sets are all the same.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_verify_same_dims(&myset1, &myset2);#
 * @return        Standard mutils error/OK code.
 * @param matset1  Pointer to the first matrix set to be checked.
 * @param matset2  Pointer to the second matrix set to be checked.
 * @see matset_verify_allsame
 * @see matset_validate_matrices
 * @see _mat_set
 */
MUTIL_LIBEXPORT mutil_errcode matset_verify_same_dims( const mat_set *matset1,
  const mat_set *matset2 );


/** Cast a matrix set to another type.  Take a matrix set with
 * matrices of one type and cast its data into another matrix set
 * with different types of matrices, rounding if necessary.
 * Implemented as a series of calls to \Ref{matuniv_cast}.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @return          Standard mutils error/OK code.
 * @usage #err_code = matset_cast(&mat, TRUE, intrp_ptr, &result);#
 *
 * @param matset    Pointer to the matrix set to cast.
 * @param clip      if TRUE and input matrix values are illegal for the
 *   output type, clip to fit; if FALSE, return an error.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param result    Pointer to output matrix set of same size.
 *
 * @see matuniv_cast
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matset_cast(
  const mat_set  *matset,
  boolean         clip,
  void           *intrp_ptr,
  mat_set        *result );


/** Assign a matrix set to another.  Take a matrix set with
 * matrices and assign its data into another matrix set
 * with different matrices of same type and size.
 * Implemented as a series of calls to \Ref{matuniv_assign}.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @return          Standard mutils error/OK code.
 * @usage #err_code = matset_set(&mat, intrp_ptr, &result);#
 *
 * @param matset    Pointer to the matrix set to copy.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param result    Pointer to output matrix set of same size.
 *
 * @see matuniv_assign
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matset_assign(
  const mat_set  *matset,
  void           *intrp_ptr,
  mat_set        *result );


/*
*************************
Operations on matrix sets
*************************
*/


/** Make a matrix set into a tiled matrix.
 * Take a one- or two-dimensional set of matrices, all of the same type and
 * size, and put its values into a single tiled matrix of the same type and
 * correct size.
 *
 * A tiled matrix is a matrix consisting of tr rows and tc columns
 * of (nr x nc)-sized non-overlapping rectangular blocks (tiles).
 * The matrix set for the untiled matrices is a set of (tr x tc)
 * matrices, and the set has dimensions (nr x nc).
 * The values in the tile at position (ti, tj)
 * within the tiled matrix correspond to the values at
 * position (ti, tj) in the matrices in the set, with
 * the value within the tile at position (i, j) in the
 * (i, j)-th matrix in the set.
 *
 * Example 1: A 4x6 matrix with 2-row, 3-column tiles:
 *
 * \begin{tabular}{|c|c|c|c|c|c|}
 * \hline 1 & 2 & 3 &  7 & 8 & 9 \\
 * \hline 4 & 5 & 6 &  10 & 11 & 12 \\
 * \hline 13 & 14 & 15  & 19 & 20 & 21 \\
 * \hline 16 & 17 & 18  & 22 & 23 & 24 \\
 * \hline \end{tabular}
 *
 * corresponds to 2x3 set of 2x2 matrices, that looks like this:
 *
 * \begin{tabular}{|c|c|c|c|c|c|c|c|}
 * \hline 1 & 7 & & 2 & 8 & & 3 & 9  \\
 * \hline 13 & 19 & & 14 & 20 & & 15 & 21 \\
 * \hline
 * \hline 4 & 10 & & 5 & 11 & & 6 & 12 \\
 * \hline 16 & 22 & & 17 & 23 & & 18 & 24 \\
 * \hline \end{tabular}
 *
 * Example 2: an (n x m) RGB color image can be stored as a tiled
 * (n x 3 m) matrix or a (1 x 3) set of (n x m) matrices, and the
 * R, G, and B values for the pixel at position (i, j) are in
 * tile (i, j) in the tiled matrix (i.e. row i, and columns 3j,
 * 3j + 1, and 3j + 2), and in the matrix set they are in the
 * first, second, and third matrices at position (i, j).  That is,
 * the matrix set has one matrix for R, one for G, and one for B,
 * and in the tiled matrix, the values are stored with the first
 * pixel's R, G, and B, and then the second pixel's R, G, and B,
 * etc.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_tile(&myset, intrp_ptr, &outmat);#
 * @return        Standard mutils error/OK code.
 * @param matset     Pointer to the matrix set to tile.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param umat       Pointer to allocated matrix for output.
 * @see matset_untile
 * @see matuniv_tile
 * @see matuniv_untile
 * @see _mat_set
 * @see _univ_mat
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matset_tile( const mat_set *matset,
  void *intrp_ptr, univ_mat *umat );



/** Make a tiled matrix into a set of matrices.
 * Take an tiled matrix and put the values,
 * untiled, into a one- or two-dimensional set of matrices.
 * The set must be pre-allocated to the correct type and size.
 * See the documentation for \Ref{matset_tile} for
 * more information on tiled matrices.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_set.h
 * @source mat\_set.c
 * @library matrix
 * @usage  #err_code = matset_untile(&mymat, intrp_ptr, &outset);#
 * @return        Standard mutils error/OK code.
 * @param umat       Pointer to matrix to untile.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param matset     Pointer to the matrix set for output.
 * @see matset_tile
 * @see matuniv_tile
 * @see matuniv_untile
 * @see _mat_set
 * @see _univ_mat
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matset_untile( const univ_mat *umat,
  void *intrp_ptr, mat_set *matset );



#ifdef __cplusplus
}
#endif

#endif /* MAT_SET_H_ */
