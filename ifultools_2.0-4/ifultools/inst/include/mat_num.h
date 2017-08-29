
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_num.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_NUM_H_
#define IN_MAT_NUM_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include "img_type.h"

/* This file contains function declarations for matrix
    numerical functions.
*/

#ifdef __cplusplus
extern "C" {
#endif


/** LU decomposition of a square matrix.
 * This function decomposes a square matrix into a product of a
 * lower-triangular matrix with ones on the diagonal (L)
 * and an upper-triangular matrix (U), and returns them both stored
 * in one matrix (with the diagonal elements of L omitted).
 * Since row permutations may be required to perform the decomposition,
 * a vector giving the index permutation is also returned.
 * See William H. Press et al., {\bf Numerical Recipes in C},
 * Second Edition, Cambridge University Press, 1992, p. 43.
 * The input and output matrices may point to the same data location.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_num.h
 * @source mat\_num.c
 * @library matrix
 * @usage  #err_code = matuniv_lu_decomposition(&in_mat, NULL, &indx, &out_mat);#
 * @return        Standard mutils error/OK code.
 * @param in_mat    Pointer to the input matrix.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param indx      Pointer to the output index column vector, of type sint32
 *    and number of rows of the input matrix, to record the
 *    row permutation of the output.
 * @param out_mat   Pointer to output matrix, of the same size and type as
 *    the input matrix.
 * @same
 *   #MUTIL_LIBEXPORT mutil_errcode matdbl_lu_decomposition(const double_mat *in_mat, void *intrp_ptr, uint32_mat *indx, double_mat *out_mat);#
 * @see matuniv_lu_solve
 * @see matuniv_inverse
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_lu_decomposition( const univ_mat *in_mat,
  void *intrp_ptr, univ_mat *indx, univ_mat *out_mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matdbl_lu_decomposition(
  const double_mat *in_mat, void *intrp_ptr, sint32_mat *indx,
  double_mat *out_mat );


/** Solve a linear matrix equation given LU decomposition.
 * This function performs forward substitution and back substitution
 * of a square matrix, after LU decomposition
 * (see \Ref{matuniv_lu_decomposition}),
 * to solve the linear equation A x = B, where A = LU is the LU
 * decomposition.
 * See William H. Press et al., {\bf Numerical Recipes in C},
 * Second Edition, Cambridge University Press, 1992, p. 43.
 * The input and output vectors may point to the same data location.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_num.h
 * @source mat\_num.c
 * @library matrix
 * @usage  #err_code = matuniv_lu_solve(&lu_mat, &in_vec, &indx, NULL, &out_vec);#
 * @return   Standard mutils error/OK code.
 * @param lu_mat    Pointer to matrix containing the LU decomposition of
 *    the input matrix A.
 * @param in_vec    Pointer to the input vector B.
 * @param indx      Pointer to the matrix containing the row permutation
 *    from the LU decomposition.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param out_vec   Pointer to the output vector, x, the solution to A x = B.
 * @same
 *   #MUTIL_LIBEXPORT mutil_errcode matdbl_lu_solve(const double_mat *lu_mat, const double_mat *in_vec, const sint32_mat *indx, void *intrp_ptr, double_mat *out_vec);#
 * @see matuniv_lu_decomposition
 * @see matuniv_inverse
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_lu_solve( const univ_mat *lu_mat,
  const univ_mat *in_vec, const univ_mat *indx, void *intrp_ptr,
  univ_mat *out_vec );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matdbl_lu_solve( const double_mat *lu_mat,
  const double_mat *in_vec, const sint32_mat *indx,
  void *intrp_ptr, double_mat *out_vec );


/** Inverse of a square matrix.
 * This routine inverts a square matrix using LU decomposition.
 * See William H. Press et al., {\bf Numerical Recipes in C},
 * Second Edition, Cambridge University Press, 1992, p. 48.
 * The input and output matrices may point to the same data location.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_num.h
 * @source mat\_num.c
 * @library matrix
 * @usage  #err_code = matuniv_inverse(&in_mat, NULL, &out_mat);#
 * @return        Standard mutils error/OK code.
 * @param in_mat    Pointer to the input square matrix.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param out_mat   Pointer for returning the matrix inverse.
 * @same
 *   #MUTIL_LIBEXPORT mutil_errcode matdbl_inverse(const double_mat *in_mat, void *intrp_ptr, double_mat *out_mat);#
 * @see matuniv_lu_decomposition
 * @see matuniv_lu_solve
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_inverse( const univ_mat *in_mat,
  void *intrp_ptr, univ_mat *out_mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matdbl_inverse( const double_mat *in_mat,
  void *intrp_ptr, double_mat *out_mat );



/** Eigenvectors and eigenvalues of an input real symmetric matrix.
 * This function computes the eigenvectors and eigenvalues of
 * an input real symmetric matrix using Jacobi rotations.
 * See William H. Press et al., {\bf Numerical Recipes in C},
 * Second Edition, Cambridge University Press, 1992, p. 467.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_num.h
 * @source mat\_num.c
 * @library matrix
 * @usage  #err_code = matuniv_eigen_jacobi(&in_mat, NULL, &evals, &evecs);#
 * @return        Standard mutils error/OK code.
 * @param in_mat        Pointer to the input matrix (must be of type
 *                      MUTIL\_DOUBLE)
 * @param intrp_ptr     Pointer for implementation of interrupt checking.
 * @param eigenvalues   Pointer to a one dimensional output matrix containing
 *                      the eigenvalues.
 *                      The length should be equal to the number of rows
 *                      or columns of the input matrix.
 * @param eigenvectors  Pointer to output matrix, of the same size and type as
 *                      the input matrix containing the eigenvectors.
 *                      The eigenvector corresponding to the k'th
 *                      element in the eigenvalue vector is in the k'th
 *                      column in this matrix.
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_eigen_jacobi( const univ_mat *in_mat,
  void *intrp_ptr, univ_mat *eigenvalues, univ_mat *eigenvectors );


/** Sort the eigenvalues and corresponding eigenvectors.
 * This function sorts the eigenvalues into descending order
 * and rearranges the columns of the eigenvectors accordingly.
 * The eigenvalues and eigenvectors can be computed by the
 * \Ref{matuniv_eigen_jacobi} function.
 * See William H. Press et al., {\bf Numerical Recipes in C},
 * Second Edition, Cambridge University Press, 1992, p. 468.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_num.h
 * @source mat\_num.c
 * @library matrix
 * @usage  #err_code = matuniv_eigen_sort(&in_eigenvalues, &in_eigenvectors, intrp_ptr, &out_eigenvalues, &out_eigenvectors);#
 * @return        Standard mutils error/OK code.
 * @param in_eigenvalues     Pointer to the input one-dimensional matrix
 *                           containing the eigenvalues (must be of type
 *                           MUTIL\_DOUBLE).
 * @param in_eigenvectors    Pointer to input matrix containing the
 *                           eigenvectors.
 * @param intrp_ptr          Pointer for implementation of interrupt checking.
 * @param out_eigenvalues    Pointer to output sorted eigenvalue matrix,
 *                           of the same size and type as
 *                           the input eigenvalue matrix.
 *                           Can point to the same location as the input.
 * @param out_eigenvectors   Pointer to output eigenvector matrix,
 *                           of the same size and type as
 *                           the input eigenvector matrix.
 *                           Can point to the same location as the input.
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_eigen_sort(
  const univ_mat *in_eigenvalues, const univ_mat *in_eigenvectors,
  void *intrp_ptr, univ_mat *out_eigenvalues, univ_mat *out_eigenvectors );


#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_NUM_H_*/
