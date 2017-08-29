
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_intp.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef MAT_INTP_H_
#define MAT_INTP_H_

#include "mat_set.h"

/* This file contains macros and function declarations for
   interpolation functions for matrices and sets of matrices
   for the mutils library.
 */

#ifdef __cplusplus
extern "C" {
#endif


/*
*****************************************
Operations for matrix interpolation
*****************************************
*/


/** Interpolate between two matrices.
 * Take two matrices of the same type and size, and output
 * a set of interpolated matrices, where each output matrix
 * in the set is a linear combination of the two input matrices,
 * with the given weights for the first matrix, and 1 minus the given
 * weights for the second matrix.  If weights are negative or
 * greater than one, extrapolation will result. For integer valued
 * output data, interpolated values will be rounded.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intp.h
 * @source mat\_intp.c
 * @library matrix
 * @usage  #err_code = matset_interpolate(&mat1, &mat2, &weights, intrp_ptr, &outset);#
 * @return        Standard mutils error/OK code.
 * @param mat1       Pointer to first input matrix.
 * @param mat2       Pointer to second input matrix, same size and
 *   type as first.
 * @param weights    Pointer to 1-column matrix of weights for the
 *   first matrix, of type MUTIL\_DOUBLE and same length as number of elements
 *   in the output matset.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param matset     Pointer to a previously allocated 1-D set of matrices of
 *   same size and type as input matrices, for output.
 * @see _mat_set
 * @see matset_interpolate_cubic_voxels
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matset_interpolate( const univ_mat *mat1,
  const univ_mat *mat2, const univ_mat *weights, void *intrp_ptr,
  mat_set *matset );


/** Interpolate between non-cubic voxels of a matrix set
 * to create cubic voxels.
 * Take a matrix set where the voxel is not a cube, and use linear
 * interpolation to output a matrix set of the same data type where
 * the voxel is a cube of unit dimensions.  For integer valued input data,
 * output interpolated values will be rounded.
 *
 * @limits  Only interpolation across two consecutive image slices is
 * supported. Therefore, it is required that row\_dimension be equal to
 * column\_dimension be equal to 1.
 * The only interpolation type supported is MUTIL\_INTERPOLATION\_LINEAR.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_intp.h
 * @source mat\_intp.c
 * @library matrix
 * @usage  #err_code = matset_interpolate_cubic_voxels(&inset, row_dim, col_dim, slice_dim, MUTIL_INTERPOLATION_LINEAR, intrp_ptr, &outset);#
 * @return           Standard mutils error/OK code.
 * @param matset_in         Pointer to input one-dimension matrix set
 *   of same size and type.
 * @param row_dimension     Height of each row of the input matrix set.
 * @param column_dimension  Width of each column of the input matrix set.
 * @param slice_dimension   Depth of each slice of the input matrix set.
 * @param interpolation     Algorithm for interpolation.
 * @param intrp_ptr         Pointer for implementation of interrupt checking.
 * @param matset_out        Pointer to an unallocated but declared matrix set.
 * @see _mat_set
 * @see matset_interpolate
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matset_interpolate_cubic_voxels(
  const mat_set                  *matset_in,
  const double                    row_dimension,
  const double                    column_dimension,
  const double                    slice_dimension,
  const mutil_interpolation_type  interpolation,
  void                           *intrp_ptr,
  mat_set                        *matset_out );


#ifdef __cplusplus
}
#endif

#endif /* MAT_INTP_H_ */
