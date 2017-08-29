
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_tile.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file */

#ifndef IN_MAT_TILE_H
#define IN_MAT_TILE_H

#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

/* This file contains function declarations for matrix and matrix set
   tiling and untiling functions.
*/

#ifdef __cplusplus
extern "C" {
#endif


/** Tile a matrix into blocks.
 * Tiles a matrix into non-overlapping blocks of neighboring elements,
 * and puts the resulting blocks into a matrix, in the form of row vectors,
 * with one block per row. The input matrix dimensions need not be an
 * exact multiple of the block dimensions. In this case up to
 * block\_nrow - 1 of the last rows and block\_ncol - 1 of the last
 * columns of the original matrix
 * will be left out of the output matrix (see example below).
 * The output matrix must be pre-allocated to have
 * the number of rows equal to the number of blocks that will fit
 * in the input matrix, and the number of columns equal to the number
 * of elements per block.
 *
 * Example: consider the 4x12 matrix, for which the individual blocks
 * are set to the same value for illustration purposes,
 *
 * \begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|}
 * \hline 1 & 1 & 1 & 2 & 2 & 2 & 3 & 3 & 3 & 4 & 4 & 4 \\
 * \hline 1 & 1 & 1 & 2 & 2 & 2 & 3 & 3 & 3 & 4 & 4 & 4 \\
 * \hline 5 & 5 & 5 & 6 & 6 & 6 & 7 & 7 & 7 & 8 & 8 & 8 \\
 * \hline 5 & 5 & 5 & 6 & 6 & 6 & 7 & 7 & 7 & 8 & 8 & 8 \\
 * \hline \end{tabular}
 *
 * The extracted 2x3 blocks will be the rows of a 8x6 matrix,
 *
 * \begin{tabular}{|c|c|c|c|c|c|}
 * \hline 1 & 1 & 1 & 1 & 1 & 1 \\
 * \hline 2 & 2 & 2 & 2 & 2 & 2 \\
 * \hline 3 & 3 & 3 & 3 & 3 & 3 \\
 * \hline 1 & 4 & 4 & 4 & 4 & 4 \\
 * \hline 5 & 5 & 5 & 5 & 5 & 5 \\
 * \hline 6 & 6 & 6 & 6 & 6 & 6 \\
 * \hline 7 & 7 & 7 & 7 & 7 & 7 \\
 * \hline 8 & 8 & 8 & 8 & 8 & 8 \\
 * \hline \end{tabular}
 *
 * If the input matrix dimensions are not exact multiples of the block
 * dimensions, some elements of the input matrix are left out of the
 * output matrix. For example, the following input
 * matrix will generate an output matrix identical to the one above.
 *
 * \begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}
 * \hline 1 & 1 & 1 & 2 & 2 & 2 & 3 & 3 & 3 & 4 & 4 & 4 & 0 & 0 \\
 * \hline 1 & 1 & 1 & 2 & 2 & 2 & 3 & 3 & 3 & 4 & 4 & 4 & 0 & 0 \\
 * \hline 5 & 5 & 5 & 6 & 6 & 6 & 7 & 7 & 7 & 8 & 8 & 8 & 0 & 0 \\
 * \hline 5 & 5 & 5 & 6 & 6 & 6 & 7 & 7 & 7 & 8 & 8 & 8 & 0 & 0 \\
 * \hline 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 * \hline \end{tabular}
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tile.h
 * @source mat\_tile.c
 * @library matrix
 * @usage  #err_code = matuniv_tile(&in_mat, block_nrow, block_ncol, intrp_ptr, &outmat);#
 *
 * @return        Standard mutils error/OK code.
 * @param in_mat     Pointer to matrix to be tiled.
 * @param block_nrow Number of rows per block.
 * @param block_ncol Number of columns per block.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param out_mat    Pointer to output matrix.
 * @same \begin{itemize}
 * \item #MUTIL_LIBEXPORT mutil_errcode matdbl_tile( const double_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, double_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode matflt_tile( const float_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, float_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode mats32_tile( const sint32_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, sint32_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode mats16_tile( const sint16_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, sint16_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode matu32_tile( const uint32_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, uint32_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode matu16_tile( const uint16_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, uint16_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode matu8_tile( const uint8_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, uint8_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode matcpx_tile( const dcomplex_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, dcomplex_mat *out_mat );#
 * \end{itemize}
 * @see matuniv_untile
 * @see matset_tile
 * @see matset_untile
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_tile( const univ_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    univ_mat *out_mat);


/* Function documented with universal type above*/
MUTIL_LIBEXPORT mutil_errcode matdbl_tile( const double_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    double_mat *out_mat );


/* Function documented with universal type above*/
MUTIL_LIBEXPORT mutil_errcode mats32_tile( const sint32_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    sint32_mat *out_mat );


/* Function documented with universal type above*/
MUTIL_LIBEXPORT mutil_errcode mats16_tile( const sint16_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    sint16_mat *out_mat );


/* Function documented with universal type above*/
MUTIL_LIBEXPORT mutil_errcode matu8_tile( const uint8_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    uint8_mat *out_mat );


/* Function documented with universal type above */
MUTIL_LIBEXPORT mutil_errcode matu16_tile( const uint16_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    uint16_mat *out_mat );


/* Function documented with universal type above */
MUTIL_LIBEXPORT mutil_errcode matu32_tile( const uint32_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    uint32_mat *out_mat );


/* Function documented with universal type above */
MUTIL_LIBEXPORT mutil_errcode matflt_tile( const float_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    float_mat *out_mat );


/* Function documented with universal type above */
MUTIL_LIBEXPORT mutil_errcode matcpx_tile( const dcomplex_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    dcomplex_mat *out_mat );


/** Untile a matrix.
 * Given a matrix, it puts each row into a block of specified dimensions
 * contained in the output matrix. This function performs the inverse
 * operation of \Ref{matuniv_tile}. For details on the structure of a tiled and
 * untiled matrix, refer to its documentation.
 * The input and output matrix must have the same number of elements. In
 * addition, the input matrix must have block\_nrow x block\_ncol columns, and
 * a number of rows equal to the number of blocks to fill into the output
 * matrix.  The blocks are filled into the output matrix from left to right,
 * and from top to bottom.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_tile.h
 * @source mat\_tile.c
 * @library matrix
 * @usage  #err_code = matuniv_untile(&in_mat, block_nrow, block_ncol, intrp_ptr, &outmat);#
 *
 * @return        Standard mutils error/OK code.
 * @param in_mat     Pointer to matrix to be untiled.
 * @param block_nrow Number of rows per block.
 * @param block_ncol Number of columns per block.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param out_mat    Pointer to output matrix.
 * @same \begin{itemize}
 * \item #MUTIL_LIBEXPORT mutil_errcode matdbl_untile( const double_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, double_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode matflt_untile( const float_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, float_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode mats32_untile( const sint32_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, sint32_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode mats16_untile( const sint16_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, sint16_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode matu32_untile( const uint32_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, uint32_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode matu16_untile( const uint16_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, uint16_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode matu8_untile( const uint8_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, uint8_mat *out_mat );#
 * \item #MUTIL_LIBEXPORT mutil_errcode matcpx_untile( const dcomplex_mat *in_mat, sint32 block_nrow, sint32 block_ncol, void *intrp_ptr, dcomplex_mat *out_mat );#
 * \end{itemize}
 * @see matuniv_tile
 * @see matset_untile
 * @see matset_tile
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_untile( const univ_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    univ_mat *out_mat);


/* Function documented with universal type above*/
MUTIL_LIBEXPORT mutil_errcode matdbl_untile( const double_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    double_mat *out_mat );


/* Function documented with universal type above*/
MUTIL_LIBEXPORT mutil_errcode mats32_untile( const sint32_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    sint32_mat *out_mat );


/* Function documented with universal type above*/
MUTIL_LIBEXPORT mutil_errcode mats16_untile( const sint16_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    sint16_mat *out_mat );


/* Function documented with universal type above*/
MUTIL_LIBEXPORT mutil_errcode matu8_untile( const uint8_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    uint8_mat *out_mat );


/* Function documented with universal type above */
MUTIL_LIBEXPORT mutil_errcode matu16_untile( const uint16_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    uint16_mat *out_mat );


/* Function documented with universal type above */
MUTIL_LIBEXPORT mutil_errcode matu32_untile( const uint32_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    uint32_mat *out_mat );


/* Function documented with universal type above */
MUTIL_LIBEXPORT mutil_errcode matflt_untile( const float_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    float_mat *out_mat );


/* Function documented with universal type above */
MUTIL_LIBEXPORT mutil_errcode matcpx_untile( const dcomplex_mat *in_mat,
    sint32 block_nrow, sint32 block_ncol, void *intrp_ptr,
    dcomplex_mat *out_mat );


#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_TILE_H */
