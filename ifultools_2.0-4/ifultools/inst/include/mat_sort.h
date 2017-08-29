
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_sort.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_SORT_H_
#define IN_MAT_SORT_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

/* This file contains declarations for functions for sorting matrices
   and table lookup which are implemented in mat_sort.c
 */

#ifdef __cplusplus
extern "C" {
#endif


/**************
 Universal matrix functions
 **************/


/** Perform a partial sort of a matrix, returning sort indices.
 * Calculate the index permutation needed to partially sort a
 * matrix, considering the matrix data as a single flat vector.
 * Subscripting the data vector in the matrix with the returned
 * indices will result in data that contains all the input matrix
 * elements, in partially sorted ascending order.
 *
 * If the second argument to the sorting function is NULL, then the
 * result is a fully sorted vector.  If the second argument is not
 * NULL, then the result is only guaranteed to have the data in the
 * index positions from the second argument in their correct positions,
 * and the rest of the data in the correct gaps between the guaranteed
 * elements, but not necessarily in the correct order within the gaps.
 *
 * The sort is performed by the well-known quicksort algorithm
 * (see, for instance, the discussion in Chapter 8 of {\bf
 * Numerical Recipes in C, Second Edition}, W. H. Press, S. A. Teukolsky,
 * W. T. Vetterling, and B. P. Flannery, Cambridge University Press,
 * 1992), choosing the partitioning element by median selection,
 * and sorting sub-arrays of length 9 or less by insertion sort.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_sort.h
 * @source mat\_sort.c
 * @library matrix
 * @usage  #err_code = matuniv_sort_index_partial(&mymat, &partials, intrp_ptr, &out);#
 * @return          Standard mutils error/OK code.
 * @param mat       Pointer to the matrix to sort.
 * @param needed    Pointer to matrix (sint32 type) of indices that must
 *    be in right position, or NULL to perform a full sort.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param index     Pointer for matrix (sint32 type) of returned sort indices.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_sort_index_partial(const double_mat *mat, const sint32_mat *needed, void *intrp_ptr, sint32_mat *index);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_sort_index_partial(const float_mat *mat, const sint32_mat *needed, void *intrp_ptr, sint32_mat *index);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_sort_index_partial(const uint8_mat *mat, const sint32_mat *needed, void *intrp_ptr, sint32_mat *index);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_sort_index_partial(const uint16_mat *mat, const sint32_mat *needed, void *intrp_ptr, sint32_mat *index);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_sort_index_partial(const sint16_mat *mat, const sint32_mat *needed, void *intrp_ptr, sint32_mat *index);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_sort_index_partial(const uint32_mat *mat, const sint32_mat *needed, void *intrp_ptr, sint32_mat *index);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_sort_index_partial(const sint32_mat *mat, const sint32_mat *needed, void *intrp_ptr, sint32_mat *index);#
 *  \end{itemize}
 * @see matuniv_sort
 * @see matuniv_permute
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_sort_index_partial( const univ_mat *mat,
  const univ_mat *needed, void *intrp_ptr, univ_mat *index );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matdbl_sort_index_partial( const double_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matflt_sort_index_partial( const float_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matu32_sort_index_partial( const uint32_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode mats32_sort_index_partial( const sint32_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matu16_sort_index_partial( const uint16_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode mats16_sort_index_partial( const sint16_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matu8_sort_index_partial( const uint8_mat *mat,
  const sint32_mat *needed, void *intrp_ptr, sint32_mat *index );


/** Permute a matrix.
 * Copy elements of a matrix to another matrix in a different order
 * given by a permutation index (e.g. the output of the
 * \Ref{matuniv_sort_index_partial} function).  That is,
 * element i in the output matrix is assigned the value of the input
 * matrix element whose index is the ith element in the permutation index,
 * treating all three matrices as flat arrays.
 *
 * Repeated indices are allowed in the permutation index, as no attempt
 * is made to assure that all elements in the input matrix are used in the
 * output.  The input and output matrices must be the same type, and
 * the index and output matrices must be the same size, with elements in
 * the index matrix positive numbers less than the number of elements in
 * the input matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_sort.h
 * @source mat\_sort.c
 * @library matrix
 * @usage  #err_code = matuniv_permute(&mymat, &index, intrp_ptr, &out);#
 * @return        Standard mutils error/OK code.
 * @param mat       Pointer to the matrix to permute.
 * @param index     Pointer to the index matrix (same size as output matrix,
 *   and of type sint32).
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param out       Pointer to returned permuted matrix (same size
 *    as index matrix, same type as input matrix).
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_permute(const double_mat *mat, const sint32_mat *index, void *intrp_ptr, double_mat *out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_permute(const float_mat *mat, const sint32_mat *index, void *intrp_ptr, double_mat *out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_permute(const uint8_mat *mat, const sint32_mat *index, void *intrp_ptr, uint8_mat *out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_permute(const uint16_mat *mat, const sint32_mat *index, void *intrp_ptr, uint16_mat *out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_permute(const sint16_mat *mat, const sint32_mat *index, void *intrp_ptr, uint16_mat *out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_permute(const uint32_mat *mat, const sint32_mat *index, void *intrp_ptr, sint32_mat *out);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_permute(const sint32_mat *mat, const sint32_mat *index, void *intrp_ptr, sint32_mat *out);#
 *  \end{itemize}
 * @see matuniv_sort
 * @see matuniv_sort_index_partial
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_permute( const univ_mat *mat,
  const univ_mat *index, void *intrp_ptr, univ_mat *out );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matdbl_permute( const double_mat *mat,
  const sint32_mat *index, void *intrp_ptr, double_mat *out );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matflt_permute( const float_mat *mat,
  const sint32_mat *index, void *intrp_ptr, float_mat *out );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matu32_permute( const uint32_mat *mat,
  const sint32_mat *index, void *intrp_ptr, uint32_mat *out );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode mats32_permute( const sint32_mat *mat,
  const sint32_mat *index, void *intrp_ptr, sint32_mat *out );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matu16_permute( const uint16_mat *mat,
  const sint32_mat *index, void *intrp_ptr, uint16_mat *out );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode mats16_permute( const sint16_mat *mat,
  const sint32_mat *index, void *intrp_ptr, sint16_mat *out );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matu8_permute( const uint8_mat *mat,
  const sint32_mat *index, void *intrp_ptr, uint8_mat *out );


/** Sort a matrix.
 * Sort a matrix, treated as a flat array, by calling the
 * \Ref{matuniv_sort_index_partial} function and then permuting the
 * matrix with the \Ref{matuniv_permute} function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_sort.h
 * @source mat\_sort.c
 * @library matrix
 * @usage  #err_code = matuniv_sort(&mymat, intrp_ptr, &out);#
 * @return        Standard mutils error/OK code.
 * @param mat       Pointer to the matrix to sort.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param out       Pointer to returned sorted matrix (same size and type
 *    as input matrix).
 * @see matuniv_sort_index_partial
 * @see matuniv_permute
 * @see matuniv_unique
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_sort( const univ_mat *mat,
  void *intrp_ptr, univ_mat *out );


/** Table lookup.
 * Takes an input table and array of indices and uses the indices to do
 * table lookup: result\_data[i,:] = table\_data[index\_data[i],:].
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_sort.h
 * @source mat\_sort.c
 * @library matrix
 * @usage  #err_code = matuniv_table_lookup(&table, &index, intrp_ptr, &result);#
 * @return         Standard mutils error/OK code.
 * @param table The lookup table, also known as a codebook.  The type
 *     must be the same as the result type.  For example an 8-bit
 *     color palette is a table which would be stored as a 256 row,
 *     3 column matrix, where the color pixels would be indexed from
 *     0 to 255.
 * @param index The indices into the table or codebook.  index must
 *     be an integer type, ranging from 0 to MATUNIV\_NELEM(table)-1.
 * @param intrp_ptr Interrupt pointer.
 * @param result The looked up (decoded) values from the table.
 *     The data type of the result matrix must be the same as the table
 *     type.  The result matrix must have the same number of rows as
 *     there are index elements and columns as there are table columns.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_table_lookup(const double_mat *table, const univ_mat index, void *intrp_ptr, double_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_table_lookup(const float_mat *table, const univ_mat index, void *intrp_ptr, float_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_table_lookup(const uint8_mat *table, const univ_mat index, void *intrp_ptr, uint8_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_table_lookup(const uint16_mat *table, const univ_mat index, void *intrp_ptr, uint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_table_lookup(const sint16_mat *table, const univ_mat index, void *intrp_ptr, sint16_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_table_lookup(const uint32_mat *table, const univ_mat index, void *intrp_ptr, uint32_mat *result);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_table_lookup(const sint32_mat *table, const univ_mat index, void *intrp_ptr, sint32_mat *result);#
 * \end{itemize}
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_table_lookup( const univ_mat *table,
  const univ_mat *index, void *intrp_ptr, univ_mat *result );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matdbl_table_lookup( const double_mat *table,
  const univ_mat *index, void *intrp_ptr, double_mat *result );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matflt_table_lookup( const float_mat  *table,
  const univ_mat *index, void *intrp_ptr, float_mat *result );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu32_table_lookup( const uint32_mat *table,
  const univ_mat *index, void *intrp_ptr, uint32_mat *result );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats32_table_lookup( const sint32_mat *table,
  const univ_mat *index, void *intrp_ptr, sint32_mat *result );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu16_table_lookup( const uint16_mat *table,
  const univ_mat *index, void *intrp_ptr, uint16_mat *result );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats16_table_lookup( const sint16_mat *table,
  const univ_mat *index, void *intrp_ptr, sint16_mat *result );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu8_table_lookup( const uint8_mat *table,
  const univ_mat *index, void *intrp_ptr, uint8_mat *result );


/** Return the unique elements of a matrix.
 * This function takes a matrix and returns a vector matrix with
 * the unique elements of the input matrix, either sorted or unsorted.
 * If the returned elements are unsorted, then the elements are returned
 * in the order in which they first appear in the matrix array, otherwise they
 * are returned in ascending order.  If the routine fails, the output is not
 * allocated.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_sort.h
 * @source mat\_sort.c
 * @library matrix
 * @usage  #err_code = matuniv_unique(&in_mat, TRUE, intrp_ptr, &out_vector);#
 * @return              Standard mutils error/OK code.
 * @param mat           Pointer to the input matrix to determine the unique
 *                        elements.
 * @param sort          If TRUE uniq\_vector contains sorted elements,
 *                      if FALSE uniq\_vector contains unsorted elements.
 * @param intrp_ptr     Pointer for implementation of interrupt checking.
 * @param unique_vector Pointer to the returned matrix containing the unique
 *                      elements.  This matrix should only be declared and
 *                      must not be allocated. If the function call is
 *                      successful, it will be returned as a row vector
 *                      ( unique\_vector->nrow = 1 ) with the same type
 *                      as mat, otherwise it will not be allocated.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_unique(const uint8_mat *mat, boolean sort, void *intrp_ptr, uint8_mat *unique_vector);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_unique(const uint16_mat *mat, boolean sort, void *intrp_ptr, uint16_mat *unique_vector);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_unique(const uint32_mat *mat, boolean sort, void *intrp_ptr, uint32_mat *unique_vector);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_unique(const sint16_mat *mat, boolean sort, void *intrp_ptr, sint16_mat *unique_vector);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_unique(const sint32_mat *mat, boolean sort, void *intrp_ptr, sint32_mat *unique_vector);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_unique(const float_mat *mat, boolean sort, void *intrp_ptr, float_mat *unique_vector);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_unique(const double_mat *mat, boolean sort, void *intrp_ptr, double_mat *unique_vector);#
 *  \end{itemize}
 * @see matuniv_sort
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_unique( const univ_mat *mat,
  boolean sort, void *intrp_ptr, univ_mat *unique_vector );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matu8_unique( const uint8_mat *mat,
  boolean sort, void *intrp_ptr, uint8_mat *unique_vector );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matu16_unique( const uint16_mat *mat,
  boolean sort, void *intrp_ptr, uint16_mat *unique_vector );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matu32_unique( const uint32_mat *mat,
  boolean sort, void *intrp_ptr, uint32_mat *unique_vector );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode mats16_unique( const sint16_mat *mat,
  boolean sort, void *intrp_ptr, sint16_mat *unique_vector );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode mats32_unique( const sint32_mat *mat,
  boolean sort, void *intrp_ptr, sint32_mat *unique_vector );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matflt_unique( const float_mat *mat,
  boolean sort, void *intrp_ptr, float_mat *unique_vector );


/* Function documented with universal matrix version, above. */
MUTIL_LIBEXPORT mutil_errcode matdbl_unique( const double_mat *mat,
  boolean sort, void *intrp_ptr, double_mat *unique_vector );


#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_SORT_H_*/
