
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_umat.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_UMAT_H
#define IN_MAT_UMAT_H

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


/*
********************************
Macros for universal matrices
********************************
*/


/** Extract the number of elements of a universal matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_umat.h
 * @source mat\_umat.h
 * @library
 * @usage #n = MATUNIV_NELEM(&my_univ_mat);#
 * @return  The total number of elements (sint32) in the matrix.
 * @param  matrix   Pointer to a universal matrix.
 * @see _univ_mat
 * @see MATUNIV_NROW
 * @see MATUNIV_NCOL
 */
#define MATUNIV_NELEM(matrix)   (\
        (!matrix) ? ((sint32)MUTIL_INVALID_LENGTH) :\
        ((matrix)->type == MUTIL_UINT8)    ?  ((matrix)->mat.u8mat.nelem) :\
        ((matrix)->type == MUTIL_SINT8)    ?  ((matrix)->mat.s8mat.nelem) :\
        ((matrix)->type == MUTIL_UINT16)   ?  ((matrix)->mat.u16mat.nelem) :\
        ((matrix)->type == MUTIL_SINT16)   ?  ((matrix)->mat.s16mat.nelem) :\
        ((matrix)->type == MUTIL_UINT32)   ?  ((matrix)->mat.u32mat.nelem) :\
        ((matrix)->type == MUTIL_SINT32)   ?  ((matrix)->mat.s32mat.nelem) :\
        ((matrix)->type == MUTIL_FLOAT)    ?  ((matrix)->mat.fltmat.nelem) :\
        ((matrix)->type == MUTIL_DOUBLE)   ?  ((matrix)->mat.dblmat.nelem) :\
        ((matrix)->mat.cpxmat.nelem)\
    )


/** Extract the number of rows of a universal matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_umat.h
 * @source mat\_umat.h
 * @library
 * @usage #n = MATUNIV_NROW(&my_univ_mat);#
 * @return   The number of rows (sint32) in the matrix.
 * @param   matrix   Pointer to a universal matrix.
 * @see _univ_mat
 * @see MATUNIV_NCOL
 * @see MATUNIV_NELEM
 */
#define MATUNIV_NROW(matrix)    (\
        (!matrix) ? ((sint32) MUTIL_INVALID_LENGTH) :\
        ((matrix)->type == MUTIL_UINT8)    ?  ((matrix)->mat.u8mat.nrow) :\
        ((matrix)->type == MUTIL_SINT8)    ?  ((matrix)->mat.s8mat.nrow) :\
        ((matrix)->type == MUTIL_UINT16)   ?  ((matrix)->mat.u16mat.nrow) :\
        ((matrix)->type == MUTIL_SINT16)   ?  ((matrix)->mat.s16mat.nrow) :\
        ((matrix)->type == MUTIL_UINT32)   ?  ((matrix)->mat.u32mat.nrow) :\
        ((matrix)->type == MUTIL_SINT32)   ?  ((matrix)->mat.s32mat.nrow) :\
        ((matrix)->type == MUTIL_FLOAT)    ?  ((matrix)->mat.fltmat.nrow) :\
        ((matrix)->type == MUTIL_DOUBLE)   ?  ((matrix)->mat.dblmat.nrow) :\
        ((matrix)->mat.cpxmat.nrow)\
    )


/** Extract the number of columns of a universal matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_umat.h
 * @source mat\_umat.h
 * @library
 * @usage #n = MATUNIV_NCOL(&my_univ_mat);#
 * @return   The number of columns (sint32) in the matrix.
 * @param   matrix   Pointer to a universal matrix.
 * @see _univ_mat
 * @see MATUNIV_NROW
 * @see MATUNIV_NELEM
 */
#define MATUNIV_NCOL(matrix)    (\
        (!matrix) ? ((sint32)MUTIL_INVALID_LENGTH) :\
        ((matrix)->type == MUTIL_UINT8)    ?  ((matrix)->mat.u8mat.ncol) :\
        ((matrix)->type == MUTIL_SINT8)    ?  ((matrix)->mat.s8mat.ncol) :\
        ((matrix)->type == MUTIL_UINT16)   ?  ((matrix)->mat.u16mat.ncol) :\
        ((matrix)->type == MUTIL_SINT16)   ?  ((matrix)->mat.s16mat.ncol) :\
        ((matrix)->type == MUTIL_UINT32)   ?  ((matrix)->mat.u32mat.ncol) :\
        ((matrix)->type == MUTIL_SINT32)   ?  ((matrix)->mat.s32mat.ncol) :\
        ((matrix)->type == MUTIL_FLOAT)    ?  ((matrix)->mat.fltmat.ncol) :\
        ((matrix)->type == MUTIL_DOUBLE)   ?  ((matrix)->mat.dblmat.ncol) :\
        ((matrix)->mat.cpxmat.ncol)\
    )


/** Extract the data pointer in a universal matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_umat.h
 * @source mat\_umat.h
 * @library
 * @usage #myvoidptr = MATUNIV_DATA(&my_univ_mat);#
 * @return   A pointer to the data in the matrix, cast to void *.
 * @param  matrix   Pointer to a universal matrix.
 * @see _univ_mat
 */
#define MATUNIV_DATA(matrix)    (\
     (!matrix) ? (void *)(matrix) :\
     ((matrix)->type==MUTIL_UINT8)    ?  (void *)((matrix)->mat.u8mat.data) :\
     ((matrix)->type==MUTIL_SINT8)    ?  (void *)((matrix)->mat.s8mat.data) :\
     ((matrix)->type==MUTIL_UINT16)   ?  (void *)((matrix)->mat.u16mat.data) :\
     ((matrix)->type==MUTIL_SINT16)   ?  (void *)((matrix)->mat.s16mat.data) :\
     ((matrix)->type==MUTIL_UINT32)   ?  (void *)((matrix)->mat.u32mat.data) :\
     ((matrix)->type==MUTIL_SINT32)   ?  (void *)((matrix)->mat.s32mat.data) :\
     ((matrix)->type==MUTIL_FLOAT)    ?  (void *)((matrix)->mat.fltmat.data) :\
     ((matrix)->type==MUTIL_DOUBLE)   ?  (void *)((matrix)->mat.dblmat.data) :\
     (void *)((matrix)->mat.cpxmat.data)\
   )


/** Extract the pointer to the matrix in a universal matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_umat.h
 * @source mat\_umat.h
 * @library
 * @usage #myvoidptr = MATUNIV_MATRIX(&my_univ_mat);#
 * @return   A pointer to the matrix, cast to void *.
 * @param  matrix   Pointer to a universal matrix.
 * @see _univ_mat
 */
#define MATUNIV_MATRIX(matrix)    ( \
     (!matrix) ? (void *)(matrix) : \
     ((matrix)->type==MUTIL_UINT8)    ?  (void *) &((matrix)->mat.u8mat) : \
     ((matrix)->type==MUTIL_SINT8)    ?  (void *) &((matrix)->mat.s8mat) : \
     ((matrix)->type==MUTIL_UINT16)   ?  (void *) &((matrix)->mat.u16mat) : \
     ((matrix)->type==MUTIL_SINT16)   ?  (void *) &((matrix)->mat.s16mat) : \
     ((matrix)->type==MUTIL_UINT32)   ?  (void *) &((matrix)->mat.u32mat) : \
     ((matrix)->type==MUTIL_SINT32)   ?  (void *) &((matrix)->mat.s32mat) : \
     ((matrix)->type==MUTIL_FLOAT)    ?  (void *) &((matrix)->mat.fltmat) : \
     ((matrix)->type==MUTIL_DOUBLE)   ?  (void *) &((matrix)->mat.dblmat) : \
     (void *) &((matrix)->mat.cpxmat) \
   )


/** Extract an element from a universal matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_umat.h
 * @source mat\_umat.h
 * @library
 * @limits Does not work with complex matrices, and does not check that
 *    the indices are valid or that the matrix is non-null.
 * @usage #mynum = MATUNIV_ELEM(&my_univ_mat, i, j);#
 * @return   The value of the element of a universal matrix, cast to double.
 * @param   matrix   Pointer to a universal matrix.
 * @param   i        Row of element to extract.
 * @param   j        Column of element to extract.
 * @see MATUNIV_ELEM_ASSIGN
 * @see _univ_mat
 */
#define MATUNIV_ELEM(matrix, i, j)    (\
        ((matrix)->type==MUTIL_UINT8)    ? \
           (double)((matrix)->mat.u8mat.data[ \
             MATANY_INDEX(&(matrix)->mat.u8mat, (i), (j))]) :\
        ((matrix)->type==MUTIL_SINT8)    ? \
           (double)((matrix)->mat.s8mat.data[ \
             MATANY_INDEX(&(matrix)->mat.s8mat, (i), (j))]) :\
        ((matrix)->type==MUTIL_UINT16)   ? \
           (double)((matrix)->mat.u16mat.data[ \
             MATANY_INDEX(&(matrix)->mat.u16mat, (i), (j))]) :\
        ((matrix)->type==MUTIL_SINT16)   ? \
           (double)((matrix)->mat.s16mat.data[ \
             MATANY_INDEX(&(matrix)->mat.s16mat, (i), (j))]) :\
        ((matrix)->type==MUTIL_UINT32)   ? \
           (double)((matrix)->mat.u32mat.data[ \
             MATANY_INDEX(&(matrix)->mat.u32mat, (i), (j))]) :\
        ((matrix)->type==MUTIL_SINT32)   ? \
           (double)((matrix)->mat.s32mat.data[ \
             MATANY_INDEX(&(matrix)->mat.s32mat, (i), (j))]) :\
        ((matrix)->type==MUTIL_FLOAT)    ? \
           (double)((matrix)->mat.fltmat.data[ \
             MATANY_INDEX(&(matrix)->mat.fltmat, (i), (j))]) :\
           (double)((matrix)->mat.dblmat.data[ \
             MATANY_INDEX(&(matrix)->mat.dblmat, (i), (j))]) \
    )


/** Assign to an element in a universal matrix.
 *  Casts the value input to fit the appropriate data type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_umat.h
 * @source mat\_umat.h
 * @library
 * @limits Does not work with complex matrices, and does not check that
 *    the indices are valid, that the matrix is non-null, or that data
 *    types are compatible for assignment.
 * @param  matrix  Pointer to a universal matrix.
 * @param  i       Row of element to assign to.
 * @param  j       Column of element to assign to.
 * @param  val     Value to assign
 * @usage #MATUNIV_ELEM_ASSIGN(&my_univ_mat, i, j, val);#
 * @see MATUNIV_ELEM
 * @see _univ_mat
 */
#define MATUNIV_ELEM_ASSIGN(matrix, i, j, val)  \
        switch ((matrix)->type) {\
          case MUTIL_UINT8: \
            (matrix)->mat.u8mat.data[ \
              MATANY_INDEX(&(matrix)->mat.u8mat, (i), (j))] = (uint8)(val);\
            break;\
          case MUTIL_SINT8: \
            (matrix)->mat.s8mat.data[ \
              MATANY_INDEX(&(matrix)->mat.s8mat, (i), (j))] = (sint8)(val);\
            break;\
          case MUTIL_UINT16: \
            (matrix)->mat.u16mat.data[ \
               MATANY_INDEX(&(matrix)->mat.u16mat, (i), (j))] = (uint16)(val);\
             break;\
          case MUTIL_SINT16: \
            (matrix)->mat.s16mat.data[ \
              MATANY_INDEX(&(matrix)->mat.s16mat, (i), (j))] = (sint16)(val);\
            break;\
          case MUTIL_UINT32: \
            (matrix)->mat.u32mat.data[ \
              MATANY_INDEX(&(matrix)->mat.u32mat, (i), (j))] = (uint32)(val);\
            break;\
          case MUTIL_SINT32: \
            (matrix)->mat.s32mat.data[ \
              MATANY_INDEX(&(matrix)->mat.s32mat, (i), (j))] = (sint32)(val);\
            break;\
          case MUTIL_FLOAT: \
            (matrix)->mat.fltmat.data[ \
              MATANY_INDEX(&(matrix)->mat.fltmat, (i), (j))] = (float)(val);\
            break;\
          case MUTIL_DOUBLE:\
          default: \
            (matrix)->mat.dblmat.data[ \
              MATANY_INDEX(&(matrix)->mat.dblmat, (i), (j))] = (double)(val);\
            break;\
        }


/** Check if the data types of two universal matrices are the same.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_umat.h
 * @source mat\_umat.h
 * @library
 * @usage #ok = MATUNIV_CHECK_TYPE(mat1_ptr, mat2_ptr);#
 * @return   TRUE (same types) or FALSE (not the same or NULL pointers).
 * @param  mat1  Pointer to the first universal matrix.
 * @param  mat2  Pointer to the second universal matrix.
 * @see _mutil_data_type
 * @see _univ_mat
 */
#define MATUNIV_CHECK_TYPE( mat1, mat2 )  \
        ( mat1 && mat2 && ( (mat1)->type == (mat2)->type) )


/** Extract an element from a universal matrix given the data index.
 * Casts the value returned to double.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_umat.h
 * @source mat\_umat.h
 * @library
 * @limits Does not work with complex matrices, and does not check that
 *    the indices are valid or that the matrix is non-null.
 * @usage #mynum = MATUNIV_ELEM_BY_INDEX(&my_univ_mat, i);#
 * @return   The value of the element of a universal matrix, cast to double.
 * @param   matrix   Pointer to a universal matrix.
 * @param   i        Index of element to extract from the data array.
 * @see MATUNIV_ELEM
 * @see MATUNIV_ELEM_ASSIGN
 * @see MATUNIV_ELEM_ASSIGN_BY_INDEX
 * @see _univ_mat
 */
#define MATUNIV_ELEM_BY_INDEX( matrix, i )    (\
        ((matrix)->type==MUTIL_UINT8)    ? \
           (double)((matrix)->mat.u8mat.data[ (i) ]) :\
        ((matrix)->type==MUTIL_SINT8)    ? \
           (double)((matrix)->mat.s8mat.data[ (i) ]) :\
        ((matrix)->type==MUTIL_UINT16)   ? \
           (double)((matrix)->mat.u16mat.data[ (i) ]) :\
        ((matrix)->type==MUTIL_SINT16)   ? \
           (double)((matrix)->mat.s16mat.data[ (i) ]) :\
        ((matrix)->type==MUTIL_UINT32)   ? \
           (double)((matrix)->mat.u32mat.data[ (i) ]) :\
        ((matrix)->type==MUTIL_SINT32)   ? \
           (double)((matrix)->mat.s32mat.data[ (i) ]) :\
        ((matrix)->type==MUTIL_FLOAT)    ? \
           (double)((matrix)->mat.fltmat.data[ (i) ]) :\
           (double)((matrix)->mat.dblmat.data[ (i) ]) \
    )


/** Assign to an element in a universal matrix given the data index.
 * Casts the value input to fit the appropriate data type.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_umat.h
 * @source mat\_umat.h
 * @library
 * @limits Does not work with complex matrices, and does not check that
 *    the indices are valid, that the matrix is non-null, or that data
 *    types are compatible for assignment.
 * @param  matrix  Pointer to a universal matrix.
 * @param  i       Index of element in the data array to assign to.
 * @param  val     Value to assign
 * @usage #MATUNIV_ELEM_ASSIGN_BY_INDEX(&my_univ_mat, i, val);#
 * @see MATUNIV_ELEM
 * @see MATUNIV_ELEM_ASSIGN
 * @see MATUNIV_ELEM_BY_INDEX
 * @see _univ_mat
 */
#define MATUNIV_ELEM_ASSIGN_BY_INDEX( matrix, i, val )\
        switch ((matrix)->type) {\
          case MUTIL_UINT8: \
            (matrix)->mat.u8mat.data[ (i) ] = (uint8)(val);\
            break;\
          case MUTIL_SINT8: \
            (matrix)->mat.s8mat.data[ (i) ] = (sint8)(val);\
            break;\
          case MUTIL_UINT16: \
            (matrix)->mat.u16mat.data[ (i) ] = (uint16)(val);\
             break;\
          case MUTIL_SINT16: \
            (matrix)->mat.s16mat.data[ (i) ] = (sint16)(val);\
            break;\
          case MUTIL_UINT32: \
            (matrix)->mat.u32mat.data[ (i) ] = (uint32)(val);\
            break;\
          case MUTIL_SINT32: \
            (matrix)->mat.s32mat.data[ (i) ] = (sint32)(val);\
            break;\
          case MUTIL_FLOAT: \
            (matrix)->mat.fltmat.data[ (i) ] = (float)(val);\
            break;\
          case MUTIL_DOUBLE:\
          default: \
            (matrix)->mat.dblmat.data[ (i) ] = (double)(val);\
            break;\
        }


#ifdef __cplusplus
}
#endif

#endif /* ifdef IN_MAT_UMAT_H */
