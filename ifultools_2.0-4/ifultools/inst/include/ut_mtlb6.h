
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_mtlb6.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

/*
 *********************************************************************
 *                                                                   *
 * File: ut_mtlb6.h                                                  *
 *                                                                   *
 * Main header file for the Matlab/MUtils wrapper library functions. *
 *                                                                   *
 * The functions declared within are used for converting Matlab data *
 * types to MUTILS data types, and vice versa. These functions may   *
 * be used to create MEX functions to be called from the Matlab      *
 * command prompt or in stand-alone executables to be run from the   *
 * DOS/Unix command prompt.                                          *
 *                                                                   *
 * NOTE: These wrapper functions have been tested for Matlab         *
 *       version 6.5, under the Windows XP OS.                       *
 *                                                                   *
 *********************************************************************
*/

#ifndef _UT_MTLB6_H
#define _UT_MTLB6_H

/* MUTILS HEADERS */
#include "cls_type.h"
#include "img_ls2.h"
#include "img_type.h"
#include "mat_iter.h"
#include "mat_type.h"
#include "pgn_type.h"
#include "ut_debug.h"
#include "ut_type.h"
#include "vid_type.h"
#include "wav_type.h"

/* MATLAB MEX API HEADER */
#include "mex.h"



/*
 ******************
 *                *
 * LIBRARY MACROS *
 *                *
 ******************
 */

/** Macro to convert MUtils matrix flat index to Matlab matrix flat index.
 * Matrices in MUtils are stored in a flat array that is row-major. Matlab
 * uses the FORTRAN convention storing matrices in column-major flat arrays.
 * This macro converts a MUtils flat index to a Matlab flat index.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.h
 */
#define MATRIX_INDEX_MUTILS_TO_MATLAB( N, ROW, COL ) \
  ( (N / COL) + (N % COL) * ROW )

/** Macro to convert Matlab matrix flat index to MUtils matrix flat index.
 * Matrices in Matlab are stored in a flat array that is column-major. MUtils
 * uses the flat array row-major storing convention.
 * This macro converts a Matlab flat index to a MUtils flat index.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.h
 */
#define MATRIX_INDEX_MATLAB_TO_MUTILS( N, ROW, COL ) \
  MATRIX_INDEX_MUTILS_TO_MATLAB( N, COL, ROW )


/*
 *
 * An enumeration of the MUTILS data types that will be supported
 * by automated mex source generation.
 *
 */
enum _mm_data_type {

  /* Void type */
  MM_VOID,

  /* Scalar types */
  MM_BOOLEAN,
  MM_DCOMPLEX,
  MM_DOUBLE,
  MM_FLOAT,
  MM_SINT32,
  MM_SINT16,
  MM_SINT8,
  MM_UINT32,
  MM_UINT16,
  MM_UINT8,
  MM_UNIV_SCALAR,

  /* General MUTILS types */
  MM_MUTIL_BEZIER_CURVE,
  MM_MUTIL_BEZIER_CURVE_LIST,
  MM_MUTIL_BOUNDARY_TYPE,
  MM_MUTIL_DATA_TYPE,
  MM_MUTIL_ERRCODE,
  MM_MUTIL_INTERPOLATION_TYPE,
  MM_MUTIL_ITERATOR,
  MM_MUTIL_ITERATOR_DIM,
  MM_MUTIL_MSG_LEVEL,
  MM_MUTIL_VECTOR,

  /* Classify types */
  MM_CLS_SEARCH_TYPE,
  MM_TSVQ_NODE,

  /* Image types */
  MM_IMG_CONNECTED_TYPE,
  MM_IMG_COOCCURRENCE_FEATURE_TYPE,
  MM_IMG_DIFFUSION_MODEL_TYPE,
  MM_IMG_DISTANCE_TYPE,
  MM_IMG_LAWS_FILTER_TYPE,
  MM_IMG_LSET_PARAMS,
  MM_IMG_MODEL_TYPE,
  MM_IMG_ORIENTATION_TYPE,
  MM_IMG_SNAKE_PARAMS,
  MM_IMG_THRESH_TYPE,
  
  /* Matrix types */
  MM_MAT_SET,
  MM_UNIV_MAT,
  MM_DCOMPLEX_MAT,
  MM_DOUBLE_MAT,
  MM_FLOAT_MAT,
  MM_SINT32_MAT,
  MM_SINT16_MAT,
  MM_SINT8_MAT,
  MM_UINT32_MAT,
  MM_UINT16_MAT,
  MM_UINT8_MAT,

  /* Polygon types */
  MM_PGN_CHAIN,
  MM_PGN_POINT_CONTAIN_TYPE,

  /* Video types */
  MM_VID_GRADIENT_TYPE,

  /* Wavelet-related types */
  MM_WAV_FDP_ESTIMATOR,
  MM_WAV_FILTER_TYPE,
  MM_WAV_HMT_MODEL,
  MM_WAV_PARAMS_GABOR_2D,
  MM_WAV_SUBBAND_DWT_2D,
  MM_WAV_TRANSFORM
};

/* Type definition of the above enum */
typedef enum _mm_data_type mm_data_type;


#ifdef __cplusplus
extern "C" {
#endif

/*
 *******************************************************************
 *                                                                 *
 * FUNCTIONS CONSTRUCTING MATLAB DATA TYPES FROM MUTILS DATA TYPES *
 *                                                                 *
 *******************************************************************
 */

/*
 ****************
 * MATRIX TYPES *
 ****************
 */

/** Convert an MUTILS matrix set to a Matlab cell array.
 * This function converts the MUTILS matrix set data structure to a
 * Matlab cell array. The dimensions of the resulting cell array, and
 * the elements of that array, will be identical to those of the matrix
 * set and its internal universal matrices, respectively. The class of
 * of the cell array matrices is determined by the type of the matrix
 * set universal matrices. See \Ref{matuniv_to_mxArray} for Matlab/MUTILS
 * class/type correspondence.
 *
 * NOTE: In the case where the matrix set is one-dimensional with N elements
 * the resulting Matlab cell array will be two-dimensional with dimensions
 * 1 x N. This is because, by convention, all Matlab cell arrays have at least
 * two dimensions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.c
 * @library wrap_mtl
 * @usage #err = matset_to_mxArray( mset_ptr, ptr_array_ptr );#
 * @return Standard mutils error/OK code.
 * @param mset_ptr  Pointer to the matrix set to be converted.
 * @param ptr_array_ptr Pointer to a Matlab mxArray pointer.
 */
MUTIL_WRAPEXPORT mutil_errcode
  matset_to_mxArray( const mat_set* mset_ptr, mxArray** ptr_array_ptr );


/** Convert an MUTILS universal matrix to a Matlab matrix.
 * This function converts an MUTILS universal matrix to a
 * Matlab matrix of the same size and data type. The data type of the
 * universal matrix determines the class of the resulting Matlab matrix.
 * The correspondence between Matlab data classes and MUTILS data types
 * is as follows:
 * 
 * /begin{tabbing}{2}
 *   'double' (complexity = mxCOMPLEX) & MUTIL_DCOMPLEX \\
 *   'double' (complexity = mxREAL) & MUTIL_DOUBLE \\
 *   'single' & MUTIL_FLOAT \\
 *   'int32' & MUTIL_SINT32 \\
 *   'int16' & MUTIL_SINT16 \\
 *   'int8' & MUTIL_SINT8 \\
 *   'uint32' & MUTIL_SINT32 \\
 *   'uint16' & MUTIL_SINT16 \\
 *   'uint8' & MUTIL_SINT8 \\
 *  \end{tabbing}
 *
 * Following
 * the MUTILS convention, universal matrices of type MUTIL_DCOMPLEX
 * will be converted to floating
 * point complex Matlab matrices with their complexity flag set to
 * mxCOMPLEX and data class 'double.' 
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.c
 * @library wrap_mtl
 * @usage #err = matuniv_to_mxArray( mset_ptr, ptr_array_ptr );#
 * @return Standard mutils error/OK code.
 * @param mat_ptr  Pointer to the universal matrix to be converted.
 * @param ptr_array_ptr Pointer to a Matlab mxArray pointer.
 * @same \begin{itemize}
 *  \item #MUTIL_WRAPEXPORT mutil_errcode matcpx_to_mxArray( const dcomplex_mat* mat_ptr, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode matdbl_to_mxArray( const double_mat* mat_ptr, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode matflt_to_mxArray( const float_mat* mat_ptr, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mats32_to_mxArray( const sint32_mat* mat_ptr, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mats16_to_mxArray( const sint16_mat* mat_ptr, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode matu32_to_mxArray( const uint32_mat* mat_ptr, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode matu16_to_mxArray( const uint16_mat* mat_ptr, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode matu8_to_mxArray( const uint8_mat* mat_ptr, mxArray** ptr_array_ptr );#
 *  \end{itemize}
 */
MUTIL_WRAPEXPORT mutil_errcode
  matuniv_to_mxArray( const univ_mat* mat_ptr, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  matcpx_to_mxArray( const dcomplex_mat* mat_ptr, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  matdbl_to_mxArray( const double_mat* mat_ptr, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  matflt_to_mxArray( const float_mat* mat_ptr, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mats32_to_mxArray( const sint32_mat* mat_ptr, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mats16_to_mxArray( const sint16_mat* mat_ptr, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  matu32_to_mxArray( const uint32_mat* mat_ptr, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  matu16_to_mxArray( const uint16_mat* mat_ptr, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  matu8_to_mxArray( const uint8_mat* mat_ptr, mxArray** ptr_array_ptr );

/*
 ****************
 * SCALAR TYPES *
 ****************
 */

/** Convert an MUTILS universal scalar to a Matlab scalar.
 * This function converts an MUTILS universal scalar to a
 * Matlab scalar of the same type. The data type of the universal
 * scalar determines the class of the resulting Matlab scalar. Following
 * the MUTILS convention, universal scalars of type MUTIL_DCOMPLEX
 * will be converted to floating
 * point complex Matlab scalars of class 'double.' See
 * \Ref{matuniv_to_mxArray} for Matlab/MUTILS class/type correspondence.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.c
 * @library wrap_mtl
 * @usage #err = scauniv_to_mxArray( scalar, ptr_array_ptr );#
 * @return Standard mutils error/OK code.
 * @param scalar  Universal scalar to be converted.
 * @param ptr_array_ptr Pointer to a Matlab mxArray pointer.
 * @same \begin{itemize}
 *  \item #MUTIL_WRAPEXPORT mutil_errcode dcomplex_to_mxArray( const dcomplex scalar, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode double_to_mxArray( const double scalar, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode float_to_mxArray( const float scalar, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode sint32_to_mxArray( const sint32 scalar, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode sint16_to_mxArray( const sint16 scalar, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode sint8_to_mxArray( const sint8 scalar, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode uint32_to_mxArray( const uint32 scalar, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode uint16_to_mxArray( const uint16 scalar, mxArray** ptr_array_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode uint8_to_mxArray( const uint8 scalar, mxArray** ptr_array_ptr );#
 *  \end{itemize}
 */
MUTIL_WRAPEXPORT mutil_errcode
  scauniv_to_mxArray( const univ_scalar scalar, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  dcomplex_to_mxArray( const dcomplex scalar, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  double_to_mxArray( const double scalar, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  float_to_mxArray( const float scalar, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  sint32_to_mxArray( const sint32 scalar, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  sint16_to_mxArray( const sint16 scalar, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  sint8_to_mxArray( const sint8 scalar, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  uint32_to_mxArray( const uint32 scalar, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  uint16_to_mxArray( const uint16 scalar, mxArray** ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  uint8_to_mxArray( const uint8 scalar, mxArray** ptr_array_ptr );


/** Convert an MUTILS boolean scalar to a Matlab logical scalar.
 * This function converts an MUTILS boolean scalar to a
 * Matlab logical scalar of numerical value 1 or 0. If the MUTILS
 * boolean scalar is non-zero, the Matlab logical scalar will be set to
 * 1; otherwise, the Matlab logical scalar will be set to zero.
 *
 * This function simply creates a double floating point Matlab scalar of
 * value 1 or 0, and sets its logical flag.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.c
 * @library wrap_mtl
 * @usage #err = boolean_to_mxArray( scalar, ptr_array_ptr );#
 * @return Standard mutils error/OK code.
 * @param scalar  Boolean scalar to be converted.
 * @param ptr_array_ptr Pointer to a Matlab mxArray pointer.
 */
MUTIL_WRAPEXPORT mutil_errcode
  boolean_to_mxArray( const boolean scalar, mxArray** ptr_array_ptr );

/*
 ********************
 * ENUMERATED TYPES *
 ********************
 */

/* General enumerated types */

MUTIL_WRAPEXPORT mutil_errcode mutil_boundary_type_to_mxArray(
  const mutil_boundary_type   enum_val,
  mxArray**                   ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode mutil_data_type_to_mxArray(
  const mutil_data_type   enum_val,
  mxArray**               ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode mutil_interpolation_type_to_mxArray(
  const mutil_interpolation_type   enum_val,
  mxArray**                        ptr_array_ptr );

/* Classify enumerated types */

MUTIL_WRAPEXPORT mutil_errcode cls_search_type_to_mxArray(
  const cls_search_type   enum_val,
  mxArray**               ptr_array_ptr );

/* Image enumerated types */

MUTIL_WRAPEXPORT mutil_errcode img_connected_type_to_mxArray(
  const img_connected_type   enum_val,
  mxArray**                  ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode img_distance_type_to_mxArray(
  const img_distance_type   enum_val,
  mxArray**                 ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode img_thresh_type_to_mxArray(
  const img_thresh_type   enum_val,
  mxArray**               ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode img_model_type_to_mxArray(
  const img_model_type   enum_val,
  mxArray**              ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode img_diffusion_model_type_to_mxArray(
  const img_diffusion_model_type   enum_ptr,
  mxArray**                        ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode img_orientation_type_to_mxArray(
  const img_orientation_type   enum_val,
  mxArray**                    ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode img_cooccurrence_feature_type_to_mxArray(
  const img_cooccurrence_feature_type   enum_val,
  mxArray**                             ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode img_laws_filter_type_to_mxArray(
  const img_laws_filter_type   enum_val,
  mxArray**                    ptr_array_ptr );


/* Polygon enumerated types */

MUTIL_WRAPEXPORT mutil_errcode pgn_point_contain_type_to_mxArray(
  const pgn_point_contain_type   enum_val,
  mxArray**                      ptr_array_ptr );


/* Video enumerated types */

MUTIL_WRAPEXPORT mutil_errcode vid_gradient_type_to_mxArray(
  const vid_gradient_type   enum_val,
  mxArray**                 ptr_array_ptr );


/* Wavelet enumerated types */

MUTIL_WRAPEXPORT mutil_errcode wav_fdp_estimator_to_mxArray(
  const wav_fdp_estimator   enum_val,
  mxArray**                 ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode wav_filter_type_to_mxArray(
  const wav_filter_type   enum_val,
  mxArray**               ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode wav_subband_dwt_2d_to_mxArray(
  const wav_subband_dwt_2d   enum_val,
  mxArray**                  ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode wav_transform_to_mxArray(
  const wav_transform   enum_val,
  mxArray**             ptr_array_ptr );

/* 
 *******************
 * STRUCTURE TYPES *
 *******************
 */

/* General structure types */
MUTIL_WRAPEXPORT mutil_errcode mutil_bezier_curve_to_mxArray(
  const mutil_bezier_curve*   struct_ptr,
  mxArray**                   ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode mutil_bezier_curve_list_to_mxArray(
  const mutil_bezier_curve_list*   struct_ptr,
  mxArray**                        ptr_array_ptr );


/* Classify structure types */
MUTIL_WRAPEXPORT mutil_errcode cls_tsvq_node_to_mxArray(
  const cls_tsvq_node*   struct_ptr,
  mxArray**              ptr_array_ptr );


/* Image structure types */
MUTIL_WRAPEXPORT mutil_errcode img_snake_params_to_mxArray(
  const img_snake_params*   struct_ptr,
  mxArray**                 ptr_array_ptr );


/* Polygon structure types */
MUTIL_WRAPEXPORT mutil_errcode pgn_chain_to_mxArray(
  const pgn_chain*   struct_ptr,
  mxArray**          ptr_array_ptr );


/* Wavelet structure types */
MUTIL_WRAPEXPORT mutil_errcode wav_hmt_model_to_mxArray(
  const wav_hmt_model*   struct_ptr,
  mxArray**              ptr_array_ptr );

MUTIL_WRAPEXPORT mutil_errcode wav_params_gabor_2d_to_mxArray(
  const wav_params_gabor_2d*   struct_ptr,
  mxArray**              ptr_array_ptr );


/*
 *******************************************************************
 *                                                                 *
 * FUNCTIONS CONSTRUCTING MUTILS DATA TYPES FROM MATLAB DATA TYPES *
 *                                                                 *
 *******************************************************************
 */

/*
 ****************
 * MATRIX TYPES *
 ****************
 */

/** Convert a Matlab cell array to a MUTILS matrix set.
 * This function converts a Matlab cell array to a MUTILS matrix set.
 * The dimensions of the resulting matrix set, and
 * the internal universal matrices, will be identical to those of the
 * cell array and its elements, respectively.
 *
 * To comply with the MUTILS convention for matrix sets, it is required
 * the elements of the cell array be of the same data class. See
 * \Ref{matuniv_to_mxArray} for Matlab/MUTILS
 * class/type correspondence.
 *
 * NOTE: In the case where the cell array dimension is of the form 1 x N
 * or N x 1 (considered to be two-dimensional in Matlab), the resulting
 * matrix set will be one-dimensional with N elements.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.c
 * @library wrap_mtl
 * @usage #err = mxArray_to_matset( array_ptr, mset );#
 * @return Standard mutils error/OK code.
 * @param array_ptr Pointer to the Matlab mxArray to be converted.
 * @param mset_ptr  Pointer to an unallocated matrix set.
 */
MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_matset( const mxArray* array_ptr, mat_set* mset_ptr );


/** Convert a Matlab matrix to a MUTILS universal matrix.
 * This function converts a Matlab matrix to a MUTILS universal matrix.
 * The type of the resulting universal matrix is determined by the class
 * of the Matlab matrix. See \Ref{matuniv_to_mxArray} for Matlab/MUTILS
 * class/type correspondence.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.c
 * @library wrap_mtl
 * @usage #err = mxArray_to_matuniv( array_ptr, mat_ptr );#
 * @return Standard mutils error/OK code.
 * @param array_ptr Pointer to the Matlab mxArray to be converted.
 * @param mat_ptr  Pointer to an unallocated universal matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_matcpx( const mxArray* array_ptr, dcomplex_mat* mat_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_matdbl( const mxArray* array_ptr, double_mat* mat_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_matflt( const mxArray* array_ptr, float_mat* mat_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_mats32( const mxArray* array_ptr, sint32_mat* mat_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_mats16( const mxArray* array_ptr, sint16_mat* mat_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_matu32( const mxArray* array_ptr, uint32_mat* mat_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_matu16( const mxArray* array_ptr, uint16_mat* mat_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_matu8( const mxArray* array_ptr, uint8_mat* mat_ptr );#
 *  \end{itemize}
 */
MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_matuniv( const mxArray* array_ptr, univ_mat* mat_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_matcpx( const mxArray* array_ptr, dcomplex_mat* mat_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_matdbl( const mxArray* array_ptr, double_mat* mat_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_matflt( const mxArray* array_ptr, float_mat* mat_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_mats32( const mxArray* array_ptr, sint32_mat* mat_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_mats16( const mxArray* array_ptr, sint16_mat* mat_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_matu32( const mxArray* array_ptr, uint32_mat* mat_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_matu16( const mxArray* array_ptr, uint16_mat* mat_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_matu8( const mxArray* array_ptr, uint8_mat* mat_ptr );

/*
 ****************
 * SCALAR TYPES *
 ****************
 */

/** Convert a Matlab scalar to a MUTILS universal scalar.
 * This function converts a Matlab scalar to a MUTILS universal scalar.
 * The type of the resulting universal scalar is determined by the class
 * of the Matlab scalar. See \Ref{matuniv_to_mxArray} for Matlab/MUTILS
 * class/type correspondence.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.c
 * @library wrap_mtl
 * @usage #err = mxArray_to_scauniv( array_ptr, sca_ptr );#
 * @return Standard mutils error/OK code.
 * @param array_ptr Pointer to the Matlab mxArray to be converted.
 * @param mat_ptr  Pointer to an universal scalar.
 * @same \begin{itemize}
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_dcomplex( const mxArray* array_ptr, dcomplex* sca_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_double( const mxArray* array_ptr, double* sca_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_float( const mxArray* array_ptr, float* sca_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_sint32( const mxArray* array_ptr, sint32* sca_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_sint16( const mxArray* array_ptr, sint16* sca_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_sint8( const mxArray* array_ptr, sint8* sca_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_uint32( const mxArray* array_ptr, uint32* sca_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_uint16( const mxArray* array_ptr, uint16* sca_ptr );#
 *  \item #MUTIL_WRAPEXPORT mutil_errcode mxArray_to_uint8( const mxArray* array_ptr, uint8* sca_ptr );#
 *  \end{itemize}
 */
MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_scauniv( const mxArray* array_ptr, univ_scalar* sca_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_dcomplex( const mxArray* array_ptr, dcomplex* sca_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_double( const mxArray* array_ptr, double* sca_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_float( const mxArray* array_ptr, float* sca_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_sint32( const mxArray* array_ptr, sint32* sca_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_sint16( const mxArray* array_ptr, sint16* sca_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_sint8( const mxArray* array_ptr, sint8* sca_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_uint32( const mxArray* array_ptr, uint32* sca_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_uint16( const mxArray* array_ptr, uint16* sca_ptr );

MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_uint8( const mxArray* array_ptr, uint8* sca_ptr );


/** Convert a Matlab logical scalar to a MUTILS boolean scalar.
 * This function converts a Matlab logical scalar to a MUTILS boolean scalar
 * of numerical value 1 or 0. The Matlab scalar must have its logical flag
 * set.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.c
 * @library wrap_mtl
 * @usage #err = mxArray_to_boolean( scalar, ptr_array_ptr );#
 * @return Standard mutils error/OK code.
 * @param array_ptr Pointer to the Matlab logical mxArray to be converted.
 * @param scalar  Pointer to a boolean scalar.
 */
MUTIL_WRAPEXPORT mutil_errcode
  mxArray_to_boolean( const mxArray* array_ptr, boolean* sca_ptr );

/*
 ********************
 * ENUMERATED TYPES *
 ********************
 */

/* General enumerated types */

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_mutil_boundary_type(
  const mxArray*         array_ptr,
  mutil_boundary_type*   enum_ptr );


MUTIL_WRAPEXPORT mutil_errcode mxArray_to_mutil_data_type(
  const mxArray*     array_ptr,
  mutil_data_type*   enum_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_mutil_interpolation_type(
  const mxArray*              array_ptr,
  mutil_interpolation_type*   enum_ptr );

/* Classify enumerated types */

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_cls_search_type(
  const mxArray*     array_ptr,
  cls_search_type*   enum_ptr );

/* Image enumerated types */

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_img_connected_type(
  const mxArray*        array_ptr,
  img_connected_type*   enum_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_img_distance_type(
  const mxArray*       array_ptr,
  img_distance_type*   enum_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_img_thresh_type(
  const mxArray*     array_ptr,
  img_thresh_type*   enum_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_img_model_type(
  const mxArray*    array_ptr,
  img_model_type*   enum_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_img_diffusion_model_type(
  const mxArray*              array_ptr,
  img_diffusion_model_type*   enum_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_img_orientation_type(
  const mxArray*          array_ptr,
  img_orientation_type*   enum_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_img_cooccurrence_feature_type(
  const mxArray*                   array_ptr,
  img_cooccurrence_feature_type*   enum_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_img_laws_filter_type(
  const mxArray*          array_ptr,
  img_laws_filter_type*   enum_ptr );

/* Polygon enumerated types */

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_pgn_point_contain_type(
  const mxArray*            array_ptr,
  pgn_point_contain_type*   enum_ptr );

/* Video enumerated types */

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_vid_gradient_type(
  const mxArray*       array_ptr,
  vid_gradient_type*   enum_ptr );

/* Wavelet enumerated types */

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_wav_fdp_estimator(
  const mxArray*       array_ptr,
  wav_fdp_estimator*   enum_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_wav_filter_type(
  const mxArray*     array_ptr,
  wav_filter_type*   enum_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_wav_subband_dwt_2d(
  const mxArray*        array_ptr,
  wav_subband_dwt_2d*   enum_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_wav_transform(
  const mxArray*   array_ptr,
  wav_transform*   enum_ptr );

/*
 *******************
 * STRUCTURE TYPES *
 *******************
 */

/* General structure types */
MUTIL_WRAPEXPORT mutil_errcode mxArray_to_mutil_bezier_curve(
    const mxArray*        array_ptr,
    mutil_bezier_curve*   struct_ptr );

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_mutil_bezier_curve_list(
    const mxArray*             array_ptr,
    mutil_bezier_curve_list*   struct_ptr );

/* Classify structure types */
MUTIL_WRAPEXPORT mutil_errcode mxArray_to_cls_tsvq_node(
    const mxArray*   array_ptr,
    cls_tsvq_node*   struct_ptr );

/* Image structure types */
MUTIL_WRAPEXPORT mutil_errcode mxArray_to_img_snake_params(
    const mxArray*      ptr_array_ptr,
    img_snake_params*   struct_ptr );

/* Polygon structure types */
MUTIL_WRAPEXPORT mutil_errcode mxArray_to_pgn_chain(
    const mxArray*   array_ptr,
    pgn_chain*       struct_ptr );


/* Wavelet structure types */

/** Convert a Matlab structure array to a MUTILS wav_hmt_model structure.
 * This function take a Matlab structure mxArray and creates an MUTILS
 * wav_hmt_model structure. The fields names in the Matlab structure,
 * and the associated names, must be in 1-to-1 correspondence with the
 * member names of the MUTILS structure.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.c
 * @library wrap_mtl
 * @usage #err = mxArray_to_wav_hmt_model( scalar, ptr_array_ptr );#
 * @return Standard mutils error/OK code.
 * @param array_ptr Pointer to the Matlab structure array.
 * @param struct_ptr Pointer to a wav_hmt_model structure.
 * @see _wav_hmt_model
 */
MUTIL_WRAPEXPORT mutil_errcode mxArray_to_wav_hmt_model(
  const mxArray*   array_ptr,
  wav_hmt_model*   struct_ptr );

/** Convert a Matlab structure array to a MUTILS wav_params_gabor_2d structure.
 * This function take a Matlab structure mxArray and creates an MUTILS
 * wav_params_gabor_2d structure. The fields names in the Matlab structure,
 * and the associated names, must be in 1-to-1 correspondence with the
 * member names of the MUTILS structure.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut_mtlb6.h
 * @source ut_mtlb6.c
 * @library wrap_mtl
 * @usage #err = mxArray_to_wav_params_gabor_2d( scalar, ptr_array_ptr );#
 * @return Standard mutils error/OK code.
 * @param array_ptr Pointer to the Matlab structure array.
 * @param struct_ptr Pointer to a wav_params_gabor_2d structure.
 * @see _wav_params_gabor_2d
 */

MUTIL_WRAPEXPORT mutil_errcode mxArray_to_wav_params_gabor_2d(
    const mxArray*         array_ptr,
    wav_params_gabor_2d*   struct_ptr );




#ifdef __cplusplus
}
#endif

#endif
/* End of ut_mtlb6.h */
