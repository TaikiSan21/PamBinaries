
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_cast.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_CAST_H_
#define IN_MAT_CAST_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

/* This file contains function declarations for casting functions
   for the mutils library.
*/


#ifdef __cplusplus
extern "C" {
#endif


/** Cast a matrix to another type.
 * Take a matrix of one type and assign its data
 * into another type of matrix, rounding if necessary.
 * Values from the input matrix that cannot be represented in the
 * type of the output matrix can be handled by clipping (truncation)
 * to the output range or by returning an error. When casting a matrix
 * of complex numbers to another of non-complex numbers, clipping means
 * dropping the imaginary part.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_cast.h
 * @source mat\_cast.c
 * @library matrix
 * @usage  #err_code = matuniv_cast(&inmat, TRUE, intrp_ptr, &outmat);#
 * @return        Standard mutils error/OK code.
 * @param mat        Pointer to the matrix to cast.
 * @param clip       If TRUE and input matrix values are illegal for
 *    the output type, clip to fit; if FALSE, return an error (argument
 *    omitted from specific matrix type functions that can never have
 *    out-of-range values).
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     Pointer to output matrix of same size.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_flt(const double_mat *mat, boolean clip, void *intrp_ptr, float_mat *result);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_cpx(const double_mat *mat, void *intrp_ptr, dcomplex_mat *result);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_u8(const double_mat *mat, boolean clip, void *intrp_ptr, uint8_mat *result);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_u16(const double_mat *mat, boolean clip, void *intrp_ptr, uint16_mat *result);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_u32(const double_mat *mat, boolean clip, void *intrp_ptr, uint32_mat *result);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_s16( const double_mat *mat, boolean clip, void *intrp_ptr, sint16_mat *result );#
 *   \item  #MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_s32(const double_mat *mat, boolean clip, void *intrp_ptr, sint32_mat *result);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_dbl(const float_mat *mat, void *intrp_ptr, double_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_cpx(const float_mat *mat, void *intrp_ptr, dcomplex_mat *result);#
 *    \item #MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_u8(const float_mat *mat, boolean clip, void *intrp_ptr, uint8_mat *result);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_u16(const float_mat *mat, boolean clip, void *intrp_ptr, uint16_mat *result);#
 *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_u32(const float_mat *mat, boolean clip, void *intrp_ptr, uint32_mat *result);#
  *   \item  #MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_s16(const float_mat *mat, boolean clip, void *intrp_ptr, sint16_mat *result);#
  *   \item  #MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_s32(const float_mat *mat, boolean clip, void *intrp_ptr, sint32_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matcpx_cast_to_dbl(const dcomplex_mat *mat, boolean clip, void *intrp_ptr, double_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_dbl(const uint8_mat *mat, void *intrp_ptr, double_mat *result );#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_flt(const uint8_mat *mat, void *intrp_ptr, float_mat *result );#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_u32(const uint8_mat *mat, void *intrp_ptr, uint32_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_u16(const uint8_mat *mat, void *intrp_ptr, uint16_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_s16(const uint8_mat *mat, void *intrp_ptr, sint16_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_s32(const uint8_mat *mat, void *intrp_ptr, sint32_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_dbl(const uint16_mat *mat, void *intrp_ptr, double_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_flt(const uint16_mat *mat, void *intrp_ptr, float_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_u8(const uint16_mat *mat, boolean clip,  void *intrp_ptr, uint8_mat *result );#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_u32(const uint16_mat *mat, void *intrp_ptr, uint32_mat *result );#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_s16(const uint16_mat *mat, boolean clip, void *intrp_ptr, sint16_mat *result );#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_s32(const uint16_mat *mat, void *intrp_ptr, sint32_mat *result );#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_dbl(const uint32_mat *mat, void *intrp_ptr, double_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_flt(const uint32_mat *mat, void *intrp_ptr, float_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_u8(const uint32_mat *mat, boolean clip, void *intrp_ptr, uint8_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_u16(const uint32_mat *mat, boolean clip, void *intrp_ptr, uint16_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_s16(const uint32_mat *mat, boolean clip, void *intrp_ptr, sint16_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_s32(const uint32_mat *mat, boolean clip, void *intrp_ptr, sint32_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_dbl(const sint16_mat *mat, void *intrp_ptr, double_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_flt(const sint16_mat *mat, void *intrp_ptr, float_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_u8(const sint16_mat *mat, boolean clip, void *intrp_ptr, uint8_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_u16(const sint16_mat *mat, boolean clip, void *intrp_ptr, uint16_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_u32(const sint16_mat *mat, boolean clip, void *intrp_ptr, uint32_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_s32(const sint16_mat *mat, void *intrp_ptr, sint32_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_dbl(const sint32_mat *mat, void *intrp_ptr, double_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_flt(const sint32_mat *mat, void *intrp_ptr, float_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_u8(const sint32_mat *mat, boolean clip, void *intrp_ptr, uint8_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_u16(const sint32_mat *mat, boolean clip, void *intrp_ptr, uint16_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_u32(const sint32_mat *mat, boolean clip, void *intrp_ptr, uint32_mat *result);#
  *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_s16(const sint32_mat *mat, boolean clip, void *intrp_ptr, sint16_mat *result);#
 *  \end{itemize}
 * @see matset_cast
 * @see matuniv_assign
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_cast( const univ_mat *mat, boolean clip,
  void *intrp_ptr, univ_mat *result );


/****************
  Specific type cast functions
  **************/


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_flt( const double_mat *mat,
  boolean clip, void *intrp_ptr, float_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_cpx( const double_mat *mat,
  void *intrp_ptr, dcomplex_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_u8( const double_mat *mat,
  boolean clip, void *intrp_ptr, uint8_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_u16( const double_mat *mat,
  boolean clip, void *intrp_ptr, uint16_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_u32( const double_mat *mat,
  boolean clip, void *intrp_ptr, uint32_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_s16( const double_mat *mat,
  boolean clip, void *intrp_ptr, sint16_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matdbl_cast_to_s32( const double_mat *mat,
  boolean clip, void *intrp_ptr, sint32_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_dbl(const float_mat *mat,
void *intrp_ptr, double_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_cpx(const float_mat *mat,
void *intrp_ptr, dcomplex_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_u8(const float_mat *mat,
boolean clip, void *intrp_ptr, uint8_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_u16(const float_mat *mat,
boolean clip, void *intrp_ptr, uint16_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_u32(const float_mat *mat,
boolean clip, void *intrp_ptr, uint32_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_s16(const float_mat *mat,
boolean clip, void *intrp_ptr, sint16_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matflt_cast_to_s32(const float_mat *mat,
boolean clip, void *intrp_ptr, sint32_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matcpx_cast_to_dbl( const dcomplex_mat *mat,
  boolean clip, void *intrp_ptr, double_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_dbl( const uint8_mat *mat,
  void *intrp_ptr, double_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_flt( const uint8_mat *mat,
  void *intrp_ptr, float_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_u16( const uint8_mat *mat,
  void *intrp_ptr, uint16_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_u32( const uint8_mat *mat,
  void *intrp_ptr, uint32_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_s16( const uint8_mat *mat,
  void *intrp_ptr, sint16_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu8_cast_to_s32( const uint8_mat *mat,
  void *intrp_ptr, sint32_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_dbl( const uint16_mat *mat,
  void *intrp_ptr, double_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_flt( const uint16_mat *mat,
  void *intrp_ptr, float_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_u8( const uint16_mat *mat,
  boolean clip, void *intrp_ptr, uint8_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_u32( const uint16_mat *mat,
  void *intrp_ptr, uint32_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_s16( const uint16_mat *mat,
  boolean clip, void *intrp_ptr, sint16_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu16_cast_to_s32( const uint16_mat *mat,
  void *intrp_ptr, sint32_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_dbl( const uint32_mat *mat,
  void *intrp_ptr, double_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_flt( const uint32_mat *mat,
  void *intrp_ptr, float_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_u8( const uint32_mat *mat,
  boolean clip, void *intrp_ptr, uint8_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_u16( const uint32_mat *mat,
  boolean clip, void *intrp_ptr, uint16_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_s16( const uint32_mat *mat,
  boolean clip, void *intrp_ptr, sint16_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode matu32_cast_to_s32( const uint32_mat *mat,
  boolean clip, void *intrp_ptr, sint32_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_dbl(const sint16_mat *mat,
  void *intrp_ptr, double_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_flt(const sint16_mat *mat,
  void *intrp_ptr, float_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_u8(const sint16_mat *mat,
  boolean clip, void *intrp_ptr, uint8_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_u16(const sint16_mat *mat,
  boolean clip, void *intrp_ptr, uint16_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_u32(const sint16_mat *mat,
  boolean clip, void *intrp_ptr, uint32_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats16_cast_to_s32(const sint16_mat *mat,
  void *intrp_ptr, sint32_mat *result);


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_dbl( const sint32_mat *mat,
  void *intrp_ptr, double_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_flt( const sint32_mat *mat,
  void *intrp_ptr, float_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_u8( const sint32_mat *mat,
  boolean clip, void *intrp_ptr, uint8_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_u16( const sint32_mat *mat,
  boolean clip, void *intrp_ptr, uint16_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_u32( const sint32_mat *mat,
  boolean clip, void *intrp_ptr, uint32_mat *result );


/* Documented with matuniv_cast, above. */
MUTIL_LIBEXPORT mutil_errcode mats32_cast_to_s16( const sint32_mat *mat,
  boolean clip, void *intrp_ptr, sint16_mat *result );


/** Assign two double matrices to a complex matrix.
 * Assign two matrices of doubles to be the real and imaginary parts
 * of a matrix of complex numbers. Either, but not both, of the pointers to the
 * double matrices may be NULL, in which case the data in the corresponding
 * real or
 * imaginary part in the complex matrix is set equal to all zeros.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_cast.h
 * @source mat\_cast.c
 * @library matrix
 * @usage  #err_code = matdbl_to_complex(&real_part, &imag_part, intrp_ptr, &outmat);#
 * @return          Standard mutils error/OK code.
 * @param real_part Pointer to the matrix containing real part.
 * @param imag_part Pointer to matrix containing imaginary part.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param result    Pointer to output complex matrix of same size as inputs.
 * @see
 * matcpx_to_double
 */
MUTIL_LIBEXPORT mutil_errcode matdbl_to_complex( const double_mat *real_part,
    const double_mat *imag_part, void *intrp_ptr, dcomplex_mat *result );


/** Separate real and imaginary parts of a complex matrix.
 * Given a matrix of complex numbers, separate the real and imaginary
 * parts into two separate matrices, each of type double.
 * Either, but not both, of the output matrix pointers may be NULL,
 * in which case only the real or the imaginary part is extracted
 * from the complex matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_cast.h
 * @source mat\_cast.c
 * @library matrix
 * @usage  #err_code = matcpx_to_double(&in_mat, intrp_ptr, &real_part, &imag_part);#
 * @return          Standard mutils error/OK code.
 * @param in_mat    Pointer to matrix of complex numbers.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param real_part Pointer to the matrix containing real part.
 *     NULL to omit output.
 * @param imag_part Pointer to matrix containing imaginary part.
 *     NULL to omit output.
 * @see matdbl_to_complex
 */
MUTIL_LIBEXPORT mutil_errcode matcpx_to_double( const dcomplex_mat *in_mat,
    void *intrp_ptr, double_mat *real_part, double_mat *imag_part );


/** Separate magnitude and phase parts of a complex matrix.
 * Given a matrix of complex numbers, separate the
 * parts into two separate matrices, each of type double.
 * Either, but not both, of the output matrix pointers may be NULL,
 * in which case only the magnitude or the phase is extracted
 * from the complex matrix.  The origin is (0,0).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_cast.h
 * @source mat\_cast.c
 * @library matrix
 * @usage  #err_code = matcpx_to_polar(&in_mat, intrp_ptr, &mag, &phase);#
 * @return          Standard mutils error/OK code.
 * @param in_mat    Pointer to matrix of complex numbers.
 * @param intrp_ptr Pointer for implementation of interrupt checking.
 * @param magnitude Pointer to the matrix containing magnitude of the
 *    complex numbers.  NULL to omit output.
 * @param phase     Pointer to matrix containing phase [0, 2$\Pi$).
 *    of the complex numbers.  NULL to omit output.
 * @see
 * matdbl_to_complex
 */
MUTIL_LIBEXPORT mutil_errcode matcpx_to_polar( const dcomplex_mat *in_mat,
    void *intrp_ptr, double_mat *magnitude, double_mat *phase );


#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_CAST_H_*/
