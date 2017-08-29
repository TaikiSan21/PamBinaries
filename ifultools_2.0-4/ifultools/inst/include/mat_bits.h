
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_bits.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_BITS_H_
#define IN_MAT_BITS_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include <stdio.h>

/* This file contains declarations for functions involving specialized
   bit manipulation, such as compaction, for the mutils library.  */


#ifdef __cplusplus
extern "C" {
#endif


/** Bit compaction of fixed length code.
 * Bit compaction of symbols from a discrete-source,
 * fixed-length alphabet of N consecutive symbols indexed from 0 to N-1.
 * This function can be used to do both compaction and decompaction.
 * This function is useful for breaking symbols from a long alphabet
 * into symbols from a short alphabet for either source or channel
 * coding.  It is also useful for combining short alphabet symbols
 * which may result from source or channel decoding to recreate
 * symbols from a longer alphabet.
 *
 * Example 1:
 * For big\_endian = TRUE, if nbits\_in=3 and nbits\_out=6,
 * and the input index is [1,4,2,5,3,0], the result, if length 3,
 * will be [12,21,24].
 * For big\_endian = TRUE, if nbits\_in=6 and nbits\_out=3,
 * and the input index is [12,21,24], the result, if length 6,
 * will be [1,4,2,5,3,0].
 *
 * Example 2:
 * For big\_endian = FALSE, if nbits\_in=3 and nbits\_out=6,
 * and the input index is [1,4,2,5,3,0], the result, if length 3,
 * will be [33,42,3].
 * For big\_endian = FALSE, if nbits\_in=6 and nbits\_out=3,
 * and the input index is [33,42,3], the result, if length 6,
 * will be [1,4,2,5,3,0].
 *
 * Example 3:
 * For big\_endian = TRUE, if nbits\_in=3 and nbits\_out=8,
 * and the input index is [1,4,2,5,3,0], the result, if length 3,
 * will be [49,86,8].
 * For big\_endian = TRUE, if nbits\_in=8 and nbits\_out=3,
 * and the input index is [49,86,8], the result: if length 4,
 * will be [1,4,2,5];  and if length 8, will be [1,4,2,5,3,0,0,0].
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_bits.h
 * @source mat\_bits.c
 * @library matrix
 *
 * @usage  #err_code = matuniv_bit_compact_fixed_length_code(&index, 3, 8, TRUE, intrp_ptr, &result);#
 *
 * @return         Standard mutils error/OK code.
 *
 * @param index      Indices from a consecutively indexed discrete alphabet
 *   source which range in value from 0 to N-1.  The matrix type must
 *   be a signed or unsigned integer type.
 * @param nbits_in   The number of bits used in each index array element.
 *   If the index type is
 *   MUTIL\_UINT8 then nbits\_in should be 1 through 8.  If the index
 *   type is MUTIL\_UINT16 then nbits\_in should be 1 through 16.
 *   If the index type is MUTIL\_SINT16 then nbits\_in should be
 *   1 through 15.  If the index type is MUTIL\_UINT32 then
 *   nbits\_in should be 1 through 32.  If the index
 *   type is MUTIL\_SINT32 then nbits\_in should be 1 through 31.
 * @param nbits_out  The number of bits put into each output array element.
 *   If the result type is
 *   MUTIL\_UINT8 then nbits\_out should be 1 through 8.  If the result
 *   type is MUTIL\_UINT16 then nbits\_out should be 1 through 16.
 *   If the result
 *   type is MUTIL\_UINT32 then nbits\_out should be 1 through 32.
 * @param big_endian Flag to indicate whether groups of nbits\_in
 *    bits should be placed in big or little endian order.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param result     The bit-compacted string. The matrix type must be
 *   an unsigned integer type. The length of the index controls the number
 *   of bits assigned.  The result pointer cannot point to the same memory
 *   as the index pointer.
 *
 * @see matuniv_bit_and
 * @see matuniv_bit_or
 * @see matuniv_bit_not
 * @see matuniv_bit_xor
 * @see matuniv_bit_and_scalar
 * @see matuniv_bit_or_scalar
 * @see matuniv_bit_xor_scalar
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_compact_fixed_length_code(
  const univ_mat *index,
  uint32          nbits_in,
  uint32          nbits_out,
  boolean         big_endian,
  void           *intrp_ptr,
  univ_mat       *result );


/** Bit compaction of prefix code.
 * Bit compaction of symbols from a discrete-source,
 * prefix code of N consecutive symbols indexed from 0 to N-1.
 * This function can be used to do only compaction.
 * Use \Ref{matuniv_bit_decompact_prefix_code} to do decompaction.
 * This function is useful for encoding symbols from a fixed-length
 * alphabet into symbols from a prefix code for source
 * coding, such as Huffman coding.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_bits.h
 * @source mat\_bits.c
 * @library matrix
 *
 * @usage  #err_code = matuniv_bit_compact_prefix_code(&indices, &table, intrp_ptr, &result);#
 *
 * @return         Standard mutils error/OK code.
 *
 * @param index         Pointer to one-column integral-type universal
 *    matrix containing input alphabet symbol indices to be encoded.
 *    The alphabet symbols are assumed to be indices from 0 to N-1.
 * @param prefix_table  Pointer to one-column universal matrix
 *    containing the table of prefix symbols.  The symbols must
 *    be greater than 1.  The prefix table symbols must be positive
 *    integers  assigned so that the binary word can be extracted
 *    as the last N bits of the index where
 *    $N=\lfloor{\log_{2}\mbox{prefix\_table}}\rfloor$.
 *    This functions does {\em not} verify if the table is a legitimate
 *    prefix code.  The data type must be MUTIL\_SINT32.
 * @param intrp_ptr     Pointer for implementation of interrupt checking.
 * @param bit_stream    Pointer to unallocated universal matrix
 *    in which to return the data.  Data will be returned as a one-row
 *    matrix of type MUTIL\_UINT8.  Any unused bits in the last byte will
 *    be padded with zeros. Calling program is responsible for freeing.
 *    NULL to omit this output.
 *
 * @see matuniv_bit_compact_fixed_length_code
 * @see matuniv_bit_decompact_prefix_code
 * @see cdcuniv_table_create_huffman
 * @see cdcuniv_table_encode_huffman
 * @see cdcuniv_table_decode_huffman
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_compact_prefix_code(
  const univ_mat *index,
  const univ_mat *prefix_table,
  void           *intrp_ptr,
  univ_mat       *bit_stream );


/** Bit decompaction of bit stream from a prefix code.
 * Bit decompaction of symbols from a discrete-source,
 * prefix code for N consecutive symbols indexed from 0 to N-1.
 * This function can be used to do only decompaction.
 * Use \Ref{matuniv_bit_compact_prefix_code} to do compaction.
 * This function is useful for decoding symbols from a prefix code
 * into symbols from a  fixed-length alphabet for source
 * decoding, such as Huffman decoding.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_bits.h
 * @source mat\_bits.c
 * @library matrix
 *
 * @usage  #err_code = matuniv_bit_decompact_prefix_code(&bit_stream, &table, intrp_ptr, &result);#
 *
 * @return         Standard mutils error/OK code.
 *
 * @param bit_stream    Pointer to one-row universal matrix
 *    of type MUTIL\_UINT8 containing encoded bit stream.
 * @param prefix_table  Pointer to one-column universal matrix
 *    containing the table of prefix symbols.  The symbols must
 *    be greater than 1.  The prefix table symbols must be positive
 *    integers  assigned so that the binary word can be extracted
 *    as the last N bits of the index where
 *    $N=\lfloor{\log_{2}\mbox{prefix\_table}}\rfloor$.
 *    This functions does {\em not} verify if the table is a legitimate
 *    prefix code.  The data type must be MUTIL\_SINT32.
 * @param intrp_ptr     Pointer for implementation of interrupt checking.
 * @param index         Pointer to one-column universal matrix that has
 *    been previously allocated to hold the desired number of output symbols.
 *    If the length is longer than the number of symbols that can be
 *    decoded from the bit stream, then the matrix will be reallocated to
 *    the correct size.  The returned alphabet indices
 *    range from 0 to N-1.  The type must be MUTIL\_SINT32.
 *
 * @see matuniv_bit_compact_fixed_length_code
 * @see matuniv_bit_compact_prefix_code
 * @see cdcuniv_table_create_huffman
 * @see cdcuniv_table_encode_huffman
 * @see cdcuniv_table_decode_huffman
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_bit_decompact_prefix_code(
  const univ_mat *bit_stream,
  const univ_mat *prefix_table,
  void           *intrp_ptr,
  univ_mat       *index );


#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_BITS_H_*/
