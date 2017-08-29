
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_io.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_IO_H_
#define IN_UT_IO_H_

#include "ut_err.h"
#include "ut_type.h"
#include "ut_plat.h"
#include <stdio.h>

/* This file contains declarations for functions used to
   read and write non-matrix data as binary files
   for the mutils library.
*/

#ifdef __cplusplus
extern "C" {
#endif


/** Write to a raw data file.
 * Write non-structured data (i. e. not a mutils matrix),
 * in raw (binary) packed format, into a file
 * that was previously opened for writing and possibly already
 * repositioned.
 * This function exists to provide a big/little endian-independent
 * fwrite for writing data from structures that are {\em not} mutils
 * There is no error or sanity checking for this
 * function.  Use \Ref{matuniv_write_raw} for matrix data.
 *
 * The values are written out using IEEE
 * formats for floating-point data (this function will return an
 * error if the machine is not using IEEE sizes for floats and
 * doubles), and using the specified number of bits for the signed
 * and unsigned integers (this function currently returns an error
 * on signed data if the machine is not using the exact number of
 * bytes specified by the type).
 *
 * If the data type has more than one byte
 * width, the bytes of each value can be written out either with the
 * high-order byte first (big-endian), or the low-order byte first
 * (little-endian).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_io.h
 * @source ut\_io.c
 * @library matrix
 * @usage  #err_code = mutil_write_raw(file, &dbl_ptr, MUTIL_DOUBLE, 1, TRUE, intrp_ptr);#
 * @return    Standard mutils error/OK code.
 * @param out_file     Pointer to file for writing.
 * @param data_ptr     Pointer to the array to be written.
 * @param type         Mutils data type.
 * @param number       Number of items of the given type to be written.
 * @param high_first   If TRUE, write high-order byte first; if FALSE, write
 *     low-order byte first.
 * @param intrp_ptr    Pointer for implementation of interrupt checking.
 * @see mutil_read_raw
 * @see matuniv_write_raw
 * @see matuniv_read_raw
 * @see scauniv_write_raw
 * @see scauniv_read_raw
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode mutil_write_raw( FILE *out_file, void *data_ptr,
  mutil_data_type type, sint32 number, boolean high_first, void *intrp_ptr );


/** Read from a raw data file.
 * Read data in a raw (binary) packed data file
 * that was previously opened for reading, and possibly
 * already repositioned, into a previously-allocated
 * array.  This function exists to provide a big/little endian-independent
 * fread for reading data into structures that are {\em not} mutils
 * matrices.  There is no error or sanity checking for this
 * function.  Use \Ref{matuniv_read_raw} for matrix data.
 *
 * The values are read using IEEE
 * formats for floating-point data (this function will return an
 * error if the machine is not using IEEE sizes for floats and
 * doubles), and using the specified number of bits for the signed
 * and unsigned integers (this function currently returns an error
 * on signed data if the machine is not using the exact number of
 * bytes specified by the type).
 *
 * If the data type has more than one byte
 * width, the bytes of each value can be read either with the
 * high-order byte first (big-endian), or the low-order byte first
 * (little-endian).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_io.h
 * @source ut\_io.c
 * @library matrix
 * @usage  #err_code = mutils_read_raw(file, MUTIL_DOUBLE, 1, TRUE, intrp_ptr, &dbl_pointer);#
 * @return    Standard mutils error/OK code.
 * @param in_file      Pointer to file to be read.
 * @param type         Mutils data type.
 * @param number       Number of items of the given type to be read.
 * @param high_first   If TRUE, read high-order byte first; if FALSE, read
 *     low-order byte first.
 * @param intrp_ptr    Pointer for implementation of interrupt checking.
 * @param data_ptr     Pointer to the array to put the data in.
 * @see mutil_write_raw
 * @see matuniv_write_raw
 * @see matuniv_read_raw
 * @see scauniv_write_raw
 * @see scauniv_read_raw
 * @see Matrix Data Types
 */
MUTIL_LIBEXPORT mutil_errcode mutil_read_raw( FILE *in_file,
  mutil_data_type type, sint32 number, boolean high_first,
  void *intrp_ptr, void *data_ptr );


/** Write a universal scalar to a raw data file.
 * Write a universal scalar in raw (binary) format, into a file
 * that was previously opened for writing and possibly already
 * repositioned.
 * This function exists to provide a big/little endian-independent
 * fwrite for writing a universal scalar to a file.
 * There is no error or sanity checking for this
 * function.
 *
 * The values are written out using IEEE
 * formats for floating-point data (this function will return an
 * error if the machine is not using IEEE sizes for floats and
 * doubles), and using the specified number of bits for the signed
 * and unsigned integers (this function currently returns an error
 * on signed data if the machine is not using the exact number of
 * bytes specified by the type).
 *
 * If the data type has more than one byte
 * width, the bytes of each value can be written out either with the
 * high-order byte first (big-endian), or the low-order byte first
 * (little-endian).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_io.h
 * @source ut\_io.c
 * @library matrix
 * @usage  #err_code = scauniv_write_raw(file, usca, TRUE, intrp_ptr);#
 * @return    Standard mutils error/OK code.
 * @param out_file     Pointer to file for writing.
 * @param usca         Universal scalar to be written.
 * @param high_first   If TRUE, write high-order byte first; if FALSE, write
 *     low-order byte first.
 * @param intrp_ptr    Pointer for implementation of interrupt checking.
 * @see scauniv_read_raw
 * @see mutil_write_raw
 * @see mutil_read_raw
 * @see matuniv_write_raw
 * @see matuniv_read_raw
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode scauniv_write_raw( FILE *out_file,
  univ_scalar usca, boolean high_first, void *intrp_ptr );

/** Read a universal scalar from a raw data file.
 * Read data in a raw (binary) data file
 * that was previously opened for reading, and possibly
 * already repositioned, into a universal scalar.
 * The type is determined by the universal scalar type.
 * This function exists to provide a big/little endian-independent
 * fread for reading data into universal scalars.
 * There is no error or sanity checking for this
 * function.
 *
 * The values are read using IEEE
 * formats for floating-point data (this function will return an
 * error if the machine is not using IEEE sizes for floats and
 * doubles), and using the specified number of bits for the signed
 * and unsigned integers (this function currently returns an error
 * on signed data if the machine is not using the exact number of
 * bytes specified by the type).
 *
 * If the data type has more than one byte
 * width, the bytes of each value can be read either with the
 * high-order byte first (big-endian), or the low-order byte first
 * (little-endian).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_io.h
 * @source ut\_io.c
 * @library matrix
 * @usage  #err_code = scauniv_read_raw(file, TRUE, intrp_ptr, &usca);#
 * @return    Standard mutils error/OK code.
 * @param in_file      Pointer to file to be read.
 * @param high_first   If TRUE, write high-order byte first; if FALSE, write
 *     low-order byte first.
 * @param intrp_ptr    Pointer for implementation of interrupt checking.
 * @param usca         Pointer to universal scalar to put the data in.
 * @see scauniv_read_raw
 * @see mutil_write_raw
 * @see mutil_read_raw
 * @see matuniv_write_raw
 * @see matuniv_read_raw
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode scauniv_read_raw( FILE *in_file,
  boolean high_first, void *intrp_ptr, univ_scalar *usca );


#ifdef __cplusplus
}
#endif

#endif /*ifndef IN_UT_IO_H_*/
