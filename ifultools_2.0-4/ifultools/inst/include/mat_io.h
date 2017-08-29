
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/mat_io.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_MAT_IO_H_
#define IN_MAT_IO_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"
#include <stdio.h>

/* This file contains declarations for functions used to
   read and write matrices as ascii and binary files
   for the mutils library.
*/


#ifdef __cplusplus
extern "C" {
#endif


/******************
 Universal Matrix functions
 *****************/


/** Pretty-print a matrix.
 * Make a pretty printout of a matrix.  The first line of
 * the printout indicates the matrix type (so that it can be parsed
 * later), and subsequent lines contain the matrix data.
 *
 * This function tries to be pretty about printing the results.
 * For example, it tries to skip to new lines to avoid ugly wraparound.
 * Since each element may occupy a different
 * number of characters, a fixed format should be supplied to get the
 * columns neatly lined up. Also, if you want to be able to read
 * the matrix later, the format string should only print a single number,
 * followed by at least one space and preceded by optional whitespace.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_io.h
 * @source mat\_io.c
 * @library matrix
 * @usage  #err_code = matuniv_print(file, &my_mat, "%10.2g", 80, intrp_ptr);#
 * @return    Standard mutils error/OK code.
 * @param out_file   Pointer to file for writing.
 * @param mat        Pointer to the matrix to be printed.
 * @param fmt        Printf-style  format string;
 *      if NULL, uses default format for each matrix type (listed below).
 * @param line_width Maximum approximate number of characters to print per
 *      line; if negative or zero, defaults to 75.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @same \begin{itemize}
 *   \item #MUTIL_LIBEXPORT mutil_errcode matdbl_print(FILE *out_file, const double_mat *mat, const char *fmt, int line_width, void *intrp_ptr);#
 *  (default format is #"%g "#)
 *   \item #MUTIL_LIBEXPORT mutil_errcode matflt_print(FILE *out_file, const float_mat *mat, const char *fmt, int line_width, void *intrp_ptr);#
 *  (default format is #"%g "#)
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu8_print(FILE *out_file, const uint8_mat *mat, const char *fmt, int line_width, void *intrp_ptr);#
 *  (default format is #"%hu "#)
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu16_print(FILE *out_file, const uint16_mat *mat, const char *fmt, int line_width, void *intrp_ptr);#
 *  (default format is #"%hu "#)
 *   \item #MUTIL_LIBEXPORT mutil_errcode matu32_print(FILE *out_file, const uint32_mat *mat, const char *fmt, int line_width, void *intrp_ptr);#
 *  (default format is #"%lu "#)
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats16_print(FILE *out_file, const sint16_mat *mat, const char *fmt, int line_width, void *intrp_ptr);#
 *  (default format is #"%hd "#)
 *   \item #MUTIL_LIBEXPORT mutil_errcode mats32_print(FILE *out_file, const sint32_mat *mat, const char *fmt, int line_width, void *intrp_ptr);#
 *  (default format is #"%ld "#)
 * \end{itemize}
 * @see matuniv_parse
 * @see matuniv_write_ascii
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_print( FILE *out_file,
  const univ_mat *mat, const char *fmt, int line_width,  void *intrp_ptr );


/** Parse a pretty-printed matrix.
 * Parse a matrix that was printed via \Ref{matuniv_print},
 * and store the result in a matrix that should not have been
 * previously allocated.
 * The calling program is responsible for freeing the matrix.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_io.h
 * @source mat\_io.c
 * @library matrix
 * @usage  #err_code = matuniv_parse(file, intrp_ptr, &new_mat);#
 * @return    Standard mutils error/OK code.
 * @param in_file    Pointer to the input file.
 * @param intrp_ptr  Pointer for implementation of interrupt checking.
 * @param mat        Pointer to an unallocated matrix.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_parse(FILE *in_file, void *intrp_ptr, double_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_parse(FILE *in_file, void *intrp_ptr, float_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_parse(FILE *in_file, void *intrp_ptr, uint8_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_parse(FILE *in_file, void *intrp_ptr, uint16_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_parse(FILE *in_file, void *intrp_ptr, uint32_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_parse(FILE *in_file, void *intrp_ptr, sint16_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_parse(FILE *in_file, void *intrp_ptr, sint32_mat *mat);#
 * \end{itemize}
 * @see matuniv_print
 * @see matuniv_read_ascii
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_parse( FILE *in_file, void *intrp_ptr,
  univ_mat *mat );


/** Write a matrix to an ASCII text file.
 * Write the data in a matrix, as ASCII text,
 * into a file  that was previously opened for writing (possibly already
 * containing a header).  The values are written as decimal numbers,
 * separated by the given string (or a single space if it is NULL),
 * with up to the given number of characters printed per line.  No
 * attempt is made to line up the numbers or to separate them into
 * rows as in the matrix.
 * The values in the matrix are written out in the order that they
 * appear in the matrix's data array (row-order).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_io.h
 * @source mat\_io.c
 * @library matrix
 * @limits Complex matrices cannot be written with this function.
 * @usage  #err_code = matuniv_write_ascii(file, &my_mat, ", ", 70, intrp_ptr);#
 * @return    Standard mutils error/OK code.
 * @param out_file      Pointer to file for writing.
 * @param mat           Pointer to the matrix to be written.
 * @param sep           Character string to use to separate the values; or
 *     NULL to use a single space.
 * @param char_per_line Maximum number of characters to print out per line.
 * @param intrp_ptr     Pointer for implementation of interrupt checking.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_write_ascii(FILE *out_file, const double_mat *mat, const char *sep, uint16 char_per_line, void *intrp_ptr);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_write_ascii(FILE *out_file, const float_mat *mat, const char *sep, uint16 char_per_line, void *intrp_ptr);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_write_ascii(FILE *out_file, const uint8_mat *mat, const char *sep, uint16 char_per_line, void *intrp_ptr);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_write_ascii(FILE *out_file, const uint16_mat *mat, const char *sep, uint16 char_per_line, void *intrp_ptr);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_write_ascii(FILE *out_file, const uint32_mat *mat, const char *sep, uint16 char_per_line, void *intrp_ptr);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_write_ascii(FILE *out_file, const sint16_mat *mat, const char *sep, uint16 char_per_line, void *intrp_ptr);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_write_ascii(FILE *out_file, const sint32_mat *mat, const char *sep, uint16 char_per_line, void *intrp_ptr);#
 * \end{itemize}
 * @see matuniv_read_ascii
 * @see matuniv_write_raw
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_write_ascii( FILE *out_file,
  const univ_mat *mat, const char *sep, uint16 char_per_line,
  void *intrp_ptr );


/** Read a matrix from an ASCII text file.
 * Read data from an ASCII file previously opened for reading,
 * and possibly already repositioned, into a previously-allocated
 * matrix.  The values are read as decimal numbers,
 * separated by whitespace and/or the given separator character,
 * and portions of lines following the given comment character
 * are ignored.  Line boundaries in the file are not used for
 * any special purpose beyond separating numbers, like any other
 * whitespace, and ending comments.  (In particular, they are not
 * used to signal ends of matrix rows or columns.)
 *
 * If the data in the file are not consistent with the number of values and
 * data type of the matrix, an error will be generated.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_io.h
 * @source mat\_io.c
 * @library matrix
 * @limits Complex matrices cannot be read with this function.
 *    If the comment character can
 *    be read as a number or whitespace, it will be ignored.
 * @usage  #err_code = matuniv_read_ascii(file, ',', 'C', intrp_ptr, &my_mat);#
 * @return    Standard mutils error/OK code.
 * @param in_file     Pointer to file for reading.
 * @param sep         Character (besides whitespace) that separates values.
 * @param comment     Character that designates rest of line is comment.
 * @param intrp_ptr   Pointer for implementation of interrupt checking.
 * @param mat         Pointer to matrix to put the data in.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_read_ascii(FILE *in_file, char sep, char comment, void *intrp_ptr, double_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_read_ascii(FILE *in_file, char sep, char comment, void *intrp_ptr, float_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_read_ascii(FILE *in_file, char sep, char comment, void *intrp_ptr, uint8_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_read_ascii(FILE *in_file, char sep, char comment, void *intrp_ptr, uint16_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_read_ascii(FILE *in_file, char sep, char comment, void *intrp_ptr , uint32_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_read_ascii(FILE *in_file, char sep, char comment, void *intrp_ptr, sint16_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_read_ascii(FILE *in_file, char sep, char comment, void *intrp_ptr , sint32_mat *mat);#
 * \end{itemize}
 * @see matuniv_write_ascii
 * @see matuniv_read_raw
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_read_ascii( FILE *in_file, char sep,
  char comment, void *intrp_ptr, univ_mat *mat);


/** Write a matrix to a raw data file.
 * Write the data in a matrix, in raw (binary) packed format, into a file
 * that was previously opened for writing (possibly already
 * containing a header).  The values are written out using IEEE
 * formats for floating-point data (this function will return an
 * error if the machine is not using IEEE sizes for floats and
 * doubles), and using the specified number of bits for the signed
 * and unsigned integers (this function currently returns an error
 * on signed data if the machine is not using the exact number of
 * bytes specified by the type).
 *
 * If the matrix contains data of more than one byte
 * width, the bytes of each value can be written out either with the
 * high-order byte first (big-endian), or the low-order byte first
 * (little-endian).
 *
 * The values in the matrix are written out in the order that they
 * appear in the matrix's data array (row-order).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_io.h
 * @source mat\_io.c
 * @library matrix
 * @usage  #err_code = matuniv_write_raw(file, &my_mat, TRUE, intrp_ptr);#
 * @return    Standard mutils error/OK code.
 * @param out_file     Pointer to file for writing.
 * @param mat          Pointer to the matrix to be written.
 * @param high_first   If TRUE, write high-order byte first; if FALSE, write
 *     low-order byte first.
 * @param intrp_ptr    Pointer for implementation of interrupt checking.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_write_raw(FILE *out_file, const double_mat *mat, boolean high_first, void *intrp_ptr);#
 *  (must be 8-byte double)
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_write_raw(FILE *out_file, const float_mat *mat, boolean high_first, void *intrp_ptr);#
 *  (must be 4-byte float)
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_write_raw(FILE *out_file, const uint8_mat *mat, void *intrp_ptr);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_write_raw(FILE *out_file, const uint16_mat *mat, boolean high_first, void *intrp_ptr);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_write_raw(FILE *out_file, const uint32_mat *mat, boolean high_first, void *intrp_ptr);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_write_raw(FILE *out_file, const sint16_mat *mat, boolean high_first, void *intrp_ptr);#
 * (must be 2-byte sint16)
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_write_raw(FILE *out_file, const sint32_mat *mat, boolean high_first, void *intrp_ptr);#
 * (must be 4-byte sint32)
 * \end{itemize}
 * @see matuniv_read_raw
 * @see matuniv_write_ascii
 * @see scauniv_write_raw
 * @see scauniv_read_raw
 * @see mutil_write_raw
 * @see mutil_read_raw
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_write_raw( FILE *out_file,
  const univ_mat *mat, boolean high_first, void *intrp_ptr );


/** Read a matrix from a raw data file.
 * Read data in a raw (binary) packed data file
 * that was previously opened for reading, and possibly
 * already repositioned past a header, into a previously-allocated
 * matrix.  The type of raw data and number of values
 * in the file are taken from the matrix header of the output matrix,
 * and the values in the file are assumed to be in row-major order.
 * The values are read using IEEE
 * formats for floating-point data (this function will return an
 * error if the machine is not using IEEE sizes for floats and
 * doubles), and using the specified number of bits for the signed
 * and unsigned integers (this function currently returns an error
 * on signed data if the machine is not using the exact number of
 * bytes specified by the type).
 *
 * If the matrix contains data of more than one byte
 * width, the bytes of each value can be read either with the
 * high-order byte first (big-endian), or the low-order byte first
 * (little-endian).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include mat\_io.h
 * @source mat\_io.c
 * @library matrix
 * @usage  #err_code = matuniv_read_raw(file, 0, TRUE, intrp_ptr, &my_mat);#
 * @return    Standard mutils error/OK code.
 * @param in_file      Pointer to file to be read.
 * @param skip         Number of bytes to skip before reading data.
 * @param high_first   If TRUE, read high-order byte first; if FALSE, read
 *     low-order byte first.
 * @param intrp_ptr    Pointer for implementation of interrupt checking.
 * @param mat          Pointer to the matrix to put the data in.
 * @same \begin{itemize}
 *  \item #MUTIL_LIBEXPORT mutil_errcode matdbl_read_raw(FILE *in_file, sint32 skip, boolean high_first, void *intrp_ptr, double_mat *mat);#
 *  (must be 8-byte double)
 *  \item #MUTIL_LIBEXPORT mutil_errcode matflt_read_raw(FILE *in_file, sint32 skip, boolean high_first, void *intrp_ptr, float_mat *mat);#
 *  (must be 4-byte float)
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu8_read_raw(FILE *in_file, sint32 skip, void *intrp_ptr, uint8_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu16_read_raw(FILE *in_file, sint32 skip, boolean high_first, void *intrp_ptr, uint16_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode matu32_read_raw(FILE *in_file, sint32 skip, boolean high_first, void *intrp_ptr, uint32_mat *mat);#
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats16_read_raw(FILE *in_file, sint32 skip, boolean high_first, void *intrp_ptr, sint16_mat *mat);#
 *  (must be 2-byte sint16)
 *  \item #MUTIL_LIBEXPORT mutil_errcode mats32_read_raw(FILE *in_file, sint32 skip, boolean high_first, void *intrp_ptr, sint32_mat *mat);#
 *  (must be 4-byte sint32)
 *  \end{itemize}
 * @see matuniv_write_raw
 * @see matuniv_read_ascii
 * @see scauniv_write_raw
 * @see scauniv_read_raw
 * @see mutil_write_raw
 * @see mutil_read_raw
 * @see Matrix Data Types
 * @see Interrupt Handling
 */
MUTIL_LIBEXPORT mutil_errcode matuniv_read_raw( FILE *in_file, sint32 skip,
  boolean high_first, void *intrp_ptr, univ_mat *mat );


/******************
 Non-universal Matrix functions
 *****************/


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matdbl_print( FILE *out_file,
  const double_mat *mat, const char *fmt, int line_width, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matflt_print(FILE *out_file,
  const float_mat *mat, const char *fmt, int line_width, void *intrp_ptr);


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu8_print( FILE *out_file,
  const uint8_mat *mat, const char *fmt, int line_width, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu16_print( FILE *out_file,
  const uint16_mat *mat, const char *fmt, int line_width, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu32_print( FILE *out_file,
  const uint32_mat *mat, const char *fmt, int line_width, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats16_print( FILE *out_file,
  const sint16_mat *mat, const char *fmt, int line_width, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats32_print( FILE *out_file,
  const sint32_mat *mat, const char *fmt, int line_width, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matdbl_parse( FILE *in_file, void *intrp_ptr,
  double_mat *mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matflt_parse(FILE *in_file, void *intrp_ptr,
  float_mat *mat);


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu8_parse( FILE *in_file, void *intrp_ptr,
  uint8_mat *mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu16_parse( FILE *in_file, void *intrp_ptr,
  uint16_mat *mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu32_parse( FILE *in_file, void *intrp_ptr,
  uint32_mat *mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats16_parse( FILE *in_file, void *intrp_ptr,
  sint16_mat *mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats32_parse( FILE *in_file, void *intrp_ptr,
  sint32_mat *mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matdbl_write_ascii( FILE *out_file,
  const double_mat *mat, const char *sep, uint16 char_per_line,
  void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matflt_write_ascii( FILE *out_file,
  const float_mat *mat,
  const char *sep, uint16 char_per_line, void *intrp_ptr );

/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu8_write_ascii( FILE *out_file,
  const uint8_mat *mat, const char *sep, uint16 char_per_line,
  void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu16_write_ascii( FILE *out_file,
  const uint16_mat *mat, const char *sep, uint16 char_per_line,
  void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu32_write_ascii( FILE *out_file,
  const uint32_mat *mat, const char *sep, uint16 char_per_line,
  void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats16_write_ascii( FILE *out_file,
  const sint16_mat *mat, const char *sep, uint16 char_per_line,
  void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats32_write_ascii( FILE *out_file,
  const sint32_mat *mat, const char *sep, uint16 char_per_line,
  void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matdbl_read_ascii( FILE *in_file, char sep,
  char comment, void *intrp_ptr, double_mat *mat);


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matflt_read_ascii(FILE *in_file, char sep,
  char comment, void *intrp_ptr, float_mat *mat);


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu8_read_ascii( FILE *in_file, char sep,
  char comment, void *intrp_ptr, uint8_mat *mat);


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu16_read_ascii( FILE *in_file, char sep,
  char comment, void *intrp_ptr, uint16_mat *mat);


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu32_read_ascii( FILE *in_file, char sep,
  char comment, void *intrp_ptr , uint32_mat *mat);


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats16_read_ascii( FILE *in_file, char sep,
  char comment, void *intrp_ptr, sint16_mat *mat);


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats32_read_ascii( FILE *in_file, char sep,
  char comment, void *intrp_ptr , sint32_mat *mat);


/* Documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matdbl_write_raw( FILE *out_file,
  const double_mat *mat, boolean high_first, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matflt_write_raw( FILE *out_file,
  const float_mat *mat, boolean high_first, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu8_write_raw( FILE *out_file,
  const uint8_mat *mat, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu16_write_raw( FILE *out_file,
  const uint16_mat *mat, boolean high_first, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu32_write_raw( FILE *out_file,
  const uint32_mat *mat, boolean high_first, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats16_write_raw( FILE *out_file,
  const sint16_mat *mat, boolean high_first, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats32_write_raw( FILE *out_file,
  const sint32_mat *mat, boolean high_first, void *intrp_ptr );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matdbl_read_raw( FILE *in_file, sint32 skip,
  boolean high_first, void *intrp_ptr, double_mat *mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matflt_read_raw( FILE *in_file, sint32 skip,
  boolean high_first, void *intrp_ptr, float_mat *mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu8_read_raw( FILE *in_file, sint32 skip,
  void *intrp_ptr, uint8_mat *mat);


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu16_read_raw( FILE *in_file, sint32 skip,
  boolean high_first, void *intrp_ptr, uint16_mat *mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode matu32_read_raw( FILE *in_file, sint32 skip,
  boolean high_first, void *intrp_ptr, uint32_mat *mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats16_read_raw( FILE *in_file, sint32 skip,
  boolean high_first, void *intrp_ptr, sint16_mat *mat );


/* Function documented with universal matrix version, above */
MUTIL_LIBEXPORT mutil_errcode mats32_read_raw( FILE *in_file, sint32 skip,
  boolean high_first, void *intrp_ptr, sint32_mat *mat );


#ifdef __cplusplus
}
#endif

#endif /* IN_MAT_IO_H_*/
