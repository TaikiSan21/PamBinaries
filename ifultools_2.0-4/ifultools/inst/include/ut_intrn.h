
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_intrn.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_INTRN_H
#define IN_UT_INTRN_H

#include "ut_err.h"
#include "ut_type.h"

#include <stdlib.h>
#include <stdio.h>

/*
   Utilities (ut_) for internal (intrn) Insightful use (ut_intrn).
   This file contains function declarations for private functions that
   are internal to the mutils library, but that are used by more than
   one c file, and therefore cannot be static functions.
*/


#ifdef __cplusplus
extern "C" {
#endif


/** Copy data to/from fixed/actual width arrays.
 * Convert data of a certain width and byte order to data of another
 * width and byte order, padding blank bytes with zeros if necessary.
 * Primarily used as preprocessing for writing raw data files,
 * and post-processing for reading raw data files.
 *
 * This internal function does not call the error message function if it
 * returns an error code, and it should never return an error code
 * unless it is called with improper inputs, NULL pointers.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrn.h
 * @source ut\_intrn.c
 * @library matrix
 * @usage  To take uint32 (4 or more bytes per number, depending on platform),
 *   big-endian data and store in exactly 4-byte little-endian
 *   array for output:
 *   #err_code = mutil_copy_bytes((uint8 *) mylongarray, nelem, 4, sizeof(uint32), 4, TRUE, FALSE, mybuf );#
 *   To take exactly 4-byte little-endian input read from a file and put it
 *     into a big-endian uint32:
 *   #err_code = mutil_copy_bytes(mybuf, nelem, 4, 4, sizeof(uint32), FALSE, TRUE, (uint8 *) mylongarray );#
 * @return    Standard mutils error/OK code.
 * @param in_buf         Array of bytes to switch.
 * @param num_nums       Number of numbers (not bytes) in array.
 * @param eff_size       Effective number of bytes per number in data.
 * @param in_size        Actual number of bytes per number in input data.
 * @param out_size       Actual number of bytes per number in output data.
 * @param in_big_endian  Flag telling whether input is big-endian or not.
 * @param out_big_endian Flag telling whether output is big-endian or not.
 * @param out_buf        Array of bytes for output, pre-allocated to
 *           out\_size * num\_nums.
 * @see MUTIL_SYSTEM_IS_BIGENDIAN
 * @see MUTIL_SYSTEM_IS_LITTLEENDIAN
 * @private
 */
extern mutil_errcode mutil_copy_bytes( uint8 *in_buf, sint32 num_nums,
  size_t eff_size, size_t in_size, size_t out_size,
  boolean in_big_endian, boolean out_big_endian, uint8 *out_buf );


/** Pack 1-bit data into a buffer.
 * Take an array of data, and convert it into 1-bit data packed
 * into an 8-bit output buffer.  The input values can be converted to
 * bits either so that non-zero values become 0 or so that non-zero
 * values become 1.
 *
 * This internal function does not call the error message function if it
 * returns an error code, and it should never return an error code
 * unless it is called with improper inputs, such as NULL pointers.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrn.h
 * @source ut\_intrn.c
 * @library matrix
 * @usage  #err = mutil_pack_bits(mydata, ndat, TRUE, buf);#
 * @return    Standard mutils error/OK code.
 * @param data     Data array to pack.
 * @param ndata    Length of data array.
 * @param onetrue  If TRUE, convert non-zero data to 1 bits; otherwise
 *     convert non-zero data to 0 bits.
 * @param outbuf   Output buffer, of length ndata/8 (rounded up).
 * @see mutil_unpack_bits
 * @private
 */
extern mutil_errcode mutil_pack_bits( uint8 *data, sint32 ndata,
  boolean onetrue, uint8 *outbuf );


/** Unpack 1-bit data from a buffer.
 * Take a buffer of 1-bit data, and convert it into an 8-bit data array.
 * The input values can be converted from
 * bits either so that 0 bits become data value 0 or so that
 * 0 bits become data value 1.
 *
 * This internal function does not call the error message function if it
 * returns an error code, and it should never return an error code
 * unless it is called with improper inputs, such as NULL pointers.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrn.h
 * @source ut\_intrn.c
 * @library matrix
 * @usage  #err = mutil_unpack_bits(buf, ndat, TRUE, outdat);#
 * @return    Standard mutils error/OK code.
 * @param inbuf    Input buffer, of length ndata/8 (rounded up).
 * @param ndata    Length of data array.
 * @param onetrue  If TRUE, convert 1 bits to 1 data; otherwise
 *     convert 1 bits to 0 data.
 * @param data     Data array to put unpacked output in.
 * @see mutil_pack_bits
 * @private
 */
extern mutil_errcode mutil_unpack_bits( uint8 *inbuf, sint32 ndata,
  boolean onetrue, uint8 *data );


/** Skip comments and separator characters in a file.
 * Read characters from an ascii file, skipping the given
 * separator character and whitespace, and also skipping comment
 * characters and the remainders of the lines containing them.
 * Return a non-error value when any other characters are encountered;
 * return an error value and generate an error message
 * if the end of the file is reached.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrn.h
 * @source ut\_intrn.c
 * @library matrix
 * @limits If the comment character can
 *    be read as a number or whitespace, it will be ignored.
 * @usage  #err_code = localfn_skip_comments(file, ',', 'C');#
 * @return    Standard mutils error/OK code.
 * @param in_file     Pointer to file for reading.
 * @param sep_char    Character (besides whitespace) that separates data.
 * @param comment_char Character that designates rest of line is comment.
 * @private
 */
extern mutil_errcode mutil_skip_comments( FILE *in_file, char sep_char,
  char comment_char );


/** Macro to free matrix with warning.
 * This macro wraps around a call to any matrix free function and if the
 * free function fails, it warns of a possible memory leak
 * by calling \Ref{MUTIL_WARN}.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrn.h
 * @source ut\_intrn.h
 * @library matrix
 * @usage  To free a universal matrix by calling the
 *   matuniv\_free function: #MUTIL_FREE_WARN(matuniv, mybuf);#
 * @param FREE_FN_PREFIX   Prefix to the free function to use,
 *         e.g., matuniv, matdbl, etc.
 * @param MAT_PTR          Pointer to the matrix to free.
 */

#define MUTIL_FREE_WARN( FREE_FN_PREFIX, MAT_PTR ) \
   if( FREE_FN_PREFIX ## _free( MAT_PTR )) \
     MUTIL_WARN( "Free failed, possible memory leak" )


/** Macro to free a buffer of memory with warning.
 * This macro calls \Ref{mutil_free} to free a buffer of
 * previously allocated memory. If the operation fails, it
 * warns of a possible memory leak by calling \Ref{MUTIL_WARN}
 * The macro casts the data argument to void* and
 * the size argument to sint32, as required by mutil\_free.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrn.h
 * @source ut\_intrn.h
 * @library matrix
 * @usage  #MUTIL_FREE_BUFFER_WARN(DATA, SIZE);#
 * @param DATA A pointer to the buffer to be freed.
 * @param SIZE Size of the buffer to be freed.
 * @see MUTIL_FREE_WARN
 */

#define MUTIL_FREE_BUFFER_WARN( DATA, SIZE ) \
   if( mutil_free( (void *) (DATA), (sint32) (SIZE) ) ) \
      MUTIL_WARN("Free failed: possible memory leak")

/** Macro to completely free a matrix set with warning.
 * This macro calls \Ref{matset_matrices_free} and
 * \Ref{matset_free} to free data and headers of a matrix set, and
 * if either fails, warns of a possible memory leak.  It is two statements, so
 * be careful to enclose in brackets.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrn.h
 * @source ut\_intrn.h
 * @library matrix
 * @usage  #MUTIL_FREEALL_MATSET_WARN(myset);#
 * @param SET A pointer to the matrix set to be freed.
 * @see MUTIL_FREE_WARN
 * @see matset_matrices_free
 * @see matset_free
 */
#define MUTIL_FREEALL_MATSET_WARN( SET ) \
   MUTIL_FREE_WARN( matset_matrices, SET ); \
   MUTIL_FREE_WARN( matset, SET )


/** Macro to implement standard casting.
 * This macro is for use in functions that use the standard
 * mutils casting conventions.  It will cast a matrix to the
 * output matrix type, if it is legal, or wrap the matrix if it
 * is already the correct type. It keeps track of whether the matrix was
 * allocated and cast or wrapped in a boolean variable.  It
 * checks to see that the matrix can be cast without data loss,
 * sets an error variable if any of the operations fail,
 * and assures that MUTIL\_ERROR has been called if there is an error.
 *
 * Note that this expands to a set of statements.  Use
 * \Ref{MUTIL_FREE_STANDARD_CASTING} to free the resulting matrix.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrn.h
 * @source ut\_intrn.h
 * @library matrix
 * @usage  #MUTIL_DO_STANDARD_CASTING(in_image, out_image->type, cast_in_image, intrp_ptr, in_image_alloc, errcode);#
 * @param  IN_MATRIX    Universal matrix to cast or wrap.
 * @param  MAT_TYPE     Data type to cast to.
 * @param  CAST_MATRIX  Universal matrix variable to use for cast/wrapped
 *    result.
 * @param  INTRP_PTR    Pointer for interrupt checking.
 * @param  DID_ALLOC    boolean variable to store whether matrix was
 *    allocated and cast (TRUE) or wrapped (FALSE).
 * @param  ERR_VAR      mutil\_errcode variable set with error code of
 *    cast, alloc, or wrap operation.
 * @see MUTIL_FREE_STANDARD_CASTING
 * @see MUTIL_CHECK_STANDARD_CASTING_SCALAR
 */
#define MUTIL_DO_STANDARD_CASTING( IN_MATRIX, MAT_TYPE,  \
     CAST_MATRIX, INTRP_PTR, DID_ALLOC, ERR_VAR ) \
  DID_ALLOC = FALSE; \
  if( (IN_MATRIX)->type == (MAT_TYPE) ) { \
    ERR_VAR = matuniv_wrap_data( &CAST_MATRIX, \
      (void *) ( MATUNIV_DATA( IN_MATRIX )), \
      MATUNIV_NROW( IN_MATRIX ), MATUNIV_NCOL( IN_MATRIX ), MAT_TYPE ); \
  } \
  else { \
    /* do not allow cast of floating to integer */ \
    if(( (IN_MATRIX)->type == MUTIL_DOUBLE || \
         (IN_MATRIX)->type == MUTIL_FLOAT || \
         (IN_MATRIX)->type == MUTIL_DCOMPLEX ) && \
       ( (MAT_TYPE) != MUTIL_DOUBLE && \
         (MAT_TYPE) != MUTIL_FLOAT && \
         (MAT_TYPE) != MUTIL_DCOMPLEX )) { \
      MUTIL_ERROR( "Incompatible types for input and output, cannot cast" ); \
      ERR_VAR = MUTIL_ERR_ILLEGAL_TYPE; \
    } \
    else { \
      /* allocate space; cast with clip = FALSE to assure no loss of data */ \
      ERR_VAR = matuniv_malloc( &CAST_MATRIX,  \
        MATUNIV_NROW( IN_MATRIX ), MATUNIV_NCOL( IN_MATRIX ), MAT_TYPE ); \
      if( ERR_VAR == MUTIL_ERR_OK ) { \
        ERR_VAR = matuniv_cast( IN_MATRIX, FALSE, INTRP_PTR, &CAST_MATRIX ); \
        if( ERR_VAR ) { \
           MUTIL_FREE_WARN( matuniv, &CAST_MATRIX ); \
        } \
        else { \
          DID_ALLOC = TRUE; \
        } \
      } /* end of if(malloc succeeded) */  \
    } /* end of if/else on types compatible */ \
  } /* end of if/else on types same */



/** Macro to check range for standard casting of scalar.
 * This macro is for use in functions that use the standard
 * mutils casting conventions.  It checks that a universal
 * scalar can be cast to the output data type without overflow or
 * loss of floating-point precision.  If not, it sets an error variable,
 * and assures that MUTIL\_ERROR has been called.  Note that
 * complex numbers can never be cast.
 *
 * Note that this expands to a set of statements.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrn.h
 * @source ut\_intrn.h
 * @library matrix
 * @usage  #MUTIL_CHECK_STANDARD_CASTING_SCALAR(scalar, out_image->type, errcode, tmpdbl);#
 * @param  USCALAR      Universal scalar to check.
 * @param  DATA_TYPE     Data type to check compatibility with.
 * @param  ERR_VAR      mutil\_errcode variable set with error code of
 *    result.
 * @param  TMP_VAR      double variable used internally.
 * @see MUTIL_DO_STANDARD_CASTING
 */
#define MUTIL_CHECK_STANDARD_CASTING_SCALAR( USCALAR, DATA_TYPE, ERR_VAR, \
  TMP_VAR )  \
  if( (USCALAR).type == MUTIL_DCOMPLEX || (DATA_TYPE) == MUTIL_DCOMPLEX ) { \
    MUTIL_ERROR( "Complex types not supported" ); \
    ERR_VAR = MUTIL_ERR_ILLEGAL_TYPE; \
  } \
  else if( (USCALAR).type == (DATA_TYPE) ) { \
    ERR_VAR = MUTIL_ERR_OK; \
  } \
  else if(( (USCALAR).type == MUTIL_DOUBLE || \
            (USCALAR).type == MUTIL_FLOAT ) && \
          ( (DATA_TYPE) != MUTIL_DOUBLE && \
            (DATA_TYPE) != MUTIL_FLOAT )) { \
    /* do not allow cast of floating to integer */ \
    MUTIL_ERROR( "Incompatible types for input and output, cannot cast" ); \
    ERR_VAR = MUTIL_ERR_ILLEGAL_TYPE; \
  } \
  else { \
    /* different types, check range */ \
    TMP_VAR = SCAUNIV_CAST( USCALAR, double ); \
    ERR_VAR = MUTIL_ERR_OK; \
    switch( DATA_TYPE ) { \
      case MUTIL_SINT32: \
        if( TMP_VAR > MUTIL_SINT32_MAX || TMP_VAR < MUTIL_SINT32_MIN ) { \
          MUTIL_ERROR( "Incompatible types for input and output, cannot cast" ); \
          ERR_VAR = MUTIL_ERR_ILLEGAL_TYPE; \
        }  \
        break; \
      case MUTIL_SINT16: \
        if( TMP_VAR > MUTIL_SINT16_MAX || TMP_VAR < MUTIL_SINT16_MIN ) { \
          MUTIL_ERROR( "Incompatible types for input and output, cannot cast" ); \
          ERR_VAR = MUTIL_ERR_ILLEGAL_TYPE; \
        }  \
        break; \
      case MUTIL_SINT8: \
        if( TMP_VAR > MUTIL_SINT8_MAX || TMP_VAR < MUTIL_SINT8_MIN ) { \
          MUTIL_ERROR( "Incompatible types for input and output, cannot cast" ); \
          ERR_VAR = MUTIL_ERR_ILLEGAL_TYPE; \
        }  \
        break; \
      case MUTIL_UINT32: \
        if( TMP_VAR > MUTIL_UINT32_MAX || TMP_VAR < 0 ) { \
          MUTIL_ERROR( "Incompatible types for input and output, cannot cast" ); \
          ERR_VAR = MUTIL_ERR_ILLEGAL_TYPE; \
        } \
        break; \
      case MUTIL_UINT16: \
        if( TMP_VAR > MUTIL_UINT16_MAX || TMP_VAR < 0 ) { \
          MUTIL_ERROR( "Incompatible types for input and output, cannot cast" ); \
          ERR_VAR = MUTIL_ERR_ILLEGAL_TYPE; \
        } \
        break; \
      case MUTIL_UINT8: \
        if( TMP_VAR > MUTIL_UINT8_MAX || TMP_VAR < 0 ) { \
          MUTIL_ERROR( "Incompatible types for input and output, cannot cast" ); \
          ERR_VAR = MUTIL_ERR_ILLEGAL_TYPE; \
        } \
        break; \
    } /* end of switch on output type */ \
  } /* end of if/else to determine compatibility */


/** Macro to free result of standard casting.
 * This macro is for use in functions that use the standard
 * mutils casting conventions which have called
 * \Ref{MUTIL_DO_STANDARD_CASTING}.  It frees the matrix,
 * warning of memory leaks on error, if it was allocated, and
 * does nothing if it was wrapped.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrn.h
 * @source ut\_intrn.h
 * @library matrix
 * @usage  #MUTIL_FREE_STANDARD_CASTING(cast_in_image, in_image_alloc);#
 * @param  CAST_MATRIX  Universal matrix variable used for cast/wrapped
 *    matrix.
 * @param  DID_ALLOC    boolean variable storing whether matrix was
 *    allocated and cast (TRUE) or wrapped (FALSE).
 * @see MUTIL_DO_STANDARD_CASTING
 */
#define MUTIL_FREE_STANDARD_CASTING( CAST_MATRIX, DID_ALLOC ) \
  if( DID_ALLOC ) { \
     MUTIL_FREE_WARN( matuniv, &(CAST_MATRIX )); \
  }


#ifdef __cplusplus
}
#endif

#endif /* ifndef IN_UT_INTRN_H*/
