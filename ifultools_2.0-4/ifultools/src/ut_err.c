
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/ut_err.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "ut_err.h"
#include "ut_debug.h"

/*
   This file contains the implementations of the function(s)
   declared in ut_err.h.
*/


/* Array of character strings for error messages */
static const char *error_messages[] = {
  /* MUTIL_ERR_OK */
  "OK, no error",

  /* MUTIL_ERR_MEM_ALLOC */
  "Memory allocation error",

  /* MUTIL_ERR_NULL_POINTER */
  "Passed in a null pointer",

  /* MUTIL_ERR_ILLEGAL_ADDRESS */
  "Illegal memory address",

  /* MUTIL_ERR_ILLEGAL_SIZE */
  "Illegal data size",

  /* MUTIL_ERR_ILLEGAL_TYPE */
  "Illegal data type",

  /* MUTIL_ERR_ILLEGAL_VALUE */
  "Illegal data value",

  /* MUTIL_ERR_OUT_OF_BOUNDS */
  "Access past data bounds",

  /* MUTIL_ERR_DIVIDE_BY_ZERO */
  "Divide by zero",

  /* MUTIL_ERR_INPUT_OUTPUT */
  "Input-output error",

  /* MUTIL_ERR_OVERFLOW */
  "Data overflow",

  /* MUTIL_ERR_INTERRUPT */
  "User interrupted calculation",

  /* MUTIL_ERR_SINGULAR_MATRIX */
  "Singular or nearly singular matrix",

  /* MUTIL_ERR_NOT_CONVERGING */
  "Iterative algorithm is not converging",

  /* MUTIL_ERR_CUM_ROUNDOFF */
  "Cumulative roundoff errors",

  /* MUTIL_ERR_SINGULARITY */
  "Singularity encountered in data",

  /* MUTIL_ERR_TREE_STRUCTURE */
  "Tree-structured algorithm has error in structure",

  /* MUTIL_ERR_ZERO_NEIGHBORS_FOUND */
  "No neighbors found in a nearest neighbor search",

  /* MUTIL_ERR_FEATURE_NOT_IMPLEMENTED */
  "Feature has not been implemented"
};

/* Function documented in ut_err.h */
/* This function just grabs error messages from the static array above.*/
/* Written by Jennifer Hodgdon */
mutil_errcode mutil_err_string(mutil_errcode msg_code, const char **message)
{
  MUTIL_ASSERT( msg_code >= 0 && msg_code <= 18 );

  /* avoid lint warning */
  (void) whatssi;

  *message = error_messages[msg_code];
  return MUTIL_ERR_OK;
}
