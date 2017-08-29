
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/ut_debug.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";

/* This is a self-documenting doc++ file. */

#include "ut_RS.h"
#include "ut_debug.h"
#include <stdio.h>

/*
   This file contains the R wrapper implementations
   of the functions in ut_debug.h, for C programs running from
   within R.
*/

#define MAX_MESSAGE_ELEMENT 2


/* default beginning values */

static mutil_msg_level print_level        = MUTIL_MSG_LEVEL_ERROR;
static boolean         print_line_numbers = FALSE;

/* Function documented in ut_debug.h */
/* This version prints error messages to stdout. */
/* Written by William Constantine */
mutil_errcode mutil_msg_print(mutil_msg_level level, const char *message,
  const char *filename, int linenum)
{
  char          linstr[50];

  if( level == MUTIL_MSG_LEVEL_NONE ) {
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if( level < print_level ) {
    return MUTIL_ERR_OK;
  }

  if( !message || !filename ) {
    return MUTIL_ERR_NULL_POINTER;
  }

  /* figure out the line number string */
  /* note that we put file name along with it */
  if( print_line_numbers ) {
    if( sprintf( linstr, "%s(%d)", filename, linenum ) < 0 ) {
      return MUTIL_ERR_INPUT_OUTPUT;
    }
  }
  else if( sprintf(linstr, "%s", "") < 0 ) {
    return MUTIL_ERR_INPUT_OUTPUT;
  }

  /* print the message */

  if( level < MUTIL_MSG_LEVEL_WARN ) {
/* R does not have an equivalent PRINT_IT macro ... */
/*      PROBLEM "%s: %s", linstr, message PRINT_IT; */
      PROBLEM "%s: %s", linstr, message WARN;
  }
  else {
      PROBLEM "%s: %s", linstr, message WARN;
  }

  return MUTIL_ERR_OK;
}


/* Function documented in ut_debug.h */
/* This version sets the local static message level variable */
/* Written by William Constantine */
void mutil_msg_set_level( mutil_msg_level level )
{
  /* avoid lint warning */
  (void) whatssi;

  print_level = level;
}


/* Function documented in ut_debug.h */
/* Written by William Constantine */
void mutil_msg_print_line( boolean show_nums )
{
  print_line_numbers = show_nums;
}


/* Function documented in ut_debug.h */
/* Written by William Constantine */
void mutil_msg_print_time( boolean show_time )
{
  (void) show_time;
}


/* Function documented in ut_debug.h */
/* This version calls R error on fatal errors, which causes
   the function to crash back to the R command prompt. */
/* Written by William Constantine */
void mutil_abort()
{
  PROBLEM "Aborting function" ERROR;
}
