
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/ut_intrp.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "ut_intrp.h"
#include "R.h"
#include "R_ext/Utils.h"

/* This file contains the R definitions of the
   interrupt function and variable from ut_intrp.h */

/* Interrupt checking.             */
/*                                 */
/* Documented in ut_intrp.h        */
/* Written by William Constantine  */
boolean mutil_interrupt( double nops, double *next_check, void *ptr )
{
  /* avoid lint warning */
  (void) whatssi;
  (void) ptr;

  R_CheckUserInterrupt();

  if( next_check && ( nops > *next_check )) {
    *next_check += mutil_interrupt_threshold;
  }

  return FALSE;
}

/* Million floating-point operations. */
/*                                    */
/* Documented in ut_intrp.h           */
/* Written by William Constantine     */
double mutil_interrupt_threshold = 1000000;
