
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_intrp.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_INTRP_H
#define IN_UT_INTRP_H

#include "ut_plat.h"
#include "ut_type.h"

/*
   This file contains a function declaration and macros for
   user-interrupt handling.  The function is user-definable to
   let different applications have different interrupt handling.

   The INTERRUPT\_ENABLE compile flag must be defined for interrupt
   handling to be used.
*/

#ifdef __cplusplus
extern "C" {
#endif

/** Check for user interrupt.
 * Check to see if the user has requested that processing be interrupted.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it. Note that there is
 * no user interrupt checking in the wrap\_std library -- this function
 * always returns FALSE.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrp.h
 * @source ut\_intrp.c
 * @library wrap
 * @usage  #if( nops > nextcheck ) stop_now = mutil_interrupt( nops, &nextcheck, intrp_ptr );#
 * @return  TRUE if user requested interrupt, and FALSE if not.
 * @param  nops  Approximate number of arithmetic operations that
 *    have occurred in this function.
 * @param  next_check  Number of operations for next time to check
 *    for interrupt; if a check occurs, this will be incremented
 *    by the interrupt threshold.
 * @param  ptr   A pointer that is passed into functions that are
 *    interrupt enabled, for implementation-specific use.
 * @see MUTIL_INTERRUPT
 * @see mutil_interrupt_threshold
 * @see Interrupt Handling
 * @see Compile-Time Options
 */
MUTIL_WRAPEXPORT boolean mutil_interrupt( double nops, double *next_check, void *ptr );


/** Threshold for user interrupts.
 * Global variable storing the approximate number of arithmetic
 * operations to do in between checking for user interrupts.
 *
 * This variable is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrp.h
 * @source ut\_intrp.c
 * @library wrap
 * @see MUTIL_INTERRUPT
 * @see MUTIL_INTERRUPT_INIT
 */
MUTIL_WRAPEXPORT double mutil_interrupt_threshold;


#ifdef INTERRUPT_ENABLE


/** Check for user interrupt.
 * Check to see if the user has requested that processing be interrupted,
 * if the number of arithmetic operations exceeds the threshold.
 * All interrupt checking in the mutils library goes through this
 * macro.
 *
 * If the INTERRUPT\_ENABLE compile flag is not used, this macro
 * is defined to be FALSE; otherwise, it calls the \Ref{mutil_interrupt}
 * function to check whether the user has requested an interrupt.
 * The macro must be used in conjunction with the
 * \Ref{MUTIL_INTERRUPT_INIT} macro (otherwise, it will not compile).
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrp.h
 * @source ut\_intrp.h
 * @library wrap
 * @usage  #if(MUTIL_INTERRUPT(i, intrp_ptr)) return MUTIL_ERR_INTERRUPT;#
 * @return TRUE if user requested interrupt, and FALSE if not.
 * @param  NUM_OPS  Approximate number of arithmetic operations that
 *            have occurred in this function.
 * @param  INTRP_PTR A pointer that is passed into functions that are
 *            interrupt enabled, for implementation-specific use.
 * @see mutil_interrupt
 * @see MUTIL_INTERRUPT_INIT
 * @see Interrupt Handling
 * @see Compile-Time Options
 */
#define MUTIL_INTERRUPT(NUM_OPS, INTRP_PTR) (\
         ((NUM_OPS) > interrupt_next_check ) ? \
            mutil_interrupt((NUM_OPS), &interrupt_next_check, (INTRP_PTR)) \
            : FALSE )


/** Initialize user interrupts.
 * Set up a variable local to a function that allows the
 * \Ref{MUTIL_INTERRUPT} macro to work correctly.  This macro must
 * be put in the variable declaration section of any function that needs
 * to check for user interrupts.
 *
 * If the INTERRUPT\_ENABLE compile flag is not used, this macro
 * is defined to be a null operation.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_intrp.h
 * @source ut\_intrp.h
 * @library wrap
 * @usage  #MUTIL_INTERRUPT_INIT(intrp_ptr);#
 * @param  INTRP_PTR  A pointer that is passed into functions that are
 *            interrupt enabled, for implementation-specific use.
 * @see MUTIL_INTERRUPT
 */
#define MUTIL_INTERRUPT_INIT(INTRP_PTR)  double interrupt_next_check = 0


#else /* ifdef INTERRUPT_ENABLE */

  /* non-interrupt-enabled versions of above macros -- see doc above */

/* this define will get rid of compile warnings about unused variable */
#define MUTIL_INTERRUPT_INIT(INTRP_PTR)  (void) (INTRP_PTR)
/* this define says user never interrupts */
#define MUTIL_INTERRUPT(NUM_OPS, INTRP_PTR) FALSE

#endif /* ifdef INTERRUPT_ENABLE */

#ifdef __cplusplus
}
#endif

#endif /* ifndef IN_UT_INTRP_H */
