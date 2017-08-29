
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/ut_rand.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";

/* This is a self-documenting doc++ file. */

/*
   This file contains the R wrapper implementations
   of the functions in ut_rand.c, for C programs 
   running from within R that use the R random 
   number generation functions.
*/

#include "ut_rand.h"
#include "ut_RS.h"
#include "ut_debug.h"
#include "R_ext/Random.h"
#include <time.h>

/* Initializes R random number generation. */
/* The rand_ptr is not used in R random    */
/* number generation.                      */
/*                                         */
/* Documented in ut_rand.h                 */
/* Written by William Constantine          */
mutil_errcode mutil_rand_begin( void **rand_ptr )
{
  MUTIL_TRACE("Start mutil_rand_begin()");

  /* avoid lint warning */
  (void) whatssi;
  (void) rand_ptr;

  GetRNGstate();

  MUTIL_TRACE("Done with mutil_rand_begin()");
  return MUTIL_ERR_OK;
}

/* Terminates R random number generation.  */
/* The rand_ptr is not used in R random    */
/* number generation.                      */
/*                                         */
/* Documented in ut_rand.h                 */
/* Written by William Constantine          */
mutil_errcode mutil_rand_end( void *rand_ptr )
{
  MUTIL_TRACE("Start mutil_rand_end()");

  /* avoid lint warning */
  (void) rand_ptr;

  PutRNGstate();

  MUTIL_TRACE("Done with mutil_rand_end()");
  return MUTIL_ERR_OK;
}

/* The random number generator is private to R;      */
/* there is no way to select the kind of RNG or      */
/* set the seed except by evaluating calls to        */
/* the R functions. Therefore, mutil_rand_set_seed() */
/* does nothing.                                     */
/*                                                   */
/* Documented in ut_rand.h                           */
/* Written by William Constantine                    */
mutil_errcode mutil_rand_set_seed( const void *seed, void *rand_ptr )
{
  MUTIL_TRACE("Start mutil_rand_set_seed()");

  /* avoid lint warning */
  (void) rand_ptr;

  MUTIL_TRACE("Done with mutil_rand_set_seed()");
  return MUTIL_ERR_OK;
}

/* Generate a uniformly distributed   */
/* pseudo-random deviate.             */
/*                                    */
/* Documented in ut_rand.h            */
/* Written by William Constantine     */
mutil_errcode mutil_rand_uniform( void *rand_ptr, double *num_out )
{
  MUTIL_TRACE("Start mutil_rand_uniform()");

  /* avoid lint warning */
  (void) rand_ptr;

  if( !num_out ) {
    MUTIL_ERROR( "NULL pointer for output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do not allow function to return endpoints */

  do {
    *num_out = unif_rand();
  } while( *num_out == 0 || *num_out == 1 );

  MUTIL_TRACE("Done with mutil_rand_uniform()");
  return MUTIL_ERR_OK;
}

/* Generate a normally distributed    */
/* pseudo-random deviate.             */
/*                                    */
/* Documented in ut_rand.h            */
/* Written by William Constantine     */
mutil_errcode mutil_rand_normal( void *rand_ptr, double *num_out )
{
  MUTIL_TRACE("Start mutil_rand_normal()");

  /* avoid lint warning */
  (void) rand_ptr;

  if( !num_out ) {
    MUTIL_ERROR( "NULL pointer for output" );
    return MUTIL_ERR_NULL_POINTER;
  }

  *num_out = norm_rand();

  MUTIL_TRACE("Done with mutil_rand_normal()");
  return MUTIL_ERR_OK;
}


