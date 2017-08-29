
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_rand.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_RAND_H
#define IN_UT_RAND_H

#include "ut_plat.h"
#include "ut_err.h"
#include "ut_type.h"

/*
   This file contains function declarations for random number
   generation.  The functions are user-definable to let different
   applications have different random number generation, and all
   random number generation in the mutils library goes through these
   functions.
*/

#ifdef __cplusplus
extern "C" {
#endif


/** Initialize random number generation.
 * This function is called to initialize the random number generator
 * pointer.  The action of this function depends on the specific
 * wrapper implementation of the random number generator, but typically
 * it would allocate space for storage of the random number seed
 * and any other information the random number generator needs, and
 * initialize the seed to a legal value (but not necessarily the same
 * value every time).
 *
 * This function would typically be called once per application using
 * random numbers.  After all usage of random numbers is finished,
 * the \Ref{mutil_rand_end} function must be called.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_rand.h
 * @source ut\_rand.c
 * @library wrap
 * @usage  #err_code = mutil_rand_begin(&rand_ptr);#
 * @return     Standard mutils error/OK code.
 * @param  rand_ptr    Pointer for returning a new pointer for
 *    implementing random number generator.
 * @see mutil_rand_end
 * @see mutil_rand_set_seed
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_rand_begin( void **rand_ptr );


/** End random number generation.
 * This function is called to perform any tasks necessary to
 * end random number generation.
 * The action of this function depends on the specific
 * wrapper implementation of the random number generator, but typically
 * it would clean up the space allocated by the \Ref{mutil_rand_begin}
 * function.  After passing a random number generator pointer to this
 * function, the pointer should no longer be used.
 *
 * This function would typically be called once per application using
 * random numbers, after all usage of random numbers is finished.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_rand.h
 * @source ut\_rand.c
 * @library wrap
 * @usage  #err_code = mutil_rand_end(rand_ptr);#
 * @return     Standard mutils error/OK code.
 * @param  rand_ptr    Initialized pointer for implementing random
 *    number generator.
 * @see mutil_rand_begin
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_rand_end( void *rand_ptr );


/** Set random seed.
 * Set the seed value for the random number generator using the
 * given input value.  The specific implementation of the random
 * number generator will determine how the input seed value is used,
 * but any implementation should ensure that a given seed value
 * passed in via this function always results in the same sequence
 * of random numbers.  Typically, the values pointed to by the
 * random number generator pointer will be changed by this function.
 *
 * The standard wrapper implementation of this function takes a
 * pointer to a signed long integer as input, and uses the long integer
 * as the seed; it must be positive.
 * The S-PLUS wrapper implementation takes a pointer
 * to a signed long integer as input, and uses the long integer
 * as the seed; it should be between 0 and 1023.
 *
 * If the input for the seed is a NULL pointer, this function will
 * set the seed to a pseudo-randomized number, based on the date and time of
 * day or some other somewhat random input, depending on the implementation.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_rand.h
 * @source ut\_rand.c
 * @library wrap
 * @usage  #err_code = mutil_rand_set_seed(&seed_val, rand_ptr);#
 * @return     Standard mutils error/OK code.
 * @param  seed        Input value for seeding random number generator.
 * @param  rand_ptr    Initialized pointer for implementing random
 *    number generator.
 * @see mutil_rand_begin
 * @see mutil_rand_uniform
 * @see mutil_rand_normal
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_rand_set_seed( const void *seed,
  void *rand_ptr );


/** Generate a uniformly-distributed random number.
 * Generate a single random number from a distribution that
 * is uniform between 0.0 and 1.0, not including the endpoints.
 * Typically, the values pointed to by the
 * random number generator pointer will be changed by this function.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_rand.h
 * @source ut\_rand.c
 * @library wrap
 * @usage  #err_code = mutil_rand_uniform(rand_ptr, &val);#
 * @return     Standard mutils error/OK code.
 * @param  rand_ptr    Initialized pointer for implementing random
 *    number generator.
 * @param  num_out     Pointer for returning random number.
 * @see mutil_rand_normal
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_rand_uniform( void *rand_ptr,
  double *num_out );


/** Generate a normally-distributed random number.
 * Generate a single random number from a normal distribution that
 * has standard deviation 1.0 and mean 0.0.
 * Typically, the values pointed to by the
 * random number generator pointer will be changed by this function.
 *
 * To convert the numbers returned to have a normal distribution
 * with different mean and standard deviation, multiply the numbers
 * output by this function by the desired standard deviation and then
 * add the desired mean.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_rand.h
 * @source ut\_rand.c
 * @library wrap
 * @usage  #err_code = mutil_rand_normal(rand_ptr, &val);#
 * @return     Standard mutils error/OK code.
 * @param  rand_ptr    Initialized pointer for implementing random
 *    number generator.
 * @param  num_out     Pointer for returning random number.
 * @see mutil_rand_uniform
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_rand_normal( void *rand_ptr,
  double *num_out );


#ifdef __cplusplus
}
#endif

#endif /* ifndef IN_UT_RAND_H*/
