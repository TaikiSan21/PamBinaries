
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/mat_rand.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */

#include "mat_assn.h"
#include "mat_rand.h"
#include "mat_sort.h"
#include "mat_stat.h"
#include "mat_type.h"

#include "ut_math.h"
#include "ut_mem.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_rand.h"

#include "ut_mem.h"

#include <time.h>

/* This file contains function definitions
   for the creation of randomized matrices
   and belongs to the mutils matrix library. */

/* Random index generation with   */
/* optional replacement.          */
/*                                */
/* Documented in mat_rand.h       */
/* Written by William Constantine */

 mutil_errcode matuniv_random_uniform_indices(
  const sint32    nrow,
  const sint32    ncol,
  const boolean   replacement,
  void           *rand_ptr,
  void           *intrp_ptr,
  univ_mat       *result )
{
  mutil_errcode  err;

  MUTIL_TRACE( "Start matuniv_random_uniform_indices()" );

  /* avoid lint warning */
  (void) whatssi;

  err = mats32_random_uniform_indices( nrow, ncol, replacement, rand_ptr,
    intrp_ptr, &( result->mat.s32mat) );
  if ( err ){
    return err;
  }

  /* force the type to be sint32 */

  result->type = MUTIL_SINT32;

  MUTIL_TRACE( "Done with matuniv_random_uniform_indices()" );

  return MUTIL_ERR_OK;
}

mutil_errcode mats32_random_uniform_indices(
  const sint32    nrow,
  const sint32    ncol,
  const boolean   replacement,
  void           *rand_ptr,
  void           *intrp_ptr,
  sint32_mat     *result )
{
  double         deviate;
  double        *pd_rand;
  mutil_errcode  err;
  sint32         i;
  sint32         nelem = nrow * ncol;
  sint32        *ps_result;
  double_mat     random;
  memlist        list;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start mats32_random_uniform_indices()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* allocate memory */

  if ( !replacement ){
    err = matdbl_malloc_register( &random, nrow, ncol, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  err = mats32_malloc_register( result, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* assign pointers */

  pd_rand   = random.data;
  ps_result = result->data;

  /* create the random deviates */

  for ( i = 0; i < nelem; i++ ){

    /* generate a uniformly distributed random deviate */

    err = mutil_rand_uniform( rand_ptr, &deviate );
    if ( err ) {
      MUTIL_ERROR( "Could not generate a uniformly-distributed "
        "random number" );
      return err;
    }

    if ( replacement ){

      *ps_result = (sint32) floor( deviate * (double) nelem );
      ps_result++;
    }
    else{

      *pd_rand = deviate;
      pd_rand++;
    }

    if ( MUTIL_INTERRUPT( (double) nelem, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }
  }

  if ( !replacement ){

    err = matdbl_sort_index_partial( &random, (sint32_mat *) NULL, intrp_ptr, result );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* free nodes corresponding to registered
     memory for the result, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free other malloced space and
     corresponding nodes in memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with mats32_random_uniform_indices()" );

  return MUTIL_ERR_OK;
}

