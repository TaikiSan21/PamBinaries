
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/ut_alloc.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";

/* This is a self-documenting doc++ file. */

/*
   This file contains the R wrapper implementations
   of the functions in ut_alloc.c, for C programs running from within
   R that use the R dynamic memory allocation functions.
*/

#include <string.h>

#include "ut_alloc.h"
#include "ut_math.h"

#include "ut_RS.h"
#include "ut_debug.h"
#include "ut_limit.h"

/* Memory allocation.              */
/* Note: zeroing also occurs       */
/*                                 */
/* Documented in ut_alloc.h        */
/* Written by William Constantine  */
mutil_errcode mutil_malloc( sint32 size, void **data )
{

  MUTIL_TRACE("Start mutil_malloc()");

  /* avoid lint warning */
  (void) whatssi;

  if( !data ){
    MUTIL_ERROR("Must pass in a non-null pointer for allocated data");
    return MUTIL_ERR_NULL_POINTER;
  }

  if( size <= 0 || size >= MUTIL_SIZE_T_MAX ){
    MUTIL_ERROR( "Size to allocate must be positive and no more than MUTIL_SIZE_T_MAX" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* S_alloc's memory cleanup is controlled by R */
  /* *data  = (void *) S_alloc( size, 1 );      */
  //*data = (void *) R_chk_calloc( (size_t) size, (size_t) 1 );
  *data = (void *) (Calloc(size, int));
  if( !*data ){
    MUTIL_ERROR("Calloc failed");
    return MUTIL_ERR_MEM_ALLOC;
  }

  MUTIL_TRACE("Done with mutil_malloc()");
  return MUTIL_ERR_OK;
}

/* Memory allocation with zeroing.       */
/*                                       */
/* Documented in ut_alloc.h              */
/* Written by William Constantine        */
mutil_errcode mutil_calloc( sint32 nobj, size_t size, void **data )
{
  mutil_errcode trouble;
  double        tmp1;
  sint32        tmp2;

  MUTIL_TRACE("Start mutil_calloc()");

  /* make sure size is OK */
  tmp1 = (double) nobj * (double) size;
  tmp2 = (sint32) nobj * (sint32) size;
  if( tmp1 != (double) tmp2 ){
    MUTIL_ERROR( "Unable to allocate -- size too large" );
    return MUTIL_ERR_MEM_ALLOC;
  }

  trouble = mutil_malloc( tmp2, data );
  if( trouble ){
    return trouble;
  }

  MUTIL_TRACE( "Done with mutil_calloc()" );
  return MUTIL_ERR_OK;
}

/* Memory re-allocation.          */
/*                                */
/* Documented in ut_alloc.h       */
/* Written by William Constantine */
mutil_errcode mutil_realloc( void **data, sint32 new_size, sint32 old_size )
{
  /* avoid lint warning */
  (void) whatssi;
  (void) old_size;

  MUTIL_TRACE("Start mutil_realloc()");

  if( new_size <= 0 || new_size >= MUTIL_SIZE_T_MAX ){
    MUTIL_ERROR( "Size to allocate must be positive and no more than MUTIL_SIZE_T_MAX" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if( !data ){
    MUTIL_ERROR("Must pass in a non-null pointer for allocated data");
    return MUTIL_ERR_NULL_POINTER;
  }

  /* S_realloc's memory cleanup is controlled by R                                      */
  /* *data = (void *) S_realloc( (char *) *data, (long) new_size, (long) old_size, 1 ); */
//  *data = (void *) R_chk_realloc( (void *) *data, (size_t) new_size );

//  #define Realloc(p,n,t) (t *) R_chk_realloc( (void *)(p), (size_t)((n) * sizeof(t)) )
  *data = (void *) (Realloc( *data, new_size, int ));

  if( !*data ){
    MUTIL_ERROR("Realloc failed");
    return MUTIL_ERR_MEM_ALLOC;
  }

  MUTIL_TRACE("Done with mutil_realloc()");
  return MUTIL_ERR_OK;
}

/* Function documented in ut_alloc.h    */
/* This version uses R's Free function  */
/* to free the data                     */
/* Written by William Constantine       */
mutil_errcode mutil_free( void *data, sint32 old_size )
{
  char *pdata;

  /* avoid lint warning */
  (void) whatssi;
  (void) old_size;

  MUTIL_TRACE("Start mutil_free()");

  if( !data ){
    MUTIL_WARN("Attempt to free NULL pointer");
    return MUTIL_ERR_OK;
  }

  pdata = (char *) data;

  Free( pdata );

  MUTIL_TRACE("Done with mutil_free()");
  return MUTIL_ERR_OK;
}
