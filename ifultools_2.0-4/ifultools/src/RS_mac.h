
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_mac.h $: $Revision: #1 $, $Date: 2008/03/21 $ */
/* This is a self-documenting doc++ file. */

#ifndef IN_RS_MAC_H_
#define IN_RS_MAC_H_

#include "ut_plat.h"
#include "ut_type.h"
#include "mat_type.h"
#include "ut_err.h"

/* This file contains function declarations the mutils fractal library. */

#ifdef __cplusplus
extern "C" {
#endif

#define MEMLIST_FREE_ON_ERROR_SPLUS( ERR_CODE, MEM_LIST, PROBLEM_STRING ) \
  if ( ERR_CODE ){ \
    MUTIL_FREE_WARN( memlist, MEM_LIST ); \
    PROBLEM #PROBLEM_STRING ERROR; \
  }

#define MEMLIST_FREE_ON_ERROR_REGISTER( ERR_CODE, MEM_LIST ) \
  if ( ERR_CODE ){ \
    MUTIL_FREE_WARN( memlist, MEM_LIST ); \
    PROBLEM "Unable to register memory with the memory manager" ERROR; \
  }

#define READ_MATSET_REGISTER( RS_POINTER, TYPE, MUTILS_PTR ) \
  err = matset_from_R( RS_POINTER, TYPE, MUTILS_PTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert (RS_POINTER) R or S-PLUS list to MUTILS matrix set (MUTILS_PTR) ); \
  err = memlist_member_register( &list, MUTILS_PTR, MEMTYPE_MATSET ); \
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list )

#define READ_MATRIX_REGISTER( RS_POINTER, MUTILS_PTR ) \
  err = mutil_R_type( RS_POINTER, &type ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to read (RS_POINTER) R or S-PLUS matrix type ); \
  err = matuniv_from_R( RS_POINTER, type, MUTILS_PTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert R or S-PLUS matrix (RS_POINTER) to an MUTILS matrix (MUTILS_PTR) ); \
  err = memlist_member_register( &list, MUTILS_PTR, MEMTYPE_MATUNIV ); \
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list )

#define CONVERT_MATRIX_AND_RETURN( MUTILS_FUNCTION, INMAT_PTR, OUTMAT_PTR ) \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Problem calling MUTILS_FUNCTION function ); \
  err = memlist_member_register( &list, INMAT_PTR, MEMTYPE_MATUNIV ); \
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list ); \
  err = matuniv_to_R( INMAT_PTR, (mutil_R_class_type) MUTIL_R_MATRIX, OUTMAT_PTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert output data to R format ); \
  MUTIL_FREE_WARN( memlist, &list ); \
  return *(OUTMAT_PTR)

#define CONVERT_MATSET_AND_RETURN( MUTILS_FUNCTION, INMATSET_PTR, OUTLIST_PTR ) \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Problem calling MUTILS_FUNCTION function ); \
  err = memlist_member_register( &list, INMATSET_PTR, MEMTYPE_MATSET ); \
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list ); \
  err = matset_to_R_list( INMATSET_PTR, OUTLIST_PTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert output data to R format ); \
  MUTIL_FREE_WARN( memlist, &list ); \
  return *(OUTLIST_PTR)

#define CHECK_CONVERSION( TYPE, IN, OUT ) \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert TYPE type argument IN to OUT )

#define SINT32_FROM_R( IN_PTR, OUTPTR ) \
  err = sint32_from_R( IN_PTR, OUTPTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert sint32 type argument IN_PTR to OUTPTR )

#define DOUBLE_FROM_R( IN_PTR, OUTPTR ) \
  err = double_from_R( IN_PTR, OUTPTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert double type argument IN_PTR to OUTPTR )

#define BOOLEAN_FROM_R( IN_PTR, OUTPTR ) \
  err = boolean_from_R( IN_PTR, OUTPTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert boolean type argument IN_PTR to OUTPTR )

#define DISTANCE_METRIC_FROM_R( IN_PTR, OUTPTR ) \
  err = fra_distance_metric_from_R( IN_PTR, OUTPTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert fra_distance_metric type argument IN_PTR to OUTPTR )

#define WAV_FDP_ESTIMATOR_FROM_R( IN_PTR, OUTPTR ) \
  err = wav_fdp_estimator_from_R( IN_PTR, OUTPTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert wav_fdp_estimator type argument IN_PTR to OUTPTR )

#define WAV_FDP_ESTIMATOR_FROM_R( IN_PTR, OUTPTR ) \
  err = wav_fdp_estimator_from_R( IN_PTR, OUTPTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert wav_fdp_estimator type argument IN_PTR to OUTPTR )

#define WAV_FILTER_TYPE_FROM_R( IN_PTR, OUTPTR ) \
  err = wav_filter_type_from_R( IN_PTR, OUTPTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert wav_filter_type type argument IN_PTR to OUTPTR )

#define WAV_TRANSFORM_FROM_R( IN_PTR, OUTPTR ) \
  err = wav_transform_from_R( IN_PTR, OUTPTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert wav_transform type argument IN_PTR to OUTPTR )

#define WAV_TRANSFORM_PEAK_FROM_R( IN_PTR, OUTPTR ) \
  err = wav_transform_peak_from_R( IN_PTR, OUTPTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert wav_transform_peak type argument IN_PTR to OUTPTR )

#define WAV_SHRINK_THRESHOLD_FROM_R( IN_PTR, OUTPTR ) \
  err = wav_shrink_threshold_from_R( IN_PTR, OUTPTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert wav_shrink_threshold type argument IN_PTR to OUTPTR )

#define WAV_SHRINK_FUNCTION_FROM_R( IN_PTR, OUTPTR ) \
  err = wav_shrink_function_from_R( IN_PTR, OUTPTR ); \
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, Unable to convert wav_shrink_function type argument IN_PTR to OUTPTR )


#ifdef __cplusplus
}
#endif

#endif /* IN_RS_MAC_H_ */
