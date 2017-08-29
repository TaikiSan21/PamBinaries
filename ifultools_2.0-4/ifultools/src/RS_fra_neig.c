
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/RS_fra_neig.c $: $Revision: #1 $, $Date: 2008/03/21 $ ";

/* This is a self-documenting doc++ file. */

/* This file contains wrapper functions, callable from R
   for the MUTILS fractal library.

   Functions wrapped:

   frauniv_neighbor_find()
   frauniv_neighbor_find_arbitrary()
*/

#include "fra_neig.h"
#include "ut_RS.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_mem.h"
#include "RS_mac.h"

/** Find the nearest neighbors of a multidimensional embedding of data in the phase space.
 * @source RS\_fra\_neig.c
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @usage #itCall( "RS_fractal_neighbor_find", embedding, n.neighbor, distance.max, distance.metric, pdmatrix, sort.distances, orbital.lag))#
 * @return                 An R ... containing ...
 * @param embedding        Pointer to an R object containing ... embedding
 * @param n.neighbor       Pointer to an R object containing ... n.neighbor
 * @param distance.max     Pointer to an R object containing ... distance.max
 * @param distance.metric  Pointer to an R object containing ... distance.metric
 * @param pdmatrix         Pointer to an R object containing ... pdmatrix
 * @param sort.distances   Pointer to an R object containing ... sort.distances
 * @param orbital.lag      Pointer to an R object containing ... orbital.lag
 * @see _fra_distance_metric
 * @see frauniv_embed_neighbors
*/
EXTERN_R SEXP RS_fractal_neighbor_find(
 SEXP pr_embedding,
 SEXP pr_n_neighbor,
 SEXP pr_distance_max,
 SEXP pr_distance_metric,
 SEXP pr_sort_distances,
 SEXP pr_orbital_lag )
{
  SEXP                  pr_ret_neighbor_distances;  
  SEXP                  pr_ret_neighbor_indices;    
  SEXP                  pr_ret_obj;                 
  SEXP                  pr_ret_original_indices;    
  boolean               sort_distances;             
  double                distance_max;               
  fra_distance_metric   distance_metric;            
  mutil_data_type       type;                       
  mutil_errcode         err;                        
  sint32                n_neighbor;                 
  sint32                orbital_lag;                
  univ_mat              embedding;                  
  univ_mat              neighbor_distances;         
  univ_mat              neighbor_indices;           
  univ_mat              original_indices;           
  void                  *VPNULL = NULL;             
  memlist               list;

  /* Avoid lint warning */
  (void) whatssi;

  /* initialize memory list */
  MEMLIST_INIT( list );

  /* Conversion of input data ... */

  /* ... pr_time_series to time_series */
  READ_MATRIX_REGISTER( pr_embedding, &embedding );

  /* ... pr_n_neighbor to n_neighbor */
  SINT32_FROM_R( pr_n_neighbor, &n_neighbor );

  /* ... pr_distance_max to distance_max */
  DOUBLE_FROM_R( pr_distance_max, &distance_max );

  /* ... pr_distance_metric to distance_metric */
  DISTANCE_METRIC_FROM_R( pr_distance_metric, &distance_metric );

  /* ... pr_sort_distances to sort_distances */
  BOOLEAN_FROM_R( pr_sort_distances, &sort_distances );

  /* ... pr_orbital_lag to orbital_lag */
  SINT32_FROM_R( pr_orbital_lag, &orbital_lag );

  /* Call the function */
  err = frauniv_neighbor_find(
    &embedding,
    n_neighbor,
    distance_max,
    distance_metric,
    (univ_mat*) NULL,
    sort_distances,
    orbital_lag,
    VPNULL,
    &original_indices,
    &neighbor_indices,
    &neighbor_distances );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Problem calling frauniv_neighbor_find() function" );
  err = memlist_member_register( &list, &original_indices, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
  err = memlist_member_register( &list, &neighbor_indices, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );
  err = memlist_member_register( &list, &neighbor_distances, MEMTYPE_MATUNIV);
  MEMLIST_FREE_ON_ERROR_REGISTER( err, &list );

  /* create the output R object */
  err = matuniv_to_R( &original_indices, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_original_indices );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );
  err = matuniv_to_R( &neighbor_indices, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_neighbor_indices );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );
  err = matuniv_to_R( &neighbor_distances, (mutil_R_class_type) MUTIL_R_MATRIX, &pr_ret_neighbor_distances );
  MEMLIST_FREE_ON_ERROR_SPLUS( err, &list, "Unable to convert output data to R format" );

  PROTECT( pr_ret_obj = allocVector( VECSXP, 3 ) );
  SET_VECTOR_ELT( pr_ret_obj, 0, pr_ret_original_indices );
  SET_VECTOR_ELT( pr_ret_obj, 1, pr_ret_neighbor_indices );
  SET_VECTOR_ELT( pr_ret_obj, 2, pr_ret_neighbor_distances );
  UNPROTECT(1);
  
  /* free registered local memory */
  MUTIL_FREE_WARN( memlist, &list );

  return pr_ret_obj;
}
