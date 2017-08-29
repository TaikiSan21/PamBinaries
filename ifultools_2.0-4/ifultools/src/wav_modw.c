
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/wav_modw.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "wav_modw.h"

#include "mat_assn.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_limit.h"
#include "ut_math.h"
#include "ut_mem.h"

/* Static functions declared here and defined at end of file */

static mutil_errcode localfn_modw_input_check(
   const univ_mat *time_series,
   const mat_set  *filters,
   sint32          n_level);

static mutil_errcode localfn_imodwt_input_check(
   const mat_set *modwt,
   const mat_set *filters );

static mutil_errcode localfn_filters_check(
   const mat_set *filters );

#undef LOCALDEF_CHECK_NULL_POINTER_MODW
#define LOCALDEF_CHECK_NULL_POINTER_MODW( DATA_PTR, DATA_TYPE, \
                                     TYPE_PREFIX )             \
   err = TYPE_PREFIX ## _validate( DATA_PTR );                 \
   if ( err ) return err;                                      \
   if ( DATA_PTR == ( DATA_TYPE * ) NULL ){                    \
     MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" ); \
     return MUTIL_ERR_NULL_POINTER;                            \
   }

#undef LOCALDEF_ILOG2
#define LOCALDEF_ILOG2( VALUE ) \
(sint32) floor( MUTIL_LOG2( (double) ( VALUE ) + MUTIL_DOUBLE_EPSILON ) )

#undef LOCALDEF_IS_EVEN
#define LOCALDEF_IS_EVEN(N) ( (N % 2) == 0 ? 1 : 0 )

/* The maximum overlap discrete wavelet transform function */
/* Documented in wav_modw.h                                */
/* Written by William Constantine                          */

mutil_errcode wavuniv_transform_maximum_overlap(
   const univ_mat *time_series,
   const mat_set  *filters,
   sint32          n_level,
   void           *intrp_ptr,
   mat_set        *result )
{
   double         sum_scaling;
   double         sum_wavelet;
   double        *pd_new_scaling_coeffs;
   double        *pd_old_scaling_coeffs;
   double        *pd_result;
   double        *pd_scaling;
   double        *pd_wavelet;
   memlist        list;
   mutil_errcode  err;
   sint32         dims;
   sint32         filter_length;
   sint32         j;
   sint32         jj;
   sint32         l;
   sint32         n_sample;
   sint32         n_scale;
   sint32         t;
   univ_mat       new_scaling_coeffs;
   univ_mat       old_scaling_coeffs;

   /* initialize interrupt pointer */

   MUTIL_INTERRUPT_INIT( intrp_ptr );

   MUTIL_TRACE( "Start wavuniv_transform_maximum_overlap()" );

   /* avoid lint warning */

   ( void ) whatssi;

   /* initialize memory list */

   MEMLIST_INIT( list );

   /* check inputs for errors */

   err = localfn_modw_input_check(
      time_series,
      filters,
      n_level);
   if ( err ) return err;

   /* obtain sizes of matrices */

   filter_length = MATUNIV_NELEM( &filters->mats[ 0 ] );
   n_sample      = MATUNIV_NELEM( time_series );

   /* allocate space for output */

   dims = n_level + 1;
   err = matset_malloc_register( result, 1, &dims, &list );
   MEMLIST_FREE_ON_ERROR( err, &list );

   err = matset_malloc_matrices( result, 1, n_sample, MUTIL_DOUBLE );
   MEMLIST_FREE_ON_ERROR( err, &list );

   /* malloc space for scaling coeffs,
      one for the new and one for the old */

   err = matuniv_malloc_register( &new_scaling_coeffs, 1, n_sample,
                              MUTIL_DOUBLE, &list );
   MEMLIST_FREE_ON_ERROR( err, &list );

   err = matuniv_malloc_register( &old_scaling_coeffs, 1, n_sample,
                              MUTIL_DOUBLE, &list );
   MEMLIST_FREE_ON_ERROR( err, &list );

   /* create pointers */

   pd_new_scaling_coeffs = new_scaling_coeffs.mat.dblmat.data;
   pd_old_scaling_coeffs = old_scaling_coeffs.mat.dblmat.data;
   pd_wavelet            = filters->mats[ 0 ].mat.dblmat.data;
   pd_scaling            = filters->mats[ 1 ].mat.dblmat.data;

   /* copy input time series into "old" scaling vector. here,
      the time series (which can be of any type other than
      complex) is "converted" to type double */

   if ( MATUNIV_NCOL( time_series) == 1 ) {
     for ( t = 0; t < n_sample; t++ ){
       pd_old_scaling_coeffs[ t ] =
         (double) MATUNIV_ELEM( time_series, t, 0 );
     }
   }
   else {
     for ( t = 0; t < n_sample; t++ ){
       pd_old_scaling_coeffs[ t ] =
         (double) MATUNIV_ELEM( time_series, 0, t );
     }
   }

   /* calculate the modwt wavelet and scaling coefficient matrices */

   for( j = 0; j < n_level; j++ ){

     pd_result = result->mats[ j ].mat.dblmat.data;

     n_scale = 1 << j ;

     for ( t = 0; t < n_sample; t++ ){

       sum_wavelet = 0.0;
       sum_scaling = 0.0;

       /* perform convolution operations */

       for ( l = 0; l < filter_length ; l++ ){

         jj = t - n_scale * l;

         while ( jj < 0 ) jj += n_sample;

         sum_wavelet += pd_old_scaling_coeffs[ jj ] * pd_wavelet[ l ];
         sum_scaling += pd_old_scaling_coeffs[ jj ] * pd_scaling[ l ];
       }

       /* store convolution results in wavelet and scaling coefficient matrices */

       pd_result[ t ]             = sum_wavelet;
       pd_new_scaling_coeffs[ t ] = sum_scaling;

     } /* end loop over time index  */

     /* store the new scaling coefficients */

     if ( j == ( n_level - 1 ) ){

       err = matuniv_assign( &new_scaling_coeffs,
                             intrp_ptr,
                             &result->mats[ n_level ] );
     }
     else{

       err = matuniv_assign( &new_scaling_coeffs,
                             intrp_ptr,
                             &old_scaling_coeffs );
     }

     MEMLIST_FREE_ON_ERROR( err, &list );

     /* Check for interrupts */

     if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
       MUTIL_ERROR( "user interrupt" );
       MUTIL_FREE_WARN( memlist, &list );
       return MUTIL_ERR_INTERRUPT;
     }

   }  /* end loop over scale index */

   /* free nodes corresponding to registered
      modwt matrix memory, but do not free
      the memory itself */

   err = memlist_member_unregister( result, &list );
   MEMLIST_FREE_ON_ERROR( err, &list );

   /* free malloced space and corresponding linked list structures */

   MUTIL_FREE_WARN( memlist, &list );

   MUTIL_TRACE( "End wavuniv_transform_maximum_overlap()" );

   return MUTIL_ERR_OK;
}

/* The maximum overlap discrete wavelet packet transform function */
/* Documented in wav_modw.h                                       */
/* Written by William Constantine                                 */
mutil_errcode wavuniv_transform_maximum_overlap_packet(
   const univ_mat *time_series,
   const mat_set  *filters,
   sint32          n_level,
   void           *intrp_ptr,
   mat_set        *result )
{
  double            convolution;
  double           *pd_filter;
  mutil_errcode     err;
  sint32            j;                   /* scale index      */
  sint32            l;                   /* filter index     */
  sint32            n;                   /* local node index */
  sint32            n_child_base_row;
  sint32            n_child_node;
  sint32            n_child_row;
  sint32            n_child_time;
  sint32            n_mod;
  sint32            n_sample;
  sint32            n_nodes;
  sint32            n_parent_base_row;
  sint32            n_parent_node;
  sint32            n_parent_row;
  sint32            n_parent_time;
  sint32            n_scale;
  sint32            filter_length;
  sint32            t;                   /* time index       */
  memlist           list;

  MUTIL_INTERRUPT_INIT( intrp_ptr );
  MEMLIST_INIT( list );

  MUTIL_TRACE( "Start wavuniv_transform_maximum_overlap_packet()" );

  /* check inputs */

  err = localfn_modw_input_check(
    time_series,
    filters,
    n_level);
  if ( err ) return err;

  /* obtain sizes of matrices */

  filter_length = MATUNIV_NELEM( &filters->mats[ 0 ] );
  n_sample      = MATUNIV_NELEM( time_series );
  /* n_nodes       = (sint32) MUTIL_POW( 2, n_level + 1 ) - 1; */
  /* nlevel is positive */
  n_nodes       = ( 2 << n_level ) - 1;

  /* allocate space for output */

  err = matset_malloc_register( result, 1, &n_nodes, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices( result, 1, n_sample, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* Copy time series into first matrix of MODWPT
     matrix set. This represents node Wmod_{0,0}. */

  if ( MATUNIV_NCOL( time_series) == 1 ) {
    for ( t = 0; t < n_sample; t++ ){
      result->mats[ 0 ].mat.dblmat.data[ t ] =
        (double) MATUNIV_ELEM( time_series, t, 0 );
    }
  }
  else {
    for ( t = 0; t < n_sample; t++ ){
      result->mats[ 0 ].mat.dblmat.data[ t ] =
        (double) MATUNIV_ELEM( time_series, 0, t );
    }
  }

  /* loop over each scale */

  for ( j = 1; j <= n_level; j++ ){

    n_scale           = 1 << ( j - 1 );
    n_child_base_row  = ( 1 << j ) - 1;
    n_parent_base_row = n_scale - 1;

    /* loop over each local node */

    for ( n = 0; n < (sint32) ( 1 << j ); n++ ){

      /* Select appropriate filter */

      n_mod = n % 4;

      if ( ( n_mod == 0 ) || ( n_mod == 3 ) ){
        pd_filter = filters->mats[ 1 ].mat.dblmat.data;
      }
      else{
        pd_filter = filters->mats[ 0 ].mat.dblmat.data;
      }

      /* set indices used to point to parent coefficients */

      n_parent_node     = n >> 1;
      n_parent_row      = n_parent_base_row + n_parent_node;

      /* set indices used to point to child coefficients */

      n_child_node     = n;
      n_child_row      = n_child_base_row + n_child_node;

      /* loop over time indices */

      for ( t = 0; t < n_sample; t++ ) {

        /* set remaining child indices */

        n_child_time = t;

        /* initialize convolution summation to zero */

        convolution = (double) 0.0;

        /* loop over each filter coefficient */

        for ( l = 0; l < filter_length; l++ ){

          /* set remaining parent indices */

          n_parent_time = ( t - n_scale * l ) % n_sample;

          while ( n_parent_time < 0 ){
            n_parent_time += n_sample;
          }

          /* perform convolution operation */

          convolution +=
            result->mats[ n_parent_row ].mat.dblmat.data[ n_parent_time ] *
            pd_filter[ l ];

          /* check for interrupts */

          if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
            MUTIL_FREE_WARN( memlist, &list );

            MUTIL_ERROR( "user interrupt" );
            return MUTIL_ERR_INTERRUPT;
          }

        }  /* end loop over filter coefficients */

        /* store results of convolution with filter in proper location */

        result->mats[ n_child_row ].mat.dblmat.data[ n_child_time ] =
          convolution;

      }  /* end loop over time indices     */
    }   /* end loop over local node index */
  }    /* end loop over each scale       */

  /* free nodes corresponding to registered
     modwpt matrix memory, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free malloced space and corresponding linked list structures */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_transform_maximum_overlap_packet()" );

  return MUTIL_ERR_OK;
}

/* The inverse maximum overlap discrete wavelet transform function */
/* Documented in wav_modw.h                                        */
/* Written by William Constantine                                  */

mutil_errcode wavuniv_transform_maximum_overlap_inverse(
   const mat_set  *modwt,
   const mat_set  *filters,
   void           *intrp_ptr,
   univ_mat       *result )
{
   double            sum_scaling;
   double           *pd_old_scaling_coeffs;
   double           *pd_result;
   double           *pd_scaling;
   double           *pd_wavelet;
   double           *pd_wavelet_coeffs;
   mutil_errcode     err;
   sint32            filter_length;
   sint32            j;
   sint32            l;
   sint32            n_level;
   sint32            n_sample;
   sint32            n_scale;
   sint32            nn;
   sint32            t;
   sint32            tt;
   univ_mat          old_scaling_coeffs;
   memlist list;

   MUTIL_INTERRUPT_INIT( intrp_ptr );
   MEMLIST_INIT( list );

   MUTIL_TRACE( "Start wavuniv_transform_maximum_overlap_inverse()" );

   /* obtain sizes of matrices */

   filter_length = MATUNIV_NELEM( &filters->mats[ 0 ] );
   n_level       = modwt->nelem - 1;
   n_sample      = MATUNIV_NELEM( &modwt->mats[ 0 ] );

   /* check inputs */

   err = localfn_imodwt_input_check( modwt, filters );
   if ( err ) return err;

   /* allocate space for output */

   err = matuniv_malloc_register( result, 1, n_sample, MUTIL_DOUBLE, &list );
   MEMLIST_FREE_ON_ERROR( err, &list );

   /* the inversion will proceed as follows:

      (1) the last row of the modwt matrix (which contains
      the scaling coefficients) will be copied to a temporary
      vector (old_scaling_coeffs)

      (2) old_scaling_coeffs will then be used in conjunction with
      the last scales detail coefficients (penultimate row of the
      modwt matrix) to form the next level approximation. the results
      will be written to the result vector as the coefficients are calculated.

      (3) Before the next level of synthesis, the current result vector
      will be copied into the old_scaling_coeffs vector.

      (4) then steps 2-4 will be repeated until full synthesis
      of the original time series.                     */

   /* create temporary storage vector for scaling coefficients */

   err = matuniv_malloc_register( &old_scaling_coeffs, 1,
                              n_sample, MUTIL_DOUBLE, &list );
   MEMLIST_FREE_ON_ERROR( err, &list );

   /* create pointers */

   pd_wavelet            = filters->mats[ 0 ].mat.dblmat.data;
   pd_scaling            = filters->mats[ 1 ].mat.dblmat.data;
   pd_old_scaling_coeffs = old_scaling_coeffs.mat.dblmat.data;
   pd_result             = result->mat.dblmat.data;

   /* copy last level scaling coefficients into temporary vector */

   err = matuniv_assign(
     &modwt->mats[ n_level ],
     intrp_ptr,
     &old_scaling_coeffs );
   MEMLIST_FREE_ON_ERROR( err, &list );

   for( j = n_level; j >= 1; j-- ){

     pd_wavelet_coeffs = modwt->mats[ j - 1 ].mat.dblmat.data;

     n_scale = 1 << (j - 1) ;

     for ( t = 0; t < n_sample; t++ ){

       sum_scaling = 0.0;

       /* perform convolution operations */

       for ( l = 0; l < filter_length ; l++ ){

         nn = ( l * n_scale ) % n_sample;

         tt = t + nn;

         while ( tt > n_sample - 1 ) tt -= n_sample;

         sum_scaling += pd_old_scaling_coeffs[ tt ] * pd_scaling[ l ] +
           pd_wavelet_coeffs[ tt ] * pd_wavelet[ l ];

       } /* end loop over filter index */


       /* store current scaling coefficient in result vector */

       pd_result[ t ] = sum_scaling;

     } /* end loop over time index  */

     /* copy the current scaling coefficients into the
        old_scaling_coeffs vector */

     err = matuniv_assign( result, intrp_ptr, &old_scaling_coeffs );
     MEMLIST_FREE_ON_ERROR( err, &list );

     /* update pointers */

     /* Check for interrupts */

     if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
       MUTIL_ERROR( "user interrupt" );
       MUTIL_FREE_WARN( memlist, &list );

       return MUTIL_ERR_INTERRUPT;
     }

   }  /* end loop over scale index */

   /* free nodes corresponding to registered
      inverse vector memory, but do not free
      the memory itself */

   err = memlist_member_unregister( result, &list );
   MEMLIST_FREE_ON_ERROR( err, &list );

   MUTIL_FREE_WARN( memlist, &list );

   MUTIL_TRACE( "Done with wavuniv_transform_maximum_overlap_inverse()" );

   return MUTIL_ERR_OK;
}

/* The multiresolution decomposition of a wavelet packet basis */
/* Documented in wav_modw.h                                   */
/* Written by William Constantine                             */

mutil_errcode wavuniv_transform_packet_detail(
  const mat_set  *transform,
  const mat_set  *filters,
  sint32          level,
  sint32          node,
  wav_transform   type,
  void           *intrp_ptr,
  univ_mat       *result )
{
  boolean        is_child;
  boolean        is_dwt;
  boolean        is_maximum_overlap;
  boolean        is_multiply_inheritable;
  boolean        is_scaling;
  double         conv_left;
  double         conv_right;
  double         convolution;
  double        *pd_filter;
  double        *pd_old_coeffs;
  double        *pd_result;
  double        *pd_extra;
  memlist        list;
  mutil_errcode  err;
  sint32         Nj;
  sint32         crystal;
  sint32         current_node;
  sint32         filter_length;
  sint32         i;
  sint32         j;
  sint32         k;
  sint32         l;
  sint32         level_max;
  sint32         m;
  sint32         n;
  sint32         n_extra;
  sint32         n_mod;
  sint32         n_sample;
  sint32         n_scale;
  sint32         nn;
  sint32         parent_node;
  sint32         t;
  sint32         tt;
  sint32         u;
  sint32_mat     add_extra;
  sint32_mat     n_coeff;
  univ_mat       old_coeffs;

  MUTIL_INTERRUPT_INIT( intrp_ptr );
  MEMLIST_INIT( list );

  MUTIL_TRACE( "Start wavuniv_transform_packet_detail()" );

  /* check inputs ... */

  /* ... transform and filters inputs */

  err = localfn_imodwt_input_check( transform, filters );
  if ( err ) return err;

  /* ... level argument */

  if ( level < 0 ){
    MUTIL_ERROR( "Specified level cannot be negative" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* ... node input */

  if ( node < 0 ){
    MUTIL_ERROR( "Node must be positive integer or equal to zero." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }


  /* check node argument to ensure it does not exceed
     maximum allowable value, calculate the maximum
     decomposition level, and locate the index of the
     transform matrix set which contains the coefficients
     specified by the local node index */

  switch ( type ){

  /*FALLTHROUGH*/
  case WAV_TRANSFORM_MODWT:
    level_max = transform->nelem - 1;
  case WAV_TRANSFORM_DWT:

    if ( node > 1 ){
      MUTIL_ERROR( "(MO)DWT node must be either 0 or 1." );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

    crystal = level - node;
    break;

  /*FALLTHROUGH*/
  case WAV_TRANSFORM_MODWPT:
  case WAV_TRANSFORM_DWPT:

    level_max = LOCALDEF_ILOG2( transform->nelem + 1 ) - 1;

    /* if ( node >= (sint32) MUTIL_POW( 2, level ) ) { */
    if ( node >= ( 1 << level ) ) {
      MUTIL_ERROR( "(MO)DWPT node must be less than 2^level." );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }

    /* crystal = (sint32) MUTIL_POW( 2, level ) - 1 + node; */
    crystal = ( 1 << level ) - 1 + node;
    break;

  default:
    MUTIL_ERROR( "Transform type is unsupported" );
    return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* in the case that the transform is a DWT, we have to take special
     measures to obtain the number of samples in the original
     time series (a sum of the total number of coefficients in
     the matrix set) and to obtain the number of levels */

  if ( type == WAV_TRANSFORM_DWT ){

    is_dwt = TRUE;

    n_sample = 0;

    for ( j = 0; j < transform->nelem; j++ ){
      n_sample += MATUNIV_NELEM( &(transform->mats[ j ]) );
    }

    level_max = LOCALDEF_ILOG2( n_sample );
  }
  else {

    is_dwt = FALSE;
    n_sample = MATUNIV_NELEM( &transform->mats[ 0 ] );
  }

  /* define a boolean which indicates whether or not the
     transform a maximum overlap transform */

  if ( ( type == WAV_TRANSFORM_MODWT ) | ( type == WAV_TRANSFORM_MODWPT ) ){
    is_maximum_overlap = TRUE;
  }
  else is_maximum_overlap = FALSE;

  /* ... level argument */

  if ( level > level_max ){
    MUTIL_ERROR( "Specified level exceeds maximum level in supplied transform" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* obtain sizes of matrices */

  filter_length = MATUNIV_NELEM( &filters->mats[ 0 ] );

  /* calculate the number of usable (non-extra) coefficients
  for each level of the DWPT or DWT */

  n_extra = 0;

  if ( !is_maximum_overlap ){

    err = mats32_malloc_register( &n_coeff, level + 1, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = mats32_malloc_register( &add_extra, level + 1, 1, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* initialize variables */

    t = n_sample;
    current_node = node;

    is_multiply_inheritable = TRUE;

    for ( j = 0; j <= level; j++ ){

      if ( LOCALDEF_IS_EVEN( t ) ){
        add_extra.data[ j ] = 0;
      }
      else{

        if ( t == 1 ){
          add_extra.data[ j ] = 0;
        }
        else{
          add_extra.data[ j ] = 1;
          n_extra++;
        }
      }

     /* check multiple inheritance:
      can a crystal inherit an extra
      coefficient multiple times
      (from different levels) as
      it is synthesized from levels
      j = J, ..., 1? We impose the
      rule that a crystal can be
      multiply inheritable if the original
      crystal was formed ONLY by
      successive filterings of
      the scaling filter */

      n_mod = current_node % 4;

      if  ( ( n_mod == 1 ) || ( n_mod == 1 ) ){
        is_multiply_inheritable = FALSE;
      }

      current_node = current_node / 2;

      /* override for now and just put in length of series */

      n_coeff.data[ j ] = t;

      /* halve the number of points */

      t = t / 2;
    }
  }

  /* if the transform is a DWT and we do
     have extra coefficients, then set a pointer
     to the last stored extra scaling coefficient.
     we point to the last one because the extra
     scaling coefficients are stored for levels
     0,..,J-1 and we are going to be inverting
     the transform crystal from levels J,..,1,
     i.e. we need to access the stored extra
     scaling coefficients in reverse order. */

  if ( is_dwt & ( n_extra > 0 ) ){

    pd_extra =
      transform->mats[ transform->nelem - 1 ].mat.dblmat.data +
      n_extra - 1;
  }

  /* create temporary storage vector for coefficients */

  err = matuniv_malloc_register( &old_coeffs, 1, n_sample, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* the inversion will proceed as follows:

   (1) locate the appropriate node in the the matrix set
   and copy it to a temporary vector (old_coeffs)

    (2) old_coeffs will then be convolved with the appropriate
    filter (either the scaling or wavelet filter) to form the
    next level approximation. the result will be written to the
    result vector as the coefficients are calculated.

     (3) Before the next level of synthesis, the current result vector
     will be copied into the old_coeffs vector.

      (4) then steps 2-4 will be repeated until a full synthesis
  is accomplished  */

  /* (1) locate the appropriate node in the the matrix set
  and copy it to a temporary vector (old_coeffs) */

  if ( is_maximum_overlap ){
    err = matuniv_assign( &transform->mats[ crystal ], intrp_ptr, &old_coeffs );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }
  else{

    for ( t = 0; t < MATUNIV_NELEM( &transform->mats[ crystal ] ); t++ ){
      old_coeffs.mat.dblmat.data[ t ] = transform->mats[ crystal ].mat.dblmat.data[ t ];
    }
  }

  /* allocate space for output */

  err = matuniv_malloc_register( result, 1, n_sample, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create pointers */

  pd_old_coeffs = old_coeffs.mat.dblmat.data;
  pd_result     = result->mat.dblmat.data;

  /* initialize local node index */

  current_node = node;

  /* (2) old_coeffs will then be convolved with the appropriate
  filter (either the scaling or wavelet filter) to form the
  next level approximation. the result will be written to the
  result vector as the coefficients are calculated.

   (3) Before the next level of synthesis, the current result vector
  will be copied into the old_coeffs vector.   */

  for( j = level; j >= 1; j-- ){

    /* select the appropriate filter */

    n_mod = current_node % 4;

    if ( ( n_mod == 0 ) || ( n_mod == 3 ) ){
      /* scaling filter */
      pd_filter = filters->mats[ 1 ].mat.dblmat.data;
      is_scaling = TRUE;

    }
    else{
      /* wavelet filter */
      pd_filter = filters->mats[ 0 ].mat.dblmat.data;
      is_scaling = FALSE;
    }

    if ( is_maximum_overlap ){

      n_scale = 1 << (j - 1);

      for ( t = 0; t < n_sample; t++ ){

        convolution = 0.0;

        /* perform convolution operations */

        for ( l = 0; l < filter_length; l++ ){

          nn = ( l * n_scale ) % n_sample;

          tt = t + nn;

          while ( tt > n_sample - 1 ) tt -= n_sample;

          convolution += pd_old_coeffs[ tt ] * pd_filter[ l ];

        } /* end loop over filter index */

        /* store current transform coefficient in result vector */

        pd_result[ t ] = convolution;

      } /* end loop over time index  */

    }
    else{

      l = -2;
      m = -1;
      Nj = n_coeff.data[j];

      /* perform convolution operations */

      for ( t = 0; t < Nj; t++ ){

        l += 2;
        m += 2;
        u = t;
        i = 1;
        k = 0;

        conv_left = pd_old_coeffs[ u ] * pd_filter[ i ];
        conv_right = pd_old_coeffs[ u ] * pd_filter[ k ];

        if ( filter_length > 2 ){

          for ( n = 1; n < ( filter_length / 2 ); n++ ){

            u++;
            if ( u >= Nj ) u = 0;
            i += 2;
            k += 2;

            conv_left += pd_old_coeffs[ u ] * pd_filter[ i ];
            conv_right  += pd_old_coeffs[ u ] * pd_filter[ k ];
          }

        }

        /* store current transform coefficients in result vector */

        pd_result[ l ] = conv_left;
        pd_result[ m ] = conv_right;

      } /* end loop over time index  */

        /* if the parent node has an extra coefficient
        and the current detail was a result of a convolution
        with a scaling filter, append the extra parent coefficient
      it to current detail sequence */

      Nj *= 2;

      if ( add_extra.data[ j - 1 ] ){  /* we must append something to current detail */

        if ( is_scaling ){ /* we can append the extra coefficient */

          if ( level == j ){
            is_child = TRUE;
          }
          else{
            is_child = FALSE;
          }

          /*          is_child = (boolean) ( ( level - j ) == 0 ); */

          if ( ( !is_child & is_multiply_inheritable ) | is_child ){

          /* the current level is separated from the original
            by at least two and we are multiply inheritable */

            if ( is_dwt ){

              pd_result[ Nj ] = *pd_extra;
              pd_extra--;
            }
            else{

              /* parent_node = (sint32) MUTIL_POW( 2, j - 1 ) - 1 +
                (sint32) floor( (double) current_node / 2.0 ); */
              /* j is >= 1 */
              parent_node = ( 1 << (j-1) ) - 1 + current_node / 2;

              pd_result[ Nj ] =
                transform->mats[ parent_node ].mat.dblmat.data[ Nj ];
            }

          }
          else pd_result[ Nj ] = (double) 0.0;
        }
        else pd_result[ Nj ] = (double) 0.0;

        Nj++;
      }

    } /* end is_maximum_overlap check  */

    /* copy the current detail coefficients into the old_coeffs vector */

    if ( j != 1 ){
      if ( is_maximum_overlap ){

        err = matuniv_assign( result, intrp_ptr, &old_coeffs );
        MEMLIST_FREE_ON_ERROR( err, &list );
      }
      else{

        for ( t = 0; t < Nj; t++ ){
          old_coeffs.mat.dblmat.data[ t ] = result->mat.dblmat.data[ t ];
        }
      }


      /* change current node to parent node */

      current_node = current_node / 2;

    }

    /* Check for interrupts */

    if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );

      return MUTIL_ERR_INTERRUPT;
    }

   }  /* end loop over level index */

   /* free nodes corresponding to registered
      inverse vector memory, but do not free
      the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_transform_packet_detail()" );

  return MUTIL_ERR_OK;
}


/* Checks the inputs for maximum overlap DWT functions */
/* Written by William Constantine */

static mutil_errcode localfn_modw_input_check(
   const univ_mat *time_series,
   const mat_set  *filters,
   sint32          n_level)
{
   mutil_errcode err;

   MUTIL_TRACE( "Start localfn_modw_input_check()" );

   /*** check time series ... ***/

   /* ... for valid matrix structure and NULL pointer */

   LOCALDEF_CHECK_NULL_POINTER_MODW( time_series, univ_mat, matuniv );

   /* ... to see if it is a vector */

   if ( ( MATUNIV_NCOL( time_series ) != 1 ) &&
     ( MATUNIV_NROW( time_series ) != 1 ) ){
     MUTIL_ERROR( "Time series matrix must be a single column or row." );
     return MUTIL_ERR_ILLEGAL_SIZE;
   }

   /* ... if number of elements is positive */
   if ( MATUNIV_NELEM( time_series ) < 1 ){
      MUTIL_ERROR( "Number of elements in time series must be positive." );
      return MUTIL_ERR_ILLEGAL_SIZE;
   }

   /* ... if the type is not dcomplex */
   if ( time_series->type == MUTIL_DCOMPLEX ){
      MUTIL_ERROR( "Complex time series are currently unsupported." );
      return MUTIL_ERR_ILLEGAL_TYPE;
   }

   /*** check number of levels ... ***/

   if ( n_level <= 0 ){
      MUTIL_ERROR( "Number of decomposition levels must be positive." );
      return MUTIL_ERR_ILLEGAL_VALUE ;
   }

   /*** check filters ***/

   err = localfn_filters_check( filters );
   if ( err ) return err;

   MUTIL_TRACE( "Done with localfn_modw_input_check()" );

   return MUTIL_ERR_OK;
}


/* Checks the inputs for the inverse MODWT function */
/* Written by William Constantine */

static mutil_errcode localfn_imodwt_input_check(
   const mat_set *modwt,
   const mat_set *filters )
{
   mutil_errcode err;

   MUTIL_TRACE( "Start localfn_imodwt_input_check()" );

   /*** check modwt ... ***/

   /* ... for valid matrix structure and NULL pointer */

   LOCALDEF_CHECK_NULL_POINTER_MODW( modwt, mat_set, matset );

   /* ... for type MUTIL_DOUBLE */

   if ( modwt->mats[ 0 ].type != MUTIL_DOUBLE ){
      MUTIL_ERROR( "Transform matrix must be of type MUTIL_DOUBLE." );
      return MUTIL_ERR_ILLEGAL_TYPE;
   }

   /* ... for proper dimension */

   /*** check filters ***/

   err = localfn_filters_check( filters );
   if ( err ) return err;

   MUTIL_TRACE( "Done with localfn_imodwt_input_check()" );

   return MUTIL_ERR_OK;
}


/* Checks the filters for maximum overlap DWT functions */
/* Written by William Constantine */

static mutil_errcode localfn_filters_check(
   const mat_set *filters )
{
   mutil_errcode err;
   sint32        filter_length_scaling;
   sint32        filter_length_wavelet;

   /*** check wavelet filter ... ***/

   /* ... for valid matrix structure and NULL pointer */

   LOCALDEF_CHECK_NULL_POINTER_MODW( filters, mat_set, matset );

   /* ... for type MUTIL_DOUBLE */

   if ( filters->mats[0].type != MUTIL_DOUBLE ){
      MUTIL_ERROR( "Wavelet filter must be of type MUTIL_DOUBLE." );
      return MUTIL_ERR_ILLEGAL_TYPE;
   }

   /* ... to see if it is a vector */

   if ( !MATANY_IS_VEC( &(filters->mats[ 0 ].mat.dblmat) ) ){
      MUTIL_ERROR( "Wavelet filter matrix must be a single column or row." );
      return MUTIL_ERR_ILLEGAL_SIZE;
   }

   /* ... if number of elements is positive */

   filter_length_wavelet = MATUNIV_NELEM( &filters->mats[ 0 ] );

   if ( MATUNIV_NELEM( &filters->mats[ 0 ] ) < 1 ){
      MUTIL_ERROR( "Wavelet filter length must be greater than one." );
      return MUTIL_ERR_ILLEGAL_SIZE;
   }

   /* ... if the length of the filter is even */

   if ( ( filter_length_wavelet % 2 ) != 0 ){
      MUTIL_ERROR( "Wavelet filter length must be even." );
      return MUTIL_ERR_ILLEGAL_VALUE;
   }

   /*** check scaling filter ... ***/

   /* ... for type MUTIL_DOUBLE */

   if ( filters->mats[0].type != MUTIL_DOUBLE ){
      MUTIL_ERROR( "Scaling filter must be of type MUTIL_DOUBLE." );
      return MUTIL_ERR_ILLEGAL_TYPE;
   }

   /* ... to see if it is a vector */

   if ( !MATANY_IS_VEC( &(filters->mats[ 1 ].mat.dblmat) ) ){
      MUTIL_ERROR( "Scaling filter matrix must be a single column or row." );
      return MUTIL_ERR_ILLEGAL_SIZE;
   }

   /* ... if number of elements is positive */

   filter_length_scaling = MATUNIV_NELEM( &filters->mats[ 1 ] );

   if ( MATUNIV_NELEM( &filters->mats[ 1 ] ) < 1 ){
      MUTIL_ERROR( "Scaling filter length must be greater than one." );
      return MUTIL_ERR_ILLEGAL_SIZE;
   }

   /* ... if the length of the filter is even */

   if ( ( filter_length_scaling % 2 ) != 0 ){
      MUTIL_ERROR( "Scaling filter length must be even." );
      return MUTIL_ERR_ILLEGAL_VALUE;
   }

   /*** check consistent filter lengths */

   if ( filter_length_scaling != filter_length_wavelet ){
      MUTIL_ERROR( "Wavelet filter and scaling filter must be the same length." );
      return MUTIL_ERR_ILLEGAL_SIZE;
   }

   return MUTIL_ERR_OK;
}
