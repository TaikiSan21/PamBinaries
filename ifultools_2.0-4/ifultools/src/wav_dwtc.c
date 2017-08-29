
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/wav_dwtc.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "wav_dwtc.h"

#include "mat_assn.h"
#include "mat_cast.h"
#include "mat_comp.h"
#include "mat_io.h"
#include "mat_sort.h"
#include "mat_set.h"
#include "mat_summ.h"
#include "mat_umat.h"
#include "mat_univ.h"
#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_limit.h"
#include "ut_math.h"
#include "ut_mem.h"
#include <math.h>

/* Static functions declared here and defined at end of file */

static mutil_errcode localfn_dwpt_extra(
  const mat_set  *dwpt,
  wav_dwpt_extra *extra );

static mutil_errcode localfn_dwtc_input_check(
   const univ_mat *time_series,
   const mat_set  *filters,
   sint32          n_level );

static mutil_errcode localfn_dwtc_inverse_input_check(
   const mat_set *dwt,
   const mat_set *filters );

static mutil_errcode localfn_filters_check(
   const mat_set *filters );

static mutil_errcode localfn_dwtc_index(
  sint32            n_level,
  sint32            n_sample,
  void             *intrp_ptr,
  univ_mat         *result );

static mutil_errcode localfn_wavelet_packet_children_synthesis(
  const univ_mat *packet_child_left,
  const univ_mat *packet_child_right,
  const mat_set  *filters,
  const sint32    oscillation_index_left,
  univ_mat       *result );

static mutil_errcode localfn_sort_dwpt_index_vectors(
  sint32_mat *flat,
  sint32_mat *level,
  sint32_mat *osc,
  void       *intrp_ptr,
  sint32_mat *isort );

static mutil_errcode localdef_convert_packet_indices(
  const sint32_mat *transform_indices,
  void             *intrp_ptr,
  sint32_mat       *flat,
  sint32_mat       *level,
  sint32_mat       *osc );

static mutil_errcode localfn_check_packet_input(
  const mat_set *packet );

/* Static macro definitions */

#undef LOCALDEF_ILOG2
#define LOCALDEF_ILOG2( VALUE ) \
(sint32) floor( MUTIL_LOG2( (double) ( VALUE ) + MUTIL_DOUBLE_EPSILON ) )

#undef LOCALDEF_IS_EVEN
#define LOCALDEF_IS_EVEN(N) ( (N % 2) == 0 ? 1 : 0 )

#undef LOCALDEF_IS_ODD
#define LOCALDEF_IS_ODD(N) ( (N % 2) == 1 ? 1 : 0 )

#undef LOCALDEF_CHECK_NULL_MATRIX_POINTER
#define LOCALDEF_CHECK_NULL_MATRIX_POINTER( DATA_PTR, DATA_TYPE, \
                                     TYPE_PREFIX )               \
   err = TYPE_PREFIX ## _validate( DATA_PTR );                   \
   if ( err ) return err;                                        \
   if ( DATA_PTR == ( DATA_TYPE * ) NULL ){                      \
     MUTIL_ERROR( "Pointer to " #DATA_PTR " matrix is NULL" );   \
     return MUTIL_ERR_NULL_POINTER;                              \
   }

#undef LOCALDEF_PARENT_CRYSTAL
#define LOCALDEF_PARENT_CRYSTAL( CHILD_LEVEL, CHILD_OSC, \
 PARENT_LEVEL, PARENT_OSC )                              \
(PARENT_LEVEL) = (CHILD_LEVEL) - 1;                      \
(PARENT_OSC)   = (sint32) (CHILD_OSC) / 2

#undef LOCALDEF_SEQUENCY2FLAT
#define LOCALDEF_SEQUENCY2FLAT( LEVEL, OSCILLATION_INDEX ) \
  ( 1 << ( LEVEL ) ) - 1 + ( OSCILLATION_INDEX )

/*
#undef LOCALDEF_FLAT2SEQUENCY
#define LOCALDEF_FLAT2SEQUENCY( FLAT_INDEX,  LEVEL, OSCILLATION_INDEX )  \
( LEVEL ) = LOCALDEF_ILOG2( ( FLAT_INDEX ) + 1 ); \
( OSCILLATION_INDEX ) = ( FLAT_INDEX ) + 1 - ( 1 << ( LEVEL ) )
*/

#undef LOCALDEF_FLAT2SEQUENCY
#define LOCALDEF_FLAT2SEQUENCY( FLAT_INDEX,  LEVEL, OSCILLATION_INDEX )  \
( LEVEL ) = (sint32) floor( (double) MUTIL_LOG2( ( FLAT_INDEX ) + 1 ) ); \
( OSCILLATION_INDEX ) = ( FLAT_INDEX ) + 1 - ( 1 << ( LEVEL ) )

#undef LOCALDEF_DEALLOCATE_CHILD
#define LOCALDEF_DEALLOCATE_CHILD( CRYSTAL )                        \
if ( memlist_member_exist( crystal[ (CRYSTAL) ], &list ) ){         \
								    \
 err = memlist_member_free( crystal[ (CRYSTAL) ], &list );	    \
 if ( err ){							    \
   MUTIL_FREE_WARN( memlist, &list );				    \
   (void) mutil_free( crystal[ (CRYSTAL) ], sizeof(univ_mat) );     \
   return err;							    \
 }								    \
 								    \
 (void) mutil_free( crystal[ (CRYSTAL) ], sizeof(univ_mat) );       \
}

/* The discrete wavelet transform function  */
/* using convolution style filtering        */
/* Documented in wav_dwtc.h                 */
/* Written by William Constantine           */

mutil_errcode wavuniv_transform_discrete_wavelet_convolution(
  const univ_mat *time_series,
  const mat_set  *filters,
  sint32          n_level,
  void           *intrp_ptr,
  mat_set        *result )
{
  double         Vt;            /* scaling coeff at level j, time t        */
  double         vsum;          /* scaling convolution summation variable  */
  double         wsum;          /* wavelet convolution summation variable  */
  double        *V;             /* pointer to scaling coefficients         */
  double        *W;             /* pointer to wavelet coefficients         */
  double        *g;             /* pointer to scaling filter               */
  double        *h;             /* pointer to wavelet filter               */
  memlist        list;          /* memory list                             */
  mutil_errcode  err;           /* MUTIL error code                        */
  sint32         n_sample;      /* length of original time series          */
  sint32         Nj = 0;        /* number of coefficients per level        */
  sint32         Njm1;          /* # of coefficients at level j-1          */
  sint32         extra = 0;     /* index for  extra scaling coeffs         */
  sint32         filter_length; /* length of the wavelet/scaling filter    */
  sint32         j;             /* index denoting decomposition level      */
  sint32         l;             /* index denoting filter tap               */
  sint32         n_level_max;   /* maximum allowable decomposition levels  */
  sint32         n_return;      /* # of matrices returned in matrix set    */
  sint32         save_extra;    /* logical flag for extra scaling coeffs   */
  sint32         t;             /* index denoting time                     */
  sint32         tt;            /* index denoting time                     */
  sint32         tmplen;
  sint32_mat     ncol;          /* number of cols in matset matrices       */
  sint32_mat     nrow;          /* number of rows in matset matrices       */
  univ_mat       Vextra;        /* vector that hold extra scaling coefs    */
  univ_mat       Vtemp;         /* vector copy of level (j-1) scaling coef */
  univ_mat       result_temp;   /* temp matrix to store transform coeffs   */

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_transform_discrete_wavelet_convolution()" );

  /* avoid lint warning */

  ( void ) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  err = localfn_dwtc_input_check(
    time_series,
    filters,
    n_level );
  if ( err ) return err;

  /* obtain sizes of matrices */

  n_sample = MATUNIV_NELEM( time_series );
  filter_length = MATUNIV_NELEM( &filters->mats[ 0 ] );

  /* calculate maximum possible number of scales */

  n_level_max = LOCALDEF_ILOG2( n_sample );

  /* verify sizes */

  if ( n_level_max < n_level ) {
    MUTIL_ERROR( "Number of decomposition levels exceeds maximum." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* allocate memory for output

     NOTE: we first allocate space for a universal matrix (result_temp)
     to contain all of the transform coefficients. once this
     if filled, we will map these coefficients into a matrix set.   */

  err = matuniv_malloc_register( &result_temp, 1, n_sample,
                             MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* allocate memory for local matrices */

  err = matuniv_malloc_register( &Vtemp, n_sample, 1,
    MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &Vextra, 1, n_level + 1,
    MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* copy input time series into temporary scaling vector */

  if ( MATUNIV_NCOL( time_series) == 1 ) {
    for ( t = 0; t < n_sample; t++ ){
      Vtemp.mat.dblmat.data[ t ] =
        (double) MATUNIV_ELEM( time_series, t, 0 );
    }
  }
  else {
    for ( t = 0; t < n_sample; t++ ){
      Vtemp.mat.dblmat.data[ t ] =
        (double) MATUNIV_ELEM( time_series, 0, t );
    }
  }


  err = matuniv_cast( time_series, FALSE, intrp_ptr, &Vtemp );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* set pointers */

  W = result_temp.mat.dblmat.data;
  h = filters->mats[ 0 ].mat.dblmat.data;
  g = filters->mats[ 1 ].mat.dblmat.data;

  MUTIL_TRACE( "Starting DWT loop over each level." );

  /* calculate the DWT coefficients */

  for ( j = 1, tmplen = 2; j <= n_level; j++, tmplen <<= 1 ){

    /* calculate the number of coefficients in current and previous level */

    /* Nj   = n_sample / ( sint32 ) MUTIL_POW( 2, j );
       Njm1 = n_sample / ( sint32 ) MUTIL_POW( 2, j - 1 ); */
    Nj   = n_sample / tmplen;
    Njm1 = n_sample / ( tmplen / 2 );

    /* set pointer to location to write scaling coefficients */

    V = &( W[ Nj ] );

    /* if the number of scaling coefficients in the previous level is
       odd then store the last scaling coefficient in the extra vector.
       increment counter as well. */

    save_extra = !( LOCALDEF_IS_EVEN( Njm1 ) );

    if ( save_extra ){
      Vextra.mat.dblmat.data[ extra++ ] =
        Vtemp.mat.dblmat.data[ Njm1 - 1 ];
    }

    /* perform convolution operations */

    for ( t = 0; t < Nj; t++ ){

      wsum = vsum = 0.0;

      for ( l = 0; l < filter_length ; l++ ){

        tt = ( 2 * t + 1 - l ) % ( 2 * Nj );

        while ( tt < 0 ) tt += ( 2 * Nj );

        Vt = Vtemp.mat.dblmat.data[ tt ];

        wsum += h[ l ] * Vt;
        vsum += g[ l ] * Vt;
      }

      W[ t ] = wsum;
      V[ t ] = vsum;
    }

    /* move pointer for storing wavelet coefficients */

    W = V;

    /* copy scaling coefficients into temporary vector */

    for ( t = 0; t < Nj; t++ ) Vtemp.mat.dblmat.data[ t ] = V[ t ];

    /* Check for interrupts */

    if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }

  }     /* end loop over decomposition level index */

  /* if extra scaling coefficients were stored add them to end
     of result_temp vector */

  for ( t = 0; t < extra; t++ ){
    result_temp.mat.dblmat.data[ n_sample - extra + t ] =
      Vextra.mat.dblmat.data[ t ];
  }

  MUTIL_TRACE( "Done with DWT loop over each level." );

  /* now map the results into a matrix set ... */

  /* ... calculate the number of matrices to return
     in the matrix set:

     n_level wavelet coefficient vectors +
     1       scaling coefficient vector  +
     (1       extra scaling coefficient vector) */

  n_return = n_level + 1;
  if ( extra > 0 ) n_return++;

  /* ... allocate memory for the dimension vectors
     used to size the matrix set and fill in the
     appropriate dimensions */

  err = mats32_malloc_register( &nrow, 1,  n_return, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &ncol, 1,  n_return, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* ... set dimension of each matrix
     in the matrix set  ... */

  /* ... the wavelet coefficients */

  for ( j = 0, tmplen = 2; j < n_level; j++, tmplen <<= 1 ){

    /* Nj = n_sample / (sint32) MUTIL_POW( 2, j + 1 ); */ ;
    /* j is positive */
    Nj = n_sample / tmplen;

    nrow.data[ j ] = (sint32) 1;
    ncol.data[ j ] = Nj;
  }

  /* ... the scaling coefficients */

  nrow.data[ n_level ] = (sint32) 1;
  ncol.data[ n_level ] = Nj;

  /* ... extra scaling coefficients */

  if ( extra > 0 ){

    nrow.data[ n_level + 1 ] = (sint32) 1;
    ncol.data[ n_level + 1 ] = extra;
  }

  /* allocate memory for matrix set output */

  err = matset_malloc_register( result, 1, &n_return, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices_arbitrary_size(
    result, &nrow, &ncol, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* now map the DWT coefficients to corresponding
     matrices in the matrix set */

  W = result_temp.mat.dblmat.data;

  for ( j = 0; j < result->nelem; j++ ){

    for ( t = 0; t < ncol.data[ j ]; t++ ){

      result->mats[ j ].mat.dblmat.data[ t ] = W[ t ];
    }

    W += ncol.data[ j ];
  }

  /* free node corresponding to registered
     output memory, but do not free the
     memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_transform_discrete_wavelet_convolution()" );

  return MUTIL_ERR_OK;
}


/* The discrete wavelet packet transform function */
/* Documented in wav_dwtc.h                       */
/* Written by William Constantine                 */
mutil_errcode wavuniv_transform_packet(
  const univ_mat *time_series,
  const mat_set  *filters,
  sint32          n_level,
  void           *intrp_ptr,
  mat_set        *result )
{
  double            convolution;
  double           *pd_result;
  double           *pd_filter;
  mutil_errcode     err;
  sint32            j; /* scale index      */
  sint32            l; /* filter index     */
  sint32            n; /* local node index */
  sint32            Nj;/* number of usable coeffs in parent node */
  sint32            Nj2;/* half of the number of usable coeffs in parent node */
  sint32            n_child_base_row;
  sint32            n_child_node;
  sint32            n_child_row;
  sint32            n_mod;
  sint32            n_sample;
  sint32            n_parent_base_row;
  sint32            n_parent_node;
  sint32            n_parent_row;
  sint32            n_scale;
  sint32            filter_length;
  sint32            t;  /* time index       */
  sint32            tt; /* time index       */
  sint32            n_level_max;
  memlist           list;
  sint32            n_osc;
  sint32            n_crystal;
  sint32_mat        nrow;
  sint32_mat        ncol;
  sint32           *ps_nrow;
  sint32           *ps_ncol;


  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MEMLIST_INIT( list );

  MUTIL_TRACE( "Start wavuniv_transform_packet()" );

  /* check inputs */

  err = localfn_dwtc_input_check(
    time_series,
    filters,
    n_level);
  if ( err ) return err;

  /* obtain sizes of matrices */

  filter_length = MATUNIV_NELEM( &filters->mats[ 0 ] );
  n_sample      = MATUNIV_NELEM( time_series );
  n_level_max   = LOCALDEF_ILOG2( n_sample );

  if ( n_level_max < n_level ) {
    MUTIL_ERROR( "Number of decomposition levels exceeds maximum." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( n_sample < 2 ){
    MUTIL_ERROR( "Number of samples in original time series must be at least two." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* allocate memory for the DWPT */

  n_crystal = ( 1 << ( n_level + 1 ) ) - 1;

  err = mats32_malloc_register( &nrow, n_crystal, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &ncol, n_crystal, 1, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  ps_nrow = nrow.data;
  ps_ncol = ncol.data;

  Nj = n_sample;

  for ( j = 0; j <= n_level; j++ ){

    for ( n = 0; n < ( 1 << j ); n++ ){

      *ps_nrow = 1;
      *ps_ncol = Nj;

      ps_nrow++;
      ps_ncol++;
    }

    /* truncation is intended in the following integer division */

    Nj /= 2;
  }

  err = matset_malloc_register( result, 1, &n_crystal, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices_arbitrary_size(
    result, &nrow, &ncol, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* calculate the DWPT ... */

  /* copy time series into DWPT crystal W(0,0).
   make sure to use MATUNIV_ELEM() macro here and cast it to
   a double since (in general) the inptu series can be of any
   class except complex */

  pd_result = result->mats[0].mat.dblmat.data;

  Nj = MATUNIV_NELEM( &( result->mats[0] ) );

  for ( t = 0; t < Nj; t++ ){

    *pd_result = (double) MATUNIV_ELEM( time_series, 0, t );
    pd_result++;
  }

  Nj      = n_sample;
  n_scale = 1;
  n_osc   = 2;

  for ( j = 1; j <= n_level; j++ ){

    if ( LOCALDEF_IS_ODD( Nj ) ) Nj--;

    Nj2 = Nj / 2;

    n_child_base_row  = n_osc - 1;
    n_parent_base_row = n_scale - 1;

    /* loop over each oscillaiton index */

    for ( n = 0; n < n_osc; n++ ){

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

      for ( t = 0; t < Nj2; t++ ) {

        /* initialize convolution summation to zero */

        convolution = (double) 0.0;

        /* loop over each filter coefficient */

        for ( l = 0; l < filter_length; l++ ){

          tt = ( 2 * t + 1 - l ) % Nj;

          while ( tt < 0 ) tt += Nj;

          /* perform convolution operation */

          convolution +=
            result->mats[ n_parent_row ].mat.dblmat.data[ tt ] * pd_filter[ l ];

          /* check for interrupts */

          if ( MUTIL_INTERRUPT( 3.0 * n_sample, intrp_ptr ) ) {
            MUTIL_FREE_WARN( memlist, &list );

            MUTIL_ERROR( "user interrupt" );
            return MUTIL_ERR_INTERRUPT;
          }

        }  /* end loop over filter coefficients */

        /* store results of convolution in proper location */

        result->mats[ n_child_row ].mat.dblmat.data[ t ] = convolution;

      }  /* end loop over time indices     */

    }   /* end loop over oscillation index */

    /* update variables */

    n_scale *= 2;
    n_osc *= 2;

    /* flooring is desired in the following integer division operation */

    Nj /= 2;

  }    /* end loop over each level       */

  /* free nodes corresponding to registered
     DWPT matrix memory, but do not free
     the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free malloced space and corresponding linked list structures */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_transform_packet()" );

  return MUTIL_ERR_OK;
}

/* The inverse discrete wavelet transform function   */
/* using convolution style filtering.                */
/* Documented in wav_dwtc.h                          */
/* Written by William Constantine                    */

mutil_errcode wavuniv_transform_discrete_wavelet_convolution_inverse(
  const mat_set  *dwt,
  const mat_set  *filters,
  void           *intrp_ptr,
  univ_mat       *result )
{
  boolean        extra_counted;
  double         num_ops = 0.0;
  double        *E;
  double        *V;
  double        *Vsynth;
  memlist        list;
  mutil_errcode  err;
  sint32         Nsynth;
  sint32         extra_count = 0;
  sint32         j;
  sint32         t;
  sint32        *n_coefs;
  sint32         n_extra   = 0;
  sint32         n_level;
  sint32         n_level_max;
  sint32         n_sample;
  sint32         n_sample_j;
  sint32         n_scaling = 0;
  sint32         n_wavelet = 0;
  sint32        *extra;
  univ_mat       Vtemp; /* vector copy of level (j-1) scaling coeffs */
  univ_mat       dwt_index;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start "
    "wavuniv_transform_discrete_wavelet_convolution_inverse()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs for errors */

  err = localfn_dwtc_inverse_input_check(
    dwt, filters );
  if ( err ) return err;

  /* calculate original sample size */

  for ( n_sample = 0, j = 0; j < dwt->nelem; j++ ){
    n_sample += MATUNIV_NELEM( &dwt->mats[ j ] );
  }

  /* calculate the number of decomposition levels:
  this is slightly tricky because of the storage system
  for extra scaling coefficients. */

  n_sample_j = n_sample;

  extra_counted = (boolean) FALSE;

  for ( j = 1; j < dwt->nelem; j++ ){

	  if ( j == dwt->nelem - 1 && !extra_counted ) break;

	  if ( LOCALDEF_IS_ODD( n_sample_j ) && !extra_counted ){

		  extra_counted = (boolean) TRUE;
		  break;
	  }

	  n_sample_j /= 2; /* integer division intended */
  }

  n_level = extra_counted ? dwt->nelem - 2 : dwt->nelem - 1;

  /* calculate maximum possible number of levels */

  n_level_max = LOCALDEF_ILOG2( n_sample );

  /* verify sizes */

  if ( n_level_max < n_level ) {
    MUTIL_ERROR( "Number of decomposition levels exceeds maximum." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* allocate memory for output */

  err = matuniv_malloc_register( result, n_sample, 1, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /*  obtain allocation information for wavelet,
      scaling, and extra scaling coefficients in
      returned 1D DWT vector  */

  err = localfn_dwtc_index(
    n_level,
    n_sample,
    intrp_ptr,
    &dwt_index);
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &dwt_index, MEMTYPE_MATUNIV );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* set some pointers to data */

  n_coefs = dwt_index.mat.s32mat.data;
  extra   = &( dwt_index.mat.s32mat.data[ n_level + 1 ]);

  /* calculate the number of wavelet, scaling,
     and extra scaling coefficients */

  for ( j = 0; j <= n_level; j++ ){

    if ( j < n_level ) n_wavelet += n_coefs[ j ];
    else n_scaling = n_coefs[ j ];

    n_extra += extra[ j ];
  }

  /* allocate temporary scaling vector. the length
     of this vector need only be the length of the
     first scale wavelet coefficients vector */

  err = matuniv_malloc_register( &Vtemp, n_coefs[ 0 ], 1, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* set remaining pointers */

  if ( n_extra > 0 ){
    E = dwt->mats[ n_level + 1 ].mat.dblmat.data;
  }

  V      = Vtemp.mat.dblmat.data;
  Vsynth = result->mat.dblmat.data;

  /* copy last level scaling coefficients into
     a temporary scaling vector*/

  for ( t = 0; t < n_scaling; t++ ){
    V[ t ] = dwt->mats[ n_level ].mat.dblmat.data[ t ];
  }

  MUTIL_TRACE( "Starting inverse DWT loop over each level." );

  /* perform the synthesis */

  for( j = n_level; j >= 1; j-- ){

    /*
       call the synthesis routine. we set the
       oscillation index of the left child
       to zero to "fool" the packet synthesis
       into always using g & h for left and right
       filters, respectively.
    */

    err = localfn_wavelet_packet_children_synthesis(
      &Vtemp,
      &( dwt->mats[ j - 1 ] ),
      filters,
      0,
      result );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* copy the scaling coefficients into the temporary vector
       and append any extra scaling coefficients */

    if ( j != 1 ){

      Nsynth = n_coefs[ j - 2 ];

      for ( t = 0; t < Nsynth - 1; t++ ){
        V[ t ] = Vsynth[ t ];
      }

      if ( extra[ j - 1 ] ){
        V[ Nsynth - 1 ] = E[ n_extra - 1 - extra_count++ ];
      }
      else{
        V[ Nsynth - 1 ] = Vsynth[ Nsynth - 1 ];
      }
    }

    /* Check for interrupts */

    num_ops += 10.0 * n_sample;
    if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }
  }  /* end loop over scale index */

  /* add final scaling coefficient if necessary */

  if ( extra[ 0 ] ){
    Vsynth[ n_sample - 1 ] = E[ n_extra - 1 - extra_count ];
  }

  MUTIL_TRACE( "Done with inverse DWT loop over each level." );

  /* free node corresponding to registered
     output memory, but do not free the
     memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with "
    "wavuniv_transform_discrete_wavelet_convolution_inverse()" );

  return MUTIL_ERR_OK;
}


/** Indexing function for the discrete wavelet transform (DWT) using
 * convolution style filtering.  Given J is the number of (convolution
 * style) DWT decomposition levels, this function returns the number
 * of wavelet coefficients for decomposition levels j = 1, ..., J and
 * the number of scaling coefficients at level J. These lengths are
 * stored in the first row of the return matrix.
 *
 * In addition to the (regular) wavelet and scaling coefficients,
 * so-called ``extra'' scaling coefficients can be stored in the
 * convolution style DWT.  Specifically, the last scaling coefficient
 * at a given level is preserved if the number of scaling coefficients
 * (in that level) is odd. For each level j = 0, ..., J-1, this
 * function also returns a logical flag indicating if an extra scaling
 * coefficient was preserved at the corresponding level (indicated by
 * a 0 or 1 in the second row of the return matrix).
 *
 * For more information about the algorithm, see D. B. Percival and
 * A. T. Walden, ``Wavelet Methods for Time Series Analysis,''
 * Cambridge University Press, 2000.
 *
 * @limits Only relevant for the convolution style DWT as defined in
 * the above reference.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwtc.h
 * @source wav\_dwtc.c
 * @library wavelets
 * @usage #err = localfn_dwtc_index( n_level, n_sample, intrp_ptr, &result);#
 * @return Standard mutils error/OK code.
 * @param  n_level    The number of DWT decomposition levels.
 * @param  n_sample   The number of points in the original time series.
 * @param  intrp_ptr  Pointer for implementation of interrupt checking.
 * @param  result     Pointer to a universal matrix which (upon return)
 *                    will contain the index information. The size of the
 *                    result matrix will be 2 x (J+1) where
 *                    J = num\_levels and will be of type MUTIL\_SINT32.
 *                    Given L( ) is the length operator, the first row
 *                    will contain [ L(W\_1) | L(W\_2) | ... | L(W\_J) |
 *                    L(V\_J) ], where W\_j are the wavelet coefficients
 *                    at level j and V\_J are the scaling coefficients at
 *                    level J = num\_levels.  The second row will contain
 *                    only 0's or 1's, with a 1 indicating that an extra
 *                    scaling coefficient was preserved at the corresponding
 *                    level. The result matrix is automatically allocated
 *                    within the function.
 *
 * @see wavuniv_transform_discrete_wavelet_convolution
 * @private
 */
static mutil_errcode localfn_dwtc_index(
  sint32            n_level,
  sint32            n_sample,
  void             *intrp_ptr,
  univ_mat         *result )
{
  double            num_ops = 0.0; /* number of operations            */
  mutil_errcode     err;           /* MUTIL error code                */
  sint32            j;             /* counting variable               */
  sint32            total;         /* variable used to count the
                                      number of recorded coefficients */
  sint32           *extra;         /* pointer to the ``state'' of the
                                      extra scaling coefficients at
				      each level                      */
  sint32           *length;        /* pointer to the number of
				      coefficients at each level      */
  sint32            tmplen;

  memlist           list;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE(
    "Start localfn_dwtc_index()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* allocate memory for output */

  err = matuniv_malloc_register( result, 2, n_level + 1,
    MUTIL_SINT32, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* set pointers */

  length = result->mat.s32mat.data;
  extra  = &( length[ n_level + 1 ] );

  /* check level zero scaling coefficients for odd length */

  total = extra[ 0 ] = LOCALDEF_IS_ODD( n_sample );

  /* loop over each scale and record coefficient vector lengths */

  for ( j = 1, tmplen = 2; j <= n_level; j++, tmplen <<= 1 ){

    /* length[ j - 1 ] = ( sint32 ) floor( n_sample / MUTIL_POW( 2, j ) ); */
    length[ j - 1 ] = n_sample / tmplen;

    /* here we force the number of extra scaling coefficients at
       the last level to be zero because in storing the DWT coefficients
       the last level scaling coefficients are appended to the regular
       wavelet coefficients REGARDLESS of whether or not the number of
       scaling coefficients in that last level is odd or even */

    if ( j < n_level ){
      extra[j] =
        ( LOCALDEF_IS_ODD( length[ j - 1 ] ) & ( length[ j - 1 ] != 1 ) );
    }
    else extra[ j ] = 0;

    total += extra[ j ] + length[ j - 1 ];

    /* check for interrupt */

    num_ops += 3.0 * n_level;
    if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }
  }

  /* finally, calculate the number of regular scaling coefficients
     at the final level */

  length[ n_level ] = n_sample - total;

   /* free node corresponding to registered
     output memory, but do not free the
     memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with "
    "localfn_dwtc_index()" );

  return MUTIL_ERR_OK;
}

/* Checks the inputs for DWT functions */
/* Written by William Constantine */

static mutil_errcode localfn_dwtc_input_check(
   const univ_mat *time_series,
   const mat_set  *filters,
   sint32          n_level )
{
   mutil_errcode err;

   MUTIL_TRACE( "Start localfn_dwtc_input_check()" );

   /*** check time series ... ***/

   /* ... for valid matrix structure and NULL pointer */

   LOCALDEF_CHECK_NULL_MATRIX_POINTER( time_series, univ_mat, matuniv );

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

   MUTIL_TRACE( "Done with localfn_dwtc_input_check()" );

   return MUTIL_ERR_OK;
}

/* Checks the inputs for DWT inverse function */
/* Written by William Constantine */

static mutil_errcode localfn_dwtc_inverse_input_check(
   const mat_set *dwt,
   const mat_set *filters )
{
   mutil_errcode err;
   sint32        j;

   MUTIL_TRACE( "Start localfn_dwtc_inverse_input_check()" );

   /*** check dwt ... ***/

   /* ... for valid matrix structure and NULL pointer */

   LOCALDEF_CHECK_NULL_MATRIX_POINTER( dwt, mat_set, matset );

   /* ... to see if its matrices are vectors */

   for ( j = 0; j < dwt->nelem; j++ ){

     if ( !MATANY_IS_VEC( &( dwt->mats[ j ].mat.dblmat ) ) ){
       MUTIL_ERROR( "All universal matrices in the DWT matrix set "
                    "  must have a single-column or single-row." );
       return MUTIL_ERR_ILLEGAL_SIZE;
     }
   }

   /* ... if the types of the matrices are double */

   if ( dwt->mats[ 0 ].type == MUTIL_DCOMPLEX ){
     MUTIL_ERROR( "All universal matrices in the DWT matrix set "
                  "  must be of type MUTIL_DOUBLE." );
      return MUTIL_ERR_ILLEGAL_TYPE;
   }

   /*** check filters ***/

   err = localfn_filters_check( filters );
   if ( err ) return err;

   MUTIL_TRACE( "Done with localfn_dwtc_inverse_input_check()" );

   return MUTIL_ERR_OK;
}

/* Checks the filters for DWT functions */
/* Written by William Constantine */

static mutil_errcode localfn_filters_check(
   const mat_set *filters )
{
   mutil_errcode err;
   sint32        filter_length_scaling;
   sint32        filter_length_wavelet;

   /*** check wavelet filter ... ***/

   /* ... for valid matrix structure and NULL pointer */

   LOCALDEF_CHECK_NULL_MATRIX_POINTER( filters, mat_set, matset );

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

/* Obtain a subset of DWPT crystals */
/* Written by William Constantine */

mutil_errcode wavuniv_transform_packet_basis(
  const mat_set     *dwpt,
  const univ_mat    *transform_indices,
  void              *intrp_ptr,
  mat_set           *result,
  wav_dwpt_extra    *extra )
{
  double_mat    *pd_dwpt;
  memlist        list;
  mutil_errcode  err;
  sint32         dims;
  sint32         i;
  sint32_mat     ncol;
  sint32_mat     nrow;
  sint32         n_valid;
  univ_mat       flat;
  univ_mat       level;
  univ_mat       osc;

  MUTIL_TRACE( "Start wavuniv_transform_packet_basis()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check inputs ... */

  /*** ... check dwpt  */

  err = localfn_check_packet_input( dwpt );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* convert transform indices and register results with memory manager */

  err = wavuniv_transform_packet_convert_indices(
    transform_indices, intrp_ptr, &flat, &level, &osc );
  MEMLIST_FREE_ON_ERROR( err, &list );

  MUTIL_FREE_WARN( matuniv, &level );
  MUTIL_FREE_WARN( matuniv, &osc );

  err = memlist_member_register( &list, &flat, MEMTYPE_MATUNIV );
  if ( err ){
    MUTIL_FREE_WARN( matuniv, &flat );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* check to make sure that all of the indices do not
     exceed that available in the DWPT */

  err = mats32_number_less_than_scalar( &( flat.mat.s32mat ), dwpt->nelem, intrp_ptr, &n_valid );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* be careful here to compare the number of valid matrices to the
     number of elements in the flattened matrix, because (1) that was
     the matrix we used in the comparison and (2) we will run into problems
     otherwise because transform_indices can be a 2-column or 2-row matrix
     (in the non-flattened case) and the number of elements in this case
     is double that of the number of elements in the corresponding flattened
     vector */

  if ( n_valid != MATUNIV_NELEM( &flat ) ){
    MUTIL_FREE_WARN( memlist, &list );
    MUTIL_ERROR( "Some transform indices exceed that available in the DWPT input" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* allocate memory for the new matrix set */

  dims = MATUNIV_NELEM( &flat );

  err = mats32_malloc_register( &nrow, 1, dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &ncol, 1, dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  for ( i = 0; i < dims; i++ ){

    pd_dwpt = &( dwpt->mats[ flat.mat.s32mat.data[ i ] ].mat.dblmat );

    nrow.data[ i ] = pd_dwpt->nrow;
    ncol.data[ i ] = pd_dwpt->ncol;
  }

  err = matset_malloc_register( result, 1, &dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices_arbitrary_size(
    result, &nrow, &ncol, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* copy subset of of original matrix set */

  for ( i = 0; i < dims; i++ ){

    err = matdbl_assign( &( dwpt->mats[ flat.mat.s32mat.data[ i ] ].mat.dblmat ),
      intrp_ptr,  &( result->mats[ i ].mat.dblmat ) );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* free node corresponding to registered
     output memory, but do not free the
     memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* obtain the extra DWPT atoms */

  err = localfn_dwpt_extra( dwpt, extra );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_transform_packet_basis()" );

  return MUTIL_ERR_OK;
}


/* Inverse wavelet packet transform subset */
/* Written by William Constantine          */

mutil_errcode wavuniv_transform_packet_inverse(
  const mat_set  *dwpt_basis,
  const wav_dwpt_extra *extra,
  const univ_mat *transform_indices,
  const mat_set  *filters,
  void           *intrp_ptr,
  univ_mat       *result )
{
  univ_mat       flat;
  univ_mat       level;
  univ_mat       osc;
  univ_mat     **crystal;
  univ_mat     **crystal_temp;
  memlist        list;
  mutil_errcode  err;
  sint32         n_crystal;
  sint32         n_sets;
  sint32         parent_length;
  sint32         index;
  sint32         ileft;
  sint32         iright;
  sint32         i;
  sint32         j;
  sint32         max_level;
  sint32         min_level;
  sint32_mat     match_index;
  sint32        *ps_osc;
  sint32        *ps_level;
  sint32        *ps_flat;
  univ_mat      *pum_child_left;
  univ_mat      *pum_child_right;
  sint32         parent_level;
  sint32         parent_osc;
  sint32_mat     isort;
  univ_mat      *parent;
  sint32         p;
  boolean        append_extra;
  boolean        any_extra;

  MUTIL_TRACE( "Start wavuniv_transform_packet_inverse()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* convert transform indices and register results with memory manager */

  err = wavuniv_transform_packet_convert_indices(
    transform_indices, intrp_ptr, &flat, &level, &osc );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &flat, MEMTYPE_MATUNIV );
  if ( err ){
    MUTIL_FREE_WARN( matuniv, &flat );
    MUTIL_FREE_WARN( matuniv, &level );
    MUTIL_FREE_WARN( matuniv, &osc );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  err = memlist_member_register( &list, &level, MEMTYPE_MATUNIV );
  if ( err ){
    MUTIL_FREE_WARN( matuniv, &level );
    MUTIL_FREE_WARN( matuniv, &osc );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  err = memlist_member_register( &list, &osc, MEMTYPE_MATUNIV );
  if ( err ){
    MUTIL_FREE_WARN( matuniv, &osc );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /*** check filters ***/

  err = localfn_filters_check( filters );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* initialize variables */

  n_crystal = MATUNIV_NELEM( &flat );

  /* if the only crystal is W(0,0), then
     simply return that crystal */

  if ( n_crystal == 1 && flat.mat.s32mat.data[0] == 0 ){

    MUTIL_FREE_WARN( memlist, &list );

    err = matuniv_malloc( result, dwpt_basis->mats[ 0 ].mat.dblmat.nrow,
      dwpt_basis->mats[ 0 ].mat.dblmat.ncol, MUTIL_DOUBLE );
    if ( err ) return err;

    err = matuniv_assign( &(dwpt_basis->mats[ 0 ]), intrp_ptr, result );
    if ( err ){
      MUTIL_FREE_WARN( matuniv, result );
      return err;
    }

    MUTIL_TRACE( "Done with wavuniv_transform_packet_inverse()" );

    return MUTIL_ERR_OK;
  }

  /* create an array of universal matrix pointers
     which coordinate with the level, osc, and flat
     crystal index vectors */

  err = mutil_malloc_register( (sint32)  n_crystal * sizeof(univ_mat *),
    (void **) &crystal, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  for ( i = 0; i < n_crystal; i++ ){

    crystal[ i ] = &( dwpt_basis->mats[ i ] );
  }

  err = mats32_range( &( level.mat.s32mat ), intrp_ptr, &min_level, &max_level );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* are there any extra DWPT atoms? */

  any_extra = (boolean) ( extra->nelem > 0 );
  append_extra = FALSE;

  /* perform inverse transform */

  for ( j = max_level; j > 0; j-- ){

    /* assign pointers */

    ps_level = level.mat.s32mat.data;
    ps_flat  = flat.mat.s32mat.data;

    /* search the level's index vector to find all crystals at current level */

    err = mats32_compare_scalar(  &( level.mat.s32mat ), MUTIL_RELATION_EQUAL, j,
      intrp_ptr, &match_index, NULL );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_register( &list, &match_index, MEMTYPE_MATS32 );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* check to see if there are any extra coefficients that
       we should tack on after synthesis for the parent level */

    parent_level = j - 1;

    if ( any_extra ){

      p = extra->levelmap.data[ parent_level ];

      append_extra = (boolean) ( p >= 0 );
    }

    /*
      perform checks on crystals:

      (1) make sure that the number of crystals is even. this
      has to be, otherwise we would not be able to produce a
      parent with only a single child.

      (2) make sure that the oscillation indices of adjacent
      crystals are siblings. otherwise, something is wrong
      because a parent crystal can only be created by its children
      crystals
    */

    if ( LOCALDEF_IS_ODD( match_index.nelem ) ){

      MUTIL_ERROR("The number of crystals found at current decomposition level " \
	"is odd. An even number of crystals are required to produce the " \
	"respective parent crystals." );

      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    n_sets = match_index.nelem / 2;

    for ( i = 0; i < n_sets; i++ ){

      /* 'ileft' is the index of the left
	 child of the current pair in the
	 flat, level, osc, and crystal vectors */

      ileft = match_index.data[ 2 * i ];
      iright = ileft + 1;

      ps_osc = osc.mat.s32mat.data;

      parent_osc = ps_osc[ ileft ] / 2 ;

      if ( ps_osc[ ileft ] + 1 != ps_osc[ iright ] ){

	MUTIL_ERROR("An adjacent crystal pair at the current decomposition " \
	  "level was found not to contain siblings. Adjacent crystal pairs " \
	  "must be siblings in order to recreate the parent." );

	MUTIL_FREE_WARN( memlist, &list );
	return MUTIL_ERR_ILLEGAL_VALUE;
      }

      /* assign pointers */

      pum_child_left  = crystal[ ileft ];
      pum_child_right = crystal[ iright ];

      /* ensure that the lengths of the children are the same */

      if ( MATUNIV_NELEM( pum_child_left ) != MATUNIV_NELEM( pum_child_right ) ){

	MUTIL_ERROR("The children crystals are not the same length");
	return MUTIL_ERR_ILLEGAL_SIZE;
      }

      /* now create the parent ... */

      parent_length = MATUNIV_NELEM( pum_child_left ) * 2;

      /* ... adjust length of parent for any extra coefficients */

      if ( append_extra ) parent_length++;

      /* allocate memory for new parent crystal ... */

      /* ... dynamically allocate space for the (parent) universal matrix header */

      err = mutil_malloc( (sint32) sizeof(univ_mat), (void **) &parent );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* ... dynamically allocate space for the (parent) universal matrix data
	 and fill in appropriate information in the header */

      err = matuniv_malloc_register( parent, 1, parent_length, MUTIL_DOUBLE, &list );
      if ( err ){
	MUTIL_FREE_WARN( memlist, &list );
	(void) mutil_free( &parent, sizeof(univ_mat) );
	return err;
      }

      /* create the parent */

      err = localfn_wavelet_packet_children_synthesis(
	pum_child_left,
	pum_child_right,
	filters,
	ps_osc[ ileft ],
	parent );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* tack on an extra coefficient in the case where the parent
	 has odd length (in the original decomposition) */

      if ( append_extra ){

	parent->mat.dblmat.data[ parent_length - 1 ] =
	  extra->atoms.data[ p + parent_osc ];
      }

      /* remove any registered memory of the children
	 via the memory manager. also remove the corresponding
	 dynamically allocated universal matrix header */

      LOCALDEF_DEALLOCATE_CHILD( ileft );
      LOCALDEF_DEALLOCATE_CHILD( iright );

      /* update the crystal, osc, level, and flat vectors
	 with the corresponding parent information */

      index = *match_index.data + i;

      ps_level[ index ] = parent_level;
      ps_osc  [ index ] = parent_osc;
      ps_flat [ index ] = LOCALDEF_SEQUENCY2FLAT( parent_level, parent_osc );
      crystal [ index ] = parent;

    } /* end loop over pairs/sets of children */

    /* shorten the index vectors and crystal vector
       by the number of children pairs which were
       processed */

    n_crystal -= n_sets;

    if ( n_crystal > 0 ){

      err = matuniv_realloc_register( &flat, 1, n_crystal, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = matuniv_realloc_register( &level, 1, n_crystal, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = matuniv_realloc_register( &osc, 1, n_crystal, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = mutil_realloc_register( (void **) &crystal, n_crystal * sizeof(univ_mat *),
	( n_crystal + n_sets ) * sizeof(univ_mat *), &list);
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    /* now rebuild the crystal vector: each element contains
       a pointer to a universal matrix representing a particular
       DWPT crystal. order the new set of crystals (now that the children
       at the current level have been converted to parents) based on
       the flattened index, from low to high. if there is only
       one crystal left, then we are done and do not need to sort. */

    if ( n_crystal > 1 ){

      /* re-sort the crystal index vectors */

      err = localfn_sort_dwpt_index_vectors(
	&( flat.mat.s32mat ),
	&( level.mat.s32mat ),
	&( osc.mat.s32mat ),
	intrp_ptr,
	&isort );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_register( &list, &isort, MEMTYPE_MATS32 );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* re-sort the array containing the pointers to DWPT crystals */

      err = mutil_malloc_register( (sint32)  n_crystal * sizeof(univ_mat *),
	(void **) &crystal_temp, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      /* fill in the temporary vector with pointers to universal
	 matrices corresponding to the sorted crystal index vectors */

      for ( i = 0; i < n_crystal; i++ ){
	crystal_temp[ i ] = crystal[ isort.data[ i ] ];
      }

      for ( i = 0; i < n_crystal; i++ ){

	crystal[ i ] = crystal_temp[ i ];
      }

      /* free memory in preparation for next loop */

      err = memlist_member_free( &isort, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_free( crystal_temp, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );

      err = memlist_member_free( &match_index, &list );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

  } /* end loop over decomposition level */

  /* at this point, 'parent' contains the synthesis so we free
     it from the memory list and assign it to 'result' */

  err = matuniv_malloc_register( result, MATUNIV_NROW( parent ),
    MATUNIV_NCOL( parent ), MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_assign( parent, intrp_ptr, result );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  /* free the 'parent' header */

  err = mutil_free( parent, sizeof(univ_mat) );
  MEMLIST_FREE_ON_ERROR( err, &list );

  MUTIL_TRACE( "Done with wavuniv_transform_packet_inverse()" );

  return MUTIL_ERR_OK;
}

/* Validate a DWPT subset         */
/* Written by William Constantine */

mutil_errcode wavuniv_transform_packet_convert_indices(
  const univ_mat *transform_indices,
  void           *intrp_ptr,
  univ_mat       *flat,
  univ_mat       *level,
  univ_mat       *osc )
{
  mutil_errcode err;

  MUTIL_TRACE( "Start wavuniv_transform_packet_convert_indices()" );

  /* check input for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_MATRIX_POINTER( transform_indices, univ_mat, matuniv );

  /* check input matrix type */

  if ( transform_indices->type != MUTIL_SINT32 ){
    MUTIL_ERROR( "Transform indices matrix must be of type MUTIL_SINT32." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* force the type of the univ_mat matrices to be MUTIL_SINT32 */

  flat->type  = MUTIL_SINT32;
  level->type = MUTIL_SINT32;
  osc->type   = MUTIL_SINT32;

  /* call the conversion routine */

  err = wavs32_transform_packet_convert_indices(
    &( transform_indices->mat.s32mat ), intrp_ptr,
    &( flat->mat.s32mat ), &( level->mat.s32mat ), &( osc->mat.s32mat ) );
  if ( err ) return( err );

  MUTIL_TRACE( "Done with wavuniv_transform_packet_convert_indices()" );

  return MUTIL_ERR_OK;
}


mutil_errcode wavs32_transform_packet_convert_indices(
  const sint32_mat *transform_indices,
  void             *intrp_ptr,
  sint32_mat       *flat,
  sint32_mat       *level,
  sint32_mat       *osc )
{
  double        bandwidth;
  memlist       list;
  mutil_errcode err;
  sint32        i;
  sint32_mat    unique;

  MUTIL_TRACE( "Start wavs32_transform_packet_convert_indices()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* convert transform indices and register results with memory manager */

  err = localdef_convert_packet_indices( transform_indices, intrp_ptr,
    flat, level, osc );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, flat, MEMTYPE_MATS32 );
  if ( err ){
    MUTIL_FREE_WARN( mats32, flat );
    MUTIL_FREE_WARN( mats32, level );
    MUTIL_FREE_WARN( mats32, osc );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  err = memlist_member_register( &list, level, MEMTYPE_MATS32 );
  if ( err ){
    MUTIL_FREE_WARN( mats32, level );
    MUTIL_FREE_WARN( mats32, osc );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  err = memlist_member_register( &list, osc, MEMTYPE_MATS32 );
  if ( err ){
    MUTIL_FREE_WARN( mats32, osc );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* perform tests ... */

  bandwidth = 0.0;

  for ( i = 0; i < osc->nelem; i++ ){

    /* ... check oscillation index range */

    if ( osc->data[ i ] < 0 || osc->data[ i ] > ( ( 1 << level->data[ i ] ) - 1 ) ){

      MUTIL_ERROR( "Oscillation index is out of range." );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    /* accumulate normalized bandwidth */

    bandwidth += 1.0 / (double) ( 1 << ( level->data[ i ] + 1 ) );
  }

  /* ... check for basis, i.e, normalized bandwidth
     spanned by the union of all crystals is 1/2 */

  if ( MUTIL_ABS( bandwidth - 0.5 ) > MUTIL_FLOAT_EPSILON ){
    MUTIL_ERROR( "Normalized bandwidth spanned by union of all crystals is not 1/2." );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* check for redundant crystals. this is easily done by
     flattening the transform indices, finding unique
     values, and checking whether the result has the same
     length as the flattened vector.
  */

  err = mats32_unique( flat, FALSE, intrp_ptr, &unique );
  MEMLIST_FREE_ON_ERROR( err, &list );

  if( unique.nelem != flat->nelem ){
    MUTIL_ERROR( "Crystal indices are not unique." );
    MUTIL_FREE_WARN( mats32, &unique );
    MUTIL_FREE_WARN( memlist, &list );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  MUTIL_FREE_WARN( mats32, &unique );

  /* free node corresponding to registered
     output memory, but do not free the
     memory itself */

  err = memlist_member_unregister( flat, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_unregister( level, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_unregister( osc, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavs32_transform_packet_convert_indices()" );

  return MUTIL_ERR_OK;
}

/* Conversion of wavelet packet indices into flat */
/* and sequency order form                         */
/* Written by William Constantine                 */
static mutil_errcode localdef_convert_packet_indices(
  const sint32_mat *transform_indices,
  void             *intrp_ptr,
  sint32_mat       *flat,
  sint32_mat       *level,
  sint32_mat       *osc )
{
  memlist        list;
  mutil_errcode  err;
  sint32         i;
  sint32         n_crystal;
  sint32         nrow;
  sint32         ncol;
  sint32         osc_start_row;
  sint32         osc_start_col;
  boolean        is_flat;
  sint32_mat     isort;

  MUTIL_TRACE( "Start localdef_convert_packet_indices()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check input for valid matrix structure and NULL pointer */

  LOCALDEF_CHECK_NULL_MATRIX_POINTER( transform_indices, sint32_mat, mats32 );

  /* if tranform_indices is a matrix ( dim = 2 x P ,where P != 1 )
     then the indices are not in flattened form. if is a vector
     whose only two elements are zero, then that means it is crystal
     W(0,0) and the indices are NOT flat. otherwise, it is flat. */

  is_flat = FALSE;

  if ( MATANY_IS_VEC( transform_indices ) ){

    if ( transform_indices->nelem == 2 ){

      if ( !( transform_indices->data[ 0 ] == 0 &&
	transform_indices->data[ 1 ] == 0 ) ){

	is_flat = TRUE;
      }
    }
    else is_flat = TRUE;
  }

  if ( is_flat ){ /* flattened indices */

    n_crystal = transform_indices->nelem;
    nrow      = transform_indices->nrow;
    ncol      = transform_indices->ncol;
  }
  else{ /* indices in {j,n} form */

    /* input matrix is not a vector, so we assume that
       we the indices are stored in either a 2 x P
       or a P x 2 matrix, where P is the number of
       crystals in the transform. check to make sure
       we have either 2 columns or 2 rows and allocate memory
       accordingly if we do */

    if ( transform_indices->nrow != 2 && transform_indices->nrow != 2 ){
      MUTIL_ERROR("The transform indices matrix must either be a vector or" \
	"a matrix with either 2 rows or 2 columns");
      return MUTIL_ERR_ILLEGAL_SIZE;
    }

    if ( transform_indices->nrow == 2 ){ /* 2 rows */

      n_crystal = transform_indices->ncol;

      nrow = 1;
      ncol = n_crystal;

      osc_start_row = 1;
      osc_start_col = 0;
    }
    else{ /* 2 columns */

      n_crystal = transform_indices->nrow;

      nrow = n_crystal;
      ncol = 1;

      osc_start_row = 0;
      osc_start_col = 1;
    }
  }

  /* allocate memory */

  err = mats32_malloc_register( flat, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( level, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( osc, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* convert data */

  if ( is_flat ){

    err = mats32_assign( transform_indices, intrp_ptr, flat );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* convert flattened transform indices to sequency j,n indexing */

    for ( i = 0; i < n_crystal; i++ ){

      LOCALDEF_FLAT2SEQUENCY( flat->data[ i ], level->data[ i ], osc->data[ i ] );
    }
  }
  else{

    err = mats32_extract( transform_indices, 0, 0, intrp_ptr, level );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = mats32_extract( transform_indices, osc_start_row, osc_start_col,
      intrp_ptr, osc );
    MEMLIST_FREE_ON_ERROR( err, &list );

    /* convert sequency {j,n} indexing to flattened transform indices */

    for ( i = 0; i < n_crystal; i++ ){

      flat->data[ i ] = LOCALDEF_SEQUENCY2FLAT( level->data[ i ], osc->data[ i ] );
    }
  }

  /* free node corresponding to registered
     output memory, but do not free the
     memory itself */

  err = memlist_member_unregister( flat, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_unregister( level, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_unregister( osc, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* finally, sort the flat vector of indices
     and re-order the level and oscillation index
     vectors accordingly */

  err = localfn_sort_dwpt_index_vectors( flat, level, osc, intrp_ptr, &isort );
  if ( err ){

    MUTIL_FREE_WARN( mats32, flat );
    MUTIL_FREE_WARN( mats32, level );
    MUTIL_FREE_WARN( mats32, osc );

    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  MUTIL_FREE_WARN( mats32, &isort );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localdef_convert_packet_indices()" );

  return MUTIL_ERR_OK;
}


/** Single stage synthesis of discrete wavelet packet children
 * crystals into a parent crystal.
 * This function reconstructs a discrete wavelet packet parent
 * crystal from its two children crystals.
 *
 * For more information about the algorithm, see D. B. Percival and
 * A. T. Walden, ``Wavelet Methods for Time Series Analysis,''
 * Cambridge University Press, 2000.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwtc.h
 * @source wav\_dwtc.c
 * @library wavelets
 * @usage #err = localfn_wavelet_packet_children_synthesis( &left, &right, &filters, osc, &result );#
 * @return Standard mutils error/OK code.
 * @param  packet_child_left  Pointer to a pre-allocated universal matrix
 *                    of type MUTIL\_DOUBLE containing the left child
 *                    (of the two children crystals)
 *                    to use to reconstruct the parent crystal.
 * @param  packet_child_right  Pointer to a pre-allocated universal matrix
 *                    of type MUTIL\_DOUBLE containing the right child
 *                    (of the two children crystals)
 *                    to use to reconstruct the parent crystal.
 * @param  filters    Pointer to a matrix set containing
 *                    two pre-allocated universal matrices
 *                    of type MUTIL\_DOUBLE and containing
 *                    (in order) the wavelet and scaling
 *                    filter coefficients. The filter matrices
 *                    must be a single-column or single-row
 *                    and each must contain the same number of elements.
 * @param  ocillation_index_left The oscillation index of the left child crystal.
 *                    This is used to control which of the filters
 *                    (scaling or wavelet) are to be used in reverting
 *                    the children crystals.
 * @param  result     Pointer to a pre-allocated universal matrix
 *                    of type MUTIL\_DOUBLE which upon return will
 *                    contain the synthesis of the parent. The number
 *                    of elements in this matrix (typically a vector)
 *                    must be long enought to hold the result (at least
 *                    teice as long as the length of one child crystal).
 * @see wavuniv_transform_packet
 * @private
 */
static mutil_errcode localfn_wavelet_packet_children_synthesis(
  const univ_mat *packet_child_left,
  const univ_mat *packet_child_right,
  const mat_set  *filters,
  const sint32    oscillation_index_left,
  univ_mat       *result )
{
  double          Wleft;
  double          Wright;
  double         *crystal_left;
  double         *crystal_right;
  double         *filter_left;
  double         *filter_right;
  double         *parent;
  mutil_errcode   err;
  sint32          n_mod;
  sint32          filter_length;
  sint32          i, k, l, m, n, u, t;
  sint32          children_length;

  MUTIL_TRACE("Start localfn_wavelet_packet_children_synthesis()" );

  /* check inputs */

  LOCALDEF_CHECK_NULL_MATRIX_POINTER( packet_child_left, univ_mat, matuniv );
  LOCALDEF_CHECK_NULL_MATRIX_POINTER( packet_child_right, univ_mat, matuniv );

  /* define the length of the children crystals.

     Note: due to a specialized and efficient storage system,
     the inverse DWT function
     sends in vectors for the parent and the left child
     that are longer than necessary for the inversion
     of the children into a parent. however, the length
     of the right child is always the length it should be
     (as is the case for the DWPT inversion rotuine). thus,
     we obtain the true length of the children crystals
     via the right child only, i.e., the line below is not
     arbitrary. */

  children_length = MATUNIV_NELEM( packet_child_right );

  if ( children_length * 2 > MATUNIV_NELEM( result ) ){
    MUTIL_ERROR( "Specified parent length in wavelet packet sysnthesis exceeds " \
      "the number of elements in supplied parent crystal" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* set pointers to left and right children */

  crystal_left  = packet_child_left->mat.dblmat.data;
  crystal_right = packet_child_right->mat.dblmat.data;

  /* initialize variables */

  l = -2;
  m = -1;

  n_mod         = oscillation_index_left % 4;
  filter_length = filters->mats[0].mat.dblmat.nelem;

  /* assign pointers */

  parent = result->mat.dblmat.data;

  if ( ( n_mod == 0 ) || ( n_mod == 3 ) ){
    filter_left  = filters->mats[ 1 ].mat.dblmat.data;
    filter_right = filters->mats[ 0 ].mat.dblmat.data;
  }
  else{
    filter_left  = filters->mats[ 0 ].mat.dblmat.data;
    filter_right = filters->mats[ 1 ].mat.dblmat.data;
  }

  /* perform synthesis */

  for ( t = 0; t < children_length; t++ ){

    l += 2;
    m += 2;
    u = t;
    i = 1;
    k = 0;

    Wleft  = crystal_left[ u ];
    Wright = crystal_right[ u ];

    parent[ l ] = filter_right[ i ] * Wright + filter_left[ i ] * Wleft;
    parent[ m ] = filter_right[ k ] * Wright + filter_left[ k ] * Wleft;

    if ( filter_length > 2 ){

      for ( n = 1; n < ( filter_length / 2 ); n++ ){

	u++;
	if ( u >= children_length ) u = 0;
	i += 2;
	k += 2;

	Wleft  = crystal_left[ u ];
	Wright = crystal_right[ u ];

	parent[ l ] += filter_right[ i ] * Wright + filter_left[ i ] * Wleft;
	parent[ m ] += filter_right[ k ] * Wright + filter_left[ k ] * Wleft;
      }
    }
  }

  MUTIL_TRACE("Done with localfn_wavelet_packet_children_synthesis()" );

  return MUTIL_ERR_OK;
}


static mutil_errcode localfn_sort_dwpt_index_vectors(
  sint32_mat *flat,
  sint32_mat *level,
  sint32_mat *osc,
  void       *intrp_ptr,
  sint32_mat *isort )
{
  mutil_errcode  err;
  memlist        list;
  sint32_mat     flat_sorted;
  sint32_mat     level_sorted;
  sint32_mat     osc_sorted;
  sint32         nrow = flat->nrow;
  sint32         ncol = flat->ncol;

  MUTIL_TRACE( "Start localfn_sort_dwpt_index_vectors()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* allocate memory for flattened index vector */

  err = mats32_malloc_register( &flat_sorted, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &level_sorted, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &osc_sorted, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( isort, nrow, ncol, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* sort the vectors */

  err = mats32_sort_index_partial( flat, NULL, intrp_ptr, isort );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_permute( flat, isort, intrp_ptr, &flat_sorted );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_permute( level, isort, intrp_ptr, &level_sorted );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_permute( osc, isort, intrp_ptr, &osc_sorted );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* now map sorted vectors back into originals */

  err = mats32_assign( &flat_sorted, intrp_ptr, flat );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_assign( &level_sorted, intrp_ptr, level );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_assign( &osc_sorted, intrp_ptr, osc );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free node corresponding to registered
     output memory, but do not free the
     memory itself */

  err = memlist_member_unregister( isort, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_sort_dwpt_index_vectors()" );

  return MUTIL_ERR_OK;
}

static mutil_errcode localfn_check_packet_input(
  const mat_set *packet )
{

  mutil_errcode err;

  /* perform wavelet packet matrix set checks ... */

  /* ... check for NULL pointer and validate matrix set */

  LOCALDEF_CHECK_NULL_MATRIX_POINTER( packet, mat_set, matset );

  /* ... check for type MUTIL_DOUBLE */

  if ( packet->mats[ 0 ].type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Transform matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* ... check for odd number of crystals */

  if ( !( LOCALDEF_IS_ODD( packet->nelem ) ) ){
    MUTIL_ERROR( "Number of universal amtrices in a wavelet packet matrix set must be odd" );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  return MUTIL_ERR_OK;
}

/** Extraction of 'extra' DWPT coefficients and a tracking vector.
 * If, for any given parent crystal in a DWPT, the number of
 * coefficients is odd, then then some sort of boundary treatment
 * must be applied to extend the crystal to an even length.
 * As this is somehwat arbitrary, another techniques, which we employ
 * in our DWPT and DWT functions, is to simply preserve the last
 * coefficient of an odd length crystal in a specal 'extra'
 * storage vector and put it back in
 * place upon synthesis. This method has the advantage of preserving
 * the energy of the original crystal with a disadvantage of having
 * to do a little bookkeeping.
 *
 * The extra coefficients and a map to those extra coefficients is
 * facilitated by the wav\_dwpt\_extra structure whose members are
 * an extra 'atoms' vector, a 'levelmap' vector, and an 'nelem'
 * member sued to keep count of the number of extra stored atoms.
 * For a J level DWPT, the levelmap will be
 * J element vector of integers such that
 * the jth element corresponds to the jth level of a DWPT.
 * If the jth element contains a negative integer it denotes
 * that the jth DWPT level does not contain any extra coefficients.
 * Conversely, if a positive integer p appears in the jth element,
 * it means that the jth DWPT level has extra coefficients (if any
 * crystal in a given level does, then they all do in that level),
 * and p denotes the index of the levelmap vector where the extra coefficients
 * for the first crystal (with an oscillation index of zero) is stored.
 * Thus, the extra coefficient for crystal W(j,n) can be found in
 * the (p + n)th element of the atoms vector. The nelem member of the structure
 * contains the total number of extra atoms, which can be zero. Since there
 * are mutils functions whihc complain about zero lenght matrices, this additional
 * member is needed.
 *
 * NOTE 1: This extraction is only needed for synthesis operations and since the
 * DWPT synthesis is based on a subset of the DWPT, this function is called from the
 * the \Ref{wavuniv_transform_packet_basis} function and both the DWPT and
 * extra atom structure is sent to the \Ref{wavuniv_transform_packet_inverse}
 * synthesis function.
 *
 * NOTE 2: By design, there are never any extra coefficients at the last level
 * of DWPT, regardless of whether it is a full decomposition or not. This is
 * why the levelmap vector is J elements long as opposed to the J+1 levels that
 * are stored in a DWPT (since the original time series at elvel 0 is also included).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include wav\_dwtc.h
 * @source wav\_dwtc.c
 * @library wavelets
 * @usage #err = localfn_dwpt_extra( &dwpt, &extra );#
 * @return Standard mutils error/OK code.
 * @param  dwpt       Pointer to a DWPT matrix set containing
 *                    $2^{n_level + 1 } - 1$ single-row universal
 *                    matrices of type MUTIL\_DOUBLE. The order of the
 *                    matrices is given by
 *                    [W\_{0,0} , W\_{1,0}, W\_{1,1}, W\_{2,0}, W\_{2,1},
 *                    W\_{2,2}, W\_{2,3}, ... , W\_{J,0}, ...  W\_{J, N\_J}]
 *                    where J  = n\_level and $N\_J = 2^J$. By definition,
 *                    W\_{0,0} is the original time series.
 * @param  extra      Pointer to a structure of type \Ref{_wav_dwpt_extra}.
 * @see _wav_dwpt_extra
 * @see wavuniv_transform_packet_basis
 * @see wavuniv_transform_packet_inverse
 * @see wavuniv_transform_packet
 * @private
 */
static mutil_errcode localfn_dwpt_extra(
  const mat_set  *dwpt,
  wav_dwpt_extra *extra )
{
  memlist       list;
  mutil_errcode err;
  boolean       save;
  sint32        j;
  sint32        n;
  sint32        nosc;
  sint32        n_crystal;
  sint32        n_level;
  sint32        cumsum;
  sint32        iflat;
  sint32        Nj;
  double       *pd_extra;

  MUTIL_TRACE( "Start localfn_dwpt_extra()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* initialize variables */

  n_crystal = dwpt->nelem;
  n_level   = LOCALDEF_ILOG2( n_crystal + 1 ) - 1;

  iflat        = 0;
  cumsum       = 0;
  extra->nelem = 0;

  /* allocate memory for extra structure */

  err = mats32_malloc_register( &( extra->levelmap ), 1, n_level, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matdbl_malloc_register( &( extra->atoms ), 1, n_crystal - ( 1 << n_level ), &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* assign pointers */

  pd_extra = extra->atoms.data;

  /* calculate map and store extra coefficients */

  for ( j = 0; j < n_level; j++ ){

    Nj = MATUNIV_NELEM( &( dwpt->mats[ iflat ] ) );

    save = LOCALDEF_IS_ODD( Nj );

    extra->levelmap.data[ j ] = save ? cumsum : -1;

    nosc = ( 1 << j );

    if ( save ){

      for ( n = 0; n < nosc; n++ ){

	*pd_extra = dwpt->mats[ iflat + n ].mat.dblmat.data[ Nj - 1 ];
	pd_extra++;
	(extra->nelem)++;
      }

      cumsum += nosc;
    }

    iflat += nosc;
  }

  if ( extra->nelem > 0 ){

    /* reallocate atoms vector to conform to new size */

    err = matdbl_realloc_register( &( extra->atoms ), 1, extra->nelem, &list );
    MEMLIST_FREE_ON_ERROR( err, &list );

    err = memlist_member_unregister( &( extra->atoms ), &list );
    if ( err ){
      MUTIL_FREE_WARN( mats32, &( extra->levelmap ) );
      MEMLIST_FREE_ON_ERROR( err, &list );
    }

    err = memlist_member_unregister( &( extra->levelmap ), &list );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* free memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with localfn_dwpt_extra()" );

  return MUTIL_ERR_OK;
}

