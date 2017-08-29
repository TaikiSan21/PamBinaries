
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/ut_kdtr.c $: $Revision: #1 $, $Date: 2008/03/21 $";
/* This is a self-documenting doc++ file */


/* This file contains definitions for functions used to find neighbors    */
/* of points in N-dimensional Euclidian space. The functions are declared */
/* in fra_neig.h.                                                         */

#include "ut_kdtr.h"

#include "mat_umat.h"
#include "mat_univ.h"

#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_limit.h"
#include "ut_mem.h"

#define FRA_PARTIAL_SORT_CONST   31


/*
 ****************************
 *                          *
 * STATIC (local) VARIABLES *
 *                          *
 ****************************
 */



/* These global variables used to point to memory share by many of the */
/* functions and use fewer arguments in recursive functions.           */

static double   *global_kdtree_data;
static sint32    global_kdtree_npoints;
static sint32    global_dim;
static sint32    global_bucketsize;



/*
 ****************************************
 *                                      *
 * STATIC (local) FUNCTION DECLARATIONS *
 *                                      *
 ****************************************
 */


static void localfn_partial_sort(
  const double   *data,
  const sint32    npoints,
  const sint32    one_or_two,
  sint32         *data_indices,
  sint32         *median_indices );

static double localfn_get_median(
  const double   *data,
  const sint32    npoints,
  sint32         *indices );

static double localfn_get_range(
  const double   *data,
  const sint32    limit,
  const sint32   *indices );

static void localfn_kdtree(
  sint32    npoints,
  sint32   *point_indices,
  sint32   *split_indices,
  double   *medians );



/*
 *********************************
 *                               *
 * LIBRARY FUNCTION DEFINITIONS  *
 *                               *
 *********************************
 */



/* Create a kd-tree structure from an */
/* set of data points.                */
/*                                    */
/* Documented in ut_kdtr.h            */
/* Written by Keith L. Davidson       */
mutil_errcode mutil_kdtree_malloc(
  mutil_kdtree     *kdt,
  const univ_mat   *points,
  const sint32      bucket_size )
{
  double         *data;
  sint32          nrow;
  sint32          ncol;
  sint32          tmp;
  sint32          i;
  sint32          j;
  mutil_errcode   err;

  MUTIL_TRACE( "Start mutil_kdtree_malloc()" );

  /* Avoid lint warnings */
  (void) whatssi;


  /* Error check input arguments */

  /* check the points matrix */
  if ( points == (univ_mat*) NULL ) {
    MUTIL_ERROR( "NULL pointer passed for input points" );
    return MUTIL_ERR_NULL_POINTER;
  }
  err = matuniv_validate( points );
  if ( err ) {
     MUTIL_ERROR( "Input points not a valid universal matrix" );
     return err;
  }
  if ( points->type != MUTIL_DOUBLE ) {
    MUTIL_ERROR( "Input points must be of type MUTIL_DOUBLE" );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  nrow = MATUNIV_NROW( points );
  ncol = MATUNIV_NCOL( points );

  /* check the bucket size */
  if ( bucket_size < 1 || bucket_size > nrow ) {
     MUTIL_ERROR( "Input bucket_size out of range" );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }


  /* Make sure the kdt is not a NULL pointer */
  if ( kdt == (mutil_kdtree*) NULL ) {
     MUTIL_ERROR( "Input kdt is a NULL pointer" );
     return MUTIL_ERR_NULL_POINTER;
  }


  /* Initialize memory list within the kd-tree stucture */
  MEMLIST_INIT( kdt->mlist );


#define LOCALMAC_KDT_EXIT_ON_ERROR \
   if ( err ) MEMLIST_FREE_ON_ERROR( err, &(kdt->mlist) )


  /* Allocate memory for the kd-tree structure fields */
  err = matdbl_malloc_register(
    &(kdt->points),
    nrow,
    ncol,
    &(kdt->mlist) );
  if ( err ) {
    MUTIL_ERROR( "Could not allocate memory for the kd-tree data" );
    return MUTIL_ERR_MEM_ALLOC;
  }
  err = matdbl_malloc_register(
    &(kdt->medians),
    nrow,
    (sint32) 1,
    &(kdt->mlist) );
  LOCALMAC_KDT_EXIT_ON_ERROR;
  err = mats32_malloc_register(
    &(kdt->point_index),
    nrow,
    (sint32) 1,
    &(kdt->mlist) );
  LOCALMAC_KDT_EXIT_ON_ERROR;
  err = mats32_malloc_register(
    &(kdt->split_index),
    nrow,
    (sint32) 1,
    &(kdt->mlist) );
  LOCALMAC_KDT_EXIT_ON_ERROR;


#undef LOCALMAC_KDT_EXIT_ON_ERROR


  /* Copy the points data into the kd-tree structure */
  data = (double*) MATUNIV_DATA( points );
  for ( i = 0; i < nrow; i++ )
    for ( tmp = i * ncol, j = 0; j < ncol; j++ )
      kdt->points.data[ j * nrow + i ] = data[ tmp + j ];


  /* Initialize the kd-tree data members */
  for ( i = 0; i < nrow; i++ ) {
    kdt->point_index.data[ i ] = i;
    kdt->split_index.data[ i ] = (sint32) -1;
    kdt->medians.data[ i ]     = (double) -MUTIL_DOUBLE_MAX;
  }
  kdt->bucket_size = bucket_size;
  kdt->dim         = ncol;
  kdt->npoints     = nrow;


  /* Set some global variables that appear in the recursive functions
     used to create the kd-tree */
  global_kdtree_data    = kdt->points.data;
  global_kdtree_npoints = kdt->npoints;
  global_dim            = kdt->dim;
  global_bucketsize     = kdt->bucket_size;


  /* Create the kd-tree (recursive function) */
  localfn_kdtree(
    kdt->npoints,
    kdt->point_index.data,
    kdt->split_index.data,
    kdt->medians.data );

  MUTIL_TRACE( "Done mutil_kdtree()" );

  return MUTIL_ERR_OK;
}


/* Free memory for a kd-tree structure. */
/*                                      */
/* Documented in ut_kdtr.h              */
/* Written by Keith L. Davidson         */
mutil_errcode mutil_kdtree_free( mutil_kdtree *kdt )
{
  mutil_errcode   err;

  MUTIL_TRACE( "Start mutil_kdtree_free()" );

  /* Check for NULL pointer */
  if ( kdt == (mutil_kdtree*) NULL )
    return MUTIL_ERR_NULL_POINTER;


  /* Check for valid kd-tree structure. NOTE: this also checks */
  /* the validity of the internal memory management list.      */
  err = mutil_kdtree_validate( kdt );
  if ( err ) {
     MUTIL_ERROR( "The kd-tree structure is invalid" );
     return err;
  }


  /* Free the internal memory */
  err = memlist_free( &(kdt->mlist) );
  if ( err ) {
    MUTIL_ERROR( "Could not free memory using internal memory"
      " management list" );
    return err;
  }

  MUTIL_TRACE( "Done with mutil_kdtree_free()" );

  return MUTIL_ERR_OK;
}


/* Validate a kd-tree structure. */
/*                               */
/* Documented in ut_kdtr.h       */
/* Written by Keith L. Davidson  */
mutil_errcode mutil_kdtree_validate( const mutil_kdtree *kdt )
{
  sint32          dim;
  sint32          npoints;
  memlist        *mlist;
  mutil_errcode   err;

  MUTIL_TRACE( "Start with mutil_kdtree_validate()" );

  /* Check for NULL pointer */
  if ( kdt == (mutil_kdtree*) NULL ) {
     MUTIL_ERROR( "Input kdt is a NULL pointer" );
     return MUTIL_ERR_NULL_POINTER;
  }


  /* Begin validating data members */

  /* mlist */
  mlist = (memlist*) &(kdt->mlist);
  if ( mlist->length == (sint32) 0 && mlist->root == (memlist_node*) NULL ) {
    /* an initialized memory list is valid according to memlist_validate(),
       but that's not enough for a kd-tree structure */
    MUTIL_ERROR( "Data member mlist initialized, but empty" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  err = memlist_validate( mlist );
  if ( err ) {
     MUTIL_ERROR( "Data member mlist invalid" );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* points */
  err = matdbl_validate( &(kdt->points) );
  if ( err ) {
     MUTIL_ERROR( "Data member points invalid" );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*  medians */
  err = matdbl_validate( &(kdt->medians) );
  if ( err ) {
     MUTIL_ERROR( "Data member medians invalid" );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*  split indices */
  err = mats32_validate( &(kdt->split_index) );
  if ( err ) {
     MUTIL_ERROR( "Data member split_index invalid" );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*  point indices */
  err = mats32_validate( &(kdt->point_index) );
  if ( err ) {
     MUTIL_ERROR( "Data member point_index invalid" );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* these are for convenience */
  dim     = kdt->points.ncol;
  npoints = kdt->points.nrow;


  /* bucket_size */
  if ( kdt->bucket_size < (sint32) 1 || kdt->bucket_size > npoints ) {
     MUTIL_ERROR( "Data member bucket_size out of range" );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* check for consistent dimensions */
  if ( (kdt->medians.nrow != npoints) ||
       (kdt->split_index.nrow != npoints ) ||
       (kdt->point_index.nrow != npoints ) ||
       (kdt->dim != dim ) ||
       (kdt->npoints != npoints) ) {
     MUTIL_ERROR( "Inconsistent dimensions in kd-tree structure" );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* make sure all internal memory is registered with the internal memory
  management list */
  if ( !memlist_member_exist( (void*) &(kdt->points), mlist ) ||
       !memlist_member_exist( (void*) &(kdt->medians), mlist ) ||
       !memlist_member_exist( (void*) &(kdt->point_index), mlist ) ||
       !memlist_member_exist( (void*) &(kdt->split_index), mlist ) ) {
     MUTIL_ERROR( "Structure memory is not registered with the internal"
        " memory manager" );
     return MUTIL_ERR_ILLEGAL_VALUE;
  }

  MUTIL_TRACE( "Done with mutil_kdtree_validate()" );

   return MUTIL_ERR_OK;
}



/*
 *******************************
 *                             *
 * STATIC FUNCTION DEFINITIONS *
 *                             *
 *******************************
 */



static void localfn_kdtree(
  sint32    npoints,
  sint32   *point_index,
  sint32   *split_index,
  double   *medians )
{
  static double   max_range;
  static double   tmp;
  static sint32   col;
  static sint32   max_range_dim;
  sint32          mid_index;

  MUTIL_TRACE( "Start localfn_kdtree()" );


  /* Using a certain range of points, find the dimension with maximum
     range. Note that the kd-tree data pointer points directly to the
     points, represented by a matrix in column-major order. The number of
     columns in the matrix is the dimension of the space the points are
     drawn from. The number of rows is the number of points in the
     kd-tree. */
  max_range = (double) 0;

  for ( col = 0; col < global_dim; col++ ) {

     /* get the range of column col */
     tmp = localfn_get_range(
        global_kdtree_data + col * global_kdtree_npoints,
        npoints,
        point_index );

     /* update maximal column index if current column range is larger
        than the previous column range */
     if ( tmp > max_range ) {
        max_range     = tmp;
        max_range_dim = col;
     }
  }


  /* Compute the middle index of the current range of points being
     considered */
  mid_index = ( npoints - 1 ) / 2;


  /* The split dimension of the current set of points is given by
     max_range_dim */
  split_index[ mid_index ] = max_range_dim;


  /* Sort the current set of points based on the median */
  medians[ mid_index ] = localfn_get_median(
    global_kdtree_data + max_range_dim * global_kdtree_npoints,
    npoints,
    point_index );



  /* If either of the two subsets of points formed from the split has
     more points than the specified bucket size, then continue
     splitting */
  mid_index++;
  if ( mid_index > global_bucketsize ) {

     /* recursive call */
     localfn_kdtree(
        mid_index,
        point_index,
        split_index,
        medians );
  }
  if ( ( npoints - mid_index ) > global_bucketsize ) {

     /* recursive call */
     localfn_kdtree(
        npoints - mid_index,
        point_index + mid_index,
        split_index + mid_index,
        medians + mid_index );
  }

  MUTIL_TRACE( "Done with localfn_kdtree()" );
}


static double localfn_get_median(
  const double   *data,
  const sint32    npoints,
  sint32         *indices )
{
/*
Finds the median of input vector data with n elements indexed by ind.
In the process of finding the median, it sorts the elements (via
the index) such that the elements in the vector are divided by
the median.  The algorithm is based upon a C version of an
ACM TOMS routine (below).
  */

  static double   ret_val;
  static sint32   kk[ 2 ];
  static sint32   k;

  MUTIL_TRACE( "Start localfn_get_median()" );

  if ( npoints == (sint32) 1 )
    return *data;

  k = ( npoints - 1 ) / 2;

  if ( ( npoints % 2 ) == (sint32) 0 ) {

    kk[ 0 ] = k;
    kk[ 1 ] = k + 1;

    localfn_partial_sort( data, npoints, (sint32) 2, indices, kk );

    ret_val = ( data[ indices[ k ] ] + data[ indices[ k + 1 ] ] ) / 2.0;

  } else {

    localfn_partial_sort( data, npoints, (sint32) 1, indices, &k );

    ret_val = data[ indices[ k ] ];
  }

  MUTIL_TRACE( "Done with localfn_get_median()" );

  return ret_val;
}



static double localfn_get_range(
  const double   *data,
  const sint32    limit,
  const sint32   *indices )
{
   static double   small;
   static double   large;
   static sint32   i;

   MUTIL_TRACE( "Start localfn_get_range()" );

   large = MUTIL_DOUBLE_MAX;
   small = -MUTIL_DOUBLE_MAX;

   for ( i = 0; i < limit; i++ ) {

      if ( data[ indices[ i ] ] < small )
        small = data[ indices[ i ] ];

      if ( data[ indices[ i ] ] > large )
        large = data[ indices[ i ] ];
   }

   MUTIL_TRACE( "Done with localfn_get_range()" );

   return ( large - small );
}




/*     ALGORITHM 410 COLLECTED ALGORITHMS FROM ACM. */
/*     ALGORITHM APPEARED IN COMM. ACM, VOL. 15, NO. 05, */
/*     P. 357. */
static void localfn_partial_sort(
  const double   *data,
  const sint32    npoints,
  const sint32    one_or_two,
  sint32         *data_indices,
  sint32         *median_indices )
{
  register double   tmp_data;
  register sint32   tmp_index;
  register sint32   i;
  register sint32   j;
  register sint32   ij;
  static sint32     indl[ FRA_PARTIAL_SORT_CONST ];
  static sint32     indu[ FRA_PARTIAL_SORT_CONST ];
  static sint32     iu[ FRA_PARTIAL_SORT_CONST ];
  static sint32     il[ FRA_PARTIAL_SORT_CONST ];
  sint32            k;
  sint32            l;
  sint32            m;
  sint32            p;
  sint32            jl;
  sint32            ju;
  const boolean     one = TRUE;

  MUTIL_TRACE( "Done with localfn_partial_sort()" );

/* PARAMETERS TO PSORT HAVE THE FOLLOWING MEANING */
/* A     ARRAY TO BE SORTED */
/* N     NUMBER OF ELEMENTS IN A */
/* IND   ARRAY OF INDICES IN ASCENDING ORDER */
/* NI    NUMBER OF ELEMENTS IN IND */
  jl        = 0;
  ju        = one_or_two - 1;
  indl[ 0 ] = 0;
  indu[ 0 ] = one_or_two - 1;

/* ARRAYS INDL, INDU KEEP ACCOUNT OF THE PORTION OF IND RELATED TO THE */
/* CURRENT SEGMENT OF DATA BEING ORDERED. */
  i = 0;
  j = npoints - 1;
  m = 0;

L5:

  if ( i < j ) {

L10:

     /* FIRST ORDER A(I),A(J),A((I+J)/2), AND USE MEDIAN TO SPLIT THE DATA */
     k  = i;
     ij = ( i + j ) / 2;
     tmp_data  = data[ data_indices[ ij ] ];

     if ( data[ data_indices[ i ] ] > tmp_data ) {

        tmp_index          = data_indices[ ij ];
        data_indices[ ij ] = data_indices[ i ];
        data_indices[ i ]  = tmp_index;
        tmp_data           = data[ data_indices[ ij ] ];
     }
     l = j;

     if ( data[ data_indices[ j ] ] < tmp_data ) {

        tmp_index          = data_indices[ ij ];
        data_indices[ ij ] = data_indices[ j ];
        data_indices[ j ]  = tmp_index;
        tmp_data           = data[ data_indices[ ij ] ];


        if ( data[ data_indices[ i ] ] > tmp_data ) {

           tmp_index          = data_indices[ ij ];
           data_indices[ ij ] = data_indices[ i ];
           data_indices[ i ]  = tmp_index;
           tmp_data           = data[ data_indices[ ij ] ];
        }
     }

     goto L40;

L30:

     tmp_index         = data_indices[ l ];
     data_indices[ l ] = data_indices[ k ];
     data_indices[ k ] = tmp_index;

L40:

     do { --l; } while ( data[ data_indices[ l ] ] > tmp_data );

     tmp_index = data_indices[ l ];

     /* SPLIT THE DATA INTO A(I TO L).LT.T, A(K TO J).GT.T */
     do { ++k; } while ( data[ data_indices[ k ] ] < tmp_data );

     if ( k <= l )
       goto L30;

     indl[ m ] = jl;
     indu[ m ] = ju;
     p         = m;
     ++m;

     /* SPLIT THE LARGER OF THE SEGMENTS */
     if ( ( l - i ) <= ( j - k ) )
       goto L60;

     il[ p ] = i;
     iu[ p ] = l;
     i       = k;

     /* SKIP ALL SEGMENTS NOT CORRESPONDING TO AN ENTRY IN IND */
     do {

        if ( jl > ju )
          goto L70;

        if ( median_indices[ jl ] >= i )
          break;

        ++jl;

     } while ( one );

     indu[ p ] = jl - 1;

     goto L80;

L60:

     il[ p ] = k;
     iu[ p ] = j;
     j       = l;

     do {

        if ( jl > ju )
          goto L70;

        if ( median_indices[ ju ] <= j )
          break;

        --ju;

     } while ( one );

     indl[ p ] = ju + 1;

     goto L80;

  } /* if ( i < j ) */

L70:

  do {

     --m;

     if ( m < 0 )
       goto L100;

     i  = il[ m ];
     j  = iu[ m ];
     jl = indl[ m ];
     ju = indu[ m ];

  } while ( jl > ju );

L80:

  if ( j - i > 10 )
    goto L10;

  if ( i == 1 )
    goto L5;

  --i;

  do {

     ++i;

     if ( i == j )
       goto L70;

     tmp_data = data[ data_indices[ i + 1 ] ];

     if ( data[ data_indices[ i ] ] > tmp_data ) {

        tmp_index = data_indices[ i + 1 ];
        k  = i;

        do {

           data_indices[ k + 1 ] = data_indices[ k ];
           --k;

        } while ( ( k >= 0 ) && ( tmp_data < data[ data_indices[ k ] ] ) );

        data_indices[ k + 1 ] = tmp_index;
     }

  } while ( one );

L100:

   MUTIL_TRACE( "Done with localfn_partial_sort()" );

   return;
}
