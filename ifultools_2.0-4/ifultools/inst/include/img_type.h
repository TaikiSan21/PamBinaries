
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/img_type.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IMG_TYPE_H_
#define IMG_TYPE_H_

/* This file contains typedefs, enums, and structs
   used in image processing functions in the mutils library.
*/

#ifdef __cplusplus
extern "C" {
#endif

/** Enum of connectedness types in images.
 * Image pixels can be considered to be either four-connected or
 * eight-connected to their neighbors for purposes of image processing.
 * With four-connectedness, only North, West, South, and East neighbors
 * of a pixel are considered to be neighbors in the image processing
 * operation.  With eight-connectedness, the NE, NW, SE, and SW neighbors of a
 * pixel are also considered to be neighbors.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @same #typedef enum _img_connected_type img_connected_type;#
 */
enum _img_connected_type
{
  /** Four-connectedness */
  MUTIL_FOUR_CONNECTED,

  /** Eight-connectedness */
  MUTIL_EIGHT_CONNECTED
};


/* See above for documentation on this enum. */
typedef enum _img_connected_type img_connected_type;


/** Enum of distance types in images.
 * There are several types of distance metrics that can be used
 * for image processing, and this enum is used as an input to
 * functions for measuring distances to choose which metric to use.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @see imguniv_transform_distance
 * @same #typedef enum _img_distance_type img_distance_type;#
 */
enum _img_distance_type
{
  /** Montanari distance: a 4-connected neighbor of a
   * pixel is one pixel away, and an 8-connected neighbor is square
   * root of two pixels away -- this is an approximation to
   * Euclidean distance that still allows for a fast algorithm.
   */
  MUTIL_DIST_MONTANARI,

  /** Manhattan (city block) distance: a 4-connected neighbor of a pixel
   * is one pixel away, and an 8-connected neighbor is two pixels away.
   */
  MUTIL_DIST_MANHATTAN,

  /** Chessboard distance: both a 4-connected neighbor and an 8-connected
   * neighbor of a pixel are one pixel away.
   */
  MUTIL_DIST_CHESSBOARD
};


/* See above for documentation on this enum. */
typedef enum _img_distance_type img_distance_type;


/** Enum of thresholding algorithms.
 * Automatic choice of thresholds for image binarization can be
 * done using several criteria.  This enum is used as an input
 * to automatic thresholding functions to choose which algorithm
 * to use.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @same #typedef enum _img_thresh_type img_thresh_type;#
 * @see imguniv_threshold_auto
 */
enum _img_thresh_type
{
  /** Kittler and Illingworth thresholding criterion.
   * This criterion for choosing a threshold value for binarizing
   * an image chooses the threshold value that
   * minimizes the Kullback information
   * distance between the assumed distribution (sum of two
   * Gaussians) and the actual pixel value histogram.  It is calculated
   * following Robert M. Haralick and Linda G.
   * Shapiro, {\bf Computer and Robot Vision, Volume 1}, Addison
   * Wesley, 1992, p. 23.  Haralick and Shapiro recommend using
   * this criterion for machine vision applications.
   */
  MUTIL_THRESH_KITTLER_ILLINGWORTH,

  /** Otsu thresholding criterion.
   * This criterion for choosing a threshold value for binarizing
   * an image chooses the threshold value that
   * minimizes the within-group variance (weighted sum of the two
   * variances) of the two groups
   * of pixel values separated by the threshold value.
   * It is calculated following Robert M. Haralick and Linda G.
   * Shapiro, {\bf Computer and Robot Vision, Volume 1}, Addison
   * Wesley, 1992, p. 20.  Haralick and Shapiro recommend using the
   * Kittler-Illingworth criterion instead.
   */
  MUTIL_THRESH_OTSU,

  /** Kittler and Illingworth thresholding criterion, modification 1.
   * This criterion for choosing a threshold value for binarizing
   * is the same as MUTIL\_THRESH\_KITTLER with the following modification:
   * a threshold value is ruled out if it's Gaussian
   * is less than 20%.  This eliminates the effect of small perturbations
   * to the histogram.
   */
  MUTIL_THRESH_KITTLER_ILLINGWORTH_MOD1

};


/* See above for documentation on this enum. */
typedef enum _img_thresh_type img_thresh_type;


/** Enum of image model types.
 * This enum is used as an input to segmentation functions to choose
 * which image model to use.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @same #typedef enum _img_model_type img_model_type;#
 */
enum _img_model_type
{
  /** Piecewise constant model, basis f = [1] */
  MUTIL_PIECEWISE_CONSTANT,

  /** Piecewise affine model, basis f = [1 x y] */
  MUTIL_PIECEWISE_AFFINE,

  /** Piecewise quadratic model, basis f = [1 x y x*x y*y] */
  MUTIL_PIECEWISE_QUADRATIC,

  /** Piecewise cubic model, basis f = [1 x y x*x y*y x*x*x y*y*y] */
  MUTIL_PIECEWISE_CUBIC
};


/* See above for documentation on this enum. */
typedef enum _img_model_type img_model_type;


/** Enum for image shape features.
 * Constants from this enum can be used as column indices in
 * matrices of features
 * returned from the imguniv\_shape\_features function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @usage For example, to find the area for connected
 *   component numbered label, if feature\_mat is returned from the function,
 *   you could use the MUTIL\_SHAPE\_AREA member of this enum:
 *   #dmatptr = &(feature_mat.mat.dblmat); area = dmatptr->data[MATANY_INDEX(dmatptr, label, MUTIL_SHAPE_AREA)]#
 * @see imguniv_shape_features
 */
enum _img_shape_features {

  /** Column (x) centroid of the object on the image */
  MUTIL_SHAPE_COL_CENT = 0,

  /** Row (y) centroid of the object on the image */
  MUTIL_SHAPE_ROW_CENT,


  /** Column (x) direction standard deviation of the object */
  MUTIL_SHAPE_COL_SIGMA,

  /** Row (y) direction standard deviation of the object */
  MUTIL_SHAPE_ROW_SIGMA,


  /** Standard deviation of the object along its principal axis */
  MUTIL_SHAPE_MX_SIGMA,

  /** Standard deviation of the object perpendicular to its
    principal axis */
  MUTIL_SHAPE_MY_SIGMA,


  /** Region area (connected components) or weighted area (grayscale) (m00) */
  MUTIL_SHAPE_AREA,

  /** Ratio of the area of the object to the area of the entire image */
  MUTIL_SHAPE_FRACTIONAL_AREA,


  /** Orientation of the principle axis from horizontal */
  MUTIL_SHAPE_THETA,

  /** Eccentricity of the object */
  MUTIL_SHAPE_ECCENTRICITY,


  /** Leftmost point of the object on the image */
  MUTIL_SHAPE_COL_MIN,

  /** Rightmost point of the object on the image */
  MUTIL_SHAPE_COL_MAX,

  /** Uppermost point of the object on the image */
  MUTIL_SHAPE_ROW_MIN,

  /** Lowest point of the object on the image */
  MUTIL_SHAPE_ROW_MAX,



  /** Standard moment 00 (moment indices are [row,column] or [y,x]) */
  MUTIL_SHAPE_M00,

  /** Standard moment 01 */
  MUTIL_SHAPE_M01,

  /** Standard moment 10 */
  MUTIL_SHAPE_M10,

  /** Standard moment 11 */
  MUTIL_SHAPE_M11,

  /** Standard moment 02 */
  MUTIL_SHAPE_M02,

  /** Standard moment 20 */
  MUTIL_SHAPE_M20,

  /** Standard moment 12 */
  MUTIL_SHAPE_M12,

  /** Standard moment 21 */
  MUTIL_SHAPE_M21,

  /** Standard moment 03 */
  MUTIL_SHAPE_M03,

  /** Standard moment 30 */
  MUTIL_SHAPE_M30,


  /** Central moment 00 */
  MUTIL_SHAPE_MU00,

  /** Central moment 01 */
  MUTIL_SHAPE_MU01,

  /** Central moment 10 */
  MUTIL_SHAPE_MU10,

  /** Central moment 11 */
  MUTIL_SHAPE_MU11,

  /** Central moment 02 */
  MUTIL_SHAPE_MU02,

  /** Central moment 20 */
  MUTIL_SHAPE_MU20,

  /** Central moment 12 */
  MUTIL_SHAPE_MU12,

  /** Central moment 21 */
  MUTIL_SHAPE_MU21,

  /** Central moment 03 */
  MUTIL_SHAPE_MU03,

  /** Central moment 30 */
  MUTIL_SHAPE_MU30,


  /** Normalized central moment 00 */
  MUTIL_SHAPE_NU00,

  /** Normalized central moment 01 */
  MUTIL_SHAPE_NU01,

  /** Normalized central moment 10 */
  MUTIL_SHAPE_NU10,

  /** Normalized central moment 11 */
  MUTIL_SHAPE_NU11,

  /** Normalized central moment 02 */
  MUTIL_SHAPE_NU02,

  /** Normalized central moment 20 */
  MUTIL_SHAPE_NU20,

  /** Normalized central moment 12 */
  MUTIL_SHAPE_NU12,

  /** Normalized central moment 21 */
  MUTIL_SHAPE_NU21,

  /** Normalized central moment 03 */
  MUTIL_SHAPE_NU03,

  /** Normalized central moment 30 */
  MUTIL_SHAPE_NU30,


  /** Invariant moment 1 */
  MUTIL_SHAPE_M1,

  /** Invariant moment 2 */
  MUTIL_SHAPE_M2,

  /** Invariant moment 3 */
  MUTIL_SHAPE_M3,

  /** Invariant moment 4 */
  MUTIL_SHAPE_M4,

  /** Invariant moment 5 */
  MUTIL_SHAPE_M5,

  /** Invariant moment 6 */
  MUTIL_SHAPE_M6,

  /** Invariant moment 7 */
  MUTIL_SHAPE_M7

};


/** Enum for joint image shape features in two dimensions.
 * Constants from this enum can be used as column indices in
 * matrices of features returned from the
 * the imguniv\_joint\_features\_2d functions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @usage For example, to find the area for connected
 *   component numbered label, if feature\_mat is returned from the function,
 *   you could use the MUTIL\_SHAPE2\_AREA member of this enum:
 *   #dmatptr = &(feature_mat.mat.dblmat); area = dmatptr->data[MATANY_INDEX(dmatptr, label, MUTIL_SHAPE2_AREA)]#
 * @see imgset_joint_features_3d
 * @see _img_joint_shape_features_3d
 */
enum _img_joint_shape_features_2d {

  /** Column (x) centroid of the object on the image */
  MUTIL_SHAPE2_COL_CENT = 0,

  /** Row (y) centroid of the object on the image */
  MUTIL_SHAPE2_ROW_CENT,

  /**************************************************/

  /** Region area */
  MUTIL_SHAPE2_AREA,

  /** Region perimeter */
  MUTIL_SHAPE2_PERIMETER,

  /**************************************************/

  /** Integrated grayscale intensity within a region  */
  MUTIL_SHAPE2_INTENSITY_SUM,

   /** Integrated squared grayscale intensity within a region  */
  MUTIL_SHAPE2_SQUARED_INTENSITY_SUM,

  /** Mean grayscale intensity within a region  */
  MUTIL_SHAPE2_INTENSITY_MEAN,

  /** Variance of grayscale intensity within a region  */
  MUTIL_SHAPE2_INTENSITY_VARIANCE,

  /** Minimum grayscale intensity within a region  */
  MUTIL_SHAPE2_INTENSITY_MIN,

  /** Maximum grayscale intensity within a region  */
  MUTIL_SHAPE2_INTENSITY_MAX,

  /**************************************************/

  /** Leftmost point of the object on the image */
  MUTIL_SHAPE2_COL_MIN,

  /** Rightmost point of the object on the image */
  MUTIL_SHAPE2_COL_MAX,

  /** Uppermost point of the object on the image */
  MUTIL_SHAPE2_ROW_MIN,

  /** Lowest point of the object on the image */
  MUTIL_SHAPE2_ROW_MAX,

  /**************************************************/

  /** Extent of object along column axis */
  MUTIL_SHAPE2_COL_EXTENT,

  /** Extent of object along row axis */
  MUTIL_SHAPE2_ROW_EXTENT,

  /**************************************************/

  /** Extent of object along the longest (first) principal axis */
  MUTIL_SHAPE2_PA1_EXTENT,

  /** Extent of object along the second principal axis */
  MUTIL_SHAPE2_PA2_EXTENT,

  /**************************************************/

  /** Row coordinate of one (min) end point along the first principal axis */
  MUTIL_SHAPE2_PA1_MIN_ROW,

  /** Column coordinate of one (min) end point along the first
      principal axis */
  MUTIL_SHAPE2_PA1_MIN_COL,


  /** Row coordinate of one (max) end point along the first principal axis */
  MUTIL_SHAPE2_PA1_MAX_ROW,

  /** Column coordinate of one (max) end point along the first
      principal axis */
  MUTIL_SHAPE2_PA1_MAX_COL,

  /**************************************************/

  /** Row coordinate of one (min) end point along the second principal axis */
  MUTIL_SHAPE2_PA2_MIN_ROW,

  /** Column coordinate of one (min) end point along the second
      principal axis */
  MUTIL_SHAPE2_PA2_MIN_COL,


  /** Row coordinate of one (max) end point along the second principal axis */
  MUTIL_SHAPE2_PA2_MAX_ROW,

  /** Column coordinate of one (max) end point along the second
      principal axis */
  MUTIL_SHAPE2_PA2_MAX_COL,

  /**************************************************/

  /** Standard 2D spatial moment 00 (moment indices are [row,column]
      or [y,x]) */
  MUTIL_SHAPE2_M00,

  /** Standard 2D spatial moment 01 */
  MUTIL_SHAPE2_M01,

  /** Standard 2D spatial moment 10 */
  MUTIL_SHAPE2_M10,

  /** Standard 2D spatial moment 11 */
  MUTIL_SHAPE2_M11,

  /** Standard 2D spatial moment 02 */
  MUTIL_SHAPE2_M02,

  /** Standard 2D spatial moment 20 */
  MUTIL_SHAPE2_M20,

  /**************************************************/

  /** End delimiter for enum */
  MUTIL_SHAPE2_END
};


/** Enum for joint image shape features in three dimensions.
 * Constants from this enum can be used as column indices in
 * matrices of features returned from the
 * \Ref{imgset_joint_features_3d} function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @usage For example, to find the area for connected
 *   component numbered label, if feature\_mat is returned from the function,
 *   you could use the MUTIL\_SHAPE3\_AREA member of this enum:
 *   #dmatptr = &(feature_mat.mat.dblmat); area = dmatptr->data[MATANY_INDEX(dmatptr, label, MUTIL_SHAPE3_AREA)]#
 * @see imgset_joint_features_3d
 * @see _img_joint_shape_features_2d
 */
enum _img_joint_shape_features_3d {

  /** Column (x) centroid of the object on the image */
  MUTIL_SHAPE3_COL_CENT = 0,

  /** Row (y) centroid of the object on the image */
  MUTIL_SHAPE3_ROW_CENT,

  /** Slice (z) centroid of the object on the image */
  MUTIL_SHAPE3_SLICE_CENT,

  /**************************************************/

  /** Region volume (for 3D) */
  MUTIL_SHAPE3_VOLUME,

  /** Region outer surface area (for 3D) */
  MUTIL_SHAPE3_SURFACE_AREA,

  /**************************************************/

  /** Integrated grayscale intensity within a region  */
  MUTIL_SHAPE3_INTENSITY_SUM,

   /** Integrated squared grayscale intensity within a region  */
  MUTIL_SHAPE3_SQUARED_INTENSITY_SUM,

  /** Mean grayscale intensity within a region  */
  MUTIL_SHAPE3_INTENSITY_MEAN,

  /** Variance of grayscale intensity within a region  */
  MUTIL_SHAPE3_INTENSITY_VARIANCE,

  /** Minimum grayscale intensity within a region  */
  MUTIL_SHAPE3_INTENSITY_MIN,

  /** Maximum grayscale intensity within a region  */
  MUTIL_SHAPE3_INTENSITY_MAX,

  /**************************************************/

  /** Leftmost point of the object on the image */
  MUTIL_SHAPE3_COL_MIN,

  /** Rightmost point of the object on the image */
  MUTIL_SHAPE3_COL_MAX,

  /** Uppermost point of the object on the image */
  MUTIL_SHAPE3_ROW_MIN,

  /** Lowest point of the object on the image */
  MUTIL_SHAPE3_ROW_MAX,

  /** Uppermost slice point of the object on the image (for 3D) */
  MUTIL_SHAPE3_SLICE_MIN,

  /** Lowest slice point of the object on the image (for 3D) */
  MUTIL_SHAPE3_SLICE_MAX,

  /**************************************************/

  /** Extent of object along column axis */
  MUTIL_SHAPE3_COL_EXTENT,

  /** Extent of object along row axis */
  MUTIL_SHAPE3_ROW_EXTENT,

  /** Extent of object along slice axis (for 3D) */
  MUTIL_SHAPE3_SLICE_EXTENT,

  /**************************************************/

  /** Extent of object along the longest (first) principal axis */
  MUTIL_SHAPE3_PA1_EXTENT,

  /** Extent of object along the second principal axis */
  MUTIL_SHAPE3_PA2_EXTENT,

  /** Extent of object along the third principal axis (for 3D) */
  MUTIL_SHAPE3_PA3_EXTENT,

  /**************************************************/

  /** Row coordinate of one (min) end point along the first principal axis */
  MUTIL_SHAPE3_PA1_MIN_ROW,

  /** Column coordinate of one (min) end point along the first
      principal axis */
  MUTIL_SHAPE3_PA1_MIN_COL,

  /** Slice coordinate of one (min) end point along the first principal axis */
  MUTIL_SHAPE3_PA1_MIN_SLICE,


  /** Row coordinate of one (max) end point along the first principal axis */
  MUTIL_SHAPE3_PA1_MAX_ROW,

  /** Column coordinate of one (max) end point along the first
      principal axis */
  MUTIL_SHAPE3_PA1_MAX_COL,

  /** Slice coordinate of one (max) end point along the first principal axis */
  MUTIL_SHAPE3_PA1_MAX_SLICE,

  /**************************************************/

  /** Row coordinate of one (min) end point along the second principal axis */
  MUTIL_SHAPE3_PA2_MIN_ROW,

  /** Column coordinate of one (min) end point along the second
      principal axis */
  MUTIL_SHAPE3_PA2_MIN_COL,

  /** Slice coordinate of one (min) end point along the second
      principal axis */
  MUTIL_SHAPE3_PA2_MIN_SLICE,


  /** Row coordinate of one (max) end point along the second principal axis */
  MUTIL_SHAPE3_PA2_MAX_ROW,

  /** Column coordinate of one (max) end point along the second
      principal axis */
  MUTIL_SHAPE3_PA2_MAX_COL,

  /** Slice coordinate of one (max) end point along the second
      principal axis */
  MUTIL_SHAPE3_PA2_MAX_SLICE,

  /**************************************************/

  /** Row coordinate of one (min) end point along the third principal axis */
  MUTIL_SHAPE3_PA3_MIN_ROW,

  /** Column coordinate of one (min) end point along the third
      principal axis */
  MUTIL_SHAPE3_PA3_MIN_COL,

  /** Slice coordinate of one (min) end point along the third principal axis */
  MUTIL_SHAPE3_PA3_MIN_SLICE,


  /** Row coordinate of one (max) end point along the third principal axis */
  MUTIL_SHAPE3_PA3_MAX_ROW,

  /** Column coordinate of one (max) end point along the third
      principal axis */
  MUTIL_SHAPE3_PA3_MAX_COL,

  /** Slice coordinate of one (max) end point along the third principal axis */
  MUTIL_SHAPE3_PA3_MAX_SLICE,

  /**************************************************/

  /** Standard 3D spatial moment 000 (moment indices are
      [row,column,slice] or [y,x,z]) */
  MUTIL_SHAPE3_M000,

  /** Standard 3D spatial moment 001 */
  MUTIL_SHAPE3_M001,

  /** Standard 3D spatial moment 010 */
  MUTIL_SHAPE3_M010,

  /** Standard 3D spatial moment 100 */
  MUTIL_SHAPE3_M100,

  /** Standard 3D spatial moment 011 */
  MUTIL_SHAPE3_M011,

  /** Standard 3D spatial moment 101 */
  MUTIL_SHAPE3_M101,

  /** Standard 3D spatial moment 110 */
  MUTIL_SHAPE3_M110,

  /** Standard 3D spatial moment 002 */
  MUTIL_SHAPE3_M002,

  /** Standard 3D spatial moment 020 */
  MUTIL_SHAPE3_M020,

  /** Standard 3D spatial moment 200 */
  MUTIL_SHAPE3_M200,

  /**************************************************/

  /** End delimiter for enum */
  MUTIL_SHAPE3_END
};


/** Struct for snakes parameters.
 * The various parameters for the snakes algorithm are stored
 * within this data structure. See the description of the
 * \Ref{imguniv_snake} function for more detailed
 * information about these parameters.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @same
 *   #typedef struct _img_snake_params img_snake_params;#
 * @see imguniv_snake
 */
struct _img_snake_params {
  /** The boundary type for the snake,
   for open curves use MUTIL\_BOUNDARY\_CONTINUE, for closed curves use
   MUTIL\_BOUNDARY\_PERIODIC; other boundary types are not handled */
  mutil_boundary_type  boundary;

  /** The strength of the elastic energy of the snake,
    value must not be less than 0. */
  double elastic_strength;

  /** The strength of the smoothing energy of the snake,
    value must not be less than 0. */
  double smoothing_strength;

  /** The strength of the image energy of the snake,
    value must not be less than 0. */
  double image_strength;

  /** The size of the step taken between each snake iteration,
    value must be greater than 0. */
  double step_size;

  /** The maximum number of snake iterations to be carried out (on output,
    the actual number of iterations carried out) ,
    value must be greater than 0. */
  sint32 num_iterations;

  /** A snake pixel's activity threshold; if it moves more than
    this value between iterations, it is considered active;
    its value must be greater than or equal to 0. */
  double activity_threshold;

  /** A threshold value for the fraction of active snake pixels
    which is used as the convergence criterion for the algorithm;
    value must be greater than or equal to 0. */
  double fraction_active_threshold;
};


/* See documentation for _img_snake_params (above) for description */
typedef struct _img_snake_params img_snake_params;


/** Structure for level set algorithm parameters.
 * The various parameter for the level set algorithms are stored
 * within this data structure.  See the description for
 * \Ref{imgset_levelset_contour_3d} function for more detailed
 * information about these parameters.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @same #typedef struct _img_levelset_params img_levelset_params;#
 * @see imgset_levelset_contour_3d
 *
*/
struct _img_levelset_params {

  /** The direction to propagate the contour; if TRUE
    the contour is propagated inwards, else it is propagated outwards. */
  boolean    inward;

  /** The strength of the length/curvature term of the energy
    functional, must not be less than 0. */
  double     length_penalty_strength;

  /** The strength of the hyperbolic/area term of the energy
    functional, must not be less than 0. */
  double     area_penalty_strength;

  /** Narrow band option; If TRUE, narrow banding is done. */
  boolean    narrow_band;

  /** The width of the narrow band to be used for calculations
    (the maximum bandwidth value would be less than or equal to
    half the bandwidth); it has no effect if narrow\_band is FALSE. */
  double     narrow_bandwidth;

  /** The size of the step taken between each iteration, value must be
    greater than 0. */
  double     step_size;

  /** The maximum number of level set propagation iterations
    (on output the actual number of iterations carried out), must
    be greater than 0. */
  sint32     num_iterations;

  /** The zero set's activity threshold; if a zero set point moves more than
      this value between iterations, it is considered active; its value must
      be greater than or equal to 0. */
  double    activity_threshold;

  /** Threshold value for the fraction of active zero-set pixels, used
    as the convergence criterion for the algorithm; value must in [0,1].
    If zero, then num\_interations will be performed. */
  double    fraction_active_threshold;
};

/* See documentation for _img_levelset_params (above) for description */
typedef struct _img_levelset_params img_levelset_params;


/** Enum of image anisotropic diffusion model types.
 * This enum is used as an input to anisotropic diffusion
 * functions to choose which diffusion model to use.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @same #typedef enum _img_diffusion_model_type img_diffusion_model_type;#
 */
enum _img_diffusion_model_type
{
  /** Model $D = \exp\{ - (\nabla{I}/K)^{2}\}$ */
  MUTIL_PERONA_MALIK_EXPONENTIAL,

  /** Model $D = \frac{1}{1 + (\nabla{I}/K)^{2}}$ */
  MUTIL_PERONA_MALIK_FRACTIONAL
};


/* See above for documentation on this enum. */
typedef enum _img_diffusion_model_type img_diffusion_model_type;


/** Enum of orientations to define spatial relationships between pixels.
 * The orientations are measured in the row-column coordinate axis in
 * clockwise direction relative to the column axis.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @same #typedef enum _img_orientation_type img_orientation_type;#
 */
enum _img_orientation_type
{
  /** 0 degree orientation. */
  MUTIL_DEGREE_0,

  /** 45 degree orientation. */
  MUTIL_DEGREE_45,

  /** 90 degree orientation. */
  MUTIL_DEGREE_90,

  /** 135 degree orientation. */
  MUTIL_DEGREE_135
};


/* This typedef is documented with _img_orientation_type, above. */
typedef enum _img_orientation_type img_orientation_type;


/** Enum of textural features that are computed from co-occurrence matrices.
 * Given a co-occurrence matrix, Haralick et al. defined the angular second
 * moment, contrast, correlation, sum of squares, inverse difference moment,
 * sum average, sum variance, sum entropy, entropy, correlation and maximal
 * correlation coefficient as features that relate to specific textural
 * characteristics of an image such as homogeneity, contrast and the
 * presence of organized structure. Detailed information about co-occurrence
 * matrices is given under
 * \Ref{imguniv_cooccurrence_texture_features_pixel_level}.
 *
 * For more information, see
 * R. M. Haralick, K. Shanmugam, and I. Dinstein, ``Textural features
 * for image classification,'' {\em IEEE Transactions on Systems, Man, and
 * Cybernetics}, vol. SMC-3, no. 6, pp. 610-621, November 1973.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @same #typedef enum _img_cooccurrence_feature_type img_cooccurrence_feature_type;#
 * @see imguniv_cooccurrence_texture_features_pixel_level
 */
enum _img_cooccurrence_feature_type
{
  /** Angular second moment. */
  MUTIL_COOC_ANGULAR_SECOND_MOMENT,

  /** Contrast. */
  MUTIL_COOC_CONTRAST,

  /** Correlation. */
  MUTIL_COOC_CORRELATION,

  /** Sum of squares. */
  MUTIL_COOC_SUM_OF_SQUARES,

  /** Inverse difference moment. */
  MUTIL_COOC_INVERSE_DIFFERENCE_MOMENT,

  /** Sum average. */
  MUTIL_COOC_SUM_AVERAGE,

  /** Sum variance. */
  MUTIL_COOC_SUM_VARIANCE,

  /** Sum entropy. */
  MUTIL_COOC_SUM_ENTROPY,

  /** Entropy. */
  MUTIL_COOC_ENTROPY,

  /** Difference variance. */
  MUTIL_COOC_DIFFERENCE_VARIANCE,

  /** Difference entropy. */
  MUTIL_COOC_DIFFERENCE_ENTROPY,

  /** Information measures of correlation 1. */
  MUTIL_COOC_INFORMATION_CORRELATION1,

  /** Information measures of correlation 2. */
  MUTIL_COOC_INFORMATION_CORRELATION2,

  /** Maximal correlation coefficient, feature not implemented. */
  MUTIL_COOC_MAXIMAL_CORRELATION_COEFFICIENT
};


/* This typedef is documented with _img_cooccurrence_feature_type, above. */
typedef enum _img_cooccurrence_feature_type img_cooccurrence_feature_type;


/** Enum of filter types for Laws' texture energy-based image texture analysis.
 * Laws described a texture-energy approach to measure the amount
 * of variation in an image within fixed-size windows. 5 1D filters
 * of length 5 are defined as
 * \begin{itemize}
 * \item LEVEL = { 1, 4, 6, 4, 1 }
 * \item EDGE = { -1, -2, 0, 2, 1 }
 * \item RIPPLE = { 1, -4, 6, -4, 1 }
 * \item SPOT = { -1, 0, 2, 0, -1 }
 * \item WAVE = { -1, 2, 0, -2, 1 }
 * \end{itemize}
 * and 25 5x5 filter masks (e.g.
 * LEVEL-LEVEL, EDGE-SPOT) are derived using their outer products.
 * These masks are used to filter the input image as the initial step
 * of the texture analysis algorithm.
 *
 * For more information, see
 * K. I. Laws, ``Rapid Texture Classification,'' {\em SPIE Image Processing
 * for Missile Guidance}, vol. 238, pp. 376-380, 1980
 * or
 * George Stockman and Linda G. Shapiro, {\em Computer Vision},
 * Prentice-Hall, 2000.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include img\_type.h
 * @source img\_type.h
 * @library image
 * @same #typedef enum _img_laws_filter_type img_laws_filter_type;#
 * @see imguniv_laws_texture_features
 */
enum _img_laws_filter_type
{
  /** LEVEL. */
  MUTIL_LAWS_LEVEL,

  /** EDGE. */
  MUTIL_LAWS_EDGE,

  /** RIPPLE. */
  MUTIL_LAWS_RIPPLE,

  /** SPOT. */
  MUTIL_LAWS_SPOT,

  /** WAVE. */
  MUTIL_LAWS_WAVE
};


/* This typedef is documented with _img_laws_filter_type, above. */
typedef enum _img_laws_filter_type img_laws_filter_type;


#ifdef __cplusplus
}
#endif


#endif /* IMG_TYPE_H */
