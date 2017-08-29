
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/sig_conv.c $: $Revision: #1 $, $Date: 2008/03/21 $  ";
/* This is a self-documenting doc++ file. */

#include "sig_conv.h"

#include "mat_any.h"
#include "mat_cast.h"
#include "mat_set.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "ut_debug.h"
#include "ut_intrp.h"
#include "ut_intrn.h"

/* This file contains function declarations for the convolution and correlation
   functions.
*/



/* Static functions declared here and defined at bottom of file */

static mutil_errcode localfn_sigset_correlate_convolve_check_arguments(
    const mat_set *in_set, const mat_set *kernel_set,
    const sint32_mat *row_step, const sint32_mat *col_step,
    const sint32_mat *row_overlap, const sint32_mat *col_overlap,
    const mat_set *out_set );



/* function documented in sig_conv.h */
/* written by Vikram Chalana */
mutil_errcode siguniv_convolve( const univ_mat *in_sig,
    const univ_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap,
    mutil_boundary_type boundary, void *intrp_ptr, univ_mat *out_sig )
{
  mutil_errcode err;
  univ_mat      cast_in_sig;
  univ_mat      cast_kernel;
  boolean       alloc_cast_in;
  boolean       alloc_cast_ker;

  MUTIL_TRACE( "Start siguniv_convolve()" );

  /* avoid lint warning */
  (void) whatssi;

  if( !in_sig || !kernel || !out_sig ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do standard casting */

  MUTIL_DO_STANDARD_CASTING( in_sig, out_sig->type, cast_in_sig,
    intrp_ptr, alloc_cast_in, err );
  if( err ) {
    return err;
  }

  MUTIL_DO_STANDARD_CASTING( kernel, out_sig->type, cast_kernel,
    intrp_ptr, alloc_cast_ker, err );
  if( err ) {
    MUTIL_FREE_STANDARD_CASTING( cast_in_sig, alloc_cast_in );
    return err;
  }


  /* finally call the appropriate function by
     Switching on the type of the output matrix */
  switch( out_sig->type ) {
     case MUTIL_DCOMPLEX:
       err = sigcpx_convolve( &(cast_in_sig.mat.cpxmat),
           &(cast_kernel.mat.cpxmat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.cpxmat) );
       break;

     case MUTIL_DOUBLE:
       err = sigdbl_convolve( &(cast_in_sig.mat.dblmat),
           &(cast_kernel.mat.dblmat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       err = sigflt_convolve( &(cast_in_sig.mat.fltmat),
           &(cast_kernel.mat.fltmat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       err = sigu8_convolve( &(cast_in_sig.mat.u8mat),
           &(cast_kernel.mat.u8mat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       err = sigu16_convolve( &(cast_in_sig.mat.u16mat),
           &(cast_kernel.mat.u16mat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       err = sigu32_convolve( &(cast_in_sig.mat.u32mat),
           &(cast_kernel.mat.u32mat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.u32mat) );
       break;

     case MUTIL_SINT32:
       err = sigs32_convolve( &(cast_in_sig.mat.s32mat),
           &(cast_kernel.mat.s32mat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.s32mat) );
       break;

    default:
      MUTIL_ERROR( "This matrix type is currently unsupported " );
      err = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrices */
  MUTIL_FREE_STANDARD_CASTING( cast_in_sig, alloc_cast_in );
  MUTIL_FREE_STANDARD_CASTING( cast_kernel, alloc_cast_ker );

  if( err ) {
    return err;
  }

  MUTIL_TRACE( "Done with siguniv_convolve()" );
  return MUTIL_ERR_OK;
}


/** Template macro for image convolution or correlation.
 * Macro that expands to the body of a non-universal signal
 * convolve or correlate function, such as sigu8\_convolve or
 * sigu8\_correlate. Complex data types are not supported by
 * this macro; they are handled by separate sigcpx\_convolve and
 * sigcpx\_correlate functions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source sig\_conv.c
 * @library signal
 * @param SIG_FN_PREFIX   Prefix for functions for this matrix type.
 * @param SIG_IN_PTR      Pointer to the input matrix (function argument).
 * @param KERNEL_PTR      Pointer to the kernel matrix (function argument).
 * @param ROW_STEP        Integer for row step (function argument)
 * @param COL_STEP        Integer for column step (function argument)
 * @param ROW_OVERLAP     Number of kernel rows initially overlapping the
 *                        input matrix (function argument).
 * @param COL_OVERLAP     Number of kernel columns initially overlapping the
 *                        input matrix (function argument).
 * @param BOUNDARY_TYPE   Boundary condition to use (function argument)
 * @param INTRP_PTR       Pointer for interrupt handling (function argument).
 * @param SIG_OUT_PTR     Pointer to the output matrix (function argument)
 * @param DATA_TYPE       Data type of matrix data.
 * @param CONVOLUTION     TRUE to do convolution, FALSE for correlation.
 * @usage Body of the imgu8\_convolve function:
 * #TMPL_SIG_CONVOLVE(u8, in_sig, kernel, row_step, col_step, row_overlap, col_overlap, boundary, intrp_ptr, out_sig, uint8);#
 * @private
 */
#define TMPL_SIG_CONVOLVE( SIG_FN_PREFIX, SIG_IN_PTR, KERNEL_PTR, \
    ROW_STEP, COL_STEP, ROW_OVERLAP, COL_OVERLAP, BOUNDARY_TYPE, INTRP_PTR, \
    SIG_OUT_PTR, DATA_TYPE, CONVOLUTION ) \
  sint32    in_nrow;        /* number of input rows */ \
  sint32    in_ncol;        /* number of input columns */ \
  sint32    out_nrow;       /* number of output rows */ \
  sint32    out_ncol;       /* number of output columns */ \
  sint32    kern_nrow;      /* number of kernel rows */ \
  sint32    kern_ncol;      /* number of kernel rows */ \
  sint32    kern_row_offset; \
  sint32    kern_col_offset; \
  sint32    row_output; \
  sint32    col_output; \
  sint32    row_kernel; \
  sint32    col_kernel; \
  sint32    row_input; \
  sint32    col_input; \
  sint32    outsig_idx; \
  \
  mutil_errcode err; \
  double        num_ops = 0; \
  \
  DATA_TYPE    *pinsig; \
  DATA_TYPE    inval; \
  DATA_TYPE    *poutsig; \
  DATA_TYPE    *pkern; \
  \
  MUTIL_INTERRUPT_INIT( intrp_ptr ); \
  \
  if( CONVOLUTION ) { \
    MUTIL_TRACE( "Start sig" #SIG_FN_PREFIX "_convolve()" ); \
  } \
  else { \
    MUTIL_TRACE( "Start sig" #SIG_FN_PREFIX "_correlate()" ); \
  } \
  \
  /* sanity checks */\
  err = mat ## SIG_FN_PREFIX ## _validate( SIG_IN_PTR );\
  if( err ) {\
    return err;\
  }\
  \
  err = mat ## SIG_FN_PREFIX ## _validate( KERNEL_PTR );\
  if( err ) {\
    return err;\
  }\
  \
  err = mat ## SIG_FN_PREFIX ## _validate( SIG_OUT_PTR );\
  if( err ) {\
     return err;\
  }\
  \
  if( (SIG_OUT_PTR)->nrow > (SIG_IN_PTR)->nrow + (KERNEL_PTR)->nrow - 1 || \
      (SIG_OUT_PTR)->ncol > (SIG_IN_PTR)->ncol + (KERNEL_PTR)->ncol - 1 ) { \
      MUTIL_ERROR( "Output dimensions too big" ); \
      return MUTIL_ERR_ILLEGAL_SIZE; \
      }  \
      \
  if( (SIG_IN_PTR)->data == (SIG_OUT_PTR)->data ) { \
    MUTIL_ERROR( "Input matrix cannot be the same as the output matrix" ); \
    return MUTIL_ERR_ILLEGAL_ADDRESS; \
  } \
  \
  if ( (KERNEL_PTR)->data == (SIG_OUT_PTR)->data ) { \
    MUTIL_ERROR( "Kernel matrix cannot be the same as the output matrix" ); \
    return MUTIL_ERR_ILLEGAL_ADDRESS; \
  } \
  \
  if( (ROW_STEP) < 1 || (COL_STEP) < 1 ) { \
    MUTIL_ERROR( "The step size must be at least 1" ); \
    return MUTIL_ERR_ILLEGAL_VALUE; \
  } \
  \
  if( (ROW_OVERLAP) < 1 || (COL_OVERLAP) < 1 ) { \
    MUTIL_ERROR( "The initial overlap must be at least 1" ); \
    return MUTIL_ERR_ILLEGAL_VALUE; \
  } \
  if( (ROW_OVERLAP) > (KERNEL_PTR)->nrow || \
      (COL_OVERLAP) > (KERNEL_PTR)->ncol ) { \
    MUTIL_ERROR( "Initial overlap must be less then the kernel dimensions" ); \
    return MUTIL_ERR_ILLEGAL_VALUE; \
  } \
  \
  /* get the sizes */ \
  in_nrow         = (SIG_IN_PTR)->nrow; \
  in_ncol         = (SIG_IN_PTR)->ncol; \
  out_nrow        = (SIG_OUT_PTR)->nrow; \
  out_ncol        = (SIG_OUT_PTR)->ncol; \
  kern_nrow       = (KERNEL_PTR)->nrow; \
  kern_ncol       = (KERNEL_PTR)->ncol; \
  kern_row_offset = kern_nrow - row_overlap; \
  kern_col_offset = kern_ncol - col_overlap; \
  \
  /* get the underlying flat arrays */ \
  pinsig  = (SIG_IN_PTR)->data; \
  poutsig = (SIG_OUT_PTR)->data; \
  pkern   = (KERNEL_PTR)->data; \
  \
  /* Loop over the output matrix coordinates */ \
  for( row_output = 0; row_output < out_nrow; row_output ++ ) { \
    outsig_idx = out_ncol * row_output; \
    for( col_output = 0; col_output < out_ncol; col_output ++ ) { \
      \
      /* Set initial value to zero */ \
      poutsig[ outsig_idx ] = 0; \
      \
      /* Loop over the kernel coordinates */ \
      for( row_kernel = 0;  row_kernel < kern_nrow; row_kernel ++ ) { \
       row_input = row_output * (ROW_STEP) - kern_row_offset + row_kernel; \
       for( col_kernel = 0;  col_kernel < kern_ncol; col_kernel ++ ) { \
         \
         /* Compute the column for this kernel elem in input matrix */ \
         \
         col_input = col_output * (COL_STEP) - kern_col_offset + col_kernel; \
         \
         switch( BOUNDARY_TYPE ) { \
           case MUTIL_BOUNDARY_ZERO: \
             if(( row_input < 0 ) || \
                ( col_input < 0 ) || \
                ( row_input > in_nrow - 1 ) || \
                ( col_input > in_ncol - 1 )) { \
                 inval = 0; \
             } \
             else { \
               inval = pinsig[ MATANY_INDEX( SIG_IN_PTR, row_input, \
                 col_input )]; \
             } \
             break; \
           \
           case MUTIL_BOUNDARY_CONTINUE: \
             if( row_input < 0 ) { \
               row_input = 0; \
             } \
             else if( row_input > in_nrow - 1 ) { \
               row_input = in_nrow - 1; \
             } \
             if( col_input < 0 ) { \
               col_input = 0; \
             } \
             else if( col_input > in_ncol - 1 ) { \
               col_input = in_ncol - 1; \
             } \
             inval = pinsig[ MATANY_INDEX( SIG_IN_PTR, row_input, \
               col_input )]; \
             break; \
           \
           case MUTIL_BOUNDARY_PERIODIC: \
             /* mod operator can give pos or neg results */ \
             /* also, note that kernel could be larger than matrix*/ \
             row_input = row_input % in_nrow; \
             if ( row_input < 0 ) { \
               row_input = in_nrow + row_input; \
             } \
             col_input = col_input % in_ncol; \
             if ( col_input < 0 ) { \
               col_input = in_ncol + col_input; \
             } \
             inval = pinsig[ MATANY_INDEX( (SIG_IN_PTR), row_input, \
               col_input ) ]; \
             break; \
             \
           case MUTIL_BOUNDARY_REFLECT: \
             /* note that kernel could be larger than matrix*/ \
             while( row_input < 0 || row_input > in_nrow - 1 ) { \
               if ( row_input < 0 ) { \
                 row_input =  -1 - row_input; \
               } \
               if ( row_input > in_nrow - 1 ) { \
                 row_input = in_nrow - 1 - (row_input - in_nrow); \
               } \
             } \
             while( col_input < 0 || col_input > in_ncol - 1 ) { \
               if ( col_input < 0 ) { \
                 col_input =  -1 - col_input; \
               } \
               if ( col_input > in_ncol - 1 ) { \
                 col_input = in_ncol - 1 - (col_input - in_ncol); \
               } \
             } \
             inval = pinsig[ MATANY_INDEX( (SIG_IN_PTR), row_input, \
               col_input ) ]; \
             break; \
           \
           default: \
             MUTIL_ERROR( "Boundary condition not supported" ); \
             return MUTIL_ERR_ILLEGAL_VALUE; \
         } /* end of switch on boundary condition */ \
         \
         if( CONVOLUTION ) { \
           /* flip and transpose the kernel */ \
           poutsig[ outsig_idx ] += \
             inval * pkern[MATANY_INDEX( (KERNEL_PTR), \
               kern_nrow - 1 - row_kernel, kern_ncol - 1 - col_kernel )]; \
         } \
         else { \
           /* correlation */ \
           poutsig[ outsig_idx ] += \
             inval * pkern[MATANY_INDEX( (KERNEL_PTR), \
               row_kernel, col_kernel )]; \
         } \
       } /* loop over kernel cols */ \
     } /* loop over kernel rows */ \
     outsig_idx++;\
   } /* loop over input matrix cols */ \
   \
   num_ops += 4.0 * kern_nrow * kern_ncol * in_ncol; \
   if( MUTIL_INTERRUPT( num_ops, INTRP_PTR )) { \
     MUTIL_ERROR( "User interrupt" ); \
     return MUTIL_ERR_INTERRUPT; \
   } \
 } /* loop over input matrix rows */ \
 \
 if( CONVOLUTION ) { \
   MUTIL_TRACE( "Done with sig" #SIG_FN_PREFIX "_convolve()" ); \
 } \
 else { \
   MUTIL_TRACE( "Done with sig" #SIG_FN_PREFIX "_correlate()" ); \
 } \
 return MUTIL_ERR_OK


 /* function documented in sig_conv.h */
 /* written by Vikram Chalana */
mutil_errcode sigdbl_convolve( const double_mat *in_sig,
    const double_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, double_mat *out_sig )
{
    boolean convolve = TRUE;

    TMPL_SIG_CONVOLVE( dbl, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, double, convolve );
}


/* function documented in sig_conv.h */
/* written by Andrea Borning */
mutil_errcode sigflt_convolve( const float_mat *in_sig,
    const float_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, float_mat *out_sig )
{
    boolean convolve = TRUE;

    TMPL_SIG_CONVOLVE( flt, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, float, convolve );
}


/* function documented in sig_conv.h */
/* written by Qin Cai */
mutil_errcode sigu8_convolve( const uint8_mat *in_sig,
    const uint8_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, uint8_mat *out_sig )
{
    boolean convolve = TRUE;

    TMPL_SIG_CONVOLVE( u8, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, uint8, convolve );
}


/* function documented in sig_conv.h */
/* written by Qin Cai */
mutil_errcode sigu16_convolve( const uint16_mat *in_sig,
    const uint16_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, uint16_mat *out_sig )
{
    boolean convolve = TRUE;

    TMPL_SIG_CONVOLVE( u16, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, uint16, convolve );
}


/* function documented in sig_conv.h */
/* written by Qin Cai */
mutil_errcode sigs32_convolve( const sint32_mat *in_sig,
    const sint32_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, sint32_mat *out_sig )
{
    boolean convolve = TRUE;

    TMPL_SIG_CONVOLVE( s32, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, sint32, convolve );
}


/* function documented in sig_conv.h */
/* written by Qin Cai */
mutil_errcode sigu32_convolve( const uint32_mat *in_sig,
    const uint32_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, uint32_mat *out_sig)
{
    boolean convolve = TRUE;

    TMPL_SIG_CONVOLVE( u32, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, uint32, convolve );
}


/* function documented in sig_conv.h */
/* written by Vikram Chalana */
mutil_errcode siguniv_correlate( const univ_mat *in_sig,
    const univ_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, univ_mat *out_sig )
{
  mutil_errcode err;
  univ_mat      cast_in_sig;
  univ_mat      cast_kernel;
  boolean       alloc_cast_in;
  boolean       alloc_cast_ker;

  MUTIL_TRACE( "Start siguniv_correlate()" );

  if( !in_sig || !kernel || !out_sig ) {
    MUTIL_ERROR( "NULL pointer for operand or result" );
    return MUTIL_ERR_NULL_POINTER;
  }

  /* do standard casting */

  MUTIL_DO_STANDARD_CASTING( in_sig, out_sig->type, cast_in_sig,
    intrp_ptr, alloc_cast_in, err );
  if( err ) {
    return err;
  }

  MUTIL_DO_STANDARD_CASTING( kernel, out_sig->type, cast_kernel,
    intrp_ptr, alloc_cast_ker, err );
  if( err ) {
    MUTIL_FREE_STANDARD_CASTING( cast_in_sig, alloc_cast_in );
    return err;
  }

  /* finally call the appropriate function by
     Switching on the type of the output matrix */
  switch( out_sig->type ) {
     case MUTIL_DCOMPLEX:
       err = sigcpx_correlate( &(cast_in_sig.mat.cpxmat),
           &(cast_kernel.mat.cpxmat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.cpxmat) );
       break;

     case MUTIL_DOUBLE:
       err = sigdbl_correlate( &(cast_in_sig.mat.dblmat),
           &(cast_kernel.mat.dblmat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.dblmat) );
       break;

     case MUTIL_FLOAT:
       err = sigflt_correlate( &(cast_in_sig.mat.fltmat),
           &(cast_kernel.mat.fltmat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.fltmat) );
       break;

     case MUTIL_UINT8:
       err = sigu8_correlate( &(cast_in_sig.mat.u8mat),
           &(cast_kernel.mat.u8mat), row_step, col_step, row_overlap,
           col_overlap,boundary, intrp_ptr, &(out_sig->mat.u8mat) );
       break;

     case MUTIL_UINT16:
       err = sigu16_correlate( &(cast_in_sig.mat.u16mat),
           &(cast_kernel.mat.u16mat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr,
          &(out_sig->mat.u16mat) );
       break;

     case MUTIL_UINT32:
       err = sigu32_correlate( &(cast_in_sig.mat.u32mat),
           &(cast_kernel.mat.u32mat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.u32mat) );
       break;

     case MUTIL_SINT32:
       err = sigs32_correlate( &(cast_in_sig.mat.s32mat),
           &(cast_kernel.mat.s32mat), row_step, col_step, row_overlap,
           col_overlap, boundary, intrp_ptr, &(out_sig->mat.s32mat) );
       break;

    default:
      MUTIL_ERROR( "This matrix type is currently unsupported" );
      err = MUTIL_ERR_ILLEGAL_TYPE;
  }

  /* finally, free the cast matrices */
  MUTIL_FREE_STANDARD_CASTING( cast_in_sig, alloc_cast_in );
  MUTIL_FREE_STANDARD_CASTING( cast_kernel, alloc_cast_ker );

  if( err ) {
    return err;
  }

  MUTIL_TRACE( "Done with siguniv_correlate()" );
  return MUTIL_ERR_OK;
}


/* function documented in sig_conv.h */
/* written by Vikram Chalana */
mutil_errcode sigdbl_correlate( const double_mat *in_sig,
    const double_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, double_mat *out_sig )
{
    boolean convolve = FALSE;

    TMPL_SIG_CONVOLVE( dbl, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, double, convolve );
}


/* function documented in sig_conv.h */
/* written by Andrea Borning */
mutil_errcode sigflt_correlate( const float_mat *in_sig,
    const float_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, float_mat *out_sig )
{
    boolean convolve = FALSE;

    TMPL_SIG_CONVOLVE( flt, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, float, convolve );
}


/* function documented in sig_conv.h */
/* written by Qin Cai */
mutil_errcode sigu8_correlate( const uint8_mat *in_sig,
    const uint8_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, uint8_mat *out_sig )
{
    boolean convolve = FALSE;

    TMPL_SIG_CONVOLVE( u8, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, uint8, convolve );
}


/* function documented in sig_conv.h */
/* written by Qin Cai */
mutil_errcode sigu16_correlate( const uint16_mat *in_sig,
    const uint16_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, uint16_mat *out_sig )
{
    boolean convolve = FALSE;

    TMPL_SIG_CONVOLVE( u16, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, uint16, convolve );
}


/* function documented in sig_conv.h */
/* written by Qin Cai */
mutil_errcode sigs32_correlate( const sint32_mat *in_sig,
    const sint32_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, sint32_mat *out_sig )
{
    boolean convolve = FALSE;

    TMPL_SIG_CONVOLVE( s32, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, sint32, convolve );
}


/* function documented in sig_conv.h */
/* written by Qin Cai */
mutil_errcode sigu32_correlate( const uint32_mat *in_sig,
    const uint32_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, uint32_mat *out_sig )
{
    boolean convolve = FALSE;

    TMPL_SIG_CONVOLVE( u32, in_sig, kernel, row_step, col_step, row_overlap,
        col_overlap, boundary, intrp_ptr, out_sig, uint32, convolve );
}


/* function documented in sig_conv.h */
/* Function written by Sean Shen in Cambridge; modified by
   Jennifer Hodgdon */
mutil_errcode sigcpx_correlate( const dcomplex_mat *in_sig,
    const dcomplex_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, dcomplex_mat *out_sig )
{
  sint32    in_nrow;        /* number of input rows */
  sint32    in_ncol;        /* number of input columns */
  sint32    kern_nrow;      /* number of kernel rows */
  sint32    kern_ncol;      /* number of kernel columns */
  sint32    kern_row_offset;
  sint32    kern_col_offset;
  sint32    out_nrow;       /* number of output rows */
  sint32    out_ncol;       /* number of output columns */

  sint32    row_output;     /* indices into matrices */
  sint32    col_output;
  sint32    row_kernel;
  sint32    col_kernel;
  sint32    row_input;
  sint32    col_input;
  sint32    outsig_idx;
  sint32    kernel_idx;
  sint32    row_input_offset;
  sint32    col_input_offset;

  mutil_errcode err;
  double        num_ops = 0;

  dcomplex    *pinsig;
  dcomplex    inval;
  dcomplex    *poutsig;
  dcomplex    *pkern;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start sigcpx_correlate()" );

  /* sanity checks */

  err = matcpx_validate( in_sig );
  if( err ) {
    return err;
  }
  err = matcpx_validate( kernel );
  if( err ) {
    return err;
  }
  err = matcpx_validate( out_sig );
  if( err ) {
    return err;
  }

  if( out_sig->nrow > in_sig->nrow + kernel->nrow - 1 ||
      out_sig->ncol > in_sig->ncol + kernel->ncol - 1 ) {
      MUTIL_ERROR( "Output dimensions too big" );
      return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if( (in_sig)->data == (out_sig)->data ) {
    MUTIL_ERROR( "Input image cannot be the same as the output image" );
    return MUTIL_ERR_ILLEGAL_ADDRESS;
  }

  if( (kernel)->data == (out_sig)->data ) {
    MUTIL_ERROR( "Kernel image cannot be the same as the output image" );
    return MUTIL_ERR_ILLEGAL_ADDRESS;
  }

  if( row_overlap < 1 || col_overlap < 1 ) {
    MUTIL_ERROR( "Overlap must be at least 1" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if( row_overlap > (kernel)->nrow ||
      col_overlap > (kernel)->ncol ) {
    MUTIL_ERROR( "Overlap must be less than the dimensions of the kernel" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* get the sizes */
  in_nrow   = in_sig->nrow;
  in_ncol   = in_sig->ncol;
  out_nrow  = out_sig->nrow;
  out_ncol  = out_sig->ncol;
  kern_nrow = kernel->nrow;
  kern_ncol = kernel->ncol;

  kern_row_offset = kern_nrow - row_overlap;
  kern_col_offset = kern_ncol - col_overlap;

  /* get the underlying flat arrays */
  pinsig  = in_sig->data;
  poutsig = out_sig->data;
  pkern   = kernel->data;

  /* Loop over the output image coordinates */
  for ( row_output = 0; row_output < out_nrow; row_output ++ ) {

      row_input_offset = ( row_output * row_step ) - kern_row_offset;
      outsig_idx       = row_output * out_ncol;

    for ( col_output = 0; col_output < out_ncol; col_output ++ ) {

      /* compute here to avoid recomputing inside inner loops */
      col_input_offset = ( col_output * col_step ) - kern_col_offset;

      /* Set initial value to zero */
      poutsig[ outsig_idx ].re = 0.0;
      poutsig[ outsig_idx ].im = 0.0;
      kernel_idx               = 0;

      /* Loop over the kernel coordinates */
      for ( row_kernel = 0;  row_kernel < kern_nrow; row_kernel ++ ) {
        row_input = row_input_offset + row_kernel;

        for ( col_kernel = 0;  col_kernel < kern_ncol; col_kernel ++ ) {
          col_input = col_input_offset + col_kernel;

           switch( boundary ) {

            case MUTIL_BOUNDARY_ZERO:

              if(( row_input < 0 ) ||
                 ( col_input < 0 ) ||
                 ( row_input > in_nrow - 1 ) ||
                 ( col_input > in_ncol - 1 )) {
                inval.re = 0.0;
                inval.im = 0.0;
              }
              else {
                inval = pinsig[ MATANY_INDEX( in_sig, row_input, col_input ) ];
              }

              break;

            case MUTIL_BOUNDARY_CONTINUE:

              if ( row_input < 0 ) {
                row_input = 0;
              }
              else if ( row_input > in_nrow - 1 ) {
                row_input = in_nrow - 1;
              }

              if ( col_input < 0 ) {
                col_input = 0;
              }
              else if ( col_input > in_ncol - 1 ) {
                col_input = in_ncol - 1;
              }

              inval = pinsig[ MATANY_INDEX( in_sig, row_input, col_input ) ];
              break;

            case MUTIL_BOUNDARY_PERIODIC:

             /* mod operator can give pos or neg results */
             /* also, note that kernel could be larger than input */

             row_input = row_input % in_nrow;
             if ( row_input < 0 ) {
               row_input = in_nrow + row_input;
             }

             col_input = col_input % in_ncol;
             if ( col_input < 0 ) {
               col_input = in_ncol + col_input;
             }

             inval = pinsig[ MATANY_INDEX( in_sig, row_input, col_input ) ];
             break;

            case MUTIL_BOUNDARY_REFLECT:

              /* note that kernel could be larger than input matrix */

              while( row_input < 0 || row_input > in_nrow - 1 ) {
                if ( row_input < 0 ) {
                  row_input =  -1 - row_input;
                }
                if ( row_input > in_nrow - 1 ) {
                  row_input = in_nrow - 1 - (row_input - in_nrow);
                }
              }

              while( col_input < 0 || col_input > in_ncol - 1 ) {
                if ( col_input < 0 ) {
                  col_input =  -1 - col_input;
                }
                if ( col_input > in_ncol - 1 ) {
                  col_input = in_ncol - 1 - (col_input - in_ncol);
                }
              }

              inval = pinsig[ MATANY_INDEX( in_sig, row_input, col_input ) ];
              break;

            default:
              MUTIL_ERROR( "Boundary condition not supported" );
              return MUTIL_ERR_ILLEGAL_VALUE;

          } /* end of switch on boundary conditions */

          /* multiply to get the value */
          poutsig[ outsig_idx ].re += inval.re * pkern[ kernel_idx ].re -
              inval.im * pkern[ kernel_idx ].im;
          poutsig[ outsig_idx ].im += inval.re * pkern[ kernel_idx ].im +
              inval.im * pkern[ kernel_idx ].re;

          /* incrementing index works because order of loops is equivalent
             to stepping through kernel matrix in row major order */
          kernel_idx++;

        } /* loop over kernel cols */
      } /* loop over kernel rows */

      outsig_idx++;

    } /* loop over output matrix cols */

    num_ops += 8.0 * kern_nrow * kern_ncol * in_ncol;
    if( MUTIL_INTERRUPT( num_ops, intrp_ptr )) {
      MUTIL_ERROR( "User interrupt" );
      return MUTIL_ERR_INTERRUPT;
    }

  } /* end of loop over output matrix rows */

  MUTIL_TRACE( "Done with sigcpx_correlate()" );
  return MUTIL_ERR_OK;

}


/* function documented in sig_conv.h */
/* written by Luca Cazzanti */
mutil_errcode sigcpx_convolve( const dcomplex_mat *in_sig,
    const dcomplex_mat *kernel, sint32 row_step, sint32 col_step,
    sint32 row_overlap, sint32 col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, dcomplex_mat *out_sig )
{
    mutil_errcode errcode;
    dcomplex_mat  flip_kernel;

    MUTIL_TRACE( "Start sigcpx_convolve()" );

    /* All error checking done in sigcpx_correlate.
       Here just check what is absolutely necessary */

    errcode = matcpx_validate( kernel );
    if( errcode ) {
        return errcode;
    }

    errcode = matcpx_malloc( &flip_kernel, kernel->nrow, kernel->ncol );
    if( errcode ) {
        MUTIL_ERROR( "Could not alllocate memory for flipped kernel" );
        return errcode;
    }

    /* flip kernel */

    errcode = matcpx_flip_up_down( kernel, intrp_ptr, &flip_kernel );
    if( errcode ) {
        MUTIL_ERROR( "Could not flip kernel up/down" );
        MUTIL_FREE_WARN( matcpx, &flip_kernel );
        return errcode;
    }

    errcode = matcpx_flip_left_right( &flip_kernel, intrp_ptr, &flip_kernel );
    if( errcode ) {
        MUTIL_ERROR( "Could not flip kernel left/right" );
        MUTIL_FREE_WARN( matcpx, &flip_kernel );
        return errcode;
    }

    /* call correlation function */

    errcode = sigcpx_correlate( in_sig, &flip_kernel, row_step, col_step,
        row_overlap, col_overlap, boundary, intrp_ptr, out_sig );

    MUTIL_FREE_WARN( matcpx, &flip_kernel );

    if( errcode ) {
        MUTIL_ERROR( "Could not calculate convolution" );
        return errcode;
    }

    MUTIL_TRACE( "sigcpx_convolve() done" );
    return MUTIL_ERR_OK;
}


/* function documented in sig_conv.h */
/* written by Luca Cazzanti */
mutil_errcode sigset_correlate( const mat_set *in_set,
    const mat_set *kernel_set, const sint32_mat *row_step,
    const sint32_mat *col_step, const sint32_mat *row_overlap,
    const sint32_mat *col_overlap,mutil_boundary_type boundary,
    void *intrp_ptr, mat_set *out_set )
{
    mutil_errcode errcode;
    sint32        idx;

    MUTIL_TRACE( "Start sigset_correlate()" );

    /* check arguments */

    errcode = localfn_sigset_correlate_convolve_check_arguments( in_set,
        kernel_set, row_step, col_step, row_overlap, row_overlap, out_set );
    if( errcode ) {
        return errcode;
    }

    /* compute */

    for( idx = 0; idx < in_set->nelem; idx++ ) {
        errcode = siguniv_correlate( &(in_set->mats[ idx ]),
            &(kernel_set->mats[ idx ]), row_step->data[ idx ],
            col_step->data[ idx ], row_overlap->data [idx ],
            col_overlap->data[ idx ], boundary, intrp_ptr,
            &(out_set->mats[ idx ]) );
        if( errcode ) {
            MUTIL_ERROR( "Could not compute correlations" );
            return errcode;
        }
    }

    MUTIL_TRACE( "sigset_correlate() done" );
    return MUTIL_ERR_OK;
}



/* function documented in sig_conv.h */
/* written by Luca Cazzanti */
mutil_errcode sigset_convolve( const mat_set *in_set,
    const mat_set *kernel_set, const sint32_mat *row_step,
    const sint32_mat *col_step, const sint32_mat *row_overlap,
    const sint32_mat *col_overlap, mutil_boundary_type boundary,
    void *intrp_ptr, mat_set *out_set )
{
    mutil_errcode errcode;
    sint32        idx;

    MUTIL_TRACE( "Start sigset_convolve()" );

    /* check arguments */

    errcode = localfn_sigset_correlate_convolve_check_arguments( in_set,
        kernel_set, row_step, col_step, row_overlap, col_overlap, out_set );
    if( errcode ) {
        return errcode;
    }

    /* compute */

    for( idx = 0; idx < in_set->nelem; idx++ ) {
        errcode = siguniv_convolve( &(in_set->mats[ idx ]),
            &(kernel_set->mats[ idx ]), row_step->data[ idx ],
            col_step->data[ idx ], row_overlap->data[ idx ],
            col_overlap->data[ idx ], boundary, intrp_ptr,
            &(out_set->mats[ idx ]) );
        if( errcode ) {
            MUTIL_ERROR( "Could not compute convolutions" );
            return errcode;
        }
    }

    MUTIL_TRACE( "sigset_convolve() done" );
    return MUTIL_ERR_OK;
}

/*
****************
Static functions
****************
*/

/** Check the arguments for matrix set correlation and convolution.
 * This function checks the arguments passed to \Ref{sigset_convolve}
 * and \Ref{sigset_correlate}.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source sig\_conv.c
 * @library signal
 * @usage #err_code = localfn_sigset_correlate_convolve_check_arguments( &in_set, &kernel_set, &row_step, &col_step, &out_set );#
 * @return            Standard mutils error/OK code
 * @param in_set      Pointer to input matrix set.
 * @param kernel_set  Pointer to matrix set containing kernels.
 * @param row_step    Pointer to matrix of sint32 containing the downsampling
 *   step along the rows for each convolution/correlation.
 * @param col_step    Pointer to matrix of sint32 containing the downsampling
 *   step along the columns for each convolution/correlation.
 * @param row_overlap Pointer to matrix of sint32 containing the row overlap
 *  for each convolution/correlation.
 * @param col_overlap Pointer to matrix of sint32 containing the column
 *   overlap for each convolution/correlation.
 * @param out_set     Pointer to output matrix set.
 * @private
 */
static mutil_errcode localfn_sigset_correlate_convolve_check_arguments(
    const mat_set
    *in_set, const mat_set *kernel_set, const sint32_mat *row_step,
    const sint32_mat *col_step, const sint32_mat *row_overlap,
    const sint32_mat *col_overlap, const mat_set *out_set )
{
    mutil_errcode errcode;

    MUTIL_TRACE( "Start localfn_sigset_correlate_convolve_check_arguments()" );

    errcode = matset_validate( in_set );
    if( errcode ) {
        return errcode;
    }
    errcode = matset_validate( kernel_set );
    if( errcode ) {
        return errcode;
    }
    errcode = matset_validate( out_set );
    if( errcode ) {
        return errcode;
    }

    errcode = mats32_validate( row_step );
    if( errcode ) {
        return errcode;
    }
    errcode = mats32_validate( col_step );
    if( errcode ) {
        return errcode;
    }

    errcode = mats32_validate( row_overlap );
    if( errcode ) {
        return errcode;
    }
    errcode = mats32_validate( col_overlap );
    if( errcode ) {
        return errcode;
    }

    if( in_set->nelem != kernel_set->nelem ) {
        MUTIL_ERROR( "Input/kernel sets number of elements mismatch" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }
    if( in_set->nelem != out_set->nelem ) {
        MUTIL_ERROR( "Input/output sets number of elements mismatch" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }
    if( row_step->nelem != col_step->nelem ) {
        MUTIL_ERROR( "Row steps/columns steps length mismatch" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }
    if( row_step->nelem != in_set->nelem ) {
        MUTIL_ERROR( "Matrix of step values/input matrix set size mismatch" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    if( row_overlap->nelem != col_overlap->nelem ) {
        MUTIL_ERROR( "Row overlaps/column overlaps length mismatch" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }
    if( row_overlap->nelem != in_set->nelem ) {
        MUTIL_ERROR( "Matrix of overlap values/input matrix set size mismatch" );
        return MUTIL_ERR_ILLEGAL_SIZE;
    }

    MUTIL_TRACE( "localfn_sigset_correlate_convolve_check_arguments() done" );
    return MUTIL_ERR_OK;
}



