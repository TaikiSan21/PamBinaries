
static char whatssi[] = "@(#) $File: //depot/Research/ifultools/pkg/ifultools/src/wav_filt.c $: $Revision: #1 $, $Date: 2008/03/21 $";

#include "wav_filt.h"
#include "wav_type.h"

#include "mat_assn.h"
#include "mat_set.h"
#include "mat_set.h"
#include "mat_type.h"
#include "mat_umat.h"
#include "mat_univ.h"

#include "sig_tran.h"

#include "ut_alloc.h"
#include "ut_debug.h"
#include "ut_intrn.h"
#include "ut_intrp.h"
#include "ut_limit.h"
#include "ut_math.h"
#include "ut_mem.h"

#include <math.h>

/******************************************************/
/* STATIC FILTER GENERATION FUNCTION DECLARATIONS     */
/******************************************************/

static mutil_errcode localfn_wavuniv_filters_daubechies_input_check(
  sint32           filter_length,
  wav_filter_type  filter_type );

static mutil_errcode localfn_filter_matrices_check(
  sint32    filter_length,
  const univ_mat *wavelet_filter,
  const univ_mat *scaling_filter );

static mutil_errcode localfn_filter_fill(
  sint32  filter_length,
  double *in,
  double *out );

/******************************************************/
/* STATIC FILTER INDEXING FUNCTION DECLARATIONS       */
/******************************************************/

static mutil_errcode localfn_wavuniv_transform_coefficient_boundaries_inputs_check(
  sint32         n_level,
  sint32         filter_length,
  sint32         n_sample,
  wav_transform  transform_type );

/******************************************************/
/* STATIC FILTER PHASE FUNCTION DECLARATIONS          */
/******************************************************/

static sint32 localfn_S_value(
  sint32 level,
  sint32 node );

/******************************************************/
/* STATIC FILTER GENERATION MACRO DEFINITIONS         */
/******************************************************/

#undef  LOCALDEF_UNSUPPORTED_FILTER_LENGTH
#define LOCALDEF_UNSUPPORTED_FILTER_LENGTH     \
MUTIL_ERROR( "Filter length is unsupported" ); \
err = MUTIL_ERR_FEATURE_NOT_IMPLEMENTED

#undef LOCALDEF_MAX
#define LOCALDEF_MAX(A,B) ( ((A) > (B)) ? A : B )

#undef LOCALDEF_MODULO
#define LOCALDEF_MODULO(A,B) ( (A) - floor( (A) / (B) ) * (B) )

#undef LOCALDEF_IS_EVEN
#define LOCALDEF_IS_EVEN(N) ( (N % 2) == 0 ? 1 : 0 )

#undef LOCALDEF_MORLET_CONSTANT
#define LOCALDEF_MORLET_CONSTANT(W)                    \
(double) MUTIL_POW( MUTIL_PI, - 1.0 / 4.0 ) /                   \
sqrt( 1.0 - 4.0 / sqrt(3.0) * exp( - (W)*(W) / 4.0 ) + \
sqrt(2.0) * exp( - (W)*(W) / 2.0 ) )


/******************************************************/
/* STATIC FILTER INDEXING MACRO DECLARATIONS          */
/******************************************************/

/* NONE */

/******************************************************/
/* STATIC FILTER PHASE MACRO DEFINITIONS              */
/******************************************************/

#define DWT_SHIFT( SHIFT, LEVEL ) \
- (sint32) ceil( (double) ( MUTIL_ABS( SHIFT ) + 1 ) / \
(double) MUTIL_POW( 2, LEVEL ) - 1.0 )


/******************************************************/
/* FILTER GENERATION FUNCTION DEFINITIONS             */
/******************************************************/

/* Generate frequency response for   */
/* continuous wavelet filters        */
/* Function documented in wav_filt.h */
/* Written by William Constantine    */

mutil_errcode wavuniv_filters_continuous(
  const wav_filter_type  filter_type,
  const double           filter_arg,
  const univ_mat        *frequency,
  void                  *intrp_ptr,
  univ_mat              *result )
{
  dcomplex     *pc_result;
  double        Omega2;
  double        ampweight;
  double       *pd_freq;
  sint32        i;
  sint32        num_freqs;
  double        std;
  double        omega0;
  double        num_ops = 0.0;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_filters_continuous() ::" );

  /* avoid lint warning */

  (void) whatssi;

  /* check input ... */

  /* filter type */

  switch( filter_type ){
    case WAV_FILTER_GAUSSIAN_I:
    case WAV_FILTER_GAUSSIAN_II:
    case WAV_FILTER_MORLET:
    case WAV_FILTER_HAAR:
      break;
    default:
      MUTIL_ERROR( "Filter type is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }
  /* frequency matrix */

  if ( frequency == (univ_mat *) NULL ){
    MUTIL_ERROR( "Pointer to frequency matrix is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if ( frequency->type != MUTIL_DOUBLE ){
    MUTIL_ERROR( "Frequency matrix must be of type MUTIL_DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( !MATANY_IS_VEC( &(frequency->mat.dblmat) ) ){
    MUTIL_ERROR( "Frequency matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* output matrix */

  if ( result == (univ_mat *) NULL ){
    MUTIL_ERROR( "Pointer to frequency response matrix is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  if ( result->type != MUTIL_DCOMPLEX ){
    MUTIL_ERROR( "Frequency response matrix must be of type MUTIL_DCOMPLEX." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }

  if ( !MATANY_IS_VEC( &(result->mat.cpxmat) ) ){
    MUTIL_ERROR( "Frequency response matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  if ( MATUNIV_NELEM( result ) != MATUNIV_NELEM( frequency ) ){
    MUTIL_ERROR( "Frequency response matrix must have the same number of elements as does the input frequency matrix." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* filter specific input */

  switch( filter_type ){

  case WAV_FILTER_GAUSSIAN_I:
  case WAV_FILTER_GAUSSIAN_II:

    std = filter_arg;

    if ( std <= 0.0 ){
      MUTIL_ERROR( "Standard deviation of Gaussian function must be a positive value." );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
    
    break;
    
  case WAV_FILTER_MORLET:
    
    omega0 = filter_arg;
    
    if ( omega0 <= 0.0 ){
      MUTIL_ERROR( "Morlet frequency shift omega0 must be a positive value." );
      return MUTIL_ERR_ILLEGAL_VALUE;
    }
    
    break;

  case WAV_FILTER_HAAR:
    break;
    
  default:
    MUTIL_ERROR( "Filter type is unsupported" );
    return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* end I/O check */

  /** start main code */

  /* define weights used to form response functions */

  switch( filter_type ){

    case WAV_FILTER_HAAR:
      ampweight = sqrt( 2.0 );
      break;

    case WAV_FILTER_GAUSSIAN_I:
      ampweight = - 2.0 * pow( std, 1.5 ) * pow( MUTIL_PI, 1.0 / 4.0 );
      break;

    case WAV_FILTER_GAUSSIAN_II:
      ampweight = sqrt( 8.0 / 3.0 ) * pow( std, 2.5 ) * pow( MUTIL_PI, 1.0 / 4.0 );
      break;

    case WAV_FILTER_MORLET:

      ampweight = sqrt( 2.0 * MUTIL_PI ) * LOCALDEF_MORLET_CONSTANT( omega0 );
      break;

    default:
      MUTIL_ERROR( "Filter type is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /* set pointers */

  pd_freq = frequency->mat.dblmat.data;
  pc_result = result->mat.cpxmat.data;

  num_freqs = (sint32) MATUNIV_NELEM( frequency );

  /* develop mother wavelet filter in the frequency domain */

  for ( i = 0; i < num_freqs; i ++ ){

    /* Check for interrupts */
    num_ops += (double) num_freqs;
    
    if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      return MUTIL_ERR_INTERRUPT;
    }
    switch( filter_type ){
      
    case WAV_FILTER_HAAR:

      pc_result[i].re = (double) 0.0;
      
      if ( MUTIL_ABS( pd_freq[i] ) < MUTIL_DOUBLE_EPSILON ){
	  pc_result[i].im = (double) 0.0;
      }
      else{
	pc_result[i].im = ampweight / pd_freq[i] * ( cos( pd_freq[i] ) - 1.0 );
      }
      break;
	
    case WAV_FILTER_GAUSSIAN_I:
      
      pc_result[i].re = (double) 0.0;
      pc_result[i].im = ampweight * pd_freq[i] * exp( - pd_freq[i] * pd_freq[i] * std * std / 2.0 );
      break;
	
    case WAV_FILTER_GAUSSIAN_II:
	
      pc_result[i].re = ampweight * pd_freq[i] * pd_freq[i] * exp( - pd_freq[i] * pd_freq[i] * std * std  / 2.0 );
      pc_result[i].im = (double) 0.0;
      break;
      
    case WAV_FILTER_MORLET:
      
      /*  morlet = C * sqrt( 2 * PI ) * exp( -W^2 / 2 ) * ( 1 - exp( ( W^2 - W0^2 ) / 4 ) )
	  where W = w + W0
      */
      
      Omega2 = pd_freq[i] + omega0;
      
      /* here we check the relation | w + w0 | > 38.0. If TRUE, then
	 exp( ( W^2 - W0^2 ) / 4 ) is huge while exp( -W^2 / 2 ) is
	 (asymptotically) zero. So, in these cases, we simply return
	 a 0 + 0i value, otherwise we will have floating point overflow
	 issues. */
      
      if ( MUTIL_ABS( Omega2 ) > 38.0 ){
	
	pc_result[i].re = (double) 0.0;
	pc_result[i].im = (double) 0.0;
      }
      else{
	
	Omega2 *= Omega2;
	
	pc_result[i].re = ampweight * exp( - Omega2 / 2.0  ) * ( 1.0 - exp( ( Omega2 - omega0 * omega0 ) / 4.0 ) ) ;
	pc_result[i].im = (double) 0.0;
      }

      break;

    default:
      MUTIL_ERROR( "Filter type is unsupported" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;

    }

  }

  MUTIL_TRACE( "Done with wavuniv_filters_continuous()" );

  return MUTIL_ERR_OK;
}

/* Generate Daubechies wavelet and   */
/* scaling filters.                  */
/* Function documented in wav_filt.h */
/* Written by William Constantine    */

mutil_errcode wavuniv_filters_daubechies(
  sint32           filter_length,
  wav_filter_type  filter_type,
  boolean          normalize,
  void            *intrp_ptr,
  mat_set         *result )
{
  double       *pd_scaling;
  double       *pd_wavelet;
  double        num_ops = 0.0;
  mutil_errcode err;
  sint32        tap;
  memlist       list;
  sint32        ndim = 1;
  sint32        dims = 2;

  /* explicitly set common filters for lengths 2, 4, and 6 */

  double daub2[] = { 0.7071067811865475, 0.7071067811865475};
  double daub4[] = { 0.4829629131445341, 0.8365163037378077,
    0.2241438680420134, -0.1294095225512603 };
  double daub6[] = { 0.3326705529500827, 0.8068915093110928,
    0.4598775021184915, -0.1350110200102546,
    -0.0854412738820267, 0.0352262918857096 };

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_filters_daubechies() ::" );

  /* avoid lint warning */

  (void) whatssi;

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check input arguments */

  err = localfn_wavuniv_filters_daubechies_input_check(
    filter_length,
    filter_type );
  if ( err ) return err;

  /* allocate space for output matrix set */

  err = matset_malloc_register( result, ndim, &dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices( result, filter_length,
    1, MUTIL_DOUBLE );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* assign pointers for typing ease in setting */

  pd_wavelet = result->mats[ 0 ].mat.dblmat.data;
  pd_scaling = result->mats[ 1 ].mat.dblmat.data;

  /* create the enum types for the different wavelets */

  switch( filter_type ){

  case WAV_FILTER_HAAR:

      err = localfn_filter_fill( filter_length, daub2, pd_scaling );
      break;

  case WAV_FILTER_EXTREMAL_PHASE:

    switch ( filter_length ){

    case 2:   /* Haar */

      err = localfn_filter_fill( filter_length, daub2, pd_scaling );
      break;

    case 4:  /* D(4) */

      err = localfn_filter_fill( filter_length, daub4, pd_scaling );
      break;

    case 6:  /* D(6) */

      err = localfn_filter_fill( filter_length, daub6, pd_scaling );
      break;

    case 8:  /* D(8) */

      pd_scaling[ 0 ]  =  0.2303778133074431;
      pd_scaling[ 1 ]  =  0.7148465705484058;
      pd_scaling[ 2 ]  =  0.6308807679358788;
      pd_scaling[ 3 ]  = -0.0279837694166834;
      pd_scaling[ 4 ]  = -0.1870348117179132;
      pd_scaling[ 5 ]  =  0.0308413818353661;
      pd_scaling[ 6 ]  =  0.0328830116666778;
      pd_scaling[ 7 ]  = -0.0105974017850021;
      break;

    case 10:  /* D(10) */

      pd_scaling[ 0 ]  =  0.1601023979741930;
      pd_scaling[ 1 ]  =  0.6038292697971898;
      pd_scaling[ 2 ]  =  0.7243085284377729;
      pd_scaling[ 3 ]  =  0.1384281459013204;
      pd_scaling[ 4 ]  = -0.2422948870663824;
      pd_scaling[ 5 ]  = -0.0322448695846381;
      pd_scaling[ 6 ]  =  0.0775714938400459;
      pd_scaling[ 7 ]  = -0.0062414902127983;
      pd_scaling[ 8 ]  = -0.0125807519990820;
      pd_scaling[ 9 ]  =  0.0033357252854738;
      break;

    case 12:  /* D(12) */

      pd_scaling[ 0 ]  =  0.1115407433501094;
      pd_scaling[ 1 ]  =  0.4946238903984530;
      pd_scaling[ 2 ]  =  0.7511339080210954;
      pd_scaling[ 3 ]  =  0.3152503517091980;
      pd_scaling[ 4 ]  = -0.2262646939654399;
      pd_scaling[ 5 ]  = -0.1297668675672624;
      pd_scaling[ 6 ]  =  0.0975016055873224;
      pd_scaling[ 7 ]  =  0.0275228655303053;
      pd_scaling[ 8 ]  = -0.0315820393174862;
      pd_scaling[ 9 ]  =  0.0005538422011614;
      pd_scaling[ 10 ] =  0.0047772575109455;
      pd_scaling[ 11 ] = -0.0010773010853085;
      break;

    case 14:  /* D(14) */

      pd_scaling[ 0 ]  =  0.0778520540850081;
      pd_scaling[ 1 ]  =  0.3965393194819136;
      pd_scaling[ 2 ]  =  0.7291320908462368;
      pd_scaling[ 3 ]  =  0.4697822874052154;
      pd_scaling[ 4 ]  = -0.1439060039285293;
      pd_scaling[ 5 ]  = -0.2240361849938538;
      pd_scaling[ 6 ]  =  0.0713092192668312;
      pd_scaling[ 7 ]  =  0.0806126091510820;
      pd_scaling[ 8 ]  = -0.0380299369350125;
      pd_scaling[ 9 ]  = -0.0165745416306664;
      pd_scaling[ 10 ] =  0.0125509985560993;
      pd_scaling[ 11 ] =  0.0004295779729214;
      pd_scaling[ 12 ] = -0.0018016407040474;
      pd_scaling[ 13 ] =  0.0003537137999745;
      break;

    case 16:  /* D(16) */

      pd_scaling[ 0 ]  =  0.0544158422431049;
      pd_scaling[ 1 ]  =  0.3128715909143031;
      pd_scaling[ 2 ]  =  0.6756307362972904;
      pd_scaling[ 3 ]  =  0.5853546836541907;
      pd_scaling[ 4 ]  = -0.0158291052563816;
      pd_scaling[ 5 ]  = -0.2840155429615702;
      pd_scaling[ 6 ]  =  0.0004724845739124;
      pd_scaling[ 7 ]  =  0.1287474266204837;
      pd_scaling[ 8 ]  = -0.0173693010018083;
      pd_scaling[ 9 ]  = -0.0440882539307952;
      pd_scaling[ 10 ] =  0.0139810279173995;
      pd_scaling[ 11 ] =  0.0087460940474061;
      pd_scaling[ 12 ] = -0.0048703529934518;
      pd_scaling[ 13 ] = -0.0003917403733770;
      pd_scaling[ 14 ] =  0.0006754494064506;
      pd_scaling[ 15 ] = -0.0001174767841248;
      break;

    case 18:  /* D(18) */

      pd_scaling[ 0 ]  =  0.0380779473638791;
      pd_scaling[ 1 ]  =  0.2438346746125939;
      pd_scaling[ 2 ]  =  0.6048231236901156;
      pd_scaling[ 3 ]  =  0.6572880780512955;
      pd_scaling[ 4 ]  =  0.1331973858249927;
      pd_scaling[ 5 ]  = -0.2932737832791761;
      pd_scaling[ 6 ]  = -0.0968407832229524;
      pd_scaling[ 7 ]  =  0.1485407493381306;
      pd_scaling[ 8 ]  =  0.0307256814793395;
      pd_scaling[ 9 ]  = -0.0676328290613302;
      pd_scaling[ 10 ] =  0.0002509471148340;
      pd_scaling[ 11 ] =  0.0223616621236805;
      pd_scaling[ 12 ] = -0.0047232047577520;
      pd_scaling[ 13 ] = -0.0042815036824636;
      pd_scaling[ 14 ] =  0.0018476468830564;
      pd_scaling[ 15 ] =  0.0002303857635232;
      pd_scaling[ 16 ] = -0.0002519631889427;
      pd_scaling[ 17 ] =  0.0000393473203163;
      break;

    case 20:  /* D(20) */

      pd_scaling[ 0 ]  =  0.0266700579005546;
      pd_scaling[ 1 ]  =  0.1881768000776863;
      pd_scaling[ 2 ]  =  0.5272011889317202;
      pd_scaling[ 3 ]  =  0.6884590394536250;
      pd_scaling[ 4 ]  =  0.2811723436606485;
      pd_scaling[ 5 ]  = -0.2498464243272283;
      pd_scaling[ 6 ]  = -0.1959462743773399;
      pd_scaling[ 7 ]  =  0.1273693403357890;
      pd_scaling[ 8 ]  =  0.0930573646035802;
      pd_scaling[ 9 ]  = -0.0713941471663697;
      pd_scaling[ 10 ] = -0.0294575368218480;
      pd_scaling[ 11 ] =  0.0332126740593703;
      pd_scaling[ 12 ] =  0.0036065535669880;
      pd_scaling[ 13 ] = -0.0107331754833036;
      pd_scaling[ 14 ] =  0.0013953517470692;
      pd_scaling[ 15 ] =  0.0019924052951930;
      pd_scaling[ 16 ] = -0.0006858566949566;
      pd_scaling[ 17 ] = -0.0001164668551285;
      pd_scaling[ 18 ] =  0.0000935886703202;
      pd_scaling[ 19 ] = -0.0000132642028945;
      break;

    default:

      LOCALDEF_UNSUPPORTED_FILTER_LENGTH;

      break;
    }

    break;

  case WAV_FILTER_LEAST_ASYMMETRIC:

    switch ( filter_length ){

    case 2:   /* LA(2) */

      err = localfn_filter_fill( filter_length, daub2, pd_scaling );
      break;

    case 4:  /* LA(4) */

      err = localfn_filter_fill( filter_length, daub4, pd_scaling );
      break;

    case 6:  /* LA(6) */

      err = localfn_filter_fill( filter_length, daub6, pd_scaling );
      break;

    case 8:  /* LA(8) */

      pd_scaling[ 0 ]  = -0.0757657147893407;
      pd_scaling[ 1 ]  = -0.0296355276459541;
      pd_scaling[ 2 ]  =  0.4976186676324578;
      pd_scaling[ 3 ]  =  0.8037387518052163;
      pd_scaling[ 4 ]  =  0.2978577956055422;
      pd_scaling[ 5 ]  = -0.0992195435769354;
      pd_scaling[ 6 ]  = -0.0126039672622612;
      pd_scaling[ 7 ]  =  0.0322231006040713;
      break;

    case 10:  /* LA(10) */

      pd_scaling[ 0 ]  =  0.0195388827353869;
      pd_scaling[ 1 ]  = -0.0211018340249298;
      pd_scaling[ 2 ]  = -0.1753280899081075;
      pd_scaling[ 3 ]  =  0.0166021057644243;
      pd_scaling[ 4 ]  =  0.6339789634569490;
      pd_scaling[ 5 ]  =  0.7234076904038076;
      pd_scaling[ 6 ]  =  0.1993975339769955;
      pd_scaling[ 7 ]  = -0.0391342493025834;
      pd_scaling[ 8 ]  =  0.0295194909260734;
      pd_scaling[ 9 ]  =  0.0273330683451645;
      break;

    case 12:  /* LA(12) */

      pd_scaling[ 0 ]  =  0.0154041093273377;
      pd_scaling[ 1 ]  =  0.0034907120843304;
      pd_scaling[ 2 ]  = -0.1179901111484105;
      pd_scaling[ 3 ]  = -0.0483117425859981;
      pd_scaling[ 4 ]  =  0.4910559419276396;
      pd_scaling[ 5 ]  =  0.7876411410287941;
      pd_scaling[ 6 ]  =  0.3379294217282401;
      pd_scaling[ 7 ]  = -0.0726375227866000;
      pd_scaling[ 8 ]  = -0.0210602925126954;
      pd_scaling[ 9 ]  =  0.0447249017707482;
      pd_scaling[ 10 ] =  0.0017677118643983;
      pd_scaling[ 11 ] = -0.0078007083247650;
      break;

    case 14:  /* LA(14) */

      pd_scaling[ 0 ]  =  0.0102681767084968;
      pd_scaling[ 1 ]  =  0.0040102448717033;
      pd_scaling[ 2 ]  = -0.1078082377036168;
      pd_scaling[ 3 ]  = -0.1400472404427030;
      pd_scaling[ 4 ]  =  0.2886296317509833;
      pd_scaling[ 5 ]  =  0.7677643170045710;
      pd_scaling[ 6 ]  =  0.5361019170907720;
      pd_scaling[ 7 ]  =  0.0174412550871099;
      pd_scaling[ 8 ]  = -0.0495528349370410;
      pd_scaling[ 9 ]  =  0.0678926935015971;
      pd_scaling[ 10 ] =  0.0305155131659062;
      pd_scaling[ 11 ] = -0.0126363034031526;
      pd_scaling[ 12 ] = -0.0010473848889657;
      pd_scaling[ 13 ] =  0.0026818145681164;
      break;

    case 16:  /* LA(16) */

      pd_scaling[ 0 ]  = -0.0033824159513594;
      pd_scaling[ 1 ]  = -0.0005421323316355;
      pd_scaling[ 2 ]  =  0.0316950878103452;
      pd_scaling[ 3 ]  =  0.0076074873252848;
      pd_scaling[ 4 ]  = -0.1432942383510542;
      pd_scaling[ 5 ]  = -0.0612733590679088;
      pd_scaling[ 6 ]  =  0.4813596512592012;
      pd_scaling[ 7 ]  =  0.7771857516997478;
      pd_scaling[ 8 ]  =  0.3644418948359564;
      pd_scaling[ 9 ]  = -0.0519458381078751;
      pd_scaling[ 10 ] = -0.0272190299168137;
      pd_scaling[ 11 ] =  0.0491371796734768;
      pd_scaling[ 12 ] =  0.0038087520140601;
      pd_scaling[ 13 ] = -0.0149522583367926;
      pd_scaling[ 14 ] = -0.0003029205145516;
      pd_scaling[ 15 ] =  0.0018899503329007;
      break;

    case 18:  /* LA(18) */

      pd_scaling[ 0 ]  =  0.0010694900326538;
      pd_scaling[ 1 ]  = -0.0004731544985879;
      pd_scaling[ 2 ]  = -0.0102640640276849;
      pd_scaling[ 3 ]  =  0.0088592674935117;
      pd_scaling[ 4 ]  =  0.0620777893027638;
      pd_scaling[ 5 ]  = -0.0182337707798257;
      pd_scaling[ 6 ]  = -0.1915508312964873;
      pd_scaling[ 7 ]  =  0.0352724880359345;
      pd_scaling[ 8 ]  =  0.6173384491413523;
      pd_scaling[ 9 ]  =  0.7178970827642257;
      pd_scaling[ 10 ] =  0.2387609146074182;
      pd_scaling[ 11 ] = -0.0545689584305765;
      pd_scaling[ 12 ] =  0.0005834627463312;
      pd_scaling[ 13 ] =  0.0302248788579895;
      pd_scaling[ 14 ] = -0.0115282102079848;
      pd_scaling[ 15 ] = -0.0132719677815332;
      pd_scaling[ 16 ] =  0.0006197808890549;
      pd_scaling[ 17 ] =  0.0014009155255716;
      break;

    case 20:  /* LA(20) */

      pd_scaling[ 0 ]  =  0.0007701598091030;
      pd_scaling[ 1 ]  =  0.0000956326707837;
      pd_scaling[ 2 ]  = -0.0086412992759401;
      pd_scaling[ 3 ]  = -0.0014653825833465;
      pd_scaling[ 4 ]  =  0.0459272392237649;
      pd_scaling[ 5 ]  =  0.0116098939129724;
      pd_scaling[ 6 ]  = -0.1594942788575307;
      pd_scaling[ 7 ]  = -0.0708805358108615;
      pd_scaling[ 8 ]  =  0.4716906668426588;
      pd_scaling[ 9 ]  =  0.7695100370143388;
      pd_scaling[ 10 ] =  0.3838267612253823;
      pd_scaling[ 11 ] = -0.0355367403054689;
      pd_scaling[ 12 ] = -0.0319900568281631;
      pd_scaling[ 13 ] =  0.0499949720791560;
      pd_scaling[ 14 ] =  0.0057649120455518;
      pd_scaling[ 15 ] = -0.0203549398039460;
      pd_scaling[ 16 ] = -0.0008043589345370;
      pd_scaling[ 17 ] =  0.0045931735836703;
      pd_scaling[ 18 ] =  0.0000570360843390;
      pd_scaling[ 19 ] = -0.0004593294205481;
      break;

    default:

      LOCALDEF_UNSUPPORTED_FILTER_LENGTH;
      break;
    }

    break;

  case WAV_FILTER_BEST_LOCALIZED:

    switch ( filter_length ){

    case 2:   /* BL(2) */

      err = localfn_filter_fill( filter_length, daub2, pd_scaling );
      break;

    case 4:  /* BL(4) */

      err = localfn_filter_fill( filter_length, daub4, pd_scaling );
      break;

    case 6:  /* BL(6) */

      err = localfn_filter_fill( filter_length, daub6, pd_scaling );
      break;

    case 14:  /* BL(14) */

      pd_scaling[ 0 ]  =  0.0120154192834842;
      pd_scaling[ 1 ]  =  0.0172133762994439;
      pd_scaling[ 2 ]  = -0.0649080035533744;
      pd_scaling[ 3 ]  = -0.0641312898189170;
      pd_scaling[ 4 ]  =  0.3602184608985549;
      pd_scaling[ 5 ]  =  0.7819215932965554;
      pd_scaling[ 6 ]  =  0.4836109156937821;
      pd_scaling[ 7 ]  = -0.0568044768822707;
      pd_scaling[ 8 ]  = -0.1010109208664125;
      pd_scaling[ 9 ]  =  0.0447423494687405;
      pd_scaling[ 10 ] =  0.0204642075778225;
      pd_scaling[ 11 ] = -0.0181266051311065;
      pd_scaling[ 12 ] = -0.0032832978473081;
      pd_scaling[ 13 ] =  0.0022918339541009;
      break;

    case 18:  /* BL(18) */

      pd_scaling[ 0 ]  =  0.0002594576266544;
      pd_scaling[ 1 ]  = -0.0006273974067728;
      pd_scaling[ 2 ]  = -0.0019161070047557;
      pd_scaling[ 3 ]  =  0.0059845525181721;
      pd_scaling[ 4 ]  =  0.0040676562965785;
      pd_scaling[ 5 ]  = -0.0295361433733604;
      pd_scaling[ 6 ]  = -0.0002189514157348;
      pd_scaling[ 7 ]  =  0.0856124017265279;
      pd_scaling[ 8 ]  = -0.0211480310688774;
      pd_scaling[ 9 ]  = -0.1432929759396520;
      pd_scaling[ 10 ] =  0.2337782900224977;
      pd_scaling[ 11 ] =  0.7374707619933686;
      pd_scaling[ 12 ] =  0.5926551374433956;
      pd_scaling[ 13 ] =  0.0805670008868546;
      pd_scaling[ 14 ] = -0.1143343069619310;
      pd_scaling[ 15 ] = -0.0348460237698368;
      pd_scaling[ 16 ] =  0.0139636362487191;
      pd_scaling[ 17 ] =  0.0057746045512475;
      break;

    case 20:  /* BL(20) */

      pd_scaling[ 0 ]  =  0.0008625782242896;
      pd_scaling[ 1 ]  =  0.0007154205305517;
      pd_scaling[ 2 ]  = -0.0070567640909701;
      pd_scaling[ 3 ]  =  0.0005956827305406;
      pd_scaling[ 4 ]  =  0.0496861265075979;
      pd_scaling[ 5 ]  =  0.0262403647054251;
      pd_scaling[ 6 ]  = -0.1215521061578162;
      pd_scaling[ 7 ]  = -0.0150192395413644;
      pd_scaling[ 8 ]  =  0.5137098728334054;
      pd_scaling[ 9 ]  =  0.7669548365010849;
      pd_scaling[ 10 ] =  0.3402160135110789;
      pd_scaling[ 11 ] = -0.0878787107378667;
      pd_scaling[ 12 ] = -0.0670899071680668;
      pd_scaling[ 13 ] =  0.0338423550064691;
      pd_scaling[ 14 ] = -0.0008687519578684;
      pd_scaling[ 15 ] = -0.0230054612862905;
      pd_scaling[ 16 ] = -0.0011404297773324;
      pd_scaling[ 17 ] =  0.0050716491945793;
      pd_scaling[ 18 ] =  0.0003401492622332;
      pd_scaling[ 19 ] = -0.0004101159165852;
      break;

    default:

      LOCALDEF_UNSUPPORTED_FILTER_LENGTH;
      break;
    }

    break;

  case WAV_FILTER_COIFLET:

    switch ( filter_length ){

    case 6:  /* C(6) */

      pd_scaling[ 0 ]  = -0.01565572813579199;
      pd_scaling[ 1 ]  = -0.07273261951252648;
      pd_scaling[ 2 ]  =  0.38486484686485783;
      pd_scaling[ 3 ]  =  0.85257202021160061;
      pd_scaling[ 4 ]  =  0.33789766245748187;
      pd_scaling[ 5 ]  = -0.07273261951252648;
      break;

    case 12:  /* C(12) */

      pd_scaling[ 0 ]  = -0.0007205494453679;
      pd_scaling[ 1 ]  = -0.0018232088707116;
      pd_scaling[ 2 ]  =  0.0056114348194211;
      pd_scaling[ 3 ]  =  0.0236801719464464;
      pd_scaling[ 4 ]  = -0.0594344186467388;
      pd_scaling[ 5 ]  = -0.0764885990786692;
      pd_scaling[ 6 ]  =  0.4170051844236707;
      pd_scaling[ 7 ]  =  0.8127236354493977;
      pd_scaling[ 8 ]  =  0.3861100668229939;
      pd_scaling[ 9 ]  = -0.0673725547222826;
      pd_scaling[ 10 ] = -0.0414649367819558;
      pd_scaling[ 11 ] =  0.0163873364635998;
      break;

    case 18:  /* C(18) */

      pd_scaling[ 0 ]  = -0.0000345997728362;
      pd_scaling[ 1 ]  = -0.0000709833031381;
      pd_scaling[ 2 ]  =  0.0004662169601129;
      pd_scaling[ 3 ]  =  0.0011175187708906;
      pd_scaling[ 4 ]  = -0.0025745176887502;
      pd_scaling[ 5 ]  = -0.0090079761366615;
      pd_scaling[ 6 ]  =  0.0158805448636158;
      pd_scaling[ 7 ]  =  0.0345550275730615;
      pd_scaling[ 8 ]  = -0.0823019271068856;
      pd_scaling[ 9 ]  = -0.0717998216193117;
      pd_scaling[ 10 ] =  0.4284834763776168;
      pd_scaling[ 11 ] =  0.7937772226256169;
      pd_scaling[ 12 ] =  0.4051769024096150;
      pd_scaling[ 13 ] = -0.0611233900026726;
      pd_scaling[ 14 ] = -0.0657719112818552;
      pd_scaling[ 15 ] =  0.0234526961418362;
      pd_scaling[ 16 ] =  0.0077825964273254;
      pd_scaling[ 17 ] = -0.0037935128644910;
      break;

    case 24:  /* C(24) */

      pd_scaling[ 0 ]  = -0.0000017849850031;
      pd_scaling[ 1 ]  = -0.0000032596802369;
      pd_scaling[ 2 ]  =  0.0000312298758654;
      pd_scaling[ 3 ]  =  0.0000623390344610;
      pd_scaling[ 4 ]  = -0.0002599745524878;
      pd_scaling[ 5 ]  = -0.0005890207562444;
      pd_scaling[ 6 ]  =  0.0012665619292991;
      pd_scaling[ 7 ]  =  0.0037514361572790;
      pd_scaling[ 8 ]  = -0.0056582866866115;
      pd_scaling[ 9 ]  = -0.0152117315279485;
      pd_scaling[ 10 ] =  0.0250822618448678;
      pd_scaling[ 11 ] =  0.0393344271233433;
      pd_scaling[ 12 ] = -0.0962204420340021;
      pd_scaling[ 13 ] = -0.0666274742634348;
      pd_scaling[ 14 ] =  0.4343860564915321;
      pd_scaling[ 15 ] =  0.7822389309206135;
      pd_scaling[ 16 ] =  0.4153084070304910;
      pd_scaling[ 17 ] = -0.0560773133167630;
      pd_scaling[ 18 ] = -0.0812666996808907;
      pd_scaling[ 19 ] =  0.0266823001560570;
      pd_scaling[ 20 ] =  0.0160689439647787;
      pd_scaling[ 21 ] = -0.0073461663276432;
      pd_scaling[ 22 ] = -0.0016294920126020;
      pd_scaling[ 23 ] =  0.0008923136685824;
      break;

    case 30:  /* C(30) */

      pd_scaling[ 0 ]  = -0.0000000951765727;
      pd_scaling[ 1 ]  = -0.0000001674428858;
      pd_scaling[ 2 ]  =  0.0000020637618516;
      pd_scaling[ 3 ]  =  0.0000037346551755;
      pd_scaling[ 4 ]  = -0.0000213150268122;
      pd_scaling[ 5 ]  = -0.0000413404322768;
      pd_scaling[ 6 ]  =  0.0001405411497166;
      pd_scaling[ 7 ]  =  0.0003022595818445;
      pd_scaling[ 8 ]  = -0.0006381313431115;
      pd_scaling[ 9 ]  = -0.0016628637021860;
      pd_scaling[ 10 ] =  0.0024333732129107;
      pd_scaling[ 11 ] =  0.0067641854487565;
      pd_scaling[ 12 ] = -0.0091642311634348;
      pd_scaling[ 13 ] = -0.0197617789446276;
      pd_scaling[ 14 ] =  0.0326835742705106;
      pd_scaling[ 15 ] =  0.0412892087544753;
      pd_scaling[ 16 ] = -0.1055742087143175;
      pd_scaling[ 17 ] = -0.0620359639693546;
      pd_scaling[ 18 ] =  0.4379916262173834;
      pd_scaling[ 19 ] =  0.7742896037334738;
      pd_scaling[ 20 ] =  0.4215662067346898;
      pd_scaling[ 21 ] = -0.0520431631816557;
      pd_scaling[ 22 ] = -0.0919200105692549;
      pd_scaling[ 23 ] =  0.0281680289738655;
      pd_scaling[ 24 ] =  0.0234081567882734;
      pd_scaling[ 25 ] = -0.0101311175209033;
      pd_scaling[ 26 ] = -0.0041593587818186;
      pd_scaling[ 27 ] =  0.0021782363583355;
      pd_scaling[ 28 ] =  0.0003585896879330;
      pd_scaling[ 29 ] = -0.0002120808398259;
      break;

    default:

      LOCALDEF_UNSUPPORTED_FILTER_LENGTH;
      break;

    }

    break;

  default:
    MUTIL_ERROR( "Filter type is unsupported" );
    err = MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
    break;
  }

  /* check for errors */

  MEMLIST_FREE_ON_ERROR( err, &list );

  /* normalize the filters */

  if ( normalize ){
    for ( tap = 0; tap < filter_length; tap++ ){
      pd_scaling[ tap ] /= sqrt( 2.0 );
    }
  }

  /* Now that the scaling coefficients have been developed,
     calculate the wavelet coefficients through a quadrature
     mirror filter relation and store the results. This QMF
     relation is given by h_k = g_{L-1-k} * (-1)^{k-L} for
     k = 0, ..., L - 1 where L is the filter_length, h is
     the wavelet filter, and g is the scaling filter. The actual
     storage of the data is done on a flattened (row major) array */

  for ( tap = 0; tap < filter_length; tap++ ){

    pd_wavelet[ tap ] =
      pd_scaling[ filter_length - 1 - tap ] *
      (double) MUTIL_POW( -1, tap - filter_length );

    /* Check for interrupts */
    num_ops += (double) 2.0 * filter_length;

    if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
      MUTIL_FREE_WARN( memlist, &list );
      MUTIL_ERROR( "user interrupt" );
      return MUTIL_ERR_INTERRUPT;
    }
  }

  /* free nodes corresponding to registered
     wavelet and scaling filter memory, but
     do not free the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free the memory list */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( ":: Done with wavuniv_filters_daubechies()" );

  return MUTIL_ERR_OK;
}

/* Wavelet and scaling filter verification */
/* Documented in wav_filt.h                */
/* Written by William Constantine          */

mutil_errcode wavuniv_filters_daubechies_verify(
  const univ_mat  *wavelet_filter,
  const univ_mat  *scaling_filter,
  boolean          normalize )
{
  double        sum_wavelet_even_shift;
  double        sum_cross_even_shift;
  double        sum_scaling_even_shift;
  double        tol = (double) MUTIL_FLOAT_EPSILON;
  double        energy_scaling;
  double        energy_true;
  double        energy_wavelet;
  double        sum_scaling;
  double        sum_scaling_true;
  double        sum_wavelet;
  double       *g;
  double       *h;
  mutil_errcode err;
  sint32        L;
  sint32        l;
  sint32        n;
  sint32        skip;

  /* obtain the size of the wavelet filter */

  L = MATUNIV_NELEM( wavelet_filter );

  /* check the filter arguments */

  err = localfn_filter_matrices_check(
    L,
    wavelet_filter,
    scaling_filter);

  if ( err ) return err;

  /* set normalization factor */

  if ( normalize ) {
    sum_scaling_true = 1.0;
    energy_true      = 0.5;
  }
  else{
    sum_scaling_true = sqrt( 2.0 );
    energy_true      = 1.0;
  }

  /* set pointers */

  g = scaling_filter->mat.dblmat.data;
  h = wavelet_filter->mat.dblmat.data;

  /* test the wavelet filter */

  sum_wavelet            = 0.0;
  sum_scaling            = 0.0;
  energy_wavelet         = 0.0;
  energy_scaling         = 0.0;
  sum_wavelet_even_shift = 0.0;
  sum_scaling_even_shift = 0.0;
  sum_cross_even_shift   = 0.0;

  for ( n = -( L / 2 ) + 1; n < ( L / 2 ); n++ ){

    skip = 2 * MUTIL_ABS( n );

    for ( l = 0; l < ( L - skip ); l++ ){

      /* Form the inner product of wavelet filter
         and the even-shifted scaling filter to test
         for even shift orthogonality. In order
         to test both forward and reverse even shifts,
         the cross orthogonality components are
         calculated by shifting the wavelet filter
         for positive n, and by shifting the
         scaling filter for negative n. This need
         not be done when testing for even shift
         orthogonality within an individual filter
         due to symmetry. In the case where n = 0
         (no shift), the even-shift summation
         variables contain the filter energy. */

      if ( n == 0 ){

        /* summation of filter coefficients */

        sum_wavelet          += h[ l ];
        sum_scaling          += g[ l ];
        sum_cross_even_shift += g[ l ] * h[ l ];

        /* filter energy */

        energy_wavelet += h[ l ] * h[ l ];
        energy_scaling += g[ l ] * g[ l ];
      }
      else if ( n > 0 ){

        sum_wavelet_even_shift += h[ l ] * h[ l + skip ];
        sum_scaling_even_shift += g[ l ] * g[ l + skip ];
        sum_cross_even_shift   += g[ l ] * h[ l + skip ];
      }
      else if ( n < 0 ){
        sum_cross_even_shift += h[ l ] * g[ l + skip ];
      }

    } /* ends loop over l, the filter tap variable */

  } /* ends loop over n, the even-shift index */

  /* perform checks ... */

  /* ... sum over coefficients */

  if ( sum_wavelet > tol ){
    MUTIL_ERROR( "Mean of wavelet filter is not zero" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( MUTIL_ABS( MUTIL_ABS( sum_scaling ) - sum_scaling_true ) > tol ){
    MUTIL_ERROR( "Scaling filter coefficients do not "
      "sum to the proper value" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* ... energy */

  if ( MUTIL_ABS( energy_wavelet - energy_true ) > tol ){
    MUTIL_ERROR( "Wavelet filter does not have the correct energy" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }
  if ( MUTIL_ABS( energy_scaling - energy_true ) > tol ){
    MUTIL_ERROR( "Scaling filter does not have the correct energy" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* ... even shift orthogonality */

  if ( MUTIL_ABS( sum_wavelet_even_shift ) > tol ){
    MUTIL_ERROR( "Wavelet filter is not even shift orthogonal" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( MUTIL_ABS( sum_scaling_even_shift ) > tol ){
    MUTIL_ERROR( "Scaling filter is not even shift orthogonal" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  if ( MUTIL_ABS( sum_cross_even_shift ) > tol ){
    MUTIL_ERROR( "Wavelet and scaling filter are not "
      "even shift orthogonal to one another" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  return MUTIL_ERR_OK;
}

/* The squared gain functions for Daubechies  */
/* wavelet and scaling filters.               */
/* Documented in wav_filt.h                   */
/* Written by William Constantine             */

mutil_errcode wavuniv_filters_daubechies_gain(
  wav_filter_type   filter_type,
  sint32            filter_length,
  sint32            num_levels,
  sint32            num_fft,
  boolean           normalize,
  void             *intrp_ptr,
  univ_mat         *gain_frequency,
  univ_mat         *gain_wavelet,
  univ_mat         *gain_scaling)
{
  dcomplex          temp;           /* temporary variable used to store complex math operations        */
  dcomplex         *G1;             /* pointer to level 1 scaling filter gain function                 */
  dcomplex         *G;              /* pointer to current level's scaling filter gain function         */
  dcomplex         *Glast;          /* pointer to previous level's scaling filter gain function        */
  dcomplex         *H1;             /* pointer to level 1 wavelet filter gain function                 */
  dcomplex         *H;              /* pointer to current level's wavelet filter gain function         */
  double            dfreq;          /* the frequency interval between Fourier coefficients             */
  double            freq;           /* current Fourier frequency                                       */
  double            num_ops = 0.0;  /* used for user interrupts (number of operations count)           */
  mat_set           filters;        /* wavelet and scaling filter matrix set                           */
  mutil_errcode     err;            /* MUTIL error code                                                */
  sint32            f;              /* counting index for looping through Fourier frequencies          */
  sint32            index;          /* index for (circular) Fourier coefficients                       */
  sint32            j;              /* counting index                                                  */
  sint32            l;              /* counting index                                                  */
  double            magnification;  /* variable related to scale of current decomposition level        */
  sint32            num_pad = LOCALDEF_MAX( filter_length, num_fft); /* total length of zero padded filters */
  univ_mat          Gone;           /* first level scaling filter gain function                        */
  univ_mat          Hone;           /* first level wavelet filter gain function                        */
  univ_mat          GoneTranspose;
  univ_mat          HoneTranspose;
  univ_mat          g;              /* vector containing scaling filter (impulse response) (complex)   */
  univ_mat          h;              /* vector containing wavelet filter (impulse response) (complex)   */
  memlist           list;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_filters_daubechies_gain()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /* check input arguments */

  err = localfn_wavuniv_filters_daubechies_input_check(
    filter_length,
    filter_type );
  if ( err ) return err;

  /* allocate memory for output */

  err = matuniv_malloc_register( gain_frequency, 1,
    num_pad, MUTIL_DOUBLE, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( gain_wavelet, num_levels,
    num_pad, MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( gain_scaling, num_levels,
    num_pad, MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* malloc space for matrices */

  err = matuniv_malloc_register( &g, num_pad, 1,
    MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &h, num_pad, 1,
    MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &Gone, num_pad, 1,
    MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &Hone, num_pad, 1,
    MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &GoneTranspose, 1, num_pad,
    MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_malloc_register( &HoneTranspose, 1, num_pad,
    MUTIL_DCOMPLEX, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* obtain filters and register with memory manager */

  err = wavuniv_filters_daubechies(
    filter_length,
    filter_type,
    normalize,
    intrp_ptr,
    &filters );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_register( &list, &filters, MEMTYPE_MATSET );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* explicitly fill the zero padded wavelet and scaling filter vectors.
     this operation zeros out the imaginary portion, zero pads the real
     portion beyond the filter length, and transposes the data in one swoop. */

  for ( l = 0; l < num_pad; l++ ){

    if ( l < filter_length ){
      h.mat.cpxmat.data[ l ].re = filters.mats[ 0 ].mat.dblmat.data[ l ];
      g.mat.cpxmat.data[ l ].re = filters.mats[ 1 ].mat.dblmat.data[ l ];
    }
    else{
      g.mat.cpxmat.data[ l ].re = 0;
      h.mat.cpxmat.data[ l ].re = 0;
    }

    g.mat.cpxmat.data[ l ].im = 0;
    h.mat.cpxmat.data[ l ].im = 0;
  }

  /* perform FFT */

  err = siguniv_transform_discrete_fourier( &g, FALSE, intrp_ptr, &Gone );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = siguniv_transform_discrete_fourier( &h, FALSE, intrp_ptr, &Hone );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* transpose the first level response vectors. because they are
     1D (row) vectors, the output of the transpose function can be written
     to the same data location as the input */

  err = matuniv_transpose( &Gone, intrp_ptr, &GoneTranspose );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_transpose( &Hone, intrp_ptr, &HoneTranspose );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* copy the first level gain functions into the appropriate
    rows in the result matrices   */

  err = matuniv_assign_submat( &GoneTranspose, 0, 0, intrp_ptr, gain_scaling );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matuniv_assign_submat( &HoneTranspose, 0, 0, intrp_ptr, gain_wavelet );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* create pointers */

  G1 = &( Gone.mat.cpxmat.data[ 0 ] );
  H1 = &( Hone.mat.cpxmat.data[ 0 ] );

  /* now create the remaining gain functions by
     iterating through each level */

  dfreq = 1.0 / ( double ) num_pad;

  MUTIL_TRACE( "Starting loop over each level." );

  for( j = 1; j < num_levels; j++ ){

    /* set pointers */

    Glast = &( gain_scaling->mat.cpxmat.data[ ( j - 1 ) * num_pad ]);
    G     = &( gain_scaling->mat.cpxmat.data[ j * num_pad ]);
    H     = &( gain_wavelet->mat.cpxmat.data[ j * num_pad ]);

    magnification = pow( 2.0, (double) j );

    for ( f = 0; f < num_pad; f++ ){

      freq  = f * dfreq;

      /* the math on the right will always be an integer result so flooring
	 it will cause no harm and will attenuate difficulties with LINUX casts */
      index =
        (sint32) floor( LOCALDEF_MODULO( freq * magnification , 1.0 ) / dfreq );

      /* fill the frequency vector */

      if ( j == 1 ){
        gain_frequency->mat.dblmat.data[ f ] = freq;
      }

      /* G_j(f) = G(2^{j-1} * f) * G_{j-1}(f) */

      MUTIL_CPX_MULT( G1[ index ], Glast[ f ], temp );

      G[ f ] = temp;

      /* H_j(f) = H( 2^{j-1} * f ) * G_{j-1}(f) */

      MUTIL_CPX_MULT( H1[ index ], Glast[ f ], temp );

      H[ f ] = temp;
    }

    /* Check for interrupt */

    num_ops += 3.0 * num_pad;
    if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      MUTIL_FREE_WARN( memlist, &list );
      return MUTIL_ERR_INTERRUPT;
    }

  }

  /* wrap things up */

  MUTIL_TRACE( "Done with loop over each level." );

  /* free nodes corresponding to registered
     output matrices memory, but
     do not free the memory itself */

  err = memlist_member_unregister( gain_frequency, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = memlist_member_unregister( gain_wavelet, &list );
  if ( err ){
    MUTIL_FREE_WARN( matuniv, gain_frequency );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  err = memlist_member_unregister( gain_scaling, &list );
  if ( err ){
    MUTIL_FREE_WARN( matuniv, gain_frequency );
    MUTIL_FREE_WARN( matuniv, gain_scaling );
    MEMLIST_FREE_ON_ERROR( err, &list );
  }

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_filters_daubechies_gain()" );

  return MUTIL_ERR_OK;
}

/******************************************************/
/* FILTER INDEXING FUNCTION DEFINITIONS               */
/******************************************************/


/* Boundary and interior wavelet coefficient */
/* identification for the DWT and MODWT.     */
/* Function documented in wav_filt.h         */
/* Written by William Constantine            */

mutil_errcode wavuniv_transform_coefficient_boundaries(
  sint32         n_level,
  sint32         filter_length,
  sint32         n_sample,
  wav_transform  transform_type,
  void          *intrp_ptr,
  mat_set       *result )
{
  double        num_ops = 0.0;
  memlist       list;
  mutil_errcode err;
  sint32        all_length_j;
  sint32        boundary_length_j;
  sint32        filter_width;
  sint32        interior_high_j;
  sint32        interior_length_j;
  sint32        interior_low_j;
  sint32        level;
  sint32        j;
  sint32        n_coeff;
  sint32        ndim = 1;
  sint32        dims = 5;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_transform_coefficient_boundaries()" );

  /* initialize memory manager */

  MEMLIST_INIT( list );

  /* check input arguments */

  err = localfn_wavuniv_transform_coefficient_boundaries_inputs_check(
    n_level,
    filter_length,
    n_sample,
    transform_type );
  if ( err ) return err;

  /* allocate memory for output */

  err = matset_malloc_register( result, ndim, &dims, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = matset_malloc_matrices( result, 1, n_level, MUTIL_SINT32 );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* initialize variables */

  n_coeff = n_sample;

  /* calculate boundary and interior wavelet
     coefficient indices */

  for( level = 1; level <= n_level; level++ ){

    switch( transform_type ) {

    case WAV_TRANSFORM_DWT:

      filter_width =
        (sint32) ceil( ( filter_length - 2 ) *
        ( 1 - MUTIL_POW( 2, - level ) ) );
      n_coeff =
        (sint32) floor( n_sample / MUTIL_POW( 2, level ) );
      break;

    case WAV_TRANSFORM_MODWT:

      filter_width =
        /* (sint32) floor( ( MUTIL_POW( 2, level ) - 1 ) *
        ( filter_length - 1 ) ); */
        ( (1 << level) - 1 ) * ( filter_length - 1 );
      break;

    default:

      filter_width = 0;
      MUTIL_ERROR( "This transform type is unsupported" );
      err = MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
      break;
    }

    /* check for errors */

    MEMLIST_FREE_ON_ERROR( err, &list );

    /* calculate results */

    boundary_length_j = MUTIL_MIN( filter_width, n_coeff );
    interior_low_j    = MUTIL_MIN( n_coeff, boundary_length_j + 1 );
    interior_high_j   = n_coeff;
    interior_length_j = n_coeff - boundary_length_j;
    all_length_j      = interior_length_j + boundary_length_j;

    /* assign a zero to the interior low and high indices in the case
    where the interior length is zero */

    if ( interior_length_j == 0 ){
      interior_low_j = interior_high_j = 0;
    }

    /* assign results ( matset index: description ):

       0: interior_low,
       1: interior_high,
       2: interior_length
       3: boundary_length
       4: all_length                              */

    j = level - 1;

    result->mats[ 0 ].mat.s32mat.data[ j ] = interior_low_j;
    result->mats[ 1 ].mat.s32mat.data[ j ] = interior_high_j;
    result->mats[ 2 ].mat.s32mat.data[ j ] = interior_length_j;
    result->mats[ 3 ].mat.s32mat.data[ j ] = boundary_length_j;
    result->mats[ 4 ].mat.s32mat.data[ j ] = all_length_j;

    /* Check for interrupts */
    num_ops += ( double ) 10.0 * n_level;

    if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
      MUTIL_ERROR( "user interrupt" );
      return MUTIL_ERR_INTERRUPT;
    }
  }

  /* check for errors */

  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free nodes corresponding to registered
     output memory, but do not free the
     memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free registered local memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_transform_coefficient_boundaries()" );

  return MUTIL_ERR_OK;
}

/******************************************************/
/* FILTER PHASE FUNCTION DEFINITIONS                  */
/******************************************************/

/* Zero phase shift factors for Daubechies */
/* symmlet and Coiflet filters.            */
/* Function documented in wav_filt.h       */
/* Written by William Constantine          */

mutil_errcode wavuniv_filters_zero_phase(
  wav_filter_type  filter_type,
  sint32           filter_length,
  sint32           n_level,
  void            *intrp_ptr,
  mat_set         *result )
{
  /* boolean         is_symmlet; */
  double          num_ops = 0.0;
  memlist         list;
  mutil_errcode   err;
  sint32          Lj;
  sint32          Sjn;
  sint32          dims;
  sint32          dwt_shift;
  sint32          j;
  sint32          n;
  sint32          n_nodes;
  sint32          shift;
  sint32          sym_case;
  sint32          tauj;
  sint32         *ps_dwpt;
  sint32         *ps_dwt_scaling;
  sint32         *ps_dwt_wavelet;
  sint32         *ps_modwpt;
  sint32         *ps_modwt_scaling;
  sint32         *ps_modwt_wavelet;
  sint32_mat      ncol;
  sint32_mat      nrow;

  MUTIL_INTERRUPT_INIT( intrp_ptr );

  MUTIL_TRACE( "Start wavuniv_filters_zero_phase()" );

  /* initialize memory list */

  MEMLIST_INIT( list );

  /*** filter_length argument (general check) ... */

  /* ... if filter_length is positive */

  if ( filter_length <= 1 ){
    MUTIL_ERROR( "Filter length must be greater than one." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if filter_length is even */

  if ( ( filter_length % 2 ) != 0 ){
    MUTIL_ERROR( "Filter length must be even." );

    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*** check filter_type and specific lengths */

  switch ( filter_type ) {
    case WAV_FILTER_LEAST_ASYMMETRIC:

      /* is_symmlet = TRUE; */

      if ( ( filter_length < 8 ) || ( filter_length > 20 ) ){
        MUTIL_ERROR( "Unsupported filter length" );
        return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
      }

      if ( ( ( filter_length / 2 ) % 2 ) == 0 ){
        sym_case = 0;
      }
      else if ( ( filter_length == 10 ) || ( filter_length == 18 ) ){
        sym_case = 1;
      }
      else if ( filter_length == 14 ){
        sym_case = 2;
      }
      else{
        MUTIL_ERROR( "This filter length is unsupported" );
        return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
      }

      break;

    case WAV_FILTER_COIFLET:

      /* is_symmlet = FALSE; */

      if ( ( filter_length != 6 ) &&
        ( filter_length != 12 ) &&
        ( filter_length != 18 ) &&
        ( filter_length != 24 ) &&
        ( filter_length != 30 ) ){

        MUTIL_ERROR( "Unsupported filter length" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }

      break;

    case WAV_FILTER_EXTREMAL_PHASE:

      if ( ( filter_length != 4 ) &&
        ( filter_length != 6 )){

        MUTIL_ERROR( "Unsupported filter length" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }

      break;

    case WAV_FILTER_HAAR:

      if (filter_length != 2){

        MUTIL_ERROR( "Unsupported filter length" );
        return MUTIL_ERR_ILLEGAL_VALUE;
      }

      break;

    default:
      MUTIL_ERROR( "Unsupported filter type" );
      return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  /*** check the n_level argument */

  if ( n_level <= 0 ){
    MUTIL_ERROR( "Number of decomposition levels must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  /* create dimension vectors for matset matrices */

  /* n_nodes = (sint32) MUTIL_POW( 2, n_level + 1 ) - 2; */
  n_nodes = ( 2 << n_level ) - 2;

  err = mats32_malloc_register( &nrow, 1, 4, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  err = mats32_malloc_register( &ncol, 1, 4, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  for ( j = 0; j < 4; j++ ){

    if ( j < 2 ){
      ncol.data[ j ] = n_level * 2;
    }
    else{
      ncol.data[ j ] = n_nodes;
    }

    nrow.data[ j ] = 1;
  }

  /* allocate space for output */

  dims = 4;
  err = matset_malloc_register( result, (sint32) 1, &dims, &list );
  if ( err ) return err;

  err = matset_malloc_matrices_arbitrary_size(
    result, &nrow, &ncol, MUTIL_SINT32 );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* initialize pointers */

  ps_dwt_wavelet   = result->mats[0].mat.s32mat.data;
  ps_dwt_scaling   = ps_dwt_wavelet + n_level;
  ps_modwt_wavelet = result->mats[1].mat.s32mat.data;
  ps_modwt_scaling = ps_modwt_wavelet + n_level;
  ps_dwpt          = result->mats[2].mat.s32mat.data;
  ps_modwpt        = result->mats[3].mat.s32mat.data;

  for( j = 1; j <= n_level; j++){

    /* initialize level dependent variables */

    /* Lj = - (sint32) ( ( floor( MUTIL_POW( 2, j ) ) - 1 ) *
      ( filter_length - 1 )  + 1 ) / 2; */
    Lj = - ( ( ( 1 << j ) - 1 ) * ( filter_length - 1 ) + 1 ) / 2;

    /* tauj = (sint32) MUTIL_POW( 2, j - 1 ); */
    tauj = 1 << ( j - 1 );

    /* calculate shifts */

    for ( n = 0; n < MUTIL_POW( 2, j ); n++ ){

      Sjn = localfn_S_value( j, n );

      if ( filter_type == WAV_FILTER_LEAST_ASYMMETRIC ){

        switch( sym_case ){
          case 0:
            shift = Lj + tauj - Sjn;
            break;
          case 1:
            shift = Lj - ( tauj - Sjn ) + 1;
            break;
          case 2:
            shift = Lj + 3 * ( tauj - Sjn ) - 1;
            break;
          default:
            MUTIL_ERROR( "Symmlet case out of range." );
            return MUTIL_ERR_ILLEGAL_VALUE;
        }
      }
      else if ( filter_type == WAV_FILTER_COIFLET ){
        shift = Lj - ( filter_length - 3 ) / 3 *
          ( tauj - Sjn -1/2 ) + 1/2;
      }
      else{ /* Extremal phase or Haar */
        shift = Lj + tauj - Sjn;
      }

      dwt_shift = DWT_SHIFT( shift, j );

      /* assign shift result to matrices */

      *ps_modwpt = shift;
      ps_modwpt++;

      *ps_dwpt = dwt_shift;
      ps_dwpt++;

      if ( n == 0 ){

        *ps_modwt_scaling = shift;
        ps_modwt_scaling++;

        *ps_dwt_scaling = dwt_shift;
        ps_dwt_scaling++;
      }
      else if ( n == 1 ){

        *ps_modwt_wavelet = shift;
        ps_modwt_wavelet++;

        *ps_dwt_wavelet = dwt_shift;
        ps_dwt_wavelet++;
      }

      /* Check for interrupts */
      num_ops += ( double ) 10.0 * n_level;

      if ( MUTIL_INTERRUPT( num_ops, intrp_ptr ) ) {
        MUTIL_ERROR( "user interrupt" );
        MUTIL_FREE_WARN( memlist, &list );
        return MUTIL_ERR_INTERRUPT;
      }
    } /* end loop over local node index */
  } /* end loop over decomposition level index */

  /* free node corresponding to registered
     output memory, but do not free the memory itself */

  err = memlist_member_unregister( result, &list );
  MEMLIST_FREE_ON_ERROR( err, &list );

  /* free all other memory */

  MUTIL_FREE_WARN( memlist, &list );

  MUTIL_TRACE( "Done with wavuniv_filters_zero_phase()" );

  return MUTIL_ERR_OK;
}


/******************************************************/
/* STATIC FILTER GENERATION FUNCTION DEFINITIONS      */
/******************************************************/


static mutil_errcode localfn_filter_fill(
  sint32  filter_length,
  double *in,
  double *out )
{
  sint8 i;

  for ( i = 0; i < filter_length;  i++ ){
    out[ i ] = in[ i ];
  }
  return MUTIL_ERR_OK;
}

static mutil_errcode localfn_wavuniv_filters_daubechies_input_check(
  sint32           filter_length,
  wav_filter_type  filter_type )
{
  /*** filter_length argument */

  /* ... positive length */

  if ( filter_length <= 1 ){
    MUTIL_ERROR( "Filter length must be greater than one" );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /* ... even filter length */

  if ( ( filter_length % 2 ) != 0 ){
    MUTIL_ERROR( "Filter length must be even." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*** filter_type argument */

  switch( filter_type ){
  case WAV_FILTER_EXTREMAL_PHASE:
  case WAV_FILTER_LEAST_ASYMMETRIC:
  case WAV_FILTER_BEST_LOCALIZED:
  case WAV_FILTER_COIFLET:
  case WAV_FILTER_HAAR:
    break;
  default:
    MUTIL_ERROR( "Filter type is unsupported" );
    return MUTIL_ERR_FEATURE_NOT_IMPLEMENTED;
  }

  return MUTIL_ERR_OK;
}

static mutil_errcode localfn_filter_matrices_check(
  sint32    filter_length,
  const univ_mat *wavelet_filter,
  const univ_mat *scaling_filter )
{

  mutil_errcode err;

  /*** wavelet filter */

  if ( wavelet_filter == (univ_mat *) NULL ){
    MUTIL_ERROR( "Pointer to wavelet filter matrix is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  err = matuniv_validate( wavelet_filter );
  if ( err ) return err;

  if ( wavelet_filter->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "Wavelet filter matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  if ( !MATANY_IS_VEC( &(wavelet_filter->mat.dblmat) ) ){
    MUTIL_ERROR( "Wavelet filter matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if ( MATUNIV_NELEM( wavelet_filter ) != filter_length ){
    MUTIL_ERROR( "Number of elements in wavelet_filter matrix "
      "must equal the specified filter length." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /*** scaling filter */

  if ( scaling_filter == (univ_mat *) NULL ){
    MUTIL_ERROR( "Pointer to scaling filter matrix is NULL" );
    return MUTIL_ERR_NULL_POINTER;
  }

  err = matuniv_validate( scaling_filter );
  if ( err ) return err;

  if ( scaling_filter->type != MUTIL_DOUBLE  ){
    MUTIL_ERROR( "scaling filter matrix must be of type DOUBLE." );
    return MUTIL_ERR_ILLEGAL_TYPE;
  }
  if ( !MATANY_IS_VEC( &(scaling_filter->mat.dblmat) ) ){
    MUTIL_ERROR( "Scaling filter matrix must have a single column or row." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }
  if ( MATUNIV_NELEM( scaling_filter ) != filter_length ){
    MUTIL_ERROR( "Number of elements in scaling filter matrix must "
      "equal the specified filter length." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  return MUTIL_ERR_OK;
}

/******************************************************/
/* STATIC FILTER INDEXING FUNCTION DEFINITIONS        */
/******************************************************/

/** Validate the SDF input for wavelet variance functions.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_var.c
 * @usage #err = localfn_wavuniv_transform_coefficient_boundaries_inputs_check( n_level, filter_length, n_sample, transform_type );#
 * @return                  Standard mutils error/OK code.
 * @param  num_levels       The number of decomposition levels.
 * @param  filter_length    The length of the wavelet filter.
 * @param  num_points       The number of points in the time series.
 * @param  transform_type   Specifies the type of wavelet transform.
 *
 * @see wavuniv_transform_coefficient_boundaries
 * @private
 */
static mutil_errcode localfn_wavuniv_transform_coefficient_boundaries_inputs_check(
  sint32         n_level,
  sint32         filter_length,
  sint32         n_sample,
  wav_transform  transform_type )
{
  MUTIL_TRACE( "Start localfn_wavuniv_coefficient_"
               "boundaries_inputs_check()" );

  /*** check the n_level argument ... ***/

  if ( n_level <= 0 ){
    MUTIL_ERROR( "Number of decomposition levels must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  /*** filter_length argument ... */

  /* ... if filter_length is positive */

  if ( filter_length <= 1 ){
    MUTIL_ERROR( "Filter length must be greater than one." );
    return MUTIL_ERR_ILLEGAL_SIZE;
  }

  /* ... if filter_length is even */

  if ( ( filter_length % 2 ) != 0 ){
    MUTIL_ERROR( "Filter length must be even." );
    return MUTIL_ERR_ILLEGAL_VALUE;
  }

  /*** check the n_sample argument ... ***/

  if ( n_sample <= 0 ){
    MUTIL_ERROR( "Number of samples in the time series must be positive." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  /*** check the transform type ... ***/

  if ( ( transform_type != (wav_transform) WAV_TRANSFORM_DWT ) &&
       ( transform_type != (wav_transform) WAV_TRANSFORM_MODWT ) ){
    MUTIL_ERROR( "Transform type must be either "
                 "WAV_TRANSFORM_DWT or WAV_TRANSFORM_MODWT." );
    return MUTIL_ERR_ILLEGAL_VALUE ;
  }

  MUTIL_TRACE( "Done with localfn_wavuniv_coefficient_"
               "boundaries_inputs_check()" );

  return MUTIL_ERR_OK;
}


/******************************************************/
/* STATIC FILTER PHASE FUNCTION DEFINITIONS           */
/******************************************************/

/** Calculate the S(j,n,1) index factor for one crystal
 * in a wavelet packet transform. The S(n,1) is
 * used to calculate the zero phase shift for Daubechies
 * (approximately) linear phase filters. The crystal
 * is specified via the natural order (level,node).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @source wav\_zero.c
 * @library wavelets
 * @usage #S = localfn_S_value( level, node );#
 * @return The S(n,1) value.
 * @param  level  The wavelet packet decomposition level.
 * @param  node  The
 * @see wavuniv_filters_zero_phase
 * @private
 */
static sint32 localfn_S_value( sint32 level, sint32 node )
{
  sint32 j;
  sint32 n = node;
  sint32 n_mod;
  sint32 sum = 0;

  for ( j = level; j >= 1; j-- ){

    n_mod = n % 4;

    if ( ( n_mod == 1 ) || ( n_mod == 2 ) ){

      /* sum += (sint32) MUTIL_POW( 2, j - 1 ); */
      sum += ( 1 << ( j - 1 ) );
    }

    n = (sint32) floor( (double) n / 2.0 );
  }

  return sum;
}
