
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_plat.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_PLAT_H_
#define IN_UT_PLAT_H_


/*
   This file contains platform-specific defines for the mutils
   library.
*/

#if defined(_IBMR2) && defined(_AIX) /*(*/
#define MUTIL_IBMRS6000
#define __unix
#endif

#if defined(CYGWIN_BUILD)
#undef __unix
#endif

#if defined(__unix)      /*(*/

#if defined(__linux__)     /*(*/
#if __BYTE_ORDER == __LITTLE_ENDIAN /*(*/
#define _LITTLE_ENDIAN
#else                        /*)(*/
#define _BIG_ENDIAN
#endif                       /*)*/

#elif defined(MUTIL_IBMRS6000) /*(*/
#include <sys/machine.h>
#if BYTE_ORDER == BIG_ENDIAN /*(*/
#define _BIG_ENDIAN
#else                        /*)(*/
#define _LITTLE_ENDIAN
#endif                       /*)*/

#elif defined(__hpux)      /*(*/
#include <machine/param.h>

#elif defined( _MIPSEB )   /*(*/
#define _BIG_ENDIAN
#else                      /*)(*/
#include <sys/isa_defs.h>
#endif                     /*)*/
/** Macro prefix for exported symbols.
 * This macro is used as a prefix for exported symbols in the main
 * (non-wrapper) mutils library.
 * The macro evaluates to the extern keyword on Unix platforms,
 * and to the \_declspec(dllexport) and \_declspec(dllimport)
 * keywords for the Windows platform.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_plat.h
 * @source ut\_plat.h
 * @library
 * @see MUTIL_WRAPEXPORT
 * @see Compile-Time Options
 */
#define MUTIL_LIBEXPORT extern


/** Macro prefix for wrapper library exported symbols.
 * This macro is used as a prefix for exported symbols in the wrapper
 * mutils library.
 * The macro evaluates to the extern keyword on Unix platforms,
 * and to the \_declspec(dllexport) and \_declspec(dllimport)
 * keywords for the Windows platform.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_plat.h
 * @source ut\_plat.h
 * @library
 * @see MUTIL_LIBEXPORT
 * @see Compile-Time Options
 */
#define MUTIL_WRAPEXPORT extern


/** Macro prefix for export of private symbols.
 * This macro is used as a prefix for exported private symbols in the main
 * (non-wrapper) mutils library for testing purposes.
 * The macro evaluates to the extern keyword on Unix platforms,
 * and to the \_declspec(dllexport) and \_declspec(dllimport)
 * keywords for the Windows platform.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_plat.h
 * @source ut\_plat.h
 * @library
 * @see MUTIL_LIBEXPORT
 * @see Compile-Time Options
 */
#ifdef DEBUG
#define PRIV_LIBEXPORT extern
#else
#define PRIV_LIBEXPORT
#endif


#else /* ifdef __unix */
      /* NOTE: we assume here that everything that is not Unix is
         Windows */

#ifdef MUTIL_STATIC

#define MUTIL_LIBEXPORT  extern
#define MUTIL_WRAPEXPORT extern

#ifdef DEBUG
#define PRIV_LIBEXPORT extern
#else
#define PRIV_LIBEXPORT
#endif

#else

#ifdef MUTIL_WIN_DLLCREATING
#define MUTIL_LIBEXPORT __declspec(dllexport)
#else
#define MUTIL_LIBEXPORT __declspec(dllimport)
#endif /* ifdef MUTIL_WIN_DLLCREATING */


#ifdef MUTIL_WIN_WRAPCREATING
#define MUTIL_WRAPEXPORT __declspec(dllexport)
#else
#define MUTIL_WRAPEXPORT __declspec(dllimport)
#endif /* ifdef MUTIL_WIN_WRAPCREATING */

#ifdef DEBUG
#define PRIV_LIBEXPORT MUTIL_LIBEXPORT
#else
#define PRIV_LIBEXPORT
#endif

#endif /* ifdef MUTIL_STATIC */

#endif /* if/else on defined __unix */


/* Determine whether system is little-endian or big-endian */
/* Use the system/compiler built-in defines if possible for
   highest reliability. */

#if defined(_BIG_ENDIAN)


/** Big-endian-ness of the platform.
 * This value is 1 on platforms that have a big-endian architecture
 * (where the most significant byte comes first), and 0 on
 * platforms that have a little-endian architecture (where the
 * least significant byte comes first).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_plat.h
 * @source ut\_plat.h
 * @library
 * @see MUTIL_SYSTEM_IS_LITTLEENDIAN
 */
#define MUTIL_SYSTEM_IS_BIGENDIAN 1


/** Little-endian-ness of the platform.
 * This value is 0 on platforms that have a big-endian architecture
 * (where the most significant byte comes first), and 1 on
 * platforms that have a little-endian architecture (where the
 * least significant byte comes first).
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_plat.h
 * @source ut\_plat.h
 * @library
 * @see MUTIL_SYSTEM_IS_BIGENDIAN
 */
#define MUTIL_SYSTEM_IS_LITTLEENDIAN 0

#elif defined(_LITTLE_ENDIAN)

/* see documentation above */
#define MUTIL_SYSTEM_IS_BIGENDIAN 0
#define MUTIL_SYSTEM_IS_LITTLEENDIAN 1

#elif defined(_WINDOWS)

/* see documentation above */
#define MUTIL_SYSTEM_IS_BIGENDIAN 0
#define MUTIL_SYSTEM_IS_LITTLEENDIAN 1

/* If MUTIL_SYSTEM_IS_BIGENDIAN is not getting defined, and
   the compiler complains, add
   machine-specific code here as further #elif directives. */

#endif /* ifdef to determine big-endian-ness */

#endif /* IN_UT_PLAT_H_*/
