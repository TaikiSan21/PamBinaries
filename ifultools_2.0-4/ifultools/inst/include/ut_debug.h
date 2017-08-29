
/* @(#) $File: //depot/Research/ifultools/pkg/ifultools/inst/include/ut_debug.h $: $Revision: #1 $, $Date: 2008/03/21 $  */
/* This is a self-documenting doc++ file. */

#ifndef IN_UT_DEBUG_H
#define IN_UT_DEBUG_H

#include "ut_plat.h"
#include "ut_err.h"
#include "ut_type.h"

/*
   This file contains macros, typedefs, and functions for debugging code
   in the mutils library.  There are macros for three levels of
   debug message printing, as well as a debug assertion function.
*/


#ifdef __cplusplus
extern "C" {
#endif


/** Enum of debug message levels for mutils library.
 * The mutils library supports three levels of severity of debugging messages.
 * Trace-level informational messages are used to trace function
 * calls, parameter values, and progress.  Warning messages are used
 * to indicate possible problems that can be handled by the functions.
 * Error messages are used to indicate fatal problems.
 *
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_debug.h
 * @source ut\_debug.h
 * @library wrap
 * @same
 *    #typedef enum _mutil_msg_level mutil_msg_level;#
 * @see mutil_msg_print
 * @see _mutil_errcode
 */
enum _mutil_msg_level
{
  /** Trace-level informational messages */
  MUTIL_MSG_LEVEL_TRACE,

  /** Warning-level messages */
  MUTIL_MSG_LEVEL_WARN,

  /** Error-level messages */
  MUTIL_MSG_LEVEL_ERROR,

  /** No messages */
  MUTIL_MSG_LEVEL_NONE
};


/* See the documentation for _mutil_msg_level (above) for a
 * description of the enum. */
typedef enum _mutil_msg_level mutil_msg_level;


/** Display a debug message.
 * The message includes the file name, and is only
 * displayed if the given level of message display is turned on.
 * If line number printing is turned on, the line number will also
 * be printed.  If time printing is turned on, the time is printed
 * as well.  Normally it is called only via the MUTIL\_TRACE, MUTIL\_WARN,
 * and MUTIL\_ERROR macros, which supply the severity, file name, and line
 * number automatically.
 *
 * All debug message displays in the mutils library go through this
 * function.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_debug.h
 * @source ut\_debug.c
 * @library wrap
 * @usage #err_code = mutil_msg_print(MUTIL_MSG_LEVEL_TRACE, "Informational message", "file.c", 32);#
 * @return    Standard mutils error/OK code.
 * @param  level    The severity level of the message.
 * @param  message  The message to be printed.
 * @param  filename The file name where the message was generated.
 * @param  linenum  The line number where the message was generated.
 * @see MUTIL_TRACE
 * @see MUTIL_WARN
 * @see MUTIL_ERROR
 * @see _mutil_msg_level
 * @see _mutil_errcode
 * @see mutil_msg_set_level
 * @see mutil_msg_print_line
 * @see mutil_msg_print_time
 */
MUTIL_WRAPEXPORT mutil_errcode mutil_msg_print( mutil_msg_level level,
  const char *message, const char *filename, int linenum);


/** Set the debug message level to be displayed.
 * Choose the lowest severity of debug messages to display: trace,
 * warning, error, or none.  Messages of severity equal to or greater
 * than the chosen severity are displayed.  If the severity level
 * chosen is none, or if the debug version of the library is not
 * used, no messages will be displayed. The default message level
 * is warning.
 *
 * In multi-threaded environments,
 * there is only one global debug message level value, but this
 * function sets it in a thread-safe manner.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_debug.h
 * @source ut\_debug.c
 * @library wrap
 * @usage #mutil_msg_set_level(MUTIL_MSG_LEVEL_TRACE);#
 * @param  level   The severity level of messages to be displayed.
 * @see MUTIL_TRACE
 * @see MUTIL_WARN
 * @see MUTIL_ERROR
 * @see _mutil_msg_level
 * @see mutil_msg_print
 */
MUTIL_WRAPEXPORT void mutil_msg_set_level( mutil_msg_level level );


/** Turn line-number printing on or off.
 * Choose whether or not line numbers are printed in debug messages.
 * Note that  if the DEBUG compiler directive is not used,
 * the debug message macros are disabled and no messages will be
 * generated unless the mutil\_msg\_print function is called directly.
 *
 * In multi-threaded environments,
 * there is only one global debug line-number printing flag, but this
 * function sets it in a thread-safe manner.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_debug.h
 * @source ut\_debug.c
 * @library wrap
 * @usage #mutil_msg_print_line(TRUE);#
 * @param   show_nums   TRUE to print numbers; FALSE to omit them.
 * @see mutil_msg_print
 */
MUTIL_WRAPEXPORT void mutil_msg_print_line( boolean show_nums );


/** Turn time printing on or off.
 * Choose whether or not the time is printed in debug
 * messages.  Note that  if the DEBUG compiler directive is not used,
 * the debug message macros are disabled and no messages will be
 * generated unless the mutil\_msg\_print function is called directly.
 *
 * In multi-threaded environments,
 * there is only one global debug time printing flag, but this
 * function sets it in a thread-safe manner.
 *
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_debug.h
 * @source ut\_debug.c
 * @library wrap
 * @usage #mutil_msg_print_time(TRUE);#
 * @param   show_time   TRUE to print time; FALSE to omit it.
 * @see mutil_msg_print
 */
MUTIL_WRAPEXPORT void mutil_msg_print_time( boolean show_time );


/** Abort program on fatal error.
 * This function is in the wrap\_std library, so that applications
 * needing their own implementation can redefine it.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_debug.h
 * @source ut\_debug.c
 * @library wrap
 * @usage #mutil_abort();#
 */
MUTIL_WRAPEXPORT void mutil_abort();


#ifdef DEBUG

/** Display a trace-level informational message.
 * Calls \Ref{mutil_msg_print} with message level set to trace
 * to display an informational message.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_debug.h; if the DEBUG compiler directive is not used,
 * then it is defined to be a null operation.
 * @source ut\_debug.h
 * @library wrap
 * @usage #MUTIL_TRACE("Informational message");#
 * @param   MSG    Constant character string giving the message.
 * @see MUTIL_WARN
 * @see MUTIL_ERROR
 * @see _mutil_msg_level
 * @see _mutil_errcode
 * @see mutil_msg_print
 * @see mutil_msg_set_level
 */
#define MUTIL_TRACE(MSG) \
     (void) mutil_msg_print(MUTIL_MSG_LEVEL_TRACE, (MSG), __FILE__,  __LINE__)


/** Display a warning-level informational message.
 * Calls \Ref{mutil_msg_print} with message level set to warning
 * to display a warning message.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_debug.h; if the DEBUG compiler directive is not used,
 * then it is defined to be a null operation.
 * @source ut\_debug.h
 * @library wrap
 * @usage #MUTIL_WARN("Warning message");#
 * @param   MSG    Constant character string giving the message.
 * @see MUTIL_TRACE
 * @see MUTIL_ERROR
 * @see _mutil_msg_level
 * @see _mutil_errcode
 * @see mutil_msg_print
 * @see mutil_msg_set_level
 */
#define MUTIL_WARN(MSG) \
     (void) mutil_msg_print(MUTIL_MSG_LEVEL_WARN, (MSG), __FILE__,  __LINE__)


/** Display an error-level informational message.
 * Calls \Ref{mutil_msg_print} with message level set to error
 * to display an error message.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_debug.h; if the DEBUG compiler directive is not used,
 * then it is defined to be a null operation.
 * @source ut\_debug.h
 * @library wrap
 * @usage #MUTIL_ERROR("Error message");#
 * @param   MSG    Constant character string giving the message.
 * @see MUTIL_TRACE
 * @see MUTIL_WARN
 * @see _mutil_msg_level
 * @see _mutil_errcode
 * @see mutil_msg_print
 * @see mutil_msg_set_level
 */
#define MUTIL_ERROR(MSG) \
     (void) mutil_msg_print(MUTIL_MSG_LEVEL_ERROR, (MSG), __FILE__,  __LINE__)


/** Debug assertion.
 * This macro checks to see if the asserted condition is true
 * or false.  If it is false, it then displays an
 * error-level message by calling \Ref{MUTIL_ERROR}, containing
 * the asserted condition, and then calls the \Ref{mutil_abort} function.
 *
 * @author Copyright (c), 1988, 2006 Insightful Corp.  All rights reserved.
 * @include ut\_debug.h; if the DEBUG compiler directive is not used,
 * then it is defined to be a null operation.
 * @source ut\_debug.h
 * @library wrap
 * @usage #MUTIL_ASSERT(i == end_of_loop);#
 * @param  CONDITION   The condition to assert.
 * @see MUTIL_ERROR
 * @see mutil_abort
 * @see _mutil_errcode
 */
#define MUTIL_ASSERT(CONDITION) \
     if(!(CONDITION))\
     { \
         MUTIL_ERROR("Debug assertion failed " #CONDITION "\n"); \
         mutil_abort(); \
     }

#else /* ifdef DEBUG */

/*
 * The non-debug versions of the macros above -- expand to
 * null-operations if DEBUG is not defined.  See documentation above.
 */
#define MUTIL_TRACE(MSG)
#define MUTIL_WARN(MSG)
#define MUTIL_ERROR(MSG)
#define MUTIL_ASSERT(CONDITION)

#endif /* ifdef DEBUG */

#ifdef __cplusplus
}
#endif

#endif /* ifndef IN_UT_DEBUG_H */
