#ifndef __GETOPT_H_
#define __GETOPT_H_
#ifndef _MSC_VER
/*	$NetBSD: getopt.h,v 1.4 2000/07/07 10:43:54 ad Exp $	*/
/*	$FreeBSD: src/include/getopt.h,v 1.6 2004/02/24 08:09:20 ache Exp $ */

/*-
 * Copyright (c) 2000 The NetBSD Foundation, Inc.
 * All rights reserved.
 *
 * This code is derived from software contributed to The NetBSD Foundation
 * by Dieter Baron and Thomas Klausner.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *        This product includes software developed by the NetBSD
 *        Foundation, Inc. and its contributors.
 * 4. Neither the name of The NetBSD Foundation nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE NETBSD FOUNDATION, INC. AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _GETOPT_H_
#define _GETOPT_H_

#include <sys/cdefs.h>
#include <unistd.h>

/*
 * GNU-like getopt_long()/getopt_long_only() with 4.4BSD optreset extension.
 * getopt() is declared here too for GNU programs.
 */
#define no_argument        0
#define required_argument  1
#define optional_argument  2

struct option {
	/* name of long option */
	const char *name;
	/*
	 * one of no_argument, required_argument, and optional_argument:
	 * whether option takes an argument
	 */
	int has_arg;
	/* if not NULL, set *flag to val when option found */
	int *flag;
	/* if flag not NULL, value to set *flag to; else return value */
	int val;
};

__BEGIN_DECLS
int	getopt_long(int, char * const *, const char *,
	const struct option *, int *);
int	getopt_long_only(int, char * const *, const char *,
	const struct option *, int *);
#ifndef _GETOPT
#define	_GETOPT
int	 getopt(int, char * const [], const char *) __DARWIN_ALIAS(getopt);

extern char *optarg;			/* getopt(3) external variables */
extern int optind, opterr, optopt;
#endif
#ifndef _OPTRESET
#define	_OPTRESET
extern int optreset;			/* getopt(3) external variable */
#endif
__END_DECLS
 
#endif /* !_GETOPT_H_ */
#else
/* Getopt for Microsoft C  
This code is a modification of the Free Software Foundation, Inc. 
Getopt library for parsing command line argument the purpose was
to provide a Microsoft Visual C friendly derivative. This code
provides functionality for both Unicode and Multibyte builds.

Date: 02/03/2011 - Ludvik Jerabek - Initial Release
Version: 1.0
Comment: Supports getopt, getopt_long, and getopt_long_only
and POSIXLY_CORRECT environment flag
License: LGPL

Revisions:

02/03/2011 - Ludvik Jerabek - Initial Release
02/20/2011 - Ludvik Jerabek - Fixed compiler warnings at Level 4
07/05/2011 - Ludvik Jerabek - Added no_argument, required_argument, optional_argument defs
08/03/2011 - Ludvik Jerabek - Fixed non-argument runtime bug which caused runtime exception
08/09/2011 - Ludvik Jerabek - Added code to export functions for DLL and LIB
02/15/2012 - Ludvik Jerabek - Fixed _GETOPT_THROW definition missing in implementation file

**DISCLAIMER**
THIS MATERIAL IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
EITHER EXPRESS OR IMPLIED, INCLUDING, BUT Not LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT. SOME JURISDICTIONS DO NOT ALLOW THE
EXCLUSION OF IMPLIED WARRANTIES, SO THE ABOVE EXCLUSION MAY NOT
APPLY TO YOU. IN NO EVENT WILL I BE LIABLE TO ANY PARTY FOR ANY
DIRECT, INDIRECT, SPECIAL OR OTHER CONSEQUENTIAL DAMAGES FOR ANY
USE OF THIS MATERIAL INCLUDING, WITHOUT LIMITATION, ANY LOST
PROFITS, BUSINESS INTERRUPTION, LOSS OF PROGRAMS OR OTHER DATA ON
YOUR INFORMATION HANDLING SYSTEM OR OTHERWISE, EVEN If WE ARE
EXPRESSLY ADVISED OF THE POSSIBILITY OF SUCH DAMAGES. 
*/

#ifdef _GETOPT_API
#undef _GETOPT_API
#endif

#if defined(EXPORTS_GETOPT) && defined(STATIC_GETOPT)
#error "The preprocessor definitions of EXPORTS_GETOPT and STATIC_GETOPT can only be used individually"
#elif defined(STATIC_GETOPT)
#pragma message("Warning static builds of getopt violate the Lesser GNU Public License")
#define _GETOPT_API
#elif defined(EXPORTS_GETOPT)
#pragma message("Exporting getopt library")
#define _GETOPT_API __declspec(dllexport)	
#else
#pragma message("Importing getopt library")
#define _GETOPT_API __declspec(dllimport)
#endif


#include <tchar.h>

// Standard GNU options
#define	null_argument		0 /*Argument Null*/
#define	no_argument		0 /*Argument Switch Only*/
#define required_argument	1 /*Argument Required*/
#define optional_argument	2 /*Argument Optional*/

// Shorter Versions of options
#define ARG_NULL 0 /*Argument Null*/
#define ARG_NONE 0 /*Argument Switch Only*/
#define ARG_REQ 1  /*Argument Required*/
#define ARG_OPT 2  /*Argument Optional*/

// Change behavior for C\C++
#ifdef __cplusplus
#define _BEGIN_EXTERN_C extern "C" {
#define _END_EXTERN_C }
#define _GETOPT_THROW throw()
#else
#define _BEGIN_EXTERN_C
#define _END_EXTERN_C
#define _GETOPT_THROW
#endif

_BEGIN_EXTERN_C

	extern _GETOPT_API TCHAR *optarg;
extern _GETOPT_API int optind;
extern _GETOPT_API int opterr;
extern _GETOPT_API int optopt;

struct option
{
	const TCHAR* name;
	int has_arg;
	int *flag;
	TCHAR val;
};

extern _GETOPT_API int getopt(int argc, TCHAR *const *argv, const TCHAR *optstring) _GETOPT_THROW;
extern _GETOPT_API int getopt_long(int ___argc, TCHAR *const *___argv, const TCHAR *__shortopts, const struct option *__longopts, int *__longind) _GETOPT_THROW;
extern _GETOPT_API int getopt_long_only(int ___argc, TCHAR *const *___argv, const TCHAR *__shortopts, const struct option *__longopts, int *__longind) _GETOPT_THROW;
_END_EXTERN_C

// Undefine so the macros are not included
#undef _BEGIN_EXTERN_C
#undef _END_EXTERN_C
#undef _GETOPT_THROW
#undef _GETOPT_API

#endif
#endif  // __GETOPT_H_
