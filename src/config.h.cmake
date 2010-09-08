/* config.h.  Generated from config.h.cmake by cmake  */


#define OPENMD_VERSION_MAJOR "${VERSION_MAJOR}"
#define OPENMD_VERSION_MINOR "${VERSION_MINOR}"
#define OPENMD_VERSION_TINY  "${VERSION_TINY}"

#define MK_STR(s) # s
#define STR_DEFINE(t, s) t = MK_STR(s)
/* The file extension used for shared modules */
#undef MODULE_EXTENSION


#include "config_f.h"
#define FC_FUNC FC_GLOBAL
#define FC_FUNC_ FC_GLOBAL_

/* Is defined if OpenMD should be compiled with single precision arithmetic. */
#cmakedefine SINGLE_PRECISION

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
#cmakedefine FC_DUMMY_MAIN

/* Define if F77 and FC dummy `main' functions are identical. */
#cmakedefine FC_DUMMY_MAIN_EQ_F77

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#cmakedefine FC_FUNC

/* As FC_FUNC, but for C identifiers containing underscores. */
#cmakedefine FC_FUNC_

/* define if fstream::open() accepts third parameter. */
#cmakedefine FSTREAM_OPEN_PROT

/* Is defined if the qhull library is available. */
#cmakedefine HAVE_QHULL

/* Define to 1 if the system has the type `clock_t'. */
#cmakedefine HAVE_CLOCK_T 1

/* Define to 1 if you have the <cmath> header file. */
#cmakedefine HAVE_CMATH 1

/* Define to 1 if you have the <conio.h> header file. */
#cmakedefine HAVE_CONIO_H 1

/* Define to 1 if you have the <ctype.h> header file. */
#cmakedefine HAVE_CTYPE_H 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#cmakedefine HAVE_DLFCN_H 1

/* define if fftw3.h exists */
#cmakedefine HAVE_FFTW3_H

/* define if fftw.h exists */
#cmakedefine HAVE_FFTW_H

/* define if dfftw.h exists */
#cmakedefine HAVE_DFFTW_H



/* Define to 1 if you have the `getsysinfo' function. */
#cmakedefine HAVE_GETSYSINFO 1

/* Define to 1 if you have the <libintl.h> header file. */
#cmakedefine HAVE_LIBINTL_H 1

/* Define to 1 if you have the `z' library (-lz). */
#cmakedefine HAVE_LIBZ 1

/* Define to 1 if you have the <limits.h> header file. */
#cmakedefine HAVE_LIMITS_H 1

/* Define to 1 if you have the <machine/hal_sysinfo.h> header file. */
#cmakedefine HAVE_MACHINE_HAL_SYSINFO_H 1



/* Define to 1 if you have the `pstat_getdynamic' function. */
#cmakedefine HAVE_PSTAT_GETDYNAMIC 1

/* Define to 1 if you have the `pstat_getstatic' function. */
#cmakedefine HAVE_PSTAT_GETSTATIC 1


/* Define to 1 if you have the `strcasecmp' function. */
#cmakedefine HAVE_STRCASECMP 1


/* Define to 1 if you have the `stricmp' function. */
#cmakedefine HAVE_STRICMP 1


/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H 1

/* Define to 1 if you have the `sysmp' function. */
#cmakedefine HAVE_SYSMP 1

/* Define to 1 if you have the <sys/param.h> header file. */
#cmakedefine HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/pstat.h> header file. */
#cmakedefine HAVE_SYS_PSTAT_H 1

/* Define to 1 if you have the <sys/select.h> header file. */
#cmakedefine HAVE_SYS_SELECT_H 1

/* Define to 1 if you have the <sys/socket.h> header file. */
#cmakedefine HAVE_SYS_SOCKET_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/sysctl.h> header file. */
#cmakedefine HAVE_SYS_SYSCTL_H 1

/* Define to 1 if you have the <sys/sysinfo.h> header file. */
#cmakedefine HAVE_SYS_SYSINFO_H 1

/* Define to 1 if you have the <sys/sysmp.h> header file. */
#cmakedefine HAVE_SYS_SYSMP_H 1

/* Define to 1 if you have the <sys/systemcfg.h> header file. */
#cmakedefine HAVE_SYS_SYSTEMCFG_H 1

/* Define to 1 if you have the <sys/table.h> header file. */
#cmakedefine HAVE_SYS_TABLE_H 1

/* Define to 1 if you have the `table' function. */
#cmakedefine HAVE_TABLE 1


/* Define if you have the _system_configuration variable. */
#cmakedefine HAVE__SYSTEM_CONFIGURATION

/* Define to the one symbol short name of this package. */
#cmakedefine PACKAGE_TARNAME

/* Define to the version of this package. */
#cmakedefine PACKAGE_VERSION

/* Path to ps program */
#cmakedefine PSCOMMAND

/* ps uses BSD-style arguments */
#define PSTYPE_IS_BSD


/* Define to the type of arg 1 for `select'. */
#cmakedefine SELECT_TYPE_ARG1

/* Define to the type of args 2, 3 and 4 for `select'. */
#cmakedefine SELECT_TYPE_ARG234

/* Define to the type of arg 5 for `select'. */
#cmakedefine SELECT_TYPE_ARG5

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#cmakedefine TIME_WITH_SYS_TIME

/* needed by DEC/Compaq/HP cxx to activate ANSI standard iostream. */
#cmakedefine __USE_STD_IOSTREAM

/* Define to empty if `const' does not conform to ANSI C. */
#cmakedefine const

/* Code compiled in debug mode */
#cmakedefine debug

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
#undef inline
#endif

/* Define to rpl_malloc if the replacement function should be used. */
#undef malloc

/* Define to rpl_realloc if the replacement function should be used. */
#undef realloc

/* Define to `unsigned int' if <sys/types.h> does not define. */
#undef size_t

#ifdef SINGLE_PRECISION
typedef float RealType;
#ifdef IS_MPI
#define MPI_REALTYPE MPI_FLOAT
#define REALTYPE FLOAT
#define REALTYPE_INT FLOAT_INT
#endif
#else
typedef double RealType;
#ifdef IS_MPI
#define MPI_REALTYPE MPI_DOUBLE
#define REALTYPE DOUBLE
#define REALTYPE_INT DOUBLE_INT
#endif
#endif
