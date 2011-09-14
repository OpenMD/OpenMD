/* config.h.  Generated from config.h.cmake by cmake  */

#define OPENMD_VERSION_MAJOR "${VERSION_MAJOR}"
#define OPENMD_VERSION_MINOR "${VERSION_MINOR}"
#define OPENMD_VERSION_TINY  "${VERSION_TINY}"

#define MK_STR(s) # s
#define STR_DEFINE(t, s) t = MK_STR(s)

/* Is defined if OpenMD should be compiled with single precision arithmetic. */
#cmakedefine SINGLE_PRECISION

/* Is defined if the qhull library is available. */
#cmakedefine HAVE_QHULL

/* Define to 1 if you have the <conio.h> header file. */
#cmakedefine HAVE_CONIO_H 1

/* define if fftw3.h exists */
#cmakedefine HAVE_FFTW3_H

/* define if fftw.h exists */
#cmakedefine HAVE_FFTW_H

/* define if dfftw.h exists */
#cmakedefine HAVE_DFFTW_H

/* Define to 1 if you have the `z' library (-lz). */
#cmakedefine HAVE_LIBZ 1

/* Define to 1 if you have the `strcasecmp' function. */
#cmakedefine HAVE_STRCASECMP 1

/* Define to 1 if you have the `stricmp' function. */
#cmakedefine HAVE_STRICMP 1

/* Define to the one symbol short name of this package. */
#cmakedefine PACKAGE_TARNAME

/* Define to the version of this package. */
#cmakedefine PACKAGE_VERSION

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
