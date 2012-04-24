/* config.h.  Generated from config.h.cmake by CMake for @PROJECT_NAME@  */

#define OPENMD_VERSION_MAJOR "${VERSION_MAJOR}"
#define OPENMD_VERSION_MINOR "${VERSION_MINOR}"
#define OPENMD_VERSION_TINY  "${VERSION_TINY}"

#define MK_STR(s) # s
#define STR_DEFINE(t, s) t = MK_STR(s)

/* Is defined if OpenMD should be compiled with single precision arithmetic. */
#cmakedefine SINGLE_PRECISION

/* Is defined if the qhull library is available. */

#cmakedefine HAVE_QHULL 1
#cmakedefine HAVE_QHULL_2011 1
#ifdef DISABLE_QHULL
#undef HAVE_QHULL
#endif

/* have <conio.h> */
#cmakedefine HAVE_CONIO_H 1

/* have symbol strncasecmp */
#cmakedefine HAVE_STRNCASECMP 1

/* define if fftw3.h exists */
#cmakedefine HAVE_FFTW3_H 1 

/* Define to 1 if you have the `z' library (-lz). */
#cmakedefine HAVE_LIBZ 1

/* Define to the one symbol short name of this package. */
#cmakedefine PACKAGE_TARNAME

/* Define to the version of this package. */
#cmakedefine PACKAGE_VERSION

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
