# Top level cmake script for OpenMD. 

cmake_minimum_required (VERSION 2.6)
project(OpenMD)


# OpenMD version number.
set (VERSION_MAJOR "2")
set (VERSION_MINOR "0")
set (VERSION_TINY "0")  

# cmake modules that are needed to build oopse
set (CMAKE_MODULE_PATH ${OpenMD_SOURCE_DIR}/cmake/modules)

enable_language(CXX  C)


if ( APPLE )
set(MPI_COMPILER /opt/local/lib/openmpi/bin/mpic++)
set(MPI_LIBRARY "MPI_LIBRARY-NOTFOUND" )
endif( APPLE )
FIND_PACKAGE(MPI)
IF ( MPI_FOUND )
  SET(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_PATH};${CMAKE_REQUIRED_INCLUDES}")
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

FIND_PACKAGE(Subversion)
IF ( Subversion_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/.svn )
  Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
  MESSAGE("Current revision is ${Project_WC_REVISION}")
  Subversion_WC_LOG(${PROJECT_SOURCE_DIR} Project)
  SET(SVN_REV ${Project_WC_REVISION})
  ADD_DEFINITIONS( -DSVN_REV=${Project_WC_REVISION} )
ENDIF()

# Optional libraries: If we can find these, we will build with them
# Look for OpenBabel libraries
find_package(OpenBabel2)
IF(OPENBABEL2_FOUND)
  SET(USE_OPENBABEL)
  include_directories(${OPENBABEL2_INCLUDE_DIR})
  message(STATUS "Openbabel libs: " ${OPENBABEL2_LIBRARIES})
#  TARGET_LINK_LIBRARIES(${OPENBABEL2_LIBRARIES})
ENDIF(OPENBABEL2_FOUND)

#Look for QHULL Libraries
find_package(QHULL)
IF(QHULL_FOUND)
  SET(HAVE_QHULL 1)
  include_directories(${QHULL_INCLUDE_DIR})
  LINK_LIBRARIES(${QHULL_LIBRARIES})
ENDIF(QHULL_FOUND)


# zlib stuff
find_package(ZLIB)
if(ZLIB_FOUND)
  SET(HAVE_LIBZ 1)
  add_definitions(-DHAVE_LIBZ)
  include_directories(${ZLIB_INCLUDE_DIR})
  LINK_LIBRARIES(${ZLIB_LIBRARIES})
endif(ZLIB_FOUND)

#FFTW
INCLUDE(${CMAKE_MODULE_PATH}/FindFFTW.cmake)
IF(FFTW_FOUND)
  SET(HAVE_LIBFFTW 1)
  SET(HAVE_FFTW3_H)
  INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
  LINK_LIBRARIES(${FFTW_LIBRARIES})
ENDIF(FFTW_FOUND)


#include checks
include(CheckIncludeFile)
CHECK_INCLUDE_FILE(fftw3.h HAVE_FFTW3_H)
CHECK_INCLUDE_FILE(fftw.h HAVE_FFWT_H)
CHECK_INCLUDE_FILE(dfftw.h HAVE_DFFTW_H)
CHECK_INCLUDE_FILE(cmath.h HAVE_CMATH_H)
check_include_file(limits.h HAVE_LIMITS_H)
check_include_file(stdlib.h HAVE_STDLIB_H)
check_include_file(string.h HAVE_STRING_H)
check_include_file("sys/ioctl.h" HAVE_SYS_IOCTL_H)
check_include_file(unistd.h HAVE_UNISTD_H)

#symbol checks
include(CheckSymbolExists)

#function checks
include(CheckFunctionExists)
check_function_exists(atexit HAVE_ATEXIT)
check_function_exists(floor HAVE_FLOOR)
check_function_exists(pow HAVE_POW)
check_function_exists(select HAVE_SELECT)
check_function_exists(sqrt HAVE_SQRT)
check_function_exists(strcasecmp HAVE_STRCASECMP)
check_function_exists(strchr HAVE_STRCHR)
check_function_exists(strdup HAVE_STRDUP)
check_function_exists(strncasecmp HAVE_STRNCASECMP)
check_function_exists(strstr HAVE_STRSTR)
check_function_exists(strtol HAVE_STRTOL)
check_function_exists(strtoull HAVE_STRTOULL)
check_function_exists(sysmp HAVE_SYSMP)
check_function_exists(table HAVE_TABLE)


include(CMakePrintSystemInformation)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.cmake"
  "${PROJECT_SOURCE_DIR}/src/config.h"
  )
include_directories (${PROJECT_SOURCE_DIR}/src)

#CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)


#Add executables for build
set (PROGRAMS openmd openmd_MPI Dump2XYZ simpleBuilder StaticProps DynamicProps
randomBuilder nanoparticleBuilder thermalizer atom2md Hydro)
# Example of how to set core libraries
#set(CORELIBS ${GLUT_LIBRARY} ${OPENGL_LIBRARY} m)

#Add subirectories for build

add_subdirectory(lib)
add_subdirectory(bin)
