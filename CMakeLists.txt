# Top level cmake script for OpenMD. 

cmake_minimum_required (VERSION 2.6)
project(OpenMD)


# OpenMD version number.
set (VERSION_MAJOR "2")
set (VERSION_MINOR "0")
set (VERSION_TINY "0")  

# cmake modules that are needed to build oopse
set (CMAKE_MODULE_PATH ${OpenMD_SOURCE_DIR}/cmake/modules)


enable_language( Fortran C  CXX)



# Figure out compiler stuff. We need to know which c, c++ and fortran
#  compiler we are using as well as the name magling stuff.
INCLUDE(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_HEADER( ${PROJECT_SOURCE_DIR}/src/config_f.h
                            MACRO_NAMESPACE "FC_"
                            SYMBOL_NAMESPACE "FC_")



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

#Figure out MPI stuff
find_package(MPI)
if(MPI_FOUND)
SET(IS_MPI 1)
#add_definitions(-DIS_MPI)
endif(MPI_FOUND)

#Add subversion Information
find_package(Subversion)
if(Subversion_FOUND)
  Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} ER)
  set(SUBVERSION_REVISION ${ER_WC_REVISION})
endif(Subversion_FOUND)
set(CPACK_PACKAGE_VERSION_PATCH "${SUBVERSION_REVISION}")

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
check_include_file (errno.h HAVE_ERRNO_H)
check_include_file (inttypes.h HAVE_INTTYPES_H)
check_include_file (limits.h HAVE_LIMITS_H)
check_include_file (math.h HAVE_MATH_H)
check_include_file (memory.h HAVE_MEMORY_H)
check_include_file (stdint.h HAVE_STDINT_H)
check_include_file (stdio.h HAVE_STDIO_H)
check_include_file (stdlib.h HAVE_STDLIB_H)
check_include_file (strings.h HAVE_STRINGS_H)
check_include_file (string.h HAVE_STRING_H)
check_include_file ("sys/param.h" HAVE_SYS_PARAM_H)
check_include_file ("sys/pstat.h" HAVE_SYS_PSTAT_H)
check_include_file ("sys/select.h" HAVE_SYS_SELECT_H)
check_include_file ("sys/socket.h" HAVE_SYS_SOCKET_H)
check_include_file ("sys/stat.h" HAVE_SYS_STAT_H)
check_include_file ("sys/sysctl.h" HAVE_SYS_SYSCTL_H)
check_include_file ("sys/sysinfo.h" HAVE_SYS_SYSINFO_H)
check_include_file ("sys/param.h" HAVE_SYS_PARAM_H)
check_include_file ("sys/sysmp.h" HAVE_SYS_SYSMP_H)
check_include_file ("sys/systemcfg.h" HAVE_SYS_SYSTEMCFG_H)
check_include_file ("sys/table.h" HAVE_SYS_TABLE_H)



find_path (HAVE_SYS_STAT_H stat.h
  PATHS ${include_locations}
  PATH_SUFFIXES sys
)

find_path (HAVE_SYS_TYPES_H types.h
  PATHS ${include_locations}
  PATH_SUFFIXES sys
)
find_path (HAVE_UNISTD_H unistd.h ${include_locations})

#symbol checks
include(CheckSymbolExists)


#function checks
include(CheckFunctionExists)
check_function_exists(floor HAVE_FLOOR)
check_function_exists (memset HAVE_MEMSET)
check_function_exists (pow HAVE_POW)
check_function_exists (sqrt HAVE_SQRT)
check_function_exists (strchr HAVE_STRCHR)
check_function_exists (strstr HAVE_STRSTR)
check_function_exists (sysmp HAVE_SYSMP)
check_function_exists (table HAVE_TABLE)



include(CMakePrintSystemInformation)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.cmake"
  "${PROJECT_SOURCE_DIR}/src/config.h"
  )
include_directories (${PROJECT_SOURCE_DIR}/src)

#CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)



#Add executables for build
set (PROGRAMS openmd Dump2XYZ simpleBuilder StaticProps DynamicProps
randomBuilder nanoparticleBuilder thermalizer atom2md Hydro)
# Example of how to set core libraries
#set(CORELIBS ${GLUT_LIBRARY} ${OPENGL_LIBRARY} m)

#foreach(program ${PROGRAMS})
#  add_executable(${program} ${program}.cpp)
#  target_link_libraries(${program} ${CORELIBS})
#endforeach(program)

#add_subdirectory(src/lattice)
#add_subdirectory(src/brains)
#add_subdirectory(src/UseTheForce)


#Add subirectories for build


#add_subdirectory(src)
add_subdirectory(lib)
add_subdirectory(bin)
