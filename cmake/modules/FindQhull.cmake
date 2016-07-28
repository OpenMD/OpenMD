###############################################################################
# Find QHULL
#
# This sets the following variables:
# QHULL_FOUND - True if QHULL was found.
# QHULL_INCLUDE_DIRS - Directories containing the QHULL include files.
# QHULL_LIBRARIES - Libraries needed to use QHULL.
# QHULL_DEFINITIONS - Compiler flags for QHULL.
# If QHULL_USE_STATIC is specified then look for static libraries ONLY else
# look for shared ones

find_program(QHULL_BINARY qhull PATHS ${QHULL_ROOT}/bin)

if(QHULL_BINARY STREQUAL "QHULL_BINARY-NOTFOUND")
    set(QHULL_FOUND FALSE)
else()
    set(QHULL_FOUND TRUE)
    execute_process(COMMAND "${QHULL_BINARY}" OUTPUT_VARIABLE _QHULL_OUT ERROR_VARIABLE _QHULL_OUT)
    string(REGEX MATCH "([0-9]+\\.[0-9]+)" QHULL_VERSION "${_QHULL_OUT}")
endif()

# qhull Versions
if (QHULL_FOUND)
  if(${QHULL_VERSION} STREQUAL "2003.1")
      set(QHULL_MAJOR_VERSION "4")
  elseif(${QHULL_VERSION} VERSION_GREATER "2009" AND ${QHULL_VERSION} VERSION_LESS "2010")
      set(QHULL_MAJOR_VERSION  "5")
  elseif(${QHULL_VERSION} VERSION_GREATER "2010" AND ${QHULL_VERSION} VERSION_LESS "2015" )
      set(QHULL_MAJOR_VERSION  "6")
  elseif(${QHULL_VERSION} VERSION_GREATER "2015")
      set(QHULL_MAJOR_VERSION  "7")
  endif()
else()
  set(QHULL_MAJOR_VERSION "7")
endif()

if(QHULL_USE_STATIC)
  set(QHULL_RELEASE_NAME qhullstatic)
  set(QHULL_DEBUG_NAME qhullstatic_d)
else(QHULL_USE_STATIC)
  # prefer reentrant, then pointer (qh_QHpointer = 1) version, then
  # static object (qh_QHpointer = 0) or unknown
  set(QHULL_RELEASE_NAME qhull_r qhull${QHULL_MAJOR_VERSION}_r qhull_r${QHULL_MAJOR_VERSION} 
                         qhull_p qhull${QHULL_MAJOR_VERSION}_p qhull_p${QHULL_MAJOR_VERSION} 
                         qhull qhull${QHULL_MAJOR_VERSION})
  set(QHULL_DEBUG_NAME qhull_rd qhull${QHULL_MAJOR_VERSION}_rd qhull_rd${QHULL_MAJOR_VERSION} 
                       qhull_pd qhull${QHULL_MAJOR_VERSION}_pd qhull_pd${QHULL_MAJOR_VERSION} 
                       qhull_d qhull${QHULL_MAJOR_VERSION}_d qhull_d${QHULL_MAJOR_VERSION})
endif(QHULL_USE_STATIC)

find_file(QHULL_HEADER
          NAMES libqhull_r/libqhull_r.h libqhull/libqhull.h qhull.h
          HINTS "${QHULL_ROOT}" "$ENV{QHULL_ROOT}" "${QHULL_INCLUDE_DIR}"
          PATHS "$ENV{PROGRAMFILES}/QHull" "$ENV{PROGRAMW6432}/QHull"
          PATH_SUFFIXES qhull src/libqhull libqhull_r libqhull include)

set(QHULL_HEADER "${QHULL_HEADER}" CACHE INTERNAL "QHull header" FORCE )

if(QHULL_HEADER)
    get_filename_component(qhull_header ${QHULL_HEADER} NAME_WE)
    if("${qhull_header}" STREQUAL "libqhull_r")
        set(HAVE_QHULL_REENTRANT_HEADERS ON)
        set(HAVE_QHULL_2011 OFF)
        get_filename_component(QHULL_INCLUDE_DIR ${QHULL_HEADER} PATH)
        get_filename_component(QHULL_INCLUDE_DIR ${QHULL_INCLUDE_DIR} PATH)
    elseif("${qhull_header}" STREQUAL "qhull")
        set(HAVE_QHULL_2011 OFF)
        set(HAVE_QHULL_REENTRANT_HEADERS OFF)
        get_filename_component(QHULL_INCLUDE_DIR ${QHULL_HEADER} PATH)
    elseif("${qhull_header}" STREQUAL "libqhull")
        set(HAVE_QHULL_2011 ON)
        set(HAVE_QHULL_REENTRANT_HEADERS OFF)
        get_filename_component(QHULL_INCLUDE_DIR ${QHULL_HEADER} PATH)
        get_filename_component(QHULL_INCLUDE_DIR ${QHULL_INCLUDE_DIR} PATH)
    endif()
else(QHULL_HEADER)
    set(QHULL_INCLUDE_DIR "QHULL_INCLUDE_DIR-NOTFOUND")
endif(QHULL_HEADER)

set(QHULL_INCLUDE_DIR "${QHULL_INCLUDE_DIR}" CACHE PATH "QHull include dir." FORCE)

find_library(QHULL_LIBRARY
             NAMES ${QHULL_RELEASE_NAME}
             HINTS "${QHULL_ROOT}" "$ENV{QHULL_ROOT}"
             PATHS "$ENV{PROGRAMFILES}/QHull" "$ENV{PROGRAMW6432}/QHull"
             PATH_SUFFIXES project build bin lib)

find_library(QHULL_LIBRARY_DEBUG
             NAMES ${QHULL_DEBUG_NAME} ${QHULL_RELEASE_NAME}
             HINTS "${QHULL_ROOT}" "$ENV{QHULL_ROOT}"
             PATHS "$ENV{PROGRAMFILES}/QHull" "$ENV{PROGRAMW6432}/QHull"
             PATH_SUFFIXES project build bin lib)

if(NOT QHULL_LIBRARY_DEBUG)
  set(QHULL_LIBRARY_DEBUG ${QHULL_LIBRARY})
endif(NOT QHULL_LIBRARY_DEBUG)

if(QHULL_LIBRARY)
    get_filename_component(qhull_library ${QHULL_LIBRARY} NAME_WE)
    if("${qhull_library}" MATCHES "libqhull[a-zA-Z0-9]*_[a-zA-Z0-9]*r[a-zA-Z0-9]*")
        set(HAVE_QHULL_REENTRANT_LIB ON)
        set(HAVE_QHULL_POINTER_LIB OFF)
    elseif("${qhull_library}" MATCHES "libqhull[a-zA-Z0-9]*_[a-zA-Z0-9]*p[a-zA-Z0-9]*")
        set(HAVE_QHULL_REENTRANT_LIB OFF)
        set(HAVE_QHULL_POINTER_LIB ON)
    else()
        # could still be pointer or reentrant version of library - we don't know
    endif()
endif(QHULL_LIBRARY)

set(QHULL_INCLUDE_DIRS ${QHULL_INCLUDE_DIR})
set(QHULL_LIBRARIES optimized ${QHULL_LIBRARY} debug ${QHULL_LIBRARY_DEBUG})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Qhull DEFAULT_MSG QHULL_LIBRARY QHULL_INCLUDE_DIR)

mark_as_advanced(QHULL_LIBRARY QHULL_LIBRARY_DEBUG QHULL_INCLUDE_DIR HAVE_QHULL_REENTRANT_HEADERS HAVE_QHULL_REENTRANT_LIB)

if(QHULL_FOUND)
  set(HAVE_QHULL ON)

  if(HAVE_QHULL_REENTRANT_HEADERS AND HAVE_QHULL_REENTRANT_LIB)
    set(HAVE_QHULL_REENTRANT ON)
  endif(HAVE_QHULL_REENTRANT_HEADERS AND HAVE_QHULL_REENTRANT_LIB)

  if (QHULL_USE_STATIC)
    add_definitions("-Dqh_QHpointer=0")
    if(MSVC)
      add_definitions("-Dqh_QHpointer_dllimport")
    endif(MSVC)
  endif(QHULL_USE_STATIC)

  if(NOT HAVE_QHULL_REENTRANT AND NOT QHULL_USE_STATIC)
    #possibly a dangerous assumption that if we don't have the reentrant headers and libs and 
    #if we also don't request static, then we want the pointer version:
    add_definitions("-Dqh_QHpointer=1")
  endif(NOT HAVE_QHULL_REENTRANT AND NOT QHULL_USE_STATIC)
  message(STATUS "QHULL found (include: ${QHULL_INCLUDE_DIRS}, lib: ${QHULL_LIBRARIES})")
endif(QHULL_FOUND)
