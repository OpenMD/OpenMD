# Find stl_reader
#
# Find the stl_reader include directory
# 
# if you need to add a custom library search path, do it via via CMAKE_PREFIX_PATH 
# 
# This module defines
#  stl_reader_INCLUDE_DIRS, where to find headers, etc.
#  stl_reader_FOUND, If false, do not try to use stl_reader.

if(stl_reader_INCLUDE_DIR)
  # in cache already or user-specified
  set(stl_reader_FOUND TRUE)
else()

  if(NOT stl_reader_INCLUDE_DIR)
  	 find_path(stl_reader_INCLUDE_DIR
	   NAMES "stl_reader.h"
	   DOC "stl_reader include dir"
	   HINTS ${stl_reader_ROOT} /usr/include )
    if(stl_reader_INCLUDE_DIR)
      message(STATUS "Found stl_reader include file at ${stl_reader_INCLUDE_DIR}")
    endif()
  endif()
  if(stl_reader_INCLUDE_DIR)
    set(stl_reader_FOUND TRUE)
    set(stl_reader_INCLUDE_DIRS ${stl_reader_INCLUDE_DIR})
    message("Setting stl_reader found ${stl_reader_FOUND}")
  endif()

  mark_as_advanced(stl_reader_FOUND stl_reader_INCLUDE_DIR)
endif()
