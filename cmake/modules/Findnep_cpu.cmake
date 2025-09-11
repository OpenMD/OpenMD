# Find nep_cpu
#
# Find the nep_cpu include directory
#
# if you need to add a custom library search path, do it via via CMAKE_PREFIX_PATH
#
# This module defines
#  nep_cpu_INCLUDE_DIRS, where to find headers, etc.
#  nep_cpu_FOUND, If false, do not try to use nep_cpu.

if(nep_cpu_INCLUDE_DIRS)
  # in cache already or user-specified
  set(nep_cpu_FOUND TRUE)
else()
  set(nep_cpu_FOUND FALSE)
  if(NOT nep_cpu_INCLUDE_DIR)
     find_path(nep_cpu_INCLUDE_DIR
               NAMES "nep.h dfdt3para.h"
               DOC "nep_cpu include dir"
               HINTS ${nep_cpu_ROOT} /usr/include)
    if(nep_cpu_INCLUDE_DIR)
      message(STATUS "Found nep_cpu include file at ${nep_cpu_INCLUDE_DIR}")
    endif()
  endif()

  if(nep_cpu_INCLUDE_DIR)
    set(nep_cpu_FOUND TRUE)
    set(nep_cpu_INCLUDE_DIRS ${nep_cpu_INCLUDE_DIR})
    message("Setting nep_cpu found ${nep_cpu_FOUND}")
  endif()

  mark_as_advanced(nep_cpu_FOUND nep_cpu_INCLUDE_DIR)
endif()
