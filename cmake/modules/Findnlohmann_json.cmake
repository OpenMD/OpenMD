# Find nlohmann_json
#
# Find the nlohmann_json include directory
# 
# if you need to add a custom library search path, do it via via CMAKE_PREFIX_PATH 
# 
# This module defines
#  nlohmann_json_INCLUDE_DIRS, where to find headers, etc.
#  nlohmann_json_FOUND, If false, do not try to use nlohmann_json.

if(nlohmann_json_INCLUDE_DIR)
  # in cache already or user-specified
  set(nlohmann_json_FOUND TRUE)
else()

  if(NOT nlohmann_json_INCLUDE_DIR)
  	 find_path(nlohmann_json_INCLUDE_DIR
	   NAMES "nlohmann/json.hpp"
	   DOC "nlohmann_json include dir"
	   HINTS ${nlohmann_json_ROOT}/include  /usr/include )
    if(nlohmann_json_INCLUDE_DIR)
      message(STATUS "Found nlohmann_json include files at ${nlohmann_json_INCLUDE_DIR}")
    endif()
  endif()
  if(nlohmann_json_INCLUDE_DIR)
    set(nlohmann_json_FOUND TRUE)
    set(nlohmann_json_INCLUDE_DIRS ${nlohmann_json_INCLUDE_DIR})
    message("Setting nlohmann_json found ${nlohmann_json_FOUND}")
  endif()

  mark_as_advanced(nlohmann_json_FOUND nlohmann_json_INCLUDE_DIR)
endif()
