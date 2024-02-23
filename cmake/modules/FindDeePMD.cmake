# FindDeePMD.cmake
# Try to find DeePMD-kit's headers and libraries
# Defines:
#
#  DeePMD_FOUND - system has DeePMD-kit
#  DeePMD_INCLUDE_DIR - the DeePMD-kit include directory
#  DeePMD_LIBRARIES - Link these to use DeePMD-kit
#  IF DeePMD_DIR is defined, will look there first

if(OPENBABEL_INCLUDE_DIR AND DeePMD_LIBRARIES)
  # in cache already or user-specified
  set(DeePMD_FOUND TRUE)

else()

  if(NOT DeePMD_INCLUDE_DIR)
      find_path(DeePMD_INCLUDE_DIR deepmd/DeepPot.h
        PATHS
        ${DeePMD_DIR}/include/deepmd
        ${DeePMD_DIR}/include
        $ENV{DeePMD_INCLUDE_DIR}/deepmd
        $ENV{DeePMD_INCLUDE_DIR}
        $ENV{DeePMD_INCLUDE_PATH}/deepmd
        $ENV{DeePMD_INCLUDE_PATH}
        $ENV{DeePMD_DIR}/include/deepmd
        $ENV{DeePMD_DIR}/include
        $ENV{DeePMD_PATH}/include/deepmd
        $ENV{DeePMD_PATH}/include
        $ENV{deepmd_root}/include/deepmd
        $ENV{deepmd_root}/include
        /usr/include/deepmd
        /usr/include
        /usr/local/include/deepmd
        /usr/local/include
        /usr/local/deepmd-kit/include/deepmd
        /usr/local/deepmd-kit/include
        /usr/local/deepmd-kit/include/deepmd
        /usr/local/deepmd-kit/include
        /opt/local/include/deepmd
        /opt/local/include
        /opt/local/deepmd-kit/include/deepmd
        /opt/local/deepmd-kit/include
        /opt/local/deepmd-kit/include/deepmd
        /opt/local/deepmd-kit/include
        /opt/homebrew/include/deepmd
        /opt/homebrew/include
        /opt/homebrew/deepmd-kit/include/deepmd
        /opt/homebrew/deepmd-kit/include
        /opt/homebrew/deepmd-kit/include/deepmd
        /opt/homebrew/deepmd-kit/include
        ~/include/deepmd-kit
        ~/include
      )
    if(DeePMD_INCLUDE_DIR)
      message(STATUS "Found DeepMD-kit include files at ${DeePMD_INCLUDE_DIR}")
    endif()
  endif()

  if(NOT DeePMD_LIBRARIES)
  find_library(DeePMD_LIBRARIES NAMES deepmd_cc deepmd_c deepmd
      HINTS
      ${DeePMD_DIR}/lib
      ${DeePMD_DIR}/windows-vc2008/build/src/Release
      ${DeePMD_DIR}/build/src/Release
      PATHS
      $ENV{DeePMD_DIR}/lib
      $ENV{DeePMD_DIR}/windows-vc2008/build/src/Release
      $ENV{DeePMD_PATH}/lib
      $ENV{DeePMD_BASE}/lib
      /usr/lib
      /usr/local/lib
      ~/lib
      $ENV{LD_LIBRARY_PATH}
    )
    if(DeePMD_LIBRARIES)
      message(STATUS "Found DeePMD library at ${DeePMD_LIBRARIES}")
    endif()
  endif()

  if(DeePMD_INCLUDE_DIR AND DeePMD_LIBRARIES)
    set(DeePMD_FOUND TRUE)
  endif()

  mark_as_advanced(DeePMD_INCLUDE_DIR DeePMD_LIBRARIES)
endif()

