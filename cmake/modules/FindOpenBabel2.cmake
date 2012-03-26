# - Try to find OpenBabel2
# Once done this will define
#
#  OPENBABEL2_FOUND - system has OpenBabel2
#  OPENBABEL2_INCLUDE_DIR - the OpenBabel2 include directory
#  OPENBABEL2_LIBRARIES - Link these to use OpenBabel2
#
# A simplified version of FindOpenBabel2.cmake which doesn't rely
# on PkgConfig and doesn't search for the executable. 
#
# Copyright (c) 2006, 2007 Carsten Niehaus, <cniehaus@gmx.de>
# Copyright (C) 2008 Marcus D. Hanwell <marcus@cryos.org>
# Copyright (C) 2012 J. Daniel Gezelter <gezelter@openscience.org>
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.



FIND_PATH(OPENBABEL2_INCLUDE_DIR openbabel-2.0/openbabel/obconversion.h 
          HINTS "${_obDir}" "${GNUWIN32_DIR}" "${OPENBABEL2_ROOT}" "$ENV{OPENBABEL2_ROOT}" "$ENV{OPENBABEL2_INCLUDE_DIR}"
          PATH_SUFFIXES include )

if(OPENBABEL2_INCLUDE_DIR)
    set(OPENBABEL2_INCLUDE_DIR ${OPENBABEL2_INCLUDE_DIR}/openbabel-2.0)
endif(OPENBABEL2_INCLUDE_DIR)


FIND_LIBRARY(OPENBABEL2_LIBRARY NAMES openbabel openbabel2 
               HINTS "${_obDir}" "${GNUWIN32_DIR}" "${OPENBABEL2_ROOT}" "$ENV{OPENBABEL2_ROOT}" "$ENV{OPENBABEL2_LIBRARIES}"
               PATH_SUFFIXES project build bin lib lib64 )

# handle the QUIETLY and REQUIRED arguments and set OPENBABEL2_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OpenBabel2 DEFAULT_MSG OPENBABEL2_LIBRARY OPENBABEL2_INCLUDE_DIR)


if(OPENBABEL2_FOUND)
  set( OPENBABEL2_LIBRARIES ${OPENBABEL2_LIBRARY} )
  set( OPENBABEL2_INCLUDE_DIRS ${OPENBABEL2_INCLUDE_DIR} )
endif()


MARK_AS_ADVANCED(OPENBABEL2_INCLUDE_DIR OPENBABEL2_LIBRARY)

