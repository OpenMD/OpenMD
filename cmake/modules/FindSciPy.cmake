# - Find the SciPy libraries
# This module finds if SciPy is installed, and sets the following variables
# indicating where it is.
#
#
#  SCIPY_FOUND               - was SciPy found
#  SCIPY_VERSION             - the version of SciPy found as a string
#  SCIPY_VERSION_MAJOR       - the major version number of SciPy
#  SCIPY_VERSION_MINOR       - the minor version number of SciPy
#  SCIPY_VERSION_PATCH       - the patch version number of SciPy
#  SCIPY_VERSION_DECIMAL     - e.g. version 1.6.1 is 10601
#  SCIPY_INCLUDE_DIRS        - path to the SciPy include files

#============================================================================
# Copyright 2013 Martin Köhler
# Copyright 2012 Continuum Analytics, Inc.
#
# MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
#============================================================================


# lint_cmake: -package/consistency


# Finding SciPy involves calling the Python interpreter
if(SciPy_FIND_REQUIRED)
    find_package(PythonInterp REQUIRED)
else()
    find_package(PythonInterp)
endif()

if(NOT PYTHONINTERP_FOUND)
    set(SCIPY_FOUND FALSE)
    return()
endif()

execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import scipy as s;  print(s.__version__); print(s.__file__);"
    RESULT_VARIABLE _SCIPY_SEARCH_SUCCESS
    OUTPUT_VARIABLE _SCIPY_VALUES_OUTPUT
    ERROR_VARIABLE _SCIPY_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT _SCIPY_SEARCH_SUCCESS MATCHES 0)
    if(SciPy_FIND_REQUIRED)
        message(FATAL_ERROR "SciPy import failure:\n${_SCIPY_ERROR_VALUE}")
    endif()
    set(SCIPY_FOUND FALSE)
    return()
endif()

# Convert the process output into a list
string(REGEX REPLACE ";" "\\\\;" _SCIPY_VALUES ${_SCIPY_VALUES_OUTPUT})
string(REGEX REPLACE "\n" ";" _SCIPY_VALUES ${_SCIPY_VALUES})
# Just in case there is unexpected output from the Python command.
list(GET _SCIPY_VALUES -2 SCIPY_VERSION)
list(GET _SCIPY_VALUES -1 SCIPY_INCLUDE_DIRS)

string(REGEX MATCH "^[0-9]+\\.[0-9]+\\.[0-9]+" _VER_CHECK "${SCIPY_VERSION}")
if("${_VER_CHECK}" STREQUAL "")
    # The output from Python was unexpected. Raise an error always
    # here, because we found SciPy, but it appears to be corrupted somehow.
    message(FATAL_ERROR
        "Requested version and include path from SciPy, got instead:\n${_SCIPY_VALUES_OUTPUT}\n")
    return()
endif()

# Make sure all directory separators are '/'
string(REGEX REPLACE "\\\\" "/" SCIPY_INCLUDE_DIRS ${SCIPY_INCLUDE_DIRS})
string(REGEX REPLACE "__init__.py[c]*" "" SCIPY_INCLUDE_DIRS ${SCIPY_INCLUDE_DIRS})

# Get the major and minor version numbers
string(REGEX REPLACE "\\." ";" _SCIPY_VERSION_LIST ${SCIPY_VERSION})
list(GET _SCIPY_VERSION_LIST 0 SCIPY_VERSION_MAJOR)
list(GET _SCIPY_VERSION_LIST 1 SCIPY_VERSION_MINOR)
list(GET _SCIPY_VERSION_LIST 2 SCIPY_VERSION_PATCH)
string(REGEX MATCH "[0-9]*" SCIPY_VERSION_PATCH ${SCIPY_VERSION_PATCH})
math(EXPR SCIPY_VERSION_DECIMAL
    "(${SCIPY_VERSION_MAJOR} * 10000) + (${SCIPY_VERSION_MINOR} * 100) + ${SCIPY_VERSION_PATCH}")

find_package_message(SCIPY
    "Found SciPy: version \"${SCIPY_VERSION}\" ${SCIPY_INCLUDE_DIRS}"
    "${SCIPY_INCLUDE_DIRS}${SCIPY_VERSION}")

set(SCIPY_FOUND TRUE)
