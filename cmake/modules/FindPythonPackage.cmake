# FIXME: Remove the code for CMake <3.12 once we have switched over completely.
# FindPython3 is a new module since CMake 3.12. It deprecates FindPythonInterp and FindPythonLibs.
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    # Use FindPythonInterp and FindPythonLibs for CMake <3.12
    find_package(PythonInterp 3.4 REQUIRED)

    # Define variables that are available in FindPython3, so there's no need to branch off in the later part.
    set(Python3_Interpreter_FOUND ${PYTHONINTERP_FOUND})
    set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
else()
    # Use FindPython3 for CMake >=3.12
    find_package(Python3 3.4 REQUIRED COMPONENTS Interpreter)
endif()

include(CMakeParseArguments)

function(FindPythonPackage)
    set(_options REQUIRED)
    set(_single_args PACKAGE_NAME MODULE_NAME VERSION_PROPERTY)

    cmake_parse_arguments("" "${_options}" "${_single_args}" "" ${ARGN})
    if(_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "FindPythonPackage called with unknown arguments ${FPP_UNPARSED_ARGUMENTS}")
    endif()

    if(NOT _PACKAGE_NAME)
        message(FATAL_ERROR "FindPythonPackage expects a PACKAGE_NAME argument")
    endif()

    if(NOT _MODULE_NAME)
        set(_MODULE_NAME ${_PACKAGE_NAME})
    endif()

    if(NOT _VERSION_PROPERTY)
        set(_VERSION_PROPERTY "__version__")
    endif()

    if(NOT Python3_Interpreter_FOUND)
        if(_REQUIRED)
            message(FATAL_ERROR "Could not find Python interpreter for required dependency ${_MODULE_NAME}")
        else()
            message(STATUS "Could not find ${_MODULE_NAME} because no Python interpreter was found")
            return()
        endif()
    endif()

    if(Python3_Interpreter_FOUND)
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import ${_MODULE_NAME}; print(${_MODULE_NAME}.__file__, ${_MODULE_NAME}.${_VERSION_PROPERTY})"
            RESULT_VARIABLE _process_status
            OUTPUT_VARIABLE _process_output
#             ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()

    set(_path "")
    set(_version "")
    if(_process_status EQUAL 0)
        string(REGEX MATCH "[^ ]+" _path ${_process_output})
        get_filename_component(_path ${_path} DIRECTORY)
        string(REGEX MATCH " .+" _version ${_process_output})
        string(STRIP ${_version} _version)
    endif()

    set(${_PACKAGE_NAME}_PATH "${_path}" PARENT_SCOPE)
    set(${_PACKAGE_NAME}_VERSION "${_version}" PARENT_SCOPE)
endfunction()
