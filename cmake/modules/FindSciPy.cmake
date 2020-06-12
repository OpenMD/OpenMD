include(FindPythonPackage)

FindPythonPackage(PACKAGE_NAME "SciPy" MODULE_NAME "scipy")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SciPy REQUIRED_VARS SciPy_PATH VERSION_VAR SciPy_VERSION)
