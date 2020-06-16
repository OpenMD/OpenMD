include(FindPythonPackage)

FindPythonPackage(PACKAGE_NAME "NumPy" MODULE_NAME "numpy")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NumPy REQUIRED_VARS NumPy_PATH VERSION_VAR NumPy_VERSION)
