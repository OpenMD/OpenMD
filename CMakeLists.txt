# Top level cmake script for OpenMD. 

cmake_minimum_required (VERSION 2.6)
project(OpenMD)


# OpenMD version number.
set (OpenMD_VERSION_MAJOR 2)
set (OpenMD_VERSION_MINOR 0)
set (OpenMD_VERSION_TINY 0)  q

# cmake modules that are needed to build oopse
set (CMAKE_MODULE_PATH ${OpenMD_SOURCE_DIR}/cmake/modules)


enable_language( Fortran C  CXX)


configure_file (
  "${PROJECT_SOURCE_DIR}/Config.h.in"
  "${PROJECT_BINARY_DIR}/Config.h"
  )

# Figure out compiler stuff. We need to know which c, c++ and fortran
  compiler we are using as well as the name magling stuff.
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_" FROM_LIBRARY blas[daxpy])




# Optional libraries: If we can find these, we will build with them
# Look for OpenBabel libraries
find_package(OpenBabel)
#Look for QHULL Libraries
find_package(QHULL)


#Figure out MPI stuff


#Add executables for build
set (PROGRAMS openmd Dump2XYZ simpleBuilder StaticProps DynamicProps
randomBuilder nanoparticleBuilder thermalizer atom2md Hydro)
# Example of how to set core libraries
#set(CORELIBS ${GLUT_LIBRARY} ${OPENGL_LIBRARY} m)

foreach(program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${CORELIBS})
endforeach(program)




#Add subirectories for build


add_subdirectory(src)
