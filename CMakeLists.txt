# Top level cmake script for OpenMD. 

cmake_minimum_required (VERSION 2.6)
project(OpenMD)

# OpenMD version number.
set (VERSION_MAJOR "2")
set (VERSION_MINOR "0")
set (VERSION_TINY "0")  

# cmake modules that are needed to build oopse
SET(CMAKE_MODULE_PATH ${OpenMD_SOURCE_DIR}/cmake/modules)
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX
    "/usr/local/openmd" CACHE PATH "OpenMD install prefix" FORCE
    )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# ---------- Setup output Directories -------------------------
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/lib
  CACHE PATH
  "Single directory for all shared libraries"
  )
# --------- Setup the Executable output Directory -------------
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/bin
  CACHE PATH
  "Single directory for all Executables."
  )
# --------- Setup the static library directory -------------
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/lib
  CACHE PATH
  "Single directory for all static libraries."
  )


enable_language(CXX  C)

if ( APPLE )
set(MPI_COMPILER /opt/local/lib/openmpi/bin/mpic++)
set(MPI_LIBRARY "MPI_LIBRARY-NOTFOUND" )
endif( APPLE )
FIND_PACKAGE(MPI)
IF ( MPI_FOUND )
  SET(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_PATH};${CMAKE_REQUIRED_INCLUDES}")
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF()

FIND_PACKAGE(Subversion)
IF ( Subversion_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/.svn )
  Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
  MESSAGE("Current revision is ${Project_WC_REVISION}")
  Subversion_WC_LOG(${PROJECT_SOURCE_DIR} Project)
  SET(SVN_REV ${Project_WC_REVISION})
  ADD_DEFINITIONS( -DSVN_REV=${Project_WC_REVISION} )
ENDIF()

# Optional libraries: If we can find these, we will build with them
# Look for OpenBabel libraries
find_package(OpenBabel2)
IF(OPENBABEL2_FOUND)
  SET(USE_OPENBABEL)
  include_directories(${OPENBABEL2_INCLUDE_DIR})
  message(STATUS "Openbabel libs: " ${OPENBABEL2_LIBRARIES})
#  TARGET_LINK_LIBRARIES(${OPENBABEL2_LIBRARIES})
ENDIF(OPENBABEL2_FOUND)

#Look for QHULL Libraries
find_package(QHULL)
IF(QHULL_FOUND)
  SET(HAVE_QHULL 1)
  include_directories(${QHULL_INCLUDE_DIR})
  LINK_LIBRARIES(${QHULL_LIBRARIES})
ENDIF(QHULL_FOUND)


# zlib stuff
find_package(ZLIB)
if(ZLIB_FOUND)
  SET(HAVE_LIBZ 1)
  add_definitions(-DHAVE_LIBZ)
  include_directories(${ZLIB_INCLUDE_DIR})
  LINK_LIBRARIES(${ZLIB_LIBRARIES})
endif(ZLIB_FOUND)

#FFTW
INCLUDE(${CMAKE_MODULE_PATH}/FindFFTW.cmake)
IF(FFTW_FOUND)
  SET(HAVE_LIBFFTW 1)
  SET(HAVE_FFTW3_H)
  INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
  LINK_LIBRARIES(${FFTW_LIBRARIES})
ENDIF(FFTW_FOUND)


#include checks
include(CheckIncludeFile)
CHECK_INCLUDE_FILE(fftw3.h HAVE_FFTW3_H)
CHECK_INCLUDE_FILE(fftw.h HAVE_FFWT_H)
CHECK_INCLUDE_FILE(dfftw.h HAVE_DFFTW_H)
CHECK_INCLUDE_FILE(cmath.h HAVE_CMATH_H)
check_include_file(limits.h HAVE_LIMITS_H)
check_include_file(stdlib.h HAVE_STDLIB_H)
check_include_file(string.h HAVE_STRING_H)
check_include_file("sys/ioctl.h" HAVE_SYS_IOCTL_H)
check_include_file(unistd.h HAVE_UNISTD_H)

#symbol checks
include(CheckSymbolExists)

#function checks
include(CheckFunctionExists)
check_function_exists(atexit HAVE_ATEXIT)
check_function_exists(floor HAVE_FLOOR)
check_function_exists(pow HAVE_POW)
check_function_exists(select HAVE_SELECT)
check_function_exists(sqrt HAVE_SQRT)
check_function_exists(strcasecmp HAVE_STRCASECMP)
check_function_exists(strchr HAVE_STRCHR)
check_function_exists(strdup HAVE_STRDUP)
check_function_exists(strncasecmp HAVE_STRNCASECMP)
check_function_exists(strstr HAVE_STRSTR)
check_function_exists(strtol HAVE_STRTOL)
check_function_exists(strtoull HAVE_STRTOULL)
check_function_exists(sysmp HAVE_SYSMP)
check_function_exists(table HAVE_TABLE)


FIND_PACKAGE(Perl)
IF(PERL_FOUND)
  SET(PERL ${PERL_EXECUTABLE})
ELSE(PERL_FOUND)
  MESSAGE(STATUS "Failed to find perl")
ENDIF(PERL_FOUND)

FIND_PACKAGE(PythonInterp)
IF(PYTHON_EXECUTABLE)
  SET(PYTHON ${PYTHON_EXECUTABLE})
ELSE(PYTHON_EXECUTABLE)
  MESSAGE(STATUS "Failed to find python")
ENDIF(PYTHON_EXECUTABLE)

SET(PERL_INSTALLDIRS "site" CACHE STRING "Perl installation locations")


include(CMakePrintSystemInformation)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.cmake"
  "${PROJECT_SOURCE_DIR}/src/config.h"
  )
include_directories (${PROJECT_SOURCE_DIR}/src)

#CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)


#Add executables for build
set (PROGRAMS 
openmd
openmd_MPI
Dump2XYZ 
simpleBuilder 
StaticProps 
DynamicProps
randomBuilder 
nanoparticleBuilder 
thermalizer 
atom2md 
Hydro)

set(SOURCE
src/antlr/ANTLRUtil.cpp
src/antlr/ASTFactory.cpp
src/antlr/ASTNULLType.cpp
src/antlr/ASTRefCount.cpp
src/antlr/BaseAST.cpp
src/antlr/BitSet.cpp
src/antlr/CharBuffer.cpp
src/antlr/CharScanner.cpp
src/antlr/CommonAST.cpp
src/antlr/CommonASTWithHiddenTokens.cpp
src/antlr/CommonHiddenStreamToken.cpp
src/antlr/CommonToken.cpp
src/antlr/InputBuffer.cpp
src/antlr/LLkParser.cpp
src/antlr/MismatchedCharException.cpp
src/antlr/MismatchedTokenException.cpp
src/antlr/NoViableAltException.cpp
src/antlr/NoViableAltForCharException.cpp
src/antlr/Parser.cpp
src/antlr/RecognitionException.cpp
src/antlr/String.cpp
src/antlr/Token.cpp
src/antlr/TokenBuffer.cpp
src/antlr/TokenRefCount.cpp
src/antlr/TokenStreamBasicFilter.cpp
src/antlr/TokenStreamHiddenTokenFilter.cpp
src/antlr/TokenStreamRewriteEngine.cpp
src/antlr/TokenStreamSelector.cpp
src/antlr/TreeParser.cpp
src/brains/BlockSnapshotManager.cpp
src/brains/DataStorage.cpp
src/brains/MoleculeCreator.cpp
src/brains/PairList.cpp
src/brains/Register.cpp
src/brains/SimSnapshotManager.cpp
src/brains/Snapshot.cpp
src/brains/Stats.cpp
src/constraints/Rattle.cpp
src/hydrodynamics/Ellipsoid.cpp
src/hydrodynamics/HydroProp.cpp
src/hydrodynamics/Sphere.cpp
src/integrators/DLM.cpp
src/integrators/Integrator.cpp
src/integrators/IntegratorFactory.cpp
src/integrators/LangevinDynamics.cpp
src/integrators/LDForceManager.cpp
src/integrators/NgammaT.cpp
src/integrators/NPAT.cpp
src/integrators/NPrT.cpp
src/integrators/NPT.cpp
src/integrators/NPTf.cpp
src/integrators/NPTi.cpp
src/integrators/NPTsz.cpp
src/integrators/NPTxyz.cpp
src/integrators/NVE.cpp
src/integrators/NVT.cpp
src/integrators/VelocityVerletIntegrator.cpp
src/io/AtomTypesSectionParser.cpp
src/io/BaseAtomTypesSectionParser.cpp
src/io/BendTypesSectionParser.cpp
src/io/BondTypesSectionParser.cpp
src/io/ChargeAtomTypesSectionParser.cpp
src/io/DirectionalAtomTypesSectionParser.cpp
src/io/EAMAtomTypesSectionParser.cpp
src/io/ForceFieldOptions.cpp
src/io/GayBerneAtomTypesSectionParser.cpp
src/io/Globals.cpp
src/io/gzstream.cpp
src/io/InversionTypesSectionParser.cpp
src/io/LennardJonesAtomTypesSectionParser.cpp
src/io/MultipoleAtomTypesSectionParser.cpp
src/io/NonBondedInteractionsSectionParser.cpp
src/io/OptionSectionParser.cpp
src/io/ParamConstraint.cpp
src/io/SCAtomTypesSectionParser.cpp
src/io/SectionParser.cpp
src/io/SectionParserManager.cpp
src/io/ShapeAtomTypesSectionParser.cpp
src/io/StickyAtomTypesSectionParser.cpp
src/io/StickyPowerAtomTypesSectionParser.cpp
src/io/TorsionTypesSectionParser.cpp
src/io/ZConsReader.cpp
src/lattice/CubicLattice.cpp
src/lattice/FCCLattice.cpp
src/lattice/Lattice.cpp
src/lattice/LatticeFactory.cpp
src/lattice/shapedLattice.cpp
src/math/ChebyshevT.cpp
src/math/ChebyshevU.cpp
src/math/CubicSpline.cpp
src/math/LegendrePolynomial.cpp
src/math/RealSphericalHarmonic.cpp
src/math/RMSD.cpp
src/math/SeqRandNumGen.cpp
src/math/SphericalHarmonic.cpp
src/math/Wigner3jm.cpp
src/mdParser/FilenameObserver.cpp
src/mdParser/MDLexer.cpp
src/mdParser/MDParser.cpp
src/mdParser/MDTreeParser.cpp
src/minimizers/MinimizerFactory.cpp
src/minimizers/MinimizerParameterSet.cpp
src/nonbonded/EAM.cpp
src/nonbonded/Electrostatic.cpp
src/nonbonded/GB.cpp
src/nonbonded/InteractionManager.cpp
src/nonbonded/LJ.cpp
src/nonbonded/MAW.cpp
src/nonbonded/Morse.cpp
src/nonbonded/RepulsivePower.cpp
src/nonbonded/SC.cpp
src/nonbonded/Sticky.cpp
src/nonbonded/SwitchingFunction.cpp
src/primitives/Atom.cpp
src/primitives/Bend.cpp
src/primitives/DirectionalAtom.cpp
src/primitives/GhostBend.cpp
src/primitives/GhostTorsion.cpp
src/primitives/Inversion.cpp
src/primitives/Molecule.cpp
src/primitives/RigidBody.cpp
src/primitives/StuntDouble.cpp
src/primitives/Torsion.cpp
src/primitives/UreyBradleyBend.cpp
src/restraints/MolecularRestraint.cpp
src/restraints/ObjectRestraint.cpp
src/selection/DistanceFinder.cpp
src/selection/HullFinder.cpp
src/selection/IndexFinder.cpp
src/selection/NameFinder.cpp
src/selection/SelectionCompiler.cpp
src/selection/SelectionEvaluator.cpp
src/selection/SelectionManager.cpp
src/selection/SelectionToken.cpp
src/selection/TokenMap.cpp
src/types/AtomStamp.cpp
src/types/AtomType.cpp
src/types/BendStamp.cpp
src/types/BondStamp.cpp
src/types/CharmmTorsionType.cpp
src/types/Component.cpp
src/types/CutoffGroupStamp.cpp
src/types/DirectionalAtomType.cpp
src/types/FragmentStamp.cpp
src/types/ImproperCosineInversionType.cpp
src/types/InversionStamp.cpp
src/types/MoleculeStamp.cpp
src/types/NonBondedInteractionType.cpp
src/types/RestraintStamp.cpp
src/types/RigidBodyStamp.cpp
src/types/ShapeAtomType.cpp
src/types/TorsionStamp.cpp
src/types/ZconsStamp.cpp
src/UseTheForce/Amber_FF.cpp
src/UseTheForce/CLAYFF.cpp
src/UseTheForce/DUFF.cpp
src/UseTheForce/EADM_FF.cpp
src/UseTheForce/EAM_FF.cpp
src/UseTheForce/ForceField.cpp
src/UseTheForce/ForceFieldFactory.cpp
src/UseTheForce/MnM_FF.cpp
src/UseTheForce/SC_FF.cpp
src/UseTheForce/SHAPES_FF.cpp
src/utils/ElementsTable.cpp
src/utils/MoLocator.cpp
src/utils/OpenMDBitSet.cpp
src/utils/PropertyMap.cpp
src/utils/StringTokenizer.cpp
src/utils/StringUtils.cpp
src/utils/Trim.cpp
src/utils/Utility.cpp
src/utils/wildcards.cpp
src/visitors/AtomNameVisitor.cpp
src/visitors/AtomVisitor.cpp
src/visitors/CompositeVisitor.cpp
src/visitors/LipidTransVisitor.cpp
src/visitors/OtherVisitor.cpp
src/visitors/ReplacementVisitor.cpp
src/visitors/RigidBodyVisitor.cpp
src/visitors/ZconsVisitor.cpp
)

set( PARALLEL_SOURCE
src/brains/ForceManager.cpp
src/brains/SimCreator.cpp
src/brains/SimInfo.cpp
src/brains/Thermo.cpp
src/constraints/ZconstraintForceManager.cpp
src/integrators/RNEMD.cpp
src/integrators/Velocitizer.cpp
src/io/DumpReader.cpp
src/io/DumpWriter.cpp
src/io/RestReader.cpp
src/io/RestWriter.cpp
src/io/StatWriter.cpp
src/io/ZConsWriter.cpp
src/io/ifstrstream.cpp
src/math/ParallelRandNumGen.cpp
src/minimizers/CGFamilyMinimizer.cpp
src/minimizers/Minimizer.cpp
src/minimizers/PRCG.cpp
src/minimizers/SDMinimizer.cpp
src/parallel/ForceDecomposition.cpp
src/parallel/ForceMatrixDecomposition.cpp
src/restraints/RestraintForceManager.cpp
src/restraints/ThermoIntegrationForceManager.cpp
src/utils/ProgressBar.cpp
src/utils/simError.c
)

IF(QHULL_FOUND)
set(QHULL_SOURCE 
src/integrators/LangevinHullDynamics.cpp 
src/math/Triangle.cpp
)
set(QHULL_PARALLEL_SOURCE 
src/integrators/LangevinHullForceManager.cpp 
src/math/ConvexHull.cpp 
src/math/AlphaHull.cpp 
)
ENDIF(QHULL_FOUND)

add_library(openmd_core STATIC ${SOURCE} ${QHULL_SOURCE} )
add_library(openmd_single STATIC ${PARALLEL_SOURCE} ${QHULL_PARALLEL_SOURCE} )

IF(MPI_FOUND)
add_library(openmd_parallel STATIC ${PARALLEL_SOURCE} ${QHULL_PARALLEL_SOURCE} )
set_target_properties(openmd_parallel PROPERTIES
COMPILE_DEFINITIONS IS_MPI
)
ENDIF(MPI_FOUND)

add_executable(openmd src/applications/openmd/openmd.cpp)
target_link_libraries(openmd openmd_single openmd_core openmd_single openmd_core)

if (MPI_FOUND)
add_executable(openmd_MPI src/applications/openmd/openmd.cpp)
set_target_properties(openmd_MPI PROPERTIES
COMPILE_DEFINITIONS IS_MPI
)
target_link_libraries(openmd_MPI openmd_parallel openmd_core openmd_parallel openmd_core)
  INSTALL(TARGETS 
        openmd_parallel openmd_MPI
  RUNTIME DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  LIBRARY DESTINATION lib PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ 
  ARCHIVE DESTINATION lib PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ 
)
ENDIF (MPI_FOUND)

set (DUMP2XYZSOURCE 
src/applications/dump2Xyz/Dump2XYZ.cpp
src/applications/dump2Xyz/Dump2XYZCmd.c
)

set (DYNAMICPROPSSOURCE 
src/applications/dynamicProps/ActionCorrFunc.cpp
src/applications/dynamicProps/CrossTimeCorrFunc.cpp
src/applications/dynamicProps/DipoleCorrFunc.cpp
src/applications/dynamicProps/DirectionalRCorrFunc.cpp
src/applications/dynamicProps/DynamicProps.cpp
src/applications/dynamicProps/EnergyCorrFunc.cpp
src/applications/dynamicProps/FrameTimeCorrFunc.cpp
src/applications/dynamicProps/LegendreCorrFunc.cpp
src/applications/dynamicProps/MomentumCorrFunc.cpp
src/applications/dynamicProps/ParticleTimeCorrFunc.cpp
src/applications/dynamicProps/RadialRCorrFunc.cpp
src/applications/dynamicProps/RCorrFunc.cpp
src/applications/dynamicProps/StressCorrFunc.cpp
src/applications/dynamicProps/SystemDipoleCorrFunc.cpp
src/applications/dynamicProps/ThetaCorrFunc.cpp
src/applications/dynamicProps/TimeCorrFunc.cpp
src/applications/dynamicProps/VCorrFunc.cpp
src/applications/dynamicProps/DynamicPropsCmd.c
)

set (HYDROSOURCE 
src/applications/hydrodynamics/AnalyticalModel.cpp
src/applications/hydrodynamics/ApproximationModel.cpp
src/applications/hydrodynamics/BeadModel.cpp
src/applications/hydrodynamics/CompositeShape.cpp
src/applications/hydrodynamics/Hydro.cpp
src/applications/hydrodynamics/HydrodynamicsModel.cpp
src/applications/hydrodynamics/HydrodynamicsModelFactory.cpp
src/applications/hydrodynamics/RoughShell.cpp
src/applications/hydrodynamics/ShapeBuilder.cpp
src/applications/hydrodynamics/HydroCmd.c
)

set (STATICPROPSSOURCE 
src/applications/staticProps/AngleR.cpp
src/applications/staticProps/BondAngleDistribution.cpp
src/applications/staticProps/BondOrderParameter.cpp
src/applications/staticProps/BOPofR.cpp
src/applications/staticProps/DensityPlot.cpp
src/applications/staticProps/GofAngle2.cpp
src/applications/staticProps/GofR.cpp
src/applications/staticProps/GofRAngle.cpp
src/applications/staticProps/GofRZ.cpp
src/applications/staticProps/GofXyz.cpp
src/applications/staticProps/GofZ.cpp
src/applications/staticProps/Hxy.cpp
src/applications/staticProps/NanoLength.cpp
src/applications/staticProps/NanoVolume.cpp
src/applications/staticProps/ObjectCount.cpp
src/applications/staticProps/P2OrderParameter.cpp
src/applications/staticProps/pAngle.cpp
src/applications/staticProps/RadialDistrFunc.cpp
src/applications/staticProps/RhoR.cpp
src/applications/staticProps/RhoZ.cpp
src/applications/staticProps/RippleOP.cpp
src/applications/staticProps/SCDOrderParameter.cpp
src/applications/staticProps/StaticProps.cpp
src/applications/staticProps/TetrahedralityParam.cpp
src/applications/staticProps/TwoDGofR.cpp
src/applications/staticProps/StaticPropsCmd.c
)

set (NANOPARTICLEBUILDERSOURCE 
src/applications/nanoparticleBuilder/nanoparticleBuilder.cpp
src/applications/nanoparticleBuilder/shapedLatticeSpherical.cpp
src/applications/nanoparticleBuilder/nanoparticleBuilderCmd.c
)

set (RANDOMBUILDERSOURCE
src/applications/randomBuilder/randomBuilder.cpp
src/applications/randomBuilder/randomBuilderCmd.c
)

set(SIMPLEBUILDERSOURCE 
src/applications/simpleBuilder/simpleBuilder.cpp
src/applications/simpleBuilder/simpleBuilderCmd.c
)

set(THERMALIZERSOURCE 
src/applications/thermalizer/thermalizer.cpp
src/applications/thermalizer/thermalizerCmd.c 
)

add_executable(Dump2XYZ ${DUMP2XYZSOURCE})
target_link_libraries(Dump2XYZ openmd_single openmd_core openmd_single openmd_core)
add_executable(DynamicProps ${DYNAMICPROPSSOURCE})
target_link_libraries(DynamicProps openmd_single openmd_core openmd_single openmd_core)
add_executable(Hydro ${HYDROSOURCE})
target_link_libraries(Hydro openmd_single openmd_core openmd_single openmd_core)
add_executable(StaticProps ${STATICPROPSSOURCE})
target_link_libraries(StaticProps openmd_single openmd_core openmd_single openmd_core)
add_executable(nanoparticleBuilder ${NANOPARTICLEBUILDERSOURCE})
target_link_libraries(nanoparticleBuilder openmd_single openmd_core openmd_single openmd_core)
add_executable(randomBuilder ${RANDOMBUILDERSOURCE})
target_link_libraries(randomBuilder openmd_single openmd_core openmd_single openmd_core)
add_executable(simpleBuilder ${SIMPLEBUILDERSOURCE})
target_link_libraries(simpleBuilder openmd_single openmd_core openmd_single openmd_core)
add_executable(thermalizer ${THERMALIZERSOURCE})
target_link_libraries(thermalizer openmd_single openmd_core openmd_single openmd_core)

if (OPENBABEL2_FOUND)
set (ATOM2MDSOURCE 
src/applications/atom2md/atom2md.cpp 
src/applications/atom2md/openmdformat.cpp
)
  add_executable(atom2md ${ATOM2MDSOURCE})
  target_link_libraries(atom2md openmd_single openmd_core openmd_single openmd_core ${OPENBABEL2_LIBRARIES})
  INSTALL(TARGETS atom2md RUNTIME DESTINATION bin
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
ENDIF (OPENBABEL2_FOUND)


set(PY_FILES
src/applications/hydrodynamics/diffExplainer
src/applications/utilities/affineScale
src/applications/utilities/dumpConverter
src/applications/utilities/md-solvator
src/applications/utilities/md2md
src/applications/utilities/mdSplit
src/applications/utilities/principalAxisCalculator
src/applications/utilities/stat2visco
src/applications/utilities/waterRotator
)

IF(PYTHON_EXECUTABLE)
   foreach(PY_FILE ${PY_FILES})
      GET_FILENAME_COMPONENT(filename "${PY_FILE}" NAME)
      CONFIGURE_FILE(${PY_FILE} 
                     "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${filename}" @ONLY)
      INSTALL(FILES "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${filename}"
              DESTINATION bin 
              PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
   endforeach(PY_FILE)
ENDIF(PYTHON_EXECUTABLE)

set(PERL_FILES
src/applications/utilities/solvator
src/applications/utilities/waterBoxer
)

IF(PERL_FOUND)
   foreach(PERL_FILE ${PERL_FILES})
      GET_FILENAME_COMPONENT(filename "${PERL_FILE}" NAME)
      configure_file(${PERL_FILE} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${filename}" @ONLY)
      INSTALL(FILES "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${filename}"
              DESTINATION bin 
              PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

   endforeach(PERL_FILE)
ENDIF(PERL_FOUND)

INSTALL(FILES AUTHORS LICENSE README DESTINATION . 
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(TARGETS 
        openmd_core 
        openmd_single
        openmd
        Dump2XYZ 
        simpleBuilder 
        StaticProps 
        DynamicProps
        randomBuilder 
        nanoparticleBuilder 
        thermalizer 
        Hydro
  RUNTIME DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  LIBRARY DESTINATION lib PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ 
  ARCHIVE DESTINATION lib PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ 
)

install(DIRECTORY forceFields/   
        DESTINATION forceFields 
        PATTERN ".svn" EXCLUDE 
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
install(DIRECTORY samples/       
        DESTINATION samples 
        PATTERN ".svn" EXCLUDE 
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
