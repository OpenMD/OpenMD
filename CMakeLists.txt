# Top level cmake script for OpenMD.
project(OpenMD CXX)
cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

IF(WIN32)
   ADD_DEFINITIONS(/D_CRT_SECURE_NO_WARNINGS)
ENDIF(WIN32)

if (NOT NO_PARALLEL)
    if(WIN32)
        #This is for finding MS-MPI.
        set(MSMPI_BIN $ENV{MSMPI_BIN})
        set(MSMPI_INC $ENV{MSMPI_INC})
        set(MSMPI_LIB32 $ENV{MSMPI_LIB32})
        set(MSMPI_LIB64 $ENV{MSMPI_LIB64})
        message(STATUS "MSMPI_BIN=${MSMPI_BIN}")
        message(STATUS "MSMPI_INC=${MSMPI_INC}")
        message(STATUS "MSMPI_LIB32=${MSMPI_LIB32}")
        message(STATUS "MSMPI_LIB64=${MSMPI_LIB64}")
        IF( "${MSMPI_BIN}" STREQUAL "" OR
                "${MSMPI_INC}" STREQUAL "" OR
                "${MSMPI_LIB32}" STREQUAL "" OR
                "${MSMPI_LIB64}" STREQUAL "")
            message(FATAL_ERROR "Can not find MSMPI install. Please install MSMPI or check the ENV variable:
            MSMPI_INC, MSMPI_LIB32, MSMPI_LIB64, MSMPI_BIN")
        endif()
        set(MPIEXEC "${MSMPI_BIN}/mpiexec.exe")
        # For building MPI programs the selected Visual Studio compiler is used, namely cl.exe.
        # So there is no need to set a specific MPI compiler.
        #set(MPI_CXX_COMPILER "${CMAKE_CXX_COMPILER}")
        set(MPI_CXX_INCLUDE_PATH "${MSMPI_INC}")
        # Make sure the correct libraries (64-bit or 32-bit) are selected.
        # Decide between 32-bit and 64-bit libraries for Microsoft's MPI
        if("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
            set(MS_MPI_ARCH_LIB ${MSMPI_LIB64})
        else()
            set(MS_MPI_ARCH_LIB ${MSMPI_LIB32})
        endif()
        set(MPI_CXX_LIBRARIES "${MS_MPI_ARCH_LIB}/msmpi.lib")
        set(MPI_C_INCLUDE_PATH "${MPI_CXX_INCLUDE_PATH}")
        set(MPI_C_LIBRARIES "{${MPI_CXX_LIBRARIES}")
        set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "Flag used by MPI to specify the number of processes for MPIEXEC;
        the next option will be the number of processes.")
        set(MPI_CXX_FOUND ON)
        set(MPI_FOUND ON)
    else()
        find_package(MPI)
    endif()
    if (NOT MPI_CXX_FOUND)
      message(STATUS "========== OpenMD Parallel Information ==========")
      message(STATUS )
      message(STATUS "No MPI compiler found. Perhaps you want to set one explicitly?")
      message(STATUS "To override the default compiler, set the environment variable")
      message(STATUS "     export CXX=/full/path/to/mpic++")
      message(STATUS "in bash or ksh or sh.  In csh or tcsh, use:")
      message(STATUS "     setenv CXX /full/path/to/mpic++")
      message(STATUS "before running the cmake command.")
      message(STATUS )
      message(STATUS "=================================================")
   endif()
   IF ( MPI_CXX_FOUND )
      INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_EXE_LINKER_FLAGS}")
      if(DEFINED $ENV{CXX})
        message(STATUS "Found MPI compiler, but using the CXX compiler set in the CXX environment variable, $ENV{CXX}")
      else()
        set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
        set(ENV{CXX} ${MPI_CXX_COMPILER})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
      endif()
   ENDIF()
endif()

include(CheckCXXCompilerFlag)
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)

if ( CMAKE_COMPILER_IS_GNUCXX )
    add_definitions("-Wall")
endif ( CMAKE_COMPILER_IS_GNUCXX )
if ( MSVC )
    add_definitions("/W3 /DNOMINMAX")
endif ( MSVC )

# OpenMD version number.
set (VERSION_MAJOR "2")
set (VERSION_MINOR "7")
set (VERSION_TINY "0")
option(SINGLE_PRECISION "Build Single precision (float) version" OFF)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if (CMAKE_HOST_UNIX)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/openmd"
      CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
  else (CMAKE_HOST_UNIX)
    GetDefaultWindowsPrefixBase(CMAKE_GENERIC_PROGRAM_FILES)
    set (CMAKE_INSTALL_PREFIX
      "${CMAKE_GENERIC_PROGRAM_FILES}/OpenMD"
      CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
    set (CMAKE_GENERIC_PROGRAM_FILES)
  endif (CMAKE_HOST_UNIX)
endif (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

ADD_DEFINITIONS( -DFRC_PATH=${CMAKE_INSTALL_PREFIX}/forceFields )

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

find_package(Git)
if(GIT_FOUND)
  include(GetGitRevisionDescription)
  get_git_head_revision(GIT_REFSPEC GIT_SHA1)
  if (${GIT_SHA1} MATCHES "GITDIR-NOTFOUND")
    set(GIT_SHA1 Release)
  endif()
  MESSAGE(STATUS "Current revision is ${GIT_SHA1}")
  macro(make_release_tarball targetname releasename)
    add_custom_target(${targetname} COMMAND ${GIT_EXECUTABLE} archive --prefix="${releasename}/" --format=tar HEAD --output "${releasename}.tar")
  endmacro()
  make_release_tarball("openmd_release" "openmd-${VERSION_MAJOR}.${VERSION_MINOR}")
else()
  SET(GIT_SHA1 Release)
endif()

ADD_DEFINITIONS(-DCMDLINE_PARSER_VERSION="\(OpenMD\ ${VERSION_MAJOR}.${VERSION_MINOR}\ Revision\ ${GIT_SHA1}\)")

macro(LONG_TODAY RESULT)
  set(MONTHS ";Jan;Feb;Mar;Apr;May;Jun;Jul;Aug;Sep;Oct;Nov;Dec")
  set(DAYS_OF_WEEK "Sun;Mon;Tue;Wed;Thu;Fri;Sat")

  string(TIMESTAMP DOW "%w")
  list(GET DAYS_OF_WEEK "${DOW}" DOW)
  string(TIMESTAMP M "%m")
  list(GET MONTHS "${M}" MON)
  string(TIMESTAMP D "%d")
  string(TIMESTAMP TIME "%H:%M:%S")
  string(TIMESTAMP Y "%Y")
  set(${RESULT} "${DOW} ${MON} ${D} ${TIME} ${Y}")
endmacro()

LONG_TODAY(BUILD_DATE)

check_include_file_cxx(conio.h      HAVE_CONIO_H)
check_cxx_symbol_exists(strncasecmp   "string.h"   HAVE_STRNCASECMP)

# Optional libraries: If we can find these, we will build with them
# Look for OpenBabel libraries
find_package(OpenBabel3)
IF(OPENBABEL3_FOUND)
  SET(USE_OPENBABEL true)
  set(OPENBABEL_INCLUDE_DIR ${OPENBABEL3_INCLUDE_DIR})
  set(OPENBABEL_LIBRARIES ${OPENBABEL3_LIBRARIES})
  include_directories(${OPENBABEL_INCLUDE_DIR})
ELSE(OPENBABEL3_FOUND)
  find_package(OpenBabel2)
  IF(OPENBABEL2_FOUND)
     SET(USE_OPENBABEL true)
     set(OPENBABEL_INCLUDE_DIR ${OPENBABEL2_INCLUDE_DIR})
     set(OPENBABEL_LIBRARIES ${OPENBABEL2_LIBRARIES})
     include_directories(${OPENBABEL_INCLUDE_DIR})
  ELSE(OPENBABEL2_FOUND)
     SET(USE_OPENBABEL false)
     MESSAGE(STATUS "No OpenBabel found - will not build atom2omd")
  ENDIF(OPENBABEL2_FOUND)
ENDIF(OPENBABEL3_FOUND)


#Look for QHULL Libraries
find_package(Qhull)
IF(QHULL_FOUND)
  SET(HAVE_QHULL 1)
  include_directories(${QHULL_INCLUDE_DIR})
  LINK_LIBRARIES(${QHULL_LIBRARIES})
ELSE(QHULL_FOUND)
  MESSAGE(STATUS "No Qhull found - will be missing some features")
ENDIF(QHULL_FOUND)

# zlib stuff
find_package(ZLIB)
if(ZLIB_FOUND)
  SET(HAVE_LIBZ 1)
  add_definitions(-DHAVE_LIBZ)
  include_directories(${ZLIB_INCLUDE_DIR})
  LINK_LIBRARIES(${ZLIB_LIBRARIES})
ELSE(ZLIB_FOUND)
  MESSAGE(STATUS "No zlib found - will be missing compressed dump files")
endif(ZLIB_FOUND)

#FFTW3
IF(SINGLE_PRECISION)
  find_package(FFTW3 COMPONENTS single)
else()
  find_package(FFTW3 COMPONENTS double)
endif(SINGLE_PRECISION)
if (FFTW3_FOUND)
  include_directories (${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} ${FFTW3_INCLUDE_DIR})
  SET(HAVE_LIBFFTW 1)
  SET(HAVE_FFTW3_H 1)
  INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})
  LINK_LIBRARIES(${FFTW3_LIBRARIES})
ELSE(FFTW3_FOUND)
  MESSAGE(STATUS "No fftw3 found - will be missing some analysis modules")
endif (FFTW3_FOUND)


# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in" "${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile" @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doc"
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif(DOXYGEN_FOUND)

FIND_PACKAGE(Perl)
IF(PERL_FOUND)
  SET(PERL ${PERL_EXECUTABLE})
ELSE(PERL_FOUND)
  MESSAGE(STATUS "Failed to find perl - some scripts will not be made")
ENDIF(PERL_FOUND)

SET(PERL_INSTALLDIRS "site" CACHE STRING "Perl installation
locations")

FIND_PACKAGE(Python3 COMPONENTS Interpreter)
IF(Python3_Interpreter_FOUND)
  SET(PYTHON ${Python3_EXECUTABLE})
ELSE(Python3_Interpreter_FOUND)
  MESSAGE(STATUS "Failed to find python - some scripts will not be made")
ENDIF(Python3_Interpreter_FOUND)

FIND_PACKAGE(NumPy)
IF(NOT NUMPY_FOUND)
  MESSAGE(STATUS "Failed to find NumPy - some scripts will not be made")
ENDIF(NOT NUMPY_FOUND)

FIND_PACKAGE(SciPy)
IF(NOT SCIPY_FOUND)
  MESSAGE(STATUS "Failed to find SciPy - some scripts will not be made")
ENDIF(NOT SCIPY_FOUND)


configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.cmake"
  "${PROJECT_BINARY_DIR}/config.h"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/src/utils/Revision.cpp.in"
  "${PROJECT_BINARY_DIR}/Revision.cpp"
)

include_directories("${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/src")

configure_file (
  "${PROJECT_SOURCE_DIR}/cmake/openmd.spec.in"
  "${PROJECT_SOURCE_DIR}/openmd.spec"
  @ONLY
  )

#Add executables for build
set (PROGRAMS
openmd
openmd_MPI
Dump2XYZ
StaticProps
DynamicProps
SequentialProps
simpleBuilder
randomBuilder
nanoparticleBuilder
nanorodBuilder
nanorod_pentBuilder
omd2omd
equationofstate
thermalizer
recenter
elasticConstants
Hydro)

set(SOURCE
"${PROJECT_BINARY_DIR}/Revision.cpp"
src/antlr/ANTLRUtil.cpp
src/antlr/ASTFactory.cpp
src/antlr/ASTNULLType.cpp
src/antlr/ASTRefCount.cpp
src/antlr/BaseAST.cpp
src/antlr/BitSet.cpp
src/antlr/CommonAST.cpp
src/antlr/CommonASTWithHiddenTokens.cpp
src/antlr/CommonHiddenStreamToken.cpp
src/antlr/CommonToken.cpp
src/antlr/LLkParser.cpp
src/antlr/MismatchedCharException.cpp
src/antlr/MismatchedTokenException.cpp
src/antlr/NoViableAltException.cpp
src/antlr/NoViableAltForCharException.cpp
src/antlr/Parser.cpp
src/antlr/RecognitionException.cpp
src/antlr/String.cpp
src/antlr/Token.cpp
src/antlr/TokenRefCount.cpp
src/antlr/TokenStreamBasicFilter.cpp
src/antlr/TokenStreamHiddenTokenFilter.cpp
src/antlr/TokenStreamRewriteEngine.cpp
src/antlr/TokenStreamSelector.cpp
src/antlr/TreeParser.cpp
src/brains/BlockSnapshotManager.cpp
src/brains/DataStorage.cpp
src/brains/ForceField.cpp
src/brains/MoleculeCreator.cpp
src/brains/PairList.cpp
src/brains/Register.cpp
src/brains/SimSnapshotManager.cpp
src/brains/Snapshot.cpp
src/brains/Stats.cpp
src/hydrodynamics/Ellipsoid.cpp
src/hydrodynamics/HydroProp.cpp
src/hydrodynamics/Sphere.cpp
src/integrators/DLM.cpp
src/flucq/FluctuatingChargeDamped.cpp
src/flucq/FluctuatingChargeLangevin.cpp
src/flucq/FluctuatingChargeParameters.cpp
src/flucq/FluctuatingChargeNVE.cpp
src/flucq/FluctuatingChargeNVT.cpp
src/integrators/Integrator.cpp
src/integrators/IntegratorFactory.cpp
src/integrators/LangevinDynamics.cpp
src/integrators/LDForceManager.cpp
src/integrators/NgammaT.cpp
src/integrators/NPA.cpp
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
src/io/FluctuatingChargeAtomTypesSectionParser.cpp
src/io/ForceFieldOptions.cpp
src/io/GayBerneAtomTypesSectionParser.cpp
src/io/Globals.cpp
src/io/InversionTypesSectionParser.cpp
src/io/LennardJonesAtomTypesSectionParser.cpp
src/io/MultipoleAtomTypesSectionParser.cpp
src/io/NonBondedInteractionsSectionParser.cpp
src/io/OptionSectionParser.cpp
src/io/ParamConstraint.cpp
src/io/PolarizableAtomTypesSectionParser.cpp
src/io/SCAtomTypesSectionParser.cpp
src/io/SectionParser.cpp
src/io/SectionParserManager.cpp
src/io/ShapeAtomTypesSectionParser.cpp
src/io/StickyAtomTypesSectionParser.cpp
src/io/StickyPowerAtomTypesSectionParser.cpp
src/io/TorsionTypesSectionParser.cpp
src/io/UFFAtomTypesSectionParser.cpp
src/io/ZConsReader.cpp
src/lattice/CubicLattice.cpp
src/lattice/FCCLattice.cpp
src/lattice/SCLattice.cpp
src/lattice/BCCLattice.cpp
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
src/optimization/OptimizationFactory.cpp
src/optimization/Armijo.cpp
src/optimization/BFGS.cpp
src/optimization/BoxObjectiveFunction.cpp
src/optimization/ConjugateGradient.cpp
src/optimization/Constraint.cpp
src/optimization/EndCriteria.cpp
src/optimization/LineSearch.cpp
src/optimization/LineSearchBasedMethod.cpp
src/optimization/SteepestDescent.cpp
src/optimization/MinimizerParameters.cpp
src/optimization/Problem.cpp
src/nonbonded/Buckingham.cpp
src/nonbonded/EAM.cpp
src/nonbonded/GB.cpp
src/nonbonded/InteractionManager.cpp
src/nonbonded/InversePowerSeries.cpp
src/nonbonded/LJ.cpp
src/nonbonded/MAW.cpp
src/nonbonded/Morse.cpp
src/nonbonded/RepulsivePower.cpp
src/nonbonded/Mie.cpp
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
src/primitives/ShortRangeInteraction.cpp
src/primitives/Torsion.cpp
src/primitives/UreyBradleyBend.cpp
src/restraints/MolecularRestraint.cpp
src/restraints/ObjectRestraint.cpp
src/selection/HullFinder.cpp
src/selection/IndexFinder.cpp
src/selection/NameFinder.cpp
src/selection/SelectionCompiler.cpp
src/selection/SelectionEvaluator.cpp
src/selection/SelectionToken.cpp
src/selection/TokenMap.cpp
src/types/AtomStamp.cpp
src/types/AtomType.cpp
src/types/BendStamp.cpp
src/types/BendTypeParser.cpp
src/types/BondStamp.cpp
src/types/BondTypeParser.cpp
src/types/CharmmTorsionType.cpp
src/types/Component.cpp
src/types/ConstraintStamp.cpp
src/types/CutoffGroupStamp.cpp
src/types/DirectionalAdapter.cpp
src/types/EAMAdapter.cpp
src/types/FixedChargeAdapter.cpp
src/types/FluctuatingChargeAdapter.cpp
src/types/FragmentStamp.cpp
src/types/GayBerneAdapter.cpp
src/types/ImproperCosineInversionType.cpp
src/types/InversionStamp.cpp
src/types/InversionTypeParser.cpp
src/types/LennardJonesAdapter.cpp
src/types/MoleculeStamp.cpp
src/types/MultipoleAdapter.cpp
src/types/NonBondedInteractionType.cpp
src/types/PolarizableAdapter.cpp
src/types/RestraintStamp.cpp
src/types/RigidBodyStamp.cpp
src/types/ShapeAtomType.cpp
src/types/StickyAdapter.cpp
src/types/SuttonChenAdapter.cpp
src/types/TorsionStamp.cpp
src/types/TorsionTypeParser.cpp
src/types/UFFAdapter.cpp
src/types/ZconsStamp.cpp
src/utils/ElementsTable.cpp
src/utils/MoLocator.cpp
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
src/rnemd/RNEMDParameters.cpp
src/clusters/Icosahedron.cpp
src/clusters/Cuboctahedron.cpp
src/clusters/Decahedron.cpp
)

set( PARALLEL_SOURCE
src/antlr/CharBuffer.cpp
src/antlr/CharScanner.cpp
src/antlr/InputBuffer.cpp
src/antlr/TokenBuffer.cpp
src/mdParser/MDLexer.cpp
src/mdParser/MDParser.cpp
src/mdParser/MDTreeParser.cpp
src/brains/ForceManager.cpp
src/brains/SimCreator.cpp
src/brains/SimInfo.cpp
src/brains/Thermo.cpp
src/brains/Velocitizer.cpp
src/constraints/ZconstraintForceManager.cpp
src/constraints/Rattle.cpp
src/constraints/Shake.cpp
src/flucq/FluctuatingChargeConstraints.cpp
src/flucq/FluctuatingChargeObjectiveFunction.cpp
src/flucq/FluctuatingChargeForces.cpp
src/flucq/FluctuatingChargePropagator.cpp
src/integrators/LangevinHullForceManager.cpp
src/rnemd/RNEMD.cpp
src/io/ConstraintWriter.cpp
src/io/DumpReader.cpp
src/io/DumpWriter.cpp
src/io/RestReader.cpp
src/io/RestWriter.cpp
src/io/StatWriter.cpp
src/io/ZConsWriter.cpp
src/io/ifstrstream.cpp
src/math/ParallelRandNumGen.cpp
src/nonbonded/Electrostatic.cpp
src/optimization/PotentialEnergyObjectiveFunction.cpp
src/parallel/ForceDecomposition.cpp
src/parallel/ForceMatrixDecomposition.cpp
src/perturbations/UniformField.cpp
src/perturbations/MagneticField.cpp
src/perturbations/UniformGradient.cpp
src/restraints/RestraintForceManager.cpp
src/restraints/ThermoIntegrationForceManager.cpp
src/selection/DistanceFinder.cpp
src/selection/SelectionManager.cpp
src/selection/SelectionSet.cpp
src/utils/ProgressBar.cpp
src/utils/simError.cpp
src/utils/OpenMDBitSet.cpp
)

IF(ZLIB_FOUND)
set(ZLIB_SOURCE
src/io/gzstream.cpp
)
ENDIF(ZLIB_FOUND)

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

IF(WIN32)
set(GETOPT_SOURCE
src/utils/wingetopt.cpp
)
ENDIF(WIN32)

add_library(openmd_core STATIC ${SOURCE} ${QHULL_SOURCE} ${ZLIB_SOURCE} )
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
src/applications/dump2Xyz/Dump2XYZCmd.cpp
)

set (OMD2OMDSOURCE
src/applications/omd2omd/omd2omd.cpp
src/applications/omd2omd/omd2omdCmd.cpp
)

set (EOSSOURCE
src/applications/eos/equationofstate.cpp
src/applications/eos/equationofstateCmd.cpp
)
set (DYNAMICPROPSSOURCE
src/applications/dynamicProps/ActionCorrFunc.cpp
src/applications/dynamicProps/CollectiveDipoleDisplacement.cpp
src/applications/dynamicProps/CurrentDensityAutoCorrFunc.cpp
src/applications/dynamicProps/DipoleCorrFunc.cpp
src/applications/dynamicProps/DirectionalRCorrFunc.cpp
src/applications/dynamicProps/Displacement.cpp
src/applications/dynamicProps/DynamicProps.cpp
src/applications/dynamicProps/FreqFlucCorrFunc.cpp
src/applications/dynamicProps/HBondJump.cpp
src/applications/dynamicProps/HBondPersistence.cpp
src/applications/dynamicProps/LegendreCorrFunc.cpp
src/applications/dynamicProps/LegendreCorrFuncZ.cpp
src/applications/dynamicProps/cOHz.cpp
src/applications/dynamicProps/MomAngMomCorrFunc.cpp
src/applications/dynamicProps/ForTorCorrFunc.cpp
src/applications/dynamicProps/ChargeKineticCorrFunc.cpp
src/applications/dynamicProps/ChargeOrientationCorrFunc.cpp
src/applications/dynamicProps/ForceAutoCorrFunc.cpp
src/applications/dynamicProps/TorForCorrFunc.cpp
src/applications/dynamicProps/TorqueAutoCorrFunc.cpp
src/applications/dynamicProps/RCorrFunc.cpp
src/applications/dynamicProps/SelectionCorrFunc.cpp
src/applications/dynamicProps/StressCorrFunc.cpp
src/applications/dynamicProps/SystemDipoleCorrFunc.cpp
src/applications/dynamicProps/ThetaCorrFunc.cpp
src/applications/dynamicProps/TimeCorrFunc.cpp
src/applications/dynamicProps/VCorrFunc.cpp
src/applications/dynamicProps/VelocityAutoOutProductCorrFunc.cpp
src/applications/dynamicProps/AngularVelocityAutoOutProductCorrFunc.cpp
src/applications/dynamicProps/VelAngularVelOutProdCorrFunc.cpp
src/applications/dynamicProps/AngularVelVelOutProdCorrFunc.cpp
src/applications/dynamicProps/WCorrFunc.cpp
src/applications/dynamicProps/BondCorrFunc.cpp
src/applications/dynamicProps/DynamicPropsCmd.cpp

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
src/applications/hydrodynamics/HydroCmd.cpp
)

set (SEQUENTIALPROPSSOURCE
src/applications/sequentialProps/CenterOfMass.cpp
src/applications/sequentialProps/equipartitionTest.cpp
src/applications/sequentialProps/COMVel.cpp
src/applications/sequentialProps/ContactAngle1.cpp
src/applications/sequentialProps/ContactAngle2.cpp
src/applications/sequentialProps/GCNSeq.cpp
src/applications/sequentialProps/SequentialAnalyzer.cpp
src/applications/sequentialProps/SequentialProps.cpp
src/applications/sequentialProps/SequentialPropsCmd.cpp
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
src/applications/staticProps/GofRAngle2.cpp
src/applications/staticProps/GofRZ.cpp
src/applications/staticProps/GofXyz.cpp
src/applications/staticProps/GofZ.cpp
src/applications/staticProps/CoordinationNumber.cpp
src/applications/staticProps/Hxy.cpp
src/applications/staticProps/HBondGeometric.cpp
src/applications/staticProps/Kirkwood.cpp
src/applications/staticProps/Field.cpp
src/applications/staticProps/MultipoleSum.cpp
src/applications/staticProps/NanoLength.cpp
src/applications/staticProps/NanoVolume.cpp
src/applications/staticProps/NitrileFrequencyMap.cpp
src/applications/staticProps/ObjectCount.cpp
src/applications/staticProps/P2OrderParameter.cpp
src/applications/staticProps/PipeDensity.cpp
src/applications/staticProps/PotDiff.cpp
src/applications/staticProps/pAngle.cpp
src/applications/staticProps/RadialDistrFunc.cpp
src/applications/staticProps/RhoR.cpp
src/applications/staticProps/RhoZ.cpp
src/applications/staticProps/RippleOP.cpp
src/applications/staticProps/RNEMDStats.cpp
src/applications/staticProps/SCDOrderParameter.cpp
src/applications/staticProps/SpatialStatistics.cpp
src/applications/staticProps/StaticAnalyser.cpp
src/applications/staticProps/StaticProps.cpp
src/applications/staticProps/SurfaceDiffusion.cpp
src/applications/staticProps/TetrahedralityHBMatrix.cpp
src/applications/staticProps/TetrahedralityParam.cpp
src/applications/staticProps/TetrahedralityParamDens.cpp
src/applications/staticProps/TetrahedralityParamZ.cpp
src/applications/staticProps/TetrahedralityParamXYZ.cpp
src/applications/staticProps/TwoDGofR.cpp
src/applications/staticProps/StaticPropsCmd.cpp
src/applications/staticProps/VelocityZ.cpp
src/applications/staticProps/DensityHistogram.cpp
src/applications/staticProps/MomentumHistogram.cpp
src/applications/staticProps/ChargeHistogram.cpp
src/applications/staticProps/CurrentDensity.cpp
src/applications/staticProps/ChargeZ.cpp
src/applications/staticProps/PositionZ.cpp
src/applications/staticProps/DipoleOrientation.cpp
src/applications/staticProps/ChargeDensityZ.cpp
src/applications/staticProps/OrderParameterProbZ.cpp

)

set (NANOPARTICLEBUILDERSOURCE
src/applications/nanoparticleBuilder/nanoparticleBuilder.cpp
src/applications/nanoparticleBuilder/shapedLatticeSpherical.cpp
src/applications/nanoparticleBuilder/nanoparticleBuilderCmd.cpp
)

set (NANORODBUILDERSOURCE
src/applications/nanoparticleBuilder/nanorodBuilder.cpp
src/applications/nanoparticleBuilder/shapedLatticeRod.cpp
src/applications/nanoparticleBuilder/shapedLatticeEllipsoid.cpp
src/applications/nanoparticleBuilder/nanorodBuilderCmd.cpp
)

set (NANOROD_PENTBUILDERSOURCE
src/applications/nanoparticleBuilder/nanorod_pentBuilder.cpp
src/applications/nanoparticleBuilder/shapedLatticePentRod.cpp
src/applications/nanoparticleBuilder/nanorod_pentBuilderCmd.cpp
)

set (ICOSAHEDRALBUILDERSOURCE
src/applications/nanoparticleBuilder/icosahedralBuilder.cpp
src/applications/nanoparticleBuilder/icosahedralBuilderCmd.cpp
)

set (RANDOMBUILDERSOURCE
src/applications/randomBuilder/randomBuilder.cpp
src/applications/randomBuilder/randomBuilderCmd.cpp
)

set(SIMPLEBUILDERSOURCE
src/applications/simpleBuilder/simpleBuilder.cpp
src/applications/simpleBuilder/simpleBuilderCmd.cpp
)

set(THERMALIZERSOURCE
src/applications/thermalizer/thermalizer.cpp
src/applications/thermalizer/thermalizerCmd.cpp
)

set(RECENTERSOURCE
src/applications/recenter/recenter.cpp
src/applications/recenter/recenterCmd.cpp
)

set(ELASTICCONSTANTSSOURCE
src/applications/elasticConstants/elasticConstants.cpp
src/applications/elasticConstants/elasticConstantsCmd.cpp
)

add_executable(Dump2XYZ ${DUMP2XYZSOURCE} ${GETOPT_SOURCE})
target_link_libraries(Dump2XYZ openmd_single openmd_core openmd_single openmd_core)
add_executable(DynamicProps ${DYNAMICPROPSSOURCE} ${GETOPT_SOURCE})
target_link_libraries(DynamicProps openmd_single openmd_core openmd_single openmd_core)
add_executable(Hydro ${HYDROSOURCE} ${GETOPT_SOURCE})
target_link_libraries(Hydro openmd_single openmd_core openmd_single openmd_core)
add_executable(omd2omd ${OMD2OMDSOURCE} ${GETOPT_SOURCE})
target_link_libraries(omd2omd openmd_single openmd_core openmd_single openmd_core)
add_executable(equationofstate ${EOSSOURCE} ${GETOPT_SOURCE})
target_link_libraries(equationofstate openmd_single openmd_core openmd_single openmd_core)
add_executable(StaticProps ${STATICPROPSSOURCE} ${GETOPT_SOURCE})
target_link_libraries(StaticProps openmd_single openmd_core openmd_single openmd_core)
add_executable(SequentialProps ${SEQUENTIALPROPSSOURCE} ${GETOPT_SOURCE})
target_link_libraries(SequentialProps openmd_single openmd_core openmd_single openmd_core)
add_executable(nanoparticleBuilder ${NANOPARTICLEBUILDERSOURCE} ${GETOPT_SOURCE})
target_link_libraries(nanoparticleBuilder openmd_single openmd_core openmd_single openmd_core openmd_single)
add_executable(nanorodBuilder ${NANORODBUILDERSOURCE} ${GETOPT_SOURCE})
target_link_libraries(nanorodBuilder openmd_single openmd_core openmd_single openmd_core openmd_single)
add_executable(nanorod_pentBuilder ${NANOROD_PENTBUILDERSOURCE} ${GETOPT_SOURCE})
target_link_libraries(nanorod_pentBuilder openmd_single openmd_core openmd_single openmd_core openmd_single)
add_executable(icosahedralBuilder ${ICOSAHEDRALBUILDERSOURCE} ${GETOPT_SOURCE})
target_link_libraries(icosahedralBuilder openmd_single openmd_core openmd_single openmd_core openmd_single)
add_executable(randomBuilder ${RANDOMBUILDERSOURCE} ${GETOPT_SOURCE})
target_link_libraries(randomBuilder openmd_single openmd_core openmd_single openmd_core openmd_single)
add_executable(simpleBuilder ${SIMPLEBUILDERSOURCE} ${GETOPT_SOURCE})
target_link_libraries(simpleBuilder openmd_single openmd_core openmd_single openmd_core openmd_single)
add_executable(thermalizer ${THERMALIZERSOURCE} ${GETOPT_SOURCE})
target_link_libraries(thermalizer openmd_single openmd_core openmd_single openmd_core openmd_single)
add_executable(recenter ${RECENTERSOURCE} ${GETOPT_SOURCE})
target_link_libraries(recenter openmd_single openmd_core openmd_single openmd_core openmd_single)
add_executable(elasticConstants ${ELASTICCONSTANTSSOURCE} ${GETOPT_SOURCE})
target_link_libraries(elasticConstants openmd_single openmd_core openmd_single openmd_core openmd_single)

if (USE_OPENBABEL)
set (ATOM2OMDSOURCE
src/applications/atom2omd/atom2omd.cpp
src/applications/atom2omd/openmdformat.cpp
)
  add_executable(atom2omd ${ATOM2OMDSOURCE} ${GETOPT_SOURCE})
  target_link_libraries(atom2omd openmd_single openmd_core openmd_single openmd_core ${OPENBABEL_LIBRARIES})
  INSTALL(TARGETS atom2omd RUNTIME DESTINATION bin
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
ENDIF (USE_OPENBABEL)


set(PYNONUMPY_FILES
src/applications/utilities/affineScale
src/applications/utilities/dumpConverter
src/applications/utilities/funcflExtractor
src/applications/utilities/omd-solvator
src/applications/utilities/omdCombined
src/applications/utilities/omdSplit
src/applications/utilities/omdShrink
src/applications/utilities/omdLast
src/applications/utilities/stat2dielectric
src/applications/utilities/stat2tension
src/applications/utilities/stat2thcond
src/applications/utilities/stat2visco
src/applications/hydrodynamics/hydroExplainer
src/applications/utilities/waterRotator
src/applications/utilities/stat2dipolecorr
)

set(PYNUMPY_FILES
src/applications/utilities/hbtetAnalyzer
src/applications/utilities/pack2omd
src/applications/utilities/principalAxisCalculator
src/applications/utilities/protonSampler
src/applications/utilities/slabBuilder
src/applications/utilities/vcorr2spectrum
src/applications/utilities/waterReplacer
src/applications/utilities/wcorr2spectrum
)

set(PYSCIPY_FILES
src/applications/utilities/lcorrzFit
src/applications/utilities/slipLength
src/applications/utilities/solLiqFricCalc
)

IF(SCIPY_FOUND)
  set(PY_FILES ${PYNONUMPY_FILES} ${PYNUMPY_FILES} ${PYSCIPY_FILES})
ELSE(SCIPY_FOUND)
  IF(Python3_NumPy_FOUND)
    set(PY_FILES ${PYNONUMPY_FILES} ${PYNUMPY_FILES})
  ELSE(Python3_NumPy_FOUND)
    set(PY_FILES ${PYNONUMPY_FILES})
  ENDIF(Python3_NumPy_FOUND)
ENDIF(SCIPY_FOUND)

IF(Python3_Interpreter_FOUND)
  IF(WIN32)
    SET(PYEXT ".py")
  ELSE(WIN32)
    SET(PYEXT "")
  ENDIF(WIN32)

  foreach(PY_FILE ${PY_FILES})
    GET_FILENAME_COMPONENT(filename "${PY_FILE}" NAME)
    CONFIGURE_FILE(${PY_FILE}
                  "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${filename}${PYEXT}" @ONLY)
    INSTALL(FILES "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${filename}${PYEXT}"
            DESTINATION bin
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endforeach(PY_FILE)
ENDIF(Python3_Interpreter_FOUND)

set(PERL_FILES
src/applications/utilities/solvator
src/applications/utilities/waterBoxer
doc/git2cl
)

IF(PERL_FOUND)
  IF(WIN32)
    SET(PLEXT ".PL")
  ELSE(WIN32)
    SET(PLEXT "")
  ENDIF(WIN32)

  foreach(PERL_FILE ${PERL_FILES})
     GET_FILENAME_COMPONENT(filename "${PERL_FILE}" NAME)
     configure_file(${PERL_FILE} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${filename}${PLEXT}" @ONLY)
     INSTALL(FILES "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${filename}${PLEXT}"
             DESTINATION bin
             PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

  endforeach(PERL_FILE)
ENDIF(PERL_FOUND)


INSTALL(FILES AUTHORS LICENSE README.md INSTALL DESTINATION .  PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
INSTALL(TARGETS
        openmd_core
        openmd_single
        openmd
        Dump2XYZ
        StaticProps
        DynamicProps
        SequentialProps
        simpleBuilder
        randomBuilder
        nanoparticleBuilder
        icosahedralBuilder
        nanorodBuilder
        nanorod_pentBuilder
        omd2omd
        equationofstate
        thermalizer
        recenter
        elasticConstants
        Hydro
  RUNTIME DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  LIBRARY DESTINATION lib PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
  ARCHIVE DESTINATION lib PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
)

install(DIRECTORY forceFields/
        DESTINATION forceFields
        DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
        PATTERN ".svn" EXCLUDE
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
install(DIRECTORY samples/
        DESTINATION samples
        DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
        PATTERN ".svn" EXCLUDE
	PATTERN "*.dump" EXCLUDE
	PATTERN "*.report" EXCLUDE
	PATTERN "*.stat" EXCLUDE
	PATTERN "*.eor" EXCLUDE
	PATTERN "*.rnemd" EXCLUDE
	PATTERN "*.fz" EXCLUDE
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

configure_file( samples/builders/runMe.in
                "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/runMe" @ONLY)
INSTALL(FILES "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/runMe"
             DESTINATION samples/builders
             PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)


INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/doc/OpenMDmanual.pdf"
        DESTINATION doc
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

message( STATUS)
message( STATUS "========== OpenMD Build Information ==========")
message( STATUS "Current revision ........... = ${GIT_SHA1}")
message( STATUS "Current date/time .......... = ${BUILD_DATE}")
message( STATUS "CMAKE_SYSTEM ............... = ${CMAKE_SYSTEM}")
message( STATUS "==============================================")
message( STATUS "CMAKE_BUILD_TYPE ........... = ${CMAKE_BUILD_TYPE}")
message( STATUS "CMAKE_INSTALL_PREFIX ....... = ${CMAKE_INSTALL_PREFIX}")
message( STATUS "Build as SINGLE_PRECISION .. = ${SINGLE_PRECISION}")
message( STATUS "CMAKE_CXX_COMPILER ......... = ${CMAKE_CXX_COMPILER}")
message( STATUS "MPI_CXX_COMPILER ........... = ${MPI_CXX_COMPILER}")
message( STATUS "MPI_CXX_INCLUDE_PATH ....... = ${MPI_CXX_INCLUDE_PATH}")
message( STATUS "MPI_CXX_LIBRARIES .......... = ${MPI_CXX_LIBRARIES}")
message( STATUS "OPENBABEL_DIR .............. = ${OPENBABEL_DIR}")
message( STATUS "OPENBABEL_INCLUDE_DIR ...... = ${OPENBABEL_INCLUDE_DIR}")
message( STATUS "OPENBABEL_LIBRARIES ........ = ${OPENBABEL_LIBRARIES}")
message( STATUS "QHULL_ROOT ................. = ${QHULL_ROOT}")
message( STATUS "QHULL_INCLUDE_DIR .......... = ${QHULL_INCLUDE_DIR}")
message( STATUS "QHULL_LIBRARIES ............ = ${QHULL_LIBRARIES}")
message( STATUS "ZLIB_ROOT .................. = ${ZLIB_ROOT}")
message( STATUS "ZLIB_INCLUDE_DIR ........... = ${ZLIB_INCLUDE_DIR}")
message( STATUS "ZLIB_LIBRARIES ............. = ${ZLIB_LIBRARIES}")
message( STATUS "FFTW3_ROOT ................. = ${FFTW3_ROOT}")
message( STATUS "FFTW3_INCLUDE_DIR .......... = ${FFTW3_INCLUDE_DIR}")
message( STATUS "FFTW3_LIBRARIES ............ = ${FFTW3_LIBRARIES}")
message( STATUS "PERL_EXECUTABLE ............ = ${PERL_EXECUTABLE}")
message( STATUS "Python3_EXECUTABLE ......... = ${Python3_EXECUTABLE}")
message( STATUS "DOXYGEN_EXECUTABLE ......... = ${DOXYGEN_EXECUTABLE}")
message( STATUS )
message( STATUS "To override these options, add -D{OPTION_NAME}=... to the cmake command" )
message( STATUS "Particularly useful defines are for:")
message( STATUS )
message( STATUS "  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} (where OpenMD will be installed)")
message( STATUS "         -DOPENBABEL_DIR=/path/to/openbabel")
message( STATUS "            -DQHULL_ROOT=/path/to/qhull")
message( STATUS "            -DFFTW3_ROOT=/path/to/fftw3")
message( STATUS )
IF(WIN32)
  message( STATUS "To build and install OpenMD, enter: ")
  message( STATUS "  msbuild /m ALL_BUILD.vcxproj /p:Configuration=${CMAKE_BUILD_TYPE}")
  message( STATUS "followed by:")
  message( STATUS "  msbuild /m INSTALL.vcxproj /p:Configuration=${CMAKE_BUILD_TYPE}")
ELSE(WIN32)
  message( STATUS "To build and install OpenMD, enter \"make\" and \"make install\"")
ENDIF(WIN32)
message( STATUS )

# build a CPack driven installer package
IF(NOT DEFINED CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS)
   SET(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS ON)
ENDIF()
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE
     "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set (CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
include (CPack)
