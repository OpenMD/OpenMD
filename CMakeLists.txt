# Top level cmake script for OpenMD. 
project(OpenMD CXX)
cmake_minimum_required(VERSION 2.8.5)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

find_package(MPI)
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
ENDIF()

include(CheckCXXCompilerFlag)
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)

# OpenMD version number.
set (VERSION_MAJOR "2")
set (VERSION_MINOR "0")
set (VERSION_TINY "0")  
option(SINGLE_PRECISION "Build Single precision (float) version" OFF)

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

FIND_PACKAGE(Subversion)
IF ( Subversion_FOUND AND EXISTS ${PROJECT_SOURCE_DIR}/.svn )
  Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
  MESSAGE(STATUS "Current revision is ${Project_WC_REVISION}")
  Subversion_WC_LOG(${PROJECT_SOURCE_DIR} Project)
  SET(SVN_REV ${Project_WC_REVISION})
  ADD_DEFINITIONS( -DSVN_REV=${Project_WC_REVISION} )
  SET(GENERATELOGS "${PROJECT_SOURCE_DIR}/doc/svn2cl")
  ADD_CUSTOM_TARGET(ChangeLog
     COMMAND ${GENERATELOGS} "--group-by-day" "--include-rev"
     "--linelen=78" "--output=${PROJECT_BINARY_DIR}/ChangeLog" "${PROJECT_SOURCE_DIR}"
  )
ENDIF()

check_include_file_cxx(conio.h      HAVE_CONIO_H)
check_cxx_symbol_exists(strncasecmp   "string.h"   HAVE_STRNCASECMP)

# Optional libraries: If we can find these, we will build with them
# Look for OpenBabel libraries
find_package(OpenBabel2)
IF(OPENBABEL2_FOUND)
  SET(USE_OPENBABEL)
  include_directories(${OPENBABEL2_INCLUDE_DIR})
  # link libraries are added only for targets that need them.
ELSE(OPENBABEL2_FOUND)
  MESSAGE(STATUS "No OpenBabel found - will not build atom2md")
ENDIF(OPENBABEL2_FOUND)

#Look for QHULL Libraries
set(QHULL_USE_STATIC 1)
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

FIND_PACKAGE(PythonInterp)
IF(PYTHON_EXECUTABLE)
  SET(PYTHON ${PYTHON_EXECUTABLE})
ELSE(PYTHON_EXECUTABLE)
  MESSAGE(STATUS "Failed to find python - some scripts will not be made")
ENDIF(PYTHON_EXECUTABLE)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.cmake"
  "${PROJECT_BINARY_DIR}/config.h"
)
include_directories("${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/src")

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
src/integrators/FluctuatingChargeNVT.cpp
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
src/io/FluctuatingChargeAtomTypesSectionParser.cpp
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
src/io/PolarizableAtomTypesSectionParser.cpp
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
src/types/DirectionalAdapter.cpp
src/types/EAMAdapter.cpp
src/types/FixedChargeAdapter.cpp
src/types/FluctuatingChargeAdapter.cpp
src/types/FragmentStamp.cpp
src/types/GayBerneAdapter.cpp
src/types/ImproperCosineInversionType.cpp
src/types/InversionStamp.cpp
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
src/utils/simError.cpp
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
src/applications/dump2Xyz/Dump2XYZCmd.cpp
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
src/applications/staticProps/StaticPropsCmd.cpp
)

set (NANOPARTICLEBUILDERSOURCE 
src/applications/nanoparticleBuilder/nanoparticleBuilder.cpp
src/applications/nanoparticleBuilder/shapedLatticeSpherical.cpp
src/applications/nanoparticleBuilder/nanoparticleBuilderCmd.cpp
)

set (NANORODBUILDERSOURCE 
src/applications/nanoparticleBuilder/nanorodBuilder.cpp
src/applications/nanoparticleBuilder/shapedLatticeRod.cpp
src/applications/nanoparticleBuilder/nanorodBuilderCmd.cpp
)

set (NANOROD_PENTBUILDERSOURCE 
src/applications/nanoparticleBuilder/nanorod_pentBuilder.cpp
src/applications/nanoparticleBuilder/nanorod_pentBuilderCmd.cpp
src/applications/nanoparticleBuilder/shapedLatticePentRod.cpp
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
add_executable(nanorodBuilder ${NANORODBUILDERSOURCE})
target_link_libraries(nanorodBuilder openmd_single openmd_core openmd_single openmd_core)
add_executable(nanorod_pentBuilder ${NANOROD_PENTBUILDERSOURCE})
target_link_libraries(nanorod_pentBuilder openmd_single openmd_core openmd_single openmd_core)
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
src/applications/utilities/waterReplacer
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

INSTALL(FILES AUTHORS LICENSE README INSTALL DESTINATION . 
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
        DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
        PATTERN ".svn" EXCLUDE
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
install(DIRECTORY samples/       
        DESTINATION samples 
        DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
        PATTERN ".svn" EXCLUDE
	PATTERN "*.dump" EXCLUDE
	PATTERN "*.stat" EXCLUDE
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/doc/OpenMDmanual.pdf"
        DESTINATION doc
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

message( STATUS)
message( STATUS "========== OpenMD Build Information ==========")
message( STATUS "Current revision ........... = ${SVN_REV}")
message( STATUS "CMAKE_SYSTEM ............... = ${CMAKE_SYSTEM}")
message( STATUS "==============================================")
message( STATUS "CMAKE_BUILD_TYPE ........... = ${CMAKE_BUILD_TYPE}")
message( STATUS "CMAKE_INSTALL_PREFIX ....... = ${CMAKE_INSTALL_PREFIX}")
message( STATUS "Build as SINGLE_PRECISION .. = ${SINGLE_PRECISION}")
message( STATUS "CMAKE_CXX_COMPILER ......... = ${CMAKE_CXX_COMPILER}")
message( STATUS "MPI_CXX_COMPILER ........... = ${MPI_CXX_COMPILER}")
message( STATUS "MPI_CXX_INCLUDE_PATH ....... = ${MPI_CXX_INCLUDE_PATH}")
message( STATUS "MPI_CXX_LIBRARIES .......... = ${MPI_CXX_LIBRARIES}")
message( STATUS "OPENBABEL2_ROOT ............ = ${OPENBABEL2_ROOT}")
message( STATUS "OPENBABEL2_INCLUDE_DIR ..... = ${OPENBABEL2_INCLUDE_DIR}")
message( STATUS "OPENBABEL2_LIBRARIES ....... = ${OPENBABEL2_LIBRARIES}")
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
message( STATUS "PYTHON_EXECUTABLE .......... = ${PYTHON_EXECUTABLE}")
message( STATUS "DOXYGEN_EXECUTABLE ......... = ${DOXYGEN_EXECUTABLE}")
message( STATUS )
message( STATUS "To override these options, add -D{OPTION_NAME}=... to the cmake command" )
message( STATUS "Particularly useful defines are for:")
message( STATUS )
message( STATUS "  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} (where OpenMD will be installed)")
message( STATUS "       -DOPENBABEL2_ROOT=/path/to/openbabel")
message( STATUS "            -DQHULL_ROOT=/path/to/qhull")
message( STATUS "            -DFFTW3_ROOT=/path/to/fftw3")
message( STATUS )
message( STATUS "To build and install OpenMD, enter \"make\" and \"make install\"")
message( STATUS )
