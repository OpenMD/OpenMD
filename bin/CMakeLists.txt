add_executable(openmd ../src/applications/openmd/openmd.cpp)
target_link_libraries(openmd openmd_core openmd_single)

if (MPI_FOUND)
add_executable(openmd_MPI ../src/applications/openmd/openmd.cpp)
set_target_properties(openmd_MPI PROPERTIES
COMPILE_DEFINITIONS IS_MPI
)
target_link_libraries(openmd_MPI openmd_core openmd_parallel)
ENDIF (MPI_FOUND)

set (DUMP2XYZSOURCE 
../src/applications/dump2Xyz/Dump2XYZ.cpp
../src/applications/dump2Xyz/Dump2XYZCmd.c
)

set (DYNAMICPROPSSOURCE 
../src/applications/dynamicProps/ActionCorrFunc.cpp
../src/applications/dynamicProps/CrossTimeCorrFunc.cpp
../src/applications/dynamicProps/DipoleCorrFunc.cpp
../src/applications/dynamicProps/DirectionalRCorrFunc.cpp
../src/applications/dynamicProps/DynamicProps.cpp
../src/applications/dynamicProps/EnergyCorrFunc.cpp
../src/applications/dynamicProps/FrameTimeCorrFunc.cpp
../src/applications/dynamicProps/LegendreCorrFunc.cpp
../src/applications/dynamicProps/MomentumCorrFunc.cpp
../src/applications/dynamicProps/ParticleTimeCorrFunc.cpp
../src/applications/dynamicProps/RadialRCorrFunc.cpp
../src/applications/dynamicProps/RCorrFunc.cpp
../src/applications/dynamicProps/StressCorrFunc.cpp
../src/applications/dynamicProps/SystemDipoleCorrFunc.cpp
../src/applications/dynamicProps/ThetaCorrFunc.cpp
../src/applications/dynamicProps/TimeCorrFunc.cpp
../src/applications/dynamicProps/VCorrFunc.cpp
../src/applications/dynamicProps/DynamicPropsCmd.c
)

set (HYDROSOURCE 
../src/applications/hydrodynamics/AnalyticalModel.cpp
../src/applications/hydrodynamics/ApproximationModel.cpp
../src/applications/hydrodynamics/BeadModel.cpp
../src/applications/hydrodynamics/CompositeShape.cpp
../src/applications/hydrodynamics/Hydro.cpp
../src/applications/hydrodynamics/HydrodynamicsModel.cpp
../src/applications/hydrodynamics/HydrodynamicsModelFactory.cpp
../src/applications/hydrodynamics/RoughShell.cpp
../src/applications/hydrodynamics/ShapeBuilder.cpp
../src/applications/hydrodynamics/HydroCmd.c
)

set (STATICPROPSSOURCE 
../src/applications/staticProps/AngleR.cpp
../src/applications/staticProps/BondAngleDistribution.cpp
../src/applications/staticProps/BondOrderParameter.cpp
../src/applications/staticProps/BOPofR.cpp
../src/applications/staticProps/DensityPlot.cpp
../src/applications/staticProps/GofAngle2.cpp
../src/applications/staticProps/GofR.cpp
../src/applications/staticProps/GofRAngle.cpp
../src/applications/staticProps/GofRZ.cpp
../src/applications/staticProps/GofXyz.cpp
../src/applications/staticProps/GofZ.cpp
../src/applications/staticProps/Hxy.cpp
../src/applications/staticProps/NanoLength.cpp
../src/applications/staticProps/NanoVolume.cpp
../src/applications/staticProps/ObjectCount.cpp
../src/applications/staticProps/P2OrderParameter.cpp
../src/applications/staticProps/pAngle.cpp
../src/applications/staticProps/RadialDistrFunc.cpp
../src/applications/staticProps/RhoR.cpp
../src/applications/staticProps/RhoZ.cpp
../src/applications/staticProps/RippleOP.cpp
../src/applications/staticProps/SCDOrderParameter.cpp
../src/applications/staticProps/StaticProps.cpp
../src/applications/staticProps/TetrahedralityParam.cpp
../src/applications/staticProps/TwoDGofR.cpp
../src/applications/staticProps/StaticPropsCmd.c
)

set (NANOPARTICLEBUILDERSOURCE 
../src/applications/nanoparticleBuilder/nanoparticleBuilder.cpp
../src/applications/nanoparticleBuilder/shapedLatticeSpherical.cpp
../src/applications/nanoparticleBuilder/nanoparticleBuilderCmd.c
)

set (RANDOMBUILDERSOURCE
../src/applications/randomBuilder/randomBuilder.cpp
../src/applications/randomBuilder/randomBuilderCmd.c
)

set(SIMPLEBUILDERSOURCE 
../src/applications/simpleBuilder/simpleBuilder.cpp
../src/applications/simpleBuilder/simpleBuilderCmd.c
)

set(THERMALIZERSOURCE 
../src/applications/thermalizer/thermalizer.cpp
../src/applications/thermalizer/thermalizerCmd.c 
)

add_executable(Dump2XYZ ${DUMP2XYZSOURCE})
target_link_libraries(Dump2XYZ openmd_core openmd_single)
add_executable(DynamicProps ${DYNAMICPROPSSOURCE})
target_link_libraries(DynamicProps openmd_core openmd_single)
add_executable(Hydro ${HYDROSOURCE})
target_link_libraries(Hydro openmd_core openmd_single)
add_executable(StaticProps ${STATICPROPSSOURCE})
target_link_libraries(StaticProps openmd_core openmd_single)
add_executable(nanoparticleBuilder ${NANOPARTICLEBUILDERSOURCE})
target_link_libraries(nanoparticleBuilder openmd_core openmd_single)
add_executable(randomBuilder ${RANDOMBUILDERSOURCE})
target_link_libraries(randomBuilder openmd_core openmd_single)
add_executable(simpleBuilder ${SIMPLEBUILDERSOURCE})
target_link_libraries(simpleBuilder openmd_core openmd_single)
add_executable(thermalizer ${THERMALIZERSOURCE})
target_link_libraries(thermalizer openmd_core openmd_single)

if (OPENBABEL2_FOUND)
set (ATOM2MDSOURCE 
../src/applications/atom2md/atom2md.cpp 
../src/applications/atom2md/openmdformat.cpp
)
  add_executable(atom2md ${ATOM2MDSOURCE})
  target_link_libraries(atom2md openmd_core openmd_single ${OPENBABEL2_LIBRARIES})
ENDIF (OPENBABEL2_FOUND)
