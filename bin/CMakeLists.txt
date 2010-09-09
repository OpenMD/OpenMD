add_executable(OpenMD ../src/applications/openmd/openmd.cpp)
target_link_libraries(OpenMD openmd)

set (DUMP2XYZSOURCE ../src/applications/dump2Xyz/Dump2XYZ.cpp ../src/applications/dump2Xyz/Dump2XYZCmd.c)

set (DYNAMICPROPSSOURCE ../src/applications/dynamicProps/TimeCorrFunc.cpp  ../src/applications/dynamicProps/CrossTimeCorrFunc.cpp  ../src/applications/dynamicProps/ParticleTimeCorrFunc.cpp  ../src/applications/dynamicProps/FrameTimeCorrFunc.cpp  ../src/applications/dynamicProps/DipoleCorrFunc.cpp  ../src/applications/dynamicProps/RCorrFunc.cpp  ../src/applications/dynamicProps/VCorrFunc.cpp  ../src/applications/dynamicProps/DynamicProps.cpp  ../src/applications/dynamicProps/DynamicPropsCmd.c  ../src/applications/dynamicProps/LegendreCorrFunc.cpp  ../src/applications/dynamicProps/RadialRCorrFunc.cpp 
 ../src/applications/dynamicProps/ThetaCorrFunc.cpp  ../src/applications/dynamicProps/EnergyCorrFunc.cpp  ../src/applications/dynamicProps/DirectionalRCorrFunc.cpp)


set (HYDROSOURCE  ../src/applications/hydrodynamics/AnalyticalModel.cpp ../src/applications/hydrodynamics/ApproximationModel.cpp ../src/applications/hydrodynamics/BeadModel.cpp ../src/applications/hydrodynamics/HydrodynamicsModel.cpp ../src/applications/hydrodynamics/RoughShell.cpp ../src/applications/hydrodynamics/CompositeShape.cpp ../src/applications/hydrodynamics/Hydro.cpp ../src/applications/hydrodynamics/HydroCmd.c ../src/applications/hydrodynamics/ShapeBuilder.cpp ../src/applications/hydrodynamics/HydrodynamicsModelFactory.cpp)

set (STATICPROPSSOURCE  ../src/applications/staticProps/GofAngle2.cpp ../src/applications/staticProps/GofRAngle.cpp ../src/applications/staticProps/GofRZ.cpp ../src/applications/staticProps/GofR.cpp
../src/applications/staticProps/GofZ.cpp  
../src/applications/staticProps/GofXyz.cpp
../src/applications/staticProps/TwoDGofR.cpp 
../src/applications/staticProps/RadialDistrFunc.cpp 
../src/applications/staticProps/P2OrderParameter.cpp 
../src/applications/staticProps/SCDOrderParameter.cpp
../src/applications/staticProps/DensityPlot.cpp
../src/applications/staticProps/StaticPropsCmd.c
../src/applications/staticProps/StaticProps.cpp
../src/applications/staticProps/RhoZ.cpp
../src/applications/staticProps/Hxy.cpp
../src/applications/staticProps/RippleOP.cpp
../src/applications/staticProps/BondOrderParameter.cpp
../src/applications/staticProps/BOPofR.cpp
../src/applications/staticProps/RhoR.cpp
../src/applications/staticProps/BondAngleDistribution.cpp
../src/applications/staticProps/pAngle.cpp
../src/applications/staticProps/NanoVolume.cpp)

set (NANOPARTICLEBUILDERSOURCE ../src/applications/nanoparticleBuilder/nanoparticleBuilderCmd.c ../src/applications/nanoparticleBuilder/nanoparticleBuilder.cpp  ../src/applications/nanoparticleBuilder/shapedLatticeSpherical.cpp)

set (RANDOMBUILDERSOURCE ../src/applications/randomBuilder/randomBuilderCmd.c ../src/applications/randomBuilder/randomBuilder.cpp)

set(SIMPLEBUILDERSOURCE ../src/applications/simpleBuilder/simpleBuilderCmd.c ../src/applications/simpleBuilder/simpleBuilder.cpp )

set(THERMALIZERSOURCE ../src/applications/thermalizer/thermalizerCmd.c ../src/applications/thermalizer/thermalizer.cpp)

add_executable(Dump2XYZ ${DUMP2XYZSOURCE})
target_link_libraries(Dump2XYZ openmd)
add_executable(DynamicProps ${DYNAMICPROPSSOURCE})
target_link_libraries(DynamicProps openmd)
add_executable(Hydro ${HYDROSOURCE})
target_link_libraries(Hydro openmd)
add_executable(StaticProps ${STATICPROPSSOURCE})
target_link_libraries(StaticProps openmd)
add_executable(nanoparticleBuilder ${NANOPARTICLEBUILDERSOURCE})
target_link_libraries(nanoparticleBuilder openmd)
add_executable(randomBuilder ${RANDOMBUILDERSOURCE})
target_link_libraries(randomBuilder openmd)
add_executable(simpleBuilder ${SIMPLEBUILDERSOURCE})
target_link_libraries(simpleBuilder openmd)
add_executable(thermalizer ${THERMALIZERSOURCE})
target_link_libraries(thermalizer openmd)

if (OPENBABEL2_FOUND)
  set (ATOM2MDSOURCE ../src/applications/atom2md/atom2md.cpp ../src/applications/atom2md/openmdformat.cpp)
  add_executable(atom2md ${ATOM2MDSOURCE})
  target_link_libraries(atom2md openmd ${OPENBABEL2_LIBRARIES})
ENDIF (OPENBABEL2_FOUND)