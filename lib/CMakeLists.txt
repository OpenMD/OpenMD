set(BRAINSSOURCE ../src/brains/PairList.cpp ../src/brains/SimInfo.cpp ../src/brains/DataStorage.cpp ../src/brains/SimSnapshotManager.cpp ../src/brains/Snapshot.cpp ../src/brains/Stats.cpp ../src/brains/Thermo.cpp ../src/brains/ForceManager.cpp ../src/brains/MoleculeCreator.cpp ../src/brains/SimCreator.cpp ../src/brains/Register.cpp ../src/brains/BlockSnapshotManager.cpp)

set(NONBONDEDSOURCE ../src/nonbonded/EAM.cpp ../src/nonbonded/GB.cpp ../src/nonbonded/LJ.cpp 
../src/nonbonded/SC.cpp ../src/nonbonded/Sticky.cpp )

set(USETHEFORCESOURCE ../src/UseTheForce/doForces.F90 ../src/UseTheForce/doForces_interface.F90 ../src/UseTheForce/ForceField.cpp ../src/UseTheForce/ForceFieldFactory.cpp ../src/UseTheForce/DUFF.cpp ../src/UseTheForce/EAM_FF.cpp ../src/UseTheForce/SHAPES_FF.cpp ../src/UseTheForce/SC_FF.cpp ../src/UseTheForce/MnM_FF.cpp ../src/UseTheForce/Amber_FF.cpp ../src/UseTheForce/CLAYFF.cpp)

set(DARKSIDESOURCE ../src/UseTheForce/DarkSide/simParallel.F90 ../src/UseTheForce/DarkSide/simParallel_interface.F90 ../src/UseTheForce/DarkSide/neighborLists.F90  ../src/UseTheForce/DarkSide/neighborLists_interface.F90 ../src/UseTheForce/DarkSide/atype.F90 ../src/UseTheForce/DarkSide/force_globals.F90 ../src/UseTheForce/DarkSide/fForceOptions.F90 ../src/UseTheForce/DarkSide/fForceOptions_interface.F90 ../src/UseTheForce/DarkSide/switcheroo.F90 ../src/UseTheForce/DarkSide/switcheroo_interface.F90 ../src/UseTheForce/DarkSide/simulation.F90 ../src/UseTheForce/DarkSide/simulation_interface.F90 ../src/UseTheForce/DarkSide/electrostatic.F90 ../src/UseTheForce/DarkSide/electrostatic_interface.F90 ../src/UseTheForce/DarkSide/shapes.F90 ../src/UseTheForce/DarkSide/shapes_interface.F90 ../src/UseTheForce/DarkSide/MetalNonMetal.F90 ../src/UseTheForce/DarkSide/MetalNonMetal_interface.F90
../src/UseTheForce/DarkSide/atype_interface.F90)

set(UTILSSOURCE ../src/utils/definitions.F90 ../src/utils/status.F90 ../src/utils/interpolation.F90 ../src/utils/fError.c ../src/utils/PropertyMap.cpp ../src/utils/OpenMDMPI.F90 ../src/utils/simError.c ../src/utils/FileLocation.cpp ../src/utils/SystemUtilities.cpp ../src/utils/Exception.cpp ../src/utils/Report.cpp ../src/utils/StringUtils.cpp ../src/utils/vector_class.F90 ../src/utils/StringTokenizer.cpp ../src/utils/Utility.cpp ../src/utils/Trim.cpp ../src/utils/OpenMDBitSet.cpp ../src/utils/MoLocator.cpp ../src/utils/wildcards.cpp ../src/utils/ElementsTable.cpp ../src/utils/ProgressBar.cpp ../src/utils/residentMem.c ../src/utils/physmem.c)

set(MATHSOURCE ../src/math/ParallelRandNumGen.cpp ../src/math/SeqRandNumGen.cpp ../src/math/RealSphericalHarmonic.cpp ../src/math/SphericalHarmonic.cpp ../src/math/ChebyshevT.cpp ../src/math/ChebyshevU.cpp ../src/math/LegendrePolynomial.cpp ../src/math/RMSD.cpp ../src/math/CubicSpline.cpp ../src/math/linearAlgebra.F90 ../src/math/Wigner3jm.F90)

set(ANTLRSOURCE ../src/antlr/ANTLRUtil.cpp    ../src/antlr/ASTFactory.cpp ../src/antlr/ASTNULLType.cpp ../src/antlr/ASTRefCount.cpp ../src/antlr/BaseAST.cpp ../src/antlr/BitSet.cpp ../src/antlr/CharBuffer.cpp ../src/antlr/CharScanner.cpp ../src/antlr/CommonAST.cpp ../src/antlr/CommonASTWithHiddenTokens.cpp ../src/antlr/CommonHiddenStreamToken.cpp ../src/antlr/CommonToken.cpp ../src/antlr/InputBuffer.cpp ../src/antlr/LLkParser.cpp ../src/antlr/MismatchedCharException.cpp ../src/antlr/MismatchedTokenException.cpp ../src/antlr/NoViableAltException.cpp ../src/antlr/NoViableAltForCharException.cpp ../src/antlr/Parser.cpp ../src/antlr/RecognitionException.cpp ../src/antlr/String.cpp ../src/antlr/TokenBuffer.cpp ../src/antlr/Token.cpp ../src/antlr/TokenRefCount.cpp ../src/antlr/TokenStreamBasicFilter.cpp ../src/antlr/TokenStreamHiddenTokenFilter.cpp ../src/antlr/TokenStreamRewriteEngine.cpp ../src/antlr/TokenStreamSelector.cpp ../src/antlr/TreeParser.cpp)

set(DEBUGSOURCE ../src/debug/Backtrace.cpp ../src/debug/Debugger.cpp ../src/debug/Pipe.cpp ../src/debug/Process.cpp)

set(CONSTRAINTSOURCE ../src/constraints/Rattle.cpp ../src/constraints/ZconstraintForceManager.cpp)

set(HYDRODYNAMICSSOURCE ../src/hydrodynamics/HydroProp.cpp ../src/hydrodynamics/Sphere.cpp ../src/hydrodynamics/Ellipsoid.cpp)

set(INTEGRATORSOURCE ../src/integrators/Integrator.cpp ../src/integrators/IntegratorFactory.cpp ../src/integrators/Velocitizer.cpp ../src/integrators/RNEMD.cpp ../src/integrators/VelocityVerletIntegrator.cpp ../src/integrators/DLM.cpp ../src/integrators/NVE.cpp ../src/integrators/NVT.cpp ../src/integrators/NPT.cpp ../src/integrators/NPTi.cpp ../src/integrators/NPTf.cpp ../src/integrators/NPTxyz.cpp ../src/integrators/NPAT.cpp ../src/integrators/NPrT.cpp ../src/integrators/NgammaT.cpp ../src/integrators/LDForceManager.cpp ../src/integrators/LangevinDynamics.cpp ../src/integrators/)

set(IOSOURCE ../src/io/AtomTypesSectionParser.cpp ../src/io/BaseAtomTypesSectionParser.cpp ../src/io/DirectionalAtomTypesSectionParser.cpp ../src/io/SectionParser.cpp ../src/io/EAMAtomTypesSectionParser.cpp ../src/io/SectionParserManager.cpp ../src/io/ChargeAtomTypesSectionParser.cpp ../src/io/MultipoleAtomTypesSectionParser.cpp ../src/io/StickyAtomTypesSectionParser.cpp ../src/io/GayBerneAtomTypesSectionParser.cpp ../src/io/TorsionTypesSectionParser.cpp ../src/io/InversionTypesSectionParser.cpp ../src/io/StickyPowerAtomTypesSectionParser.cpp ../src/io/ShapeAtomTypesSectionParser.cpp ../src/io/LennardJonesAtomTypesSectionParser.cpp ../src/io/SCAtomTypesSectionParser.cpp ../src/io/MetalNonMetalInteractionsSectionParser.cpp ../src/io/BondTypesSectionParser.cpp ../src/io/BendTypesSectionParser.cpp ../src/io/Globals.cpp ../src/io/DumpReader.cpp ../src/io/DumpWriter.cpp ../src/io/StatWriter.cpp ../src/io/ZConsReader.cpp ../src/io/ZConsWriter.cpp ../src/io/RestReader.cpp ../src/io/RestWriter.cpp ../src/io/gzstream.cpp ../src/io/ParamConstraint.cpp ../src/io/OptionSectionParser.cpp ../src/io/ForceFieldOptions.cpp)

set(LATTICESOURCE ../src/lattice/Lattice.cpp ../src/lattice/CubicLattice.cpp ../src/lattice/FCCLattice.cpp ../src/lattice/LatticeFactory.cpp ../src/lattice/shapedLattice.cpp)

set(MDPARSERSOURCE ../src/mdParser/FilenameObserver.cpp ../src/mdParser/MDParser.cpp ../src/mdParser/MDLexer.cpp ../src/mdParser/MDTreeParser.cpp)

set(MINIMIZERSOURCE ../src/minimizers/MinimizerParameterSet.cpp ../src/minimizers/Minimizer.cpp ../src/minimizers/MinimizerFactory.cpp ../src/minimizers/CGFamilyMinimizer.cpp ../src/minimizers/PRCG.cpp ../src/minimizers/SDMinimizer.cpp)

set(PRIMITIVESOURCE   ../src/primitives/StuntDouble.cpp ../src/primitives/Atom.cpp ../src/primitives/DirectionalAtom.cpp ../src/primitives/RigidBody.cpp ../src/primitives/Molecule.cpp ../src/primitives/Bend.cpp ../src/primitives/GhostBend.cpp ../src/primitives/Torsion.cpp ../src/primitives/Inversion.cpp ../src/primitives/GhostTorsion.cpp ../src/primitives/UreyBradleyBend.cpp )

#THIS PROFILEING CODE SHOULD BE REMOVED
set(PROFILINGSOURCE ../src/profiling/mdProfile.cpp ../src/profiling/timing.F90)

set(RESTRAINTSOURCE  ../src/restraints/ObjectRestraint.cpp ../src/restraints/MolecularRestraint.cpp ../src/restraints/RestraintForceManager.cpp ../src/restraints/ThermoIntegrationForceManager.cpp)

set(SELECTIONSOURCE   ../src/selection/SelectionToken.cpp ../src/selection/TokenMap.cpp ../src/selection/SelectionCompiler.cpp ../src/selection/SelectionEvaluator.cpp ../src/selection/NameFinder.cpp ../src/selection/DistanceFinder.cpp  ../src/selection/HullFinder.cpp ../src/selection/IndexFinder.cpp ../src/selection/SelectionManager.cpp)

set(TYPESSOURCE  ../src/types/AtomStamp.cpp ../src/types/AtomType.cpp ../src/types/BendStamp.cpp  ../src/types/BondStamp.cpp ../src/types/Component.cpp ../src/types/CutoffGroupStamp.cpp ../src/types/MoleculeStamp.cpp ../src/types/RigidBodyStamp.cpp ../src/types/TorsionStamp.cpp ../src/types/InversionStamp.cpp ../src/types/ZconsStamp.cpp ../src/types/RestraintStamp.cpp ../src/types/AtomType.cpp ../src/types/DirectionalAtomType.cpp ../src/types/ShapeAtomType.cpp ../src/types/CharmmTorsionType.cpp ../src/types/ImproperCosineInversionType.cpp ../src/types/NonBondedInteractionType.cpp ../src/types/FragmentStamp.cpp)

set(VISITORSOURCE ../src/visitors/AtomVisitor.cpp ../src/visitors/ReplacementVisitor.cpp ../src/visitors/CompositeVisitor.cpp ../src/visitors/OtherVisitor.cpp ../src/visitors/RigidBodyVisitor.cpp ../src/visitors/ZconsVisitor.cpp ../src/visitors/LipidTransVisitor.cpp ../src/visitors/AtomNameVisitor.cpp)

IF(QHULL_FOUND)
  set(QHULLSRC ../src/integrators/SMIPDForceManager.cpp ../src/integrators/SMIPDynamics.cpp ../src/math/ConvexHull.cpp ../src/math/AlphaHull.cpp ../src/math/Triangle.cpp)
ENDIF(QHULL_FOUND)

add_library(openmd SHARED ${UTILSSOURCE} ${MATHSOURCE} ${ANTLRSOURCE} ${DEBUGSOURCE} ${CONSTRAINTSOURCE} ${HYDRODYNAMICSSOURCE} ${INTEGRATORSOURCE} ${BRAINSSOURCE} ${IOSOURCE} ${USETHEFORCESOURCE} ${DARKSIDESOURCE} ${NONBONDEDSOURCE} ${LATTICESOURCE} ${LATTICESOURCE} ${MDPARSERSOURCE} ${PRIMITIVESOURCE} ${MINIMIZERSOURCE} ${PROFILINGSOURCE} ${RESTRAINTSOURCE} ${SELECTIONSOURCE} ${TYPESSOURCE} ${VISITORSOURCE} ${QHULLSRC} )
