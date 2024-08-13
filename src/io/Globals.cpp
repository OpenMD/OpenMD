/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "io/Globals.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "io/ParamConstraint.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {
  Globals::Globals() {
    flucQpars_     = new FluctuatingChargeParameters();
    rnemdPars_     = new RNEMD::RNEMDParameters();
    lightPars_     = new Perturbations::LightParameters();
    minimizerPars_ = new MinimizerParameters();

    DefineParameter(ForceField, "forceField");

    DefineOptionalParameter(TargetTemp, "targetTemp");
    DefineOptionalParameter(Ensemble, "ensemble");
    DefineOptionalParameter(Dt, "dt");
    DefineOptionalParameter(RunTime, "runTime");
    DefineOptionalParameter(FinalConfig, "finalConfig");
    DefineOptionalParameter(SampleTime, "sampleTime");
    DefineOptionalParameter(ResetTime, "resetTime");
    DefineOptionalParameter(StatusTime, "statusTime");
    DefineOptionalParameter(CutoffRadius, "cutoffRadius");
    DefineOptionalParameter(SwitchingRadius, "switchingRadius");
    DefineOptionalParameter(TempSet, "tempSet");
    DefineOptionalParameter(ThermalTime, "thermalTime");
    DefineOptionalParameter(TargetPressure, "targetPressure");
    DefineOptionalParameter(TauThermostat, "tauThermostat");
    DefineOptionalParameter(TauBarostat, "tauBarostat");
    DefineOptionalParameter(LangevinPistonDrag, "langevinPistonDrag");
    DefineOptionalParameter(ZconsTime, "zconsTime");
    DefineOptionalParameter(ZconsTol, "zconsTol");
    DefineOptionalParameter(ZconsForcePolicy, "zconsForcePolicy");
    DefineOptionalParameter(Seed, "seed");
    DefineOptionalParameter(ZconsGap, "zconsGap");
    DefineOptionalParameter(ZconsFixtime, "zconsFixtime");
    DefineOptionalParameter(ZconsUsingSMD, "zconsUsingSMD");
    DefineOptionalParameter(ThermodynamicIntegrationLambda,
                            "thermodynamicIntegrationLambda");
    DefineOptionalParameter(ThermodynamicIntegrationK,
                            "thermodynamicIntegrationK");
    DefineOptionalParameter(ForceFieldVariant, "forceFieldVariant");
    DefineOptionalParameter(ForceFieldFileName, "forceFieldFileName");
    DefineOptionalParameter(DampingAlpha, "dampingAlpha");
    DefineOptionalParameter(SurfaceTension, "surfaceTension");
    DefineOptionalParameter(PrintPressureTensor, "printPressureTensor");
    DefineOptionalParameter(PrintVirialTensor, "printVirialTensor");
    DefineOptionalParameter(ElectricField, "electricField");
    DefineOptionalParameter(UniformField, "uniformField");
    // magnetic field optional parameter added
    DefineOptionalParameter(MagneticField, "magneticField");
    DefineOptionalParameter(UniformGradientStrength, "uniformGradientStrength");
    DefineOptionalParameter(UniformGradientDirection1,
                            "uniformGradientDirection1");
    DefineOptionalParameter(UniformGradientDirection2,
                            "uniformGradientDirection2");
    // DefineOptionalParameter(PeriodicField, "periodicField");
    // DefineOptionalParameter(PeriodicFieldStrength, "periodicFieldStrength");

    DefineOptionalParameter(TaggedAtomPair, "taggedAtomPair");
    DefineOptionalParameter(PrintTaggedPairDistance, "printTaggedPairDistance");
    DefineOptionalParameter(SwitchingFunctionType, "switchingFunctionType");
    DefineOptionalParameter(HydroPropFile, "HydroPropFile");
    DefineOptionalParameter(Viscosity, "viscosity");
    DefineOptionalParameter(BeadSize, "beadSize");
    DefineOptionalParameter(FrozenBufferRadius, "frozenBufferRadius");
    DefineOptionalParameter(LangevinBufferRadius, "langevinBufferRadius");
    DefineOptionalParameter(NeighborListNeighbors, "NeighborListNeighbors");
    DefineOptionalParameter(UseMultipleTemperatureMethod,
                            "useMultipleTemperatureMethod");
    DefineOptionalParameter(ElectrostaticSummationMethod,
                            "electrostaticSummationMethod");
    DefineOptionalParameter(MTM_Ce, "MTM_Ce");
    DefineOptionalParameter(MTM_G, "MTM_G");
    DefineOptionalParameter(MTM_Io, "MTM_Io");
    DefineOptionalParameter(MTM_Sigma, "MTM_Sigma");
    DefineOptionalParameter(MTM_R, "MTM_R");
    DefineOptionalParameter(Alpha, "alpha");
    DefineOptionalParameter(ConstraintTime, "constraintTime");

    DefineOptionalParameter(PotentialSelection, "potentialSelection");

    DefineOptionalParameterWithDefaultValue(SkipPairLoop, "skipPairLoop",
                                            false);
    DefineOptionalParameterWithDefaultValue(
        UsePeriodicBoundaryConditions, "usePeriodicBoundaryConditions", true);
    DefineOptionalParameterWithDefaultValue(ConserveLinearMomentum,
                                            "conserveLinearMomentum", true);
    DefineOptionalParameterWithDefaultValue(ConserveAngularMomentum,
                                            "conserveAngularMomentum", true);
    DefineOptionalParameterWithDefaultValue(UseAtomicVirial, "useAtomicVirial",
                                            true);
    DefineOptionalParameterWithDefaultValue(UseLongRangeCorrections,
                                            "useLongRangeCorrections", true);
    DefineOptionalParameterWithDefaultValue(UseInitalTime, "useInitialTime",
                                            false);
    DefineOptionalParameterWithDefaultValue(
        UseIntialExtendedSystemState, "useInitialExtendedSystemState", false);
    DefineOptionalParameterWithDefaultValue(OrthoBoxTolerance,
                                            "orthoBoxTolerance", 1E-6);
    DefineOptionalParameterWithDefaultValue(CutoffMethod, "cutoffMethod",
                                            "SHIFTED_FORCE");
    DefineOptionalParameterWithDefaultValue(
        ElectrostaticScreeningMethod, "electrostaticScreeningMethod", "DAMPED");
    DefineOptionalParameter(UseSurfaceTerm, "useSurfaceTerm");
    DefineOptionalParameter(UseSlabGeometry, "useSlabGeometry");
    DefineOptionalParameterWithDefaultValue(Dielectric, "dielectric", 80.0);
    DefineOptionalParameterWithDefaultValue(CompressDumpFile,
                                            "compressDumpFile", false);
    DefineOptionalParameterWithDefaultValue(PrintHeatFlux, "printHeatFlux",
                                            false);
    DefineOptionalParameterWithDefaultValue(OutputForceVector,
                                            "outputForceVector", false);
    DefineOptionalParameterWithDefaultValue(OutputParticlePotential,
                                            "outputParticlePotential", false);
    DefineOptionalParameterWithDefaultValue(OutputElectricField,
                                            "outputElectricField", false);
    DefineOptionalParameterWithDefaultValue(OutputFluctuatingCharges,
                                            "outputFluctuatingCharges", false);
    DefineOptionalParameterWithDefaultValue(OutputSitePotential,
                                            "outputSitePotential", false);
    DefineOptionalParameterWithDefaultValue(OutputDensity, "outputDensity",
                                            false);
    DefineOptionalParameterWithDefaultValue(SkinThickness, "skinThickness",
                                            1.0);
    DefineOptionalParameterWithDefaultValue(
        StatFileFormat, "statFileFormat",
        "TIME|TOTAL_ENERGY|POTENTIAL_ENERGY|KINETIC_ENERGY|TEMPERATURE|"
        "PRESSURE|"
        "VOLUME|"
        "CONSERVED_QUANTITY");
    DefineOptionalParameterWithDefaultValue(StatFilePrecision,
                                            "statFilePrecision", 8);
    DefineOptionalParameterWithDefaultValue(UseSphericalBoundaryConditions,
                                            "useSphericalBoundaryConditions",
                                            false);
    DefineOptionalParameterWithDefaultValue(AccumulateBoxDipole,
                                            "accumulateBoxDipole", false);
    DefineOptionalParameterWithDefaultValue(AccumulateBoxQuadrupole,
                                            "accumulateBoxQuadrupole", false);
    DefineOptionalParameterWithDefaultValue(UseRestraints, "useRestraints",
                                            false);
    DefineOptionalParameterWithDefaultValue(Restraint_file, "Restraint_file",
                                            "idealCrystal.in");
    DefineOptionalParameterWithDefaultValue(
        UseThermodynamicIntegration, "useThermodynamicIntegration", false);
    DefineOptionalParameterWithDefaultValue(HULL_Method, "HULL_Method",
                                            "Convex");

    DefineOptionalParameterWithDefaultValue(PrivilegedAxis, "privilegedAxis",
                                            "z");

    deprecatedKeywords_.insert("nComponents");
    deprecatedKeywords_.insert("nZconstraints");
    deprecatedKeywords_.insert("initialConfig");
    deprecatedKeywords_.insert("thermIntDistSpringConst");
    deprecatedKeywords_.insert("thermIntThetaSpringConst");
    deprecatedKeywords_.insert("thermIntOmegaSpringConst");
    deprecatedKeywords_.insert("useSolidThermInt");
    deprecatedKeywords_.insert("useLiquidThermInt");
    deprecatedKeywords_.insert("minimizerMaxIter");
    deprecatedKeywords_.insert("minimizerWriteFreq");
    deprecatedKeywords_.insert("minimizerStepSize");
    deprecatedKeywords_.insert("minimizerFTol");
    deprecatedKeywords_.insert("minimizerGTol");
    deprecatedKeywords_.insert("minimizerLSTol");
    deprecatedKeywords_.insert("minimizerLSMaxIter");
    deprecatedKeywords_.insert("electricField");
    deprecatedKeywords_.insert("cutoffPolicy");
    deprecatedKeywords_.insert("beadSize");
  }

  Globals::~Globals() {
    Utils::deletePointers(moleculeStamps_);

    Utils::deletePointers(components_);
    Utils::deletePointers(zconstraints_);
    Utils::deletePointers(restraints_);

    delete flucQpars_;
    delete rnemdPars_;
    delete lightPars_;
    delete minimizerPars_;
  }

  void Globals::validate() {
    DataHolder::validate();

    CheckParameter(ForceField, isNotEmpty());
    CheckParameter(TargetTemp, isPositive());
    CheckParameter(
        Ensemble,
        isEqualIgnoreCase("NVE") || isEqualIgnoreCase("NVT") ||
            isEqualIgnoreCase("NPTi") || isEqualIgnoreCase("NPTf") ||
            isEqualIgnoreCase("NPTxyz") || isEqualIgnoreCase("NPTsz") ||
            isEqualIgnoreCase("NPAT") || isEqualIgnoreCase("NPA") ||
            isEqualIgnoreCase("LANGEVINDYNAMICS") || isEqualIgnoreCase("LD") ||
            isEqualIgnoreCase("NPRT") || isEqualIgnoreCase("NPGT") ||
            isEqualIgnoreCase("NGammaT") || isEqualIgnoreCase("NGT") ||
            isEqualIgnoreCase("LANGEVINHULL") || isEqualIgnoreCase("LHULL") ||
            isEqualIgnoreCase("SMIPD") || isEqualIgnoreCase("LANGEVINPISTON") ||
            isEqualIgnoreCase("SPF"));
    CheckParameter(Dt, isPositive());
    CheckParameter(RunTime, isPositive());
    CheckParameter(FinalConfig, isNotEmpty());
    CheckParameter(SampleTime, isNonNegative());
    CheckParameter(ResetTime, isNonNegative());
    CheckParameter(StatusTime, isNonNegative());
    CheckParameter(CutoffRadius, isPositive());
    CheckParameter(SwitchingRadius, isNonNegative());
    CheckParameter(Dielectric, isPositive());
    CheckParameter(ThermalTime, isNonNegative());
    CheckParameter(TauThermostat, isPositive());
    CheckParameter(TauBarostat, isPositive());
    CheckParameter(ZconsTime, isPositive());
    CheckParameter(ZconsTol, isPositive());
    CheckParameter(Seed, isPositive());
    CheckParameter(ZconsGap, isPositive());
    CheckParameter(ZconsFixtime, isPositive());
    CheckParameter(ThermodynamicIntegrationLambda, isNonNegative());
    CheckParameter(ThermodynamicIntegrationK, isPositive());
    CheckParameter(ForceFieldVariant, isNotEmpty());
    CheckParameter(ForceFieldFileName, isNotEmpty());
    CheckParameter(CutoffMethod, isEqualIgnoreCase("HARD") ||
                                     isEqualIgnoreCase("SWITCHED") ||
                                     isEqualIgnoreCase("SHIFTED_POTENTIAL") ||
                                     isEqualIgnoreCase("SHIFTED_FORCE") ||
                                     isEqualIgnoreCase("TAYLOR_SHIFTED") ||
                                     isEqualIgnoreCase("EWALD_FULL"));
    CheckParameter(ElectrostaticSummationMethod,
                   isEqualIgnoreCase("NONE") || isEqualIgnoreCase("HARD") ||
                       isEqualIgnoreCase("SWITCHED") ||
                       isEqualIgnoreCase("SHIFTED_POTENTIAL") ||
                       isEqualIgnoreCase("SHIFTED_FORCE") ||
                       isEqualIgnoreCase("REACTION_FIELD") ||
                       isEqualIgnoreCase("TAYLOR_SHIFTED") ||
                       isEqualIgnoreCase("EWALD_FULL"));
    CheckParameter(
        ElectrostaticScreeningMethod,
        isEqualIgnoreCase("UNDAMPED") || isEqualIgnoreCase("DAMPED"));
    CheckParameter(SwitchingFunctionType,
                   isEqualIgnoreCase("CUBIC") ||
                       isEqualIgnoreCase("FIFTH_ORDER_POLYNOMIAL"));
    CheckParameter(OrthoBoxTolerance, isPositive());
    CheckParameter(DampingAlpha, isNonNegative());
    CheckParameter(SkinThickness, isPositive());
    CheckParameter(Viscosity, isNonNegative());
    CheckParameter(BeadSize, isPositive());
    CheckParameter(FrozenBufferRadius, isPositive());
    CheckParameter(LangevinBufferRadius, isPositive());
    CheckParameter(NeighborListNeighbors, isPositive());
    CheckParameter(HULL_Method, isEqualIgnoreCase("Convex") ||
                                    isEqualIgnoreCase("AlphaShape"));
    CheckParameter(Alpha, isPositive());
    CheckParameter(StatFilePrecision, isPositive());
    CheckParameter(PrivilegedAxis, isEqualIgnoreCase("x") ||
                                       isEqualIgnoreCase("y") ||
                                       isEqualIgnoreCase("z"));

    for (std::vector<Component*>::iterator i = components_.begin();
         i != components_.end(); ++i) {
      if (!(*i)->findMoleculeStamp(moleculeStamps_)) {
        std::ostringstream oss;
        oss << "Globals Error: can not find molecule stamp for component "
            << (*i)->getType() << std::endl;
        throw OpenMDException(oss.str());
      }
    }
  }

  bool Globals::addComponent(Component* comp) {
    components_.push_back(comp);
    return true;
  }

  bool Globals::addZConsStamp(ZConsStamp* zcons) {
    zconstraints_.push_back(zcons);
    return true;
  }

  bool Globals::addRestraintStamp(RestraintStamp* rest) {
    restraints_.push_back(rest);
    return true;
  }

  bool Globals::addFluctuatingChargeParameters(
      FluctuatingChargeParameters* fqp) {
    if (flucQpars_ != NULL) delete flucQpars_;

    flucQpars_ = fqp;
    return true;
  }

  bool Globals::addRNEMDParameters(RNEMD::RNEMDParameters* rnemdPars) {
    if (rnemdPars_ != NULL) delete rnemdPars_;

    rnemdPars_ = rnemdPars;
    return true;
  }

  bool Globals::addLightParameters(Perturbations::LightParameters* lightPars) {
    if (lightPars_ != NULL) delete lightPars_;

    lightPars_ = lightPars;
    return true;
  }

  bool Globals::addMinimizerParameters(MinimizerParameters* miniPars) {
    if (minimizerPars_ != NULL) delete minimizerPars_;

    minimizerPars_ = miniPars;
    return true;
  }

  bool Globals::addMoleculeStamp(MoleculeStamp* molStamp) {
    std::string molStampName = molStamp->getName();
    std::map<std::string, MoleculeStamp*>::iterator i;
    bool ret = false;
    i        = moleculeStamps_.find(molStampName);
    if (i == moleculeStamps_.end()) {
      moleculeStamps_.insert(std::map<std::string, MoleculeStamp*>::value_type(
          molStampName, molStamp));
      ret = true;
    } else {
      std::ostringstream oss;
      oss << "Globals Error: Molecule Stamp " << molStamp->getName()
          << "appears multiple times\n";
      throw OpenMDException(oss.str());
    }
    return ret;
  }

  bool Globals::addFragmentStamp(FragmentStamp* fragStamp) {
    std::string fragStampName = fragStamp->getName();
    std::map<std::string, FragmentStamp*>::iterator i;
    bool ret = false;
    i        = fragmentStamps_.find(fragStampName);
    if (i == fragmentStamps_.end()) {
      fragmentStamps_.insert(std::map<std::string, FragmentStamp*>::value_type(
          fragStampName, fragStamp));
      ret = true;
    } else {
      std::ostringstream oss;
      oss << "Globals Error: Fragment Stamp " << fragStamp->getName()
          << "appears multiple times\n";
      throw OpenMDException(oss.str());
    }
    return ret;
  }
}  // namespace OpenMD
