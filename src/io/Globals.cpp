/*
 * Copyright (c) 2005, 2010 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *                                                                      
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "io/Globals.hpp"
#include "io/ParamConstraint.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {
  Globals::Globals() {
  
    flucQpars_ = new FluctuatingChargeParameters();
    rnemdPars_ = new RNEMDParameters();
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
    DefineOptionalParameter(ElectricField, "electricField");

    DefineOptionalParameter(TaggedAtomPair, "taggedAtomPair");
    DefineOptionalParameter(PrintTaggedPairDistance, "printTaggedPairDistance");
    DefineOptionalParameter(SwitchingFunctionType, "switchingFunctionType");
    DefineOptionalParameter(HydroPropFile, "HydroPropFile");
    DefineOptionalParameter(Viscosity, "viscosity");
    DefineOptionalParameter(BeadSize, "beadSize");
    DefineOptionalParameter(FrozenBufferRadius, "frozenBufferRadius");
    DefineOptionalParameter(LangevinBufferRadius, "langevinBufferRadius");
    DefineOptionalParameter(NeighborListNeighbors,"NeighborListNeighbors");
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

  
    DefineOptionalParameterWithDefaultValue(UsePeriodicBoundaryConditions, 
                                            "usePeriodicBoundaryConditions", 
                                            true);
    DefineOptionalParameterWithDefaultValue(UseAtomicVirial, "useAtomicVirial", 
                                            true);
    DefineOptionalParameterWithDefaultValue(UseLongRangeCorrections, 
                                            "useLongRangeCorrections", true);
    DefineOptionalParameterWithDefaultValue(UseInitalTime, "useInitialTime", 
                                            false);
    DefineOptionalParameterWithDefaultValue(UseIntialExtendedSystemState, 
                                            "useInitialExtendedSystemState", 
                                            false);
    DefineOptionalParameterWithDefaultValue(OrthoBoxTolerance, 
                                            "orthoBoxTolerance", 1E-6);  
    DefineOptionalParameterWithDefaultValue(CutoffMethod, "cutoffMethod", 
                                            "SHIFTED_FORCE");
    DefineOptionalParameterWithDefaultValue(ElectrostaticScreeningMethod, 
                                            "electrostaticScreeningMethod", 
                                            "DAMPED");
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
    DefineOptionalParameterWithDefaultValue(SkinThickness, "skinThickness", 
                                            1.0);
    DefineOptionalParameterWithDefaultValue(StatFileFormat, 
                                            "statFileFormat", 
                                            "TIME|TOTAL_ENERGY|POTENTIAL_ENERGY|KINETIC_ENERGY|TEMPERATURE|PRESSURE|VOLUME|CONSERVED_QUANTITY");    
    DefineOptionalParameterWithDefaultValue(UseSphericalBoundaryConditions, 
                                            "useSphericalBoundaryConditions", 
                                            false);
    DefineOptionalParameterWithDefaultValue(AccumulateBoxDipole, 
                                            "accumulateBoxDipole", false);
    DefineOptionalParameterWithDefaultValue(UseRestraints, "useRestraints", 
                                            false);
    DefineOptionalParameterWithDefaultValue(Restraint_file, "Restraint_file", 
                                            "idealCrystal.in");
    DefineOptionalParameterWithDefaultValue(UseThermodynamicIntegration, 
                                            "useThermodynamicIntegration", 
                                            false);
    DefineOptionalParameterWithDefaultValue(HULL_Method,"HULL_Method","Convex");

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

    
  }

  Globals::~Globals() {
    MemoryUtils::deletePointers(components_);
    MemoryUtils::deletePointers(zconstraints_);
    MemoryUtils::deletePointers(restraints_);
  }

  void Globals::validate() {
    DataHolder::validate();

    CheckParameter(ForceField, isNotEmpty());
    CheckParameter(TargetTemp, isPositive());
    CheckParameter(Ensemble, isEqualIgnoreCase("NVE") || 
                   isEqualIgnoreCase("NVT") || isEqualIgnoreCase("NPTi") || 
                   isEqualIgnoreCase("NPTf") || isEqualIgnoreCase("NPTxyz") || 
                   isEqualIgnoreCase("NPTsz") || isEqualIgnoreCase("NPAT")  || 
                   isEqualIgnoreCase("LANGEVINDYNAMICS") || 
                   isEqualIgnoreCase("LD") || isEqualIgnoreCase("NPRT") || 
                   isEqualIgnoreCase("NPGT") || isEqualIgnoreCase("NGammaT") || 
                   isEqualIgnoreCase("NGT") || 
                   isEqualIgnoreCase("LANGEVINHULL") || 
                   isEqualIgnoreCase("LHULL") || isEqualIgnoreCase("SMIPD"));
    CheckParameter(Dt, isPositive());
    CheckParameter(RunTime, isPositive());
    CheckParameter(FinalConfig, isNotEmpty());
    CheckParameter(SampleTime, isNonNegative());
    CheckParameter(ResetTime, isNonNegative());
    CheckParameter(StatusTime, isNonNegative());
    CheckParameter(CutoffRadius, isPositive());
    CheckParameter(SwitchingRadius, isNonNegative());
    CheckParameter(Dielectric, isPositive());
    CheckParameter(ThermalTime,  isNonNegative());
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
                   isEqualIgnoreCase("TAYLOR_SHIFTED"));
    CheckParameter(CutoffPolicy, isEqualIgnoreCase("MIX") || 
                   isEqualIgnoreCase("MAX") || 
                   isEqualIgnoreCase("TRADITIONAL"));
    CheckParameter(ElectrostaticSummationMethod, isEqualIgnoreCase("NONE") || 
                   isEqualIgnoreCase("HARD") || isEqualIgnoreCase("SWITCHED") || 
                   isEqualIgnoreCase("SHIFTED_POTENTIAL") || 
                   isEqualIgnoreCase("SHIFTED_FORCE") || 
                   isEqualIgnoreCase("REACTION_FIELD") || 
                   isEqualIgnoreCase("TAYLOR_SHIFTED"));
    CheckParameter(ElectrostaticScreeningMethod, 
                   isEqualIgnoreCase("UNDAMPED") || isEqualIgnoreCase("DAMPED")); 
    CheckParameter(SwitchingFunctionType, isEqualIgnoreCase("CUBIC") || 
                   isEqualIgnoreCase("FIFTH_ORDER_POLYNOMIAL"));
    CheckParameter(OrthoBoxTolerance, isPositive());  
    CheckParameter(DampingAlpha,isNonNegative());
    CheckParameter(SkinThickness, isPositive());
    CheckParameter(Viscosity, isNonNegative());
    CheckParameter(BeadSize, isPositive());
    CheckParameter(FrozenBufferRadius, isPositive());
    CheckParameter(LangevinBufferRadius, isPositive());
    CheckParameter(NeighborListNeighbors, isPositive());
    CheckParameter(HULL_Method, isEqualIgnoreCase("Convex") || 
                   isEqualIgnoreCase("AlphaShape")); 
    CheckParameter(Alpha, isPositive()); 
  
    for(std::vector<Component*>::iterator i = components_.begin(); 
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

  bool Globals::addFluctuatingChargeParameters(FluctuatingChargeParameters* fqp) {
    if (flucQpars_ != NULL)
      delete flucQpars_;
    
    flucQpars_ = fqp;
    return true;
  }

  bool Globals::addRNEMDParameters(RNEMDParameters* rnemdPars) {
    if (rnemdPars_ != NULL)
      delete rnemdPars_;
    
    rnemdPars_ = rnemdPars;
    return true;
  }

  bool Globals::addMinimizerParameters(MinimizerParameters* miniPars) {
    if (minimizerPars_ != NULL)
      delete minimizerPars_;
    
    minimizerPars_ = miniPars;
    return true;
  }

  bool Globals::addMoleculeStamp(MoleculeStamp* molStamp) {
    std::string molStampName = molStamp->getName();
    std::map<std::string, MoleculeStamp*>::iterator i;
    bool ret = false;
    i = moleculeStamps_.find(molStampName);
    if (i == moleculeStamps_.end()) {
      moleculeStamps_.insert(std::map<std::string, MoleculeStamp*>::value_type(molStampName, molStamp));
      ret = true;
    } else {
      std::ostringstream oss;
      oss << "Globals Error: Molecule Stamp " << molStamp->getName() 
          << "appears multiple times\n";
      throw OpenMDException(oss.str());  
    }
    return ret;
  }
}
