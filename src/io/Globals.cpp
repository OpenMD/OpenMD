/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "io/Globals.hpp"
#include "io/ParamConstraint.hpp"
#include "utils/MemoryUtils.hpp"
#include "utils/simError.h"

namespace oopse {
Globals::Globals() {
  DefineParameter(ForceField, "forceField")
 
  DefineOptionalParameter(TargetTemp, "targetTemp");
  DefineOptionalParameter(Ensemble, "ensemble");
  DefineOptionalParameter(Dt, "dt");
  DefineOptionalParameter(RunTime, "runTime");
  //DefineOptionalParameter(InitialConfig, "initialConfig");
  DefineOptionalParameter(FinalConfig, "finalConfig");
  DefineOptionalParameter(SampleTime, "sampleTime");
  DefineOptionalParameter(ResetTime, "resetTime");
  DefineOptionalParameter(StatusTime, "statusTime");
  DefineOptionalParameter(CutoffRadius, "cutoffRadius");
  DefineOptionalParameter(SwitchingRadius, "switchingRadius");
  DefineOptionalParameter(Dielectric, "dielectric");
  DefineOptionalParameter(TempSet, "tempSet");
  DefineOptionalParameter(ThermalTime, "thermalTime");
  DefineOptionalParameter(TargetPressure, "targetPressure");
  DefineOptionalParameter(TauThermostat, "tauThermostat");
  DefineOptionalParameter(TauBarostat, "tauBarostat");
  DefineOptionalParameter(ZconsTime, "zconsTime");
  DefineOptionalParameter(ZconsTol, "zconsTol");
  DefineOptionalParameter(ZconsForcePolicy, "zconsForcePolicy");
  DefineOptionalParameter(Seed, "seed");
  DefineOptionalParameter(Minimizer, "minimizer");
  DefineOptionalParameter(MinimizerMaxIter,"minimizerMaxIter");
  DefineOptionalParameter(MinimizerWriteFrq, "minimizerWriteFrq");
  DefineOptionalParameter(MinimizerStepSize, "minimizerStepSize");
  DefineOptionalParameter(MinimizerFTol, "minimizerFTol");
  DefineOptionalParameter(MinimizerGTol, "minimizerGTol");
  DefineOptionalParameter(MinimizerLSTol, "minimizerLSTol");
  DefineOptionalParameter(MinimizerLSMaxIter, "minimizerLSMaxIter");
  DefineOptionalParameter(ZconsGap, "zconsGap");
  DefineOptionalParameter(ZconsFixtime, "zconsFixtime");
  DefineOptionalParameter(ZconsUsingSMD, "zconsUsingSMD");
  DefineOptionalParameter(ThermodynamicIntegrationLambda, "thermodynamicIntegrationLambda");
  DefineOptionalParameter(ThermodynamicIntegrationK, "thermodynamicIntegrationK");
  DefineOptionalParameter(ForceFieldVariant, "forceFieldVariant");
  DefineOptionalParameter(ForceFieldFileName, "forceFieldFileName");
  DefineOptionalParameter(ThermIntDistSpringConst, "thermIntDistSpringConst");
  DefineOptionalParameter(ThermIntThetaSpringConst, "thermIntThetaSpringConst");
  DefineOptionalParameter(ThermIntOmegaSpringConst, "thermIntOmegaSpringConst");
  DefineOptionalParameter(SurfaceTension, "surfaceTension");
  DefineOptionalParameter(PrintPressureTensor, "printPressureTensor");
  DefineOptionalParameter(ElectrostaticSummationMethod, "electrostaticSummationMethod");
  DefineOptionalParameter(ElectrostaticScreeningMethod, "electrostaticScreeningMethod");
  DefineOptionalParameter(CutoffPolicy, "cutoffPolicy");
  DefineOptionalParameter(SwitchingFunctionType, "switchingFunctionType");
  DefineOptionalParameter(HydroPropFile, "HydroPropFile");
  DefineOptionalParameter(Viscosity, "viscosity");
  DefineOptionalParameter(BeadSize, "beadSize");
  DefineOptionalParameter(FrozenBufferRadius, "frozenBufferRadius");
  DefineOptionalParameter(LangevinBufferRadius, "langevinBufferRadius");
  
  DefineOptionalParameterWithDefaultValue(UsePeriodicBoundaryConditions, "usePeriodicBoundaryConditions", true);
  DefineOptionalParameterWithDefaultValue(UseInitalTime, "useInitialTime", false);
  DefineOptionalParameterWithDefaultValue(UseIntialExtendedSystemState, "useInitialExtendedSystemState", false);
  DefineOptionalParameterWithDefaultValue(OrthoBoxTolerance, "orthoBoxTolerance", 1E-6);  
  DefineOptionalParameterWithDefaultValue(UseSolidThermInt, "useSolidThermInt", false);
  DefineOptionalParameterWithDefaultValue(UseLiquidThermInt, "useLiquidThermInt", false);
  DefineOptionalParameterWithDefaultValue(ThermIntDistSpringConst, "thermIntDistSpringConst", 6.0);
  DefineOptionalParameterWithDefaultValue(ThermIntThetaSpringConst, "thermIntThetaSpringConst", 7.5);
  DefineOptionalParameterWithDefaultValue(ThermIntOmegaSpringConst, "thermIntOmegaSpringConst", 13.5);
  DefineOptionalParameter(DampingAlpha, "dampingAlpha");
  DefineOptionalParameterWithDefaultValue(CompressDumpFile, "compressDumpFile", 0);
  DefineOptionalParameterWithDefaultValue(OutputForceVector, "outputForceVector", 0);
  DefineOptionalParameterWithDefaultValue(SkinThickness, "skinThickness", 1.0);
  DefineOptionalParameterWithDefaultValue(StatFileFormat, "statFileFormat", "TIME|TOTAL_ENERGY|POTENTIAL_ENERGY|KINETIC_ENERGY|TEMPERATURE|PRESSURE|VOLUME|CONSERVED_QUANTITY");    
  DefineOptionalParameterWithDefaultValue(UseSphericalBoundaryConditions, "useSphericalBoundaryConditions", false);
  DefineOptionalParameterWithDefaultValue(AccumulateBoxDipole, "accumulateBoxDipole", false);


    deprecatedKeywords_.insert("nComponents");
    deprecatedKeywords_.insert("nZconstraints");
    deprecatedKeywords_.insert("initialConfig");
    
}

Globals::~Globals() {
    MemoryUtils::deletePointers(components_);
    MemoryUtils::deletePointers(zconstraints_);
}

void Globals::validate() {
  DataHolder::validate();

  CheckParameter(ForceField, isNotEmpty());
  CheckParameter(TargetTemp, isPositive());
  CheckParameter(Ensemble, isEqualIgnoreCase("NVE") || isEqualIgnoreCase("NVT") || isEqualIgnoreCase("NPTi") || isEqualIgnoreCase("NPTf") || isEqualIgnoreCase("NPTxyz") || isEqualIgnoreCase("NPAT")  || isEqualIgnoreCase("LANGEVINDYNAMICS") || isEqualIgnoreCase("LD") || isEqualIgnoreCase("NPRT") || isEqualIgnoreCase("NPGT") || isEqualIgnoreCase("NGammaT") || isEqualIgnoreCase("NGT"));
  CheckParameter(Dt, isPositive());
  CheckParameter(RunTime, isPositive());
  //CheckParameter(InitialConfig, isNotEmpty());
  CheckParameter(FinalConfig, isNotEmpty());
  CheckParameter(SampleTime, isNonNegative());
  CheckParameter(ResetTime, isNonNegative());
  CheckParameter(StatusTime, isNonNegative());
  CheckParameter(CutoffRadius, isPositive());
  CheckParameter(SwitchingRadius, isNonNegative());
  CheckParameter(Dielectric, isPositive());
  CheckParameter(ThermalTime,  isNonNegative());
  //  CheckParameter(TargetPressure,  isPositive());
  CheckParameter(TauThermostat, isPositive());
  CheckParameter(TauBarostat, isPositive());
  CheckParameter(ZconsTime, isPositive());
  CheckParameter(ZconsTol, isPositive());
  CheckParameter(Seed, isPositive());
  CheckParameter(Minimizer, isEqualIgnoreCase("SD") || isEqualIgnoreCase("CG"));
  CheckParameter(MinimizerMaxIter, isPositive());
  CheckParameter(MinimizerWriteFrq, isPositive());
  CheckParameter(MinimizerStepSize, isPositive());
  CheckParameter(MinimizerFTol, isPositive());
  CheckParameter(MinimizerGTol, isPositive());
  CheckParameter(MinimizerLSTol, isPositive());
  CheckParameter(MinimizerLSMaxIter, isPositive());
  CheckParameter(ZconsGap, isPositive());
  CheckParameter(ZconsFixtime, isPositive());
  CheckParameter(ThermodynamicIntegrationLambda, isNonNegative());
  CheckParameter(ThermodynamicIntegrationK, isPositive());
  CheckParameter(ForceFieldVariant, isNotEmpty());
  CheckParameter(ForceFieldFileName, isNotEmpty());
  CheckParameter(ThermIntDistSpringConst, isPositive());
  CheckParameter(ThermIntThetaSpringConst, isPositive());
  CheckParameter(ThermIntOmegaSpringConst, isPositive());
  //  CheckParameter(SurfaceTension, isNonNegative());
  CheckParameter(ElectrostaticSummationMethod, isEqualIgnoreCase("NONE") || isEqualIgnoreCase("SHIFTED_POTENTIAL") || isEqualIgnoreCase("SHIFTED_FORCE") || isEqualIgnoreCase("REACTION_FIELD"));
  CheckParameter(ElectrostaticScreeningMethod, isEqualIgnoreCase("UNDAMPED") || isEqualIgnoreCase("DAMPED")); 
  CheckParameter(CutoffPolicy, isEqualIgnoreCase("MIX") || isEqualIgnoreCase("MAX") || isEqualIgnoreCase("TRADITIONAL"));
  CheckParameter(SwitchingFunctionType, isEqualIgnoreCase("CUBIC") || isEqualIgnoreCase("FIFTH_ORDER_POLYNOMIAL"));
  //CheckParameter(StatFileFormat,);     
  CheckParameter(OrthoBoxTolerance, isPositive());  
  CheckParameter(ThermIntDistSpringConst, isPositive());
  CheckParameter(ThermIntThetaSpringConst, isPositive());
  CheckParameter(ThermIntOmegaSpringConst, isPositive());
  CheckParameter(DampingAlpha,isNonNegative());
  CheckParameter(SkinThickness, isPositive());
  CheckParameter(Viscosity, isNonNegative());
  CheckParameter(BeadSize, isPositive());
  CheckParameter(FrozenBufferRadius, isPositive());
  CheckParameter(LangevinBufferRadius, isPositive());
  for(std::vector<Component*>::iterator i = components_.begin(); i != components_.end(); ++i) {
    if (!(*i)->findMoleculeStamp(moleculeStamps_)) {
        std::ostringstream oss;
        oss << "Globals Error: can not find molecule stamp for component " << (*i)->getType() << std::endl;
        throw OOPSEException(oss.str());           
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
        oss << "Globals Error: Molecule Stamp " << molStamp->getName() << "appears multiple times\n";
        throw OOPSEException(oss.str());  
    }
    return ret;
}
 

}
