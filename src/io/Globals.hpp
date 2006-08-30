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
 
#ifndef IO_GLOBALS_HPP
#define IO_GLOBALS_HPP

#include <iostream>

#include <stdlib.h>
#include <vector>
#include <string>
#include <map>

#include "types/Component.hpp"
#include "types/ZconsStamp.hpp"
#include "types/MoleculeStamp.hpp"
#include "utils/ParameterManager.hpp"

namespace oopse {
class Globals : public DataHolder {
  public:
    Globals();
    virtual ~Globals();
    
  DeclareParameter(ForceField, std::string);
  DeclareParameter(TargetTemp, RealType);
  DeclareParameter(Ensemble, std::string);
  DeclareParameter(Dt, RealType);
  DeclareParameter(RunTime, RealType);
  //DeclareParameter(InitialConfig, std::string);
  DeclareParameter(FinalConfig, std::string);
  DeclareParameter(SampleTime, RealType);
  DeclareParameter(ResetTime, RealType);
  DeclareParameter(StatusTime, RealType);
  DeclareParameter(CutoffRadius, RealType);
  DeclareParameter(SwitchingRadius, RealType);
  DeclareParameter(Dielectric, RealType);
  DeclareParameter(TempSet, bool);
  DeclareParameter(ThermalTime, RealType);
  DeclareParameter(UsePeriodicBoundaryConditions, bool);
  DeclareParameter(TargetPressure, RealType);
  DeclareParameter(TauThermostat, RealType);
  DeclareParameter(TauBarostat, RealType);
  DeclareParameter(ZconsTime, RealType);
  DeclareParameter(ZconsTol, RealType);
  DeclareParameter(ZconsForcePolicy, std::string);
  DeclareParameter(Seed, int);
  DeclareParameter(UseInitalTime, bool);
  DeclareParameter(UseIntialExtendedSystemState, bool);
  DeclareParameter(OrthoBoxTolerance, RealType);
  DeclareParameter(Minimizer, std::string);
  DeclareParameter(MinimizerMaxIter, RealType);
  DeclareParameter(MinimizerWriteFrq, int);
  DeclareParameter(MinimizerStepSize, RealType);
  DeclareParameter(MinimizerFTol, RealType);
  DeclareParameter(MinimizerGTol, RealType);
  DeclareParameter(MinimizerLSTol, RealType);
  DeclareParameter(MinimizerLSMaxIter, int);
  DeclareParameter(ZconsGap, RealType);
  DeclareParameter(ZconsFixtime, RealType);
  DeclareParameter(ZconsUsingSMD, bool);
  DeclareParameter(UseSolidThermInt, bool);
  DeclareParameter(UseLiquidThermInt, bool);
  DeclareParameter(ThermodynamicIntegrationLambda, RealType);
  DeclareParameter(ThermodynamicIntegrationK, RealType);
  DeclareParameter(ForceFieldVariant, std::string);
  DeclareParameter(ForceFieldFileName, std::string);
  DeclareParameter(ThermIntDistSpringConst, RealType);
  DeclareParameter(ThermIntThetaSpringConst, RealType);
  DeclareParameter(ThermIntOmegaSpringConst, RealType);
  DeclareParameter(SurfaceTension, RealType);
  DeclareParameter(PrintPressureTensor, bool);
  DeclareParameter(ElectrostaticSummationMethod, std::string);
  DeclareParameter(ElectrostaticScreeningMethod, std::string);
  DeclareParameter(DampingAlpha, RealType);
  DeclareParameter(CutoffPolicy, std::string);
  DeclareParameter(SwitchingFunctionType, std::string);
  DeclareParameter(CompressDumpFile, bool);
  DeclareParameter(OutputForceVector, bool);
  DeclareParameter(SkinThickness, RealType);
  DeclareParameter(StatFileFormat, std::string);    
  DeclareParameter(HydroPropFile, std::string);
  DeclareParameter(Viscosity, RealType);
  DeclareParameter(BeadSize, RealType);  
  DeclareParameter(UseSphericalBoundaryConditions, bool);
  DeclareParameter(FrozenBufferRadius, RealType);
  DeclareParameter(LangevinBufferRadius, RealType);
  DeclareParameter(AccumulateBoxDipole, bool);
  
  public:
    bool addComponent(Component* comp);
    bool addZConsStamp(ZConsStamp* zcons);
    bool addMoleculeStamp(MoleculeStamp* molStamp);
    int getNComponents() {return components_.size();}
    std::vector<Component*> getComponents() {return components_;}
    Component* getComponentAt(int index) {return components_.at(index);}    

    int getNZconsStamps() {return zconstraints_.size();}
    std::vector<ZConsStamp*> getZconsStamps() {return zconstraints_;}
    ZConsStamp* getZconsStampAt(int index) {return zconstraints_.at(index);}    

    virtual void validate();
  private:
    
    std::vector<Component*> components_;
    std::vector<ZConsStamp*> zconstraints_;    
    std::map<std::string, MoleculeStamp*> moleculeStamps_;

};
}
#endif

