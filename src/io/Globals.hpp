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
  DeclareParameter(TargetTemp, double);
  DeclareParameter(Ensemble, std::string);
  DeclareParameter(Dt, double);
  DeclareParameter(RunTime, double);
  DeclareParameter(InitialConfig, std::string);
  DeclareParameter(FinalConfig, std::string);
  DeclareParameter(SampleTime, double);
  DeclareParameter(ResetTime, double);
  DeclareParameter(StatusTime, double);
  DeclareParameter(CutoffRadius, double);
  DeclareParameter(SwitchingRadius, double);
  DeclareParameter(Dielectric, double);
  DeclareParameter(TempSet, bool);
  DeclareParameter(ThermalTime, double);
  DeclareParameter(UsePeriodicBoundaryConditions, bool);
  DeclareParameter(TargetPressure, double);
  DeclareParameter(TauThermostat, double);
  DeclareParameter(TauBarostat, double);
  DeclareParameter(ZconsTime, double);
  DeclareParameter(ZconsTol, double);
  DeclareParameter(ZconsForcePolicy, std::string);
  DeclareParameter(Seed, int);
  DeclareParameter(UseInitalTime, bool);
  DeclareParameter(UseIntialExtendedSystemState, bool);
  DeclareParameter(OrthoBoxTolerance, double);
  DeclareParameter(Minimizer, std::string);
  DeclareParameter(MinimizerMaxIter, double);
  DeclareParameter(MinimizerWriteFrq, int);
  DeclareParameter(MinimizerStepSize, double);
  DeclareParameter(MinimizerFTol, double);
  DeclareParameter(MinimizerGTol, double);
  DeclareParameter(MinimizerLSTol, double);
  DeclareParameter(MinimizerLSMaxIter, int);
  DeclareParameter(ZconsGap, double);
  DeclareParameter(ZconsFixtime, double);
  DeclareParameter(ZconsUsingSMD, bool);
  DeclareParameter(UseSolidThermInt, bool);
  DeclareParameter(UseLiquidThermInt, bool);
  DeclareParameter(ThermodynamicIntegrationLambda, double);
  DeclareParameter(ThermodynamicIntegrationK, double);
  DeclareParameter(ForceFieldVariant, std::string);
  DeclareParameter(ForceFieldFileName, std::string);
  DeclareParameter(ThermIntDistSpringConst, double);
  DeclareParameter(ThermIntThetaSpringConst, double);
  DeclareParameter(ThermIntOmegaSpringConst, double);
  DeclareParameter(SurfaceTension, double);
  DeclareParameter(PrintPressureTensor, bool);
  DeclareParameter(ElectrostaticSummationMethod, std::string);
  DeclareParameter(ElectrostaticScreeningMethod, std::string);
  DeclareParameter(DampingAlpha, double);
  DeclareParameter(CutoffPolicy, std::string);
  DeclareParameter(SwitchingFunctionType, std::string);
  DeclareParameter(CompressDumpFile, bool);
  DeclareParameter(OutputForceVector, bool);
  DeclareParameter(SkinThickness, double);
  DeclareParameter(StatFileFormat, std::string);    

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

