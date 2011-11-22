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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
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
#include "types/RestraintStamp.hpp"
#include "types/MoleculeStamp.hpp"
#include "utils/ParameterManager.hpp"

namespace OpenMD {
  class Globals : public DataHolder {
  public:
    typedef std::pair<int, int> intPair;

    Globals();
    virtual ~Globals();
    
    DeclareParameter(ForceField, std::string);
    DeclareParameter(TargetTemp, RealType);
    DeclareParameter(Ensemble, std::string);
    DeclareParameter(Dt, RealType);
    DeclareParameter(RunTime, RealType);
    DeclareParameter(FinalConfig, std::string);
    DeclareParameter(SampleTime, RealType);
    DeclareParameter(ResetTime, RealType);
    DeclareParameter(StatusTime, RealType);
    DeclareParameter(CutoffRadius, RealType);
    DeclareParameter(SwitchingRadius, RealType);
    DeclareParameter(TempSet, bool);
    DeclareParameter(ThermalTime, RealType);
    DeclareParameter(UsePeriodicBoundaryConditions, bool);
    DeclareParameter(TargetPressure, RealType);
    DeclareParameter(UseAtomicVirial, bool);
    DeclareParameter(TauThermostat, RealType);
    DeclareParameter(TauBarostat, RealType);
    DeclareParameter(ZconsTime, RealType);
    DeclareParameter(ZconsTol, RealType);
    DeclareParameter(ZconsForcePolicy, std::string);
    DeclareParameter(Seed, unsigned long int);
    DeclareParameter(UseInitalTime, bool);
    DeclareParameter(UseIntialExtendedSystemState, bool);
    DeclareParameter(OrthoBoxTolerance, RealType);
    DeclareParameter(Minimizer, std::string);
    DeclareParameter(MinimizerMaxIter, RealType);
    DeclareParameter(MinimizerWriteFreq, int);
    DeclareParameter(MinimizerStepSize, RealType);
    DeclareParameter(MinimizerFTol, RealType);
    DeclareParameter(MinimizerGTol, RealType);
    DeclareParameter(MinimizerLSTol, RealType);
    DeclareParameter(MinimizerLSMaxIter, int);
    DeclareParameter(ZconsGap, RealType);
    DeclareParameter(ZconsFixtime, RealType);
    DeclareParameter(ZconsUsingSMD, bool);
    DeclareParameter(UseThermodynamicIntegration, bool);
    DeclareParameter(ThermodynamicIntegrationLambda, RealType);
    DeclareParameter(ThermodynamicIntegrationK, RealType);
    DeclareParameter(ForceFieldVariant, std::string);
    DeclareParameter(ForceFieldFileName, std::string);
    DeclareParameter(SurfaceTension, RealType);
    DeclareParameter(PrintPressureTensor, bool);
    DeclareParameter(TaggedAtomPair, intPair);
    DeclareParameter(PrintTaggedPairDistance, bool);
    DeclareParameter(ElectrostaticSummationMethod, std::string);
    DeclareParameter(ElectrostaticScreeningMethod, std::string);
    DeclareParameter(DampingAlpha, RealType);
    DeclareParameter(Dielectric, RealType);
    DeclareParameter(CutoffMethod, std::string);
    DeclareParameter(CutoffPolicy, std::string);
    DeclareParameter(SwitchingFunctionType, std::string);
    DeclareParameter(CompressDumpFile, bool);
    DeclareParameter(OutputForceVector, bool);
    DeclareParameter(OutputParticlePotential, bool);
    DeclareParameter(SkinThickness, RealType);
    DeclareParameter(StatFileFormat, std::string);    
    DeclareParameter(HydroPropFile, std::string);
    DeclareParameter(Viscosity, RealType);
    DeclareParameter(BeadSize, RealType);  
    DeclareParameter(UseSphericalBoundaryConditions, bool);
    DeclareParameter(FrozenBufferRadius, RealType);
    DeclareParameter(LangevinBufferRadius, RealType);
    DeclareParameter(AccumulateBoxDipole, bool);
    DeclareParameter(NeighborListNeighbors, int);
    DeclareParameter(UseMultipleTemperatureMethod, bool);
    DeclareParameter(MTM_Ce, RealType);
    DeclareParameter(MTM_G, RealType);
    DeclareParameter(MTM_Io, RealType);
    DeclareParameter(MTM_Sigma, RealType);    
    DeclareParameter(MTM_R, RealType);    
    DeclareParameter(UseRNEMD, bool);
    DeclareParameter(RNEMD_exchangeTime, RealType);
    DeclareParameter(RNEMD_nBins, int);
    DeclareParameter(RNEMD_logWidth, int);
    DeclareParameter(RNEMD_exchangeType, std::string);
    DeclareParameter(RNEMD_objectSelection, std::string);
    DeclareParameter(RNEMD_targetFlux, RealType);
    DeclareParameter(RNEMD_binShift, bool);
    DeclareParameter(RNEMD_outputDimensionalTemperature, bool);
    DeclareParameter(UseRestraints, bool);
    DeclareParameter(Restraint_file, std::string);
    DeclareParameter(HULL_Method, std::string);
    DeclareParameter(Alpha, RealType);
    DeclareAlterableParameter(MDfileVersion, int);
    
  public:
    bool addComponent(Component* comp);
    bool addZConsStamp(ZConsStamp* zcons);
    bool addRestraintStamp(RestraintStamp* rest);
    bool addMoleculeStamp(MoleculeStamp* molStamp);
    int getNComponents() {return components_.size();}
    std::vector<Component*> getComponents() {return components_;}
    Component* getComponentAt(int index) {return components_.at(index);}    
    
    int getNZconsStamps() {return zconstraints_.size();}
    std::vector<ZConsStamp*> getZconsStamps() {return zconstraints_;}
    ZConsStamp* getZconsStampAt(int index) {return zconstraints_.at(index);}    

    int getNRestraintStamps() {return restraints_.size();}
    std::vector<RestraintStamp*> getRestraintStamps() {return restraints_;}
    RestraintStamp* getRestraintStampAt(int index) {return restraints_.at(index);}    
    
    virtual void validate();
  private:
    
    std::vector<Component*> components_;
    std::vector<ZConsStamp*> zconstraints_;    
    std::vector<RestraintStamp*> restraints_;    
    std::map<std::string, MoleculeStamp*> moleculeStamps_;
    std::pair<int, int> taggedAtomPair_;
};
}
#endif
