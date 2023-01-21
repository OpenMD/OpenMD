/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#ifndef IO_GLOBALS_HPP
#define IO_GLOBALS_HPP

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "flucq/FluctuatingChargeParameters.hpp"
#include "optimization/MinimizerParameters.hpp"
#include "rnemd/RNEMDParameters.hpp"
#include "types/Component.hpp"
#include "types/FragmentStamp.hpp"
#include "types/MoleculeStamp.hpp"
#include "types/RestraintStamp.hpp"
#include "types/ZconsStamp.hpp"
#include "utils/ParameterManager.hpp"

namespace OpenMD {
  class Globals : public DataHolder {
  public:
    using intPair = std::pair<int, int>;

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
    DeclareParameter(UseLongRangeCorrections, bool);
    DeclareParameter(TauThermostat, RealType);
    DeclareParameter(TauBarostat, RealType);
    DeclareParameter(LangevinPistonDrag, RealType);
    DeclareParameter(ZconsTime, RealType);
    DeclareParameter(ZconsTol, RealType);
    DeclareParameter(ZconsForcePolicy, std::string);
    DeclareParameter(Seed, unsigned long int);
    DeclareParameter(UseInitalTime, bool);
    DeclareParameter(UseIntialExtendedSystemState, bool);
    DeclareParameter(OrthoBoxTolerance, RealType);
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
    DeclareParameter(PrintVirialTensor, bool);
    DeclareParameter(PrintHeatFlux, bool);
    DeclareParameter(TaggedAtomPair, intPair);
    DeclareParameter(PrintTaggedPairDistance, bool);
    DeclareParameter(ElectrostaticSummationMethod, std::string);
    DeclareParameter(ElectrostaticScreeningMethod, std::string);
    DeclareParameter(UseSurfaceTerm, bool);
    DeclareParameter(UseSlabGeometry, bool);
    DeclareParameter(DampingAlpha, RealType);
    DeclareParameter(Dielectric, RealType);
    DeclareParameter(CutoffMethod, std::string);
    DeclareParameter(SwitchingFunctionType, std::string);
    DeclareParameter(CompressDumpFile, bool);
    DeclareParameter(OutputForceVector, bool);
    DeclareParameter(OutputParticlePotential, bool);
    DeclareParameter(OutputElectricField, bool);
    DeclareParameter(OutputFluctuatingCharges, bool);
    DeclareParameter(OutputSitePotential, bool);
    DeclareParameter(OutputDensity, bool);
    DeclareParameter(SkinThickness, RealType);
    DeclareParameter(StatFileFormat, std::string);
    DeclareParameter(StatFilePrecision, int);
    DeclareParameter(HydroPropFile, std::string);
    DeclareParameter(Viscosity, RealType);
    DeclareParameter(BeadSize, RealType);
    DeclareParameter(UseSphericalBoundaryConditions, bool);
    DeclareParameter(FrozenBufferRadius, RealType);
    DeclareParameter(LangevinBufferRadius, RealType);
    DeclareParameter(AccumulateBoxDipole, bool);
    DeclareParameter(AccumulateBoxQuadrupole, bool);
    DeclareParameter(NeighborListNeighbors, int);
    DeclareParameter(UseMultipleTemperatureMethod, bool);
    DeclareParameter(MTM_Ce, RealType);
    DeclareParameter(MTM_G, RealType);
    DeclareParameter(MTM_Io, RealType);
    DeclareParameter(MTM_Sigma, RealType);
    DeclareParameter(MTM_R, RealType);
    DeclareParameter(UseRestraints, bool);
    DeclareParameter(Restraint_file, std::string);
    DeclareParameter(HULL_Method, std::string);
    DeclareParameter(Alpha, RealType);
    DeclareAlterableParameter(MDfileVersion, int);
    DeclareParameter(UniformField, std::vector<RealType>);
    DeclareParameter(MagneticField, std::vector<RealType>);
    DeclareParameter(UniformGradientStrength, RealType);
    DeclareParameter(UniformGradientDirection1, std::vector<RealType>);
    DeclareParameter(UniformGradientDirection2, std::vector<RealType>);
    DeclareParameter(ElectricField, std::vector<RealType>);
    DeclareParameter(ConstraintTime, RealType);
    DeclareParameter(PotentialSelection, std::string);
    DeclareParameter(PrivilegedAxis, std::string);

  public:
    bool addComponent(Component* comp);
    bool addZConsStamp(ZConsStamp* zcons);
    bool addRestraintStamp(RestraintStamp* rest);
    bool addMoleculeStamp(MoleculeStamp* molStamp);
    bool addFragmentStamp(FragmentStamp* fragStamp);
    size_t getNComponents() { return components_.size(); }
    std::vector<Component*> getComponents() { return components_; }
    Component* getComponentAt(int index) { return components_.at(index); }

    size_t getNZconsStamps() { return zconstraints_.size(); }
    std::vector<ZConsStamp*> getZconsStamps() { return zconstraints_; }
    ZConsStamp* getZconsStampAt(int index) { return zconstraints_.at(index); }

    size_t getNRestraintStamps() { return restraints_.size(); }
    std::vector<RestraintStamp*> getRestraintStamps() { return restraints_; }
    RestraintStamp* getRestraintStampAt(int index) {
      return restraints_.at(index);
    }

    bool addFluctuatingChargeParameters(FluctuatingChargeParameters* flucqPars);
    FluctuatingChargeParameters* getFluctuatingChargeParameters() {
      return flucQpars_;
    }

    bool addRNEMDParameters(RNEMD::RNEMDParameters* rnemdPars);
    RNEMD::RNEMDParameters* getRNEMDParameters() { return rnemdPars_; }

    bool addMinimizerParameters(MinimizerParameters* miniPars);
    MinimizerParameters* getMinimizerParameters() { return minimizerPars_; }

    virtual void validate();

  private:
    std::vector<Component*> components_;
    std::vector<ZConsStamp*> zconstraints_;
    std::vector<RestraintStamp*> restraints_;
    std::map<std::string, MoleculeStamp*> moleculeStamps_;
    std::map<std::string, FragmentStamp*> fragmentStamps_;
    std::pair<int, int> taggedAtomPair_;
    FluctuatingChargeParameters* flucQpars_;
    RNEMD::RNEMDParameters* rnemdPars_;
    MinimizerParameters* minimizerPars_;
  };
}  // namespace OpenMD

#endif
