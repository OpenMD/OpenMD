/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 
#include "FluctuatingChargeNVT.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
#include "utils/PhysicalConstants.hpp"

namespace OpenMD {

  FluctuatingChargeNVT::FluctuatingChargeNVT(SimInfo* info) : 
    FluctuatingChargePropagator(info), chiTolerance_ (1e-6), 
    maxIterNum_(4), thermo(info), 
    currentSnapshot_(info->getSnapshotManager()->getCurrentSnapshot()) {
    
    if (hasFlucQ_) {     
      if (info_->getSimParams()->haveDt()) {
        dt_ = info_->getSimParams()->getDt();
        dt2_ = dt_ * 0.5;
      } else {
        sprintf(painCave.errMsg,
                "FluctuatingChargeNVT Error: dt is not set\n");
        painCave.isFatal = 1;
        simError();
      }
      
      if (!info_->getSimParams()->getUseIntialExtendedSystemState()) {
        currentSnapshot_->setChiElectronic(0.0);
        currentSnapshot_->setIntegralOfChiElectronicDt(0.0);
      }
      
      if (!fqParams_->haveTargetTemp()) {
        sprintf(painCave.errMsg, "You can't use the FluctuatingChargeNVT "
                "propagator without a flucQ.targetTemp!\n");
        painCave.isFatal = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      } else {
        targetTemp_ = fqParams_->getTargetTemp();
      }
      
      // We must set tauThermostat.
      
      if (!fqParams_->haveTauThermostat()) {
        sprintf(painCave.errMsg, "If you use the FluctuatingChargeNVT\n"
                "\tpropagator, you must set flucQ.tauThermostat .\n");
        
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
      } else {
        tauThermostat_ = fqParams_->getTauThermostat();
      }
      updateSizes();
    }
  }

  void FluctuatingChargeNVT::initialize() {
    FluctuatingChargePropagator::initialize();
  }


  void FluctuatingChargeNVT::moveA() {

    if (!hasFlucQ_) return;

    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator  j;
    Molecule* mol;
    Atom* atom;
    RealType cvel, cpos, cfrc, cmass;

    RealType chi = currentSnapshot_->getChiElectronic();
    RealType integralOfChidt = currentSnapshot_->getIntegralOfChiElectronicDt();
    RealType instTemp = thermo.getElectronicTemperature();


    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {
      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        
        cvel = atom->getFlucQVel();
        cpos = atom->getFlucQPos();
        cfrc = atom->getFlucQFrc();
        cmass = atom->getChargeMass();
        
	cerr << atom->getType() << "\n";
	cerr << "Before\n";
	cerr << "cvel: " << cvel << "\tcpos: " << cpos << "\tcfrc: " << cfrc << "\tcmass: " << cmass << "\n";
        
	// velocity half step
        cvel += dt2_ * cfrc / cmass - dt2_*chi*cvel;                    
        // position whole step
        cpos += dt_ * cvel;
	
	cerr << "After\n";
	cerr << "cvel: " << cvel << "\tcpos: " << cpos << "\n\n";        

        atom->setFlucQVel(cvel);
        atom->setFlucQPos(cpos);
      }
    }
    
    chi += dt2_ * (instTemp / targetTemp_ - 1.0) / 
      (tauThermostat_ * tauThermostat_);

    integralOfChidt += chi * dt2_;
    cerr << "Move A instTemp: " << instTemp << "\n";
    currentSnapshot_->setChiElectronic(chi);
    currentSnapshot_->setIntegralOfChiElectronicDt(integralOfChidt);
       
  }

  void FluctuatingChargeNVT::updateSizes() {
    if (!hasFlucQ_) return;
    oldVel_.resize(info_->getNFluctuatingCharges());
  }

  void FluctuatingChargeNVT::moveB() {
    if (!hasFlucQ_) return;
    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator  j;
    Molecule* mol;
    Atom* atom;
    RealType instTemp;
    RealType chi = currentSnapshot_->getChiElectronic();
    RealType oldChi = chi;
    RealType prevChi;
    RealType integralOfChidt = currentSnapshot_->getIntegralOfChiElectronicDt();
    int index;
    RealType cfrc, cvel, cmass;

    index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {
      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        
        oldVel_[index] = atom->getFlucQVel();
        ++index;
      }
    }
        
    // do the iteration:
    
    for(int k = 0; k < maxIterNum_; k++) {
      index = 0;
      instTemp = thermo.getElectronicTemperature();
      cerr << "MoveB instTemp: " << instTemp << "\n";
      // evolve chi another half step using the temperature at t + dt/2
      prevChi = chi;
      chi = oldChi + dt2_ * (instTemp / targetTemp_ - 1.0) / 
        (tauThermostat_ * tauThermostat_);

      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {
        for (atom = mol->beginFluctuatingCharge(j);  atom != NULL;
             atom = mol->nextFluctuatingCharge(j)) {

          cfrc = atom->getFlucQFrc();
          cvel =atom->getFlucQVel();
          cmass = atom->getChargeMass();
          
          // velocity half step
          cvel = oldVel_[index] + dt2_ * cfrc / cmass - dt2_*chi*oldVel_[index];
          
          atom->setFlucQVel(cvel);      
          ++index;          
        }
      }
      if (fabs(prevChi - chi) <= chiTolerance_)
        break;
    }
    integralOfChidt += dt2_ * chi;
    currentSnapshot_->setChiElectronic(chi);
    currentSnapshot_->setIntegralOfChiElectronicDt(integralOfChidt);
  }
  
  void FluctuatingChargeNVT::resetPropagator() {
    if (!hasFlucQ_) return;
    currentSnapshot_->setChiElectronic(0.0);
    currentSnapshot_->setIntegralOfChiElectronicDt(0.0);
  }
  
  RealType FluctuatingChargeNVT::calcConservedQuantity() {
    if (!hasFlucQ_) return 0.0;
    RealType chi = currentSnapshot_->getChiElectronic();
    RealType integralOfChidt = currentSnapshot_->getIntegralOfChiElectronicDt();
    RealType fkBT = info_->getNFluctuatingCharges() * 
      PhysicalConstants::kB *targetTemp_;

    RealType thermostat_kinetic = fkBT * tauThermostat_ * tauThermostat_ * 
      chi * chi / (2.0 * PhysicalConstants::energyConvert);

    RealType thermostat_potential = fkBT * integralOfChidt / 
      PhysicalConstants::energyConvert;

    return thermostat_kinetic + thermostat_potential;
  }
}
