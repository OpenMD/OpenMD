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
 
#include "integrators/NVT.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
#include "utils/PhysicalConstants.hpp"

namespace OpenMD {

  NVT::NVT(SimInfo* info) : VelocityVerletIntegrator(info), chiTolerance_ (1e-6), maxIterNum_(4) {

    Globals* simParams = info_->getSimParams();

    if (!simParams->getUseIntialExtendedSystemState()) {
      Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
      snap->setThermostat(make_pair(0.0, 0.0));
    }
    
    if (!simParams->haveTargetTemp()) {
      sprintf(painCave.errMsg, "You can't use the NVT integrator without a targetTemp_!\n");
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    } else {
      targetTemp_ = simParams->getTargetTemp();
    }

    // We must set tauThermostat.

    if (!simParams->haveTauThermostat()) {
      sprintf(painCave.errMsg, "If you use the constant temperature\n"
	      "\tintegrator, you must set tauThermostat.\n");

      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    } else {
      tauThermostat_ = simParams->getTauThermostat();
    }

    updateSizes();
  }

  void NVT::doUpdateSizes() {
    oldVel_.resize(info_->getNIntegrableObjects());
    oldJi_.resize(info_->getNIntegrableObjects());
  }

  void NVT::moveA() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* sd;
    Vector3d Tb;
    Vector3d ji;
    RealType mass;
    Vector3d vel;
    Vector3d pos;
    Vector3d frc;

    pair<RealType, RealType> thermostat = snap->getThermostat();

    // We need the temperature at time = t for the chi update below:

    RealType instTemp = thermo.getTemperature();

    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {

      for (sd = mol->beginIntegrableObject(j); sd != NULL;
	   sd = mol->nextIntegrableObject(j)) {

        vel = sd->getVel();
        pos = sd->getPos();
        frc = sd->getFrc();

        mass = sd->getMass();

        // velocity half step (use chi from previous step here):
        vel += dt2 *PhysicalConstants::energyConvert/mass*frc 
          - dt2*thermostat.first*vel;
        
        // position whole step
        pos += dt * vel;

        sd->setVel(vel);
        sd->setPos(pos);

        if (sd->isDirectional()) {

	  //convert the torque to body frame
	  Tb = sd->lab2Body(sd->getTrq());

	  // get the angular momentum, and propagate a half step

	  ji = sd->getJ();

	  ji += dt2*PhysicalConstants::energyConvert*Tb 
            - dt2*thermostat.first *ji;

	  rotAlgo_->rotate(sd, ji, dt);

	  sd->setJ(ji);
        }
      }

    }
    
    flucQ_->moveA();
    rattle_->constraintA();

    // Finally, evolve chi a half step (just like a velocity) using
    // temperature at time t, not time t+dt/2

    thermostat.first += dt2 * (instTemp / targetTemp_ - 1.0) 
      / (tauThermostat_ * tauThermostat_);
    thermostat.second += thermostat.first * dt2;

    snap->setThermostat(thermostat);
  }

  void NVT::moveB() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* sd;
    
    Vector3d Tb;
    Vector3d ji;    
    Vector3d vel;
    Vector3d frc;
    RealType mass;
    RealType instTemp;
    int index;
    // Set things up for the iteration:

    pair<RealType, RealType> thermostat = snap->getThermostat();
    RealType oldChi = thermostat.first;
    RealType  prevChi;

    index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {

      for (sd = mol->beginIntegrableObject(j); sd != NULL;
	   sd = mol->nextIntegrableObject(j)) {

	oldVel_[index] = sd->getVel();
        
        if (sd->isDirectional()) 
          oldJi_[index] = sd->getJ();                
        
	++index;    
      }           
    }

    // do the iteration:

    for(int k = 0; k < maxIterNum_; k++) {
      index = 0;
      instTemp = thermo.getTemperature();

      // evolve chi another half step using the temperature at t + dt/2

      prevChi = thermostat.first;
      thermostat.first = oldChi + dt2 * (instTemp / targetTemp_ - 1.0) 
        / (tauThermostat_ * tauThermostat_);

      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {
        
	for (sd = mol->beginIntegrableObject(j); sd != NULL;
	     sd = mol->nextIntegrableObject(j)) {

	  frc = sd->getFrc();
	  vel = sd->getVel();

	  mass = sd->getMass();

	  // velocity half step

	  vel = oldVel_[index] 
            + dt2/mass*PhysicalConstants::energyConvert * frc 
            - dt2*thermostat.first*oldVel_[index];
            
	  sd->setVel(vel);

	  if (sd->isDirectional()) {

	    // get and convert the torque to body frame

	    Tb =  sd->lab2Body(sd->getTrq());

	    ji = oldJi_[index] + dt2*PhysicalConstants::energyConvert*Tb 
              - dt2*thermostat.first *oldJi_[index];

	    sd->setJ(ji);
	  }


	  ++index;
	}
      }
    
      rattle_->constraintB();

      if (fabs(prevChi - thermostat.first) <= chiTolerance_)
	break;

    }

    flucQ_->moveB();

    thermostat.second += dt2 * thermostat.first;
    snap->setThermostat(thermostat);
  }

  void NVT::resetIntegrator() {
    snap->setThermostat(make_pair(0.0, 0.0));
  }
  
  RealType NVT::calcConservedQuantity() {

    pair<RealType, RealType> thermostat = snap->getThermostat();
    RealType conservedQuantity;
    RealType fkBT;
    RealType Energy;
    RealType thermostat_kinetic;
    RealType thermostat_potential;
    
    fkBT = info_->getNdf() *PhysicalConstants::kB *targetTemp_;

    Energy = thermo.getTotalEnergy();

    thermostat_kinetic = fkBT * tauThermostat_ * tauThermostat_ * thermostat.first * thermostat.first / (2.0 * PhysicalConstants::energyConvert);

    thermostat_potential = fkBT * thermostat.second / PhysicalConstants::energyConvert;

    conservedQuantity = Energy + thermostat_kinetic + thermostat_potential;

    return conservedQuantity;
  }


}//end namespace OpenMD
