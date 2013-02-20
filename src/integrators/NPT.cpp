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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include <math.h>

#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "integrators/NPT.hpp"
#include "math/SquareMatrix3.hpp"
#include "primitives/Molecule.hpp"
#include "utils/PhysicalConstants.hpp"
#include "utils/simError.h"

// Basic isotropic thermostating and barostating via the Melchionna
// modification of the Hoover algorithm:
//
//    Melchionna, S., Ciccotti, G., and Holian, B. L., 1993,
//       Molec. Phys., 78, 533.
//
//           and
//
//    Hoover, W. G., 1986, Phys. Rev. A, 34, 2499.

namespace OpenMD {

  NPT::NPT(SimInfo* info) :
    VelocityVerletIntegrator(info), chiTolerance(1e-6), etaTolerance(1e-6), maxIterNum_(4) {

      Globals* simParams = info_->getSimParams();
    
      if (!simParams->getUseIntialExtendedSystemState()) {
        Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
        currSnapshot->setThermostat(make_pair(0.0, 0.0));
        currSnapshot->setBarostat(Mat3x3d(0.0));
      }
    
      if (!simParams->haveTargetTemp()) {
        sprintf(painCave.errMsg, "You can't use the NVT integrator without a targetTemp!\n");
        painCave.isFatal = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      } else {
        targetTemp = simParams->getTargetTemp();
      }

      // We must set tauThermostat
      if (!simParams->haveTauThermostat()) {
        sprintf(painCave.errMsg, "If you use the constant temperature\n"
		"\tintegrator, you must set tauThermostat.\n");

        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
      } else {
        tauThermostat = simParams->getTauThermostat();
      }

      if (!simParams->haveTargetPressure()) {
        sprintf(painCave.errMsg, "NPT error: You can't use the NPT integrator\n"
		"   without a targetPressure!\n");

        painCave.isFatal = 1;
        simError();
      } else {
        targetPressure = simParams->getTargetPressure();
      }
    
      if (!simParams->haveTauBarostat()) {
        sprintf(painCave.errMsg,
                "If you use the NPT integrator, you must set tauBarostat.\n");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
      } else {
        tauBarostat = simParams->getTauBarostat();
      }
    
      tt2 = tauThermostat * tauThermostat;
      tb2 = tauBarostat * tauBarostat;

      updateSizes();
    }

  NPT::~NPT() {
  }

  void NPT::doUpdateSizes() {

    oldPos.resize(info_->getNIntegrableObjects());
    oldVel.resize(info_->getNIntegrableObjects());
    oldJi.resize(info_->getNIntegrableObjects());

  }

  void NPT::moveA() {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* sd;
    Vector3d Tb, ji;
    RealType mass;
    Vector3d vel;
    Vector3d pos;
    Vector3d frc;
    Vector3d sc;
    int index;

    thermostat = snap->getThermostat();
    loadEta();
    
    instaTemp =thermo.getTemperature();
    press = thermo.getPressureTensor();
    instaPress = PhysicalConstants::pressureConvert* (press(0, 0) + press(1, 1) + press(2, 2)) / 3.0;
    instaVol =thermo.getVolume();

    Vector3d  COM = thermo.getCom();

    //evolve velocity half step

    calcVelScale();

    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {

      for (sd = mol->beginIntegrableObject(j); sd != NULL;
	   sd = mol->nextIntegrableObject(j)) {
                
	vel = sd->getVel();
	frc = sd->getFrc();

	mass = sd->getMass();

	getVelScaleA(sc, vel);

	// velocity half step  (use chi from previous step here):

	vel += dt2*PhysicalConstants::energyConvert/mass* frc - dt2*sc;
	sd->setVel(vel);

	if (sd->isDirectional()) {

	  // get and convert the torque to body frame

	  Tb = sd->lab2Body(sd->getTrq());

	  // get the angular momentum, and propagate a half step

	  ji = sd->getJ();

	  ji += dt2*PhysicalConstants::energyConvert * Tb 
            - dt2*thermostat.first* ji;
                
	  rotAlgo_->rotate(sd, ji, dt);

	  sd->setJ(ji);
	}
            
      }
    }
    // evolve chi and eta  half step

    thermostat.first += dt2 * (instaTemp / targetTemp - 1.0) / tt2;
    
    evolveEtaA();

    //calculate the integral of chidt
    thermostat.second += dt2 * thermostat.first;
    
    flucQ_->moveA();


    index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {

      for (sd = mol->beginIntegrableObject(j); sd != NULL;
	   sd = mol->nextIntegrableObject(j)) {

	oldPos[index++] = sd->getPos();            

      }
    }
    
    //the first estimation of r(t+dt) is equal to  r(t)

    for(int k = 0; k < maxIterNum_; k++) {
      index = 0;
      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {

	for (sd = mol->beginIntegrableObject(j); sd != NULL;
	     sd = mol->nextIntegrableObject(j)) {

	  vel = sd->getVel();
	  pos = sd->getPos();

	  this->getPosScale(pos, COM, index, sc);

	  pos = oldPos[index] + dt * (vel + sc);
	  sd->setPos(pos);     

	  ++index;
	}
      }

      rattle_->constraintA();
    }

    // Scale the box after all the positions have been moved:

    this->scaleSimBox();

    snap->setThermostat(thermostat);

    saveEta();
  }

  void NPT::moveB(void) {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator  j;
    Molecule* mol;
    StuntDouble* sd;
    int index;
    Vector3d Tb;
    Vector3d ji;
    Vector3d sc;
    Vector3d vel;
    Vector3d frc;
    RealType mass;

    thermostat = snap->getThermostat();
    RealType oldChi  = thermostat.first;
    RealType prevChi;

    loadEta();
    
    //save velocity and angular momentum
    index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {

      for (sd = mol->beginIntegrableObject(j); sd != NULL;
	   sd = mol->nextIntegrableObject(j)) {
                
	oldVel[index] = sd->getVel();

        if (sd->isDirectional())
	   oldJi[index] = sd->getJ();

	++index;
      }
    }

    // do the iteration:
    instaVol =thermo.getVolume();

    for(int k = 0; k < maxIterNum_; k++) {
      instaTemp =thermo.getTemperature();
      instaPress =thermo.getPressure();

      // evolve chi another half step using the temperature at t + dt/2
      prevChi = thermostat.first;
      thermostat.first = oldChi + dt2 * (instaTemp / targetTemp - 1.0) / tt2;

      //evolve eta
      this->evolveEtaB();
      this->calcVelScale();

      index = 0;
      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {

	for (sd = mol->beginIntegrableObject(j); sd != NULL;
	     sd = mol->nextIntegrableObject(j)) {            

	  frc = sd->getFrc();
	  vel = sd->getVel();

	  mass = sd->getMass();

	  getVelScaleB(sc, index);

	  // velocity half step
	  vel = oldVel[index] 
            + dt2*PhysicalConstants::energyConvert/mass* frc 
            - dt2*sc;

	  sd->setVel(vel);

	  if (sd->isDirectional()) {
	    // get and convert the torque to body frame
	    Tb = sd->lab2Body(sd->getTrq());

	    ji = oldJi[index] 
              + dt2*PhysicalConstants::energyConvert*Tb 
              - dt2*thermostat.first*oldJi[index];

	    sd->setJ(ji);
	  }

	  ++index;
	}
      }
        
      rattle_->constraintB();

      if ((fabs(prevChi - thermostat.first) <= chiTolerance) && 
          this->etaConverged())
	break;
    }

    //calculate integral of chidt
    thermostat.second += dt2 * thermostat.first;

    snap->setThermostat(thermostat);

    flucQ_->moveB();
    saveEta();
  }

  void NPT::resetIntegrator(){
    snap->setThermostat(make_pair(0.0, 0.0));
    resetEta();
  }

  void NPT::resetEta() {
    Mat3x3d etaMat(0.0);
    snap->setBarostat(etaMat);    
  }
}
