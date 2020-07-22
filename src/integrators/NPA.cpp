/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */
 
#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "integrators/IntegratorCreator.hpp"
#include "integrators/NPA.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {
  
  void NPA::moveA() {
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

    loadEta();
    
    instaTemp =thermo.getTemperature();
    press = thermo.getPressureTensor();
    instaPress = Constants::pressureConvert* (press(0, 0) + 
                                                      press(1, 1) + 
                                                      press(2, 2)) / 3.0;
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
        
	vel += dt2*Constants::energyConvert/mass* frc - dt2*sc;
	sd->setVel(vel);
        
	if (sd->isDirectional()) {
          
	  // get and convert the torque to body frame

	  Tb = sd->lab2Body(sd->getTrq());
          
	  // get the angular momentum, and propagate a half step
          
	  ji = sd->getJ();
          
	  ji += dt2*Constants::energyConvert * Tb 
            - dt2*thermostat.first* ji;
          
	  rotAlgo_->rotate(sd, ji, dt);
          
	  sd->setJ(ji);
	}            
      }
    }
    // evolve eta a half step
    
    evolveEtaA();    
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
    
    saveEta();
  }

  void NPA::moveB(void) {
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
    
    instaVol = thermo.getVolume();
    instaTemp = thermo.getTemperature();
    instaPress = thermo.getPressure();
    
    //evolve eta
    this->evolveEtaB();
    this->calcVelScale();

    index = 0;
    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {
      
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {            
        
        frc = sd->getFrc();
        mass = sd->getMass();

        getVelScaleB(sc, index);
        
        // velocity half step
        vel = oldVel[index] 
          + dt2*Constants::energyConvert/mass* frc 
          - dt2*sc;
        
        sd->setVel(vel);
        
        if (sd->isDirectional()) {
          // get and convert the torque to body frame
          Tb = sd->lab2Body(sd->getTrq());
          
          ji = oldJi[index] 
            + dt2*Constants::energyConvert*Tb 
            - dt2*thermostat.first*oldJi[index];
          
          sd->setJ(ji);
        }

        ++index;
      }
    }
        
    rattle_->constraintB();

    flucQ_->moveB();
    saveEta();
  }

  void NPA::evolveEtaA() {

    eta(2,2) += dt2 *  instaVol * (press(2, 2) - targetPressure/Constants::pressureConvert) / (NkBT*tb2);
    oldEta = eta;  
  }

  void NPA::evolveEtaB() {

    prevEta = eta;
    eta(2,2) = oldEta(2, 2) + dt2 *  instaVol *
	    (press(2, 2) - targetPressure/Constants::pressureConvert) / (NkBT*tb2);
  }

  void NPA::calcVelScale(){

    for (int i = 0; i < 3; i++ ) {
      for (int j = 0; j < 3; j++ ) {
	vScale(i, j) = eta(i, j);
      }
    }
  }

  void NPA::getVelScaleA(Vector3d& sc, const Vector3d& vel){
    sc = vScale * vel;
  }

  void NPA::getVelScaleB(Vector3d& sc, int index ) {
    sc = vScale * oldVel[index];
  }

  void NPA::getPosScale(const Vector3d& pos, const Vector3d& COM, int index, 
                        Vector3d& sc) {

    Vector3d rj = (oldPos[index] + pos)/(RealType)2.0 -COM;
    sc = eta * rj;
  }

  void NPA::scaleSimBox(){
    Mat3x3d scaleMat;
    
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        scaleMat(i, j) = 0.0;
        if(i==j) {
          scaleMat(i, j) = 1.0;
        }
      }
    }
    
    scaleMat(2, 2) = exp(dt*eta(2, 2));
    Mat3x3d hmat = snap->getHmat();
    hmat = hmat *scaleMat;
    snap->setHmat(hmat);
  }

  bool NPA::etaConverged() {
    int i;
    RealType diffEta, sumEta;

    sumEta = 0;
    for(i = 0; i < 3; i++) {
      sumEta += pow(prevEta(i, i) - eta(i, i), 2);
    }
    
    diffEta = sqrt( sumEta / 3.0 );

    return ( diffEta <= etaTolerance );
  }

  RealType NPA::calcConservedQuantity(){

    thermostat = snap->getThermostat();
    loadEta();
    
    // We need NkBT a lot, so just set it here: This is the RAW number
    // of integrableObjects, so no subtraction or addition of constraints or
    // orientational degrees of freedom:
    NkBT = info_->getNGlobalIntegrableObjects()*Constants::kB *targetTemp;

    // fkBT is used because the thermostat operates on more degrees of freedom
    // than the barostat (when there are particles with orientational degrees
    // of freedom).  
    fkBT = info_->getNdf()*Constants::kB *targetTemp;    
    
    RealType conservedQuantity;
    RealType totalEnergy;
    RealType thermostat_kinetic;
    RealType thermostat_potential;
    RealType barostat_kinetic;
    RealType barostat_potential;
    RealType trEta;

    totalEnergy = thermo.getTotalEnergy();

    thermostat_kinetic = 0.0;
    thermostat_potential = 0.0;

    SquareMatrix<RealType, 3> tmp = eta.transpose() * eta;
    trEta = tmp.trace();
    
    barostat_kinetic = NkBT * tb2 * trEta /(2.0 * Constants::energyConvert);

    barostat_potential = (targetPressure * thermo.getVolume() / Constants::pressureConvert) /Constants::energyConvert;

    conservedQuantity = totalEnergy + thermostat_kinetic + thermostat_potential +
      barostat_kinetic + barostat_potential;

    return conservedQuantity;

  }

  void NPA::loadEta() {
    eta= snap->getBarostat();

    //if (!eta.isDiagonal()) {
    //    sprintf( painCave.errMsg,
    //             "NPA error: the diagonal elements of eta matrix are not the same or etaMat is not a diagonal matrix");
    //    painCave.isFatal = 1;
    //    simError();
    //}
  }

  void NPA::saveEta() {
    snap->setBarostat(eta);
  }

}

