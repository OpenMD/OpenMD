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
 
#include <math.h>
#include <iostream>

#ifdef IS_MPI
#include <mpi.h>
#endif //is_mpi

#include "brains/Thermo.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
#include "utils/PhysicalConstants.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "math/ConvexHull.hpp"
#include "math/AlphaHull.hpp"

using namespace std;
namespace OpenMD {

  RealType Thermo::getTranslationalKinetic() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasTranslationalKineticEnergy) {
      SimInfo::MoleculeIterator miter;
      vector<StuntDouble*>::iterator iiter;
      Molecule* mol;
      StuntDouble* sd;    
      Vector3d vel;
      RealType mass;
      RealType kinetic(0.0);
      
      for (mol = info_->beginMolecule(miter); mol != NULL; 
           mol = info_->nextMolecule(miter)) {
        
        for (sd = mol->beginIntegrableObject(iiter); sd != NULL; 
             sd = mol->nextIntegrableObject(iiter)) {
          
          mass = sd->getMass();
          vel = sd->getVel();
          
          kinetic += mass * (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
          
        }
      }
      
#ifdef IS_MPI
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &kinetic, 1, MPI::REALTYPE, 
                                MPI::SUM);
#endif
      
      kinetic = kinetic * 0.5 / PhysicalConstants::energyConvert;
      
      
      snap->setTranslationalKineticEnergy(kinetic);
    }
    return snap->getTranslationalKineticEnergy();
  }

  RealType Thermo::getRotationalKinetic() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasRotationalKineticEnergy) {
      SimInfo::MoleculeIterator miter;
      vector<StuntDouble*>::iterator iiter;
      Molecule* mol;
      StuntDouble* sd;    
      Vector3d angMom;
      Mat3x3d I;
      int i, j, k;
      RealType kinetic(0.0);
      
      for (mol = info_->beginMolecule(miter); mol != NULL; 
           mol = info_->nextMolecule(miter)) {
        
        for (sd = mol->beginIntegrableObject(iiter); sd != NULL; 
             sd = mol->nextIntegrableObject(iiter)) {
          
          if (sd->isDirectional()) {
            angMom = sd->getJ();
            I = sd->getI();
            
            if (sd->isLinear()) {
              i = sd->linearAxis();
              j = (i + 1) % 3;
              k = (i + 2) % 3;
              kinetic += angMom[j] * angMom[j] / I(j, j) 
                + angMom[k] * angMom[k] / I(k, k);
            } else {                        
              kinetic += angMom[0]*angMom[0]/I(0, 0) 
                + angMom[1]*angMom[1]/I(1, 1) 
                + angMom[2]*angMom[2]/I(2, 2);
            }
          }           
        }
      }
      
#ifdef IS_MPI
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &kinetic, 1, MPI::REALTYPE, 
                                MPI::SUM);
#endif
      
      kinetic = kinetic * 0.5 / PhysicalConstants::energyConvert;
           
      snap->setRotationalKineticEnergy(kinetic);
    }
    return snap->getRotationalKineticEnergy();
  }

      

  RealType Thermo::getKinetic() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasKineticEnergy) {
      RealType ke = getTranslationalKinetic() + getRotationalKinetic();
      snap->setKineticEnergy(ke);
    }
    return snap->getKineticEnergy();
  }

  RealType Thermo::getPotential() {

    // ForceManager computes the potential and stores it in the
    // Snapshot.  All we have to do is report it.

    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    return snap->getPotentialEnergy();
  }

  RealType Thermo::getTotalEnergy() {

    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasTotalEnergy) {
      snap->setTotalEnergy(this->getKinetic() + this->getPotential());
    }

    return snap->getTotalEnergy();
  }

  RealType Thermo::getTemperature() {

    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasTemperature) {

      RealType temperature = ( 2.0 * this->getKinetic() ) 
        / (info_->getNdf()* PhysicalConstants::kb );

      snap->setTemperature(temperature);
    }
    
    return snap->getTemperature();
  }

  RealType Thermo::getElectronicTemperature() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasElectronicTemperature) {
      
      SimInfo::MoleculeIterator miter;
      vector<Atom*>::iterator iiter;
      Molecule* mol;
      Atom* atom;    
      RealType cvel;
      RealType cmass;
      RealType kinetic(0.0);
      RealType eTemp;
      
      for (mol = info_->beginMolecule(miter); mol != NULL; 
           mol = info_->nextMolecule(miter)) {
        
        for (atom = mol->beginFluctuatingCharge(iiter); atom != NULL; 
             atom = mol->nextFluctuatingCharge(iiter)) {
          
          cmass = atom->getChargeMass();
          cvel = atom->getFlucQVel();
          
          kinetic += cmass * cvel * cvel;
          
        }
      }
    
#ifdef IS_MPI
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &kinetic, 1, MPI::REALTYPE, 
                                MPI::SUM);
#endif

      kinetic *= 0.5;
      eTemp =  (2.0 * kinetic) / 
        (info_->getNFluctuatingCharges() * PhysicalConstants::kb );
     
      snap->setElectronicTemperature(eTemp);
    }

    return snap->getElectronicTemperature();
  }


  RealType Thermo::getVolume() { 
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    return snap->getVolume();
  }

  RealType Thermo::getPressure() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasPressure) {
      // Relies on the calculation of the full molecular pressure tensor
      
      Mat3x3d tensor;
      RealType pressure;
      
      tensor = getPressureTensor();
      
      pressure = PhysicalConstants::pressureConvert * 
        (tensor(0, 0) + tensor(1, 1) + tensor(2, 2)) / 3.0;
      
      snap->setPressure(pressure);
    }
    
    return snap->getPressure();    
  }

  Mat3x3d Thermo::getPressureTensor() {
    // returns pressure tensor in units amu*fs^-2*Ang^-1
    // routine derived via viral theorem description in:
    // Paci, E. and Marchi, M. J.Phys.Chem. 1996, 100, 4314-4322
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasPressureTensor) {

      Mat3x3d pressureTensor;
      Mat3x3d p_tens(0.0);
      RealType mass;
      Vector3d vcom;
      
      SimInfo::MoleculeIterator i;
      vector<StuntDouble*>::iterator j;
      Molecule* mol;
      StuntDouble* sd;    
      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {
        
        for (sd = mol->beginIntegrableObject(j); sd != NULL; 
             sd = mol->nextIntegrableObject(j)) {
          
          mass = sd->getMass();
          vcom = sd->getVel();
          p_tens += mass * outProduct(vcom, vcom);         
        }
      }
      
#ifdef IS_MPI
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, p_tens.getArrayPointer(), 9, 
                                MPI::REALTYPE, MPI::SUM);
#endif
      
      RealType volume = this->getVolume();
      Mat3x3d stressTensor = snap->getStressTensor();
      
      pressureTensor =  (p_tens + 
                         PhysicalConstants::energyConvert * stressTensor)/volume;
      
      snap->setPressureTensor(pressureTensor);
    }
    return snap->getPressureTensor();
  }




  Vector3d Thermo::getSystemDipole() {
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasSystemDipole) {
      SimInfo::MoleculeIterator miter;
      vector<Atom*>::iterator aiter;
      Molecule* mol;
      Atom* atom;
      RealType charge;
      RealType moment(0.0);
      Vector3d ri(0.0);
      Vector3d dipoleVector(0.0);
      Vector3d nPos(0.0);
      Vector3d pPos(0.0);
      RealType nChg(0.0);
      RealType pChg(0.0);
      int nCount = 0;
      int pCount = 0;
      
      RealType chargeToC = 1.60217733e-19;
      RealType angstromToM = 1.0e-10;
      RealType debyeToCm = 3.33564095198e-30;
      
      for (mol = info_->beginMolecule(miter); mol != NULL; 
           mol = info_->nextMolecule(miter)) {
        
        for (atom = mol->beginAtom(aiter); atom != NULL;
             atom = mol->nextAtom(aiter)) {
          
          charge = 0.0;
          
          FixedChargeAdapter fca = FixedChargeAdapter(atom->getAtomType());
          if ( fca.isFixedCharge() ) {
            charge = fca.getCharge();
          }
          
          FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atom->getAtomType());
          if ( fqa.isFluctuatingCharge() ) {
            charge += atom->getFlucQPos();
          }
          
          charge *= chargeToC;
          
          ri = atom->getPos();
          snap->wrapVector(ri);
          ri *= angstromToM;
          
          if (charge < 0.0) {
            nPos += ri;
            nChg -= charge;
            nCount++;
          } else if (charge > 0.0) {
            pPos += ri;
            pChg += charge;
            pCount++;
          }
          
          MultipoleAdapter ma = MultipoleAdapter(atom->getAtomType());
          if (ma.isDipole() ) {
            Vector3d u_i = atom->getElectroFrame().getColumn(2);
            moment = ma.getDipoleMoment();
            moment *= debyeToCm;
            dipoleVector += u_i * moment;
          }
        }
      }
      
      
#ifdef IS_MPI
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &pChg, 1, MPI::REALTYPE, 
                                MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &nChg, 1, MPI::REALTYPE, 
                                MPI::SUM);

      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &pCount, 1, MPI::INTEGER, 
                                MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &nCount, 1, MPI::INTEGER, 
                                MPI::SUM);

      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, pPos.getArrayPointer(), 3, 
                                MPI::REALTYPE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, nPos.getArrayPointer(), 3, 
                                MPI::REALTYPE, MPI::SUM);

      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, dipoleVector.getArrayPointer(),
                                3, MPI::REALTYPE, MPI::SUM);
#endif
      
      // first load the accumulated dipole moment (if dipoles were present)
      Vector3d boxDipole = dipoleVector;
      // now include the dipole moment due to charges
      // use the lesser of the positive and negative charge totals
      RealType chg_value = nChg <= pChg ? nChg : pChg;
      
      // find the average positions
      if (pCount > 0 && nCount > 0 ) {
        pPos /= pCount;
        nPos /= nCount;
      }
      
      // dipole is from the negative to the positive (physics notation)
      boxDipole += (pPos - nPos) * chg_value;
      snap->setSystemDipole(boxDipole);
    }

    return snap->getSystemDipole();
  }

  // Returns the Heat Flux Vector for the system
  Vector3d Thermo::getHeatFlux(){
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    SimInfo::MoleculeIterator miter;
    vector<StuntDouble*>::iterator iiter;
    Molecule* mol;
    StuntDouble* sd;    
    RigidBody::AtomIterator ai; 
    Atom* atom;       
    Vector3d vel;
    Vector3d angMom;
    Mat3x3d I;
    int i;
    int j;
    int k;
    RealType mass;

    Vector3d x_a;
    RealType kinetic;
    RealType potential;
    RealType eatom;
    RealType AvgE_a_ = 0;
    // Convective portion of the heat flux
    Vector3d heatFluxJc = V3Zero;

    /* Calculate convective portion of the heat flux */
    for (mol = info_->beginMolecule(miter); mol != NULL;
         mol = info_->nextMolecule(miter)) {
      
      for (sd = mol->beginIntegrableObject(iiter); 
           sd != NULL; 
	   sd = mol->nextIntegrableObject(iiter)) {
        
	mass = sd->getMass();
	vel = sd->getVel();

	kinetic = mass * (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
        
	if (sd->isDirectional()) {
	  angMom = sd->getJ();
	  I = sd->getI();

	  if (sd->isLinear()) {
	    i = sd->linearAxis();
	    j = (i + 1) % 3;
	    k = (i + 2) % 3;
	    kinetic += angMom[j] * angMom[j] / I(j, j) 
              + angMom[k] * angMom[k] / I(k, k);
	  } else {                        
	    kinetic += angMom[0]*angMom[0]/I(0, 0) 
              + angMom[1]*angMom[1]/I(1, 1) 
	      + angMom[2]*angMom[2]/I(2, 2);
	  }
	}

        potential = 0.0;

        if (sd->isRigidBody()) {
          RigidBody* rb = dynamic_cast<RigidBody*>(sd);
          for (atom = rb->beginAtom(ai); atom != NULL; 
               atom = rb->nextAtom(ai)) {
            potential +=  atom->getParticlePot();
          }          
        } else {
          potential = sd->getParticlePot();
        }

        potential *= PhysicalConstants::energyConvert; // amu A^2/fs^2
        // The potential may not be a 1/2 factor
        eatom = (kinetic + potential)/2.0;  // amu A^2/fs^2
        heatFluxJc[0] += eatom*vel[0]; // amu A^3/fs^3
        heatFluxJc[1] += eatom*vel[1]; // amu A^3/fs^3
        heatFluxJc[2] += eatom*vel[2]; // amu A^3/fs^3
      }
    }

    /* The J_v vector is reduced in the forceManager so everyone has
     *  the global Jv. Jc is computed over the local atoms and must be
     *  reduced among all processors.
     */
#ifdef IS_MPI
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &heatFluxJc[0], 3, MPI::REALTYPE, 
                              MPI::SUM);
#endif
    
    // (kcal/mol * A/fs) * conversion => (amu A^3)/fs^3

    Vector3d heatFluxJv = currSnapshot->getConductiveHeatFlux() * 
      PhysicalConstants::energyConvert;
        
    // Correct for the fact the flux is 1/V (Jc + Jv)
    return (heatFluxJv + heatFluxJc) / this->getVolume(); // amu / fs^3 
  }


  Vector3d Thermo::getComVel(){ 
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasCOMvel) {

      SimInfo::MoleculeIterator i;
      Molecule* mol;
      
      Vector3d comVel(0.0);
      RealType totalMass(0.0);
      
      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {
        RealType mass = mol->getMass();
        totalMass += mass;
        comVel += mass * mol->getComVel();
      }  
      
#ifdef IS_MPI
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &totalMass, 1, MPI::REALTYPE, 
                                MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, comVel.getArrayPointer(), 3, 
                                MPI::REALTYPE, MPI::SUM);
#endif
      
      comVel /= totalMass;
      snap->setCOMvel(comVel);
    }
    return snap->getCOMvel();
  }

  Vector3d Thermo::getCom(){ 
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!snap->hasCOM) {
      
      SimInfo::MoleculeIterator i;
      Molecule* mol;
      
      Vector3d com(0.0);
      RealType totalMass(0.0);
      
      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {
        RealType mass = mol->getMass();
        totalMass += mass;
        com += mass * mol->getCom();
      }  
      
#ifdef IS_MPI
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &totalMass, 1, MPI::REALTYPE, 
                                MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, com.getArrayPointer(), 3, 
                                MPI::REALTYPE, MPI::SUM);
#endif
      
      com /= totalMass;
      snap->setCOM(com);
    }
    return snap->getCOM();
  }        

  /** 
   * Returns center of mass and center of mass velocity in one
   * function call.
   */   
  void Thermo::getComAll(Vector3d &com, Vector3d &comVel){ 
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!(snap->hasCOM && snap->hasCOMvel)) {

      SimInfo::MoleculeIterator i;
      Molecule* mol;
      
      RealType totalMass(0.0);
      
      com = 0.0;
      comVel = 0.0;
      
      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {
        RealType mass = mol->getMass();
        totalMass += mass;
        com += mass * mol->getCom();
        comVel += mass * mol->getComVel();           
      }  
      
#ifdef IS_MPI
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &totalMass, 1, MPI::REALTYPE, 
                                MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, com.getArrayPointer(), 3, 
                                MPI::REALTYPE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, comVel.getArrayPointer(), 3, 
                                MPI::REALTYPE, MPI::SUM);
#endif
      
      com /= totalMass;
      comVel /= totalMass;
      snap->setCOM(com);
      snap->setCOMvel(comVel);
    }    
    com = snap->getCOM();
    comVel = snap->getCOMvel();
    return;
  }        
   
  /**
   * Return intertia tensor for entire system and angular momentum
   * Vector.
   *
   *
   *
   *    [  Ixx -Ixy  -Ixz ]
   * I =| -Iyx  Iyy  -Iyz |
   *    [ -Izx -Iyz   Izz ]
   */
  void Thermo::getInertiaTensor(Mat3x3d &inertiaTensor, 
                                Vector3d &angularMomentum){

    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    
    if (!(snap->hasInertiaTensor && snap->hasCOMw)) {
       
      RealType xx = 0.0;
      RealType yy = 0.0;
      RealType zz = 0.0;
      RealType xy = 0.0;
      RealType xz = 0.0;
      RealType yz = 0.0;
      Vector3d com(0.0);
      Vector3d comVel(0.0);
      
      getComAll(com, comVel);
      
      SimInfo::MoleculeIterator i;
      Molecule* mol;
      
      Vector3d thisq(0.0);
      Vector3d thisv(0.0);
      
      RealType thisMass = 0.0;
      
      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {
        
        thisq = mol->getCom()-com;
        thisv = mol->getComVel()-comVel;
        thisMass = mol->getMass();
        // Compute moment of intertia coefficients.
        xx += thisq[0]*thisq[0]*thisMass;
        yy += thisq[1]*thisq[1]*thisMass;
        zz += thisq[2]*thisq[2]*thisMass;
        
        // compute products of intertia
        xy += thisq[0]*thisq[1]*thisMass;
        xz += thisq[0]*thisq[2]*thisMass;
        yz += thisq[1]*thisq[2]*thisMass;
        
        angularMomentum += cross( thisq, thisv ) * thisMass;            
      }
      
      inertiaTensor(0,0) = yy + zz;
      inertiaTensor(0,1) = -xy;
      inertiaTensor(0,2) = -xz;
      inertiaTensor(1,0) = -xy;
      inertiaTensor(1,1) = xx + zz;
      inertiaTensor(1,2) = -yz;
      inertiaTensor(2,0) = -xz;
      inertiaTensor(2,1) = -yz;
      inertiaTensor(2,2) = xx + yy;
      
#ifdef IS_MPI
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, inertiaTensor.getArrayPointer(),
                                9, MPI::REALTYPE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, 
                                angularMomentum.getArrayPointer(), 3,
                                MPI::REALTYPE, MPI::SUM);
#endif
      
      snap->setCOMw(angularMomentum);
      snap->setInertiaTensor(inertiaTensor);
    }
    
    angularMomentum = snap->getCOMw();
    inertiaTensor = snap->getInertiaTensor();
    
    return;
  }

  // Returns the angular momentum of the system
  Vector3d Thermo::getAngularMomentum(){
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    
    if (!snap->hasCOMw) {
      
      Vector3d com(0.0);
      Vector3d comVel(0.0);
      Vector3d angularMomentum(0.0);
      
      getComAll(com, comVel);
      
      SimInfo::MoleculeIterator i;
      Molecule* mol;
      
      Vector3d thisr(0.0);
      Vector3d thisp(0.0);
      
      RealType thisMass;
      
      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {
        thisMass = mol->getMass(); 
        thisr = mol->getCom() - com;
        thisp = (mol->getComVel() - comVel) * thisMass;
        
        angularMomentum += cross( thisr, thisp );      
      }  
      
#ifdef IS_MPI
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, 
                                angularMomentum.getArrayPointer(), 3,
                                MPI::REALTYPE, MPI::SUM);
#endif
      
      snap->setCOMw(angularMomentum);
    }
    
    return snap->getCOMw();
  }
  
  
  /**
   * Returns the Volume of the system based on a ellipsoid with
   * semi-axes based on the radius of gyration V=4/3*Pi*R_1*R_2*R_3
   * where R_i are related to the principle inertia moments
   *  R_i = sqrt(C*I_i/N), this reduces to 
   *  V = 4/3*Pi*(C/N)^3/2*sqrt(det(I)). 
   * See S.E. Baltazar et. al. Comp. Mat. Sci. 37 (2006) 526-536.
   */
  RealType Thermo::getGyrationalVolume(){
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    
    if (!snap->hasGyrationalVolume) {
      
      Mat3x3d intTensor;
      RealType det;
      Vector3d dummyAngMom; 
      RealType sysconstants;
      RealType geomCnst;
      RealType volume;
      
      geomCnst = 3.0/2.0;
      /* Get the inertial tensor and angular momentum for free*/
      getInertiaTensor(intTensor, dummyAngMom);
      
      det = intTensor.determinant();
      sysconstants = geomCnst / (RealType)(info_->getNGlobalIntegrableObjects());
      volume = 4.0/3.0*NumericConstant::PI*pow(sysconstants,geomCnst)*sqrt(det);

      snap->setGyrationalVolume(volume);
    } 
    return snap->getGyrationalVolume();
  }
  
  void Thermo::getGyrationalVolume(RealType &volume, RealType &detI){
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (!(snap->hasInertiaTensor && snap->hasGyrationalVolume)) {
    
      Mat3x3d intTensor;
      Vector3d dummyAngMom; 
      RealType sysconstants;
      RealType geomCnst;
      
      geomCnst = 3.0/2.0;
      /* Get the inertia tensor and angular momentum for free*/
      this->getInertiaTensor(intTensor, dummyAngMom);
      
      detI = intTensor.determinant();
      sysconstants = geomCnst/(RealType)(info_->getNGlobalIntegrableObjects());
      volume = 4.0/3.0*NumericConstant::PI*pow(sysconstants,geomCnst)*sqrt(detI);
      snap->setGyrationalVolume(volume);
    } else {
      volume = snap->getGyrationalVolume();
      detI = snap->getInertiaTensor().determinant();
    }
    return;
  }
  
  RealType Thermo::getTaggedAtomPairDistance(){
    Snapshot* currSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    Globals* simParams = info_->getSimParams();
    
    if (simParams->haveTaggedAtomPair() && 
        simParams->havePrintTaggedPairDistance()) {
      if ( simParams->getPrintTaggedPairDistance()) {
        
        pair<int, int> tap = simParams->getTaggedAtomPair();
        Vector3d pos1, pos2, rab;
        
#ifdef IS_MPI        
	int mol1 = info_->getGlobalMolMembership(tap.first);
	int mol2 = info_->getGlobalMolMembership(tap.second);

        int proc1 = info_->getMolToProc(mol1);
        int proc2 = info_->getMolToProc(mol2);

	RealType data[3];
        if (proc1 == worldRank) {
          StuntDouble* sd1 = info_->getIOIndexToIntegrableObject(tap.first);
          pos1 = sd1->getPos();
          data[0] = pos1.x();
          data[1] = pos1.y();
          data[2] = pos1.z();          
          MPI_Bcast(data, 3, MPI_REALTYPE, proc1, MPI_COMM_WORLD);
        } else {
          MPI_Bcast(data, 3, MPI_REALTYPE, proc1, MPI_COMM_WORLD);
          pos1 = Vector3d(data);
        }

        if (proc2 == worldRank) {
          StuntDouble* sd2 = info_->getIOIndexToIntegrableObject(tap.second);
          pos2 = sd2->getPos();
          data[0] = pos2.x();
          data[1] = pos2.y();
          data[2] = pos2.z();          
          MPI_Bcast(data, 3, MPI_REALTYPE, proc2, MPI_COMM_WORLD);
        } else {
          MPI_Bcast(data, 3, MPI_REALTYPE, proc2, MPI_COMM_WORLD);
          pos2 = Vector3d(data);
        }
#else
        StuntDouble* at1 = info_->getIOIndexToIntegrableObject(tap.first);
        StuntDouble* at2 = info_->getIOIndexToIntegrableObject(tap.second);
        pos1 = at1->getPos();
        pos2 = at2->getPos();
#endif        
        rab = pos2 - pos1;
        currSnapshot->wrapVector(rab);
        return rab.length();
      }
      return 0.0;     
    }
    return 0.0;
  }

  RealType Thermo::getHullVolume(){
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    
    if (!snap->hasHullVolume) {

      Hull* surfaceMesh_;

      Globals* simParams = info_->getSimParams();
      const std::string ht = simParams->getHULL_Method();
      
      if (ht == "Convex") {
        surfaceMesh_ = new ConvexHull();
      } else if (ht == "AlphaShape") {
        surfaceMesh_ = new AlphaHull(simParams->getAlpha());
      } else {
        return 0.0;
      }
      
      // Build a vector of stunt doubles to determine if they are
      // surface atoms
      std::vector<StuntDouble*> localSites_;
      Molecule* mol;
      StuntDouble* sd;
      SimInfo::MoleculeIterator i;
      Molecule::IntegrableObjectIterator  j;
      
      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {          
        for (sd = mol->beginIntegrableObject(j); 
             sd != NULL;
             sd = mol->nextIntegrableObject(j)) {   
          localSites_.push_back(sd);
        }
      }   
      
      // Compute surface Mesh
      surfaceMesh_->computeHull(localSites_);
      snap->setHullVolume(surfaceMesh_->getVolume());
    }
    return snap->getHullVolume();
  }  
}
