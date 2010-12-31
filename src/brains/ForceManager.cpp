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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */
 
/**
 * @file ForceManager.cpp
 * @author tlin
 * @date 11/09/2004
 * @time 10:39am
 * @version 1.0
 */

#include "brains/ForceManager.hpp"
#include "primitives/Molecule.hpp"
#include "UseTheForce/doForces_interface.h"
#define __OPENMD_C
#include "UseTheForce/DarkSide/fInteractionMap.h"
#include "utils/simError.h"
#include "primitives/Bond.hpp"
#include "primitives/Bend.hpp"
#include "primitives/Torsion.hpp"
#include "primitives/Inversion.hpp"

namespace OpenMD {
  
  ForceManager::ForceManager(SimInfo * info) : info_(info), 
                                               NBforcesInitialized_(false) {
  }
 
  void ForceManager::calcForces() {
    

    if (!info_->isFortranInitialized()) {
      info_->update();
      nbiMan_->setSimInfo(info_);
      nbiMan_->initialize();    
      info_->setupFortran();
    }
    
    preCalculation();
    
    calcShortRangeInteraction();

    calcLongRangeInteraction();

    postCalculation();
    
  }
  
  void ForceManager::preCalculation() {
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::AtomIterator ai;
    Atom* atom;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;
    
    // forces are zeroed here, before any are accumulated.
    // NOTE: do not rezero the forces in Fortran.
    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
      for(atom = mol->beginAtom(ai); atom != NULL; atom = mol->nextAtom(ai)) {
	atom->zeroForcesAndTorques();
      }
          
      //change the positions of atoms which belong to the rigidbodies
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
	rb->zeroForcesAndTorques();
      }        
          
    }
    
    // Zero out the stress tensor
    tau *= 0.0;
    
  }
  
  void ForceManager::calcShortRangeInteraction() {
    Molecule* mol;
    RigidBody* rb;
    Bond* bond;
    Bend* bend;
    Torsion* torsion;
    Inversion* inversion;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Molecule::BondIterator bondIter;;
    Molecule::BendIterator  bendIter;
    Molecule::TorsionIterator  torsionIter;
    Molecule::InversionIterator  inversionIter;
    RealType bondPotential = 0.0;
    RealType bendPotential = 0.0;
    RealType torsionPotential = 0.0;
    RealType inversionPotential = 0.0;

    //calculate short range interactions    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {

      //change the positions of atoms which belong to the rigidbodies
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) {
        rb->updateAtoms();
      }

      for (bond = mol->beginBond(bondIter); bond != NULL; 
           bond = mol->nextBond(bondIter)) {
        bond->calcForce();
        bondPotential += bond->getPotential();
      }

      for (bend = mol->beginBend(bendIter); bend != NULL; 
           bend = mol->nextBend(bendIter)) {
        
        RealType angle;
        bend->calcForce(angle);
        RealType currBendPot = bend->getPotential();          
         
        bendPotential += bend->getPotential();
        std::map<Bend*, BendDataSet>::iterator i = bendDataSets.find(bend);
        if (i == bendDataSets.end()) {
          BendDataSet dataSet;
          dataSet.prev.angle = dataSet.curr.angle = angle;
          dataSet.prev.potential = dataSet.curr.potential = currBendPot;
          dataSet.deltaV = 0.0;
          bendDataSets.insert(std::map<Bend*, BendDataSet>::value_type(bend, dataSet));
        }else {
          i->second.prev.angle = i->second.curr.angle;
          i->second.prev.potential = i->second.curr.potential;
          i->second.curr.angle = angle;
          i->second.curr.potential = currBendPot;
          i->second.deltaV =  fabs(i->second.curr.potential -  
                                   i->second.prev.potential);
        }
      }
      
      for (torsion = mol->beginTorsion(torsionIter); torsion != NULL; 
           torsion = mol->nextTorsion(torsionIter)) {
        RealType angle;
        torsion->calcForce(angle);
        RealType currTorsionPot = torsion->getPotential();
        torsionPotential += torsion->getPotential();
        std::map<Torsion*, TorsionDataSet>::iterator i = torsionDataSets.find(torsion);
        if (i == torsionDataSets.end()) {
          TorsionDataSet dataSet;
          dataSet.prev.angle = dataSet.curr.angle = angle;
          dataSet.prev.potential = dataSet.curr.potential = currTorsionPot;
          dataSet.deltaV = 0.0;
          torsionDataSets.insert(std::map<Torsion*, TorsionDataSet>::value_type(torsion, dataSet));
        }else {
          i->second.prev.angle = i->second.curr.angle;
          i->second.prev.potential = i->second.curr.potential;
          i->second.curr.angle = angle;
          i->second.curr.potential = currTorsionPot;
          i->second.deltaV =  fabs(i->second.curr.potential -  
                                   i->second.prev.potential);
        }      
      }      

      for (inversion = mol->beginInversion(inversionIter); 
	   inversion != NULL; 
           inversion = mol->nextInversion(inversionIter)) {
        RealType angle;
        inversion->calcForce(angle);
        RealType currInversionPot = inversion->getPotential();
        inversionPotential += inversion->getPotential();
        std::map<Inversion*, InversionDataSet>::iterator i = inversionDataSets.find(inversion);
        if (i == inversionDataSets.end()) {
          InversionDataSet dataSet;
          dataSet.prev.angle = dataSet.curr.angle = angle;
          dataSet.prev.potential = dataSet.curr.potential = currInversionPot;
          dataSet.deltaV = 0.0;
          inversionDataSets.insert(std::map<Inversion*, InversionDataSet>::value_type(inversion, dataSet));
        }else {
          i->second.prev.angle = i->second.curr.angle;
          i->second.prev.potential = i->second.curr.potential;
          i->second.curr.angle = angle;
          i->second.curr.potential = currInversionPot;
          i->second.deltaV =  fabs(i->second.curr.potential -  
                                   i->second.prev.potential);
        }      
      }      
    }
    
    RealType  shortRangePotential = bondPotential + bendPotential + 
      torsionPotential +  inversionPotential;    
    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    curSnapshot->statData[Stats::SHORT_RANGE_POTENTIAL] = shortRangePotential;
    curSnapshot->statData[Stats::BOND_POTENTIAL] = bondPotential;
    curSnapshot->statData[Stats::BEND_POTENTIAL] = bendPotential;
    curSnapshot->statData[Stats::DIHEDRAL_POTENTIAL] = torsionPotential;
    curSnapshot->statData[Stats::INVERSION_POTENTIAL] = inversionPotential;
    
  }
  
  void ForceManager::calcLongRangeInteraction() {
    Snapshot* curSnapshot;
    DataStorage* config;
    RealType* frc;
    RealType* pos;
    RealType* trq;
    RealType* A;
    RealType* electroFrame;
    RealType* rc;
    RealType* particlePot;
    
    //get current snapshot from SimInfo
    curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    
    //get array pointers
    config = &(curSnapshot->atomData);
    frc = config->getArrayPointer(DataStorage::dslForce);
    pos = config->getArrayPointer(DataStorage::dslPosition);
    trq = config->getArrayPointer(DataStorage::dslTorque);
    A   = config->getArrayPointer(DataStorage::dslAmat);
    electroFrame = config->getArrayPointer(DataStorage::dslElectroFrame);
    particlePot = config->getArrayPointer(DataStorage::dslParticlePot);

    //calculate the center of mass of cutoff group
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::CutoffGroupIterator ci;
    CutoffGroup* cg;
    Vector3d com;
    std::vector<Vector3d> rcGroup;
    
    if(info_->getNCutoffGroups() > 0){
      
      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {
        for(cg = mol->beginCutoffGroup(ci); cg != NULL; 
            cg = mol->nextCutoffGroup(ci)) {
	  cg->getCOM(com);
	  rcGroup.push_back(com);
        }
      }// end for (mol)
       
      rc = rcGroup[0].getArrayPointer();
    } else {
      // center of mass of the group is the same as position of the atom  
      // if cutoff group does not exist
      rc = pos;
    }
    
    //initialize data before passing to fortran
    RealType longRangePotential[LR_POT_TYPES];
    RealType lrPot = 0.0;
    int isError = 0;

    for (int i=0; i<LR_POT_TYPES;i++){
      longRangePotential[i]=0.0; //Initialize array
    }
    
    doForceLoop(pos,
                rc,
                A,
                electroFrame,
                frc,
                trq,
	        tau.getArrayPointer(),
                longRangePotential, 
                particlePot,
                &isError );
    
    if( isError ){
      sprintf( painCave.errMsg,
	       "Error returned from the fortran force calculation.\n" );
      painCave.isFatal = 1;
      simError();
    }
    for (int i=0; i<LR_POT_TYPES;i++){
      lrPot += longRangePotential[i]; //Quick hack
    }
        
    //store the tau and long range potential    
    curSnapshot->statData[Stats::LONG_RANGE_POTENTIAL] = lrPot;
    curSnapshot->statData[Stats::VANDERWAALS_POTENTIAL] = longRangePotential[VDW_POT];
    curSnapshot->statData[Stats::ELECTROSTATIC_POTENTIAL] = longRangePotential[ELECTROSTATIC_POT];
  }

  
  void ForceManager::postCalculation() {
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::RigidBodyIterator rbIter;
    RigidBody* rb;
    Snapshot* curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();
    
    // collect the atomic forces onto rigid bodies
    
    for (mol = info_->beginMolecule(mi); mol != NULL; 
         mol = info_->nextMolecule(mi)) {
      for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
           rb = mol->nextRigidBody(rbIter)) { 
        Mat3x3d rbTau = rb->calcForcesAndTorquesAndVirial();
        tau += rbTau;
      }
    }
    
#ifdef IS_MPI
    Mat3x3d tmpTau(tau);
    MPI_Allreduce(tmpTau.getArrayPointer(), tau.getArrayPointer(), 
                  9, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif
    curSnapshot->statData.setTau(tau);
  }

} //end namespace OpenMD
