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
 
#include "constraints/Rattle.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
namespace OpenMD {

  Rattle::Rattle(SimInfo* info) : info_(info), maxConsIteration_(10), consTolerance_(1.0e-6) {
    
    if (info_->getSimParams()->haveDt()) {
      dt_ = info_->getSimParams()->getDt();
    } else {
      sprintf(painCave.errMsg,
	      "Integrator Error: dt is not set\n");
      painCave.isFatal = 1;
      simError();
    }    
    
    currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
  }

  void Rattle::constraintA() {
    if (info_->getNConstraints() > 0) {
      doConstraint(&Rattle::constraintPairA);
    }
  }
  void Rattle::constraintB() {
    if (info_->getNConstraints() > 0) {
      doConstraint(&Rattle::constraintPairB);
    }
  }

  void Rattle::doConstraint(ConstraintPairFuncPtr func) {
    Molecule* mol;
    SimInfo::MoleculeIterator mi;
    ConstraintElem* consElem;
    Molecule::ConstraintElemIterator cei;
    ConstraintPair* consPair;
    Molecule::ConstraintPairIterator cpi;
    
    for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
      for (consElem = mol->beginConstraintElem(cei); consElem != NULL; consElem = mol->nextConstraintElem(cei)) {
	consElem->setMoved(true);
	consElem->setMoving(false);
      }
    }
    
    //main loop of constraint algorithm
    bool done = false;
    int iteration = 0;
    while(!done && iteration < maxConsIteration_){
      done = true;

      //loop over every constraint pair

      for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
	for (consPair = mol->beginConstraintPair(cpi); consPair != NULL; consPair = mol->nextConstraintPair(cpi)) {


	  //dispatch constraint algorithm
	  if(consPair->isMoved()) {
	    int exeStatus = (this->*func)(consPair);

	    switch(exeStatus){
	    case consFail:
	      sprintf(painCave.errMsg,
		      "Constraint failure in Rattle::constrainA, Constraint Fail\n");
	      painCave.isFatal = 1;
	      simError();                             
                            
	      break;
	    case consSuccess:
	      //constrain the pair by moving two elements
	      done = false;
	      consPair->getConsElem1()->setMoving(true);
	      consPair->getConsElem2()->setMoving(true);
	      break;
	    case consAlready:
	      //current pair is already constrained, do not need to move the elements
	      break;
	    default:          
	      sprintf(painCave.errMsg, "ConstraintAlgorithm::doConstrain() Error: unrecognized status");
	      painCave.isFatal = 1;
	      simError();                           
	      break;
	    }      
	  }
	}
      }//end for(iter->first())


      for (mol = info_->beginMolecule(mi); mol != NULL; mol = info_->nextMolecule(mi)) {
	for (consElem = mol->beginConstraintElem(cei); consElem != NULL; consElem = mol->nextConstraintElem(cei)) {
	  consElem->setMoved(consElem->getMoving());
	  consElem->setMoving(false);
	}
      }

      iteration++;
    }//end while

    if (!done){
      sprintf(painCave.errMsg,
              "Constraint failure in Rattle::constrainA, too many iterations: %d\n",
              iteration);
      painCave.isFatal = 1;
      simError();    
    }
  }

  int Rattle::constraintPairA(ConstraintPair* consPair){
    ConstraintElem* consElem1 = consPair->getConsElem1();
    ConstraintElem* consElem2 = consPair->getConsElem2();

    Vector3d posA = consElem1->getPos();
    Vector3d posB = consElem2->getPos();

    Vector3d pab = posA -posB;   

    //periodic boundary condition

    currentSnapshot_->wrapVector(pab);

    RealType pabsq = pab.lengthSquare();

    RealType rabsq = consPair->getConsDistSquare();
    RealType diffsq = rabsq - pabsq;

    // the original rattle code from alan tidesley
    if (fabs(diffsq) > (consTolerance_ * rabsq * 2)){
    
      Vector3d oldPosA = consElem1->getPrevPos();
      Vector3d oldPosB = consElem2->getPrevPos();      

      Vector3d rab = oldPosA - oldPosB;    

      currentSnapshot_->wrapVector(rab);

      RealType rpab = dot(rab, pab);
      RealType rpabsq = rpab * rpab;

      if (rpabsq < (rabsq * -diffsq)){
	return consFail;
      }

      RealType rma = 1.0 / consElem1->getMass();
      RealType rmb = 1.0 / consElem2->getMass();

      RealType gab = diffsq / (2.0 * (rma + rmb) * rpab);

      Vector3d delta = rab * gab;

      //set atom1's position
      posA += rma * delta;    
      consElem1->setPos(posA);

      //set atom2's position
      posB -= rmb * delta;
      consElem2->setPos(posB);

      delta /= dt_;
    
      //set atom1's velocity
      Vector3d velA = consElem1->getVel();
      velA += rma * delta;
      consElem1->setVel(velA);

      //set atom2's velocity
      Vector3d velB = consElem2->getVel();
      velB -= rmb * delta;
      consElem2->setVel(velB);

      return consSuccess;
    }
    else
      return consAlready;
  
  }


  int Rattle::constraintPairB(ConstraintPair* consPair){
    ConstraintElem* consElem1 = consPair->getConsElem1();
    ConstraintElem* consElem2 = consPair->getConsElem2();

  
    Vector3d velA = consElem1->getVel();
    Vector3d velB = consElem2->getVel();

    Vector3d dv = velA - velB;

    Vector3d posA = consElem1->getPos();
    Vector3d posB = consElem2->getPos();

    Vector3d rab = posA - posB;

    currentSnapshot_->wrapVector(rab);

    RealType rma = 1.0 / consElem1->getMass();
    RealType rmb = 1.0 / consElem2->getMass();

    RealType rvab = dot(rab, dv);

    RealType gab = -rvab / ((rma + rmb) * consPair->getConsDistSquare());

    if (fabs(gab) > consTolerance_){
      Vector3d delta = rab * gab;
      
      velA += rma * delta;
      consElem1->setVel(velA);
      
      velB -= rmb * delta;
      consElem2->setVel(velB);

      return consSuccess;
    }
    else
      return consAlready;

  }

}
