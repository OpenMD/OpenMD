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
 
#include "constraints/Rattle.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"
namespace oopse {

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

  double pabsq = pab.lengthSquare();

  double rabsq = consPair->getConsDistSquare();
  double diffsq = rabsq - pabsq;

  // the original rattle code from alan tidesley
  if (fabs(diffsq) > (consTolerance_ * rabsq * 2)){
    
    Vector3d oldPosA = consElem1->getPrevPos();
    Vector3d oldPosB = consElem2->getPrevPos();      

    Vector3d rab = oldPosA - oldPosB;    

    currentSnapshot_->wrapVector(rab);

    double rpab = dot(rab, pab);
    double rpabsq = rpab * rpab;

    if (rpabsq < (rabsq * -diffsq)){
      return consFail;
    }

    double rma = 1.0 / consElem1->getMass();
    double rmb = 1.0 / consElem2->getMass();

    double gab = diffsq / (2.0 * (rma + rmb) * rpab);

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

    double rma = 1.0 / consElem1->getMass();
    double rmb = 1.0 / consElem2->getMass();

    double rvab = dot(rab, dv);

    double gab = -rvab / ((rma + rmb) * consPair->getConsDistSquare());

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
