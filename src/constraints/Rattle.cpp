/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
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
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "constraints/Rattle.hpp"

#include <cmath>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "primitives/Molecule.hpp"
#include "utils/simError.h"

namespace OpenMD {

  Rattle::Rattle(SimInfo* info) :
      info_(info), maxConsIteration_(10), consTolerance_(1.0e-6),
      doRattle_(false), currConstraintTime_(0.0) {
    if (info_->getNGlobalConstraints() > 0) doRattle_ = true;

    if (!doRattle_) return;

    Globals* simParams = info_->getSimParams();

    if (simParams->haveDt()) {
      dt_ = simParams->getDt();
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Rattle Error: dt is not set\n");
      painCave.isFatal = 1;
      simError();
    }

    currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    if (simParams->haveConstraintTime()) {
      constraintTime_ = simParams->getConstraintTime();
    } else {
      constraintTime_ = simParams->getStatusTime();
    }

    constraintOutputFile_ =
        getPrefix(info_->getFinalConfigFileName()) + ".constraintForces";

    // create ConstraintWriter
    constraintWriter_ =
        new ConstraintWriter(info_, constraintOutputFile_.c_str());

    if (!constraintWriter_) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Failed to create ConstraintWriter\n");
      painCave.isFatal = 1;
      simError();
    }
  }

  Rattle::~Rattle() { delete constraintWriter_; }

  void Rattle::constraintA() {
    if (!doRattle_) return;
    doConstraint(&Rattle::constraintPairA);
  }
  void Rattle::constraintB() {
    if (!doRattle_) return;
    doConstraint(&Rattle::constraintPairB);

    if (currentSnapshot_->getTime() >= currConstraintTime_) {
      Molecule* mol;
      SimInfo::MoleculeIterator mi;
      ConstraintPair* consPair;
      Molecule::ConstraintPairIterator cpi;
      std::list<ConstraintPair*> constraints;
      for (mol = info_->beginMolecule(mi); mol != NULL;
           mol = info_->nextMolecule(mi)) {
        for (consPair = mol->beginConstraintPair(cpi); consPair != NULL;
             consPair = mol->nextConstraintPair(cpi)) {
          constraints.push_back(consPair);
        }
      }
      constraintWriter_->writeConstraintForces(constraints);
      currConstraintTime_ += constraintTime_;
    }
  }

  void Rattle::doConstraint(ConstraintPairFuncPtr func) {
    if (!doRattle_) return;

    Molecule* mol;
    SimInfo::MoleculeIterator mi;
    ConstraintElem* consElem;
    Molecule::ConstraintElemIterator cei;
    ConstraintPair* consPair;
    Molecule::ConstraintPairIterator cpi;

    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (consElem = mol->beginConstraintElem(cei); consElem != NULL;
           consElem = mol->nextConstraintElem(cei)) {
        consElem->setMoved(true);
        consElem->setMoving(false);
      }
      for (consPair = mol->beginConstraintPair(cpi); consPair != NULL;
           consPair = mol->nextConstraintPair(cpi)) {
        consPair->resetConstraintForce();
      }
    }

    // main loop of constraint algorithm
    int done      = 0;
    int iteration = 0;
    while (!done && iteration < maxConsIteration_) {
      done = 1;

      // loop over every constraint pair

      for (mol = info_->beginMolecule(mi); mol != NULL;
           mol = info_->nextMolecule(mi)) {
        for (consPair = mol->beginConstraintPair(cpi); consPair != NULL;
             consPair = mol->nextConstraintPair(cpi)) {
          // dispatch constraint algorithm
          if (consPair->isMoved()) {
            int exeStatus = (this->*func)(consPair);

            switch (exeStatus) {
            case consFail:
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "Constraint failure in Rattle::constrainA, "
                       "Constraint Fail\n");
              painCave.isFatal = 1;
              simError();

              break;
            case consSuccess:
              // constrain the pair by moving two elements
              done = 0;
              consPair->getConsElem1()->setMoving(true);
              consPair->getConsElem2()->setMoving(true);
              break;
            case consAlready:
              // current pair is already constrained, do not need to
              // move the elements
              break;
            default:
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "ConstraintAlgorithm::doConstraint() "
                       "Error: unrecognized status");
              painCave.isFatal = 1;
              simError();
              break;
            }
          }
        }
      }  // end for(iter->first())

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
#endif

      errorCheckPoint();

      for (mol = info_->beginMolecule(mi); mol != NULL;
           mol = info_->nextMolecule(mi)) {
        for (consElem = mol->beginConstraintElem(cei); consElem != NULL;
             consElem = mol->nextConstraintElem(cei)) {
          consElem->setMoved(consElem->getMoving());
          consElem->setMoving(false);
        }
      }
      iteration++;
    }  // end while

    if (!done) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Constraint failure in Rattle::constrainA, "
               "too many iterations: %d\n",
               iteration);
      painCave.isFatal = 1;
      simError();
    }

    errorCheckPoint();
  }

  int Rattle::constraintPairA(ConstraintPair* consPair) {
    ConstraintElem* consElem1 = consPair->getConsElem1();
    ConstraintElem* consElem2 = consPair->getConsElem2();

    Vector3d posA = consElem1->getPos();
    Vector3d posB = consElem2->getPos();

    Vector3d pab = posA - posB;

    // periodic boundary condition

    currentSnapshot_->wrapVector(pab);

    RealType pabsq = pab.lengthSquare();

    RealType rabsq  = consPair->getConsDistSquare();
    RealType diffsq = rabsq - pabsq;

    // the original rattle code from alan tidesley
    if (fabs(diffsq) > (consTolerance_ * rabsq * 2.0)) {
      Vector3d oldPosA = consElem1->getPrevPos();
      Vector3d oldPosB = consElem2->getPrevPos();

      Vector3d rab = oldPosA - oldPosB;

      currentSnapshot_->wrapVector(rab);

      RealType rpab   = dot(rab, pab);
      RealType rpabsq = rpab * rpab;

      if (rpabsq < (rabsq * -diffsq)) { return consFail; }

      RealType rma = 1.0 / consElem1->getMass();
      RealType rmb = 1.0 / consElem2->getMass();

      RealType gab = diffsq / (2.0 * (rma + rmb) * rpab);

      Vector3d delta = rab * gab;

      // set atom1's position
      posA += rma * delta;
      consElem1->setPos(posA);

      // set atom2's position
      posB -= rmb * delta;
      consElem2->setPos(posB);

      delta /= dt_;

      // set atom1's velocity
      Vector3d velA = consElem1->getVel();
      velA += rma * delta;
      consElem1->setVel(velA);

      // set atom2's velocity
      Vector3d velB = consElem2->getVel();
      velB -= rmb * delta;
      consElem2->setVel(velB);

      // report the constraint force back to the constraint pair:
      Vector3d fcons = 2.0 * delta / dt_;
      RealType proj  = copysign(fcons.length(), dot(fcons, rab));

      consPair->addConstraintForce(proj);
      return consSuccess;
    } else {
      return consAlready;
    }
  }

  int Rattle::constraintPairB(ConstraintPair* consPair) {
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

    if (fabs(gab) > consTolerance_) {
      Vector3d delta = rab * gab;
      velA += rma * delta;
      consElem1->setVel(velA);

      velB -= rmb * delta;
      consElem2->setVel(velB);

      // report the constraint force back to the constraint pair:

      Vector3d fcons = 2.0 * delta / dt_;
      RealType proj  = copysign(fcons.length(), dot(fcons, rab));

      consPair->addConstraintForce(proj);
      return consSuccess;
    } else {
      return consAlready;
    }
  }

}  // namespace OpenMD
