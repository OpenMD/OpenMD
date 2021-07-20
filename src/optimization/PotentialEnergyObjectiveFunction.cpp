/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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

#include "optimization/PotentialEnergyObjectiveFunction.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

namespace OpenMD {

  PotentialEnergyObjectiveFunction::PotentialEnergyObjectiveFunction(
      SimInfo* info, ForceManager* forceMan) :
      info_(info),
      forceMan_(forceMan), thermo(info), hasFlucQ_(false) {
    shake_ = new Shake(info_);

    if (info_->usesFluctuatingCharges()) {
      if (info_->getNFluctuatingCharges() > 0) {
        hasFlucQ_      = true;
        fqConstraints_ = new FluctuatingChargeConstraints(info_);
        bool cr        = info_->getSimParams()
                      ->getFluctuatingChargeParameters()
                      ->getConstrainRegions();
        fqConstraints_->setConstrainRegions(cr);
      }
    }
  }

  RealType PotentialEnergyObjectiveFunction::value(
      const DynamicVector<RealType>& x) {
    setCoor(x);
    shake_->constraintR();
    forceMan_->calcForces();
    if (hasFlucQ_) fqConstraints_->applyConstraints();
    shake_->constraintF();
    return thermo.getPotential();
  }

  void PotentialEnergyObjectiveFunction::gradient(
      DynamicVector<RealType>& grad, const DynamicVector<RealType>& x) {
    setCoor(x);
    shake_->constraintR();
    forceMan_->calcForces();
    if (hasFlucQ_) fqConstraints_->applyConstraints();
    shake_->constraintF();
    getGrad(grad);
  }

  RealType PotentialEnergyObjectiveFunction::valueAndGradient(
      DynamicVector<RealType>& grad, const DynamicVector<RealType>& x) {
    setCoor(x);
    shake_->constraintR();
    forceMan_->calcForces();
    if (hasFlucQ_) fqConstraints_->applyConstraints();
    shake_->constraintF();
    getGrad(grad);
    return thermo.getPotential();
  }

  void PotentialEnergyObjectiveFunction::setCoor(
      const DynamicVector<RealType>& x) const {
    Vector3d position;
    Vector3d eulerAngle;
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule::AtomIterator ai;
    Molecule* mol;
    StuntDouble* sd;
    Atom* atom;

    info_->getSnapshotManager()->advance();

    int index;
#ifdef IS_MPI
    index = displacements_[myrank_];
#else
    index = 0;
#endif

    info_->getSnapshotManager()->advance();

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        position[0] = x[index++];
        position[1] = x[index++];
        position[2] = x[index++];

        sd->setPos(position);

        if (sd->isDirectional()) {
          eulerAngle[0] = x[index++];
          eulerAngle[1] = x[index++];
          eulerAngle[2] = x[index++];

          sd->setEuler(eulerAngle);

          if (sd->isRigidBody()) {
            RigidBody* rb = static_cast<RigidBody*>(sd);
            rb->updateAtoms();
          }
        }
      }
    }

    if (hasFlucQ_) {
      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (atom = mol->beginFluctuatingCharge(ai); atom != NULL;
             atom = mol->nextFluctuatingCharge(ai)) {
          atom->setFlucQPos(x[index++]);
        }
      }
    }
  }

  void PotentialEnergyObjectiveFunction::getGrad(
      DynamicVector<RealType>& grad) {
    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule::AtomIterator ai;
    Molecule* mol;
    StuntDouble* sd;
    Atom* atom;
    std::vector<RealType> myGrad;

    int index;
#ifdef IS_MPI
    index = displacements_[myrank_];
    grad.setZero();
#else
    index = 0;
#endif

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        myGrad = sd->getGrad();

        for (size_t k = 0; k < myGrad.size(); ++k) {
          grad[index++] = myGrad[k];
        }
      }
    }

    if (hasFlucQ_) {
      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (atom = mol->beginFluctuatingCharge(ai); atom != NULL;
             atom = mol->nextFluctuatingCharge(ai)) {
          grad[index++] = -atom->getFlucQFrc();
        }
      }
    }
#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &grad[0], ndf_, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif
  }

  DynamicVector<RealType> PotentialEnergyObjectiveFunction::setInitialCoords() {
#ifdef IS_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nproc_);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
    std::vector<int> onProc(nproc_, 0);

    displacements_.clear();
    displacements_.resize(nproc_, 0);
#endif

    SimInfo::MoleculeIterator i;
    Molecule::IntegrableObjectIterator j;
    Molecule::AtomIterator ai;
    Molecule* mol;
    StuntDouble* sd;
    Atom* atom;

    Vector3d pos;
    Vector3d eulerAngle;

    ndf_ = 0;

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        ndf_ += 3;

        if (sd->isDirectional()) { ndf_ += 3; }
      }
    }

    if (hasFlucQ_) {
      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (atom = mol->beginFluctuatingCharge(ai); atom != NULL;
             atom = mol->nextFluctuatingCharge(ai)) {
          ndf_++;
        }
      }
    }

#ifdef IS_MPI
    MPI_Allgather(&ndf_, 1, MPI_INT, &onProc[0], 1, MPI_INT, MPI_COMM_WORLD);

    ndf_ = 0;
    for (int iproc = 0; iproc < nproc_; iproc++) {
      ndf_ += onProc[iproc];
    }

    displacements_[0] = 0;
    for (int iproc = 1; iproc < nproc_; iproc++) {
      displacements_[iproc] = displacements_[iproc - 1] + onProc[iproc - 1];
    }
#endif

    DynamicVector<RealType> xinit(ndf_, 0.0);

    int index;
#ifdef IS_MPI
    index = displacements_[myrank_];
#else
    index = 0;
#endif

    for (mol = info_->beginMolecule(i); mol != NULL;
         mol = info_->nextMolecule(i)) {
      for (sd = mol->beginIntegrableObject(j); sd != NULL;
           sd = mol->nextIntegrableObject(j)) {
        pos            = sd->getPos();
        xinit[index++] = pos[0];
        xinit[index++] = pos[1];
        xinit[index++] = pos[2];

        if (sd->isDirectional()) {
          eulerAngle     = sd->getEuler();
          xinit[index++] = eulerAngle[0];
          xinit[index++] = eulerAngle[1];
          xinit[index++] = eulerAngle[2];
        }
      }
    }

    if (hasFlucQ_) {
      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (atom = mol->beginFluctuatingCharge(ai); atom != NULL;
             atom = mol->nextFluctuatingCharge(ai)) {
          xinit[index++] = atom->getFlucQPos();
        }
      }
    }

    return xinit;
  }
}  // namespace OpenMD
