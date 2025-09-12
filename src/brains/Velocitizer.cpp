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

#include "brains/Velocitizer.hpp"

#include <memory>
#include <random>

#include "brains/Thermo.hpp"
#include "flucq/FluctuatingChargeConstraints.hpp"
#include "math/SquareMatrix3.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/Constants.hpp"
#include "utils/RandNumGen.hpp"

namespace OpenMD {

  Velocitizer::Velocitizer(SimInfo* info) :
      info_(info), thermo_(info), globals_(info->getSimParams()),
      randNumGen_(info->getRandomNumberGenerator()) {}

  void Velocitizer::scale(RealType lambda) {
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule* mol;
    StuntDouble* sd;
    Vector3d v, j;

    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ioi); sd != NULL;
           sd = mol->nextIntegrableObject(ioi)) {
        v = sd->getVel();
        v *= lambda;
        sd->setVel(v);

        if (sd->isDirectional()) {
          j = sd->getJ();
          j *= lambda;
          sd->setJ(j);
        }
      }
    }

    // We've modified velocities, so clear any derived properties in
    // the Snapshot:
    info_->getSnapshotManager()->getCurrentSnapshot()->clearDerivedProperties();

    removeComDrift();

    // Remove angular drift if we are not using periodic boundary
    // conditions:
    if (!globals_->getUsePeriodicBoundaryConditions()) removeAngularDrift();
  }

  void Velocitizer::randomize(RealType temperature) {
    Vector3d v;
    Vector3d j;
    Mat3x3d I;
    int l, m, n;
    Vector3d vdrift;
    RealType vbar;
    RealType jbar;
    RealType av2;
    RealType kebar;

    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule* mol;
    StuntDouble* sd;

    std::normal_distribution<RealType> normalDistribution {0.0, 1.0};

    kebar = Constants::kB * temperature * info_->getNdfRaw() /
            (2.0 * info_->getNdf());
    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ioi); sd != NULL;
           sd = mol->nextIntegrableObject(ioi)) {
        // uses equipartition theory to solve for vbar in angstrom/fs

        av2  = 2.0 * kebar / sd->getMass();
        vbar = sqrt(av2);

        // picks random velocities from a gaussian distribution
        // centered on vbar

        for (int k = 0; k < 3; k++) {
          v[k] = vbar * normalDistribution(*randNumGen_);
        }
        sd->setVel(v);

        if (sd->isDirectional()) {
          I = sd->getI();

          if (sd->isLinear()) {
            l = sd->linearAxis();
            m = (l + 1) % 3;
            n = (l + 2) % 3;

            j[l] = 0.0;
            jbar = sqrt(2.0 * kebar * I(m, m));
            j[m] = jbar * normalDistribution(*randNumGen_);
            jbar = sqrt(2.0 * kebar * I(n, n));
            j[n] = jbar * normalDistribution(*randNumGen_);
          } else {
            for (int k = 0; k < 3; k++) {
              jbar = sqrt(2.0 * kebar * I(k, k));
              j[k] = jbar * normalDistribution(*randNumGen_);
            }
          }

          sd->setJ(j);
        }
      }
    }

    // We've modified velocities, so clear any derived properties in
    // the Snapshot:
    info_->getSnapshotManager()->getCurrentSnapshot()->clearDerivedProperties();
    
    removeComDrift();

    // Remove angular drift if we are not using periodic boundary
    // conditions:
    if (!globals_->getUsePeriodicBoundaryConditions()) removeAngularDrift();

  }

  void Velocitizer::randomizeChargeVelocity(RealType temperature) {
    RealType aw2;
    RealType kebar;
    RealType wbar;

    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule* mol;
    StuntDouble* sd;
    FluctuatingChargeParameters* fqParams;
    FluctuatingChargeConstraints* fqConstraints;

    std::normal_distribution<RealType> normalDistribution {0.0, 1.0};

    Globals* simParams = info_->getSimParams();
    fqParams           = simParams->getFluctuatingChargeParameters();

    fqConstraints = new FluctuatingChargeConstraints(info_);
    fqConstraints->setConstrainRegions(fqParams->getConstrainRegions());

    int nConstrain =
        fqConstraints
            ->getNumberOfFlucQConstraints();  // no of constraints in charge
    int dfRaw = fqConstraints->getNumberOfFlucQAtoms();  // no of FlucQ freedom
    int dfActual = dfRaw - nConstrain;
    kebar        = dfRaw * Constants::kb * temperature / (2 * dfActual);

    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ioi); sd != NULL;
           sd = mol->nextIntegrableObject(ioi)) {
        if (sd->isAtom()) {
          Atom* atom                   = static_cast<Atom*>(sd);
          AtomType* atomType           = atom->getAtomType();
          FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
          if (fqa.isFluctuatingCharge()) {
            // uses equipartition theory to solve for vbar in angstrom/fs

            aw2  = 2.0 * kebar / atom->getChargeMass();
            wbar = sqrt(aw2);

            // picks random velocities from a gaussian distribution
            // centered on vbar
            atom->setFlucQVel(wbar * normalDistribution(*randNumGen_));
          }
        }

        // randomization of the charge velocities for atoms in the rigidbody

        if (sd->isRigidBody()) {
          RigidBody* rigidbody = static_cast<RigidBody*>(sd);
          vector<Atom*> atomList;
          atomList = rigidbody->getAtoms();

          for (size_t i = 0; i < atomList.size(); ++i) {
            Atom* atom                   = atomList[i];
            AtomType* atomType           = atom->getAtomType();
            FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
            if (fqa.isFluctuatingCharge()) {
              // uses equipartition theory to solve for vbar in angstrom/fs
              aw2  = 2.0 * kebar / atom->getChargeMass();
              wbar = sqrt(aw2);
              // picks random velocities from a gaussian distribution
              // centered on vbar
              atom->setFlucQVel(wbar * normalDistribution(*randNumGen_));
            }
          }
        }
      }
    }
    fqConstraints->applyConstraintsOnChargeVelocities();
  }

  void Velocitizer::removeComDrift() {
    // Get the Center of Mass drift velocity.
    Vector3d vdrift = thermo_.getComVel();

    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule* mol;
    StuntDouble* sd;

    //  Corrects for the center of mass drift.
    // sums all the momentum and divides by total mass.
    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ioi); sd != NULL;
           sd = mol->nextIntegrableObject(ioi)) {
        sd->setVel(sd->getVel() - vdrift);
      }
    }

    // We've modified velocities, so clear any derived properties in
    // the Snapshot:
    info_->getSnapshotManager()->getCurrentSnapshot()->clearDerivedProperties();
  }

  void Velocitizer::removeAngularDrift() {
    // Get the Center of Mass drift velocity.
    Vector3d vdrift;
    Vector3d com;

    thermo_.getComAll(com, vdrift);

    Mat3x3d inertiaTensor;
    Vector3d angularMomentum;
    Vector3d omega;

    thermo_.getInertiaTensor(inertiaTensor, angularMomentum);

    // We now need the inverse of the inertia tensor.
    inertiaTensor = inertiaTensor.inverse();
    omega         = inertiaTensor * angularMomentum;

    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Molecule* mol;
    StuntDouble* sd;
    Vector3d tempComPos;

    // Corrects for the center of mass angular drift by summing all
    // the angular momentum and dividing by the total mass.

    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ioi); sd != NULL;
           sd = mol->nextIntegrableObject(ioi)) {
        tempComPos = sd->getPos() - com;
        sd->setVel((sd->getVel() - vdrift) - cross(omega, tempComPos));
      }
    }

    // We've modified velocities, so clear any derived properties in
    // the Snapshot:
    info_->getSnapshotManager()->getCurrentSnapshot()->clearDerivedProperties();
  }
}  // namespace OpenMD
