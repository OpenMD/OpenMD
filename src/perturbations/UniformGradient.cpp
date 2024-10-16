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

#include "perturbations/UniformGradient.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/ForceModifier.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "primitives/Molecule.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/Constants.hpp"

namespace OpenMD {
  UniformGradient::UniformGradient(SimInfo* info) :
      ForceModifier {info}, initialized {false}, doUniformGradient {false},
      doParticlePot {false} {
    simParams = info_->getSimParams();
  }

  void UniformGradient::initialize() {
    bool haveA = false;
    bool haveB = false;
    bool haveG = false;

    if (simParams->haveUniformGradientDirection1()) {
      std::vector<RealType> d1 = simParams->getUniformGradientDirection1();
      if (d1.size() != 3) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "uniformGradientDirection1: Incorrect number of parameters\n"
                 "\tspecified. There should be 3 parameters, but %zu were\n"
                 "\tspecified.\n",
                 d1.size());
        painCave.isFatal = 1;
        simError();
      }
      a_.x() = d1[0];
      a_.y() = d1[1];
      a_.z() = d1[2];

      a_.normalize();
      haveA = true;
    }

    if (simParams->haveUniformGradientDirection2()) {
      std::vector<RealType> d2 = simParams->getUniformGradientDirection2();
      if (d2.size() != 3) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "uniformGradientDirection2: Incorrect number of parameters\n"
                 "\tspecified. There should be 3 parameters, but %zu were\n"
                 "\tspecified.\n",
                 d2.size());
        painCave.isFatal = 1;
        simError();
      }
      b_.x() = d2[0];
      b_.y() = d2[1];
      b_.z() = d2[2];

      b_.normalize();
      haveB = true;
    }

    if (simParams->haveUniformGradientStrength()) {
      g_    = simParams->getUniformGradientStrength();
      haveG = true;
    }

    if (haveA && haveB && haveG) {
      doUniformGradient = true;
      cpsi_             = dot(a_, b_);

      Grad_(0, 0) = 2.0 * (a_.x() * b_.x() - cpsi_ / 3.0);
      Grad_(0, 1) = a_.x() * b_.y() + a_.y() * b_.x();
      Grad_(0, 2) = a_.x() * b_.z() + a_.z() * b_.x();
      Grad_(1, 0) = Grad_(0, 1);
      Grad_(1, 1) = 2.0 * (a_.y() * b_.y() - cpsi_ / 3.0);
      Grad_(1, 2) = a_.y() * b_.z() + a_.z() * b_.y();
      Grad_(2, 0) = Grad_(0, 2);
      Grad_(2, 1) = Grad_(1, 2);
      Grad_(2, 2) = 2.0 * (a_.z() * b_.z() - cpsi_ / 3.0);

      Grad_ *= g_ / 2.0;

    } else {
      if (!haveA) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "UniformGradient: uniformGradientDirection1 not specified.\n");
        painCave.isFatal = 1;
        simError();
      }
      if (!haveB) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "UniformGradient: uniformGradientDirection2 not specified.\n");
        painCave.isFatal = 1;
        simError();
      }
      if (!haveG) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "UniformGradient: uniformGradientStrength not specified.\n");
        painCave.isFatal = 1;
        simError();
      }
    }

    int storageLayout_ = info_->getSnapshotManager()->getAtomStorageLayout();
    if (storageLayout_ & DataStorage::dslParticlePot) doParticlePot = true;
    initialized = true;
  }

  void UniformGradient::modifyForces() {
    if (!initialized) initialize();

    SimInfo::MoleculeIterator i;
    Molecule::AtomIterator j;
    Molecule* mol;
    Atom* atom;
    AtomType* atype;
    potVec longRangePotential(0.0);

    RealType C;
    Vector3d D;
    Mat3x3d Q;

    RealType U;
    RealType fPot;
    Vector3d t;
    Vector3d f;

    Vector3d r;
    Vector3d EF;

    bool isCharge;

    if (doUniformGradient) {
      U    = 0.0;
      fPot = 0.0;

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (atom = mol->beginAtom(j); atom != NULL; atom = mol->nextAtom(j)) {
          isCharge = false;
          C        = 0.0;

          atype = atom->getAtomType();

          r = atom->getPos();
          info_->getSnapshotManager()->getCurrentSnapshot()->wrapVector(r);

          EF = Grad_ * r;

          atom->addElectricField(EF * Constants::chargeFieldConvert);

          FixedChargeAdapter fca = FixedChargeAdapter(atype);
          if (fca.isFixedCharge()) {
            isCharge = true;
            C        = fca.getCharge();
          }

          FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
          if (fqa.isFluctuatingCharge()) {
            isCharge = true;
            C += atom->getFlucQPos();
            atom->addFlucQFrc(dot(r, EF) * Constants::chargeFieldConvert);
          }

          if (isCharge) {
            f = EF * C * Constants::chargeFieldConvert;
            atom->addFrc(f);

            U = -dot(r, f);
            if (doParticlePot) { atom->addParticlePot(U); }
            fPot += U;
          }

          MultipoleAdapter ma = MultipoleAdapter(atype);
          if (ma.isDipole()) {
            D = atom->getDipole() * Constants::dipoleFieldConvert;

            f = D * Grad_;
            atom->addFrc(f);

            t = cross(D, EF);
            atom->addTrq(t);

            U = -dot(D, EF);
            if (doParticlePot) { atom->addParticlePot(U); }
            fPot += U;
          }

          if (ma.isQuadrupole()) {
            Q = atom->getQuadrupole() * Constants::dipoleFieldConvert;

            t = 2.0 * mCross(Q, Grad_);
            atom->addTrq(t);

            U = -doubleDot(Q, Grad_);
            if (doParticlePot) { atom->addParticlePot(U); }
            fPot += U;
          }
        }
      }

#ifdef IS_MPI
      MPI_Allreduce(MPI_IN_PLACE, &fPot, 1, MPI_REALTYPE, MPI_SUM,
                    MPI_COMM_WORLD);
#endif

      Snapshot* snap     = info_->getSnapshotManager()->getCurrentSnapshot();
      longRangePotential = snap->getLongRangePotentials();
      longRangePotential[ELECTROSTATIC_FAMILY] += fPot;
      snap->setLongRangePotentials(longRangePotential);
    }
  }
}  // namespace OpenMD
