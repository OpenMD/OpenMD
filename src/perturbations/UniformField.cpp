/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#include "perturbations/UniformField.hpp"

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

  UniformField::UniformField(SimInfo* info) :
      ForceModifier {info}, initialized {false}, doUniformField {false},
      doParticlePot {false} {
    simParams = info_->getSimParams();
  }

  void UniformField::initialize() {
    std::vector<RealType> ef;

    if (simParams->haveElectricField()) {
      doUniformField = true;
      ef             = simParams->getElectricField();
    }
    if (simParams->haveUniformField()) {
      doUniformField = true;
      ef             = simParams->getUniformField();
    }
    if (ef.size() != 3) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "UniformField: Incorrect number of parameters specified.\n"
               "\tthere should be 3 parameters, but %lu were specified.\n",
               ef.size());
      painCave.isFatal = 1;
      simError();
    }
    EF.x() = ef[0];
    EF.y() = ef[1];
    EF.z() = ef[2];

    int storageLayout_ = info_->getSnapshotManager()->getAtomStorageLayout();
    if (storageLayout_ & DataStorage::dslParticlePot) doParticlePot = true;
    initialized = true;
  }

  void UniformField::modifyForces() {
    if (!initialized) initialize();

    SimInfo::MoleculeIterator i;
    Molecule::AtomIterator j;
    Molecule* mol;
    Atom* atom;
    AtomType* atype;
    potVec longRangePotential(0.0);

    RealType C;
    Vector3d D;
    RealType U;
    RealType fPot;
    Vector3d t;
    Vector3d f;
    Vector3d r;

    bool isCharge;

    if (doUniformField) {
      U    = 0.0;
      fPot = 0.0;

      for (mol = info_->beginMolecule(i); mol != NULL;
           mol = info_->nextMolecule(i)) {
        for (atom = mol->beginAtom(j); atom != NULL; atom = mol->nextAtom(j)) {
          isCharge = false;
          C        = 0.0;

          atype = atom->getAtomType();

          // ad-hoc choice of the origin for potential calculation and
          // fluctuating charge force:

          r = atom->getPos();
          info_->getSnapshotManager()->getCurrentSnapshot()->wrapVector(r);

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

            t = cross(D, EF);
            atom->addTrq(t);

            U = -dot(D, EF);

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
