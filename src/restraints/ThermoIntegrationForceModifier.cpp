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

#include "restraints/ThermoIntegrationForceModifier.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/SimInfo.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD {

  ThermoIntegrationForceModifier::ThermoIntegrationForceModifier(
      SimInfo* info) :
      RestraintForceModifier {info} {
    currSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    simParam_     = info_->getSimParams();

    RealType tIntLambda, tIntK;

    if (simParam_->haveThermodynamicIntegrationLambda()) {
      tIntLambda = simParam_->getThermodynamicIntegrationLambda();
    } else {
      tIntLambda = 1.0;
      sprintf(painCave.errMsg,
              "ThermoIntegration error: the transformation parameter\n"
              "\t(lambda) was not specified. OpenMD will use a default\n"
              "\tvalue of %f. To set lambda, use the \n"
              "\tthermodynamicIntegrationLambda variable.\n",
              tIntLambda);
      painCave.isFatal = 0;
      simError();
    }

    if (simParam_->haveThermodynamicIntegrationK()) {
      tIntK = simParam_->getThermodynamicIntegrationK();
    } else {
      tIntK = 1.0;
      sprintf(painCave.errMsg,
              "ThermoIntegration Warning: the tranformation parameter\n"
              "\texponent (k) was not specified. OpenMD will use a default\n"
              "\tvalue of %f. To set k, use the thermodynamicIntegrationK\n"
              "\tvariable.\n",
              tIntK);
      painCave.isFatal = 0;
      simError();
    }

    // build the scaling factor used to modulate the forces and torques
    factor_ = pow(tIntLambda, tIntK);
  }

  void ThermoIntegrationForceModifier::modifyForces() {
    Snapshot* curSnapshot;
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::IntegrableObjectIterator ii;
    StuntDouble* sd;
    Vector3d frc;
    Vector3d trq;
    Mat3x3d tempTau;

    // now scale forces and torques of all the sds
    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ii); sd != NULL;
           sd = mol->nextIntegrableObject(ii)) {
        frc = sd->getFrc();
        frc *= factor_;
        sd->setFrc(frc);

        if (sd->isDirectional()) {
          trq = sd->getTrq();
          trq *= factor_;
          sd->setTrq(trq);
        }
      }
    }

    curSnapshot = info_->getSnapshotManager()->getCurrentSnapshot();

    // set rawPotential to be the unscaled potential energy
    RealType rawPot = curSnapshot->getPotentialEnergy();
    curSnapshot->setRawPotential(rawPot);

    RealType rp = curSnapshot->getRestraintPotential();

    // scale the potential and update the snapshot
    rawPot *= factor_;
    curSnapshot->setPotentialEnergy(rawPot);

    // scale the virial tensor
    tempTau = curSnapshot->getVirialTensor();
    tempTau *= factor_;
    curSnapshot->setVirialTensor(tempTau);

    // now, on to the applied restraining potentials (if needed):
    RealType scaledRestPot(0.0);
    RealType restPot(0.0);

    if (simParam_->getUseRestraints()) {
      // do restraints from RestraintForceModifier:
      scaledRestPot = doRestraints(1.0 - factor_);
      restPot       = getUnscaledPotential();
    }

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &scaledRestPot, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &restPot, 1, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    // give the final values to the snapshot
    curSnapshot->setPotentialEnergy(rawPot + scaledRestPot);
    curSnapshot->setRestraintPotential(rp + restPot);
  }
}  // namespace OpenMD
