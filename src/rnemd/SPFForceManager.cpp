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

#include "rnemd/SPFForceManager.hpp"

#include <config.h>

#include <vector>

#include "brains/ForceManager.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Snapshot.hpp"
#include "brains/Thermo.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"

namespace OpenMD::RNEMD {

  SPFForceManager::SPFForceManager(SimInfo* info) :
      ForceManager {info}, lambda_ {}, potentialSource_ {}, potentialSink_ {} {
    thermo_          = std::make_unique<Thermo>(info);
    currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

    k_ = info_->getSimParams()->getRNEMDParameters()->getSPFScalingPower();

    int nAtoms        = info_->getNAtoms();
    int nRigidBodies  = info_->getNRigidBodies();
    int nCutoffGroups = info_->getNCutoffGroups();
    int storageLayout = info_->getSnapshotManager()->getStorageLayout();

    bool usePBC = info_->getSimParams()->getUsePeriodicBoundaryConditions();

    temporarySourceSnapshot_ = new Snapshot(nAtoms, nRigidBodies, nCutoffGroups,
                                            storageLayout, usePBC);

    temporarySinkSnapshot_ = new Snapshot(nAtoms, nRigidBodies, nCutoffGroups,
                                          storageLayout, usePBC);
  }

  SPFForceManager::~SPFForceManager() {
    delete temporarySourceSnapshot_;
    delete temporarySinkSnapshot_;
  }

  void SPFForceManager::calcForces() {
    // Current snapshot with selected molecule in source slab
    ForceManager::calcForces();
    potentialSource_ = currentSnapshot_->getPotentialEnergy();

    if (hasSelectedMolecule_) {
      Vector3d prevSourceCom {}, currentSourceCom {}, delta {};

      // Only the processor with the selected molecule should do this step:
      if (selectedMolecule_) {
        prevSourceCom    = selectedMolecule_->getPrevCom();
        currentSourceCom = selectedMolecule_->getCom();

        delta = currentSourceCom - prevSourceCom;
      }

      if (temporarySourceSnapshot_ && currentSnapshot_) {
        *temporarySourceSnapshot_ = *currentSnapshot_;
        currentSnapshot_->clearDerivedProperties();
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Either temporarySourceSnapshot or currentSnapshot "
                 "has a null value.\n");
        painCave.isFatal  = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      }

      // Only the processor with the selected molecule should do this step:
      if (selectedMolecule_) {
        currentSinkCom_ += delta;
        selectedMolecule_->setCom(currentSinkCom_);
      }

      // Current snapshot with selected molecule in sink slab
      ForceManager::calcForces();
      potentialSink_ = currentSnapshot_->getPotentialEnergy();

      if (temporarySinkSnapshot_ && currentSnapshot_) {
        *temporarySinkSnapshot_ = *currentSnapshot_;
        currentSnapshot_->clearDerivedProperties();
      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "Either temporarySinkSnapshot or currentSnapshot "
                 "has a null value.\n");
        painCave.isFatal  = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      }

      combineForcesAndTorques();
      updatePotentials();
      updateVirialTensor();

      // Only the processor with the selected molecule should do this step:
      if (selectedMolecule_) { selectedMolecule_->setCom(currentSourceCom); }
    } else {
      *temporarySourceSnapshot_ = *currentSnapshot_;
      *temporarySinkSnapshot_   = *currentSnapshot_;
      potentialSink_            = currentSnapshot_->getPotentialEnergy();
    }
  }

  void SPFForceManager::setSelectedMolecule(Molecule* selectedMolecule,
                                            Vector3d newCom) {
    if (selectedMolecule) {
      selectedMolecule_ = selectedMolecule;
      currentSinkCom_   = newCom;
    } else {
      selectedMolecule_ = nullptr;
    }
  }

  bool SPFForceManager::updateLambda(RealType& particleTarget,
                                     RealType& deltaLambda) {
    bool updateSelectedMolecule {false};

    if (hasSelectedMolecule_) {
      lambda_ += std::fabs(particleTarget);

      if (f_lambda(lambda_ + std::fabs(particleTarget)) > 1.0 &&
          f_lambda(lambda_) < 1.0) {
        deltaLambda = particleTarget -
                      (f_lambda(lambda_ + std::fabs(particleTarget)) - 1.0);
      } else {
        deltaLambda = particleTarget;
      }

      currentSnapshot_->clearDerivedProperties();

      combineForcesAndTorques();
      updatePotentials();
      updateVirialTensor();

      if (f_lambda(lambda_) > 1.0 ||
          std::fabs(f_lambda(lambda_) - 1.0) < 1e-6) {
        lambda_ = 0.0;

        // Only the processor with the selected molecule should do this step:
        if (selectedMolecule_) {
          selectedMolecule_->setCom(currentSinkCom_);
          selectedMolecule_ = nullptr;
        }

        hasSelectedMolecule_   = false;
        updateSelectedMolecule = true;
      }
    }

    return updateSelectedMolecule;
  }

  void SPFForceManager::combineForcesAndTorques() {
    // Calculate lambda-averaged forces on all atoms and potentials:
    SimInfo::MoleculeIterator mi;
    Molecule* mol;
    Molecule::IntegrableObjectIterator ii;
    StuntDouble* sd;

    RealType result = f_lambda(lambda_);

    // Now scale forces and torques of all the sds
    for (mol = info_->beginMolecule(mi); mol != NULL;
         mol = info_->nextMolecule(mi)) {
      for (sd = mol->beginIntegrableObject(ii); sd != NULL;
           sd = mol->nextIntegrableObject(ii)) {
        sd->combineForcesAndTorques(temporarySourceSnapshot_,
                                    temporarySinkSnapshot_, 1.0 - result,
                                    result);
      }
    }
  }

  void SPFForceManager::updatePotentials() {
    updateLongRangePotentials();
    updateShortRangePotentials();
    updateSelfPotentials();
    updateExcludedPotentials();
    updateRestraintPotentials();
    if (doPotentialSelection_) updateSelectionPotentials();
  }

  void SPFForceManager::updateLongRangePotentials() {
    potVec longRangePotentials =
        linearCombination(temporarySourceSnapshot_->getLongRangePotentials(),
                          temporarySinkSnapshot_->getLongRangePotentials());
    currentSnapshot_->setLongRangePotentials(longRangePotentials);

    RealType reciprocalPotential =
        linearCombination(temporarySourceSnapshot_->getReciprocalPotential(),
                          temporarySinkSnapshot_->getReciprocalPotential());
    currentSnapshot_->setReciprocalPotential(reciprocalPotential);

    RealType surfacePotential =
        linearCombination(temporarySourceSnapshot_->getSurfacePotential(),
                          temporarySinkSnapshot_->getSurfacePotential());
    currentSnapshot_->setSurfacePotential(surfacePotential);
  }

  void SPFForceManager::updateShortRangePotentials() {
    RealType bondPotential =
        linearCombination(temporarySourceSnapshot_->getBondPotential(),
                          temporarySinkSnapshot_->getBondPotential());
    currentSnapshot_->setBondPotential(bondPotential);

    RealType bendPotential =
        linearCombination(temporarySourceSnapshot_->getBendPotential(),
                          temporarySinkSnapshot_->getBendPotential());
    currentSnapshot_->setBendPotential(bendPotential);

    RealType torsionPotential =
        linearCombination(temporarySourceSnapshot_->getTorsionPotential(),
                          temporarySinkSnapshot_->getTorsionPotential());
    currentSnapshot_->setTorsionPotential(torsionPotential);

    RealType inversionPotential =
        linearCombination(temporarySourceSnapshot_->getInversionPotential(),
                          temporarySinkSnapshot_->getInversionPotential());
    currentSnapshot_->setInversionPotential(inversionPotential);
  }

  void SPFForceManager::updateSelfPotentials() {
    potVec selfPotentials =
        linearCombination(temporarySourceSnapshot_->getSelfPotentials(),
                          temporarySinkSnapshot_->getSelfPotentials());
    currentSnapshot_->setSelfPotentials(selfPotentials);
  }

  void SPFForceManager::updateExcludedPotentials() {
    potVec excludedPotentials =
        linearCombination(temporarySourceSnapshot_->getExcludedPotentials(),
                          temporarySinkSnapshot_->getExcludedPotentials());
    currentSnapshot_->setExcludedPotentials(excludedPotentials);
  }

  void SPFForceManager::updateRestraintPotentials() {
    RealType restraintPotential =
        linearCombination(temporarySourceSnapshot_->getRestraintPotential(),
                          temporarySinkSnapshot_->getRestraintPotential());
    currentSnapshot_->setRestraintPotential(restraintPotential);
  }

  void SPFForceManager::updateSelectionPotentials() {
    potVec selectionPotentials =
        linearCombination(temporarySourceSnapshot_->getSelectionPotentials(),
                          temporarySinkSnapshot_->getSelectionPotentials());
    currentSnapshot_->setSelectionPotentials(selectionPotentials);
  }

  void SPFForceManager::updateVirialTensor() {
    Mat3x3d virialTensor =
        linearCombination(temporarySourceSnapshot_->getVirialTensor(),
                          temporarySinkSnapshot_->getVirialTensor());
    currentSnapshot_->setVirialTensor(virialTensor);

    Mat3x3d pressureTensor =
        linearCombination(thermo_->getPressureTensor(temporarySourceSnapshot_),
                          thermo_->getPressureTensor(temporarySinkSnapshot_));
    currentSnapshot_->setPressureTensor(pressureTensor);

    RealType pressure =
        linearCombination(thermo_->getPressure(temporarySourceSnapshot_),
                          thermo_->getPressure(temporarySinkSnapshot_));
    currentSnapshot_->setPressure(pressure);
  }
}  // namespace OpenMD::RNEMD
