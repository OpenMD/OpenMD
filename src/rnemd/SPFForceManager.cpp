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

#include "rnemd/SPFForceManager.hpp"

namespace OpenMD {
  namespace RNEMD {

    SPFForceManager::SPFForceManager(SimInfo* info) :
        ForceManager(info), lambda_ {0.0} {
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      int nAtoms        = info_->getNAtoms();
      int nRigidBodies  = info_->getNRigidBodies();
      int nCutoffGroups = info_->getNCutoffGroups();
      int storageLayout = info_->getSnapshotManager()->getStorageLayout();

      bool usePBC = info_->getSimParams()->getUsePeriodicBoundaryConditions();

      ghostSnapshot_  = new Snapshot(nAtoms, nRigidBodies, nCutoffGroups,
                                    storageLayout, usePBC);
      *ghostSnapshot_ = *currentSnapshot_;
      ghostSnapshot_->clearDerivedProperties();
    }

    SPFForceManager::~SPFForceManager() { delete ghostSnapshot_; }

    void SPFForceManager::calcForces() {
      // current is at t+dt,  previous is at t, ghost is at t
      // we want forces from current(t+dt) and ghost(t+dt)
      int index {};
      StuntDouble* sd;
      Molecule::IntegrableObjectIterator j;

      std::vector<Vector3d> deltas;
      std::vector<Vector3d> previous;

      ForceManager::calcForces();  // current snapshot
      potentialA_ = currentSnapshot_->getPotentialEnergy();

      if (selectedMolecule_) {
        int nAtoms        = info_->getNAtoms();
        int nRigidBodies  = info_->getNRigidBodies();
        int nCutoffGroups = info_->getNCutoffGroups();
        int storageLayout = info_->getSnapshotManager()->getStorageLayout();

        bool usePBC = info_->getSimParams()->getUsePeriodicBoundaryConditions();

        Snapshot* tempSnapshot = new Snapshot(
            nAtoms, nRigidBodies, nCutoffGroups, storageLayout, usePBC);

        // save delta from propagated current and previous:
        for (sd = selectedMolecule_->beginIntegrableObject(j); sd != NULL;
             sd = selectedMolecule_->nextIntegrableObject(j)) {
          deltas.push_back(sd->getPos() - sd->getPrevPos());
        }

        *tempSnapshot     = *currentSnapshot_;
        *currentSnapshot_ = *ghostSnapshot_;

        currentSnapshot_->clearDerivedProperties();

        // get position of tagged molecule from previous ghost
        for (sd = selectedMolecule_->beginIntegrableObject(j); sd != NULL;
             sd = selectedMolecule_->nextIntegrableObject(j)) {
          previous.push_back(sd->getPos());
        }

        // propagate the delta in position into the ghost
        for (sd = selectedMolecule_->beginIntegrableObject(j); sd != NULL;
             sd = selectedMolecule_->nextIntegrableObject(j)) {
          sd->setPos(previous[index] + deltas[index]);
          index++;
        }

        ForceManager::calcForces();  // ghost snapshot
        potentialB_ = currentSnapshot_->getPotentialEnergy();

        *ghostSnapshot_   = *currentSnapshot_;
        *currentSnapshot_ = *tempSnapshot;

        currentSnapshot_->clearDerivedProperties();

        // calculate lambda-averaged forces on all atoms and potentials
        SimInfo::MoleculeIterator mi;
        Molecule* mol;
        Molecule::IntegrableObjectIterator ii;

        // now scale forces and torques of all the sds
        for (mol = info_->beginMolecule(mi); mol != NULL;
             mol = info_->nextMolecule(mi)) {
          for (sd = mol->beginIntegrableObject(ii); sd != NULL;
               sd = mol->nextIntegrableObject(ii)) {
            sd->combineForcesAndTorques(currentSnapshot_, ghostSnapshot_,
                                        1 - lambda_, lambda_);
          }
        }

        updatePotentials();
        updateVirialTensor();

        delete tempSnapshot;
      } else {
        potentialB_ = currentSnapshot_->getPotentialEnergy();
      }
    }

    void SPFForceManager::setSelectedMolecule(Molecule* selectedMolecule,
                                              Vector3d newPosition) {
      Snapshot tempSnapshot = *currentSnapshot_;
      *currentSnapshot_     = *ghostSnapshot_;

      if (selectedMolecule) {
        // tagged molecule should be pushed from rnemd
        Vector3d delta = newPosition - selectedMolecule->getCom();
        selectedMolecule->moveCom(delta);
        selectedMolecule_ = selectedMolecule;
      } else {
        selectedMolecule_ = NULL;
      }
      *ghostSnapshot_   = *currentSnapshot_;
      *currentSnapshot_ = tempSnapshot;
    }

    void SPFForceManager::updatePotentials() {
      updateLongRangePotentials();
      updateShortRangePotentials();
      updateSelfPotentials();
      updateExcludedPotentials();
      updateRestraintPotentials();
      if (doPotentialSelection_) updateSelectionPotentials();

      // Defer potentialEnergy calculation until it is needed
      currentSnapshot_->hasPotentialEnergy = false;
    }

    void SPFForceManager::updateLongRangePotentials() {
      potVec longRangePotentials =
          linearCombination(currentSnapshot_->getLongRangePotentials(),
                            ghostSnapshot_->getLongRangePotentials());
      currentSnapshot_->setLongRangePotentials(longRangePotentials);

      RealType reciprocalPotential =
          linearCombination(currentSnapshot_->getReciprocalPotential(),
                            ghostSnapshot_->getReciprocalPotential());
      currentSnapshot_->setReciprocalPotential(reciprocalPotential);

      RealType surfacePotential =
          linearCombination(currentSnapshot_->getSurfacePotential(),
                            ghostSnapshot_->getSurfacePotential());
      currentSnapshot_->setSurfacePotential(surfacePotential);

      // Defer longRangePotential calculation until it is needed
      currentSnapshot_->hasLongRangePotential = false;
    }

    void SPFForceManager::updateShortRangePotentials() {
      RealType bondPotential =
          linearCombination(currentSnapshot_->getBondPotential(),
                            ghostSnapshot_->getBondPotential());
      currentSnapshot_->setBondPotential(bondPotential);

      RealType bendPotential =
          linearCombination(currentSnapshot_->getBendPotential(),
                            ghostSnapshot_->getBendPotential());
      currentSnapshot_->setBendPotential(bendPotential);

      RealType torsionPotential =
          linearCombination(currentSnapshot_->getTorsionPotential(),
                            ghostSnapshot_->getTorsionPotential());
      currentSnapshot_->setTorsionPotential(torsionPotential);

      RealType inversionPotential =
          linearCombination(currentSnapshot_->getInversionPotential(),
                            ghostSnapshot_->getInversionPotential());
      currentSnapshot_->setInversionPotential(inversionPotential);

      // Defer shortRangePotential calculation until it is needed
      currentSnapshot_->hasShortRangePotential = false;
    }

    void SPFForceManager::updateSelfPotentials() {
      potVec selfPotentials =
          linearCombination(currentSnapshot_->getSelfPotentials(),
                            ghostSnapshot_->getSelfPotentials());
      currentSnapshot_->setSelfPotentials(selfPotentials);

      // Defer selfPotential calculation until it is needed
      currentSnapshot_->hasSelfPotential = false;
    }

    void SPFForceManager::updateExcludedPotentials() {
      potVec excludedPotentials =
          linearCombination(currentSnapshot_->getExcludedPotentials(),
                            ghostSnapshot_->getExcludedPotentials());
      currentSnapshot_->setExcludedPotentials(excludedPotentials);

      currentSnapshot_->hasExcludedPotential = false;
    }

    void SPFForceManager::updateRestraintPotentials() {
      RealType restraintPotential =
          linearCombination(currentSnapshot_->getRestraintPotential(),
                            ghostSnapshot_->getRestraintPotential());
      currentSnapshot_->setRestraintPotential(restraintPotential);
    }

    void SPFForceManager::updateSelectionPotentials() {
      potVec selectionPotentials =
          linearCombination(currentSnapshot_->getSelectionPotentials(),
                            ghostSnapshot_->getSelectionPotentials());
      currentSnapshot_->setSelectionPotentials(selectionPotentials);
    }

    void SPFForceManager::updateVirialTensor() {
      Mat3x3d virialTensor =
          linearCombination(currentSnapshot_->getVirialTensor(),
                            ghostSnapshot_->getVirialTensor());
      currentSnapshot_->setVirialTensor(virialTensor);

      // Defer pressure and pressureTensor calculations until they are needed
      currentSnapshot_->hasPressureTensor = false;
      currentSnapshot_->hasPressure       = false;
    }
  }  // namespace RNEMD
}  // namespace OpenMD