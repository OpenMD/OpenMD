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

/**
 * @file StuntDouble.cpp
 * @author    tlin
 * @date  10/22/2004
 * @version 1.0
 */

#include "primitives/StuntDouble.hpp"

#include <cassert>
#include <memory>

namespace OpenMD {

  StuntDouble::StuntDouble(ObjectType objType, DataStoragePointer storage) :
      objType_(objType), storage_(storage), snapshotMan_(NULL), linear_(false),
      linearAxis_(-1), globalIndex_(-1), globalIntegrableObjectIndex_(-1),
      localIndex_(-1) {}

  void StuntDouble::zeroForcesAndTorques() {
    int sl = (snapshotMan_->getCurrentSnapshot()->*storage_).getStorageLayout();

    if (sl & DataStorage::dslForce) setFrc(V3Zero);
    if (sl & DataStorage::dslTorque) setTrq(V3Zero);
    if (sl & DataStorage::dslParticlePot) setParticlePot(0.0);
    if (sl & DataStorage::dslFlucQForce) setFlucQFrc(0.0);
    if (sl & DataStorage::dslElectricField) setElectricField(V3Zero);
    if (sl & DataStorage::dslSitePotential) setSitePotential(0.0);
  }

  void StuntDouble::combineForcesAndTorques(Snapshot* snapA, Snapshot* snapB,
                                            RealType multA, RealType multB) {
    int sl = (snapshotMan_->getCurrentSnapshot()->*storage_).getStorageLayout();

    assert(sl == (snapA->*storage_).getStorageLayout());
    assert(sl == (snapB->*storage_).getStorageLayout());

    if (sl & DataStorage::dslForce) {
      Vector3d forceA = (snapA->*storage_).force[localIndex_];
      Vector3d forceB = (snapB->*storage_).force[localIndex_];

      ((snapshotMan_->getCurrentSnapshot())->*storage_).force[localIndex_] =
          multA * forceA + multB * forceB;
    }

    if (sl & DataStorage::dslTorque) {
      Vector3d torqueA = (snapA->*storage_).torque[localIndex_];
      Vector3d torqueB = (snapB->*storage_).torque[localIndex_];

      ((snapshotMan_->getCurrentSnapshot())->*storage_).torque[localIndex_] =
          multA * torqueA + multB * torqueB;
    }

    if (sl & DataStorage::dslParticlePot) {
      RealType particlePotA = (snapA->*storage_).particlePot[localIndex_];
      RealType particlePotB = (snapB->*storage_).particlePot[localIndex_];

      ((snapshotMan_->getCurrentSnapshot())->*storage_)
          .particlePot[localIndex_] =
          multA * particlePotA + multB * particlePotB;
    }

    if (sl & DataStorage::dslFlucQForce) {
      RealType flucQFrcA = (snapA->*storage_).flucQFrc[localIndex_];
      RealType flucQFrcB = (snapB->*storage_).flucQFrc[localIndex_];

      ((snapshotMan_->getCurrentSnapshot())->*storage_).flucQFrc[localIndex_] =
          multA * flucQFrcA + multB * flucQFrcB;
    }

    if (sl & DataStorage::dslElectricField) {
      Vector3d electricFieldA = (snapA->*storage_).electricField[localIndex_];
      Vector3d electricFieldB = (snapB->*storage_).electricField[localIndex_];

      ((snapshotMan_->getCurrentSnapshot())->*storage_)
          .electricField[localIndex_] =
          multA * electricFieldA + multB * electricFieldB;
    }

    if (sl & DataStorage::dslSitePotential) {
      RealType sitePotentialA = (snapA->*storage_).sitePotential[localIndex_];
      RealType sitePotentialB = (snapB->*storage_).sitePotential[localIndex_];

      ((snapshotMan_->getCurrentSnapshot())->*storage_)
          .sitePotential[localIndex_] =
          multA * sitePotentialA + multB * sitePotentialB;
    }
  }

  void StuntDouble::addProperty(std::shared_ptr<GenericData> genData) {
    properties_.addProperty(genData);
  }

  void StuntDouble::removeProperty(const std::string& propName) {
    properties_.removeProperty(propName);
  }

  std::vector<std::string> StuntDouble::getPropertyNames() {
    return properties_.getPropertyNames();
  }

  std::vector<std::shared_ptr<GenericData>> StuntDouble::getProperties() {
    return properties_.getProperties();
  }

  std::shared_ptr<GenericData> StuntDouble::getPropertyByName(
      const std::string& propName) {
    return properties_.getPropertyByName(propName);
  }

}  // namespace OpenMD
