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

#ifndef OPENMD_RNEMD_SPFFORCEMANAGER_HPP
#define OPENMD_RNEMD_SPFFORCEMANAGER_HPP

#include <config.h>

#include <cmath>

#include "brains/ForceManager.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Snapshot.hpp"
#include "math/Vector3.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD::RNEMD {

  class SPFForceManager : public ForceManager {
  public:
    SPFForceManager(SimInfo* info);

    void setSelectedMolecule(Molecule* selectedMolecule, Vector3d newCom);
    bool updateLambda(RealType& particleTarget, RealType& deltaLambda);

    RealType getScaledDeltaU(RealType d_lambda) const {
      return -(f_lambda(lambda_ + d_lambda) - f_lambda(lambda_)) *
             (potentialSink_ - potentialSource_);
    }

    Molecule* getSelectedMolecule() const { return selectedMolecule_; }
    Snapshot* getTemporarySourceSnapshot() const {
      return temporarySourceSnapshot_;
    }
    Snapshot* getTemporarySinkSnapshot() const {
      return temporarySinkSnapshot_;
    }

    void combineForcesAndTorques();
    void updatePotentials();
    void updateVirialTensor();

    RealType f_lambda(RealType lambda) const { return std::pow(lambda, k_); }

  private:
    void calcForces() override;

    void updateLongRangePotentials();
    void updateShortRangePotentials();
    void updateSelfPotentials();
    void updateExcludedPotentials();
    void updateRestraintPotentials();
    void updateSelectionPotentials();

    template<typename T>
    T linearCombination(T quantityA, T quantityB) {
      RealType result = f_lambda(lambda_);

      return T {(1.0 - result) * quantityA + result * quantityB};
    }

    Snapshot* currentSnapshot_ {nullptr};
    Snapshot* temporarySourceSnapshot_ {nullptr};
    Snapshot* temporarySinkSnapshot_ {nullptr};

    Molecule* selectedMolecule_ {nullptr};
    Vector3d currentSinkCom_ {};

    RealType lambda_ {};
    const unsigned int k_ {3};

    RealType potentialSource_ {};
    RealType potentialSink_ {};
  };
}  // namespace OpenMD::RNEMD

#endif  // OPENMD_RNEMD_SPFFORCEMANAGER_HPP
