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

#ifndef OPENMD_RNEMD_SPFFORCEMANAGER_HPP
#define OPENMD_RNEMD_SPFFORCEMANAGER_HPP

#include <config.h>

#include <cmath>
#include <memory>

#include "brains/ForceManager.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Snapshot.hpp"
#include "brains/Thermo.hpp"
#include "math/Vector3.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD::RNEMD {

  class SPFMethod;

  class SPFForceManager : public ForceManager {
  public:
    SPFForceManager(SimInfo* info);
    ~SPFForceManager();

    void setSelectedMolecule(Molecule* selectedMolecule);
    void setDeltaLambda(RealType spfTarget);
    RealType getScaledDeltaU();

    bool getHasSelectedMolecule() const { return hasSelectedMolecule_; }

    void setHasSelectedMolecule(bool hasSelectedMolecule) {
      hasSelectedMolecule_ = hasSelectedMolecule;
    }

    Molecule* getSelectedMolecule() { return selectedMolecule_; }
    Snapshot getTemporarySourceSnapshot() { return *temporarySourceSnapshot_; }
    Snapshot getTemporarySinkSnapshot() { return *temporarySinkSnapshot_; }

    void updateSPFState() {
      combineForcesAndTorques();
      updatePotentials();
      updateVirialTensor();
    }

    SPFMethod* spfRNEMD_;

  private:
    void calcForces() override;

    void combineForcesAndTorques();
    void updatePotentials();
    void updateLongRangePotentials();
    void updateShortRangePotentials();
    void updateSelfPotentials();
    void updateExcludedPotentials();
    void updateRestraintPotentials();
    void updateSelectionPotentials();
    void updateVirialTensor();

    RealType f_lambda(RealType lambda) const {
      RealType result = std::pow(lambda, k_);

      return std::clamp(result, 0.0, 1.0);
    }

    template<typename T>
    T linearCombination(T quantityA, T quantityB) {
      RealType result = f_lambda(currentSnapshot_->getSPFData()->lambda);

      return T {(1.0 - result) * quantityA + result * quantityB};
    }

    std::unique_ptr<Thermo> thermo_ {nullptr};

    Snapshot* currentSnapshot_ {nullptr};
    Snapshot* temporarySourceSnapshot_ {nullptr};
    Snapshot* temporarySinkSnapshot_ {nullptr};

    // to preserve the Verlet Neighbor lists in source and sink snapshots:
    std::vector<int> sourceNeighborList_;
    std::vector<int> sourcePoint_;
    std::vector<Vector3d> sourceSavedPositions_;

    std::vector<int> sinkNeighborList_;
    std::vector<int> sinkPoint_;
    std::vector<Vector3d> sinkSavedPositions_;

    bool hasSelectedMolecule_ {};

    Molecule* selectedMolecule_ {nullptr};
    int k_ {};

    RealType deltaLambda_ {};

    RealType potentialSource_ {};
    RealType potentialSink_ {};
  };
}  // namespace OpenMD::RNEMD

#endif  // OPENMD_RNEMD_SPFFORCEMANAGER_HPP
