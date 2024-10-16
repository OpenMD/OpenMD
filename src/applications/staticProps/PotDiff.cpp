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

#include "applications/staticProps/PotDiff.hpp"

#include <algorithm>
#include <functional>

#include "brains/ForceManager.hpp"
#include "brains/SimSnapshotManager.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  PotDiff::PotDiff(SimInfo* info, const std::string& filename,
                   const std::string& sele) :
      StaticAnalyser(info, filename, 1),
      selectionScript_(sele), seleMan_(info), evaluator_(info) {
    StuntDouble* sd;
    int i;

    setOutputName(getPrefix(filename) + ".potDiff");

    // The PotDiff is computed by negating the charge on the atom type
    // using fluctuating charge values.  If we don't have any
    // fluctuating charges in the simulation, we need to expand
    // storage to hold them.
    int atomStorageLayout        = info_->getAtomStorageLayout();
    int rigidBodyStorageLayout   = info->getRigidBodyStorageLayout();
    int cutoffGroupStorageLayout = info->getCutoffGroupStorageLayout();

    atomStorageLayout |= DataStorage::dslFlucQPosition;
    atomStorageLayout |= DataStorage::dslFlucQVelocity;
    atomStorageLayout |= DataStorage::dslFlucQForce;

    info_->setAtomStorageLayout(atomStorageLayout);
    info_->setSnapshotManager(new SimSnapshotManager(info_, atomStorageLayout,
                                                     rigidBodyStorageLayout,
                                                     cutoffGroupStorageLayout));

    // now we have to figure out which AtomTypes to convert to fluctuating
    // charges
    evaluator_.loadScriptString(sele);
    seleMan_.setSelectionSet(evaluator_.evaluate());
    for (sd = seleMan_.beginSelected(i); sd != NULL;
         sd = seleMan_.nextSelected(i)) {
      AtomType* at                 = static_cast<Atom*>(sd)->getAtomType();
      FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(at);
      if (fqa.isFluctuatingCharge()) {
        selectionWasFlucQ_.push_back(true);
      } else {
        selectionWasFlucQ_.push_back(false);
        // make a fictitious fluctuating charge with an unphysical
        // charge mass and slaterN, but we need to zero out the
        // electronegativity and hardness to remove the self
        // contribution:
        fqa.makeFluctuatingCharge(1.0e9, 0.0, 0.0, 1);
        sd->setFlucQPos(0.0);
      }
    }
    info_->getSnapshotManager()->advance();
  }

  void PotDiff::process() {
    StuntDouble* sd;
    int j;

    diff_.clear();
    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();

    // We'll need the force manager to compute the potential
    ForceManager* forceMan = new ForceManager(info_);

    // We'll need thermo to report the potential
    Thermo* thermo = new Thermo(info_);

    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for (sd = seleMan_.beginSelected(j); sd != NULL;
           sd = seleMan_.nextSelected(j)) {
        if (!selectionWasFlucQ_[j]) { sd->setFlucQPos(0.0); }
      }

      forceMan->calcForces();
      RealType pot1 = thermo->getPotential();

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      for (sd = seleMan_.beginSelected(j); sd != NULL;
           sd = seleMan_.nextSelected(j)) {
        AtomType* at = static_cast<Atom*>(sd)->getAtomType();

        FixedChargeAdapter fca       = FixedChargeAdapter(at);
        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(at);

        RealType charge = 0.0;

        if (fca.isFixedCharge()) charge += fca.getCharge();
        if (fqa.isFluctuatingCharge()) charge += sd->getFlucQPos();

        sd->setFlucQPos(-charge);
      }

      currentSnapshot_->clearDerivedProperties();
      forceMan->calcForces();
      RealType pot2 = thermo->getPotential();
      RealType diff = pot2 - pot1;

      data_.add(diff);
      diff_.push_back(diff);
      times_.push_back(currentSnapshot_->getTime());

      info_->getSnapshotManager()->advance();
    }

    writeDiff();
  }

  void PotDiff::writeDiff() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);
    if (ofs.is_open()) {
      RealType mu    = data_.getAverage();
      RealType sigma = data_.getStdDev();
      RealType m95   = data_.get95percentConfidenceInterval();

      ofs << "#potDiff\n";
      ofs << "#selection: (" << selectionScript_ << ")\n";
      ofs << "# <diff> = " << mu << "\n";
      ofs << "# StdDev = " << sigma << "\n";
      ofs << "# 95% confidence interval = " << m95 << "\n";
      ofs << "# t\tdiff[t]\n";
      for (unsigned int i = 0; i < diff_.size(); ++i) {
        ofs << times_[i] << "\t" << diff_[i] << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "PotDiff: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    ofs.close();
  }

}  // namespace OpenMD
