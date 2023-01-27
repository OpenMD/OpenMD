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

/*
 * Computes the charge density distribution along preferred axis for the
 * selected atom Created by Cody R. Drisko on 06/14/14.
 */

#include "applications/staticProps/ChargeZ.hpp"

#include <algorithm>
#include <fstream>

#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  ChargeZ::ChargeZ(SimInfo* info, const std::string& filename,
                   const std::string& sele, int nzbins, int axis) :
      StaticAnalyser(info, filename, nzbins),
      selectionScript_(sele), evaluator_(info), seleMan_(info), thermo_(info),
      axis_(axis) {
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // fixed number of bins

    sliceSDLists_.resize(nBins_);
    sliceSDCount_.resize(nBins_);
    std::fill(sliceSDCount_.begin(), sliceSDCount_.end(), 0);

    chargeZ_.resize(nBins_);

    switch (axis_) {
    case 0:
      axisLabel_ = "x";
      break;
    case 1:
      axisLabel_ = "y";
      break;
    case 2:
    default:
      axisLabel_ = "z";
      break;
    }

    setOutputName(getPrefix(filename) + ".ChargeZ");
  }

  void ChargeZ::process() {
    StuntDouble* sd;
    int ii;

    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames / step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for (unsigned int i = 0; i < nBins_; i++) {
        sliceSDLists_[i].clear();
      }

      Mat3x3d hmat = currentSnapshot_->getHmat();
      zBox_.push_back(hmat(axis_, axis_));

      RealType halfBoxZ_ = hmat(axis_, axis_) / 2.0;
      RealType area      = 0.0;
      switch (axis_) {
      case 0:
        area = currentSnapshot_->getYZarea();
        break;
      case 1:
        area = currentSnapshot_->getXZarea();
        break;
      case 2:
      default:
        area = currentSnapshot_->getXYarea();
        break;
      }

      areas_.push_back(area);

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      // wrap the stuntdoubles into a cell
      for (sd = seleMan_.beginSelected(ii); sd != NULL;
           sd = seleMan_.nextSelected(ii)) {
        Vector3d pos = sd->getPos();
        if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(pos);
        sd->setPos(pos);
      }

      // determine which atom belongs to which slice
      for (sd = seleMan_.beginSelected(ii); sd != NULL;
           sd = seleMan_.nextSelected(ii)) {
        Vector3d pos = sd->getPos();
        // shift molecules by half a box to have bins start at 0
        int binNo = int(nBins_ * (halfBoxZ_ + pos[axis_]) / hmat(axis_, axis_));
        sliceSDLists_[binNo].push_back(sd);
        sliceSDCount_[binNo]++;
      }

      // loop over the slices to calculate the charge
      for (unsigned int i = 0; i < nBins_; i++) {
        RealType binC = 0;
        for (unsigned int k = 0; k < sliceSDLists_[i].size(); ++k) {
          RealType q = 0.0;
          Atom* atom = static_cast<Atom*>(sliceSDLists_[i][k]);

          AtomType* atomType = atom->getAtomType();

          if (sliceSDLists_[i][k]->isAtom()) {
            FixedChargeAdapter fca = FixedChargeAdapter(atomType);
            if (fca.isFixedCharge()) { q += fca.getCharge(); }

            FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
            if (fqa.isFluctuatingCharge()) { q += atom->getFlucQPos(); }
          }

          binC += q;
        }
        chargeZ_[i] += binC;
        // Units of (e / Ang^2 / fs)
      }
    }

    writeChargeZ();
  }

  void ChargeZ::writeChargeZ() {
    // compute average box length:
    std::vector<RealType>::iterator j;
    RealType zSum = 0.0;
    for (j = zBox_.begin(); j != zBox_.end(); ++j) {
      zSum += *j;
    }
    RealType zAve = zSum / zBox_.size();

    RealType areaSum = 0.0;
    for (j = areas_.begin(); j != areas_.end(); ++j) {
      areaSum += *j;
    }
    RealType areaAve = areaSum / areas_.size();

    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#ChargeZ "
                << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#" << axisLabel_ << "\tcharge\n";
      RealType binCharge;
      for (unsigned int i = 0; i < chargeZ_.size(); ++i) {
        RealType z = zAve * (i + 0.5) / chargeZ_.size();

        RealType volSlice = areaAve * zAve / zBox_.size();

        binCharge = chargeZ_[i] / (volSlice * nProcessed_);

        rdfStream << z << "\t" << binCharge << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "ChargeZ: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }
}  // namespace OpenMD
