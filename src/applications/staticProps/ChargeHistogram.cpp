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

/*
 * Calculates average charge profile for selected atom.
 * Created by Hemanta Bhattarai on 05/10/19.
 * @author  Hemanta Bhattarai
 */

#include "applications/staticProps/ChargeHistogram.hpp"

#include <algorithm>
#include <fstream>
#include <numeric>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  ChargeHistogram::ChargeHistogram(SimInfo* info, const std::string& filename,
                                   const std::string& sele, int nbins) :
      StaticAnalyser(info, filename, nbins),
      selectionScript_(sele), evaluator_(info), seleMan_(info), nBins_(nbins) {
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    setOutputName(getPrefix(filename) + ".Chargehist");
  }

  void ChargeHistogram::process() {
    StuntDouble* sd;
    int ii;

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames / step_;
    vector<RealType> charge;

    nProcessed_ = nFrames / step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for (sd = seleMan_.beginSelected(ii); sd != NULL;
           sd = seleMan_.nextSelected(ii)) {
        RealType q = 0.0;
        Atom* atom = static_cast<Atom*>(sd);

        AtomType* atomType = atom->getAtomType();

        FixedChargeAdapter fca = FixedChargeAdapter(atomType);
        if (fca.isFixedCharge()) { q += fca.getCharge(); }

        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
        if (fqa.isFluctuatingCharge()) { q += atom->getFlucQPos(); }

        charge.push_back(q);
      }
    }

    if (charge.empty()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Selected atom not found.\n");
      painCave.isFatal = 1;
      simError();
    }

    std::sort(charge.begin(), charge.end());

    RealType min = charge.front();
    RealType max = charge.back();

    RealType delta_charge = (max - min) / (nBins_);

    if (delta_charge == 0) {
      bincenter_.push_back(min);
      histList_.push_back(1);
    } else {
      // fill the center for histogram
      for (int j = 0; j < nBins_ + 3; ++j) {
        bincenter_.push_back(min + (j - 1) * delta_charge);
        histList_.push_back(0);
      }
      // filling up the histogram whith the densities
      int bin_center_pos = 0;
      vector<RealType>::iterator index;
      RealType charge_length = static_cast<RealType>(charge.size());

      bool hist_update;
      for (index = charge.begin(); index < charge.end(); index++) {
        hist_update = true;
        while (hist_update) {
          if (*index >= bincenter_[bin_center_pos] &&
              *index < bincenter_[bin_center_pos + 1]) {
            histList_[bin_center_pos] += 1.0 / charge_length;
            hist_update = false;
          } else {
            bin_center_pos++;
            hist_update = true;
          }
        }
      }
    }
    writeCharge();
  }

  void ChargeHistogram::writeCharge() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#Charges for selection\n";
      rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#"
                << "Bin_center"
                << "\tcount\n";
      for (unsigned int i = 0; i < histList_.size(); ++i) {
        rdfStream << bincenter_[i] << "\t" << histList_[i] << "\n";
      }
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "ChargeHistogram: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    rdfStream.close();
  }
}  // namespace OpenMD
