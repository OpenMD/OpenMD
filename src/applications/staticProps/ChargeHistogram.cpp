/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
                                 const std::string& sele, int nbins)
    : StaticAnalyser(info, filename, nbins),
      selectionScript_(sele),
      evaluator_(info),
      seleMan_(info),
      nBins_(nbins) {
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
      if (fca.isFixedCharge()) {
        q += fca.getCharge();
      }

      FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
      if (fqa.isFluctuatingCharge()) {
        q += atom->getFlucQPos();
      }

      charge.push_back(q);
    }
  }

  if (charge.empty()) {
    sprintf(painCave.errMsg, "Selected atom not found.\n");
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
    sprintf(painCave.errMsg, "ChargeHistogram: unable to open %s\n",
            outputFilename_.c_str());
    painCave.isFatal = 1;
    simError();
  }
  rdfStream.close();
}
}  // namespace OpenMD
