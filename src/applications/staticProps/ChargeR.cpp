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

#include "applications/staticProps/ChargeR.hpp"

#include <algorithm>
#include <fstream>

#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  ChargeR::ChargeR(SimInfo* info, const std::string& filename,
                   const std::string& sele, RealType len, int nrbins) :
      StaticAnalyser(info, filename, nrbins),
      selectionScript_(sele), evaluator_(info), seleMan_(info), thermo_(info),
      len_(len) {
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    deltaR_ = len_ / nBins_;

    // fixed number of bins

    sliceSDLists_.resize(nBins_);
    sliceSDCount_.resize(nBins_);
    std::fill(sliceSDCount_.begin(), sliceSDCount_.end(), 0);

    chargeR_.resize(nBins_);
    setOutputName(getPrefix(filename) + ".ChargeR");
    std::stringstream params;
    params << " len = " << len_ << ", nrbins = " << nBins_;
    const std::string paramString = params.str();
    setParameterString(paramString);
  }

  void ChargeR::process() {
    StuntDouble* sd;
    int ii;

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames / step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for (unsigned int i = 0; i < nBins_; i++) {
        sliceSDLists_[i].clear();
      }

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      // determine which atom belongs to which slice
      for (sd = seleMan_.beginSelected(ii); sd != NULL;
           sd = seleMan_.nextSelected(ii)) {
        Vector3d pos      = sd->getPos();
        RealType distance = pos.length();

        if (distance < len_) {
          int binNo = int(distance / deltaR_);
          sliceSDLists_[binNo].push_back(sd);
          sliceSDCount_[binNo]++;
        }
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
        chargeR_[i] += binC;
        // Units of (e / Ang^2 / fs)
      }
    }

    writeChargeR();
  }

  void ChargeR::writeChargeR() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#ChargeR "
                << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "# r "
                << "\tcharge\n";
      RealType binCharge;
      for (unsigned int i = 0; i < chargeR_.size(); ++i) {
        RealType rLower = i * deltaR_;
        RealType rUpper = rLower + deltaR_;
        RealType volShell =
            (4.0 * Constants::PI) * (pow(rUpper, 3) - pow(rLower, 3)) / 3.0;

        RealType r = deltaR_ * (i + 0.5);

        binCharge = chargeR_[i] / (volShell * nProcessed_);

        rdfStream << r << "\t" << binCharge << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "ChargeR: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }
}  // namespace OpenMD
