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

#include "applications/dynamicProps/CurrentDensityAutoCorrFunc.hpp"

#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/Revision.hpp"

using namespace std;
namespace OpenMD {
  CurrentDensityAutoCorrFunc::CurrentDensityAutoCorrFunc(SimInfo* info,
                                                         const string& filename,
                                                         const string& sele1,
                                                         const string& sele2) :
      SystemACF<RealType>(info, filename, sele1, sele2) {
    setCorrFuncType("Current Density Auto Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".currentDensityCorr");

    AtomTypeSet osTypes = seleMan1_.getSelectedAtomTypes();
    std::copy(osTypes.begin(), osTypes.end(), std::back_inserter(outputTypes_));

    Jc_.resize(nFrames_, V3Zero);
    JcCount_.resize(nFrames_, 0);

    typeJc_.resize(nFrames_);
    typeCounts_.resize(nFrames_);
    myHistogram_.resize(nTimeBins_);

    for (int i = 0; i < nFrames_; ++i) {
      typeJc_[i].resize(outputTypes_.size(), V3Zero);
      typeCounts_[i].resize(outputTypes_.size(), 0);
    }
    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      myHistogram_[i].resize(outputTypes_.size() + 1, 0.0);
    }
    // We'll need thermo to compute the volume:
    thermo_ = new Thermo(info_);
  }

  void CurrentDensityAutoCorrFunc::computeProperty1(int frame) {
    StuntDouble* sd1;
    AtomType* atype;
    std::vector<AtomType*>::iterator at;
    int i;

    for (sd1 = seleMan1_.beginSelected(i); sd1 != NULL;
         sd1 = seleMan1_.nextSelected(i)) {
      Vector3d v = sd1->getVel();
      RealType q = 0.0;
      int typeIndex(-1);

      if (sd1->isAtom()) {
        atype                  = static_cast<Atom*>(sd1)->getAtomType();
        FixedChargeAdapter fca = FixedChargeAdapter(atype);
        if (fca.isFixedCharge()) q = fca.getCharge();
        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
        if (fqa.isFluctuatingCharge()) q += sd1->getFlucQPos();

        typeIndex = -1;
        at        = std::find(outputTypes_.begin(), outputTypes_.end(), atype);
        if (at != outputTypes_.end()) {
          typeIndex = std::distance(outputTypes_.begin(), at);
        }
        if (typeIndex != -1) {
          typeCounts_[frame][typeIndex]++;
          typeJc_[frame][typeIndex] += q * v;
        }
      }
      JcCount_[frame]++;
      Jc_[frame] += q * v;
    }

    RealType vol = thermo_->getVolume();

    Jc_[frame] /= (vol * Constants::currentDensityConvert);
    for (unsigned int j = 0; j < outputTypes_.size(); j++) {
      typeJc_[frame][j] /= (vol * Constants::currentDensityConvert);
    }
  }

  void CurrentDensityAutoCorrFunc::correlateFrames(int frame1, int frame2,
                                                   int timeBin) {
    RealType corrVal(0.0);
    corrVal = dot(Jc_[frame1], Jc_[frame2]);
    myHistogram_[timeBin][0] += corrVal;

    for (unsigned int j = 0; j < outputTypes_.size(); j++) {
      corrVal = dot(typeJc_[frame1][j], typeJc_[frame2][j]);
      myHistogram_[timeBin][j + 1] += corrVal;
    }

    count_[timeBin]++;
  }

  void CurrentDensityAutoCorrFunc::postCorrelate() {
    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      for (unsigned int j = 0; j < outputTypes_.size() + 1; j++) {
        if (count_[i] > 0) {
          myHistogram_[i][j] /= count_[i];
        } else {
          myHistogram_[i][j] = 0.0;
        }
      }
    }
  }

  void CurrentDensityAutoCorrFunc::writeCorrelate() {
    ofstream ofs(outputFilename_.c_str());

    if (ofs.is_open()) {
      Revision r;

      ofs << "# " << getCorrFuncType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";
      ofs << "# units = Amps^2 m^-4\n";
      if (!labelString_.empty())
        ofs << "#time\t" << labelString_ << "\n";
      else
        ofs << "#time\tcorrVal\t(";

      for (unsigned int j = 0; j < outputTypes_.size(); j++) {
        ofs << outputTypes_[j]->getName() << "\t";
      }

      ofs << ")\n";

      for (unsigned int i = 0; i < nTimeBins_; ++i) {
        ofs << times_[i] - times_[0] << "\t";
        for (unsigned int j = 0; j < outputTypes_.size() + 1; j++) {
          ofs << myHistogram_[i][j] << '\t';
        }
        ofs << '\n';
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "CurrentDensityAutoCorrFunc::writeCorrelate Error: failed to "
               "open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }
}  // namespace OpenMD
