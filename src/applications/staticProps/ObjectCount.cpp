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

#include "applications/staticProps/ObjectCount.hpp"

#include <algorithm>
#include <functional>
#include <iomanip>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

namespace OpenMD {

  ObjectCount::ObjectCount(SimInfo* info, const std::string& filename,
                           const std::string& sele) :
      StaticAnalyser(info, filename, 1),
      selectionScript_(sele), seleMan_(info), evaluator_(info) {
    setOutputName(getPrefix(filename) + ".counts");

    evaluator_.loadScriptString(sele);

    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
  }

  void ObjectCount::process() {
    counts_.clear();
    counts_.resize(10, 0);
    DumpReader reader(info_, dumpFilename_);
    int nFrames             = reader.getNFrames();
    unsigned long int nsum  = 0;
    unsigned long int n2sum = 0;

    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      unsigned int count = seleMan_.getSelectionCount();

      if (counts_.size() < count + 1) { counts_.resize(count + 1, 0); }

      counts_[count]++;

      nsum += count;
      n2sum += count * count;
    }

    int nProcessed = nFrames / step_;

    nAvg  = RealType(nsum) / RealType(nProcessed);
    n2Avg = RealType(n2sum) / RealType(nProcessed);
    sDev  = sqrt(n2Avg - nAvg * nAvg);
    writeCounts();
  }

  void ObjectCount::writeCounts() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);
    if (ofs.is_open()) {
      ofs << "#counts\n";
      ofs << "#selection: (" << selectionScript_ << ")\n";
      ofs << "# <N> = " << std::fixed << std::setw(11) << std::setprecision(6)
          << nAvg << "\n";
      ofs << "# <N^2> = " << std::fixed << std::setw(11) << std::setprecision(6)
          << n2Avg << "\n";
      ofs << "# sqrt(<N^2> - <N>^2)  = " << sDev << "\n";
      ofs << "# N\tcounts[N]\n";
      for (unsigned int i = 0; i < counts_.size(); ++i) {
        ofs << i << "\t" << counts_[i] << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "ObjectCount: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    ofs.close();
  }

  MoleculeCount::MoleculeCount(SimInfo* info, const std::string& filename,
                               const std::string& sele) :
      ObjectCount(info, filename, sele) {
    setOutputName(getPrefix(filename) + ".mcounts");

    evaluator_.loadScriptString(sele);

    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }
  }

  void MoleculeCount::process() {
    counts_.clear();
    counts_.resize(10, 0);
    DumpReader reader(info_, dumpFilename_);
    int nFrames             = reader.getNFrames();
    unsigned long int nsum  = 0;
    unsigned long int n2sum = 0;

    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      unsigned int count = seleMan_.getMoleculeSelectionCount();

      if (counts_.size() < count + 1) { counts_.resize(count + 1, 0); }

      counts_[count]++;

      nsum += count;
      n2sum += count * count;
    }

    int nProcessed = nFrames / step_;

    nAvg  = RealType(nsum) / RealType(nProcessed);
    n2Avg = RealType(n2sum) / RealType(nProcessed);
    sDev  = sqrt(n2Avg - nAvg * nAvg);
    writeCounts();
  }

}  // namespace OpenMD
