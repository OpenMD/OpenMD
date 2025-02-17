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

#include "applications/sequentialProps/GCNSeq.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>

#include "io/DumpReader.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {

  GCNSeq::GCNSeq(SimInfo* info, const std::string& filename,
                 const std::string& sele1, const std::string& sele2,
                 RealType rCut, int bins) :
      SequentialAnalyzer(info, filename, sele1, sele2),
      rCut_(rCut), bins_(bins) {
    setSequenceType("Generalized Coordination Number Distribution");
    setOutputName(getPrefix(filename) + ".gcnSeq");

    nnMax_           = 12;
    RealType binMax_ = nnMax_ * 1.5;
    delta_           = binMax_ / bins_;
    usePBC_          = info->getSimParams()->getUsePeriodicBoundaryConditions();

    std::stringstream params;
    params << " rcut = " << rCut_ << ", nbins = " << bins_
           << ", max neighbors = " << nnMax_;
    const std::string paramString = params.str();
    setParameterString(paramString);
  }

  void GCNSeq::doFrame(int istep) {
    SelectionManager common(info_);

    std::vector<std::vector<int>> listNN;
    std::vector<int> globalToLocal;

    StuntDouble* sd1;
    StuntDouble* sd2;

    int iterator1;
    int iterator2;
    unsigned int mapIndex1(0);
    unsigned int mapIndex2(0);
    unsigned int tempIndex(0);
    unsigned int whichBin(0);
    RealType gcn(0.0);
    Vector3d pos1;
    Vector3d pos2;
    Vector3d diff;
    RealType distance;

    // First have to calculate lists of nearest neighbors (listNN_):

    selectionCount1_ = seleMan1_.getSelectionCount();
    selectionCount2_ = seleMan2_.getSelectionCount();

    // We need a common selection set:
    common          = seleMan1_ | seleMan2_;
    int commonCount = common.getSelectionCount();

    globalToLocal.clear();
    globalToLocal.resize(
        info_->getNGlobalAtoms() + info_->getNGlobalRigidBodies(), -1);
    for (unsigned int i = 0; i < listNN.size(); i++)
      listNN.at(i).clear();
    listNN.clear();
    listNN.resize(commonCount);
    std::vector<RealType> histo;
    histo.resize(bins_, 0.0);

    mapIndex1 = 0;
    for (sd1 = common.beginSelected(iterator1); sd1 != NULL;
         sd1 = common.nextSelected(iterator1)) {
      globalToLocal.at(sd1->getGlobalIndex()) = mapIndex1;

      pos1 = sd1->getPos();

      mapIndex2 = 0;
      for (sd2 = common.beginSelected(iterator2); sd2 != NULL;
           sd2 = common.nextSelected(iterator2)) {
        if (mapIndex1 < mapIndex2) {
          pos2 = sd2->getPos();
          diff = pos2 - pos1;
          if (usePBC_) currentSnapshot_->wrapVector(diff);
          distance = diff.length();
          if (distance < rCut_) {
            listNN.at(mapIndex1).push_back(mapIndex2);
            listNN.at(mapIndex2).push_back(mapIndex1);
          }
        }
        mapIndex2++;
      }
      mapIndex1++;
    }

    // Fill up the histogram with gcn values
    for (sd1 = seleMan1_.beginSelected(iterator1); sd1 != NULL;
         sd1 = seleMan1_.nextSelected(iterator1)) {
      mapIndex1 = globalToLocal.at(sd1->getGlobalIndex());
      gcn       = 0.0;
      for (unsigned int i = 0; i < listNN.at(mapIndex1).size(); i++) {
        // tempIndex is the index of one of i's nearest neighbors
        tempIndex = listNN.at(mapIndex1).at(i);
        gcn += listNN.at(tempIndex).size();
      }

      gcn      = gcn / nnMax_;
      whichBin = int(gcn / delta_);
      if (whichBin < histo.size()) {
        histo[whichBin] += 1;
      } else {
        cerr << "In frame " << istep << ", object " << sd1->getGlobalIndex()
             << " has GCN value = " << gcn << "\n";
      }
    }

    for (unsigned int n = 0; n < histo.size(); n++) {
      if (selectionCount1_ > 0)
        histo[n] /= RealType(selectionCount1_);
      else
        histo[n] = 0.0;
    }

    count_.push_back(selectionCount1_);
    histogram_.push_back(histo);
  }

  void GCNSeq::writeSequence() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);

    if (ofs.is_open()) {
      Revision r;
      RealType binValue(0.0);

      ofs << "# " << getSequenceType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      ofs << "#time\tvalue\n";

      for (unsigned int i = 0; i < times_.size(); ++i) {
        ofs << "#Frame " << i << "\n";
        ofs << "#Selection 1 Count: " << count_[i] << "\n";

        for (unsigned int n = 0; n < histogram_[i].size(); n++) {
          binValue = n * delta_;
          ofs << binValue << "\t" << histogram_[i][n] << "\n";
        }
      }
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "GCN::writeSequence Error: failed to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }
}  // namespace OpenMD
