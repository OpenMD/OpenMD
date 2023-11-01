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

#include "applications/dynamicProps/RotAngleDisplacement.hpp"

#include <sstream>

#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {
  RotAngleDisplacement::RotAngleDisplacement(SimInfo* info,
                                             const std::string& filename,
                                             const std::string& sele1,
                                             const std::string& sele2) :
      MoleculeACF<Vector3d>(info, filename, sele1, sele2) {
    setCorrFuncType("Rotational Angle Displacement Function");
    setOutputName(getPrefix(dumpFilename_) + ".rotAngDisp");

    if (!uniqueSelections_) { seleMan2_ = seleMan1_; }

    rotMats_.resize(nTimeBins_);
    histogram_.resize(nTimeBins_);
    counts_.resize(nTimeBins_);
    std::fill(histogram_.begin(), histogram_.end(), 0.0);
    std::fill(counts_.begin(), counts_.end(), 0);
  }

  void RotAngleDisplacement::computeFrame(int frame) {
    MoleculeACF<Vector3d>::computeFrame(frame);
  }

  int RotAngleDisplacement::computeProperty1(int frame, Molecule* mol) {
    RotMat3x3d A = mol->getRigidBodyAt(0)->getA();
    rotMats_[frame].push_back(A);
    return rotMats_[frame].size() - 1;
  }

  Vector3d RotAngleDisplacement::calcCorrVal(int frame1, int frame2, int id1,
                                             int id2) {
    RotMat3x3d A1 = rotMats_[frame1][id1];
    RotMat3x3d A2 = rotMats_[frame2][id2];

    RotMat3x3d A21 = A1.transpose() * A2;
    Vector3d rpy   = A21.toRPY();

    return rpy;
  }

  void RotAngleDisplacement::correlateFrames(int frame1, int frame2,
                                             int timeBin) {
    std::vector<int> s1;
    std::vector<int> s2;

    std::vector<int>::iterator i1;
    std::vector<int>::iterator i2;

    Vector3d corrVal(0.0);

    s1 = sele1ToIndex_[frame1];

    if (uniqueSelections_)
      s2 = sele2ToIndex_[frame2];
    else
      s2 = sele1ToIndex_[frame2];

    for (i1 = s1.begin(), i2 = s2.begin(); i1 != s1.end() && i2 != s2.end();
         ++i1, ++i2) {
      // If the selections are dynamic, they might not have the
      // same objects in both frames, so we need to roll either of
      // the selections until we have the same object to
      // correlate.

      while (i1 != s1.end() && *i1 < *i2) {
        ++i1;
      }

      while (i2 != s2.end() && *i2 < *i1) {
        ++i2;
      }

      if (i1 == s1.end() || i2 == s2.end()) break;

      corrVal = calcCorrVal(frame1, frame2, i1 - s1.begin(), i2 - s2.begin());

      histogram_[timeBin] += corrVal;
      counts_[timeBin]++;
    }
  }

  void RotAngleDisplacement::postCorrelate() {
    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      if (counts_[i] > 0) { histogram_[i] /= counts_[i]; }
    }
  }

  void RotAngleDisplacement::validateSelection(SelectionManager&) {
    Molecule* mol;
    int i;
    for (mol = seleMan1_.beginSelectedMolecule(i); mol != NULL;
         mol = seleMan1_.nextSelectedMolecule(i)) {
      if (mol->getNRigidBodies() < 1) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "RotAngleDisplacement::validateSelection Error: "
                 "at least one selected molecule does not have a rigid body\n");
        painCave.isFatal = 1;
        simError();
      }
    }
    if (seleMan1_.getMoleculeSelectionCount() < 1) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RotAngleDisplacement::validateSelection Error: "
               "There needs to be at least one selected molecule.\n");
      painCave.isFatal = 1;
      simError();
    }
  }

  void RotAngleDisplacement::writeCorrelate() {
    std::string Anglefile = getOutputFileName();
    std::ofstream ofs1(Anglefile.c_str());

    if (ofs1.is_open()) {
      Revision r;

      ofs1 << "# " << getCorrFuncType() << "\n";
      ofs1 << "# OpenMD " << r.getFullRevision() << "\n";
      ofs1 << "# " << r.getBuildDate() << "\n";
      ofs1 << "# selection script1: \"" << selectionScript1_;
      ofs1 << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";

      ofs1 << "#time\troll\tpitch\tyaw\n";

      for (unsigned int i = 0; i < nTimeBins_; ++i) {
        ofs1 << times_[i] - times_[0];

        ofs1 << "\t" << histogram_[i][0];
        ofs1 << "\t" << histogram_[i][1];
        ofs1 << "\t" << histogram_[i][2] << "\n";
      }

    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "RotAngleDisplacement::writeCorrelate Error: failed to open %s\n",
          Anglefile.c_str());
      painCave.isFatal = 1;
      simError();
    }
    ofs1.close();
  }
}  // namespace OpenMD
