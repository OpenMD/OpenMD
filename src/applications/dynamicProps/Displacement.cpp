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

#include "applications/dynamicProps/Displacement.hpp"

#include <sstream>

#include "utils/Revision.hpp"

namespace OpenMD {
  Displacement::Displacement(SimInfo* info, const std::string& filename,
                             const std::string& sele1,
                             const std::string& sele2) :
      ObjectACF<Vector3d>(info, filename, sele1, sele2) {
    setCorrFuncType("Displacement");
    setOutputName(getPrefix(dumpFilename_) + ".disp");

    positions_.resize(nFrames_);
  }

  DisplacementZ::DisplacementZ(SimInfo* info, const std::string& filename,
                               const std::string& sele1,
                               const std::string& sele2, int nZbins, int axis) :
      ObjectACF<Vector3d>(info, filename, sele1, sele2) {
    setCorrFuncType("Displacement binned by Z");
    setOutputName(getPrefix(dumpFilename_) + ".dispZ");

    positions_.resize(nFrames_);
    zBins_.resize(nFrames_);
    nZBins_ = nZbins;
    axis_   = axis;

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

    std::stringstream params;
    params << " nzbins = " << nZBins_;
    const std::string paramString = params.str();
    setParameterString(paramString);

    histograms_.resize(nTimeBins_);
    counts_.resize(nTimeBins_);
    for (unsigned int i = 0; i < nTimeBins_; i++) {
      histograms_[i].resize(nZBins_);
      counts_[i].resize(nZBins_);
      std::fill(histograms_[i].begin(), histograms_[i].end(), Vector3d(0.0));
      std::fill(counts_[i].begin(), counts_[i].end(), 0);
    }
  }

  int Displacement::computeProperty1(int frame, StuntDouble* sd) {
    positions_[frame].push_back(sd->getPos());
    return positions_[frame].size() - 1;
  }

  Vector3d Displacement::calcCorrVal(int frame1, int frame2, int id1, int id2) {
    return positions_[frame2][id2] - positions_[frame1][id1];
  }

  void DisplacementZ::computeFrame(int istep) {
    hmat_     = currentSnapshot_->getHmat();
    halfBoxZ_ = hmat_(axis_, axis_) / 2.0;

    StuntDouble* sd;

    int isd1, isd2;
    unsigned int index;

    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    if (uniqueSelections_ && evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }

    for (sd = seleMan1_.beginSelected(isd1); sd != NULL;
         sd = seleMan1_.nextSelected(isd1)) {
      index = computeProperty1(istep, sd);
      if (index == sele1ToIndex_[istep].size()) {
        sele1ToIndex_[istep].push_back(sd->getGlobalIndex());
      } else {
        sele1ToIndex_[istep].resize(index + 1);
        sele1ToIndex_[istep][index] = sd->getGlobalIndex();
      }
    }

    if (uniqueSelections_) {
      for (sd = seleMan2_.beginSelected(isd2); sd != NULL;
           sd = seleMan2_.nextSelected(isd2)) {
        index = computeProperty1(istep, sd);

        if (index == sele2ToIndex_[istep].size()) {
          sele2ToIndex_[istep].push_back(sd->getGlobalIndex());
        } else {
          sele2ToIndex_[istep].resize(index + 1);
          sele2ToIndex_[istep][index] = sd->getGlobalIndex();
        }
      }
    }
  }

  int DisplacementZ::computeProperty1(int frame, StuntDouble* sd) {
    Vector3d pos = sd->getPos();
    // we need the raw (not wrapped) positions for RMSD:
    positions_[frame].push_back(sd->getPos());

    if (info_->getSimParams()->getUsePeriodicBoundaryConditions()) {
      currentSnapshot_->wrapVector(pos);
    }
    int zBin = int(nZBins_ * (halfBoxZ_ + pos[axis_]) / hmat_(axis_, axis_));
    zBins_[frame].push_back(zBin);

    return positions_[frame].size() - 1;
  }

  void DisplacementZ::correlateFrames(int frame1, int frame2, int timeBin) {
    std::vector<int> s1;
    std::vector<int> s2;

    std::vector<int>::iterator i1;
    std::vector<int>::iterator i2;

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

      calcCorrValImpl(frame1, frame2, i1 - s1.begin(), i2 - s2.begin(),
                      timeBin);
    }
  }

  Vector3d DisplacementZ::calcCorrValImpl(int frame1, int frame2, int id1,
                                          int id2, int timeBin) {
    int zBin1 = zBins_[frame1][id1];
    int zBin2 = zBins_[frame2][id2];

    if (zBin1 == zBin2) {
      Vector3d diff = positions_[frame2][id2] - positions_[frame1][id1];
      histograms_[timeBin][zBin1] += diff;
      counts_[timeBin][zBin1]++;
    }
    return Vector3d(0.0);
  }

  void DisplacementZ::postCorrelate() {
    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      for (unsigned int j = 0; j < nZBins_; ++j) {
        if (counts_[i][j] > 0) {
          histograms_[i][j] /= counts_[i][j];
        } else {
          histograms_[i][j] = Vector3d(0.0);
        }
      }
    }
  }
  void DisplacementZ::writeCorrelate() {
    std::ofstream ofs(getOutputFileName().c_str());

    if (ofs.is_open()) {
      Revision r;

      ofs << "# " << getCorrFuncType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      ofs << "# privilegedAxis computed as " << axisLabel_ << " axis \n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      ofs << "#time\tcorrVal\n";

      for (unsigned int i = 0; i < nTimeBins_; ++i) {
        ofs << times_[i] - times_[0] << "\t";

        for (unsigned int j = 0; j < nZBins_; ++j) {
          for (int k = 0; k < 3; k++) {
            ofs << histograms_[i][j](k) << '\t';
          }
        }

        ofs << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DisplacementZ::writeCorrelate Error: fail to open %s\n",
               getOutputFileName().c_str());
      painCave.isFatal = 1;
      simError();
    }
    ofs.close();
  }
}  // namespace OpenMD
