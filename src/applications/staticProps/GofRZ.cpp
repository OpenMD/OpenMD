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

#include "applications/staticProps/GofRZ.hpp"

#include <algorithm>
#include <fstream>

#include "utils/simError.h"

namespace OpenMD {

  GofRZ::GofRZ(SimInfo* info, const std::string& filename,
               const std::string& sele1, const std::string& sele2, RealType len,
               RealType zlen, int nrbins, int nZBins, int axis) :
      RadialDistrFunc(info, filename, sele1, sele2, nrbins),
      len_(len), zLen_(zlen), nZBins_(nZBins), axis_(axis) {
    setOutputName(getPrefix(filename) + ".gofrz");

    deltaR_ = len_ / (double)nBins_;
    deltaZ_ = zLen_ / (double)nZBins_;

    histogram_.resize(nBins_);
    avgGofr_.resize(nBins_);
    for (unsigned int i = 0; i < nBins_; ++i) {
      histogram_[i].resize(nZBins_);
      avgGofr_[i].resize(nZBins_);
    }

    // Compute complementary axes to the privileged axis
    xaxis_ = (axis_ + 1) % 3;
    yaxis_ = (axis_ + 2) % 3;

    // Set the axis label for the privileged axis
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
  }

  void GofRZ::preProcess() {
    for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
      std::fill(avgGofr_[i].begin(), avgGofr_[i].end(), 0);
    }
  }

  void GofRZ::initializeHistogram() {
    npairs_ = 0;
    for (unsigned int i = 0; i < histogram_.size(); ++i) {
      std::fill(histogram_[i].begin(), histogram_[i].end(), 0);
    }
  }

  void GofRZ::processHistogram() {
    int nPairs = getNPairs();
    RealType volume =
        info_->getSnapshotManager()->getCurrentSnapshot()->getVolume();
    RealType pairDensity = nPairs / volume * 2.0;

    for (unsigned int i = 0; i < histogram_.size(); ++i) {
      RealType rLower = i * deltaR_;
      RealType rUpper = rLower + deltaR_;
      RealType volSlice =
          Constants::PI * deltaZ_ * ((rUpper * rUpper) - (rLower * rLower));
      RealType nIdeal = volSlice * pairDensity;

      for (unsigned int j = 0; j < histogram_[i].size(); ++j) {
        avgGofr_[i][j] += histogram_[i][j] / nIdeal;
      }
    }
  }

  void GofRZ::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {
    if (sd1 == sd2) { return; }
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12  = pos2 - pos1;
    if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(r12);

    RealType distance = sqrt(pow(r12[xaxis_], 2) + pow(r12[yaxis_], 2));

    int whichRBin = int(distance / deltaR_);

    if (distance <= len_) {
      RealType Z = fabs(r12[axis_]);

      if (Z <= zLen_) {
        int whichZBin = int(Z / deltaZ_);

        ++histogram_[whichRBin][whichZBin];
        ++npairs_;
      }
    }
  }

  void GofRZ::writeRdf() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#radial distribution function\n";
      rdfStream << "#selection1: (" << selectionScript1_ << ")\t";
      rdfStream << "selection2: (" << selectionScript2_ << ")\n";
      rdfStream << "#nBins = " << nBins_ << "\t maxLen = " << len_
                << "deltaR = " << deltaR_ << "\n";
      rdfStream << "#n" << axisLabel_ << "Bins =" << nZBins_ << "\t delta"
                << axisLabel_ << " = " << deltaZ_ << "\n";
      for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
        // RealType r = deltaR_ * (i + 0.5);

        for (unsigned int j = 0; j < avgGofr_[i].size(); ++j) {
          // RealType z = deltaZ_ * (j + 0.5);
          rdfStream << avgGofr_[i][j] / nProcessed_ << "\t";
        }

        rdfStream << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "GofRZ: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }

}  // namespace OpenMD
