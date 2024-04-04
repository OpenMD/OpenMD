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

#include "applications/staticProps/TranslationalOrderParamZ.hpp"

#include <algorithm>
#include <fstream>
#include <vector>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {
  TranslationalOrderParamZ::TranslationalOrderParamZ(
      SimInfo* info, const std::string& filename, const std::string& sele1,
      const std::string& sele2, double rCut, int nrbins, int nzbins,
      RealType len, RealType zlen, int axis) :
      RadialDistrFunc(info, filename, sele1, sele2, nrbins),
      rCut_(rCut), nZBins_(nzbins), len_(len), zLen_(zlen), axis_(axis) {
    setOutputName(getPrefix(filename) + ".Tz");

    deltaR_ = len_ / (double)nBins_;
    deltaZ_ = zLen_ / (double)nZBins_;

    histogram_.resize(nBins_);
    avgGofr_.resize(nBins_);
    for (unsigned int i = 0; i < nBins_; ++i) {
      histogram_[i].resize(nZBins_);
      avgGofr_[i].resize(nZBins_);
    }
    Tz_.resize(nZBins_);

    // Set up cutoff radius:
    rCut_ = rCut;

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

  void TranslationalOrderParamZ::preProcess() {
    for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
      std::fill(avgGofr_[i].begin(), avgGofr_[i].end(), 0);
    }
    std::fill(Tz_.begin(), Tz_.end(), 0);
  }

  void TranslationalOrderParamZ::initializeHistogram() {
    for (unsigned int i = 0; i < histogram_.size(); ++i) {
      std::fill(histogram_[i].begin(), histogram_[i].end(), 0);
    }
    Mat3x3d hmat = currentSnapshot_->getHmat();
    zBox_.push_back(hmat(axis_, axis_));
  }

  void TranslationalOrderParamZ::collectHistogram(StuntDouble* sd1,
                                                  StuntDouble* sd2) {
    if (sd1 == sd2) { return; }

    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();
    RealType boxZ = zBox_.back();

    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12  = pos2 - pos1;
    if (usePeriodicBoundaryConditions_) {
      currentSnapshot_->wrapVector(r12);
      currentSnapshot_->wrapVector(pos1);
      currentSnapshot_->wrapVector(pos2);
    }

    RealType distance = r12.length();

    if (distance < len_) {
      int whichBin = int(distance / deltaR_);
      int zBin1    = int(nZBins_ * (0.5 * boxZ + pos1[axis_]) / boxZ);
      int zBin2    = int(nZBins_ * (0.5 * boxZ + pos2[axis_]) / boxZ);

      histogram_[whichBin][zBin1] += 1;
      histogram_[whichBin][zBin2] += 1;
    }
  }

  void TranslationalOrderParamZ::processHistogram() {
    int nPairs = getNPairs();
    RealType volume =
        info_->getSnapshotManager()->getCurrentSnapshot()->getVolume();
    RealType pairDensity = 2 * nPairs / volume;

    for (unsigned int i = 0; i < histogram_.size(); ++i) {
      RealType rLower = i * deltaR_;
      RealType rUpper = rLower + deltaR_;
      RealType volSlice =
          4.0 * Constants::PI * (pow(rUpper, 3) - pow(rLower, 3)) / 3.0;
      RealType nIdeal = volSlice * pairDensity / nZBins_;

      for (unsigned int j = 0; j < histogram_[i].size(); ++j) {
        avgGofr_[i][j] += histogram_[i][j] / nIdeal;
      }
    }
  }

  void TranslationalOrderParamZ::postProcess() {
    for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
      for (unsigned int j = 0; j < avgGofr_[i].size(); ++j) {
        avgGofr_[i][j] /= nProcessed_;
      }
    }

    for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
      RealType rLower = i * deltaR_;
      RealType rUpper = rLower + deltaR_;

      for (unsigned int j = 0; j < avgGofr_[i].size(); ++j) {
        if (rUpper < rCut_) Tz_[j] += std::fabs(avgGofr_[i][j] - 1.0) * deltaR_;
      }
    }

    // normalize by cutoff radius
    for (unsigned int j = 0; j < Tz_.size(); ++j) {
      Tz_[j] /= rCut_;
    }
  }

  void TranslationalOrderParamZ::writeRdf() {
    // compute average box length:

    RealType zSum = 0.0;
    for (std::vector<RealType>::iterator j = zBox_.begin(); j != zBox_.end();
         ++j) {
      zSum += *j;
    }
    RealType zAve = zSum / zBox_.size();

    std::ofstream tZstream(outputFilename_.c_str());
    if (tZstream.is_open()) {
      tZstream << "#Translational Order Parameters (" << axisLabel_ << ")\n";

      tZstream << "#nFrames:\t" << zBox_.size() << "\n";
      tZstream << "#selection 1: (" << selectionScript1_ << ")\n";
      tZstream << "#selection 2: (" << selectionScript2_ << ")\n";
      tZstream << "#" << axisLabel_ << "\tT\n";

      // for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
      //   RealType r = i * deltaR_;
      //   tZstream << r << "\t" ;
      //     for (unsigned int j = 0; j < avgGofr_[i].size(); ++j) {
      //       tZstream << "\t" << avgGofr_[i][j];
      //     }
      //   tZstream << "\n";
      // }

      for (unsigned int i = 0; i < Tz_.size(); ++i) {
        RealType z = zAve * (i + 0.5) / Tz_.size();
        tZstream << z << "\t" << Tz_[i] << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "TranslationalOrderParamZ: unable to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    tZstream.close();
  }
}  // namespace OpenMD
