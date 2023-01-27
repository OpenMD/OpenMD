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

#include "applications/staticProps/KirkwoodBuff.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {

  KirkwoodBuff::KirkwoodBuff(SimInfo* info, const std::string& filename,
                             const std::string& sele1, const std::string& sele2,
                             RealType len, unsigned int nrbins) :
      MultiComponentRDF {info, filename, sele1, sele2, nrbins},
      len_ {len}, meanVol_ {0.0} {
    setAnalysisType("Kirkwood-Buff Integrals");
    setOutputName(getPrefix(filename) + ".kirkwood-buff");

    // Bins are set to the inner radius of a spherical shell.
    deltaR_ = len_ / (nBins_ - 1);

    histograms_.resize(MaxPairs);
    gofrs_.resize(MaxPairs);
    gCorr_.resize(MaxPairs);
    G_.resize(MaxPairs);

    for (std::size_t i {}; i < MaxPairs; ++i) {
      histograms_[i].resize(nBins_);
      gofrs_[i].resize(nBins_);
      gCorr_[i].resize(nBins_);
      G_[i].resize(nBins_);
    }

    std::stringstream params;
    params << " len = " << len_ << ", nrbins = " << nBins_;
    const std::string paramString = params.str();
    setParameterString(paramString);
  }

  void KirkwoodBuff::initializeHistograms() {
    for (auto& pair : histograms_)
      std::fill(pair.begin(), pair.end(), 0);
  }

  void KirkwoodBuff::collectHistograms(StuntDouble* sd1, StuntDouble* sd2,
                                       int pairIndex) {
    if (sd1 == sd2) { return; }

    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12  = pos2 - pos1;
    if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(r12);

    RealType distance = r12.length();

    // Bins are set to the inner radius of a spherical shell.
    if (distance < (len_ + deltaR_)) {
      int whichBin = static_cast<int>(distance / deltaR_);
      histograms_[pairIndex][whichBin] += 1;
    }
  }

  void KirkwoodBuff::processHistograms() {
    std::vector<int> nPairs = getNPairs();
    RealType volume =
        info_->getSnapshotManager()->getCurrentSnapshot()->getVolume();

    meanVol_ += volume;

    for (std::size_t i {}; i < histograms_.size(); ++i) {
      for (std::size_t j {}; j < histograms_[i].size(); ++j) {
        RealType rLower   = j * deltaR_;
        RealType rUpper   = rLower + deltaR_;
        RealType volSlice = 4.0 * Constants::PI *
                            (std::pow(rUpper, 3) - std::pow(rLower, 3)) / 3.0;
        RealType pairDensity = nPairs[i] / volume;
        RealType nIdeal      = volSlice * pairDensity;

        gofrs_[i][j] += histograms_[i][j] / nIdeal;
      }
    }
  }

  void KirkwoodBuff::postProcess() {
    for (std::size_t i {}; i < gofrs_.size(); ++i) {
      for (std::size_t j {}; j < gofrs_[i].size(); ++j) {
        gofrs_[i][j] /= nProcessed_;
      }
    }  // g(r) is the uncorrected g(r)
    meanVol_ /= nProcessed_;

    std::vector<int> Ns(MaxPairs);
    Ns[OneOne] = getNSelected1();
    Ns[OneTwo] = getNSelected2();
    Ns[TwoTwo] = getNSelected2();

    std::vector<int> kd(MaxPairs);
    kd[OneOne] = 1;
    kd[OneTwo] = 0;
    kd[TwoTwo] = 1;

    std::vector<RealType> rho(MaxPairs);
    rho[OneOne] = Ns[OneOne] / meanVol_;
    rho[OneTwo] = Ns[OneTwo] / meanVol_;
    rho[TwoTwo] = Ns[TwoTwo] / meanVol_;

    std::vector<std::vector<RealType>> deltaN;
    deltaN.resize(MaxPairs);
    for (auto& elem : deltaN)
      elem.resize(nBins_);

    for (std::size_t i {}; i < deltaN.size(); ++i) {
      deltaN[i][0] = 0.0;
      gCorr_[i][0] = 0.0;
      G_[i][0]     = 0.0;

      for (std::size_t j {1}; j < deltaN[i].size(); ++j) {
        RealType r = deltaR_ * j;
        RealType V = 4.0 * Constants::PI * r * r * r / 3.0;
        RealType x = r / len_;
        RealType w =
            4.0 * Constants::PI * r * r * (1 - 3 * x / 2 + std::pow(x, 3) / 2);

        deltaN[i][j] += deltaN[i][j - 1] + 4.0 * Constants::PI * r * r *
                                               rho[i] * (gofrs_[i][j] - 1) *
                                               deltaR_;
        gCorr_[i][j] = gofrs_[i][j] * (Ns[i] * (1 - V / meanVol_)) /
                       (Ns[i] * (1 - V / meanVol_) - deltaN[i][j] - kd[i]);

        G_[i][j] += G_[i][j - 1] + (gCorr_[i][j] - 1) * w * deltaR_;
      }
    }
  }

  void KirkwoodBuff::writeRdf() {
    std::ofstream ofs(outputFilename_.c_str());
    if (ofs.is_open()) {
      Revision r;
      ofs << "# " << getAnalysisType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      std::vector<std::string> labels {"g11",  "g12", "g22", "gC11", "gC12",
                                       "gC22", "G11", "G12", "G22"};

      ofs << "# r";

      for (const auto& label : labels)
        ofs << "\t" << std::setw(15) << label;

      ofs << '\n';

      for (unsigned int j = 0; j < nBins_; ++j) {
        RealType r = deltaR_ * j;
        ofs << r;
        for (unsigned int i = 0; i < MaxPairs; ++i) {
          ofs << "\t" << std::setw(15) << gofrs_[i][j];
        }
        for (unsigned int i = 0; i < MaxPairs; ++i) {
          ofs << "\t" << std::setw(15) << gCorr_[i][j];
        }
        for (unsigned int i = 0; i < MaxPairs; ++i) {
          ofs << "\t" << std::setw(15) << G_[i][j];
        }
        ofs << "\n";
      }
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "KirkwoodBuff: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    ofs.close();
  }
}  // namespace OpenMD
