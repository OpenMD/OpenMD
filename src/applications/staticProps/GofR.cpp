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

#include "applications/staticProps/GofR.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>

#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {

  GofR::GofR(SimInfo* info, const std::string& filename,
             const std::string& sele1, const std::string& sele2, RealType len,
             int nrbins) :
      RadialDistrFunc(info, filename, sele1, sele2, nrbins),
      len_(len) {
    setAnalysisType("Radial Distribution Function");
    setOutputName(getPrefix(filename) + ".gofr");

    deltaR_ = len_ / nBins_;

    histogram_.resize(nBins_);
    avgGofr_.resize(nBins_);
    sumGofr1_.resize(nBins_);
    sumGofr2_.resize(nBins_);
    std::stringstream params;
    params << " len = " << len_ << ", nrbins = " << nBins_;
    const std::string paramString = params.str();
    setParameterString(paramString);
  }

  void GofR::preProcess() {
    std::fill(avgGofr_.begin(), avgGofr_.end(), 0.0);
    std::fill(sumGofr1_.begin(), sumGofr1_.end(), 0.0);
    std::fill(sumGofr2_.begin(), sumGofr2_.end(), 0.0);
  }

  void GofR::initializeHistogram() {
    std::fill(histogram_.begin(), histogram_.end(), 0);
  }

  void GofR::processHistogram() {
    int nPairs = getNPairs();
    RealType volume =
        info_->getSnapshotManager()->getCurrentSnapshot()->getVolume();
    RealType pairDensity  = nPairs / volume * 2.0;
    RealType pairConstant = (4.0 * Constants::PI * pairDensity) / 3.0;

    for (unsigned int i = 0; i < histogram_.size(); ++i) {
      RealType rLower = i * deltaR_;
      RealType rUpper = rLower + deltaR_;
      RealType volSlice =
          (rUpper * rUpper * rUpper) - (rLower * rLower * rLower);
      RealType nIdeal = volSlice * pairConstant;

      avgGofr_[i] += histogram_[i] / nIdeal;
    }
  }

  void GofR::postProcess() {
    int nSelected1 = getNSelected1();
    int nSelected2 = getNSelected2();
    RealType volume =
        info_->getSnapshotManager()->getCurrentSnapshot()->getVolume();
    RealType constant = (4.0 * Constants::PI) / (3.0 * volume);

    RealType sum = 0.0;
    for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
      RealType rLower = i * deltaR_;
      RealType rUpper = rLower + deltaR_;

      sum += avgGofr_[i] * constant * (pow(rUpper, 3) - pow(rLower, 3));
      sumGofr1_[i] = nSelected1 * sum;
      sumGofr2_[i] = nSelected2 * sum;
    }
  }

  void GofR::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {
    if (sd1 == sd2) { return; }

    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12  = pos2 - pos1;
    if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(r12);

    RealType distance = r12.length();

    if (distance < len_) {
      int whichBin = int(distance / deltaR_);
      histogram_[whichBin] += 2;
    }
  }

  void GofR::writeRdf() {
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

      ofs << "#r\tcorrValue\n";
      for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
        RealType r = deltaR_ * (i + 0.5);
        ofs << r << "\t" << avgGofr_[i] / nProcessed_  << "\n";
      }
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "GofR: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    ofs.close();
  }
}  // namespace OpenMD
