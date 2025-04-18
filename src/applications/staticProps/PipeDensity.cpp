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

/* Calculates PipeDensity, rho(axis2,axis3) in the box */

#include "applications/staticProps/PipeDensity.hpp"

#include <algorithm>
#include <fstream>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

namespace OpenMD {

  PipeDensity::PipeDensity(SimInfo* info, const std::string& filename,
                           const std::string& sele, int nbins, int nbins2,
                           int axis) :
      StaticAnalyser(info, filename, nbins2),
      selectionScript_(sele), evaluator_(info), seleMan_(info), nBins2_(nbins),
      axis_(axis) {
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // fixed number of bins

    sliceSDLists_.resize(nBins2_);
    density_.resize(nBins2_);
    for (unsigned int i = 0; i < nBins2_; ++i) {
      sliceSDLists_[i].resize(nBins_);
      density_[i].resize(nBins_);
    }

    // Compute complementary axes to the privileged axis
    axis1_ = (axis_ + 1) % 3;
    axis2_ = (axis_ + 2) % 3;

    // Set the axis labels for the non-privileged axes
    switch (axis_) {
    case 0:
      axisLabel1_ = "y";
      axisLabel2_ = "z";
      break;
    case 1:
      axisLabel1_ = "z";
      axisLabel2_ = "x";
      break;
    case 2:
    default:
      axisLabel1_ = "x";
      axisLabel2_ = "y";
      break;
    }

    setOutputName(getPrefix(filename) + ".PipeDensity");
  }

  void PipeDensity::process() {
    StuntDouble* sd;
    int ii;

    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames / step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for (unsigned int i = 0; i < nBins2_; i++) {
        for (unsigned int j = 0; j < nBins_; j++) {
          sliceSDLists_[i][j].clear();
        }
      }

      RealType sliceVolume = currentSnapshot_->getVolume() / (nBins2_ * nBins_);
      Mat3x3d hmat         = currentSnapshot_->getHmat();

      RealType halfBox1_ = hmat(axis1_, axis1_) / 2.0;
      RealType halfBox2_ = hmat(axis2_, axis2_) / 2.0;

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      // wrap the stuntdoubles into a cell
      for (sd = seleMan_.beginSelected(ii); sd != NULL;
           sd = seleMan_.nextSelected(ii)) {
        Vector3d pos = sd->getPos();
        if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(pos);
        sd->setPos(pos);
      }

      // determine which atom belongs to which slice
      for (sd = seleMan_.beginSelected(ii); sd != NULL;
           sd = seleMan_.nextSelected(ii)) {
        Vector3d pos = sd->getPos();
        // shift molecules by half a box to have bins start at 0
        int binNo1 =
            int(nBins2_ * (halfBox1_ + pos[axis1_]) / hmat(axis1_, axis1_));
        int binNo2 =
            int(nBins_ * (halfBox2_ + pos[axis2_]) / hmat(axis2_, axis2_));
        sliceSDLists_[binNo1][binNo2].push_back(sd);
      }

      // loop over the slices to calculate the densities
      for (unsigned int i = 0; i < nBins2_; i++) {
        for (unsigned int j = 0; j < nBins_; j++) {
          RealType totalMass = 0;
          for (unsigned int k = 0; k < sliceSDLists_[i][j].size(); ++k) {
            totalMass += sliceSDLists_[i][j][k]->getMass();
          }
          density_[i][j] += totalMass / sliceVolume;
        }
      }
    }

    writeDensity();
  }

  void PipeDensity::writeDensity() {
    std::ofstream rdfStream(outputFilename_.c_str());

    if (rdfStream.is_open()) {
      rdfStream << "#PipeDensity\n";
      rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#density (" << axisLabel1_ << "," << axisLabel2_ << ")\n";
      for (unsigned int i = 0; i < density_.size(); ++i) {
        for (unsigned int j = 0; j < density_[i].size(); ++j) {
          rdfStream << Constants::densityConvert * density_[i][j] / nProcessed_;
          rdfStream << "\t";
        }
        rdfStream << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "PipeDensity: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }
}  // namespace OpenMD
