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

/* Calculates Velocity Map, V(axis2,axis3) in the box */

#include "applications/staticProps/VelocityZ.hpp"

#include <algorithm>
#include <fstream>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

namespace OpenMD {

  VelocityZ::VelocityZ(SimInfo* info, const std::string& filename,
                       const std::string& sele, int nbins1, int nbins2,
                       int axis1, int axis2) :
      StaticAnalyser(info, filename, nbins1),
      selectionScript_(sele), evaluator_(info), seleMan_(info), nBins2_(nbins2),
      axis1_(axis1), axis2_(axis2) {
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    /* fixed number of bins:
     nBins_ is along the primary axis, nBins2_ along the secondary axis
     there are nBins2_ perpendicular to each slab of nBins_
  */

    sliceSDLists_.resize(nBins2_);
    velocity_.resize(nBins2_);
    for (unsigned int i = 0; i < nBins2_; ++i) {
      sliceSDLists_[i].resize(nBins_);
      velocity_[i].resize(nBins_);
    }

    // Compute complementary axis to the two privileged axis
    axis3_ = (3 - axis1_ - axis2_);

    // Set the axis labels for the non-privileged axes
    switch (axis1_) {
    case 0:
      axisLabel1_ = "x";
      if (axis2_ == 1)
        axisLabel2_ = "y";
      else if (axis2_ == 2)
        axisLabel2_ = "z";
      break;
    case 1:
      axisLabel1_ = "y";
      if (axis2_ == 0)
        axisLabel2_ = "x";
      else if (axis2_ == 2)
        axisLabel2_ = "z";
      break;
    case 2:
    default:
      axisLabel1_ = "z";
      if (axis2_ == 0)
        axisLabel2_ = "x";
      else if (axis2_ == 1)
        axisLabel2_ = "y";
      break;
    }

    setOutputName(getPrefix(filename) + ".VelocityZ");
  }

  void VelocityZ::process() {
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

      Mat3x3d hmat = currentSnapshot_->getHmat();

      zBox_.push_back(hmat(axis2_, axis2_));

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
            int(nBins_ * (halfBox1_ + pos[axis1_]) / hmat(axis1_, axis1_));
        int binNo2 =
            int(nBins2_ * (halfBox2_ + pos[axis2_]) / hmat(axis2_, axis2_));
        sliceSDLists_[binNo2][binNo1].push_back(sd);
      }

      // loop over the slices to calculate the velocities
      for (unsigned int i = 0; i < nBins2_; i++) {
        for (unsigned int j = 0; j < nBins_; j++) {
          RealType totalVelocity = 0;
          for (unsigned int k = 0; k < sliceSDLists_[i][j].size(); ++k) {
            totalVelocity += sliceSDLists_[i][j][k]->getVel()[axis3_];
          }

          if (sliceSDLists_[i][j].size() > 0)
            velocity_[i][j] += totalVelocity / sliceSDLists_[i][j].size();
        }
      }
    }

    writeVelocity();
  }

  void VelocityZ::writeVelocity() {
    // compute average box length:

    RealType zSum = 0.0;
    for (std::vector<RealType>::iterator j = zBox_.begin(); j != zBox_.end();
         ++j) {
      zSum += *j;
    }
    RealType zAve = zSum / zBox_.size();

    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#VelocityZ\n";
      rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#velocity (" << axisLabel1_ << "," << axisLabel2_ << ")\n";

      for (unsigned int i = 0; i < velocity_.size(); ++i) {
        RealType z = zAve * (i + 0.5) / velocity_.size();
        rdfStream << z << "\t";
        for (unsigned int j = 0; j < velocity_[i].size(); ++j) {
          rdfStream << velocity_[i][j] / nProcessed_;
          rdfStream << "\t";
        }
        rdfStream << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "VelocityZ: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }
}  // namespace OpenMD
