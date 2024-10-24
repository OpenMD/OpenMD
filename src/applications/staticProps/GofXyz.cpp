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

#include "applications/staticProps/GofXyz.hpp"

#include <algorithm>
#include <fstream>

#include "primitives/Molecule.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  GofXyz::GofXyz(SimInfo* info, const std::string& filename,
                 const std::string& sele1, const std::string& sele2,
                 const std::string& sele3, RealType len, int nrbins) :
      RadialDistrFunc(info, filename, sele1, sele2, nrbins),
      len_(len), halfLen_(len / 2), evaluator3_(info), seleMan3_(info) {
    setOutputName(getPrefix(filename) + ".gxyz");

    evaluator3_.loadScriptString(sele3);
    if (!evaluator3_.isDynamic()) {
      seleMan3_.setSelectionSet(evaluator3_.evaluate());
    }

    deltaR_ = len_ / nBins_;

    histogram_.resize(nBins_);
    for (unsigned int i = 0; i < nBins_; ++i) {
      histogram_[i].resize(nBins_);
      for (unsigned int j = 0; j < nBins_; ++j) {
        histogram_[i][j].resize(nBins_);
      }
    }
  }

  void GofXyz::preProcess() {
    for (unsigned int i = 0; i < nBins_; ++i) {
      histogram_[i].resize(nBins_);
      for (unsigned int j = 0; j < nBins_; ++j) {
        std::fill(histogram_[i][j].begin(), histogram_[i][j].end(), 0);
      }
    }
  }

  void GofXyz::initializeHistogram() {
    // Calculate the center of mass of the molecule of selected
    // StuntDouble in selection1

    if (!evaluator3_.isDynamic()) {
      seleMan3_.setSelectionSet(evaluator3_.evaluate());
    }

    assert(seleMan1_.getSelectionCount() == seleMan3_.getSelectionCount());
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    // The Dipole direction of selection3 and position of selection3 will
    // be used to determine the y-z plane
    // v1 = s3 -s1,
    // z = origin.dipole
    // x = v1 X z
    // y = z X x
    rotMats_.clear();

    int i;
    int j;
    StuntDouble* sd1;
    StuntDouble* sd3;

    for (sd1 = seleMan1_.beginSelected(i), sd3 = seleMan3_.beginSelected(j);
         sd1 != NULL || sd3 != NULL;
         sd1 = seleMan1_.nextSelected(i), sd3 = seleMan3_.nextSelected(j)) {
      Vector3d r3 = sd3->getPos();
      Vector3d r1 = sd1->getPos();
      Vector3d v1 = r3 - r1;
      if (usePeriodicBoundaryConditions_)
        info_->getSnapshotManager()->getCurrentSnapshot()->wrapVector(v1);

      AtomType* atype1     = static_cast<Atom*>(sd1)->getAtomType();
      MultipoleAdapter ma1 = MultipoleAdapter(atype1);

      Vector3d zaxis;
      if (ma1.isDipole())
        zaxis = sd1->getDipole();
      else
        zaxis = sd1->getA().transpose() * V3Z;

      Vector3d xaxis = cross(v1, zaxis);
      Vector3d yaxis = cross(zaxis, xaxis);

      xaxis.normalize();
      yaxis.normalize();
      zaxis.normalize();

      RotMat3x3d rotMat;
      rotMat.setRow(0, xaxis);
      rotMat.setRow(1, yaxis);
      rotMat.setRow(2, zaxis);

      rotMats_.insert(
          std::map<int, RotMat3x3d>::value_type(sd1->getGlobalIndex(), rotMat));
    }
  }

  void GofXyz::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12  = pos2 - pos1;
    if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(r12);

    std::map<int, RotMat3x3d>::iterator i =
        rotMats_.find(sd1->getGlobalIndex());
    assert(i != rotMats_.end());

    Vector3d newR12 = i->second * r12;
    // x, y and z's possible values range -halfLen_ to halfLen_
    int xbin = int((newR12.x() + halfLen_) / deltaR_);
    int ybin = int((newR12.y() + halfLen_) / deltaR_);
    int zbin = int((newR12.z() + halfLen_) / deltaR_);

    if (xbin < int(nBins_) && xbin >= 0 && ybin < int(nBins_) && ybin >= 0 &&
        zbin < int(nBins_) && zbin >= 0) {
      ++histogram_[xbin][ybin][zbin];
    }
  }

  void GofXyz::writeRdf() {
    std::ofstream rdfStream(outputFilename_.c_str(), std::ios::binary);
    if (rdfStream.is_open()) {
      // rdfStream << "#g(x, y, z)\n";
      // rdfStream << "#selection1: (" << selectionScript1_ << ")\t";
      // rdfStream << "selection2: (" << selectionScript2_ << ")\n";
      // rdfStream << "#nRBins = " << nBins_ << "\t maxLen = "
      //          << len_ << "deltaR = " << deltaR_ <<"\n";
      for (unsigned int i = 0; i < histogram_.size(); ++i) {
        for (unsigned int j = 0; j < histogram_[i].size(); ++j) {
          for (unsigned int k = 0; k < histogram_[i][j].size(); ++k) {
            rdfStream.write(reinterpret_cast<char*>(&histogram_[i][j][k]),
                            sizeof(histogram_[i][j][k]));
          }
        }
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "GofXyz: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }

}  // namespace OpenMD
