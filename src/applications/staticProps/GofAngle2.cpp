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

#include "applications/staticProps/GofAngle2.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>

#include "primitives/Atom.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {

  GofAngle2::GofAngle2(SimInfo* info, const std::string& filename,
                       const std::string& sele1, const std::string& sele2,
                       int nangleBins) :
      RadialDistrFunc(info, filename, sele1, sele2, nangleBins),
      doSele3_(false), seleMan3_(info), evaluator3_(info) {
    setAnalysisType("Radial Distribution Function");
    setOutputName(getPrefix(filename) + ".gto");

    deltaCosAngle_ = 2.0 / nBins_;

    std::stringstream params;
    params << " nAngleBins = " << nBins_
           << ", deltaCosAngle = " << deltaCosAngle_;
    const std::string paramString = params.str();
    setParameterString(paramString);

    histogram_.resize(nBins_);
    avgGofr_.resize(nBins_);
    for (unsigned int i = 0; i < nBins_; ++i) {
      histogram_[i].resize(nBins_);
      avgGofr_[i].resize(nBins_);
    }
  }

  GofAngle2::GofAngle2(SimInfo* info, const std::string& filename,
                       const std::string& sele1, const std::string& sele2,
                       const std::string& sele3, int nangleBins) :
      RadialDistrFunc(info, filename, sele1, sele2, nangleBins),
      doSele3_(true), seleMan3_(info), evaluator3_(info),
      selectionScript3_(sele3) {
    setOutputName(getPrefix(filename) + ".gto");

    deltaCosAngle_ = 2.0 / nBins_;

    histogram_.resize(nBins_);
    avgGofr_.resize(nBins_);
    for (unsigned int i = 0; i < nBins_; ++i) {
      histogram_[i].resize(nBins_);
      avgGofr_[i].resize(nBins_);
    }
    evaluator3_.loadScriptString(sele3);
    if (!evaluator3_.isDynamic()) {
      seleMan3_.setSelectionSet(evaluator3_.evaluate());
    }
  }

  void GofAngle2::processNonOverlapping(SelectionManager& sman1,
                                        SelectionManager& sman2) {
    StuntDouble* sd1;
    StuntDouble* sd2;
    StuntDouble* sd3;
    int i;
    int j;
    int k;

    // This is the same as a non-overlapping pairwise loop structure:
    // for (int i = 0;  i < ni ; ++i ) {
    //   for (int j = 0; j < nj; ++j) {}
    // }

    if (doSele3_) {
      if (evaluator3_.isDynamic()) {
        seleMan3_.setSelectionSet(evaluator3_.evaluate());
      }
      if (sman1.getSelectionCount() != seleMan3_.getSelectionCount()) {
        RadialDistrFunc::processNonOverlapping(sman1, sman2);
      }

      for (sd1 = sman1.beginSelected(i), sd3 = seleMan3_.beginSelected(k);
           sd1 != NULL && sd3 != NULL;
           sd1 = sman1.nextSelected(i), sd3 = seleMan3_.nextSelected(k)) {
        for (sd2 = sman2.beginSelected(j); sd2 != NULL;
             sd2 = sman2.nextSelected(j)) {
          collectHistogram(sd1, sd2, sd3);
        }
      }
    } else {
      RadialDistrFunc::processNonOverlapping(sman1, sman2);
    }
  }

  void GofAngle2::processOverlapping(SelectionManager& sman) {
    StuntDouble* sd1;
    StuntDouble* sd2;
    StuntDouble* sd3;
    int i;
    int j;
    int k;

    // This is the same as a pairwise loop structure:
    // for (int i = 0;  i < n-1 ; ++i ) {
    //   for (int j = i + 1; j < n; ++j) {}
    // }

    if (doSele3_) {
      if (evaluator3_.isDynamic()) {
        seleMan3_.setSelectionSet(evaluator3_.evaluate());
      }
      if (sman.getSelectionCount() != seleMan3_.getSelectionCount()) {
        RadialDistrFunc::processOverlapping(sman);
      }
      for (sd1 = sman.beginSelected(i), sd3 = seleMan3_.beginSelected(k);
           sd1 != NULL && sd3 != NULL;
           sd1 = sman.nextSelected(i), sd3 = seleMan3_.nextSelected(k)) {
        for (j = i, sd2 = sman.nextSelected(j); sd2 != NULL;
             sd2 = sman.nextSelected(j)) {
          collectHistogram(sd1, sd2, sd3);
        }
      }
    } else {
      RadialDistrFunc::processOverlapping(sman);
    }
  }

  void GofAngle2::preProcess() {
    for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
      std::fill(avgGofr_[i].begin(), avgGofr_[i].end(), 0);
    }
  }

  void GofAngle2::initializeHistogram() {
    npairs_ = 0;
    for (unsigned int i = 0; i < histogram_.size(); ++i)
      std::fill(histogram_[i].begin(), histogram_[i].end(), 0);
  }

  void GofAngle2::processHistogram() {
    // std::for_each(avgGofr_.begin(), avgGofr_.end(),
    // std::plus<std::vector<int>>)
  }

  void GofAngle2::collectHistogram(StuntDouble* sd1, StuntDouble* sd2) {
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    if (sd1 == sd2) { return; }

    Vector3d pos1 = sd1->getPos();
    Vector3d pos2 = sd2->getPos();
    Vector3d r12  = pos1 - pos2;
    if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(r12);

    AtomType* atype1     = static_cast<Atom*>(sd1)->getAtomType();
    AtomType* atype2     = static_cast<Atom*>(sd2)->getAtomType();
    MultipoleAdapter ma1 = MultipoleAdapter(atype1);
    MultipoleAdapter ma2 = MultipoleAdapter(atype2);

    if (!sd1->isDirectional()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "GofAngle2: attempted to use a non-directional object: %s\n",
               sd1->getType().c_str());
      painCave.isFatal = 1;
      simError();
    }

    if (!sd2->isDirectional()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "GofAngle2: attempted to use a non-directional object: %s\n",
               sd2->getType().c_str());
      painCave.isFatal = 1;
      simError();
    }

    Vector3d dipole1, dipole2;
    if (ma1.isDipole())
      dipole1 = sd1->getDipole();
    else
      dipole1 = sd1->getA().transpose() * V3Z;

    if (ma2.isDipole())
      dipole2 = sd2->getDipole();
    else
      dipole2 = sd2->getA().transpose() * V3Z;

    r12.normalize();
    dipole1.normalize();
    dipole2.normalize();

    RealType cosAngle1 = dot(r12, dipole1);
    RealType cosAngle2 = dot(dipole1, dipole2);

    RealType halfBin = (nBins_ - 1) * 0.5;
    int angleBin1    = int(halfBin * (cosAngle1 + 1.0));
    int angleBin2    = int(halfBin * (cosAngle2 + 1.0));

    ++histogram_[angleBin1][angleBin2];
    ++npairs_;
  }

  void GofAngle2::collectHistogram(StuntDouble* sd1, StuntDouble* sd2,
                                   StuntDouble* sd3) {
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    if (sd1 == sd2) { return; }

    Vector3d p1 = sd1->getPos();
    Vector3d p3 = sd3->getPos();

    Vector3d c   = 0.5 * (p1 + p3);
    Vector3d r13 = p3 - p1;

    Vector3d r12 = sd2->getPos() - c;

    if (usePeriodicBoundaryConditions_) {
      currentSnapshot_->wrapVector(r12);
      currentSnapshot_->wrapVector(r13);
    }
    r12.normalize();
    r13.normalize();

    if (!sd2->isDirectional()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "GofAngle2: attempted to use a non-directional object: %s\n",
               sd2->getType().c_str());
      painCave.isFatal = 1;
      simError();
    }

    AtomType* atype2     = static_cast<Atom*>(sd2)->getAtomType();
    MultipoleAdapter ma2 = MultipoleAdapter(atype2);

    Vector3d dipole2;
    if (ma2.isDipole())
      dipole2 = sd2->getDipole();
    else
      dipole2 = sd2->getA().transpose() * V3Z;

    dipole2.normalize();

    RealType cosAngle1 = dot(r12, r13);
    RealType cosAngle2 = dot(r13, dipole2);

    RealType halfBin = (nBins_ - 1) * 0.5;
    int angleBin1    = int(halfBin * (cosAngle1 + 1.0));
    int angleBin2    = int(halfBin * (cosAngle2 + 1.0));

    ++histogram_[angleBin1][angleBin2];
    ++npairs_;
  }

  void GofAngle2::writeRdf() {
    std::ofstream ofs(outputFilename_.c_str());
    if (ofs.is_open()) {
      Revision r;
      ofs << "# " << getAnalysisType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"";
      if (doSele3_) {
        ofs << "\tselection script3: \"" << selectionScript3_ << "\"\n";
      } else {
        ofs << "\n";
      }

      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      for (unsigned int i = 0; i < avgGofr_.size(); ++i) {
        // RealType cosAngle1 = -1.0 + (i + 0.5)*deltaCosAngle_;

        for (unsigned int j = 0; j < avgGofr_[i].size(); ++j) {
          // RealType cosAngle2 = -1.0 + (j + 0.5)*deltaCosAngle_;
          ofs << avgGofr_[i][j] / nProcessed_ << "\t";
        }
        ofs << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "GofAngle2: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }

}  // namespace OpenMD
