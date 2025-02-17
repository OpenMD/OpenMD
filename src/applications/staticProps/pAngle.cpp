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

/* Calculates Rho(theta) */

#include "applications/staticProps/pAngle.hpp"

#include <algorithm>
#include <fstream>

#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

namespace OpenMD {

  pAngle::pAngle(SimInfo* info, const std::string& filename,
                 const std::string& sele1, int nthetabins) :
      StaticAnalyser(info, filename, nthetabins),
      doVect_(true), doOffset_(false), selectionScript1_(sele1),
      seleMan1_(info), seleMan2_(info), evaluator1_(info), evaluator2_(info),
      nThetaBins_(nthetabins) {
    setOutputName(getPrefix(filename) + ".pAngle");

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    count_.resize(nThetaBins_);
    histogram_.resize(nThetaBins_);
  }

  pAngle::pAngle(SimInfo* info, const std::string& filename,
                 const std::string& sele1, const std::string& sele2,
                 int nthetabins) :
      StaticAnalyser(info, filename, nthetabins),
      doVect_(false), doOffset_(false), selectionScript1_(sele1),
      selectionScript2_(sele2), seleMan1_(info), seleMan2_(info),
      evaluator1_(info), evaluator2_(info), nThetaBins_(nthetabins) {
    setOutputName(getPrefix(filename) + ".pAngle");

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    evaluator2_.loadScriptString(sele2);
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }

    count_.resize(nThetaBins_);
    histogram_.resize(nThetaBins_);
  }

  pAngle::pAngle(SimInfo* info, const std::string& filename,
                 const std::string& sele1, int seleOffset, int nthetabins) :
      StaticAnalyser(info, filename, nthetabins),
      doVect_(false), doOffset_(true), doOffset2_(false),
      selectionScript1_(sele1), seleMan1_(info), seleMan2_(info),
      evaluator1_(info), evaluator2_(info), seleOffset_(seleOffset),
      nThetaBins_(nthetabins) {
    setOutputName(getPrefix(filename) + ".pAngle");

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    count_.resize(nThetaBins_);
    histogram_.resize(nThetaBins_);
  }

  pAngle::pAngle(SimInfo* info, const std::string& filename,
                 const std::string& sele1, int seleOffset, int seleOffset2,
                 int nthetabins) :
      StaticAnalyser(info, filename, nthetabins),
      doVect_(false), doOffset_(true), doOffset2_(true),
      selectionScript1_(sele1), seleMan1_(info), seleMan2_(info),
      evaluator1_(info), evaluator2_(info), seleOffset_(seleOffset),
      seleOffset2_(seleOffset2), nThetaBins_(nthetabins) {
    setOutputName(getPrefix(filename) + ".pAngle");

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    count_.resize(nThetaBins_);
    histogram_.resize(nThetaBins_);
  }

  void pAngle::process() {
    StuntDouble* sd1;
    StuntDouble* sd2;
    int ii;
    int jj;
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    Thermo thermo(info_);
    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();

    nProcessed_ = nFrames / step_;

    std::fill(histogram_.begin(), histogram_.end(), 0.0);
    std::fill(count_.begin(), count_.end(), 0);

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      Vector3d CenterOfMass = thermo.getCom();

      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }

      if (doVect_) {
        for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL;
             sd1 = seleMan1_.nextSelected(ii)) {
          Vector3d pos = sd1->getPos();

          Vector3d r1 = CenterOfMass - pos;
          // only do this if the stunt double actually has a vector associated
          // with it
          if (sd1->isDirectional()) {
            Vector3d vec = sd1->getA().transpose() * V3Z;
            vec.normalize();
            r1.normalize();
            RealType cosangle = dot(r1, vec);

            int binNo = int(nThetaBins_ * (1.0 + cosangle) / 2.0);
            count_[binNo]++;
          }
        }
      } else {
        if (doOffset_) {
          for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL;
               sd1 = seleMan1_.nextSelected(ii)) {
            // This will require careful rewriting if StaticProps is
            // ever parallelized.  For an example, see
            // Thermo::getTaggedAtomPairDistance
            Vector3d r1;

            if (doOffset2_) {
              int sd1Aind       = sd1->getGlobalIndex() + seleOffset2_;
              StuntDouble* sd1A = info_->getIOIndexToIntegrableObject(sd1Aind);
              r1                = CenterOfMass - sd1A->getPos();
            } else {
              r1 = CenterOfMass - sd1->getPos();
            }

            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(r1);

            int sd2Index = sd1->getGlobalIndex() + seleOffset_;
            sd2          = info_->getIOIndexToIntegrableObject(sd2Index);

            Vector3d r2 = CenterOfMass - sd2->getPos();
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(r1);

            Vector3d rc = 0.5 * (r1 + r2);
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(rc);

            Vector3d vec = r1 - r2;
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);

            rc.normalize();
            vec.normalize();
            RealType cosangle = dot(rc, vec);
            int binNo         = int(nThetaBins_ * (1.0 + cosangle) / 2.0);
            count_[binNo]++;
          }
        } else {
          if (evaluator2_.isDynamic()) {
            seleMan2_.setSelectionSet(evaluator2_.evaluate());
          }

          if (seleMan1_.getSelectionCount() != seleMan2_.getSelectionCount()) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "In frame %d, the number of selected StuntDoubles are\n"
                     "\tnot the same in --sele1 and sele2\n",
                     istep);
            painCave.severity = OPENMD_INFO;
            painCave.isFatal  = 0;
            simError();
          }

          for (sd1                             = seleMan1_.beginSelected(ii),
              sd2                              = seleMan2_.beginSelected(jj);
               sd1 != NULL && sd2 != NULL; sd1 = seleMan1_.nextSelected(ii),
              sd2                              = seleMan2_.nextSelected(jj)) {
            Vector3d r1 = CenterOfMass - sd1->getPos();
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(r1);

            Vector3d r2 = CenterOfMass - sd2->getPos();
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(r1);

            Vector3d rc = 0.5 * (r1 + r2);
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(rc);

            Vector3d vec = r1 - r2;
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);

            rc.normalize();
            vec.normalize();
            RealType cosangle = dot(rc, vec);
            int binNo         = int(nThetaBins_ * (1.0 + cosangle) / 2.0);
            count_[binNo]++;
          }
        }
      }
    }

    processHistogram();
    writeProbs();
  }

  void pAngle::processHistogram() {
    int atot = 0;
    for (unsigned int i = 0; i < count_.size(); ++i)
      atot += count_[i];

    for (unsigned int i = 0; i < count_.size(); ++i) {
      histogram_[i] = double(count_[i] / double(atot));
    }
  }

  void pAngle::writeProbs() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#pAngle\n";
      rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "#selection1: (" << selectionScript1_ << ")";
      if (!doVect_) {
        rdfStream << "\tselection2: (" << selectionScript2_ << ")";
      }
      rdfStream << "\n";
      rdfStream << "#cos(theta)\tp(cos(theta))\n";
      RealType dct = 2.0 / histogram_.size();
      for (unsigned int i = 0; i < histogram_.size(); ++i) {
        RealType ct = -1.0 + (2.0 * i + 1) / (histogram_.size());
        rdfStream << ct << "\t" << histogram_[i] / dct << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "pAngle: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }

}  // namespace OpenMD
