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

/* Calculates Angle(R) for DirectionalAtoms*/

#include "applications/staticProps/P2R.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

namespace OpenMD {

  P2R::P2R(SimInfo* info, const std::string& filename, const std::string& sele1,
           unsigned int nbins) :
      StaticAnalyser(info, filename, nbins),
      doVect_(true), doOffset_(false), selectionScript1_(sele1),
      seleMan1_(info), seleMan2_(info), evaluator1_(info), evaluator2_(info) {
    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    setAnalysisType(
        "2nd order Legendre Polynomial Correlation using r as reference axis");
    setOutputName(getPrefix(filename) + ".P2R");
    std::stringstream params;
    const std::string paramString = params.str();
    setParameterString(paramString);
  }

  P2R::P2R(SimInfo* info, const std::string& filename, const std::string& sele1,
           const std::string& sele2, unsigned int nbins) :
      StaticAnalyser(info, filename, nbins),
      doVect_(false), doOffset_(false), selectionScript1_(sele1),
      selectionScript2_(sele2), seleMan1_(info), seleMan2_(info),
      evaluator1_(info), evaluator2_(info) {
    setOutputName(getPrefix(filename) + ".P2R");

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    evaluator2_.loadScriptString(sele2);
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }

    setAnalysisType(
        "2nd order Legendre Polynomial Correlation using r as reference axis");
    setOutputName(getPrefix(filename) + ".P2R");
    std::stringstream params;
    const std::string paramString = params.str();
    setParameterString(paramString);
  }

  P2R::P2R(SimInfo* info, const std::string& filename, const std::string& sele1,
           int seleOffset, unsigned int nbins) :
      StaticAnalyser(info, filename, nbins),
      doVect_(false), doOffset_(true), doOffset2_(false),
      selectionScript1_(sele1), seleMan1_(info), seleMan2_(info),
      evaluator1_(info), evaluator2_(info), seleOffset_(seleOffset) {
    setOutputName(getPrefix(filename) + ".P2R");

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    setAnalysisType(
        "2nd order Legendre Polynomial Correlation using r as reference axis");
    setOutputName(getPrefix(filename) + ".P2R");
    std::stringstream params;
    const std::string paramString = params.str();
    setParameterString(paramString);
  }

  P2R::P2R(SimInfo* info, const std::string& filename, const std::string& sele1,
           int seleOffset, int seleOffset2, unsigned int nbins) :
      StaticAnalyser(info, filename, nbins),
      doVect_(false), doOffset_(true), doOffset2_(true),
      selectionScript1_(sele1), seleMan1_(info), seleMan2_(info),
      evaluator1_(info), evaluator2_(info), seleOffset_(seleOffset),
      seleOffset2_(seleOffset2) {
    setOutputName(getPrefix(filename) + ".P2R");

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    setAnalysisType(
        "2nd order Legendre Polynomial Correlation using r as reference axis");
    setOutputName(getPrefix(filename) + ".P2R");
    std::stringstream params;
    const std::string paramString = params.str();
    setParameterString(paramString);
  }

  void P2R::process() {
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
    P2_ = 0.0;
    count_ =  0;

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

            P2_ += 0.5 * (3.0 * cosangle * cosangle - 1.0);
            count_++;
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
            P2_ += 0.5 * (3.0 * cosangle * cosangle - 1.0);
            count_++;
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
            P2_ += 0.5 * (3.0 * cosangle * cosangle - 1.0);
            count_++;
          }
        }
      }
    }
    processHistogram();
    writeP2R();
  }

  void P2R::processHistogram() {
    if (count_ > 0) P2_ /= count_;
  }

  void P2R::writeP2R() {
    std::ofstream ofs(outputFilename_.c_str());
    if (ofs.is_open()) {
      Revision rev;

      ofs << "# " << getAnalysisType() << "\n";
      ofs << "# OpenMD " << rev.getFullRevision() << "\n";
      ofs << "# " << rev.getBuildDate() << "\n";
      ofs << "#nFrames:\t" << nProcessed_ << "\n";
      ofs << "#selection1: (" << selectionScript1_ << ")";
      if (!doVect_) { ofs << "\tselection2: (" << selectionScript2_ << ")"; }
      ofs << "\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";

      ofs << "#P2 correlation:\n";
      ofs << P2_ << "\n";

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "P2R: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }

}  // namespace OpenMD
