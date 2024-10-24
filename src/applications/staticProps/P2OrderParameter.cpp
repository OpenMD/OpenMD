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

#include "applications/staticProps/P2OrderParameter.hpp"

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {

  P2OrderParameter::P2OrderParameter(SimInfo* info, const string& filename,
                                     const string& sele1) :
      StaticAnalyser(info, filename, 1),
      doVect_(true), doOffset_(false), selectionScript1_(sele1),
      seleMan1_(info), seleMan2_(info), evaluator1_(info), evaluator2_(info) {
    setOutputName(getPrefix(filename) + ".p2");

    evaluator1_.loadScriptString(sele1);
  }

  P2OrderParameter::P2OrderParameter(SimInfo* info, const string& filename,
                                     const string& sele1, const string& sele2) :
      StaticAnalyser(info, filename, 1),
      doVect_(false), doOffset_(false), selectionScript1_(sele1),
      selectionScript2_(sele2), seleMan1_(info), seleMan2_(info),
      evaluator1_(info), evaluator2_(info) {
    setOutputName(getPrefix(filename) + ".p2");

    evaluator1_.loadScriptString(sele1);
    evaluator2_.loadScriptString(sele2);
  }

  P2OrderParameter::P2OrderParameter(SimInfo* info, const string& filename,
                                     const string& sele1, int seleOffset) :
      StaticAnalyser(info, filename, 1),
      doVect_(false), doOffset_(true), selectionScript1_(sele1),
      seleMan1_(info), seleMan2_(info), evaluator1_(info), evaluator2_(info),
      seleOffset_(seleOffset) {
    setOutputName(getPrefix(filename) + ".p2");

    evaluator1_.loadScriptString(sele1);
  }

  void P2OrderParameter::process() {
    StuntDouble* sd1;
    StuntDouble* sd2;
    int ii;
    int jj;
    int vecCount;
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();

    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      Mat3x3d orderTensor(0.0);
      vecCount = 0;

      seleMan1_.setSelectionSet(evaluator1_.evaluate());

      if (doVect_) {
        for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL;
             sd1 = seleMan1_.nextSelected(ii)) {
          if (sd1->isDirectional()) {
            Vector3d vec = sd1->getA().transpose() * V3Z;

            vec.normalize();
            orderTensor += outProduct(vec, vec);
            vecCount++;
          }
        }

        orderTensor /= vecCount;

      } else {
        if (doOffset_) {
          for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL;
               sd1 = seleMan1_.nextSelected(ii)) {
            // This will require careful rewriting if StaticProps is
            // ever parallelized.  For an example, see
            // Thermo::getTaggedAtomPairDistance

            int sd2Index = sd1->getGlobalIndex() + seleOffset_;
            sd2          = info_->getIOIndexToIntegrableObject(sd2Index);

            Vector3d vec = sd1->getPos() - sd2->getPos();

            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);

            vec.normalize();

            orderTensor += outProduct(vec, vec);
            vecCount++;
          }

          orderTensor /= vecCount;
        } else {
          seleMan2_.setSelectionSet(evaluator2_.evaluate());

          if (seleMan1_.getSelectionCount() != seleMan2_.getSelectionCount()) {
            snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                     "In frame %d, the number of selected StuntDoubles are\n"
                     "\tnot the same in --sele1 and sele2\n",
                     i);
            painCave.severity = OPENMD_INFO;
            painCave.isFatal  = 0;
            simError();
          }

          for (sd1                             = seleMan1_.beginSelected(ii),
              sd2                              = seleMan2_.beginSelected(jj);
               sd1 != NULL && sd2 != NULL; sd1 = seleMan1_.nextSelected(ii),
              sd2                              = seleMan2_.nextSelected(jj)) {
            Vector3d vec = sd1->getPos() - sd2->getPos();

            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);

            vec.normalize();

            orderTensor += outProduct(vec, vec);
            vecCount++;
          }

          orderTensor /= vecCount;
        }
      }

      if (vecCount == 0) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "In frame %d, the number of selected vectors was zero.\n"
                 "\tThis will not give a meaningful order parameter.",
                 i);
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal  = 1;
        simError();
      }

      orderTensor -= (RealType)(1.0 / 3.0) * Mat3x3d::identity();

      Vector3d eigenvalues;
      Mat3x3d eigenvectors;

      Mat3x3d::diagonalize(orderTensor, eigenvalues, eigenvectors);

      int which(-1);
      RealType maxEval = 0.0;
      for (int k = 0; k < 3; k++) {
        if (fabs(eigenvalues[k]) > maxEval) {
          which   = k;
          maxEval = fabs(eigenvalues[k]);
        }
      }
      RealType p2 = 1.5 * maxEval;

      // the eigen vector is already normalized in SquareMatrix3::diagonalize
      Vector3d director = eigenvectors.getColumn(which);
      // if (director[2] < 0) { director.negate(); }

      RealType angle = 0.0;
      vecCount       = 0;

      if (doVect_) {
        for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL;
             sd1 = seleMan1_.nextSelected(ii)) {
          if (sd1->isDirectional()) {
            Vector3d vec = sd1->getA().transpose() * V3Z;
            vec.normalize();
            angle += acos(dot(vec, director));
            vecCount++;
          }
        }
        angle = angle / (vecCount * Constants::PI) * 180.0;

      } else {
        if (doOffset_) {
          for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL;
               sd1 = seleMan1_.nextSelected(ii)) {
            // This will require careful rewriting if StaticProps is
            // ever parallelized.  For an example, see
            // Thermo::getTaggedAtomPairDistance

            int sd2Index = sd1->getGlobalIndex() + seleOffset_;
            sd2          = info_->getIOIndexToIntegrableObject(sd2Index);

            Vector3d vec = sd1->getPos() - sd2->getPos();
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);
            vec.normalize();
            angle += acos(dot(vec, director));
            vecCount++;
          }
          angle = angle / (vecCount * Constants::PI) * 180.0;

        } else {
          for (sd1                             = seleMan1_.beginSelected(ii),
              sd2                              = seleMan2_.beginSelected(jj);
               sd1 != NULL && sd2 != NULL; sd1 = seleMan1_.nextSelected(ii),
              sd2                              = seleMan2_.nextSelected(jj)) {
            Vector3d vec = sd1->getPos() - sd2->getPos();
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);
            vec.normalize();
            angle += acos(dot(vec, director));
            vecCount++;
          }
          angle = angle / (vecCount * Constants::PI) * 180.0;
        }
      }

      OrderParam param;
      param.p2       = p2;
      param.director = director;
      param.angle    = angle;

      orderParams_.push_back(param);
    }

    writeP2();
  }

  void P2OrderParameter::writeP2() {
    ofstream os(getOutputFileName().c_str());
    os << "#P2 Order parameter\n";
    os << "#selection1: (" << selectionScript1_ << ")\t";
    if (!doVect_) { os << "selection2: (" << selectionScript2_ << ")\n"; }
    os << "#p2\tdirector_x\tdirector_y\tdirector_z\tangle(degree)\n";

    RealType p2Sum {};
    RealType angleSum {};
    Vector3d directorSum(0.0);

    for (size_t i = 0; i < orderParams_.size(); ++i) {
      p2Sum += orderParams_[i].p2;
      directorSum += orderParams_[i].director;
      angleSum += orderParams_[i].angle;
    }

    p2Sum /= orderParams_.size();
    directorSum /= orderParams_.size();
    angleSum /= orderParams_.size();

    os << p2Sum << "\t" << directorSum[0] << "\t" << directorSum[1] << "\t"
       << directorSum[2] << "\t" << angleSum << "\n";
  }

}  // namespace OpenMD
