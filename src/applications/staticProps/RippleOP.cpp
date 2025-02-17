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

#include "applications/staticProps/RippleOP.hpp"

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  RippleOP::RippleOP(SimInfo* info, const std::string& filename,
                     const std::string& sele1, const std::string& sele2) :
      StaticAnalyser(info, filename, 1),
      selectionScript1_(sele1), selectionScript2_(sele2), seleMan1_(info),
      seleMan2_(info), evaluator1_(info), evaluator2_(info) {
    setOutputName(getPrefix(filename) + ".rp2");

    evaluator1_.loadScriptString(sele1);
    evaluator2_.loadScriptString(sele2);

    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "--sele1 must be static selection\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "--sele2 must be static selection\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    if (seleMan1_.getSelectionCount() != seleMan2_.getSelectionCount()) {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "The number of selected Stuntdoubles are not the same in --sele1 "
          "and sele2\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    int i;
    int j;
    StuntDouble* sd1;
    StuntDouble* sd2;
    for (sd1 = seleMan1_.beginSelected(i), sd2 = seleMan2_.beginSelected(j);
         sd1 != NULL && sd2 != NULL;
         sd1 = seleMan1_.nextSelected(i), sd2 = seleMan2_.nextSelected(j)) {
      sdPairs_.push_back(std::make_pair(sd1, sd2));
    }
  }

  void RippleOP::process() {
    StuntDouble* j1;
    StuntDouble* j2;
    StuntDouble* sd3;
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();

    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      int nMolecules   = info_->getNGlobalMolecules();
      int i1;
      int nUpper    = 0;
      int nLower    = 0;
      int nTail     = 0;
      RealType sumZ = 0.0;

      for (sd3 = seleMan2_.beginSelected(i1); sd3 != NULL;
           sd3 = seleMan2_.nextSelected(i1)) {
        Vector3d pos1 = sd3->getPos();
        if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(pos1);
        sd3->setPos(pos1);
      }

      for (sd3 = seleMan2_.beginSelected(i1); sd3 != NULL;
           sd3 = seleMan2_.nextSelected(i1)) {
        Vector3d pos1 = sd3->getPos();
        sumZ += pos1.z();
      }
      RealType avgZ = sumZ / (RealType)nMolecules;

      Mat3x3d orderTensorHeadUpper(0.0);
      Mat3x3d orderTensorTail(0.0);
      Mat3x3d orderTensorHeadLower(0.0);
      for (j1 = seleMan1_.beginSelected(i1); j1 != NULL;
           j1 = seleMan1_.nextSelected(i1)) {
        Vector3d pos = j1->getPos();
        if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(pos);
        Vector3d vecHeadUpper;
        if (pos.z() >= avgZ) {
          AtomType* atype1     = static_cast<Atom*>(j1)->getAtomType();
          MultipoleAdapter ma1 = MultipoleAdapter(atype1);
          if (ma1.isDipole())
            vecHeadUpper = j1->getDipole();
          else
            vecHeadUpper = j1->getA().transpose() * V3Z;
          nUpper++;
        }
        Vector3d vecHeadLower;
        if (pos.z() <= avgZ) {
          AtomType* atype1     = static_cast<Atom*>(j1)->getAtomType();
          MultipoleAdapter ma1 = MultipoleAdapter(atype1);
          if (ma1.isDipole())
            vecHeadLower = j1->getDipole();
          else
            vecHeadLower = j1->getA().transpose() * V3Z;
          nLower++;
        }
        orderTensorHeadUpper += outProduct(vecHeadUpper, vecHeadUpper);
        orderTensorHeadLower += outProduct(vecHeadLower, vecHeadLower);
      }
      for (j2 = seleMan2_.beginSelected(i1); j2 != NULL;
           j2 = seleMan2_.nextSelected(i1)) {
        // The lab frame vector corresponding to the body-fixed
        // z-axis is simply the second column of A.transpose()
        // or, identically, the second row of A itself.

        Vector3d vecTail = j2->getA().transpose() * V3Z;
        orderTensorTail += outProduct(vecTail, vecTail);
        nTail++;
      }

      orderTensorHeadUpper /= (RealType)nUpper;
      orderTensorHeadLower /= (RealType)nLower;
      orderTensorHeadUpper -= (RealType)(1.0 / 3.0) * Mat3x3d::identity();
      orderTensorHeadLower -= (RealType)(1.0 / 3.0) * Mat3x3d::identity();

      orderTensorTail /= (RealType)nTail;
      orderTensorTail -= (RealType)(1.0 / 3.0) * Mat3x3d::identity();

      Vector3d eigenvaluesHeadUpper, eigenvaluesHeadLower, eigenvaluesTail;
      Mat3x3d eigenvectorsHeadUpper, eigenvectorsHeadLower, eigenvectorsTail;
      Mat3x3d::diagonalize(orderTensorHeadUpper, eigenvaluesHeadUpper,
                           eigenvectorsHeadUpper);
      Mat3x3d::diagonalize(orderTensorHeadLower, eigenvaluesHeadLower,
                           eigenvectorsHeadLower);
      Mat3x3d::diagonalize(orderTensorTail, eigenvaluesTail, eigenvectorsTail);

      int whichUpper(-1), whichLower(-1), whichTail(-1);
      RealType maxEvalUpper = 0.0;
      RealType maxEvalLower = 0.0;
      RealType maxEvalTail  = 0.0;
      for (int k = 0; k < 3; k++) {
        if (fabs(eigenvaluesHeadUpper[k]) > maxEvalUpper) {
          whichUpper   = k;
          maxEvalUpper = fabs(eigenvaluesHeadUpper[k]);
        }
      }
      RealType p2HeadUpper = 1.5 * maxEvalUpper;
      for (int k = 0; k < 3; k++) {
        if (fabs(eigenvaluesHeadLower[k]) > maxEvalLower) {
          whichLower   = k;
          maxEvalLower = fabs(eigenvaluesHeadLower[k]);
        }
      }
      RealType p2HeadLower = 1.5 * maxEvalLower;
      for (int k = 0; k < 3; k++) {
        if (fabs(eigenvaluesTail[k]) > maxEvalTail) {
          whichTail   = k;
          maxEvalTail = fabs(eigenvaluesTail[k]);
        }
      }
      RealType p2Tail = 1.5 * maxEvalTail;

      // the eigenvector is already normalized in SquareMatrix3::diagonalize
      Vector3d directorHeadUpper = eigenvectorsHeadUpper.getColumn(whichUpper);
      if (directorHeadUpper[0] < 0) { directorHeadUpper.negate(); }
      Vector3d directorHeadLower = eigenvectorsHeadLower.getColumn(whichLower);
      if (directorHeadLower[0] < 0) { directorHeadLower.negate(); }
      Vector3d directorTail = eigenvectorsTail.getColumn(whichTail);
      if (directorTail[0] < 0) { directorTail.negate(); }

      OrderParam paramHeadUpper, paramHeadLower, paramTail;
      paramHeadUpper.p2       = p2HeadUpper;
      paramHeadUpper.director = directorHeadUpper;
      paramHeadLower.p2       = p2HeadLower;
      paramHeadLower.director = directorHeadLower;
      paramTail.p2            = p2Tail;
      paramTail.director      = directorTail;

      orderParamsHeadUpper_.push_back(paramHeadUpper);
      orderParamsHeadLower_.push_back(paramHeadLower);
      orderParamsTail_.push_back(paramTail);
    }

    OrderParam sumOPHeadUpper, errsumOPHeadUpper;
    OrderParam sumOPHeadLower, errsumOPHeadLower;
    OrderParam sumOPTail, errsumOPTail;

    sumOPHeadUpper.p2    = 0.0;
    errsumOPHeadUpper.p2 = 0.0;
    sumOPHeadLower.p2    = 0.0;
    errsumOPHeadLower.p2 = 0.0;
    for (std::size_t i = 0; i < orderParamsHeadUpper_.size(); ++i) {
      sumOPHeadUpper.p2 += orderParamsHeadUpper_[i].p2;
      sumOPHeadUpper.director[0] += orderParamsHeadUpper_[i].director[0];
      sumOPHeadUpper.director[1] += orderParamsHeadUpper_[i].director[1];
      sumOPHeadUpper.director[2] += orderParamsHeadUpper_[i].director[2];
    }

    avgOPHeadUpper.p2 =
        sumOPHeadUpper.p2 / (RealType)orderParamsHeadUpper_.size();
    avgOPHeadUpper.director[0] =
        sumOPHeadUpper.director[0] / (RealType)orderParamsHeadUpper_.size();
    avgOPHeadUpper.director[1] =
        sumOPHeadUpper.director[1] / (RealType)orderParamsHeadUpper_.size();
    avgOPHeadUpper.director[2] =
        sumOPHeadUpper.director[2] / (RealType)orderParamsHeadUpper_.size();
    for (std::size_t i = 0; i < orderParamsHeadUpper_.size(); ++i) {
      errsumOPHeadUpper.p2 +=
          pow((orderParamsHeadUpper_[i].p2 - avgOPHeadUpper.p2), 2);
    }
    errOPHeadUpper.p2 =
        sqrt(errsumOPHeadUpper.p2 / (RealType)orderParamsHeadUpper_.size());
    for (std::size_t i = 0; i < orderParamsHeadLower_.size(); ++i) {
      sumOPHeadLower.p2 += orderParamsHeadLower_[i].p2;
      sumOPHeadLower.director[0] += orderParamsHeadLower_[i].director[0];
      sumOPHeadLower.director[1] += orderParamsHeadLower_[i].director[1];
      sumOPHeadLower.director[2] += orderParamsHeadLower_[i].director[2];
    }

    avgOPHeadLower.p2 =
        sumOPHeadLower.p2 / (RealType)orderParamsHeadLower_.size();
    avgOPHeadLower.director[0] =
        sumOPHeadLower.director[0] / (RealType)orderParamsHeadLower_.size();
    avgOPHeadLower.director[1] =
        sumOPHeadLower.director[1] / (RealType)orderParamsHeadLower_.size();
    avgOPHeadLower.director[2] =
        sumOPHeadLower.director[2] / (RealType)orderParamsHeadLower_.size();
    for (std::size_t i = 0; i < orderParamsHeadLower_.size(); ++i) {
      errsumOPHeadLower.p2 +=
          pow((orderParamsHeadLower_[i].p2 - avgOPHeadLower.p2), 2);
    }
    errOPHeadLower.p2 =
        sqrt(errsumOPHeadLower.p2 / (RealType)orderParamsHeadLower_.size());

    sumOPTail.p2    = 0.0;
    errsumOPTail.p2 = 0.0;
    for (std::size_t i = 0; i < orderParamsTail_.size(); ++i) {
      sumOPTail.p2 += orderParamsTail_[i].p2;
      sumOPTail.director[0] += orderParamsTail_[i].director[0];
      sumOPTail.director[1] += orderParamsTail_[i].director[1];
      sumOPTail.director[2] += orderParamsTail_[i].director[2];
    }

    avgOPTail.p2 = sumOPTail.p2 / (RealType)orderParamsTail_.size();
    avgOPTail.director[0] =
        sumOPTail.director[0] / (RealType)orderParamsTail_.size();
    avgOPTail.director[1] =
        sumOPTail.director[1] / (RealType)orderParamsTail_.size();
    avgOPTail.director[2] =
        sumOPTail.director[2] / (RealType)orderParamsTail_.size();
    for (std::size_t i = 0; i < orderParamsTail_.size(); ++i) {
      errsumOPTail.p2 += pow((orderParamsTail_[i].p2 - avgOPTail.p2), 2);
    }
    errOPTail.p2 = sqrt(errsumOPTail.p2 / (RealType)orderParamsTail_.size());

    writeP2();
  }

  void RippleOP::writeP2() {
    std::ofstream os(getOutputFileName().c_str());
    os << "#selection1: (" << selectionScript1_ << ")\n";
    os << "#p2\terror\tdirector_x\tdirector_y\tdiretor_z\n";

    os << avgOPHeadUpper.p2 << "\t" << errOPHeadUpper.p2 << "\t"
       << avgOPHeadUpper.director[0] << "\t" << avgOPHeadUpper.director[1]
       << "\t" << avgOPHeadUpper.director[2] << "\n";

    os << avgOPHeadLower.p2 << "\t" << errOPHeadLower.p2 << "\t"
       << avgOPHeadLower.director[0] << "\t" << avgOPHeadLower.director[1]
       << "\t" << avgOPHeadLower.director[2] << "\n";

    os << "selection2: (" << selectionScript2_ << ")\n";
    os << "#p2\terror\tdirector_x\tdirector_y\tdiretor_z\n";

    os << avgOPTail.p2 << "\t" << errOPTail.p2 << "\t" << avgOPTail.director[0]
       << "\t" << avgOPTail.director[1] << "\t" << avgOPTail.director[2]
       << "\n";
  }

}  // namespace OpenMD
