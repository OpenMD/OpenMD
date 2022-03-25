/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#include "applications/staticProps/SCDOrderParameter.hpp"

#include <string>
#include <tuple>
#include <vector>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

namespace OpenMD {

  SCDElem::SCDElem(SimInfo* info, const std::string& sele1,
                   const std::string& sele2, const std::string& sele3) :
      sele1_(sele1),
      sele2_(sele2), sele3_(sele3) {
    usePeriodicBoundaryConditions_ =
        info->getSimParams()->getUsePeriodicBoundaryConditions();

    SelectionManager seleMan1_(info);
    SelectionManager seleMan2_(info);
    SelectionManager seleMan3_(info);
    SelectionEvaluator evaluator1_(info);
    SelectionEvaluator evaluator2_(info);
    SelectionEvaluator evaluator3_(info);

    evaluator1_.loadScriptString(sele1_);
    evaluator2_.loadScriptString(sele2_);
    evaluator3_.loadScriptString(sele3_);

    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    } else {
      sprintf(painCave.errMsg, "dynamic selection is not allowed\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    } else {
      sprintf(painCave.errMsg, "dynamic selection is not allowed\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    if (!evaluator3_.isDynamic()) {
      seleMan3_.setSelectionSet(evaluator3_.evaluate());
    } else {
      sprintf(painCave.errMsg, "dynamic selection is not allowed\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    int nselected1 = seleMan1_.getSelectionCount();
    int nselected2 = seleMan2_.getSelectionCount();
    int nselected3 = seleMan3_.getSelectionCount();

    if (nselected1 != nselected2 || nselected1 != nselected3) {
      sprintf(painCave.errMsg,
              "The number of selected Stuntdoubles must be the same\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    int i;
    int j;
    int k;
    StuntDouble* sd1;
    StuntDouble* sd2;
    StuntDouble* sd3;
    for (sd1 = seleMan1_.beginSelected(i), sd2 = seleMan2_.beginSelected(j),
        sd3 = seleMan3_.beginSelected(k);
         sd1 != NULL && sd2 != NULL && sd3 != NULL;
         sd1 = seleMan1_.nextSelected(i), sd2 = seleMan2_.nextSelected(j),
        sd3 = seleMan3_.nextSelected(k)) {
      tuples_.push_back(std::make_tuple(sd1, sd2, sd3));
    }
  }

  RealType SCDElem::calcSCD(Snapshot* snapshot) {
    Vector3d normal(0.0, 0.0, 1.0);
    RealType scd = 0.0;

    for (auto& [first, second, third] : tuples_) {
      // Egberts B. and Berendsen H.J.C, J.Chem.Phys. 89(6), 3718-3732, 1988
      Vector3d zAxis = third->getPos() - first->getPos();
      if (usePeriodicBoundaryConditions_) snapshot->wrapVector(zAxis);
      Vector3d v12 = second->getPos() - first->getPos();
      if (usePeriodicBoundaryConditions_) snapshot->wrapVector(v12);
      Vector3d xAxis = cross(v12, zAxis);
      Vector3d yAxis = cross(zAxis, xAxis);

      xAxis.normalize();
      yAxis.normalize();
      zAxis.normalize();
      RealType cosThetaX = dot(xAxis, normal);
      RealType sxx       = 0.5 * (3 * cosThetaX * cosThetaX - 1.0);
      RealType cosThetaY = dot(yAxis, normal);
      RealType syy       = 0.5 * (3 * cosThetaY * cosThetaY - 1.0);
      scd += 2.0 / 3.0 * sxx + 1.0 / 3.0 * syy;
    }

    scd /= tuples_.size();

    return scd;
  }

  SCDOrderParameter::SCDOrderParameter(SimInfo* info,
                                       const std::string& filename,
                                       const std::string& sele1,
                                       const std::string& sele2,
                                       const std::string& sele3) :
      StaticAnalyser(info, filename, 1) {
    setOutputName(getPrefix(filename) + ".scd");

    scdElems_.push_back(SCDElem(info, sele1, sele2, sele3));
    scdParam_.resize(scdElems_.size());
    std::fill(scdParam_.begin(), scdParam_.end(), 0.0);
  }

  SCDOrderParameter::SCDOrderParameter(SimInfo* info,
                                       const std::string& filename,
                                       const std::string& molname,
                                       int beginIndex, int endIndex) :
      StaticAnalyser(info, filename, 1) {
    setOutputName(getPrefix(filename) + ".scd");

    assert(beginIndex >= 0 && endIndex >= 0 && beginIndex <= endIndex - 2);
    for (int i = beginIndex; i <= endIndex - 2; ++i) {
      std::string selectionTemplate = "select " + molname + ".";
      std::string sele1             = selectionTemplate + OpenMD_itoa(i);
      std::string sele2             = selectionTemplate + OpenMD_itoa(i + 1);
      std::string sele3             = selectionTemplate + OpenMD_itoa(i + 2);

      scdElems_.push_back(SCDElem(info, sele1, sele2, sele3));
    }

    scdParam_.resize(scdElems_.size());
    std::fill(scdParam_.begin(), scdParam_.end(), 0.0);
  }

  void SCDOrderParameter::process() {
    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();

    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for (std::size_t j = 0; j < scdElems_.size(); ++j) {
        scdParam_[j] += scdElems_[j].calcSCD(currentSnapshot_);
      }
    }

    int nProcessed = nFrames / step_;
    for (std::size_t j = 0; j < scdElems_.size(); ++j) {
      scdParam_[j] /= nProcessed;
    }

    writeSCD();
  }

  void SCDOrderParameter::writeSCD() {
    std::ofstream os(getOutputFileName().c_str());
    os << "#scd parameter\n";
    for (std::size_t i = 0; i < scdElems_.size(); ++i) {
      os << "#[column " << i + 1 << "]\t"
         << "sele1: \"" << scdElems_[i].getSelection1() << "\",\t"
         << "sele2: \"" << scdElems_[i].getSelection2() << "\",\t"
         << "sele3: \"" << scdElems_[i].getSelection3() << "\"\n";
    }

    for (std::size_t i = 0; i < scdElems_.size(); ++i) {
      os << scdParam_[i] << "\t";
    }
    os << std::endl;
  }

}  // namespace OpenMD
