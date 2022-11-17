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

#include "applications/staticProps/DensityPlot.hpp"

#include <algorithm>
#include <functional>
#include <memory>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "types/LennardJonesAdapter.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

using namespace std;

namespace OpenMD {

  DensityPlot::DensityPlot(SimInfo* info, const std::string& filename,
                           const std::string& sele, const std::string& cmSele,
                           RealType len, int nrbins) :
      StaticAnalyser(info, filename, nrbins),
      len_(len), halfLen_(len / 2), nRBins_(nrbins), selectionScript_(sele),
      seleMan_(info), evaluator_(info), cmSelectionScript_(cmSele),
      cmSeleMan_(info), cmEvaluator_(info) {
    setOutputName(getPrefix(filename) + ".density");

    deltaR_ = len_ / nRBins_;
    histogram_.resize(nRBins_);
    density_.resize(nRBins_);

    std::fill(histogram_.begin(), histogram_.end(), 0);

    evaluator_.loadScriptString(sele);

    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    cmEvaluator_.loadScriptString(cmSele);
    if (!cmEvaluator_.isDynamic()) {
      cmSeleMan_.setSelectionSet(cmEvaluator_.evaluate());
    }
  }

  void DensityPlot::process() {
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      if (cmEvaluator_.isDynamic()) {
        cmSeleMan_.setSelectionSet(cmEvaluator_.evaluate());
      }

      Vector3d origin = calcNewOrigin();

      Mat3x3d hmat        = currentSnapshot_->getHmat();
      RealType slabVolume = deltaR_ * hmat(0, 0) * hmat(1, 1);
      int k;
      for (StuntDouble* sd = seleMan_.beginSelected(k); sd != NULL;
           sd              = seleMan_.nextSelected(k)) {
        if (!sd->isAtom()) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Can not calculate electron density if it is not atom\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        Atom* atom = static_cast<Atom*>(sd);
        std::shared_ptr<GenericData> data =
            atom->getAtomType()->getPropertyByName("nelectron");
        if (data == nullptr) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Can not find Parameters for nelectron\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        std::shared_ptr<DoubleGenericData> doubleData =
            std::dynamic_pointer_cast<DoubleGenericData>(data);
        if (doubleData == nullptr) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Can not cast GenericData to DoubleGenericData\n");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        RealType nelectron      = doubleData->getData();
        LennardJonesAdapter lja = LennardJonesAdapter(atom->getAtomType());
        RealType sigma          = lja.getSigma() * 0.5;
        RealType sigma2         = sigma * sigma;

        Vector3d pos = sd->getPos() - origin;
        for (int j = 0; j < nRBins_; ++j) {
          Vector3d tmp(pos);
          RealType zdist = j * deltaR_ - halfLen_;
          tmp[2] += zdist;
          if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(tmp);

          RealType wrappedZdist = tmp.z() + halfLen_;
          if (wrappedZdist < 0.0 || wrappedZdist > len_) { continue; }

          int which = int(wrappedZdist / deltaR_);
          density_[which] +=
              nelectron * exp(-zdist * zdist / (sigma2 * 2.0)) /
              (slabVolume * sqrt(2 * Constants::PI * sigma * sigma));
        }
      }
    }

    int nProcessed = nFrames / step_;
    std::transform(
        density_.begin(), density_.end(), density_.begin(),
        std::bind(std::divides<RealType>(), std::placeholders::_1, nProcessed));
    writeDensity();
  }

  Vector3d DensityPlot::calcNewOrigin() {
    int i;
    Vector3d newOrigin(0.0);
    RealType totalMass = 0.0;
    for (StuntDouble* sd = seleMan_.beginSelected(i); sd != NULL;
         sd              = seleMan_.nextSelected(i)) {
      RealType mass = sd->getMass();
      totalMass += mass;
      newOrigin += sd->getPos() * mass;
    }
    newOrigin /= totalMass;
    return newOrigin;
  }

  void DensityPlot::writeDensity() {
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);
    if (ofs.is_open()) {
      ofs << "#g(x, y, z)\n";
      ofs << "#selection: (" << selectionScript_ << ")\n";
      ofs << "#cmSelection:(" << cmSelectionScript_ << ")\n";
      ofs << "#nRBins = " << nRBins_ << "\t maxLen = " << len_
          << "\tdeltaR = " << deltaR_ << "\n";
      for (unsigned int i = 0; i < histogram_.size(); ++i) {
        ofs << i * deltaR_ - halfLen_ << "\t" << density_[i] << std::endl;
      }
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "DensityPlot: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }

}  // namespace OpenMD
