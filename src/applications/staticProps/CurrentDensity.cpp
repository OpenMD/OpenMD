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

/*
 * Computes the current density for the selected atom
 * Created by Cody R. Drisko on 06/14/19.
 */

#include "applications/staticProps/CurrentDensity.hpp"

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {

  CurrentDensity::CurrentDensity(SimInfo* info, const std::string& filename,
                                 const std::string& sele, int nbins, int axis) :
      StaticAnalyser(info, filename, nbins),
      selectionScript_(sele), evaluator_(info), seleMan_(info), thermo_(info),
      axis_(axis) {
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // fixed number of bins
    sliceSDLists_.resize(nBins_);
    currentDensity_.resize(nBins_);

    switch (axis_) {
    case 0:
      axisLabel_ = "x";
      break;
    case 1:
      axisLabel_ = "y";
      break;
    case 2:
    default:
      axisLabel_ = "z";
      break;
    }

    setOutputName(getPrefix(filename) + ".Jc");
  }

  void CurrentDensity::process() {
    StuntDouble* sd;
    int ii;

    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames            = reader.getNFrames();
    nProcessed_            = nFrames / step_;
    overallCurrentDensity_ = 0;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
      Vector3d COMvel  = thermo_.getComVel();

      for (unsigned int i = 0; i < nBins_; i++) {
        sliceSDLists_[i].clear();
      }

      RealType boxVolume   = currentSnapshot_->getVolume();
      RealType sliceVolume = boxVolume / nBins_;
      Mat3x3d hmat         = currentSnapshot_->getHmat();
      zBox_.push_back(hmat(axis_, axis_));

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      // determine which atom belongs to which slice
      for (sd = seleMan_.beginSelected(ii); sd != NULL;
           sd = seleMan_.nextSelected(ii)) {
        int binNo;
        Vector3d pos = sd->getPos();

        if (usePeriodicBoundaryConditions_) {
          currentSnapshot_->wrapVector(pos);
          binNo =
              int(nBins_ * (pos[axis_] / hmat(axis_, axis_) + 0.5)) % nBins_;
          sliceSDLists_[binNo].push_back(sd);
        }
      }

      // loop over the slices to calculate the densities
      for (unsigned int i = 0; i < nBins_; i++) {
        RealType binJc = 0;

        for (unsigned int j = 0; j < sliceSDLists_[i].size(); ++j) {
          RealType q = 0.0;
          Atom* atom = static_cast<Atom*>(sliceSDLists_[i][j]);

          AtomType* atomType = atom->getAtomType();

          if (sliceSDLists_[i][j]->isAtom()) {
            FixedChargeAdapter fca = FixedChargeAdapter(atomType);
            if (fca.isFixedCharge()) q = fca.getCharge();

            FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
            if (fqa.isFluctuatingCharge()) q += atom->getFlucQPos();

            Vector3d vel = sliceSDLists_[i][j]->getVel();
            binJc += q * (vel[axis_] - COMvel[axis_]);
          }
        }

        // Units of (e / Ang^2 / fs)
        currentDensity_[i] += binJc / sliceVolume;
        overallCurrentDensity_ += binJc / boxVolume;
      }
    }

    writeCurrentDensity();
  }

  void CurrentDensity::writeCurrentDensity() {
    // compute average box length:
    std::vector<RealType>::iterator j;
    RealType zSum = 0.0;
    for (j = zBox_.begin(); j != zBox_.end(); ++j) {
      zSum += *j;
    }
    RealType zAve = zSum / zBox_.size();

    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#Current Density = "
                << overallCurrentDensity_ / nProcessed_
                << " e / Ang^2 / fs.\n";
      rdfStream << "#J_c(" << axisLabel_ << ")\n";
      rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "#" << axisLabel_ << "\tcurrent density\n";

      for (unsigned int i = 0; i < currentDensity_.size(); ++i) {
        RealType z = zAve * (i + 0.5) / currentDensity_.size();
        rdfStream << z << "\t" << currentDensity_[i] / nProcessed_ << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "CurrentDensity: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }
}  // namespace OpenMD
