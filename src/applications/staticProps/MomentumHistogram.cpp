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

/*
 * Calculates the momentum  profile for selected atom.
 * Created by Hemanta Bhattarai on 04/30/19.
 * @author  Hemanta Bhattarai
 */

#include "applications/staticProps/MomentumHistogram.hpp"

#include <algorithm>
#include <fstream>
#include <numeric>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

namespace OpenMD {

  MomentumHistogram::MomentumHistogram(SimInfo* info,
                                       const std::string& filename,
                                       const std::string& sele, int nbins,
                                       int momentum_type,
                                       int momentum_component) :
      StaticAnalyser(info, filename, nbins),
      selectionScript_(sele), evaluator_(info), seleMan_(info), nBins_(nbins),
      mom_type_(momentum_type), mom_comp_(momentum_component) {
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    switch (mom_type_) {
    case 1:
      momentumLabel_ = "Angular Momentum: J";
      break;
    case 0:
    default:
      momentumLabel_ = "Linear Momentum: P";
      break;
    }

    switch (mom_comp_) {
    case 0:
      componentLabel_ = "x";
      break;
    case 1:
      componentLabel_ = "y";
      break;
    case 2:
    default:
      componentLabel_ = "z";
      break;
    }

    setOutputName(getPrefix(filename) + ".MomentumHistogram");
  }

  void MomentumHistogram::process() {
    StuntDouble* sd;
    int ii;

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames / step_;
    vector<RealType> momentum;

    nProcessed_ = nFrames / step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for (sd = seleMan_.beginSelected(ii); sd != NULL;
           sd = seleMan_.nextSelected(ii)) {
        switch (mom_type_) {
        case 0: {
          Vector3d linMom = sd->getVel() * sd->getMass();
          momentum.push_back(linMom[mom_comp_]);
          break;
        }
        case 1:
        default: {
          if (sd->isDirectional()) {
            Vector3d angMom = sd->getJ();
            momentum.push_back(angMom[mom_comp_]);
          }
        } break;
        }
      }
    }

    if (momentum.empty()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Momentum for selected atom not found.\n");
      painCave.isFatal = 1;
      simError();
    }

    std::sort(momentum.begin(), momentum.end());

    RealType min = momentum.front();
    RealType max = momentum.back();

    RealType delta_momentum = (max - min) / (nBins_);

    if (delta_momentum == 0) {
      bincenter_.push_back(min);
      histList_.push_back(momentum.size());
    } else {
      // fill the center for histogram
      for (int j = 0; j < nBins_ + 3; ++j) {
        bincenter_.push_back(min + (j - 1) * delta_momentum);
        histList_.push_back(0);
      }

      // filling up the histogram whith the densities
      int bin_center_pos = 0;
      vector<RealType>::iterator index;
      RealType momentum_length = static_cast<RealType>(momentum.size());

      bool hist_update;
      for (index = momentum.begin(); index < momentum.end(); index++) {
        hist_update = true;
        while (hist_update) {
          if (*index >= bincenter_[bin_center_pos] &&
              *index < bincenter_[bin_center_pos + 1]) {
            histList_[bin_center_pos] +=
                1.0 / (momentum_length * delta_momentum);
            hist_update = false;
          } else {
            bin_center_pos++;
            hist_update = true;
          }
        }
      }
    }
    write();
  }

  void MomentumHistogram::write() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#" << momentumLabel_ << componentLabel_
                << " distribtution \n";
      rdfStream << "# nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "# selection: (" << selectionScript_ << ")\n";
      rdfStream << "# "
                << "Bin_center"
                << "\tcount\n";
      for (unsigned int i = 0; i < histList_.size(); ++i) {
        rdfStream << bincenter_[i] << "\t" << histList_[i] << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "MomentumHistogram: unable to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    rdfStream.close();
  }
}  // namespace OpenMD
