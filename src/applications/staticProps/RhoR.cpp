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

/* Calculates Rho(R) for nanoparticle with radius R*/
#include "applications/staticProps/RhoR.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>

#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  RhoR::RhoR(SimInfo* info, const std::string& filename,
             const std::string& sele, RealType len, int nrbins,
             RealType particleR) :
      StaticAnalyser(info, filename, nrbins),
      selectionScript_(sele), evaluator_(info), seleMan_(info), len_(len) {
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    deltaR_ = len_ / nBins_;

    histogram_.resize(nBins_);
    avgRhoR_.resize(nBins_);
    particleR_ = particleR;
    setOutputName(getPrefix(filename) + ".RhoR");
  }

  void RhoR::process() {
    Thermo thermo(info_);
    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames / step_;

    std::fill(avgRhoR_.begin(), avgRhoR_.end(), 0.0);
    std::fill(histogram_.begin(), histogram_.end(), 0);

    for (int istep = 0; istep < nFrames; istep += step_) {
      int i;
      StuntDouble* sd;
      reader.readFrame(istep);
      currentSnapshot_      = info_->getSnapshotManager()->getCurrentSnapshot();
      Vector3d CenterOfMass = thermo.getCom();

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      // determine which atom belongs to which slice
      for (sd = seleMan_.beginSelected(i); sd != NULL;
           sd = seleMan_.nextSelected(i)) {
        Vector3d pos = sd->getPos();
        Vector3d r12 = CenterOfMass - pos;

        RealType distance = r12.length();

        if (distance < len_) {
          int whichBin = int(distance / deltaR_);
          histogram_[whichBin] += 1;
        }
      }
    }

    processHistogram();
    writeRhoR();
  }

  void RhoR::processHistogram() {
    RealType particleDensity = 3.0 * info_->getNGlobalMolecules() /
                               (4.0 * Constants::PI * pow(particleR_, 3));
    RealType pairConstant = (4.0 * Constants::PI * particleDensity) / 3.0;

    for (unsigned int i = 0; i < histogram_.size(); ++i) {
      RealType rLower = i * deltaR_;
      RealType rUpper = rLower + deltaR_;
      RealType volSlice =
          (rUpper * rUpper * rUpper) - (rLower * rLower * rLower);
      RealType nIdeal = volSlice * pairConstant;

      avgRhoR_[i] += histogram_[i] / nIdeal;
    }
  }

  void RhoR::writeRhoR() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#radial density function rho(r)\n";
      rdfStream << "#r\tcorrValue\n";
      for (unsigned int i = 0; i < avgRhoR_.size(); ++i) {
        RealType r = deltaR_ * (i + 0.5);
        rdfStream << r << "\t" << avgRhoR_[i] / nProcessed_ << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "RhoR: unable to open %s\n", outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }

}  // namespace OpenMD
