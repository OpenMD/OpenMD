/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
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
      sprintf(painCave.errMsg, "RhoR: unable to open %s\n",
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }

}  // namespace OpenMD
