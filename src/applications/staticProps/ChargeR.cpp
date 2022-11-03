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

#include "applications/staticProps/ChargeR.hpp"

#include <algorithm>
#include <fstream>

#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {
  
  ChargeR::ChargeR(SimInfo* info, const std::string& filename,
                   const std::string& sele, RealType len, int nrbins) :
    StaticAnalyser(info, filename, nrbins),
    selectionScript_(sele), evaluator_(info), seleMan_(info), thermo_(info),
    len_(len) {
    
    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    deltaR_ = len_ / nBins_;
    
    // fixed number of bins

    sliceSDLists_.resize(nBins_);
    sliceSDCount_.resize(nBins_);
    std::fill(sliceSDCount_.begin(), sliceSDCount_.end(), 0);

    chargeR_.resize(nBins_);
    setOutputName(getPrefix(filename) + ".ChargeR");
    std::stringstream params;
    params << " len = " << len_ << ", nrbins = " << nBins_;
    const std::string paramString = params.str();
    setParameterString(paramString);

  }

  void ChargeR::process() {
    StuntDouble* sd;
    int ii;

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames / step_;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      for (unsigned int i = 0; i < nBins_; i++) {
        sliceSDLists_[i].clear();
      }

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      // determine which atom belongs to which slice
      for (sd = seleMan_.beginSelected(ii); sd != NULL;
           sd = seleMan_.nextSelected(ii)) {
        Vector3d pos = sd->getPos();
        RealType distance = pos.length();
        
        if (distance < len_) {
          int binNo = int(distance / deltaR_);
          sliceSDLists_[binNo].push_back(sd);
          sliceSDCount_[binNo]++;
        }
      }

      // loop over the slices to calculate the charge
      for (unsigned int i = 0; i < nBins_; i++) {
        RealType binC = 0;
        for (unsigned int k = 0; k < sliceSDLists_[i].size(); ++k) {
          RealType q = 0.0;
          Atom* atom = static_cast<Atom*>(sliceSDLists_[i][k]);
          
          AtomType* atomType = atom->getAtomType();

          if (sliceSDLists_[i][k]->isAtom()) {
            FixedChargeAdapter fca = FixedChargeAdapter(atomType);
            if (fca.isFixedCharge()) { q += fca.getCharge(); }

            FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
            if (fqa.isFluctuatingCharge()) { q += atom->getFlucQPos(); }
          }

          binC += q;
        }
        chargeR_[i] += binC;
        // Units of (e / Ang^2 / fs)
      }
    }

    writeChargeR();
  }

  void ChargeR::writeChargeR() {
    
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#ChargeR "
                << "\n";
      rdfStream << "#selection: (" << selectionScript_ << ")\n";
      rdfStream << "# r " << "\tcharge\n";
      RealType binCharge;
      for (unsigned int i = 0; i < chargeR_.size(); ++i) {
        RealType rLower = i * deltaR_;
        RealType rUpper = rLower + deltaR_;
        RealType volShell = (4.0 * Constants::PI) *
          (pow(rUpper, 3) - pow(rLower, 3)) / 3.0;        
        
        RealType r = deltaR_ * (i + 0.5);

        binCharge = chargeR_[i] / (volShell * nProcessed_);

        rdfStream << r << "\t" << binCharge << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH, "ChargeR: unable to open %s\n",
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }
}  // namespace OpenMD
