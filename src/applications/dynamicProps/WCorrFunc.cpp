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

#include "applications/dynamicProps/WCorrFunc.hpp"

#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"

namespace OpenMD {
  WCorrFunc::WCorrFunc(SimInfo* info, const std::string& filename,
                       const std::string& sele1, const std::string& sele2) :
      ObjectACF<RealType>(info, filename, sele1, sele2,
                          DataStorage::dslFlucQVelocity) {
    setCorrFuncType(" Charge Velocity Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".wcorr");
    setLabelString("<w(0)w(t)>");
    charge_velocities_.resize(nFrames_);
  }

  int WCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    charge_velocities_[frame].push_back(sd->getFlucQVel());
    return charge_velocities_[frame].size() - 1;
  }

  RealType WCorrFunc::calcCorrVal(int frame1, int frame2, int id1, int id2) {
    RealType v2 =
        charge_velocities_[frame1][id1] * charge_velocities_[frame2][id2];
    return v2;
  }

  void WCorrFunc::validateSelection(SelectionManager& seleMan) {
    StuntDouble* sd;
    int i;

    for (sd = seleMan.beginSelected(i); sd != NULL;
         sd = seleMan.nextSelected(i)) {
      Atom* atom                   = static_cast<Atom*>(sd);
      AtomType* atomType           = atom->getAtomType();
      FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);

      if (!fqa.isFluctuatingCharge()) {
        sprintf(painCave.errMsg,
                "WCorrFunc::validateSelection Error: selection "
                "%d (%s)\n"
                "\t is not a fluq object\n",
                sd->getGlobalIndex(), sd->getType().c_str());
        painCave.isFatal = 1;
        simError();
      }
    }
  }

}  // namespace OpenMD
