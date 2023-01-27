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
      ObjectACF<RealType>(info, filename, sele1, sele2) {
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
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
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
