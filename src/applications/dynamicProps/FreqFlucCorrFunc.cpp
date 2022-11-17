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

#include "applications/dynamicProps/FreqFlucCorrFunc.hpp"

#include <string>

#include "applications/dynamicProps/TimeCorrFunc.hpp"
#include "brains/DataStorage.hpp"
#include "brains/SimInfo.hpp"
#include "math/Vector3.hpp"
#include "primitives/Atom.hpp"
#include "primitives/StuntDouble.hpp"
#include "types/MultipoleAdapter.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {

  FreqFlucCorrFunc::FreqFlucCorrFunc(SimInfo* info, const std::string& filename,
                                     const std::string& sele1,
                                     const std::string& sele2) :
      ObjectACF<RealType>(info, filename, sele1, sele2,
                          DataStorage::dslElectricField | DataStorage::dslAmat |
                              DataStorage::dslDipole) {
    setCorrFuncType("FreqFluc Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".ffcorr");
  }

  int FreqFlucCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    Vector3d e = sd->getElectricField();
    Vector3d u;
    if (sd->isRigidBody()) {
      u = sd->getA().getRow(2);
    } else {
      AtomType* at        = static_cast<Atom*>(sd)->getAtomType();
      MultipoleAdapter ma = MultipoleAdapter(at);

      if (ma.isDipole()) { u = sd->getDipole(); }
    }

    RealType uedot = dot(u, e);
    ueStats_.add(uedot);

    ue_[frame].push_back(uedot);
    return ue_[frame].size() - 1;
  }

  RealType FreqFlucCorrFunc::calcCorrVal(int frame1, int frame2, int id1,
                                         int id2) {
    RealType mean = ueStats_.getAverage();

    return (ue_[frame1][id1] - mean) * (ue_[frame2][id2] - mean);
  }

  void FreqFlucCorrFunc::validateSelection(SelectionManager&) {
    StuntDouble* sd;
    int i;
    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {
      if (!sd->isRigidBody()) {
        AtomType* at        = static_cast<Atom*>(sd)->getAtomType();
        MultipoleAdapter ma = MultipoleAdapter(at);

        if (!ma.isDipole()) {
          snprintf(
              painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
              "FreqFlucCorrFunc::validateSelection Error: selection is not a "
              "RigidBody "
              "or does\n"
              "\tnot have a dipole\n");
          painCave.isFatal = 1;
          simError();
        }
      }
    }
  }
}  // namespace OpenMD
