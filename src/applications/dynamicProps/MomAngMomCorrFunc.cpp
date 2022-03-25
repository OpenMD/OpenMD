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

#include "applications/dynamicProps/MomAngMomCorrFunc.hpp"

namespace OpenMD {
  MomAngMomCorrFunc::MomAngMomCorrFunc(SimInfo* info,
                                       const std::string& filename,
                                       const std::string& sele1,
                                       const std::string& sele2) :
      ObjectCCF<RealType>(info, filename, sele1, sele2,
                          DataStorage::dslVelocity | DataStorage::dslAmat |
                              DataStorage::dslAngularMomentum) {
    setCorrFuncType("Momentum - Angular Momentum Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".pjcorr");

    momenta_.resize(nFrames_);
    js_.resize(nFrames_);
  }

  int MomAngMomCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    momenta_[frame].push_back(sd->getMass() * sd->getVel());
    return momenta_[frame].size() - 1;
  }

  int MomAngMomCorrFunc::computeProperty2(int frame, StuntDouble* sd) {
    // The lab frame vector corresponding to the body-fixed
    // z-axis is simply the second column of A.transpose()
    // or, identically, the second row of A itself.
    Vector3d u = sd->getA().getRow(2);
    Vector3d j = sd->getJ();

    js_[frame].push_back(cross(j, u));
    return js_[frame].size() - 1;
  }

  RealType MomAngMomCorrFunc::calcCorrVal(int frame1, int frame2, int id1,
                                          int id2) {
    RealType pj = dot(momenta_[frame1][id1], js_[frame2][id2]);
    return pj;
  }

  void MomAngMomCorrFunc::validateSelection(SelectionManager& seleMan) {
    StuntDouble* sd;
    int i;

    for (sd = seleMan.beginSelected(i); sd != NULL;
         sd = seleMan.nextSelected(i)) {
      if (!sd->isDirectional()) {
        sprintf(painCave.errMsg,
                "MomAngMomCorrFunc::validateSelection Error: selection "
                "%d (%s)\n"
                "\t is not a Directional object\n",
                sd->getGlobalIndex(), sd->getType().c_str());
        painCave.isFatal = 1;
        simError();
      }
    }
  }

}  // namespace OpenMD
