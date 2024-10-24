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

#include "applications/dynamicProps/DirectionalRCorrFunc.hpp"

#include "utils/Revision.hpp"

namespace OpenMD {
  DirectionalRCorrFunc::DirectionalRCorrFunc(SimInfo* info,
                                             const std::string& filename,
                                             const std::string& sele1,
                                             const std::string& sele2) :
      ObjectACF<Vector3d>(info, filename, sele1, sele2) {
    setCorrFuncType("DirectionalRCorrFunc");
    setOutputName(getPrefix(dumpFilename_) + ".drcorr");
    setLabelString("r2\trparallel\trperpendicular");
    positions_.resize(nFrames_);
    rotMats_.resize(nFrames_);
  }

  int DirectionalRCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    positions_[frame].push_back(sd->getPos());
    rotMats_[frame].push_back(sd->getA());
    return positions_[frame].size() - 1;
  }

  Vector3d DirectionalRCorrFunc::calcCorrVal(int frame1, int frame2, int id1,
                                             int id2) {
    Vector3d diff = positions_[frame2][id2] - positions_[frame1][id1];

    // The lab frame vector corresponding to the body-fixed
    // z-axis is simply the second column of A.transpose()
    // or, identically, the second row of A itself.

    Vector3d u1  = rotMats_[frame1][id1].getRow(2);
    RealType u1l = u1.length();

    RealType rsq    = diff.lengthSquare();
    RealType rpar   = dot(diff, u1) / u1l;
    RealType rpar2  = rpar * rpar;
    RealType rperp2 = rsq - rpar2;

    return Vector3d(rsq, rpar2, rperp2);
  }

  void DirectionalRCorrFunc::validateSelection(SelectionManager&) {
    StuntDouble* sd;
    int i;
    for (sd = seleMan1_.beginSelected(i); sd != NULL;
         sd = seleMan1_.nextSelected(i)) {
      if (!sd->isDirectional()) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "DirectionalRCorrFunc::validateSelection Error: "
                 "at least one of the selected objects is not Directional\n");
        painCave.isFatal = 1;
        simError();
      }
    }
  }
}  // namespace OpenMD
