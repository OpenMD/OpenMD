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

#include "applications/dynamicProps/ForceAutoCorrFunc.hpp"

#include "math/SquareMatrix3.hpp"
#include "utils/Revision.hpp"

namespace OpenMD {
  ForceAutoCorrFunc::ForceAutoCorrFunc(SimInfo* info,
                                       const std::string& filename,
                                       const std::string& sele1,
                                       const std::string& sele2) :
      ObjectACF<Mat3x3d>(info, filename, sele1, sele2) {
    setCorrFuncType("Force - Force Auto Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".facorr");

    forces_.resize(nFrames_);
    sumForces_   = V3Zero;
    forcesCount_ = 0;
  }

  int ForceAutoCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    if (sd->isDirectional()) {
      Mat3x3d A    = sd->getA();
      Vector3d f   = sd->getFrc();
      propertyTemp = A * f;
    } else {
      Vector3d f   = sd->getFrc();
      propertyTemp = f;
    }
    forces_[frame].push_back(propertyTemp);
    sumForces_ += propertyTemp;
    forcesCount_++;
    return forces_[frame].size() - 1;
  }

  Mat3x3d ForceAutoCorrFunc::calcCorrVal(int frame1, int frame2, int id1,
                                         int id2) {
    return outProduct(forces_[frame1][id1], forces_[frame2][id2]);
  }

  void ForceAutoCorrFunc::postCorrelate() {
    // Gets the average of the forces_
    sumForces_ /= RealType(forcesCount_);

    Mat3x3d correlationOfAverages_ = outProduct(sumForces_, sumForces_);
    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
        histogram_[i] /= RealType(count_[i]);

        // The outerProduct correlation of the averages is subtracted
        // from the correlation value
        histogram_[i] -= correlationOfAverages_;
      } else {
        histogram_[i] = M3Zero;
      }
    }
  }
}  // namespace OpenMD
