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

#include "applications/dynamicProps/AngularVelVelOutProdCorrFunc.hpp"

#include "math/SquareMatrix3.hpp"

namespace OpenMD {
  AngularVelVelOutProdCorrFunc::AngularVelVelOutProdCorrFunc(
      SimInfo* info, const std::string& filename, const std::string& sele1,
      const std::string& sele2) :
      ObjectCCF<Mat3x3d>(info, filename, sele1, sele2) {
    setCorrFuncType(
        "Angular Velocity - Velocity Outer Product Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".wvOutProdcorr");

    velocity_.resize(nFrames_);
    angularVelocity_.resize(nFrames_);

    sumVelocity_        = V3Zero;
    sumAngularVelocity_ = V3Zero;

    velocityCount_        = 0;
    angularVelocityCount_ = 0;

    propertyTemp = V3Zero;
  }

  int AngularVelVelOutProdCorrFunc::computeProperty1(int frame,
                                                     StuntDouble* sd) {
    if (sd->isDirectional()) {
      Mat3x3d momentInertia = sd->getI();
      Vector3d angMom       = sd->getJ();
      Vector3d omega        = momentInertia.inverse() * angMom;
      propertyTemp          = omega;
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "The selection contains non-directional entities. Your selection should include\
 Directional atoms and/or Rigid Bodies.\n");
      painCave.isFatal = 1;
      simError();
    }

    angularVelocity_[frame].push_back(propertyTemp);
    sumAngularVelocity_ += propertyTemp;
    angularVelocityCount_++;
    return angularVelocity_[frame].size() - 1;
  }

  int AngularVelVelOutProdCorrFunc::computeProperty2(int frame,
                                                     StuntDouble* sd) {
    Vector3d vel = sd->getVel();
    propertyTemp = vel;

    velocity_[frame].push_back(propertyTemp);
    sumVelocity_ += propertyTemp;
    velocityCount_++;
    return velocity_[frame].size() - 1;
  }

  Mat3x3d AngularVelVelOutProdCorrFunc::calcCorrVal(int frame1, int frame2,
                                                    int id1, int id2) {
    Mat3x3d tmpMat_1;
    tmpMat_1 =
        outProduct(angularVelocity_[frame1][id1], velocity_[frame2][id2]);
    Mat3x3d tmpMat_2;
    tmpMat_2 =
        outProduct(angularVelocity_[frame2][id2], velocity_[frame1][id1]);
    Mat3x3d tmpMat_3;
    tmpMat_3 = 0.5 * (tmpMat_1 + tmpMat_2);
    return tmpMat_3;
  }

  void AngularVelVelOutProdCorrFunc::postCorrelate() {
    // Gets the average of the angular velocities
    sumAngularVelocity_ /= RealType(angularVelocityCount_);

    // Gets the average of the velocities
    sumVelocity_ /= RealType(velocityCount_);

    Mat3x3d correlationOfAverages_ =
        outProduct(sumAngularVelocity_, sumVelocity_);

    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
        histogram_[i] /= RealType(count_[i]);

        // The outerProduct correlation of the averages is subtracted
        // from the correlation value:
        histogram_[i] -= correlationOfAverages_;
      } else {
        histogram_[i] = M3Zero;
      }
    }
  }

}  // namespace OpenMD
