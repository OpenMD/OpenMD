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

#include "applications/dynamicProps/ChargeOrientationCorrFunc.hpp"

#include "math/SquareMatrix3.hpp"
#include "primitives/Molecule.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"

namespace OpenMD {
  ChargeOrientationCorrFunc::ChargeOrientationCorrFunc(
      SimInfo* info, const std::string& filename, const std::string& sele1,
      const std::string& sele2, const RealType dipoleX, const RealType dipoleY,
      const RealType dipoleZ, const RealType cutOff, const int axis) :
      ObjectCCF<RealType>(info, filename, sele1, sele2),
      axis_(axis) {
    setCorrFuncType(
        "Charge - Orientation Order Parameter Cross Correlation Function");
    setOutputName(getPrefix(dumpFilename_) + ".QScorr");

    charges_.resize(nFrames_);
    CosTheta_.resize(nFrames_);

    sumCharge_     = 0;
    sumCosTheta_   = 0;
    chargeCount_   = 0;
    CosThetaCount_ = 0;

    dipoleVector_ = Vector3d(dipoleX, dipoleY, dipoleZ);
    dipoleVector_.normalize();

    switch (axis_) {
    case 0:
      axisLabel_ = "x";
      refAxis_   = Vector3d(1, 0, 0);
      break;
    case 1:
      axisLabel_ = "y";
      refAxis_   = Vector3d(0, 1, 0);
      break;
    case 2:
    default:
      axisLabel_ = "z";
      refAxis_   = Vector3d(0, 0, 1);
      break;
    }
  }

  void ChargeOrientationCorrFunc::validateSelection(SelectionManager& seleMan) {
    StuntDouble* sd;
    int i;

    for (sd = seleMan.beginSelected(i); sd != NULL;
         sd = seleMan.nextSelected(i)) {
      Atom* atom                   = static_cast<Atom*>(sd);
      AtomType* atomType           = atom->getAtomType();
      FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);

      if (!sd->isDirectional() && !fqa.isFluctuatingCharge()) {
        snprintf(
            painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
            "ChargeOrientationCorrFunc::validateSelection Error: selection "
            "%d (%s)\n"
            "\t is not a Directional object\n",
            sd->getGlobalIndex(), sd->getType().c_str());
        painCave.isFatal = 1;
        simError();
      }
    }
  }

  int ChargeOrientationCorrFunc::computeProperty1(int frame, StuntDouble* sd) {
    RealType q = 0.0;
    Atom* atom = static_cast<Atom*>(sd);

    AtomType* atomType = atom->getAtomType();

    FixedChargeAdapter fca = FixedChargeAdapter(atomType);
    if (fca.isFixedCharge()) { q += fca.getCharge(); }

    FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
    if (fqa.isFluctuatingCharge()) { q += atom->getFlucQPos(); }

    propertyTemp = q;
    charges_[frame].push_back(propertyTemp);
    sumCharge_ += propertyTemp;
    chargeCount_++;
    return charges_[frame].size() - 1;
  }

  int ChargeOrientationCorrFunc::computeProperty2(int frame, StuntDouble* sd) {
    SquareMatrix3<RealType> rotMat;
    Vector3d rotatedDipoleVector;
    RealType ctheta(0.0);

    rotMat              = sd->getA();
    rotatedDipoleVector = rotMat * dipoleVector_;
    rotatedDipoleVector.normalize();
    ctheta = dot(rotatedDipoleVector, refAxis_);

    propertyTemp = ctheta;
    CosTheta_[frame].push_back(propertyTemp);
    sumCosTheta_ += propertyTemp;
    CosThetaCount_++;
    return CosTheta_[frame].size() - 1;
  }

  RealType ChargeOrientationCorrFunc::calcCorrVal(int frame1, int frame2,
                                                  int id1, int id2) {
    return charges_[frame1][id1] * CosTheta_[frame2][id2];
  }

  void ChargeOrientationCorrFunc::postCorrelate() {
    // gets the average of the charges
    sumCharge_ /= RealType(chargeCount_);

    // gets the average of the CosTheta
    sumCosTheta_ /= RealType(CosThetaCount_);

    RealType correlationOfAverages_ = sumCharge_ * sumCosTheta_;
    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
        histogram_[i] /= RealType(count_[i]);

        // The  correlation of the averages is subtracted
        // from the correlation value:
        histogram_[i] -= correlationOfAverages_;
      } else {
        histogram_[i] = 0;
      }
    }
  }
}  // namespace OpenMD
