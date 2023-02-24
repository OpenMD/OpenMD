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

#include "hydrodynamics/HydroProp.hpp"

#include "math/CholeskyDecomposition.hpp"
#include "utils/StringTokenizer.hpp"
#include "utils/simError.h"

namespace OpenMD {

  HydroProp::HydroProp() : hasCOR_(false), hasXi_(false), hasS_(false) {}

  HydroProp::HydroProp(Vector3d cor, Mat6x6d Xi) :
      cor_(cor), Xi_(Xi), hasCOR_(true), hasXi_(true), hasS_(false) {}

  void HydroProp::complete() {
    if (hasXi_) {
      CholeskyDecomposition(Xi_, S_);
      hasS_ = true;
    } else {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "HydroProp was asked to complete without a Resistance Tensor.\n");
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }
  }

  Mat6x6d HydroProp::getS() {
    if (!hasS_) { complete(); }
    return S_;
  }

  Mat3x3d HydroProp::getXitt() {
    Mat3x3d Xitt;
    Xi_.getSubMatrix(0, 0, Xitt);
    return Xitt;
  }
  Mat3x3d HydroProp::getXirt() {
    Mat3x3d Xirt;
    Xi_.getSubMatrix(0, 3, Xirt);
    return Xirt;
  }
  Mat3x3d HydroProp::getXitr() {
    Mat3x3d Xitr;
    Xi_.getSubMatrix(3, 0, Xitr);
    return Xitr;
  }
  Mat3x3d HydroProp::getXirr() {
    Mat3x3d Xirr;
    Xi_.getSubMatrix(3, 3, Xirr);
    return Xirr;
  }

  Mat6x6d HydroProp::getDiffusionTensor(RealType temperature) {
    Mat6x6d XiCopy = Xi_;
    Mat6x6d D;
    invertMatrix(XiCopy, D);
    RealType kt = Constants::kb * temperature;  // in kcal mol^-1
    D *= kt;  // now in angstroms^2 fs^-1  (at least for Trans-trans)
    return D;
  }

  Mat6x6d HydroProp::getResistanceTensorAtPos(Vector3d pos) {
    // Vector from reference point to center of resistance  = cor_
    // Vector from reference point to new location = pos
    // Vector from center of resistance to new location = pos - cor_
    Vector3d cp = pos - cor_;
    Mat3x3d U;
    U.setupSkewMat(cp);

    Mat3x3d Xitt = getXitt();
    Mat3x3d Xitr = getXitr();
    Mat3x3d Xirr = getXirr();

    Mat3x3d Xipostt;
    Mat3x3d Xiposrr;
    Mat3x3d Xipostr;

    // Resistance tensors at the new location
    Xipostt = Xitt;
    Xipostr = (Xitr - U * Xitt);
    Xiposrr = Xirr - U * Xitt * U + Xitr * U - U * Xitr.transpose();

    Mat6x6d Xipos;
    Xipos.setSubMatrix(0, 0, Xipostt);
    Xipos.setSubMatrix(0, 3, Xipostr.transpose());
    Xipos.setSubMatrix(3, 0, Xipostr);
    Xipos.setSubMatrix(3, 3, Xiposrr);
    return Xipos;
  }

  Mat6x6d HydroProp::getDiffusionTensorAtPos(Vector3d pos,
                                             RealType temperature) {
    // Vector from reference point to center of resistance  = cor_
    // Vector from reference point to new location = pos
    // Vector from center of resistance to new location = pos - cor_

    Vector3d cp = pos - cor_;
    Mat3x3d U;
    U.setupSkewMat(cp);

    Mat6x6d D = getDiffusionTensor(temperature);

    Mat3x3d Dtt;
    Mat3x3d Dtr;
    Mat3x3d Drr;

    D.getSubMatrix(0, 0, Dtt);
    D.getSubMatrix(3, 0, Dtr);
    D.getSubMatrix(3, 3, Drr);

    // calculate Diffusion Tensor at new location
    Mat3x3d Dpostt;  // translational diffusion tensor at new location
    Mat3x3d Dpostr;  // rotational diffusion tensor at new location
    Mat3x3d Dposrr;  // translation-rotation coupling diffusion tensor
                     // at new location

    Dpostt = Dtt - U * Drr * U + Dtr.transpose() * U - U * Dtr;
    Dposrr = Drr;
    Dpostr = Dtr + Drr * U;

    Mat6x6d Dpos;
    Dpos.setSubMatrix(0, 0, Dpostt);
    Dpos.setSubMatrix(0, 3, Dpostr.transpose());
    Dpos.setSubMatrix(3, 0, Dpostr);
    Dpos.setSubMatrix(3, 3, Dposrr);
    return Dpos;
  }

  Vector3d HydroProp::getCenterOfDiffusion(RealType temperature) {
    // First get the Diffusion tensor at the origin of the coordinate system:

    Vector3d origin(0.0);
    Mat6x6d Do = getDiffusionTensorAtPos(origin, temperature);

    Mat3x3d Dotr;
    Mat3x3d Dorr;

    Do.getSubMatrix(3, 0, Dotr);
    Do.getSubMatrix(3, 3, Dorr);

    Mat3x3d tmp;
    Mat3x3d tmpInv;
    Vector3d tmpVec;

    // Find center of diffusion
    tmp(0, 0) = Dorr(1, 1) + Dorr(2, 2);
    tmp(0, 1) = -Dorr(0, 1);
    tmp(0, 2) = -Dorr(0, 2);
    tmp(1, 0) = -Dorr(0, 1);
    tmp(1, 1) = Dorr(0, 0) + Dorr(2, 2);
    tmp(1, 2) = -Dorr(1, 2);
    tmp(2, 0) = -Dorr(0, 2);
    tmp(2, 1) = -Dorr(1, 2);
    tmp(2, 2) = Dorr(1, 1) + Dorr(0, 0);

    // Vector3d tmpVec;
    tmpVec[0] = Dotr(1, 2) - Dotr(2, 1);
    tmpVec[1] = Dotr(2, 0) - Dotr(0, 2);
    tmpVec[2] = Dotr(0, 1) - Dotr(1, 0);

    // invert tmp Matrix
    invertMatrix(tmp, tmpInv);

    // center of difussion
    Vector3d cod = tmpInv * tmpVec;
    return cod;
  }

  Mat3x3d HydroProp::getPitchMatrix() {
    Mat3x3d P;
    P = - getXitt().inverse() * getXirt();
    return P;
  }

  RealType HydroProp::getScalarPitch() {
    Mat3x3d P = getPitchMatrix();
    Vector3d evals;
    Mat3x3d evects;
    Mat3x3d::diagonalize(P, evals, evects);
    RealType pScalar(0.0);

    for (int i = 0; i < 3; i++) {
      pScalar += pow(evals[i], 2);
    }
    pScalar /= 3.0;

    return sqrt(pScalar);
  }

  void HydroProp::pitchAxes(Mat3x3d& pitchAxes, Vector3d& pitches,
                            RealType& pitchScalar) {
    Mat3x3d P = getPitchMatrix();

    Mat3x3d::diagonalize(P, pitches, pitchAxes);

    pitchScalar = 0.0;
    for (int i = 0; i < 3; i++) {
      pitchScalar += pow(pitches[i], 2);
    }
    pitchScalar /= 3.0;

    pitchScalar = sqrt(pitchScalar);
  }

  Vector3d HydroProp::getCenterOfPitch() {
    Mat3x3d P = getPitchMatrix();
    Vector3d cop;
    cop[0] = 0.5 * (P(1, 2) - P(2, 1));
    cop[1] = 0.5 * (P(2, 0) - P(0, 2));
    cop[2] = 0.5 * (P(0, 1) - P(1, 0));
    // Assume Xi_ was stored at center of resistance, so center of
    // pitch is returned relative to origin if center of resistance is
    // not the origin:
    return (cor_ + cop);
  }

}  // namespace OpenMD
