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

#include "hydrodynamics/HydroProp.hpp"

#include "math/CholeskyDecomposition.hpp"
#include "utils/StringTokenizer.hpp"

namespace OpenMD {

  HydroProp::HydroProp() : hasCOR(false), hasXi(false) {}

  HydroProp::HydroProp(Vector3d cor, Mat6x6d Xi, Mat6x6d D) :
      cor_(cor), Xi_(Xi), D_(D), hasCOR(true), hasXi(true) {}

  HydroProp::HydroProp(const std::string& frictionLine) :
      hasCOR(false), hasXi(false) {
    StringTokenizer tokenizer(frictionLine);
    if (tokenizer.countTokens() >= 40) {
      name_   = tokenizer.nextToken();
      cor_[0] = tokenizer.nextTokenAsDouble();
      cor_[1] = tokenizer.nextTokenAsDouble();
      cor_[2] = tokenizer.nextTokenAsDouble();

      hasCOR = true;
      Mat3x3d Xitt(0.0);
      Mat3x3d Xirt(0.0);
      Mat3x3d Xitr(0.0);
      Mat3x3d Xirr(0.0);

      Xitt(0, 0) = tokenizer.nextTokenAsDouble();
      Xitt(0, 1) = tokenizer.nextTokenAsDouble();
      Xitt(0, 2) = tokenizer.nextTokenAsDouble();
      Xitt(1, 0) = tokenizer.nextTokenAsDouble();
      Xitt(1, 1) = tokenizer.nextTokenAsDouble();
      Xitt(1, 2) = tokenizer.nextTokenAsDouble();
      Xitt(2, 0) = tokenizer.nextTokenAsDouble();
      Xitt(2, 1) = tokenizer.nextTokenAsDouble();
      Xitt(2, 2) = tokenizer.nextTokenAsDouble();

      Xirt(0, 0) = tokenizer.nextTokenAsDouble();
      Xirt(0, 1) = tokenizer.nextTokenAsDouble();
      Xirt(0, 2) = tokenizer.nextTokenAsDouble();
      Xirt(1, 0) = tokenizer.nextTokenAsDouble();
      Xirt(1, 1) = tokenizer.nextTokenAsDouble();
      Xirt(1, 2) = tokenizer.nextTokenAsDouble();
      Xirt(2, 0) = tokenizer.nextTokenAsDouble();
      Xirt(2, 1) = tokenizer.nextTokenAsDouble();
      Xirt(2, 2) = tokenizer.nextTokenAsDouble();

      Xitr(0, 0) = tokenizer.nextTokenAsDouble();
      Xitr(0, 1) = tokenizer.nextTokenAsDouble();
      Xitr(0, 2) = tokenizer.nextTokenAsDouble();
      Xitr(1, 0) = tokenizer.nextTokenAsDouble();
      Xitr(1, 1) = tokenizer.nextTokenAsDouble();
      Xitr(1, 2) = tokenizer.nextTokenAsDouble();
      Xitr(2, 0) = tokenizer.nextTokenAsDouble();
      Xitr(2, 1) = tokenizer.nextTokenAsDouble();
      Xitr(2, 2) = tokenizer.nextTokenAsDouble();

      Xirr(0, 0) = tokenizer.nextTokenAsDouble();
      Xirr(0, 1) = tokenizer.nextTokenAsDouble();
      Xirr(0, 2) = tokenizer.nextTokenAsDouble();
      Xirr(1, 0) = tokenizer.nextTokenAsDouble();
      Xirr(1, 1) = tokenizer.nextTokenAsDouble();
      Xirr(1, 2) = tokenizer.nextTokenAsDouble();
      Xirr(2, 0) = tokenizer.nextTokenAsDouble();
      Xirr(2, 1) = tokenizer.nextTokenAsDouble();
      Xirr(2, 2) = tokenizer.nextTokenAsDouble();

      Xi_.setSubMatrix(0, 0, Xitt);
      Xi_.setSubMatrix(0, 3, Xirt);
      Xi_.setSubMatrix(3, 0, Xitr);
      Xi_.setSubMatrix(3, 3, Xirr);

      hasXi = true;
      complete();
    }
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

  Mat3x3d HydroProp::getPitchMatrix() {
    Mat3x3d P;
    P = Constants::TWO_PI * getXitt().inverse() * getXirt();
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

  void HydroProp::complete() {
    if (hasXi) { CholeskyDecomposition(Xi_, S_); }
  }

}  // namespace OpenMD
