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

      Xitt_(0, 0) = tokenizer.nextTokenAsDouble();
      Xitt_(0, 1) = tokenizer.nextTokenAsDouble();
      Xitt_(0, 2) = tokenizer.nextTokenAsDouble();
      Xitt_(1, 0) = tokenizer.nextTokenAsDouble();
      Xitt_(1, 1) = tokenizer.nextTokenAsDouble();
      Xitt_(1, 2) = tokenizer.nextTokenAsDouble();
      Xitt_(2, 0) = tokenizer.nextTokenAsDouble();
      Xitt_(2, 1) = tokenizer.nextTokenAsDouble();
      Xitt_(2, 2) = tokenizer.nextTokenAsDouble();

      Xirt_(0, 0) = tokenizer.nextTokenAsDouble();
      Xirt_(0, 1) = tokenizer.nextTokenAsDouble();
      Xirt_(0, 2) = tokenizer.nextTokenAsDouble();
      Xirt_(1, 0) = tokenizer.nextTokenAsDouble();
      Xirt_(1, 1) = tokenizer.nextTokenAsDouble();
      Xirt_(1, 2) = tokenizer.nextTokenAsDouble();
      Xirt_(2, 0) = tokenizer.nextTokenAsDouble();
      Xirt_(2, 1) = tokenizer.nextTokenAsDouble();
      Xirt_(2, 2) = tokenizer.nextTokenAsDouble();

      Xitr_(0, 0) = tokenizer.nextTokenAsDouble();
      Xitr_(0, 1) = tokenizer.nextTokenAsDouble();
      Xitr_(0, 2) = tokenizer.nextTokenAsDouble();
      Xitr_(1, 0) = tokenizer.nextTokenAsDouble();
      Xitr_(1, 1) = tokenizer.nextTokenAsDouble();
      Xitr_(1, 2) = tokenizer.nextTokenAsDouble();
      Xitr_(2, 0) = tokenizer.nextTokenAsDouble();
      Xitr_(2, 1) = tokenizer.nextTokenAsDouble();
      Xitr_(2, 2) = tokenizer.nextTokenAsDouble();

      Xirr_(0, 0) = tokenizer.nextTokenAsDouble();
      Xirr_(0, 1) = tokenizer.nextTokenAsDouble();
      Xirr_(0, 2) = tokenizer.nextTokenAsDouble();
      Xirr_(1, 0) = tokenizer.nextTokenAsDouble();
      Xirr_(1, 1) = tokenizer.nextTokenAsDouble();
      Xirr_(1, 2) = tokenizer.nextTokenAsDouble();
      Xirr_(2, 0) = tokenizer.nextTokenAsDouble();
      Xirr_(2, 1) = tokenizer.nextTokenAsDouble();
      Xirr_(2, 2) = tokenizer.nextTokenAsDouble();

      Xi_.setSubMatrix(0, 0, Xitt_);
      Xi_.setSubMatrix(0, 3, Xirt_);
      Xi_.setSubMatrix(3, 0, Xitr_);
      Xi_.setSubMatrix(3, 3, Xirr_);

      hasXi = true;

      CholeskyDecomposition(Xi_, S_);
    }
  }

  void HydroProp::complete() {
    if (hasXi) {
      for (int i1 = 0; i1 < 3; i1++) {
        for (int j1 = 0; j1 < 3; j1++) {
          Xitt_(i1, j1) = Xi_(i1, j1);
          Xirt_(i1, j1) = Xi_(i1, j1 + 3);
          Xitr_(i1, j1) = Xi_(i1 + 3, j1);
          Xirr_(i1, j1) = Xi_(i1 + 3, j1 + 3);
        }
      }
      CholeskyDecomposition(Xi_, S_);
    }
  }

}  // namespace OpenMD
