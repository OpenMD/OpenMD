/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
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
