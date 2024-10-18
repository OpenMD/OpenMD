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

#ifndef APPLICATIONS_STATICPROPS_TETRAHEDRALITYPARAMR_HPP
#define APPLICATIONS_STATICPROPS_TETRAHEDRALITYPARAMR_HPP

#include "applications/staticProps/StaticAnalyser.hpp"
#include "math/Vector3.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"

namespace OpenMD {

  /**
   * @class TetrahedralityParamR
   * @brief Tetrahedrality ParameterR
   *
   * Computes local tetrahedral order parameter Q as introduced in:
   *
   *  "A new order parameter for tetrahedral configurations," by P.-L. Chau and
   *    A.J. Hardwick, Mol. Phys. 93, pp. 511-518 (1998).
   *
   *
   * Note that we use a rescaled version of the tetrahedral order
   * parameter 'Q' such that a perfectly tetrahedral configuration has a Q value
   * of 1 and an ideal gas configuration has a Q value of 0. This rescaled
   * version of the tetrahedrality parameter was first introduced in:
   *
   *  "Relationship between structural order and the anomalies of liquid water,"
   *    by J.R. Errington and P.G. Debenedetti, Nature 409, pp. 318-321 (2001).
   *
   *
   * Characterization of the spatial correlations of the the local
   * order parameter Q are done according to the procedure outlined
   * in:
   *
   *   "Space-time correlations in the orientational order parameter and the
   *    orientational entropy of water," by P. Kumar, S.V. Buldyrev, and
   *    H.E. Stanley, arXiv:0807.4699v1 [cond-mat.soft] 29 Jul 2008.
   *
   */
  class TetrahedralityParamR : public StaticAnalyser {
  public:
    TetrahedralityParamR(SimInfo* info, const std::string& filename,
                         const std::string& sele1, const std::string& sele2,
                         const std::string& sele3, RealType rCut, RealType len,
                         int nrbins);

    virtual void process();

  private:
    void writeQr();

    Snapshot* currentSnapshot_;
    std::string selectionScript1_;
    std::string selectionScript2_;
    std::string selectionScript3_;
    SelectionManager seleMan1_;
    SelectionManager seleMan2_;
    SelectionManager seleMan3_;
    SelectionEvaluator evaluator1_;
    SelectionEvaluator evaluator2_;
    SelectionEvaluator evaluator3_;
    RealType len_;
    RealType rCut_;
    RealType deltaR_;
    int nBins_;
    std::vector<RealType> sliceQ_;
    std::vector<int> sliceCount_;
  };
}  // namespace OpenMD

#endif
