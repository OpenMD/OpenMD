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

#ifndef APPLICATIONS_STATICPROPS_TETRAHEDRALITYPARAMXYZ_HPP
#define APPLICATIONS_STATICPROPS_TETRAHEDRALITYPARAMXYZ_HPP

#include "applications/staticProps/StaticAnalyser.hpp"
#include "math/Vector3.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"

namespace OpenMD {

  /**
   * @class TetrahedralityParamXYZ
   * @brief Tetrahedrality Parameter XYZ
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
  class TetrahedralityParamXYZ : public StaticAnalyser {
  public:
    TetrahedralityParamXYZ(SimInfo* info, const std::string& filename,
                           const std::string& sele1, const std::string& sele2,
                           RealType rCut, RealType voxelSize,
                           RealType gaussWidth);

    virtual void process();

  private:
    void writeQxyz();
    void writeVTKGrid();

    Snapshot* currentSnapshot_;
    std::string selectionScript1_;
    std::string selectionScript2_;
    SelectionManager seleMan1_;
    SelectionManager seleMan2_;
    SelectionEvaluator evaluator1_;
    SelectionEvaluator evaluator2_;

    RealType rCut_;
    RealType voxelSize_;
    RealType gaussWidth_;

    Vector3i nBins_;
    std::vector<std::vector<std::vector<RealType>>> count_;
    std::vector<std::vector<std::vector<RealType>>> hist_;
  };
}  // namespace OpenMD

#endif
