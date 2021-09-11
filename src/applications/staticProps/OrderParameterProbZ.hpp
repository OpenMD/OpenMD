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

/*
 */

#ifndef APPLICATIONS_STATICPROPS_ORDERPARAMETERPROBZ_HPP
#define APPLICATIONS_STATICPROPS_ORDERPARAMETERPROBZ_HPP

#include <string>
#include <vector>

#include "applications/staticProps/StaticAnalyser.hpp"
#include "brains/SimInfo.hpp"
#include "brains/Thermo.hpp"
#include "math/Vector3.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"

namespace OpenMD {

  class OrderParameterProbZ : public StaticAnalyser {
  public:
    OrderParameterProbZ(SimInfo* info, const std::string& filename,
                        const std::string& sele, const RealType dipoleX,
                        const RealType dipoleY, const RealType dipoleZ,
                        int nbins, int axis = 2);

    virtual void process();

  private:
    virtual void writeOrderCount();

    Snapshot* currentSnapshot_;
    int totalCount_;
    std::string selectionScript_;
    SelectionEvaluator evaluator_;
    SelectionManager seleMan_;
    Thermo thermo_;

    std::vector<RealType> sliceSDCount_;
    std::vector<RealType> flucSliceSDCount_;
    Vector3d refAxis_, dipoleVector_;

    std::vector<RealType> Count_;
    int nbins_, axis_;
    std::string axisLabel_;
  };
}  // namespace OpenMD

#endif
