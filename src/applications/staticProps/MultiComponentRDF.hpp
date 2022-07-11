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

#ifndef APPLICATIONS_STATICPROPS_RADIALDISTRFUNC_HPP
#define APPLICATIONS_STATICPROPS_RADIALDISTRFUNC_HPP

#include <string>
#include <unordered_map>

#include "applications/staticProps/StaticAnalyser.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"
#include "utils/Constants.hpp"

namespace OpenMD {

  /**
   * @class MultiComponentRDF
   * @brief Multi-Component Radial Distribution Function
   */
  class MultiComponentRDF : public StaticAnalyser {
  public:
    MultiComponentRDF(SimInfo* info, const std::string& filename,
                      const std::string& sele1, const std::string& sele2,
                      unsigned int nbins);

    void process();

  protected:
    virtual void preProcess() {}
    virtual void postProcess() {}

    void processSelections(SelectionManager& sman1, SelectionManager& sman2);

    virtual void processNonOverlapping(SelectionManager& sman1,
                                       SelectionManager& sman2);
    virtual void processOverlapping(SelectionManager& sman);

    std::unordered_map<int, int> getNPairs() const { return nPairs_; }
    int getNSelected1() const { return nSelected1_; }
    int getNSelected2() const { return nSelected2_; }

    Snapshot* currentSnapshot_;

    std::string selectionScript1_;
    std::string selectionScript2_;
    int nProcessed_;
    SelectionEvaluator evaluator1_;
    SelectionEvaluator evaluator2_;

    SelectionManager seleMan1_;
    SelectionManager seleMan2_;
    SelectionManager sele1_minus_common_;
    SelectionManager sele2_minus_common_;
    SelectionManager common_;

  private:
    virtual void initializeHistogram() {}
    virtual void collectHistogram(StuntDouble* sd1, StuntDouble* sd2) = 0;
    virtual void processHistogram() {}

    virtual void validateSelection1(SelectionManager&) {}
    virtual void validateSelection2(SelectionManager&) {}
    virtual void writeRdf() = 0;

    std::unordered_map<int, int> nPairs_;
    int nSelected1_;
    int nSelected2_;
  };
}  // namespace OpenMD

#endif
