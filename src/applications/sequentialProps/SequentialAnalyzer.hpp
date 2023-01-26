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

#ifndef APPLICATIONS_SEQUENTIALPROPS_SEQUENTIALANALYZER_HPP
#define APPLICATIONS_SEQUENTIALPROPS_SEQUENTIALANALYZER_HPP

#include <string>
#include <vector>

#include "brains/BlockSnapshotManager.hpp"
#include "brains/SimInfo.hpp"
#include "primitives/StuntDouble.hpp"
#include "selection/SelectionEvaluator.hpp"
#include "selection/SelectionManager.hpp"

namespace OpenMD {

  /**
   * @class SequentialAnalyzer SequentialAnalyzer.hpp
   * "applications/sequentialProps/SequentialAnalyzer"
   * @brief Base class for Sequence Analyzer
   */

  class SequentialAnalyzer {
  public:
    SequentialAnalyzer(SimInfo* info, const std::string& filename,
                       const std::string& sele1, const std::string& sele2);

    virtual ~SequentialAnalyzer() = default;
    virtual void doSequence();

    void setOutputName(const std::string& filename) {
      outputFilename_ = filename;
    }

    const std::string& getOutputFileName() const { return outputFilename_; }

    void setStep(int step) {
      assert(step > 0);
      step_ = step;
    }

    int getStep() { return step_; }

    const std::string& getSequenceType() const { return sequenceType_; }

    void setSequenceType(const std::string& type) { sequenceType_ = type; }

    void setParameterString(const std::string& params) {
      paramString_ = params;
    }

  protected:
    virtual void preSequence() {}
    virtual void postSequence() {}
    virtual void writeSequence();
    virtual void doFrame(int frame) = 0;

    SimInfo* info_ {nullptr};
    Snapshot* currentSnapshot_;
    std::string dumpFilename_;

    SelectionManager seleMan1_;
    std::string selectionScript1_;
    SelectionEvaluator evaluator1_;

    SelectionManager seleMan2_;
    std::string selectionScript2_;
    SelectionEvaluator evaluator2_;

    int step_;

    std::string outputFilename_;
    int frame_;
    // int storageLayout_;
    std::vector<RealType> times_;
    std::vector<RealType> values_;
    std::string sequenceType_;
    std::string paramString_;
  };
}  // namespace OpenMD

#endif
