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

#ifndef APPLICATIONS_STATICPROPS_GCNSEQ_HPP
#define APPLICATIONS_STATICPROPS_GCNSEQ_HPP

#include <string>
#include <vector>

#include "applications/sequentialProps/SequentialAnalyzer.hpp"

namespace OpenMD {

  /**
   * @class GCNSeq
   * @brief Generalized Coordinate Number Sequence
   *
   * Computes a frame-by-frame sequence of distributions of
   * generalized coordinate numbers as described in:
   *
   *   "Finding optimal surface sites on heterogeneous catalysts by
   *    counting nearest neighbors," by F. Calle-Vallejo et al.,
   *    Science 350(6257) pp. 185-189 (2015).
   *    http://dx.doi.org/10.1126/science.aab3501
   *
   * Note that extra parameters mussed be declared:
   *
   *   rCut = cutoff radius for finding lists of nearest neighbors
   *   sele1 = selection of StuntDoubles used for the GCN distribution
   *   sele2 = selection of StuntDoubles used for nearest neighbor computation
   */
  class GCNSeq : public SequentialAnalyzer {
  public:
    GCNSeq(SimInfo* info, const std::string& filename, const std::string& sele1,
           const std::string& sele2, RealType rCut, int bins);

    virtual void doFrame(int istep);
    virtual void writeSequence();

  private:
    RealType rCut_;
    int bins_;

    bool usePBC_;
    int nnMax_;
    RealType delta_;
    int selectionCount1_;
    int selectionCount2_;

    std::vector<int> count_;
    std::vector<std::vector<RealType>> histogram_;
  };
}  // namespace OpenMD

#endif
