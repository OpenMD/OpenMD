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

#include "MultiComponentRDF.hpp"

#include <algorithm>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"

namespace OpenMD {

  MultiComponentRDF::MultiComponentRDF(SimInfo* info,
                                       const std::string& filename,
                                       const std::string& sele1,
                                       const std::string& sele2,
                                       unsigned int nbins) :
      StaticAnalyser(info, filename, nbins),
      selectionScript1_(sele1), selectionScript2_(sele2), evaluator1_(info),
      evaluator2_(info), seleMan1_(info), seleMan2_(info),
      sele1_minus_common_(info), sele2_minus_common_(info), common_(info) {
    nPairs_.resize(MaxPairs);

    evaluator1_.loadScriptString(sele1);
    evaluator2_.loadScriptString(sele2);

    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
      validateSelection1(seleMan1_);
    }
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
      validateSelection2(seleMan2_);
    }

    // If all selections are static, we can precompute the number
    // of real pairs.
    if (!evaluator1_.isDynamic() && !evaluator2_.isDynamic()) { calcNPairs(); }
  }

  void MultiComponentRDF::calcNPairs() {
    common_             = seleMan1_ & seleMan2_;
    sele1_minus_common_ = seleMan1_ - common_;
    sele2_minus_common_ = seleMan2_ - common_;
    nSelected1_         = seleMan1_.getSelectionCount();
    nSelected2_         = seleMan2_.getSelectionCount();
    int nIntersect      = common_.getSelectionCount();

    nPairs_[OneOne] = nSelected1_ * (nSelected1_ - 1) / 2;
    nPairs_[OneTwo] =
        nSelected1_ * nSelected2_ - (nIntersect + 1) * nIntersect / 2;
    nPairs_[TwoTwo] = nSelected2_ * (nSelected2_ - 1) / 2;
  }

  void MultiComponentRDF::process() {
    preProcess();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames / step_;

    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      initializeHistograms();

      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
        validateSelection1(seleMan1_);
      }
      if (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
        validateSelection2(seleMan2_);
      }

      // Selections may overlap, and we need a bit of logic to deal
      // with this.
      //
      // |     s1    |
      // | s1 -c | c |
      //         | c | s2 - c |
      //         |    s2      |
      //
      // s1 : Set of StuntDoubles in selection1
      // s2 : Set of StuntDoubles in selection2
      // c  : Intersection of selection1 and selection2
      //
      // When we loop over the pairs, we can divide the looping into 3
      // stages:
      //
      // Stage 1 :     [s1-c]      [s2]
      // Stage 2 :     [c]         [s2 - c]
      // Stage 3 :     [c]         [c]
      // Stages 1 and 2 are completely non-overlapping.
      // Stage 3 is completely overlapping.

      if (evaluator1_.isDynamic() || evaluator2_.isDynamic()) { calcNPairs(); }

      processNonOverlapping(sele1_minus_common_, seleMan2_, OneTwo);
      processNonOverlapping(common_, sele2_minus_common_, OneTwo);

      processOverlapping(common_, OneTwo);
      processOverlapping(seleMan1_, OneOne);
      processOverlapping(seleMan2_, TwoTwo);

      processHistograms();
    }

    postProcess();

    writeRdf();
  }

  void MultiComponentRDF::processNonOverlapping(SelectionManager& sman1,
                                                SelectionManager& sman2,
                                                int pairIndex) {
    StuntDouble* sd1;
    StuntDouble* sd2;
    int i;
    int j;

    // This is the same as a non-overlapping pairwise loop structure:
    // for (int i = 0;  i < ni ; ++i ) {
    //   for (int j = 0; j < nj; ++j) {}
    // }

    for (sd1 = sman1.beginSelected(i); sd1 != NULL;
         sd1 = sman1.nextSelected(i)) {
      for (sd2 = sman2.beginSelected(j); sd2 != NULL;
           sd2 = sman2.nextSelected(j)) {
        collectHistograms(sd1, sd2, pairIndex);
      }
    }
  }

  void MultiComponentRDF::processOverlapping(SelectionManager& sman,
                                             int pairIndex) {
    StuntDouble* sd1;
    StuntDouble* sd2;
    int i;
    int j;

    // This is the same as a pairwise loop structure:
    // for (int i = 0;  i < n-1 ; ++i ) {
    //   for (int j = i + 1; j < n; ++j) {}
    // }

    for (sd1 = sman.beginSelected(i); sd1 != NULL; sd1 = sman.nextSelected(i)) {
      for (j = i, sd2 = sman.nextSelected(j); sd2 != NULL;
           sd2 = sman.nextSelected(j)) {
        collectHistograms(sd1, sd2, pairIndex);
      }
    }
  }
}  // namespace OpenMD
