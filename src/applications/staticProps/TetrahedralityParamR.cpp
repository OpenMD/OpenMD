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

#include "applications/staticProps/TetrahedralityParamR.hpp"

#include <algorithm>
#include <fstream>
#include <vector>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Revision.hpp"
#include "utils/simError.h"

#define HONKING_LARGE_VALUE 1.0e10

using namespace std;
namespace OpenMD {
  TetrahedralityParamR::TetrahedralityParamR(
      SimInfo* info, const std::string& filename, const std::string& sele1,
      const std::string& sele2, const std::string& sele3, RealType rCut,
      RealType len, int nrbins) :
      StaticAnalyser(info, filename, nrbins),
      selectionScript1_(sele1), selectionScript2_(sele2),
      selectionScript3_(sele3), seleMan1_(info), seleMan2_(info),
      seleMan3_(info), evaluator1_(info), evaluator2_(info), evaluator3_(info),
      len_(len), nBins_(nrbins) {
    setAnalysisType("Tetrahedrality Parameter(r)");

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    evaluator2_.loadScriptString(sele2);
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }

    evaluator3_.loadScriptString(sele3);
    if (!evaluator3_.isDynamic()) {
      seleMan3_.setSelectionSet(evaluator3_.evaluate());
    }

    // Set up cutoff radius:
    rCut_ = rCut;

    deltaR_ = len_ / nBins_;

    // fixed number of bins
    sliceQ_.resize(nBins_);
    sliceCount_.resize(nBins_);
    std::fill(sliceQ_.begin(), sliceQ_.end(), 0.0);
    std::fill(sliceCount_.begin(), sliceCount_.end(), 0);

    setOutputName(getPrefix(filename) + ".Qr");
  }

  void TetrahedralityParamR::process() {
    StuntDouble* sd;
    StuntDouble* sd2;
    StuntDouble* sd3;
    StuntDouble* sdi;
    StuntDouble* sdj;
    int myIndex;
    Vector3d vec;
    Vector3d ri, rj, rk, rik, rkj;
    RealType r;
    RealType cospsi;
    RealType Qk;
    std::vector<std::pair<RealType, StuntDouble*>> myNeighbors;
    int isd1;
    int isd2;
    int isd3;
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }

      if (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
      }

      if (evaluator3_.isDynamic()) {
        seleMan3_.setSelectionSet(evaluator3_.evaluate());
      }

      // outer loop is over the selected StuntDoubles:
      for (sd = seleMan1_.beginSelected(isd1); sd != NULL;
           sd = seleMan1_.nextSelected(isd1)) {
        myIndex = sd->getGlobalIndex();

        Qk = 1.0;
        myNeighbors.clear();

        for (sd2 = seleMan2_.beginSelected(isd2); sd2 != NULL;
             sd2 = seleMan2_.nextSelected(isd2)) {
          if (sd2->getGlobalIndex() != myIndex) {
            vec = sd->getPos() - sd2->getPos();

            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);

            r = vec.length();

            // Check to see if neighbor is in bond cutoff

            if (r < rCut_) { myNeighbors.push_back(std::make_pair(r, sd2)); }
          }
        }

        // Sort the vector using predicate and std::sort
        std::sort(myNeighbors.begin(), myNeighbors.end());

        // Use only the 4 closest neighbors to do the rest of the work:

        int nbors = myNeighbors.size() > 4 ? 4 : myNeighbors.size();
        int nang  = int(0.5 * (nbors * (nbors - 1)));

        rk = sd->getPos();

        for (int i = 0; i < nbors - 1; i++) {
          sdi = myNeighbors[i].second;
          ri  = sdi->getPos();
          rik = rk - ri;
          if (usePeriodicBoundaryConditions_) currentSnapshot_->wrapVector(rik);

          rik.normalize();

          for (int j = i + 1; j < nbors; j++) {
            sdj = myNeighbors[j].second;
            rj  = sdj->getPos();
            rkj = rk - rj;
            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(rkj);
            rkj.normalize();

            cospsi = dot(rik, rkj);

            // Calculates scaled Qk for each molecule using calculated
            // angles from 4 or fewer nearest neighbors.
            Qk -= (pow(cospsi + 1.0 / 3.0, 2) * 2.25 / nang);
          }
        }

        if (nang > 0) {
          RealType shortest = HONKING_LARGE_VALUE;

          // loop over selection 3 to find closest atom in selection 3:
          for (sd3 = seleMan3_.beginSelected(isd3); sd3 != NULL;
               sd3 = seleMan3_.nextSelected(isd3)) {
            vec = sd->getPos() - sd3->getPos();

            if (usePeriodicBoundaryConditions_)
              currentSnapshot_->wrapVector(vec);

            r = vec.length();

            if (r < shortest) shortest = r;
          }

          int whichBin = int(shortest / deltaR_);
          if (whichBin < int(nBins_)) {
            sliceQ_[whichBin] += Qk;
            sliceCount_[whichBin] += 1;
          }
        }
      }
    }
    writeQr();
  }

  void TetrahedralityParamR::writeQr() {
    Revision rev;
    std::ofstream qRstream(outputFilename_.c_str());
    if (qRstream.is_open()) {
      qRstream << "# " << getAnalysisType() << "\n";
      qRstream << "# OpenMD " << rev.getFullRevision() << "\n";
      qRstream << "# " << rev.getBuildDate() << "\n";
      qRstream << "#selection 1: (" << selectionScript1_ << ")\n";
      qRstream << "#selection 2: (" << selectionScript2_ << ")\n";
      qRstream << "#selection 3: (" << selectionScript3_ << ")\n";
      if (!paramString_.empty())
        qRstream << "# parameters: " << paramString_ << "\n";

      qRstream << "#distance"
               << "\tQk\n";
      for (unsigned int i = 0; i < sliceQ_.size(); ++i) {
        RealType Rval = (i + 0.5) * deltaR_;
        if (sliceCount_[i] != 0) {
          qRstream << Rval << "\t" << sliceQ_[i] / (RealType)sliceCount_[i]
                   << "\n";
        }
      }

    } else {
      sprintf(painCave.errMsg, "TetrahedralityParamR: unable to open %s\n",
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    qRstream.close();
  }
}  // namespace OpenMD
