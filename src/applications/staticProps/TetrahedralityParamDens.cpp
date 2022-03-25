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

#include "applications/staticProps/TetrahedralityParamDens.hpp"

#include <algorithm>
#include <fstream>
#include <vector>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {
  TetrahedralityParamDens::TetrahedralityParamDens(SimInfo* info,
                                                   const std::string& filename,
                                                   const std::string& sele1,
                                                   const std::string& sele2,
                                                   double rCut, int ndensbins) :
      StaticAnalyser(info, filename, ndensbins),
      selectionScript1_(sele1), selectionScript2_(sele2), seleMan1_(info),
      seleMan2_(info), evaluator1_(info), evaluator2_(info), rCut_(rCut) {
    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    evaluator2_.loadScriptString(sele2);
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }

    // Q can take values of 0 to 1.
    MinQ_   = -3.0;
    MaxQ_   = 1.1;
    deltaQ_ = (MaxQ_ - MinQ_) / nBins_;

    // zeroing the number of molecules investigated counter
    count_ = 0;

    // fixed number of bins
    sliceCount_.resize(nBins_);
    std::fill(sliceCount_.begin(), sliceCount_.end(), 0);

    setOutputName(getPrefix(filename) + ".Qdens");
  }

  void TetrahedralityParamDens::process() {
    StuntDouble* sd;
    StuntDouble* sd2;
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
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      // Mat3x3d hmat = currentSnapshot_->getHmat();
      // zBox_.push_back(hmat(2,2));

      // RealType halfBoxZ_ = hmat(2,2) / 2.0;

      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }

      if (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
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

            if (r < rCut_) {
              myNeighbors.push_back(std::make_pair(r, sd2));
              if (myIndex == 513) {
                // std::cerr<< "sd2 index = " << sd2->getGlobalIndex() << "\n";
              }
            }
          }
        }

        // Sort the vector using predicate and std::sort
        std::sort(myNeighbors.begin(), myNeighbors.end());

        // Use only the 4 closest neighbors to do the rest of the work:

        int nbors = myNeighbors.size() > 4 ? 4 : myNeighbors.size();
        int nang  = int(0.5 * (nbors * (nbors - 1)));
        if (nang < 4) {
          // std::cerr << "nbors = " << nbors << " nang = " << nang << "\n";
        }
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
            /*if (Qk < 0 ) {
            std::cerr << "(pow(cospsi + 1.0 / 3.0, 2) * 2.25 / nang)" <<
            (pow(cospsi + 1.0 / 3.0, 2) * 2.25 / nang) << "\n"; std::cerr << "Qk
            = " << Qk << " nang = " << nang <<  " nbors = " << nbors << "\n";
            std::cerr << "myIndex = " << myIndex <<
            "\n";
            }*/
          }
        }

        // Based on the molecule's q value, populate the histogram.
        if (nang > 0) {
          int binNo = int((Qk - MinQ_) / deltaQ_);

          if (binNo < int(sliceCount_.size()) && binNo >= 0) {
            sliceCount_[binNo] += 1;
          } else {
            std::cerr << "binNo = " << binNo << " Qk = " << Qk << "\n";
            // std::cerr << "nbors = " << nbors << "\n";
          }

          // std::cout << "count = " << count_ << "\n";
          count_++;
        }
      }
    }
    writeQdens();
  }

  void TetrahedralityParamDens::writeQdens() {
    std::ofstream qdensstream(outputFilename_.c_str());
    if (qdensstream.is_open()) {
      qdensstream << "#Tetrahedrality Parameters \n";
      qdensstream << "#nMolecules:\t" << count_ << " \n";
      qdensstream << "#selection 1: (" << selectionScript1_ << ")\n";
      qdensstream << "#selection 2: (" << selectionScript2_ << ")\n";
      qdensstream << "#Qk\tfractional probability \n";

      for (unsigned int i = 0; i < sliceCount_.size(); ++i) {
        RealType q = MinQ_ + (i + 0.5) * deltaQ_;
        if (count_ != 0) {
          qdensstream << q << "\t" << sliceCount_[i] / RealType(count_) << "\n";
        }
      }

    } else {
      sprintf(painCave.errMsg, "TetrahedralityParamDens: unable to open %s\n",
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    qdensstream.close();
  }
}  // namespace OpenMD
