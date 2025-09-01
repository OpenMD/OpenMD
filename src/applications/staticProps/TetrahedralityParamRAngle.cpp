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

#include "applications/staticProps/TetrahedralityParamRAngle.hpp"

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
  TetrahedralityParamRAngle::TetrahedralityParamRAngle(
      SimInfo* info, const std::string& filename, const std::string& sele1,
      const std::string& sele2, const std::string& sele3, const int seleOffset3,
      RealType rCut, RealType len, int nrbins, int nangleBins) :
      StaticAnalyser(info, filename, nrbins),
      selectionScript1_(sele1), selectionScript2_(sele2),
      selectionScript3_(sele3), seleMan1_(info), seleMan2_(info),
      seleMan3_(info), evaluator1_(info), evaluator2_(info), evaluator3_(info),
      seleOffset3_(seleOffset3), len_(len), nBins_(nrbins), nAngleBins_(nangleBins) {
    
    setAnalysisType("Tetrahedrality Parameter(r, cos(theta))");

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
    deltaCosAngle_ = 2.0 / (double)nAngleBins_;

    // fixed number of bins
    sliceQ_.resize(nBins_);
    sliceCount_.resize(nBins_);
    for (unsigned int i = 0; i < nBins_; ++i) {
      sliceQ_[i].resize(nAngleBins_);
      sliceCount_[i].resize(nAngleBins_);
      std::fill(sliceQ_[i].begin(), sliceQ_[i].end(), 0.0);
      std::fill(sliceCount_[i].begin(), sliceCount_[i].end(), 0);
    }

    setOutputName(getPrefix(filename) + ".Qrangle");
  }

  void TetrahedralityParamRAngle::process() {
    StuntDouble* sd;
    StuntDouble* sd2;
    StuntDouble* sd3;
    StuntDouble* sdi;
    StuntDouble* sdj;
    int myIndex;
    Vector3d vec;
    Vector3d ri, rj, rk, rik, rkj;
    Vector3d r34;
    RealType r;
    RealType cospsi;
    RealType Qk;
    std::vector<std::pair<RealType, StuntDouble*>> myNeighbors;
    int isd1;
    int isd2;
    int isd3;
    int whichBin;
    int whichThetaBin;
    RealType halfBin;
    RealType cosAngle;
    
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
	    
            if (usePeriodicBoundaryConditions_) {
              currentSnapshot_->wrapVector(vec);
	    }

            r = vec.length();
	    
            if (r < shortest) {
	      shortest = r;
	      vec.normalize();
	      whichBin = int(r / deltaR_);
	      
	      // Find the angle that vec makes with the 3->4 vector and figure out the bin:
	      int sd4ind       = sd3->getGlobalIndex() + seleOffset3_;
	      StuntDouble* sd4 = info_->getIOIndexToIntegrableObject(sd4ind);	     
	      r34 = sd4->getPos() - sd3->getPos();
	      if (usePeriodicBoundaryConditions_) {
		currentSnapshot_->wrapVector(r34);
	      }
	      r34.normalize();	      
	      cosAngle = dot(r34, vec);
	      
	      halfBin  = (nAngleBins_ - 1) * 0.5;
	      whichThetaBin = int(halfBin * (cosAngle + 1.0));
	    }	      
          }

	  if (whichBin < int(nBins_)) {
	    sliceQ_[whichBin][whichThetaBin] += Qk;
	    sliceCount_[whichBin][whichThetaBin] += 1;
          }
        }
      }
    }
    writeQrAngle();
  }

  void TetrahedralityParamRAngle::writeQrAngle() {
    Revision rev;
    std::ofstream qRstream(outputFilename_.c_str());
    if (qRstream.is_open()) {
      qRstream << "# " << getAnalysisType() << "\n";
      qRstream << "# OpenMD " << rev.getFullRevision() << "\n";
      qRstream << "# " << rev.getBuildDate() << "\n";
      qRstream << "#selection 1: (" << selectionScript1_ << ")\n";
      qRstream << "#selection 2: (" << selectionScript2_ << ")\n";
      qRstream << "#selection 3: (" << selectionScript3_ << ")\n";
      qRstream << "#selection 3 offset: " << seleOffset3_ << "\n";
      if (!paramString_.empty())
        qRstream << "# parameters: " << paramString_ << "\n";

      for (unsigned int i = 0; i < sliceQ_.size(); ++i) {
        for (unsigned int j = 0; j < sliceQ_[i].size(); ++j) {	
	  if (sliceCount_[i][j] != 0) {
	    qRstream << sliceQ_[i][j] / (RealType)sliceCount_[i][j] << "\t";
	  } else {
	    qRstream << 0.0 << "\t";
	  }
	}
	qRstream << "\n";
      }
      
    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "TetrahedralityParamRAngle: unable to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }
    qRstream.close();
  }
}  // namespace OpenMD
