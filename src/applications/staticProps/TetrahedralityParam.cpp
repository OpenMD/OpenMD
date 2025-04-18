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

#include "applications/staticProps/TetrahedralityParam.hpp"

#include <vector>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

namespace OpenMD {

  TetrahedralityParam::TetrahedralityParam(SimInfo* info,
                                           const std::string& filename,
                                           const std::string& sele, double rCut,
                                           int nbins) :
      StaticAnalyser(info, filename, nbins),
      selectionScript_(sele), seleMan_(info), evaluator_(info) {
    setOutputName(getPrefix(filename) + ".q");

    evaluator_.loadScriptString(sele);
    if (!evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    // Set up cutoff radius:

    rCut_ = rCut;

    Q_histogram_.resize(nBins_);

    // Q can take values from 0 to 1

    MinQ_   = 0.0;
    MaxQ_   = 1.1;
    deltaQ_ = (MaxQ_ - MinQ_) / nBins_;
  }

  void TetrahedralityParam::initializeHistogram() {
    std::fill(Q_histogram_.begin(), Q_histogram_.end(), 0);
  }

  void TetrahedralityParam::process() {
    Molecule* mol;
    StuntDouble* sd;
    StuntDouble* sd2;
    StuntDouble* sdi;
    StuntDouble* sdj;
    int myIndex;
    SimInfo::MoleculeIterator mi;
    Molecule::IntegrableObjectIterator ioi;
    Vector3d vec;
    Vector3d ri, rj, rk, rik, rkj, dposition, tposition;
    RealType r;
    RealType cospsi;
    RealType Qk;
    std::vector<std::pair<RealType, StuntDouble*>> myNeighbors;
    int isd;
    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    DumpReader reader(info_, dumpFilename_);
    int nFrames   = reader.getNFrames();
    frameCounter_ = 0;

    Distorted_.clear();
    Tetrahedral_.clear();

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      frameCounter_++;
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      // outer loop is over the selected StuntDoubles:

      for (sd = seleMan_.beginSelected(isd); sd != NULL;
           sd = seleMan_.nextSelected(isd)) {
        myIndex = sd->getGlobalIndex();
        Qk      = 1.0;

        myNeighbors.clear();

        // inner loop is over all StuntDoubles in the system:

        for (mol = info_->beginMolecule(mi); mol != NULL;
             mol = info_->nextMolecule(mi)) {
          for (sd2 = mol->beginIntegrableObject(ioi); sd2 != NULL;
               sd2 = mol->nextIntegrableObject(ioi)) {
            if (sd2->getGlobalIndex() != myIndex) {
              vec = sd->getPos() - sd2->getPos();

              if (usePeriodicBoundaryConditions_)
                currentSnapshot_->wrapVector(vec);

              r = vec.length();

              // Check to see if neighbor is in bond cutoff

              if (r < rCut_) { myNeighbors.push_back(std::make_pair(r, sd2)); }
            }
          }
        }

        // Sort the vector using predicate and std::sort
        std::sort(myNeighbors.begin(), myNeighbors.end());

        // std::cerr << myNeighbors.size() <<  " neighbors within "
        //          << rCut_  << " A" << " \n";

        // Use only the 4 closest neighbors to do the rest of the work:

        int nbors = myNeighbors.size() > 4 ? 4 : myNeighbors.size();
        int nang  = int(0.5 * (nbors * (nbors - 1)));

        rk = sd->getPos();
        // std::cerr<<nbors<<endl;
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

            // std::cerr << "cos(psi) = " << cospsi << " \n";

            // Calculates scaled Qk for each molecule using calculated
            // angles from 4 or fewer nearest neighbors.
            Qk = Qk - (pow(cospsi + 1.0 / 3.0, 2) * 2.25 / nang);
            // std::cerr<<Qk<<"\t"<<nang<<endl;
          }
        }
        // std::cerr<<nang<<endl;
        if (nang > 0) {
          collectHistogram(Qk);

          // Saves positions of StuntDoubles & neighbors with distorted
          // coordination (low Qk value)
          if ((Qk < 0.55) && (Qk > 0.45)) {
            // std::cerr<<Distorted_.size()<<endl;
            Distorted_.push_back(sd);
            // std::cerr<<Distorted_.size()<<endl;
            dposition = sd->getPos();
            // std::cerr << "distorted position \t" << dposition << "\n";
          }

          // Saves positions of StuntDoubles & neighbors with
          // tetrahedral coordination (high Qk value)
          if (Qk > 0.05) {
            Tetrahedral_.push_back(sd);

            tposition = sd->getPos();
            // std::cerr << "tetrahedral position \t" << tposition << "\n";
          }

          // std::cerr<<Tetrahedral_.size()<<endl;
        }
      }
    }

    writeOrderParameter();
    std::cerr << "number of distorted StuntDoubles = " << Distorted_.size()
              << "\n";
    std::cerr << "number of tetrahedral StuntDoubles = " << Tetrahedral_.size()
              << "\n";
  }

  void TetrahedralityParam::collectHistogram(RealType Qk) {
    if (Qk > MinQ_ && Qk < MaxQ_) {
      int whichBin = int((Qk - MinQ_) / deltaQ_);
      Q_histogram_[whichBin] += 1;
    }
  }

  void TetrahedralityParam::writeOrderParameter() {
    int nSelected = 0;

    for (int i = 0; i < nBins_; ++i) {
      nSelected = nSelected + int(Q_histogram_[i] * deltaQ_);
    }

    std::ofstream osq((getOutputFileName() + "Q").c_str());

    if (osq.is_open()) {
      osq << "# Tetrahedrality Parameters\n";
      osq << "# selection: (" << selectionScript_ << ")\n";
      osq << "# \n";
      // Normalize by number of frames and write it out:
      for (int i = 0; i < nBins_; ++i) {
        RealType Qval = MinQ_ + (i + 0.5) * deltaQ_;
        osq << Qval;
        osq << "\t" << (RealType)(Q_histogram_[i] / deltaQ_) / nSelected;
        osq << "\n";
      }

      osq.close();

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "TetrahedralityParam: unable to open %s\n",
               (getOutputFileName() + "q").c_str());
      painCave.isFatal = 1;
      simError();
    }

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();

    if (nFrames == 1) {
      std::vector<StuntDouble*>::iterator iter;
      std::ofstream osd((getOutputFileName() + "dxyz").c_str());

      if (osd.is_open()) {
        osd << Distorted_.size() << "\n";
        osd << "\n";

        for (iter = Distorted_.begin(); iter != Distorted_.end(); ++iter) {
          Vector3d position;
          position = (*iter)->getPos();
          osd << "O  "
              << "\t";
          for (unsigned int z = 0; z < position.size(); z++) {
            osd << position[z] << "  "
                << "\t";
          }
          osd << "\n";
        }
        osd.close();
      }

      std::ofstream ost((getOutputFileName() + "txyz").c_str());

      if (ost.is_open()) {
        ost << Tetrahedral_.size() << "\n";
        ost << "\n";

        for (iter = Tetrahedral_.begin(); iter != Tetrahedral_.end(); ++iter) {
          Vector3d position;
          position = (*iter)->getPos();

          ost << "O  "
              << "\t";

          for (unsigned int z = 0; z < position.size(); z++) {
            ost << position[z] << "  "
                << "\t";
          }
          ost << "\n";
        }
        ost.close();
      }
    }
  }
}  // namespace OpenMD
