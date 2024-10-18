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

#include "applications/staticProps/HBondGeometric.hpp"

#include <vector>

#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/Constants.hpp"
#include "utils/simError.h"

namespace OpenMD {

  HBondGeometric::HBondGeometric(SimInfo* info, const std::string& filename,
                                 const std::string& sele1,
                                 const std::string& sele2, double rCut,
                                 double thetaCut, int nbins) :
      StaticAnalyser(info, filename, nbins),
      selectionScript1_(sele1), seleMan1_(info), evaluator1_(info),
      selectionScript2_(sele2), seleMan2_(info), evaluator2_(info) {
    setOutputName(getPrefix(filename) + ".hbg");

    ff_ = info_->getForceField();

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    evaluator2_.loadScriptString(sele2);
    if (!evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }

    // Set up cutoff values:
    nBins_    = nbins;
    rCut_     = rCut;
    thetaCut_ = thetaCut;

    nHBonds_.resize(nBins_);
    nDonor_.resize(nBins_);
    nAcceptor_.resize(nBins_);

    initializeHistogram();
  }

  void HBondGeometric::initializeHistogram() {
    std::fill(nHBonds_.begin(), nHBonds_.end(), 0);
    std::fill(nDonor_.begin(), nDonor_.end(), 0);
    std::fill(nAcceptor_.begin(), nAcceptor_.end(), 0);
    nSelected_ = 0;
  }

  void HBondGeometric::process() {
    Molecule* mol1;
    Molecule* mol2;
    Molecule::HBondDonor* hbd1;
    Molecule::HBondDonor* hbd2;
    std::vector<Molecule::HBondDonor*>::iterator hbdi;
    std::vector<Molecule::HBondDonor*>::iterator hbdj;
    std::vector<Atom*>::iterator hbai;
    std::vector<Atom*>::iterator hbaj;
    Atom* hba1;
    Atom* hba2;
    Vector3d dPos;
    Vector3d aPos;
    Vector3d hPos;
    Vector3d DH;
    Vector3d DA;
    RealType DAdist, DHdist, theta, ctheta;
    int ii, jj;
    int nHB, nA, nD;

    DumpReader reader(info_, dumpFilename_);
    int nFrames   = reader.getNFrames();
    frameCounter_ = 0;

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      frameCounter_++;
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }
      if (evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
      }

      for (mol1 = seleMan1_.beginSelectedMolecule(ii); mol1 != NULL;
           mol1 = seleMan1_.nextSelectedMolecule(ii)) {
        // We're collecting statistics on the molecules in selection 1:
        nHB = 0;
        nA  = 0;
        nD  = 0;

        for (mol2 = seleMan2_.beginSelectedMolecule(jj); mol2 != NULL;
             mol2 = seleMan2_.nextSelectedMolecule(jj)) {
          // loop over the possible donors in molecule 1:
          for (hbd1 = mol1->beginHBondDonor(hbdi); hbd1 != NULL;
               hbd1 = mol1->nextHBondDonor(hbdi)) {
            dPos = hbd1->donorAtom->getPos();
            hPos = hbd1->donatedHydrogen->getPos();
            DH   = hPos - dPos;
            currentSnapshot_->wrapVector(DH);
            DHdist = DH.length();

            // loop over the possible acceptors in molecule 2:
            for (hba2 = mol2->beginHBondAcceptor(hbaj); hba2 != NULL;
                 hba2 = mol2->nextHBondAcceptor(hbaj)) {
              aPos = hba2->getPos();
              DA   = aPos - dPos;
              currentSnapshot_->wrapVector(DA);
              DAdist = DA.length();

              // Distance criteria: are the donor and acceptor atoms
              // close enough?
              if (DAdist < rCut_) {
                ctheta = dot(DH, DA) / (DHdist * DAdist);
                theta  = acos(ctheta) * 180.0 / Constants::PI;

                // Angle criteria: are the D-H and D-A and vectors close?
                if (theta < thetaCut_) {
                  // molecule 1 is a Hbond donor:
                  nHB++;
                  nD++;
                }
              }
            }
          }

          // now loop over the possible acceptors in molecule 1:
          for (hba1 = mol1->beginHBondAcceptor(hbai); hba1 != NULL;
               hba1 = mol1->nextHBondAcceptor(hbai)) {
            aPos = hba1->getPos();

            // loop over the possible donors in molecule 2:
            for (hbd2 = mol2->beginHBondDonor(hbdj); hbd2 != NULL;
                 hbd2 = mol2->nextHBondDonor(hbdj)) {
              dPos = hbd2->donorAtom->getPos();

              DA = aPos - dPos;
              currentSnapshot_->wrapVector(DA);
              DAdist = DA.length();

              // Distance criteria: are the donor and acceptor atoms
              // close enough?
              if (DAdist < rCut_) {
                hPos = hbd2->donatedHydrogen->getPos();
                DH   = hPos - dPos;
                currentSnapshot_->wrapVector(DH);
                DHdist = DH.length();
                ctheta = dot(DH, DA) / (DHdist * DAdist);
                theta  = acos(ctheta) * 180.0 / Constants::PI;
                // Angle criteria: are the D-H and D-A and vectors close?
                if (theta < thetaCut_) {
                  // molecule 1 is a Hbond acceptor:
                  nHB++;
                  nA++;
                }
              }
            }
          }
        }
        collectHistogram(nHB, nA, nD);
      }
    }
    writeHistogram();
  }

  void HBondGeometric::collectHistogram(int nHB, int nA, int nD) {
    nHBonds_[nHB] += 1;
    nAcceptor_[nA] += 1;
    nDonor_[nD] += 1;
    nSelected_++;
  }

  void HBondGeometric::writeHistogram() {
    std::ofstream osq(getOutputFileName().c_str());

    if (osq.is_open()) {
      osq << "# HydrogenBonding Statistics\n";
      osq << "# selection1: (" << selectionScript1_ << ")"
          << "\tselection2: (" << selectionScript2_ << ")\n";
      osq << "# molecules in selection1: " << nSelected_ << "\n";
      osq << "# "
             "nHBonds\tnAcceptor\tnDonor\tp(nHBonds)\tp(nAcceptor)\tp(nDonor)"
             "\n";
      // Normalize by number of frames and write it out:
      for (int i = 0; i < nBins_; ++i) {
        osq << i;
        osq << "\t" << nHBonds_[i];
        osq << "\t" << nAcceptor_[i];
        osq << "\t" << nDonor_[i];
        osq << "\t" << (RealType)(nHBonds_[i]) / nSelected_;
        osq << "\t" << (RealType)(nAcceptor_[i]) / nSelected_;
        osq << "\t" << (RealType)(nDonor_[i]) / nSelected_;
        osq << "\n";
      }
      osq.close();

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "HBondGeometric: unable to open %s\n",
               (getOutputFileName() + "q").c_str());
      painCave.isFatal = 1;
      simError();
    }
  }
}  // namespace OpenMD
