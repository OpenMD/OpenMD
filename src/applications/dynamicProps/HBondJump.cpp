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

#include "applications/dynamicProps/HBondJump.hpp"

#include <algorithm>

#include "utils/Constants.hpp"
#include "utils/Revision.hpp"

namespace OpenMD {

  HBondJump::HBondJump(SimInfo* info, const std::string& filename,
                       const std::string& sele1, const std::string& sele2,
                       double OOcut, double thetaCut, double OHcut) :
      TimeCorrFunc<RealType>(info, filename, sele1, sele2),
      OOCut_(OOcut), thetaCut_(thetaCut), OHCut_(OHcut),
      sele1_minus_common_(info), sele2_minus_common_(info), common_(info) {
    setCorrFuncType("HBondJump");
    setOutputName(getPrefix(dumpFilename_) + ".jump");

    std::stringstream params;
    params << " OOcut = " << OOCut_ << ", thetacut = " << thetaCut_
           << ", OHcut = " << OHCut_;
    const std::string paramString = params.str();
    setParameterString(paramString);

    if (!uniqueSelections_) { seleMan2_ = seleMan1_; }
    evaluator1_.loadScriptString(selectionScript1_);
    // if selection is static, we only need to evaluate it once
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
      validateSelection(seleMan1_);
    }

    if (uniqueSelections_) {
      evaluator2_.loadScriptString(selectionScript2_);
      if (!evaluator2_.isDynamic()) {
        seleMan2_.setSelectionSet(evaluator2_.evaluate());
        validateSelection(seleMan2_);
      }
    }

    if (!evaluator1_.isDynamic() && !evaluator2_.isDynamic()) {
      // If all selections are static, we can pre-set the selections:
      common_             = seleMan1_ & seleMan2_;
      sele1_minus_common_ = seleMan1_ - common_;
      sele2_minus_common_ = seleMan2_ - common_;
    }

    // nFrames_ is initialized in MultipassCorrFunc:
    GIDtoH_.resize(nFrames_);
    hydrogen_.resize(nFrames_);
    acceptor_.resize(nFrames_);
    lastAcceptor_.resize(nFrames_);
    acceptorStartFrame_.resize(nFrames_);
    selected_.resize(nFrames_);
  }

  void HBondJump::computeFrame(int istep) {
    // Map of atomic global IDs to HBond donor hydrogens:
    GIDtoH_[istep].resize(info_->getNGlobalAtoms(), -1);

    // Find all of the HBonds in this frame
    findHBonds(istep);

    if (!uniqueSelections_) { seleMan2_ = seleMan1_; }
    if (evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    if (uniqueSelections_ && evaluator2_.isDynamic()) {
      seleMan2_.setSelectionSet(evaluator2_.evaluate());
    }
    if (evaluator1_.isDynamic() || evaluator2_.isDynamic()) {
      common_             = seleMan1_ & seleMan2_;
      sele1_minus_common_ = seleMan1_ - common_;
      sele2_minus_common_ = seleMan2_ - common_;
    }

    // Label the found HBonds as selected:
    processNonOverlapping(istep, sele1_minus_common_, seleMan2_);
    processNonOverlapping(istep, common_, sele2_minus_common_);
    processOverlapping(istep, common_);
  }

  void HBondJump::correlation() {
    std::vector<int> s1;
    std::vector<int>::iterator i1;
    RealType corrVal;
    int index1, index2, count, gid, aInd1, aInd2;

    for (int i = 0; i < nFrames_; ++i) {
      RealType time1 = times_[i];
      s1             = hydrogen_[i];

      for (int j = i; j < nFrames_; ++j) {
        // Perform a sanity check on the actual configuration times to
        // make sure the configurations are spaced the same amount the
        // sample time said they were spaced:

        RealType time2 = times_[j];

        if (fabs((time2 - time1) - (j - i) * deltaTime_) > 1.0e-4) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "HBondJump::correlation Error: sampleTime (%f)\n"
                   "\tin %s does not match actual time-spacing between\n"
                   "\tconfigurations %d (t = %f) and %d (t = %f).\n",
                   deltaTime_, dumpFilename_.c_str(), i, time1, j, time2);
          painCave.isFatal = 1;
          simError();
        }

        int timeBin = int((time2 - time1) / deltaTime_ + 0.5);

        corrVal = 0.0;
        count   = 0;

        // loop over the Hydrogens found in frame i:

        for (i1 = s1.begin(); i1 != s1.end(); ++i1) {
          // gid is the global ID of Hydrogen index1 in frame i:
          gid    = *i1;
          index1 = GIDtoH_[i][gid];

          // find matching hydrogen in frame j:
          index2 = GIDtoH_[j][gid];

          if (selected_[i][index1]) {
            count++;

            if (acceptor_[i][index1] == -1) {
              aInd1 = lastAcceptor_[i][index1];
            } else {
              aInd1 = acceptor_[i][index1];
            }

            if (acceptor_[j][index2] == -1) {
              aInd2 = lastAcceptor_[j][index2];
            } else {
              aInd2 = acceptor_[j][index2];
            }

            // aInd1 = acceptor_[i][index1];
            // aInd2 = acceptor_[j][index2];

            if (aInd1 != aInd2) {
              // different acceptor so nA(0) . nB(t) = 1
              corrVal += 1;
            } else {
              // same acceptor, but we need to look at the start frames
              // for these H-bonds to make sure it is the same H-bond:
              if (acceptorStartFrame_[i][index1] !=
                  acceptorStartFrame_[j][index2]) {
                // different start frame, so this is considered a
                // different H-bond:
                corrVal += 1;
              } else {
                // same start frame, so this is considered the same H-bond:
                corrVal += 0;
              }
            }
          }
        }
        histogram_[timeBin] += corrVal;
        count_[timeBin] += count;
      }
    }
  }

  void HBondJump::postCorrelate() {
    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      if (count_[i] > 0) {
        histogram_[i] /= count_[i];
      } else {
        histogram_[i] = 0;
      }
      histogram_[i] = 1.0 - histogram_[i];
    }
  }

  void HBondJump::processNonOverlapping(int frame, SelectionManager& sman1,
                                        SelectionManager& sman2) {
    Molecule* mol1;
    Molecule* mol2;
    int i, j;
    std::vector<Molecule::HBondDonor*>::iterator hbdi;
    Molecule::HBondDonor* hbd;
    std::vector<Atom*>::iterator hbai;
    Atom* hba;
    int hInd, index, aInd;

    // This is the same as a non-overlapping pairwise loop structure:
    // for (int i = 0;  i < ni ; ++i ) {
    //   for (int j = 0; j < nj; ++j) {}
    // }

    for (mol1 = sman1.beginSelectedMolecule(i); mol1 != NULL;
         mol1 = sman1.nextSelectedMolecule(i)) {
      for (mol2 = sman2.beginSelectedMolecule(j); mol2 != NULL;
           mol2 = sman2.nextSelectedMolecule(j)) {
        // loop over the possible donors in molecule 1:
        for (hbd = mol1->beginHBondDonor(hbdi); hbd != NULL;
             hbd = mol1->nextHBondDonor(hbdi)) {
          hInd  = hbd->donatedHydrogen->getGlobalIndex();
          index = GIDtoH_[frame][hInd];
          aInd  = acceptor_[frame][index];

          for (hba = mol2->beginHBondAcceptor(hbai); hba != NULL;
               hba = mol2->nextHBondAcceptor(hbai)) {
            if (hba->getGlobalIndex() == aInd) {
              selected_[frame][index] = true;
            }
          }
        }

        // loop over the possible donors in molecule 2:
        for (hbd = mol2->beginHBondDonor(hbdi); hbd != NULL;
             hbd = mol2->nextHBondDonor(hbdi)) {
          hInd  = hbd->donatedHydrogen->getGlobalIndex();
          index = GIDtoH_[frame][hInd];
          aInd  = acceptor_[frame][index];

          for (hba = mol1->beginHBondAcceptor(hbai); hba != NULL;
               hba = mol1->nextHBondAcceptor(hbai)) {
            if (hba->getGlobalIndex() == aInd) {
              selected_[frame][index] = true;
            }
          }
        }
      }
    }
  }

  void HBondJump::processOverlapping(int frame, SelectionManager& sman) {
    Molecule* mol1;
    Molecule* mol2;
    int i, j;

    std::vector<Molecule::HBondDonor*>::iterator hbdi;
    Molecule::HBondDonor* hbd;
    std::vector<Atom*>::iterator hbai;
    Atom* hba;
    int hInd, index, aInd;

    // This is the same as a pairwise loop structure:
    // for (int i = 0;  i < n-1 ; ++i ) {
    //   for (int j = i + 1; j < n; ++j) {}
    // }

    for (mol1 = sman.beginSelectedMolecule(i); mol1 != NULL;
         mol1 = sman.nextSelectedMolecule(i)) {
      // loop over the possible donors in molecule 1:
      for (hbd = mol1->beginHBondDonor(hbdi); hbd != NULL;
           hbd = mol1->nextHBondDonor(hbdi)) {
        hInd  = hbd->donatedHydrogen->getGlobalIndex();
        index = GIDtoH_[frame][hInd];
        aInd  = acceptor_[frame][index];

        for (j = i, mol2 = sman.nextSelectedMolecule(j); mol2 != NULL;
             mol2 = sman.nextSelectedMolecule(j)) {
          for (hba = mol2->beginHBondAcceptor(hbai); hba != NULL;
               hba = mol2->nextHBondAcceptor(hbai)) {
            if (hba->getGlobalIndex() == aInd) {
              selected_[frame][index] = true;
            }
          }
        }
      }

      for (hba = mol1->beginHBondAcceptor(hbai); hba != NULL;
           hba = mol1->nextHBondAcceptor(hbai)) {
        aInd = hba->getGlobalIndex();

        for (j = i, mol2 = sman.nextSelectedMolecule(j); mol2 != NULL;
             mol2 = sman.nextSelectedMolecule(j)) {
          for (hbd = mol2->beginHBondDonor(hbdi); hbd != NULL;
               hbd = mol2->nextHBondDonor(hbdi)) {
            hInd  = hbd->donatedHydrogen->getGlobalIndex();
            index = GIDtoH_[frame][hInd];

            if (acceptor_[frame][index] == aInd) {
              selected_[frame][index] = true;
            }
          }
        }
      }
    }
  }

  void HBondJump::findHBonds(int frame) {
    Molecule* mol1;
    Molecule* mol2;
    SimInfo::MoleculeIterator mi, mj;
    std::vector<Molecule::HBondDonor*>::iterator hbdi;
    Molecule::HBondDonor* hbd;
    std::vector<Atom*>::iterator hbai;
    Atom* hba;
    Vector3d dPos, hPos, aPos;
    int hInd, index, aInd;

    // Register all the possible HBond donor hydrogens:
    for (mol1 = info_->beginMolecule(mi); mol1 != NULL;
         mol1 = info_->nextMolecule(mi)) {
      for (hbd = mol1->beginHBondDonor(hbdi); hbd != NULL;
           hbd = mol1->nextHBondDonor(hbdi)) {
        hInd  = hbd->donatedHydrogen->getGlobalIndex();
        index = registerHydrogen(frame, hInd);
      }
    }

    for (mol1 = info_->beginMolecule(mi); mol1 != NULL;
         mol1 = info_->nextMolecule(mi)) {
      for (hbd = mol1->beginHBondDonor(hbdi); hbd != NULL;
           hbd = mol1->nextHBondDonor(hbdi)) {
        hInd  = hbd->donatedHydrogen->getGlobalIndex();
        index = GIDtoH_[frame][hInd];

        dPos = hbd->donorAtom->getPos();
        hPos = hbd->donatedHydrogen->getPos();

        for (mj = mi, mol2 = info_->beginMolecule(mj); mol2 != NULL;
             mol2 = info_->nextMolecule(mj)) {
          for (hba = mol2->beginHBondAcceptor(hbai); hba != NULL;
               hba = mol2->nextHBondAcceptor(hbai)) {
            aPos = hba->getPos();

            if (isHBond(dPos, hPos, aPos)) {
              aInd = hba->getGlobalIndex();
              registerHydrogenBond(frame, index, hInd, aInd);
            }
          }
        }
      }

      for (hba = mol1->beginHBondAcceptor(hbai); hba != NULL;
           hba = mol1->nextHBondAcceptor(hbai)) {
        aPos = hba->getPos();

        for (mj = mi, mol2 = info_->beginMolecule(mj); mol2 != NULL;
             mol2 = info_->nextMolecule(mj)) {
          for (hbd = mol2->beginHBondDonor(hbdi); hbd != NULL;
               hbd = mol2->nextHBondDonor(hbdi)) {
            hInd = hbd->donatedHydrogen->getGlobalIndex();
            // no need to register, just look up the index:
            index = GIDtoH_[frame][hInd];

            dPos = hbd->donorAtom->getPos();
            hPos = hbd->donatedHydrogen->getPos();

            if (isHBond(dPos, hPos, aPos)) {
              aInd = hba->getGlobalIndex();
              registerHydrogenBond(frame, index, hInd, aInd);
            }
          }
        }
      }
    }
  }

  bool HBondJump::isHBond(Vector3d donorPos, Vector3d hydrogenPos,
                          Vector3d acceptorPos) {
    Vector3d DA = acceptorPos - donorPos;
    currentSnapshot_->wrapVector(DA);
    RealType DAdist = DA.length();

    // Distance criteria: are the donor and acceptor atoms
    // close enough?
    if (DAdist < OOCut_) {
      Vector3d DH = hydrogenPos - donorPos;
      currentSnapshot_->wrapVector(DH);
      RealType DHdist = DH.length();

      Vector3d HA = acceptorPos - hydrogenPos;
      currentSnapshot_->wrapVector(HA);
      RealType HAdist = HA.length();

      RealType ctheta = dot(DH, DA) / (DHdist * DAdist);
      RealType theta  = acos(ctheta) * 180.0 / Constants::PI;

      // Angle criteria: are the D-H and D-A and vectors close?
      if (theta < thetaCut_ && HAdist < OHCut_) { return true; }
    }
    return false;
  }

  int HBondJump::registerHydrogen(int frame, int hIndex) {
    int index;

    // If this hydrogen wasn't already registered, register it:
    if (GIDtoH_[frame][hIndex] == -1) {
      index                  = hydrogen_[frame].size();
      GIDtoH_[frame][hIndex] = index;
      hydrogen_[frame].push_back(hIndex);
      acceptor_[frame].push_back(-1);
      selected_[frame].push_back(false);

      if (frame == 0) {
        lastAcceptor_[frame].push_back(-1);
        acceptorStartFrame_[frame].push_back(frame);
      } else {
        // Copy the last acceptor.
        int prevIndex = GIDtoH_[frame - 1][hIndex];
        lastAcceptor_[frame].push_back(lastAcceptor_[frame - 1][prevIndex]);
        acceptorStartFrame_[frame].push_back(
            acceptorStartFrame_[frame - 1][prevIndex]);
      }
    } else {
      // This hydrogen was already registered.  Just return the index:
      index = GIDtoH_[frame][hIndex];
    }
    return index;
  }

  void HBondJump::registerHydrogenBond(int frame, int index, int hIndex,
                                       int acceptorIndex) {
    acceptor_[frame][index]     = acceptorIndex;
    lastAcceptor_[frame][index] = acceptorIndex;

    if (frame == 0) {
      acceptorStartFrame_[frame][index] = frame;
    } else {
      int prevIndex = GIDtoH_[frame - 1][hIndex];
      if (acceptorIndex != lastAcceptor_[frame - 1][prevIndex]) {
        acceptorStartFrame_[frame][index] = frame;
      } else {
        acceptorStartFrame_[frame][index] =
            acceptorStartFrame_[frame - 1][prevIndex];
      }
    }
  }

  HBondJumpZ::HBondJumpZ(SimInfo* info, const std::string& filename,
                         const std::string& sele1, const std::string& sele2,
                         double OOcut, double thetaCut, double OHcut,
                         int nZBins, int axis) :
      HBondJump(info, filename, sele1, sele2, OOcut, thetaCut, OHcut),
      nZBins_(nZBins), axis_(axis) {
    setCorrFuncType("HBondJumpZ");
    setOutputName(getPrefix(dumpFilename_) + ".jumpZ");

    switch (axis_) {
    case 0:
      axisLabel_ = "x";
      break;
    case 1:
      axisLabel_ = "y";
      break;
    case 2:
    default:
      axisLabel_ = "z";
      break;
    }

    zbin_.resize(nFrames_);
    histogram_.resize(nTimeBins_);
    counts_.resize(nTimeBins_);
    for (unsigned int i = 0; i < nTimeBins_; i++) {
      histogram_[i].resize(nZBins_);
      std::fill(histogram_[i].begin(), histogram_[i].end(), 0.0);
      counts_[i].resize(nZBins_);
      std::fill(counts_[i].begin(), counts_[i].end(), 0);
    }
  }

  void HBondJumpZ::findHBonds(int frame) {
    Molecule* mol1;
    Molecule* mol2;
    SimInfo::MoleculeIterator mi, mj;
    std::vector<Molecule::HBondDonor*>::iterator hbdi;
    Molecule::HBondDonor* hbd;
    std::vector<Atom*>::iterator hbai;
    Atom* hba;
    Vector3d dPos, hPos, aPos, pos;
    int hInd, index, aInd, zBin;

    Mat3x3d hmat       = currentSnapshot_->getHmat();
    RealType halfBoxZ_ = hmat(axis_, axis_) / 2.0;

    // Register all the possible HBond donor hydrogens:
    for (mol1 = info_->beginMolecule(mi); mol1 != NULL;
         mol1 = info_->nextMolecule(mi)) {
      for (hbd = mol1->beginHBondDonor(hbdi); hbd != NULL;
           hbd = mol1->nextHBondDonor(hbdi)) {
        hInd  = hbd->donatedHydrogen->getGlobalIndex();
        index = registerHydrogen(frame, hInd);
      }
    }

    for (mol1 = info_->beginMolecule(mi); mol1 != NULL;
         mol1 = info_->nextMolecule(mi)) {
      for (hbd = mol1->beginHBondDonor(hbdi); hbd != NULL;
           hbd = mol1->nextHBondDonor(hbdi)) {
        hInd  = hbd->donatedHydrogen->getGlobalIndex();
        index = GIDtoH_[frame][hInd];

        dPos = hbd->donorAtom->getPos();
        hPos = hbd->donatedHydrogen->getPos();

        for (mj = mi, mol2 = info_->beginMolecule(mj); mol2 != NULL;
             mol2 = info_->nextMolecule(mj)) {
          for (hba = mol2->beginHBondAcceptor(hbai); hba != NULL;
               hba = mol2->nextHBondAcceptor(hbai)) {
            aPos = hba->getPos();

            if (isHBond(dPos, hPos, aPos)) {
              aInd = hba->getGlobalIndex();
              registerHydrogenBond(frame, index, hInd, aInd);
              pos = hPos;
              if (info_->getSimParams()->getUsePeriodicBoundaryConditions())
                currentSnapshot_->wrapVector(pos);
              zBin =
                  int(nZBins_ * (halfBoxZ_ + pos[axis_]) / hmat(axis_, axis_));
              zbin_[frame][index] = zBin;
            }
          }
        }
      }

      for (hba = mol1->beginHBondAcceptor(hbai); hba != NULL;
           hba = mol1->nextHBondAcceptor(hbai)) {
        aPos = hba->getPos();

        for (mj = mi, mol2 = info_->beginMolecule(mj); mol2 != NULL;
             mol2 = info_->nextMolecule(mj)) {
          for (hbd = mol2->beginHBondDonor(hbdi); hbd != NULL;
               hbd = mol2->nextHBondDonor(hbdi)) {
            hInd = hbd->donatedHydrogen->getGlobalIndex();
            // no need to register, just look up the index:
            index = GIDtoH_[frame][hInd];

            dPos = hbd->donorAtom->getPos();
            hPos = hbd->donatedHydrogen->getPos();

            if (isHBond(dPos, hPos, aPos)) {
              aInd = hba->getGlobalIndex();
              registerHydrogenBond(frame, index, hInd, aInd);
              pos = hPos;
              if (info_->getSimParams()->getUsePeriodicBoundaryConditions())
                currentSnapshot_->wrapVector(pos);
              zBin =
                  int(nZBins_ * (halfBoxZ_ + pos[axis_]) / hmat(axis_, axis_));
              zbin_[frame][index] = zBin;
            }
          }
        }
      }
    }
  }

  int HBondJumpZ::registerHydrogen(int frame, int hIndex) {
    int index;

    // If this hydrogen wasn't already registered, register it:
    if (GIDtoH_[frame][hIndex] == -1) {
      index                  = hydrogen_[frame].size();
      GIDtoH_[frame][hIndex] = index;
      hydrogen_[frame].push_back(hIndex);
      acceptor_[frame].push_back(-1);
      selected_[frame].push_back(false);
      zbin_[frame].push_back(-1);

      if (frame == 0) {
        lastAcceptor_[frame].push_back(-1);
        acceptorStartFrame_[frame].push_back(frame);
      } else {
        // Copy the last acceptor.
        int prevIndex = GIDtoH_[frame - 1][hIndex];
        lastAcceptor_[frame].push_back(lastAcceptor_[frame - 1][prevIndex]);
        acceptorStartFrame_[frame].push_back(
            acceptorStartFrame_[frame - 1][prevIndex]);
      }
    } else {
      // This hydrogen was already registered.  Just return the index:
      index = GIDtoH_[frame][hIndex];
    }
    return index;
  }

  void HBondJumpZ::correlation() {
    std::vector<int> s1;
    std::vector<int>::iterator i1;
    int index1, index2, gid, aInd1, aInd2, zBin;

    for (int i = 0; i < nFrames_; ++i) {
      RealType time1 = times_[i];
      s1             = hydrogen_[i];

      for (int j = i; j < nFrames_; ++j) {
        // Perform a sanity check on the actual configuration times to
        // make sure the configurations are spaced the same amount the
        // sample time said they were spaced:

        RealType time2 = times_[j];

        if (fabs((time2 - time1) - (j - i) * deltaTime_) > 1.0e-4) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "HBondJump::correlation Error: sampleTime (%f)\n"
                   "\tin %s does not match actual time-spacing between\n"
                   "\tconfigurations %d (t = %f) and %d (t = %f).\n",
                   deltaTime_, dumpFilename_.c_str(), i, time1, j, time2);
          painCave.isFatal = 1;
          simError();
        }

        int timeBin = int((time2 - time1) / deltaTime_ + 0.5);

        // loop over the Hydrogens found in frame i:

        for (i1 = s1.begin(); i1 != s1.end(); ++i1) {
          // gid is the global ID of Hydrogen index1 in frame i:
          gid    = *i1;
          index1 = GIDtoH_[i][gid];

          // find matching hydrogen in frame j:
          index2 = GIDtoH_[j][gid];

          if (selected_[i][index1]) {
            zBin = zbin_[i][index1];
            counts_[timeBin][zBin]++;

            if (acceptor_[i][index1] == -1) {
              aInd1 = lastAcceptor_[i][index1];
            } else {
              aInd1 = acceptor_[i][index1];
            }

            if (acceptor_[j][index2] == -1) {
              aInd2 = lastAcceptor_[j][index2];
            } else {
              aInd2 = acceptor_[j][index2];
            }

            // aInd1 = acceptor_[i][index1];
            // aInd2 = acceptor_[j][index2];

            if (aInd1 != aInd2) {
              // different acceptor so nA(0) . nB(t) = 1
              histogram_[timeBin][zBin] += 1;
            } else {
              // same acceptor, but we need to look at the start frames
              // for these H-bonds to make sure it is the same H-bond:
              if (acceptorStartFrame_[i][index1] !=
                  acceptorStartFrame_[j][index2]) {
                // different start frame, so this is considered a
                // different H-bond:
                histogram_[timeBin][zBin] += 1;
              } else {
                // same start frame, so this is considered the same H-bond:
                histogram_[timeBin][zBin] += 0;
              }
            }
          }
        }
      }
    }
  }

  void HBondJumpZ::postCorrelate() {
    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      for (unsigned int j = 0; j < nZBins_; ++j) {
        if (counts_[i][j] > 0) {
          histogram_[i][j] /= counts_[i][j];
        } else {
          histogram_[i][j] = 0;
        }
        histogram_[i][j] = 1.0 - histogram_[i][j];
      }
    }
  }

  void HBondJumpZ::writeCorrelate() {
    ofstream ofs(outputFilename_.c_str());

    if (ofs.is_open()) {
      Revision r;

      ofs << "# " << getCorrFuncType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      ofs << "# privilegedAxis computed as " << axisLabel_ << " axis \n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";
      ofs << "#time\tcorrVal\n";

      for (unsigned int i = 0; i < nTimeBins_; ++i) {
        ofs << times_[i] - times_[0];

        for (unsigned int j = 0; j < nZBins_; ++j) {
          ofs << "\t" << histogram_[i][j];
        }
        ofs << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "HBondJumpZ::writeCorrelate Error: fail to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }

  HBondJumpR::HBondJumpR(SimInfo* info, const std::string& filename,
                         const std::string& sele1, const std::string& sele2,
                         const std::string& sele3, RealType OOcut,
                         RealType thetaCut, RealType OHcut, RealType len,
                         int nRBins) :
      HBondJump(info, filename, sele1, sele2, OOcut, thetaCut, OHcut),
      len_(len), nRBins_(nRBins), seleMan3_(info_), selectionScript3_(sele3),
      evaluator3_(info_) {
    setCorrFuncType("HBondJumpR");
    setOutputName(getPrefix(dumpFilename_) + ".jumpR");

    rbin_.resize(nFrames_);
    histogram_.resize(nTimeBins_);
    counts_.resize(nTimeBins_);
    for (unsigned int i = 0; i < nTimeBins_; i++) {
      histogram_[i].resize(nRBins_);
      counts_[i].resize(nRBins_);
    }

    evaluator3_.loadScriptString(selectionScript3_);
    if (!evaluator3_.isDynamic()) {
      seleMan3_.setSelectionSet(evaluator3_.evaluate());
    }

    deltaR_ = len_ / nRBins_;
  }

  void HBondJumpR::findHBonds(int frame) {
    Molecule* mol1;
    Molecule* mol2;
    SimInfo::MoleculeIterator mi, mj;
    std::vector<Molecule::HBondDonor*>::iterator hbdi;
    Molecule::HBondDonor* hbd;
    std::vector<Atom*>::iterator hbai;
    Atom* hba;
    Vector3d dPos, hPos, aPos, pos;
    int hInd, index, aInd;
    StuntDouble* sd3;
    Vector3d vec;
    RealType r;
    int isd3;

    if (evaluator3_.isDynamic()) {
      seleMan3_.setSelectionSet(evaluator3_.evaluate());
    }

    bool usePeriodicBoundaryConditions_ =
        info_->getSimParams()->getUsePeriodicBoundaryConditions();

    // Register all the possible HBond donor hydrogens:
    for (mol1 = info_->beginMolecule(mi); mol1 != NULL;
         mol1 = info_->nextMolecule(mi)) {
      for (hbd = mol1->beginHBondDonor(hbdi); hbd != NULL;
           hbd = mol1->nextHBondDonor(hbdi)) {
        hInd  = hbd->donatedHydrogen->getGlobalIndex();
        index = registerHydrogen(frame, hInd);
      }
    }

    for (mol1 = info_->beginMolecule(mi); mol1 != NULL;
         mol1 = info_->nextMolecule(mi)) {
      for (hbd = mol1->beginHBondDonor(hbdi); hbd != NULL;
           hbd = mol1->nextHBondDonor(hbdi)) {
        hInd  = hbd->donatedHydrogen->getGlobalIndex();
        index = GIDtoH_[frame][hInd];

        dPos = hbd->donorAtom->getPos();
        hPos = hbd->donatedHydrogen->getPos();

        for (mj = mi, mol2 = info_->beginMolecule(mj); mol2 != NULL;
             mol2 = info_->nextMolecule(mj)) {
          for (hba = mol2->beginHBondAcceptor(hbai); hba != NULL;
               hba = mol2->nextHBondAcceptor(hbai)) {
            aPos = hba->getPos();

            if (isHBond(dPos, hPos, aPos)) {
              aInd = hba->getGlobalIndex();
              registerHydrogenBond(frame, index, hInd, aInd);
              pos = hPos;

              RealType shortest = HONKING_LARGE_VALUE;

              // loop over selection 3 to find closest atom in selection 3:
              for (sd3 = seleMan3_.beginSelected(isd3); sd3 != NULL;
                   sd3 = seleMan3_.nextSelected(isd3)) {
                vec = pos - sd3->getPos();

                if (usePeriodicBoundaryConditions_)
                  currentSnapshot_->wrapVector(vec);

                r = vec.length();

                if (r < shortest) shortest = r;
              }

              int whichBin = int(shortest / deltaR_);
              if (whichBin < int(nRBins_)) { rbin_[frame][index] = whichBin; }
            }
          }
        }
      }

      for (hba = mol1->beginHBondAcceptor(hbai); hba != NULL;
           hba = mol1->nextHBondAcceptor(hbai)) {
        aPos = hba->getPos();

        for (mj = mi, mol2 = info_->beginMolecule(mj); mol2 != NULL;
             mol2 = info_->nextMolecule(mj)) {
          for (hbd = mol2->beginHBondDonor(hbdi); hbd != NULL;
               hbd = mol2->nextHBondDonor(hbdi)) {
            hInd = hbd->donatedHydrogen->getGlobalIndex();
            // no need to register, just look up the index:
            index = GIDtoH_[frame][hInd];

            dPos = hbd->donorAtom->getPos();
            hPos = hbd->donatedHydrogen->getPos();

            if (isHBond(dPos, hPos, aPos)) {
              aInd = hba->getGlobalIndex();
              registerHydrogenBond(frame, index, hInd, aInd);
              pos               = hPos;
              RealType shortest = HONKING_LARGE_VALUE;

              // loop over selection 3 to find closest atom in selection 3:
              for (sd3 = seleMan3_.beginSelected(isd3); sd3 != NULL;
                   sd3 = seleMan3_.nextSelected(isd3)) {
                vec = pos - sd3->getPos();

                if (usePeriodicBoundaryConditions_)
                  currentSnapshot_->wrapVector(vec);

                r = vec.length();

                if (r < shortest) shortest = r;
              }

              int whichBin = int(shortest / deltaR_);
              if (whichBin < int(nRBins_)) { rbin_[frame][index] = whichBin; }
            }
          }
        }
      }
    }
  }

  int HBondJumpR::registerHydrogen(int frame, int hIndex) {
    int index;

    // If this hydrogen wasn't already registered, register it:
    if (GIDtoH_[frame][hIndex] == -1) {
      index = hydrogen_[frame].size();

      GIDtoH_[frame][hIndex] = index;
      hydrogen_[frame].push_back(hIndex);
      acceptor_[frame].push_back(-1);
      selected_[frame].push_back(false);
      rbin_[frame].push_back(-1);

      if (frame == 0) {
        lastAcceptor_[frame].push_back(-1);
        acceptorStartFrame_[frame].push_back(frame);
      } else {
        // Copy the last acceptor.
        int prevIndex = GIDtoH_[frame - 1][hIndex];
        lastAcceptor_[frame].push_back(lastAcceptor_[frame - 1][prevIndex]);
        acceptorStartFrame_[frame].push_back(
            acceptorStartFrame_[frame - 1][prevIndex]);
      }
    } else {
      // This hydrogen was already registered.  Just return the index:
      index = GIDtoH_[frame][hIndex];
    }
    return index;
  }

  void HBondJumpR::correlation() {
    std::vector<int> s1;
    std::vector<int>::iterator i1;
    int index1, index2, gid, aInd1, aInd2;

    for (int i = 0; i < nFrames_; ++i) {
      RealType time1 = times_[i];
      s1             = hydrogen_[i];

      for (int j = i; j < nFrames_; ++j) {
        // Perform a sanity check on the actual configuration times to
        // make sure the configurations are spaced the same amount the
        // sample time said they were spaced:

        RealType time2 = times_[j];

        if (fabs((time2 - time1) - (j - i) * deltaTime_) > 1.0e-4) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "HBondJump::correlation Error: sampleTime (%f)\n"
                   "\tin %s does not match actual time-spacing between\n"
                   "\tconfigurations %d (t = %f) and %d (t = %f).\n",
                   deltaTime_, dumpFilename_.c_str(), i, time1, j, time2);
          painCave.isFatal = 1;
          simError();
        }

        int timeBin = int((time2 - time1) / deltaTime_ + 0.5);

        // loop over the Hydrogens found in frame i:

        for (i1 = s1.begin(); i1 != s1.end(); ++i1) {
          // gid is the global ID of Hydrogen index1 in frame i:
          gid    = *i1;
          index1 = GIDtoH_[i][gid];

          // find matching hydrogen in frame j:
          index2 = GIDtoH_[j][gid];

          if (selected_[i][index1]) {
            if (int rBin {rbin_[i][index1]}; rBin != -1) {
              counts_[timeBin][rBin]++;

              if (acceptor_[i][index1] == -1) {
                aInd1 = lastAcceptor_[i][index1];
              } else {
                aInd1 = acceptor_[i][index1];
              }

              if (acceptor_[j][index2] == -1) {
                aInd2 = lastAcceptor_[j][index2];
              } else {
                aInd2 = acceptor_[j][index2];
              }

              // aInd1 = acceptor_[i][index1];
              // aInd2 = acceptor_[j][index2];

              if (aInd1 != aInd2) {
                // different acceptor so nA(0) . nB(t) = 1
                histogram_[timeBin][rBin] += 1;
              } else {
                // same acceptor, but we need to look at the start frames
                // for these H-bonds to make sure it is the same H-bond:
                if (acceptorStartFrame_[i][index1] !=
                    acceptorStartFrame_[j][index2]) {
                  // different start frame, so this is considered a
                  // different H-bond:
                  histogram_[timeBin][rBin] += 1;
                } else {
                  // same start frame, so this is considered the same H-bond:
                  histogram_[timeBin][rBin] += 0;
                }
              }
            }
          }
        }
      }
    }
  }

  void HBondJumpR::postCorrelate() {
    for (unsigned int i = 0; i < nTimeBins_; ++i) {
      for (unsigned int j = 0; j < nRBins_; ++j) {
        if (counts_[i][j] > 0) {
          histogram_[i][j] /= counts_[i][j];
        } else {
          histogram_[i][j] = 0;
        }
        histogram_[i][j] = 1.0 - histogram_[i][j];
      }
    }
  }

  void HBondJumpR::writeCorrelate() {
    ofstream ofs(outputFilename_.c_str());

    if (ofs.is_open()) {
      Revision r;

      ofs << "# " << getCorrFuncType() << "\n";
      ofs << "# OpenMD " << r.getFullRevision() << "\n";
      ofs << "# " << r.getBuildDate() << "\n";
      ofs << "# selection script1: \"" << selectionScript1_;
      ofs << "\"\tselection script2: \"" << selectionScript2_ << "\"\n";
      if (!paramString_.empty())
        ofs << "# parameters: " << paramString_ << "\n";
      ofs << "#time\tcorrVal\n";

      for (unsigned int i = 0; i < nTimeBins_; ++i) {
        ofs << times_[i] - times_[0];

        for (unsigned int j = 0; j < nRBins_; ++j) {
          ofs << "\t" << histogram_[i][j];
        }
        ofs << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "HBondJumpR::writeCorrelate Error: fail to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    ofs.close();
  }

}  // namespace OpenMD
