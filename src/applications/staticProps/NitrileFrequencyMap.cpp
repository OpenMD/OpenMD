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

#include "applications/staticProps/NitrileFrequencyMap.hpp"

#include <algorithm>
#include <fstream>

#include "brains/Thermo.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

namespace OpenMD {

  NitrileFrequencyMap::NitrileFrequencyMap(SimInfo* info,
                                           const std::string& filename,
                                           const std::string& sele1,
                                           int nbins) :
      StaticAnalyser(info, filename, nbins),
      info_(info), selectionScript1_(sele1), seleMan1_(info_),
      evaluator1_(info_) {
    setOutputName(getPrefix(filename) + ".freqs");

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }

    count_.resize(nBins_);
    histogram_.resize(nBins_);

    freqs_.resize(info_->getNGlobalMolecules());

    minFreq_ = -50;
    maxFreq_ = 50;

    // Values from Choi et. al. "Nitrile and thiocyanate IR probes:
    // Quantum chemistry calculation studies and multivariate
    // least-square ﬁtting analysis," J. Chem. Phys. 128, 134506 (2008).
    //
    // These map site electrostatic potentials onto frequency shifts
    // in the same energy units that one computes the total potential.

    frequencyMap_["CN"]     = 0.0801;
    frequencyMap_["NC"]     = 0.00521;
    frequencyMap_["RCHar3"] = -0.00182;
    frequencyMap_["SigmaN"] = 0.00157;
    frequencyMap_["PiN"]    = -0.00167;
    frequencyMap_["PiC"]    = -0.00896;

    ForceField* forceField_ = info_->getForceField();
    AtomTypeSet atypes      = info_->getSimulatedAtomTypes();
    PairList* excludes      = info_->getExcludedInteractions();
    int nAtoms =
        info->getSnapshotManager()->getCurrentSnapshot()->getNumberOfAtoms();

    RealType rcut;
    if (info_->getSimParams()->haveCutoffRadius()) {
      rcut = info_->getSimParams()->getCutoffRadius();
    } else {
      rcut = 12.0;
    }

    EF_ = V3Zero;

    std::vector<RealType> ef;
    bool efSpec = false;

    if (info_->getSimParams()->haveElectricField()) {
      efSpec = true;
      ef     = info_->getSimParams()->getElectricField();
    }
    if (info_->getSimParams()->haveUniformField()) {
      efSpec = true;
      ef     = info_->getSimParams()->getUniformField();
    }
    if (efSpec) {
      if (ef.size() != 3) {
        snprintf(
            painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
            "NitrileFrequencyMap: Incorrect number of parameters specified for "
            "uniformField.\n"
            "\tthere should be 3 parameters, but %zu were specified.\n",
            ef.size());
        painCave.isFatal = 1;
        simError();
      }
      EF_.x() = ef[0];
      EF_.y() = ef[1];
      EF_.z() = ef[2];
    }

    excludesForAtom.clear();
    excludesForAtom.resize(nAtoms);

    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < nAtoms; j++) {
        if (excludes->hasPair(i, j)) excludesForAtom[i].push_back(j);
      }
    }

    electrostatic_ = new Electrostatic();
    electrostatic_->setSimInfo(info_);
    electrostatic_->setForceField(forceField_);
    electrostatic_->setSimulatedAtomTypes(atypes);
    electrostatic_->setCutoffRadius(rcut);
  }

  bool NitrileFrequencyMap::excludeAtomPair(int atom1, int atom2) {
    for (vector<int>::iterator i = excludesForAtom[atom1].begin();
         i != excludesForAtom[atom1].end(); ++i) {
      if ((*i) == atom2) return true;
    }

    return false;
  }

  void NitrileFrequencyMap::process() {
    Molecule* mol;
    Atom* atom;
    AtomType* atype;
    SimInfo::MoleculeIterator mi;
    Molecule::AtomIterator ai2;
    Atom* atom2;
    StuntDouble* sd1;
    int ii, sdID, molID, sdID2;
    RealType li(0.0);
    RealType sPot, s1, s2;
    RealType freqShift;
    std::string name;
    map<string, RealType>::iterator fi;
    bool excluded;
    const RealType chrgToKcal = 23.0609;

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();

    nProcessed_ = nFrames / step_;

    std::fill(histogram_.begin(), histogram_.end(), 0.0);
    std::fill(count_.begin(), count_.end(), 0);

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      std::fill(freqs_.begin(), freqs_.end(), 0.0);

      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }

      for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL;
           sd1 = seleMan1_.nextSelected(ii)) {
        sdID  = sd1->getGlobalIndex();
        molID = info_->getGlobalMolMembership(sdID);
        mol   = info_->getMoleculeByGlobalIndex(molID);

        Vector3d CNcentroid = mol->getRigidBodyAt(2)->getPos();
        Vector3d ra         = sd1->getPos();

        atom  = dynamic_cast<Atom*>(sd1);
        atype = atom->getAtomType();
        name  = atype->getName();
        fi    = frequencyMap_.find(name);
        if (fi != frequencyMap_.end()) {
          li = (*fi).second;
        } else {
          // throw error
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "NitrileFrequencyMap::process: Unknown atype requested.\n"
                   "\t(Selection specified %s .)\n",
                   name.c_str());
          painCave.isFatal = 1;
          simError();
        }

        sPot = sd1->getSitePotential();

        // Subtract out the contribution from every other site on this
        // molecule:
        for (atom2 = mol->beginAtom(ai2); atom2 != NULL;
             atom2 = mol->nextAtom(ai2)) {
          sdID2 = atom2->getGlobalIndex();
          if (sdID == sdID2) {
            excluded = true;
          } else {
            excluded = excludeAtomPair(sdID, sdID2);
          }

          electrostatic_->getSitePotentials(atom, atom2, excluded, s1, s2);

          sPot -= s1;
        }

        // Add the contribution from the electric field:

        sPot += dot(EF_, ra - CNcentroid) * chrgToKcal;

        freqShift = sPot * li;

        // convert the kcal/mol energies to wavenumbers:
        freqShift *= 349.757;

        freqs_[molID] += freqShift;
      }

      for (int i = 0; i < info_->getNGlobalMolecules(); ++i) {
        int binNo =
            int(nBins_ * (freqs_[i] - minFreq_) / (maxFreq_ - minFreq_));

        count_[binNo]++;
      }
    }

    processHistogram();
    writeProbs();
  }

  void NitrileFrequencyMap::processHistogram() {
    int atot = 0;
    for (unsigned int i = 0; i < count_.size(); ++i)
      atot += count_[i];

    for (unsigned int i = 0; i < count_.size(); ++i) {
      histogram_[i] = double(count_[i] / double(atot));
    }
  }

  void NitrileFrequencyMap::writeProbs() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#NitrileFrequencyMap\n";
      rdfStream << "#nFrames:\t" << nProcessed_ << "\n";
      rdfStream << "#selection1: (" << selectionScript1_ << ")";
      rdfStream << "\n";
      rdfStream << "#nu\tp(nu))\n";
      for (unsigned int i = 0; i < histogram_.size(); ++i) {
        RealType freq = minFreq_ + (RealType)(i) * (maxFreq_ - minFreq_) /
                                       (RealType)histogram_.size();
        rdfStream << freq << "\t" << histogram_[i] << "\n";
      }

    } else {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "NitrileFrequencyMap: unable to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }

}  // namespace OpenMD
