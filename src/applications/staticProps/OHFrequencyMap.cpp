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

#include "applications/staticProps/OHFrequencyMap.hpp"

#include <algorithm>
#include <fstream>

#include "brains/ForceManager.hpp"
#include "brains/SimSnapshotManager.hpp"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "utils/simError.h"

namespace OpenMD {

  OHFrequencyMap::OHFrequencyMap(SimInfo* info, const std::string& filename,
				 const std::string& sele1, int nbins) :
    StaticAnalyser(info, filename, nbins),
    selectionScript1_(sele1), seleMan1_(info_), evaluator1_(info_) {
    
    setOutputName(getPrefix(filename) + ".OHfreqs");
    
    dumpHasElectricFields_ = info_->getSimParams()->getOutputElectricField();
    
    // If we don't have electric fields stored in the dump file, we
    // need to expand storage to hold them during a force calculation.
    
    if(!dumpHasElectricFields_) {
      int atomStorageLayout        = info_->getAtomStorageLayout();
      int rigidBodyStorageLayout   = info->getRigidBodyStorageLayout();
      int cutoffGroupStorageLayout = info->getCutoffGroupStorageLayout();
      
      atomStorageLayout |= DataStorage::dslElectricField;
      rigidBodyStorageLayout |= DataStorage::dslElectricField;
      info_->setAtomStorageLayout(atomStorageLayout);
      info_->setRigidBodyStorageLayout(rigidBodyStorageLayout);
      info_->setCutoffGroupStorageLayout(cutoffGroupStorageLayout);
      
      info_->setSnapshotManager(new SimSnapshotManager(info_, atomStorageLayout,
						       rigidBodyStorageLayout,
						       cutoffGroupStorageLayout));
    }

    evaluator1_.loadScriptString(sele1);
    if (!evaluator1_.isDynamic()) {
      seleMan1_.setSelectionSet(evaluator1_.evaluate());
    }
    
    count_.resize(nBins_);
    histogram_.resize(nBins_);

    minFreq_ = 2900;
    maxFreq_ = 3900;

    // Values from "Combined electronic structure/molecular dynamics
    // approach for ultrafast infrared spectroscopy of dilute HOD in
    // liquid H2O and D2O," by S. A. Corcelli, C. P. Lawrence, and
    // J. L. Skinner, J. Chem. Phys. 120, 8107–8117 (2004),
    // https://doi.org/10.1063/1.1683072
    //
    // These map site electric fields (in atomic units) onto
    // frequencies in wavenumbers.

    frequencyMap_["H_TIP4P"]     = std::make_pair(3832.0, -12141.0);
    frequencyMap_["H_SPCE"]      = std::make_pair(3806.0, -10792.0);
    // Assuming OH frequency map for TIP4P-Ice is the same as TIP4P.
    // This may not be a great assumption.
    frequencyMap_["H_TIP4P-Ice"] = std::make_pair(3832.0, -12141.0);

    ForceField* forceField_ = info_->getForceField();
    AtomTypeSet atypes      = info_->getSimulatedAtomTypes();
    int nAtoms =
        info->getSnapshotManager()->getCurrentSnapshot()->getNumberOfAtoms();

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
		 "OHFrequencyMap: Incorrect number of parameters specified for "
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
  }

  void OHFrequencyMap::process() {
    ForceManager* forceMan;
    Molecule* mol;
    Atom* atom;
    AtomType* atype;
    SimInfo::MoleculeIterator mi;
    Atom* atom2;
    StuntDouble* sd1;
    int ii, sdID, molID;
    RealType a0(0.0), a1(0.0);
    Vector3d sE(0.0);
    RealType freq;
    std::string name;
    map<string, std::pair<RealType, RealType>>::iterator fi;
    const RealType chrgToKcal = 23.060548;
    const RealType kcalToAU = 0.0008432975573; // fields in atomic units

    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();
    nProcessed_ = nFrames / step_;

    if (!dumpHasElectricFields_) {
      // We'll need the force manager to compute the fields
      forceMan = new ForceManager(info_);
    }
    
    std::fill(histogram_.begin(), histogram_.end(), 0.0);
    std::fill(count_.begin(), count_.end(), 0);

    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      if (!dumpHasElectricFields_) {
	forceMan->setDoElectricField(true);
	forceMan->calcForces();
      }
      
      if (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }

      for (sd1 = seleMan1_.beginSelected(ii); sd1 != NULL;
           sd1 = seleMan1_.nextSelected(ii)) {
	
        sdID  = sd1->getGlobalIndex();
        molID = info_->getGlobalMolMembership(sdID);
        mol   = info_->getMoleculeByGlobalIndex(molID);

	Vector3d Opos = mol->getAtomAt(0)->getPos();
        Vector3d Hpos   = sd1->getPos();

	Vector3d rOH = Hpos - Opos;
	rOH.normalize();

	if (sd1->isAtom()) {
	  atom  = dynamic_cast<Atom*>(sd1);
	  atype = atom->getAtomType();
	  name  = atype->getName();

	  fi    = frequencyMap_.find(name);
	  if (fi != frequencyMap_.end()) {
	    a0 = (*fi).second.first;
	    a1 = (*fi).second.second;
	    
	    sE = sd1->getElectricField(); // in kcal mol^-1 angstrom^-1 e^-1
	    
	    // Add the applied electric field (in Volts / angstrom);
	    
	    sE += EF_ * chrgToKcal;
	    
	    // convert to atomic units (hartree electrons^-1 bohr^-1):
	    
	    sE *= kcalToAU;
	    
	    freq = a0 + a1 * dot(rOH, sE);
	  
	    int binNo =
	      int(nBins_ * (freq - minFreq_) / (maxFreq_ - minFreq_));
	    if (binNo >= 0 && binNo < nBins_) {
	      count_[binNo]++;
	    } else {
	      std::cerr << "freq = " << freq
			<< " is outside histogram range: (" << minFreq_ << ", "
			<< maxFreq_ << ")\n";
	    }
	  }
	}
      }
    }

    processHistogram();
    writeProbs();
  }

  void OHFrequencyMap::processHistogram() {
    int atot = 0;
    for (unsigned int i = 0; i < count_.size(); ++i)
      atot += count_[i];

    for (unsigned int i = 0; i < count_.size(); ++i) {
      histogram_[i] = double(count_[i] / double(atot));
    }
  }

  void OHFrequencyMap::writeProbs() {
    std::ofstream rdfStream(outputFilename_.c_str());
    if (rdfStream.is_open()) {
      rdfStream << "#OHFrequencyMap\n";
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
               "OHFrequencyMap: unable to open %s\n",
               outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();
    }

    rdfStream.close();
  }

}  // namespace OpenMD
