/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4] Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [4] , Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011). *
 */

#include <algorithm>
#include <fstream>
#include "applications/staticProps/NitrileFrequencyMap.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "brains/Thermo.hpp"

namespace OpenMD {
  
  NitrileFrequencyMap::NitrileFrequencyMap(SimInfo* info, 
                                           const std::string& filename, 
                                           const std::string& sele1, 
                                           int nbins)
    : StaticAnalyser(info, filename), selectionScript1_(sele1), 
      evaluator1_(info), seleMan1_(info), nBins_(nbins), info_(info) {
    
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
    
    frequencyMap_["CN"] = 0.0801;
    frequencyMap_["NC"] = 0.00521;
    frequencyMap_["RCHar3"] = -0.00182;
    frequencyMap_["SigmaN"] = 0.00157;
    frequencyMap_["PiN"] = -0.00167;
    frequencyMap_["PiC"] = -0.00896;

    ForceField* forceField_ = info_->getForceField();
    set<AtomType*> atypes = info_->getSimulatedAtomTypes();
    PairList* excludes = info_->getExcludedInteractions();
    int nAtoms = info->getSnapshotManager()->getCurrentSnapshot()->getNumberOfAtoms();

    RealType rcut;
    if (info_->getSimParams()->haveCutoffRadius()) {
      rcut = info_->getSimParams()->getCutoffRadius();
    } else {
      rcut = 12.0;
    }

    EF_ = V3Zero;

    if (info_->getSimParams()->haveElectricField())
      EF_ = info_->getSimParams()->getElectricField();

    excludesForAtom.clear();
    excludesForAtom.resize(nAtoms);

    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < nAtoms; j++) {
        if (excludes->hasPair(i, j)) 
          excludesForAtom[i].push_back(j);              
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
      if ( (*i) == atom2 ) return true;
    }

    return false;
  }

  void NitrileFrequencyMap::process() {
    Molecule* mol;
    RigidBody* rb;
    Atom* atom;
    AtomType* atype;
    SimInfo::MoleculeIterator mi;
    Molecule::RigidBodyIterator rbIter;
    Molecule::AtomIterator ai2;
    Atom* atom2;
    StuntDouble* sd1;
    int ii, sdID, molID, sdID2;
    RealType li;
    RealType sPot, s1, s2;
    RealType freqShift;
    std::string name;
    map<string,RealType>::iterator fi;
    bool excluded;
    const RealType chrgToKcal = 23.0609;

    DumpReader reader(info_, dumpFilename_);    
    int nFrames = reader.getNFrames();

    nProcessed_ = nFrames/step_;

    std::fill(histogram_.begin(), histogram_.end(), 0.0);
    std::fill(count_.begin(), count_.end(), 0);
    
    for (int istep = 0; istep < nFrames; istep += step_) {
      reader.readFrame(istep);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();

      std::fill(freqs_.begin(), freqs_.end(), 0.0);

      for (mol = info_->beginMolecule(mi); mol != NULL; 
           mol = info_->nextMolecule(mi)) {
        //change the positions of atoms which belong to the rigidbodies
        for (rb = mol->beginRigidBody(rbIter); rb != NULL; 
             rb = mol->nextRigidBody(rbIter)) {
          rb->updateAtoms();
        }
      }     
      
      if  (evaluator1_.isDynamic()) {
        seleMan1_.setSelectionSet(evaluator1_.evaluate());
      }
     
      for (sd1 = seleMan1_.beginSelected(ii);
           sd1 != NULL;
           sd1 = seleMan1_.nextSelected(ii)) {
        
        sdID = sd1->getGlobalIndex();
        molID = info_->getGlobalMolMembership(sdID);
        mol = info_->getMoleculeByGlobalIndex(molID);

        Vector3d CNcentroid = mol->getRigidBodyAt(2)->getPos();
        Vector3d ra = sd1->getPos();

        atom = dynamic_cast<Atom *>(sd1);
        atype = atom->getAtomType();
        name = atype->getName();
        fi = frequencyMap_.find(name);
        if ( fi != frequencyMap_.end() ) {
          li = (*fi).second;
        } else {
          // throw error
          sprintf( painCave.errMsg,
                   "NitrileFrequencyMap::process: Unknown atype requested.\n"
                   "\t(Selection specified %s .)\n",
                   name.c_str() );
          painCave.isFatal = 1;
          simError();
        }

        sPot = sd1->getSitePotential();
        
        // Subtract out the contribution from every other site on this molecule:
        for(atom2 = mol->beginAtom(ai2); atom2 != NULL; 
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

        sPot += dot(EF_, ra - CNcentroid) * chrgToKcal ;

        freqShift = sPot * li;

        // convert the kcal/mol energies to wavenumbers:
        freqShift *= 349.757;

        freqs_[molID] += freqShift;
      }

      for (int i = 0; i < info_->getNGlobalMolecules(); ++i) {
        int binNo = int(nBins_ * (freqs_[i] - minFreq_)/(maxFreq_-minFreq_));

        count_[binNo]++;
      }
    }
      
    processHistogram();
    writeProbs();
    
  }
  
  void NitrileFrequencyMap::processHistogram() {
    
    int atot = 0;
    for(unsigned int i = 0; i < count_.size(); ++i) 
      atot += count_[i];
    
    for(unsigned int i = 0; i < count_.size(); ++i) {
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
        RealType freq = minFreq_ + (RealType)(i)*(maxFreq_-minFreq_) /
          (RealType)histogram_.size();
        rdfStream << freq << "\t" << histogram_[i] << "\n";
      }
      
    } else {
      
      sprintf(painCave.errMsg, "NitrileFrequencyMap: unable to open %s\n", 
              outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }
    
    rdfStream.close();
  }
  
}

