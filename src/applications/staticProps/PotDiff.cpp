/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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

#include <algorithm>
#include <functional>
#include "applications/staticProps/PotDiff.hpp"
#include "brains/ForceManager.hpp"
#include "brains/SimSnapshotManager.hpp"
#include "utils/simError.h"
#include "io/DumpReader.hpp"
#include "primitives/Molecule.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"

namespace OpenMD {

  PotDiff::PotDiff(SimInfo* info, const std::string& filename, 
                   const std::string& sele)
    : StaticAnalyser(info, filename, 1), selectionScript_(sele), 
      seleMan_(info), evaluator_(info) {
    
    StuntDouble* sd;
    int i;
    
    setOutputName(getPrefix(filename) + ".potDiff");

    // The PotDiff is computed by negating the charge on the atom type
    // using fluctuating charge values.  If we don't have any
    // fluctuating charges in the simulation, we need to expand
    // storage to hold them.
    int storageLayout = info_->getStorageLayout();
    storageLayout |= DataStorage::dslFlucQPosition;
    storageLayout |= DataStorage::dslFlucQVelocity;
    storageLayout |= DataStorage::dslFlucQForce;
    info_->setStorageLayout(storageLayout);
    info_->setSnapshotManager(new SimSnapshotManager(info_, storageLayout));

    // now we have to figure out which AtomTypes to convert to fluctuating
    // charges
    evaluator_.loadScriptString(sele);    
    seleMan_.setSelectionSet(evaluator_.evaluate());
    for (sd = seleMan_.beginSelected(i); sd != NULL;
         sd = seleMan_.nextSelected(i)) {      
      AtomType* at = static_cast<Atom*>(sd)->getAtomType();
      FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(at);
      if (fqa.isFluctuatingCharge()) {
        selectionWasFlucQ_.push_back(true);
      } else {
        selectionWasFlucQ_.push_back(false);                
        // make a fictitious fluctuating charge with an unphysical
        // charge mass and slaterN, but we need to zero out the
        // electronegativity and hardness to remove the self
        // contribution:
        fqa.makeFluctuatingCharge(1.0e9, 0.0, 0.0, 1);
        sd->setFlucQPos(0.0);
      }
    }
    info_->getSnapshotManager()->advance();
  }
  
  void PotDiff::process() {
    StuntDouble* sd;
    int j;
  
    diff_.clear();
    DumpReader reader(info_, dumpFilename_);
    int nFrames = reader.getNFrames();

    // We'll need the force manager to compute the potential
    
    ForceManager* forceMan = new ForceManager(info_);

    // We'll need thermo to report the potential
    
    Thermo* thermo =  new Thermo(info_);

    for (int i = 0; i < nFrames; i += step_) {
      reader.readFrame(i);
      currentSnapshot_ = info_->getSnapshotManager()->getCurrentSnapshot();
    
      for (sd = seleMan_.beginSelected(j); sd != NULL;
           sd = seleMan_.nextSelected(j)) {
        if (!selectionWasFlucQ_[j])  {
          sd->setFlucQPos(0.0);
        }
      }
            
      forceMan->calcForces();
      RealType pot1 = thermo->getPotential();    

      if (evaluator_.isDynamic()) {
        seleMan_.setSelectionSet(evaluator_.evaluate());
      }

      for (sd = seleMan_.beginSelected(j); sd != NULL;
           sd = seleMan_.nextSelected(j)) {

        AtomType* at = static_cast<Atom*>(sd)->getAtomType();

        FixedChargeAdapter fca = FixedChargeAdapter(at);
        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(at);

        RealType charge = 0.0;
        
        if (fca.isFixedCharge()) charge += fca.getCharge();
        if (fqa.isFluctuatingCharge()) charge += sd->getFlucQPos();

        sd->setFlucQPos(-charge);
      }

      currentSnapshot_->clearDerivedProperties();
      forceMan->calcForces();
      RealType pot2 = thermo->getPotential();
      RealType diff = pot2-pot1;
      
      data_.add(diff);
      diff_.push_back(diff);
      times_.push_back(currentSnapshot_->getTime());

      info_->getSnapshotManager()->advance();
    }
   
    writeDiff();   
  }
  
  void PotDiff::writeDiff() {

    RealType mu, sigma, m95;
    std::ofstream ofs(outputFilename_.c_str(), std::ios::binary);
    if (ofs.is_open()) {

      data_.getAverage(mu);
      data_.getStdDev(sigma);
      data_.get95percentConfidenceInterval(m95);
      
      ofs << "#potDiff\n";
      ofs << "#selection: (" << selectionScript_ << ")\n";
      ofs << "# <diff> = "<< mu << "\n";
      ofs << "# StdDev = " << sigma << "\n";
      ofs << "# 95% confidence interval = " << m95 << "\n";
      ofs << "# t\tdiff[t]\n";
      for (unsigned int i = 0; i < diff_.size(); ++i) {
        ofs << times_[i] << "\t" << diff_[i] << "\n";
      }
      
    } else {
      
      sprintf(painCave.errMsg, "PotDiff: unable to open %s\n", 
	      outputFilename_.c_str());
      painCave.isFatal = 1;
      simError();  
    }
    ofs.close();
  }
  
}


