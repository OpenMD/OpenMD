/*
 * Copyright (c) 2012 The University of Notre Dame. All Rights Reserved.
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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#include "FluctuatingChargeConstraints.hpp"
#include "primitives/Molecule.hpp"

#ifdef IS_MPI
#include <mpi.h>
#endif

namespace OpenMD {

  FluctuatingChargeConstraints::FluctuatingChargeConstraints(SimInfo* info) : 
    info_(info), constrainRegions_(false), hasFlucQ_(false), initialized_(false) {
  }

  void FluctuatingChargeConstraints::initialize(){
    if(info_->usesFluctuatingCharges()){
      if(info_->getNFluctuatingCharges() > 0){
	hasFlucQ_ = true;
      }
    }
    initialized_ = true;
  }


  void FluctuatingChargeConstraints::setConstrainRegions(bool cr) {
    constrainRegions_ = cr;

    if (!initialized_) initialize();

    regionKeys_.clear();
    regionForce_.clear();
    regionCharges_.clear();

    if (constrainRegions_) {
      SimInfo::MoleculeIterator i;
      Molecule* mol;
      int reg;
      std::set<int> regions;      

      for (mol = info_->beginMolecule(i); mol != NULL; 
           mol = info_->nextMolecule(i)) {
        reg = mol->getRegion();
        if (reg >= 0) regions.insert(reg);
      }
      // resize the keys vector to the largest found value for regions.
      regionKeys_.resize( *(regions.end()) );
      int which = 0;
      for (std::set<int>::iterator r=regions.begin(); r!=regions.end(); ++r) {
	regionKeys_[ (*r) ] = which;
	which++;
      }
      regionForce_.resize( regionKeys_.size() );
      regionCharges_.resize( regionKeys_.size() );
    }
  }


  void FluctuatingChargeConstraints::applyConstraints() {
    if (!initialized_) initialize();
    if (!hasFlucQ_) return;

    SimInfo::MoleculeIterator i;
    Molecule::FluctuatingChargeIterator  j;
    Molecule* mol;
    Atom* atom;
    
    RealType totalFrc, totalMolFrc, regionFrc, constrainedFrc;
    // accumulate the total system fluctuating charge forces
    totalFrc = 0.0;
    if (constrainRegions_) {
      std::fill(regionForce_.begin(), regionForce_.end(), 0.0); 
      std::fill(regionCharges_.begin(), regionCharges_.end(), 0); 
    }

    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {

      int region = mol->getRegion();
          
      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {

        RealType frc = atom->getFlucQFrc();
        totalFrc += frc;
        if (constrainRegions_) {
	  if (region >= 0) {	    
	    regionForce_[regionKeys_[region]] += frc;
	    regionCharges_[regionKeys_[region]] += 1;
	  }
        }
      }
    }

#ifdef IS_MPI
    // in parallel, we need to add up the contributions from all
    // processors:
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &totalFrc, 1, MPI::REALTYPE, 
                              MPI::SUM);

    if (constrainRegions_) {
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &regionForce_[0], 
                                regionForce_.size(), MPI::REALTYPE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &regionCharges_[0], 
                                regionCharges_.size(), MPI::INT, MPI::SUM);
    }

#endif

    // divide by the total number of fluctuating charges:
    totalFrc /= info_->getNFluctuatingCharges();
    
    // do the same in the regions:
    if (constrainRegions_) {
      for (int i = 0; i < regionForce_.size(); ++i)  {
        regionForce_[ i ] /= regionCharges_[ i ];
      }
    }

    for (mol = info_->beginMolecule(i); mol != NULL; 
         mol = info_->nextMolecule(i)) {     
      
      if (constrainRegions_) {
        int region = mol->getRegion();
	if (region >= 0) 
	  regionFrc = regionForce_[regionKeys_[region]];
	else
	  regionFrc = 0.0;
      } else {
        regionFrc = 0.0;
      }

      totalMolFrc = 0.0;
      
      // molecular constraints can be done with a second loop.
      if (mol->constrainTotalCharge()) {
        for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
             atom = mol->nextFluctuatingCharge(j)) {
          totalMolFrc += atom->getFlucQFrc();
        }
        totalMolFrc /= mol->getNFluctuatingCharges();
      }

      for (atom = mol->beginFluctuatingCharge(j); atom != NULL;
           atom = mol->nextFluctuatingCharge(j)) {
        //constrainedFrc = atom->getFlucQFrc() - totalFrc - totalMolFrc;

        constrainedFrc = atom->getFlucQFrc() - totalMolFrc;

        if (constrainRegions_) 
          constrainedFrc -= regionFrc;
        
        atom->setFlucQFrc(constrainedFrc);
      }      
    }
  }
}
