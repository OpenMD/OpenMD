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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Vardeman & Gezelter, in progress (2009).                        
 */
 
#include "nonbonded/InteractionManager.hpp"

namespace OpenMD {

  InteractionManager::InteractionManager() {

    initialized_ = false;
        
    lj_ = new LJ();
    gb_ = new GB();
    sticky_ = new Sticky();
    morse_ = new Morse();
    eam_ = new EAM();
    sc_ = new SC();
    electrostatic_ = new Electrostatic();
    maw_ = new MAW();
  }

  void InteractionManager::initialize() {
    
    ForceField* forceField_ = info_->getForceField();
    
    lj_->setForceField(forceField_);
    gb_->setForceField(forceField_);
    sticky_->setForceField(forceField_);
    eam_->setForceField(forceField_);
    sc_->setForceField(forceField_);
    morse_->setForceField(forceField_);
    electrostatic_->setSimInfo(info_);
    electrostatic_->setForceField(forceField_);
    maw_->setForceField(forceField_);

    ForceField::AtomTypeContainer* atomTypes = forceField_->getAtomTypes();
    ForceField::AtomTypeContainer::MapTypeIterator i1, i2;
    AtomType* atype1;
    AtomType* atype2;
    pair<AtomType*, AtomType*> key;
    pair<set<NonBondedInteraction*>::iterator, bool> ret;
    
    for (atype1 = atomTypes->beginType(i1); atype1 != NULL; 
         atype1 = atomTypes->nextType(i1)) {
      
      // add it to the map:
      AtomTypeProperties atp = atype1->getATP();    
      
      pair<map<int,AtomType*>::iterator,bool> ret;    
      ret = typeMap_.insert( pair<int, AtomType*>(atp.ident, atype1) );
      if (ret.second == false) {
        sprintf( painCave.errMsg,
                 "InteractionManager already had a previous entry with ident %d\n",
                 atp.ident);
        painCave.severity = OPENMD_INFO;
        painCave.isFatal = 0;
        simError();                 
      }
    }
    
    // Now, iterate over all known types and add to the interaction map:
    
    map<int, AtomType*>::iterator it1, it2;
    for (it1 = typeMap_.begin(); it1 != typeMap_.end(); ++it1) {
      atype1 = (*it1).second;

      for( it2 = typeMap_.begin(); it2 != typeMap_.end(); ++it2) {        
        atype2 = (*it2).second;
        
        bool vdwExplicit = false;
        bool metExplicit = false;
        bool hbExplicit = false;
                       
        key = make_pair(atype1, atype2);
        
        if (atype1->isLennardJones() && atype2->isLennardJones()) {
          interactions_[key].insert(lj_);
        }
        if (atype1->isElectrostatic() && atype2->isElectrostatic() ) {
          interactions_[key].insert(electrostatic_);
        }
        if (atype1->isSticky() && atype2->isSticky() ) {
          interactions_[key].insert(sticky_);
        }
        if (atype1->isStickyPower() && atype2->isStickyPower() ) {
          interactions_[key].insert(sticky_);
        }
        if (atype1->isEAM() && atype2->isEAM() ) {
          interactions_[key].insert(eam_);
        }
        if (atype1->isSC() && atype2->isSC() ) {
          interactions_[key].insert(sc_);
        }
        if (atype1->isGayBerne() && atype2->isGayBerne() ) {
          interactions_[key].insert(gb_);
        }
        if ((atype1->isGayBerne() && atype2->isLennardJones())
            || (atype1->isLennardJones() && atype2->isGayBerne())) {
          interactions_[key].insert(gb_);
        } 
        
        // look for an explicitly-set non-bonded interaction type using the 
        // two atom types.
        NonBondedInteractionType* nbiType = forceField_->getNonBondedInteractionType(atype1->getName(), atype2->getName());
        
        if (nbiType != NULL) {

          if (nbiType->isLennardJones()) {
            // We found an explicit Lennard-Jones interaction.  
            // override all other vdw entries for this pair of atom types:
            set<NonBondedInteraction*>::iterator it;
            for (it = interactions_[key].begin(); 
                 it != interactions_[key].end(); ++it) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == VANDERWAALS_FAMILY) interactions_[key].erase(*it);
            }
            interactions_[key].insert(lj_);
            vdwExplicit = true;
          }
          
          if (nbiType->isMorse()) {
            if (vdwExplicit) {
              sprintf( painCave.errMsg,
                       "InteractionManager::initialize found more than one "
                       "explicit \n"
                       "\tvan der Waals interaction for atom types %s - %s\n",
                       atype1->getName().c_str(), atype2->getName().c_str());
              painCave.severity = OPENMD_ERROR;
              painCave.isFatal = 1;
              simError();
            }
            // We found an explicit Morse interaction.  
            // override all other vdw entries for this pair of atom types:
            set<NonBondedInteraction*>::iterator it;
            for (it = interactions_[key].begin(); 
                 it != interactions_[key].end(); ++it) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == VANDERWAALS_FAMILY) interactions_[key].erase(*it);
            }
            interactions_[key].insert(morse_);
            vdwExplicit = true;
          }
          
          if (nbiType->isEAM()) {
            // We found an explicit EAM interaction.  
            // override all other metallic entries for this pair of atom types:
            set<NonBondedInteraction*>::iterator it;
            for (it = interactions_[key].begin(); 
                 it != interactions_[key].end(); ++it) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == METALLIC_FAMILY) interactions_[key].erase(*it);
            }
            interactions_[key].insert(eam_);
            metExplicit = true;
          }
          
          if (nbiType->isSC()) {
            if (metExplicit) {
              sprintf( painCave.errMsg,
                       "InteractionManager::initialize found more than one "
                       "explicit\n"
                       "\tmetallic interaction for atom types %s - %s\n",
                       atype1->getName().c_str(), atype2->getName().c_str());
              painCave.severity = OPENMD_ERROR;
              painCave.isFatal = 1;
              simError();
            }
            // We found an explicit Sutton-Chen interaction.  
            // override all other metallic entries for this pair of atom types:
            set<NonBondedInteraction*>::iterator it;
            for (it = interactions_[key].begin(); 
                 it != interactions_[key].end(); ++it) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == METALLIC_FAMILY) interactions_[key].erase(*it);
            }
            interactions_[key].insert(sc_);
            metExplicit = true;
          }
          
          if (nbiType->isMAW()) {
            if (vdwExplicit) {
              sprintf( painCave.errMsg,
                       "InteractionManager::initialize found more than one "
                       "explicit\n"
                       "\tvan der Waals interaction for atom types %s - %s\n",
                       atype1->getName().c_str(), atype2->getName().c_str());
              painCave.severity = OPENMD_ERROR;
              painCave.isFatal = 1;
              simError();
            }
            // We found an explicit MAW interaction.  
            // override all other vdw entries for this pair of atom types:
            set<NonBondedInteraction*>::iterator it;
            for (it = interactions_[key].begin(); 
                 it != interactions_[key].end(); ++it) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == VANDERWAALS_FAMILY) interactions_[key].erase(*it);
            }
            interactions_[key].insert(maw_);
            vdwExplicit = true;
          }        
        }
      }
    }
    
    
    // Make sure every pair of atom types in this simulation has a
    // non-bonded interaction.  If not, just inform the user.

    set<AtomType*> simTypes = info_->getSimulatedAtomTypes();
    set<AtomType*>::iterator it, jt;

    for (it = simTypes.begin(); it != simTypes.end(); ++it) {
      atype1 = (*it);
      for (jt = it; jt != simTypes.end(); ++jt) {
        atype2 = (*jt);
        key = make_pair(atype1, atype2);
        
        if (interactions_[key].size() == 0) {
          sprintf( painCave.errMsg,
                   "InteractionManager could not find a matching non-bonded\n"
                   "\tinteraction for atom types %s - %s\n"
                   "\tProceeding without this interaction.\n",
                   atype1->getName().c_str(), atype2->getName().c_str());
          painCave.severity = OPENMD_INFO;
          painCave.isFatal = 0;
          simError();
        }
      }
    }

    initialized_ = true;
  }

  void InteractionManager::setCutoffRadius(RealType rcut) {
    
    electrostatic_->setCutoffRadius(rcut);
    eam_->setCutoffRadius(rcut);
  }

  void InteractionManager::setSwitchingRadius(RealType rswitch) {
    electrostatic_->setSwitchingRadius(rswitch);
  }
  
  void InteractionManager::doPrePair(InteractionData idat){
    
    if (!initialized_) initialize();
	 
    // excluded interaction, so just return
    if (idat.excluded) return;

    set<NonBondedInteraction*>::iterator it;

    for (it = interactions_[ idat.atypes ].begin(); 
         it != interactions_[ idat.atypes ].end(); ++it){
      if ((*it)->getFamily() == METALLIC_FAMILY) {
        dynamic_cast<MetallicInteraction*>(*it)->calcDensity(idat);
      }
    }
    
    return;    
  }
  
  void InteractionManager::doPreForce(SelfData sdat){

    if (!initialized_) initialize();
    
    pair<AtomType*, AtomType*> key = make_pair(sdat.atype, sdat.atype);
    set<NonBondedInteraction*>::iterator it;
    
    for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it){
      if ((*it)->getFamily() == METALLIC_FAMILY) {
        dynamic_cast<MetallicInteraction*>(*it)->calcFunctional(sdat);
      }
    }
    
    return;    
  }

  void InteractionManager::doPair(InteractionData idat){
    
    if (!initialized_) initialize();

    set<NonBondedInteraction*>::iterator it;

    for (it = interactions_[ idat.atypes ].begin(); 
         it != interactions_[ idat.atypes ].end(); ++it) {

      // electrostatics still has to worry about indirect
      // contributions from excluded pairs of atoms:

      if (!idat.excluded || (*it)->getFamily() == ELECTROSTATIC_FAMILY) {
        (*it)->calcForce(idat);
      }
    }
    
    return;    
  }

  void InteractionManager::doSelfCorrection(SelfData sdat){

    if (!initialized_) initialize();
    
    pair<AtomType*, AtomType*> key = make_pair(sdat.atype, sdat.atype);
    set<NonBondedInteraction*>::iterator it;

    for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it){
      if ((*it)->getFamily() == ELECTROSTATIC_FAMILY) {
        dynamic_cast<ElectrostaticInteraction*>(*it)->calcSelfCorrection(sdat);
      }
    }
      
    return;    
  }

  RealType InteractionManager::getSuggestedCutoffRadius(int *atid) {
    if (!initialized_) initialize();
    
    AtomType* atype = typeMap_[*atid];

    pair<AtomType*, AtomType*> key = make_pair(atype, atype);
    set<NonBondedInteraction*>::iterator it;
    RealType cutoff = 0.0;
    
    for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it)
      cutoff = max(cutoff, (*it)->getSuggestedCutoffRadius(key));   
    return cutoff;    
  }

  RealType InteractionManager::getSuggestedCutoffRadius(AtomType* atype) {
    if (!initialized_) initialize();
    
    pair<AtomType*, AtomType*> key = make_pair(atype, atype);
    set<NonBondedInteraction*>::iterator it;
    RealType cutoff = 0.0;
    
    for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it)
      cutoff = max(cutoff, (*it)->getSuggestedCutoffRadius(key));   
    return cutoff;    
  }
} //end namespace OpenMD
