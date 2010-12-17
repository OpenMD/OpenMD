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
 
  bool InteractionManager::initialized_ = false;
  ForceField* InteractionManager::forceField_ = NULL;  
  InteractionManager* InteractionManager::_instance = NULL;
  map<int, AtomType*> InteractionManager::typeMap_;
  map<pair<AtomType*, AtomType*>, set<NonBondedInteraction*> > InteractionManager::interactions_;

  LJ* InteractionManager::lj_ = new LJ();
  GB* InteractionManager::gb_ = new GB();
  Sticky* InteractionManager::sticky_ = new Sticky();
  Morse* InteractionManager::morse_ = new Morse();
  EAM* InteractionManager::eam_ = new EAM();
  SC* InteractionManager::sc_ = new SC();
  Electrostatic* InteractionManager::electrostatic_ = new Electrostatic();
 
  InteractionManager* InteractionManager::Instance() {
    if (!_instance) {
      _instance = new InteractionManager();
    }
    return _instance;
  }

  void InteractionManager::initialize() {
    
    lj_->setForceField(forceField_);
    gb_->setForceField(forceField_);
    sticky_->setForceField(forceField_);
    eam_->setForceField(forceField_);
    sc_->setForceField(forceField_);
    morse_->setForceField(forceField_);
    electrostatic_->setForceField(forceField_);

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

        if (nbiType->isLennardJones()) {
          // We found an explicit Lennard-Jones interaction.  
          // override all other vdw entries for this pair of atom types:
          set<NonBondedInteraction*>::iterator it;
          for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it) {
            InteractionFamily ifam = (*it)->getFamily();
            if (ifam == VANDERWAALS_FAMILY) interactions_[key].erase(*it);
          }
          interactions_[key].insert(lj_);
          vdwExplicit = true;
        }

        if (nbiType->isMorse()) {
          if (vdwExplicit) {
            sprintf( painCave.errMsg,
                     "InteractionManager::initialize found more than one explicit\n"
                     "\tvan der Waals interaction for atom types %s - %s\n",
                     atype1->getName().c_str(), atype2->getName().c_str());
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();
          }
          // We found an explicit Morse interaction.  
          // override all other vdw entries for this pair of atom types:
          set<NonBondedInteraction*>::iterator it;
          for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it) {
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
          for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it) {
            InteractionFamily ifam = (*it)->getFamily();
            if (ifam == METALLIC_FAMILY) interactions_[key].erase(*it);
          }
          interactions_[key].insert(eam_);
          metExplicit = true;
        }

        if (nbiType->isSC()) {
          if (metExplicit) {
            sprintf( painCave.errMsg,
                     "InteractionManager::initialize found more than one explicit\n"
                     "\tmetallic interaction for atom types %s - %s\n",
                     atype1->getName().c_str(), atype2->getName().c_str());
            painCave.severity = OPENMD_ERROR;
            painCave.isFatal = 1;
            simError();
          }
          // We found an explicit Sutton-Chen interaction.  
          // override all other metallic entries for this pair of atom types:
          set<NonBondedInteraction*>::iterator it;
          for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it) {
            InteractionFamily ifam = (*it)->getFamily();
            if (ifam == METALLIC_FAMILY) interactions_[key].erase(*it);
          }
          interactions_[key].insert(sc_);
          metExplicit = true;
        }
      }
    }
    
    // make sure every pair of atom types has a non-bonded interaction:
    for (atype1 = atomTypes->beginType(i1); atype1 != NULL; 
         atype1 = atomTypes->nextType(i1)) {
      for (atype2 = atomTypes->beginType(i2); atype2 != NULL; 
           atype2 = atomTypes->nextType(i2)) {
        key = make_pair(atype1, atype2);
        
        if (interactions_[key].size() == 0) {
          sprintf( painCave.errMsg,
                   "InteractionManager unable to find an appropriate non-bonded\n"
                   "\tinteraction for atom types %s - %s\n",
                   atype1->getName().c_str(), atype2->getName().c_str());
          painCave.severity = OPENMD_INFO;
          painCave.isFatal = 1;
          simError();
        }
      }
    }
  } 
  
  void InteractionManager::doPrePair(int *atid1, int *atid2, RealType *rij, RealType *rho_i_at_j, RealType *rho_j_at_i){
    
    if (!initialized_) initialize();
	  
    DensityData ddat;

    ddat.atype1 = typeMap_[*atid1];
    ddat.atype2 = typeMap_[*atid2];
    ddat.rij = *rij;
    ddat.rho_i_at_j = *rho_i_at_j;
    ddat.rho_j_at_i = *rho_j_at_i;

    pair<AtomType*, AtomType*> key = make_pair(ddat.atype1, ddat.atype2);
    set<NonBondedInteraction*>::iterator it;

    for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it){
      if ((*it)->getFamily() == METALLIC_FAMILY) {
        dynamic_cast<MetallicInteraction*>(*it)->calcDensity(ddat);
      }
    }
    
    return;    
  }
  
  void InteractionManager::doPreForce(int *atid, RealType *rho, RealType *frho, RealType *dfrhodrho){

    if (!initialized_) initialize();
	  
    FunctionalData fdat;

    fdat.atype = typeMap_[*atid];
    fdat.rho = *rho;
    fdat.frho = *frho;
    fdat.dfrhodrho = *dfrhodrho;

    pair<AtomType*, AtomType*> key = make_pair(fdat.atype, fdat.atype);
    set<NonBondedInteraction*>::iterator it;
    
    for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it){
      if ((*it)->getFamily() == METALLIC_FAMILY) {
        dynamic_cast<MetallicInteraction*>(*it)->calcFunctional(fdat);
      }
    }
    
    return;    
  }

  void InteractionManager::doPair(int *atid1, int *atid2, RealType *d, RealType *r, RealType *r2, RealType *rcut, RealType *sw, RealType *vdwMult,RealType *electroMult, RealType *pot, RealType *vpair, RealType *f1, RealType *eFrame1, RealType *eFrame2, RealType *A1, RealType *A2, RealType *t1, RealType *t2, RealType *rho1, RealType *rho2, RealType *dfrho1, RealType *dfrho2, RealType *fshift1, RealType *fshift2){
    
    if (!initialized_) initialize();
    
    InteractionData idat;
    
    idat.atype1 = typeMap_[*atid1];
    idat.atype2 = typeMap_[*atid2];
    idat.d = Vector3d(d);
    idat.rij = *r;
    idat.r2 = *r2;
    idat.rcut = *rcut;
    idat.sw = *sw;
    idat.vdwMult = *vdwMult;
    idat.electroMult = *electroMult;
    idat.pot = *pot;
    idat.vpair = *vpair;
    idat.f1 = Vector3d(f1);
    idat.eFrame1 = Mat3x3d(eFrame1);
    idat.eFrame2 = Mat3x3d(eFrame2);
    idat.A1 = RotMat3x3d(A1);
    idat.A2 = RotMat3x3d(A2);
    idat.t1 = Vector3d(t1);
    idat.t2 = Vector3d(t2);
    idat.rho1 = *rho1;
    idat.rho2 = *rho2;
    idat.dfrho1 = *dfrho1;
    idat.dfrho2 = *dfrho2;
    idat.fshift1 = *fshift1;
    idat.fshift2 = *fshift2;

    pair<AtomType*, AtomType*> key = make_pair(idat.atype1, idat.atype2);
    set<NonBondedInteraction*>::iterator it;

    for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it)
      (*it)->calcForce(idat);
    
    f1[0] = idat.f1.x();
    f1[1] = idat.f1.y();
    f1[2] = idat.f1.z();
    
    t1[0] = idat.t1.x();
    t1[1] = idat.t1.y();
    t1[2] = idat.t1.z();
    
    t2[0] = idat.t2.x();
    t2[1] = idat.t2.y();
    t2[2] = idat.t2.z();

    return;    
  }

  void InteractionManager::doSkipCorrection(int *atid1, int *atid2, RealType *d, RealType *r, RealType *skippedCharge1, RealType *skippedCharge2, RealType *sw, RealType *electroMult, RealType *pot, RealType *vpair, RealType *f1, RealType *eFrame1, RealType *eFrame2, RealType *t1, RealType *t2){

    if (!initialized_) initialize();
    
    SkipCorrectionData skdat;
    
    skdat.atype1 = typeMap_[*atid1];
    skdat.atype2 = typeMap_[*atid2];
    skdat.d = Vector3d(d);
    skdat.rij = *r;
    skdat.skippedCharge1 = *skippedCharge1;
    skdat.skippedCharge2 = *skippedCharge2;
    skdat.sw = *sw;
    skdat.electroMult = *electroMult;
    skdat.pot = *pot;
    skdat.vpair = *vpair;
    skdat.f1 = Vector3d(f1);
    skdat.eFrame1 = Mat3x3d(eFrame1);
    skdat.eFrame2 = Mat3x3d(eFrame2);
    skdat.t1 = Vector3d(t1);
    skdat.t2 = Vector3d(t2);

    pair<AtomType*, AtomType*> key = make_pair(skdat.atype1, skdat.atype2);
    set<NonBondedInteraction*>::iterator it;

    for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it){
      if ((*it)->getFamily() == ELECTROSTATIC_FAMILY) {
        dynamic_cast<ElectrostaticInteraction*>(*it)->calcSkipCorrection(skdat);
      }
    }
    
    f1[0] = skdat.f1.x();
    f1[1] = skdat.f1.y();
    f1[2] = skdat.f1.z();
    
    t1[0] = skdat.t1.x();
    t1[1] = skdat.t1.y();
    t1[2] = skdat.t1.z();
    
    t2[0] = skdat.t2.x();
    t2[1] = skdat.t2.y();
    t2[2] = skdat.t2.z();

    return;    
  }

  void InteractionManager::doSelfCorrection(int *atid, RealType *eFrame, RealType *skippedCharge, RealType *pot, RealType *t){

    if (!initialized_) initialize();
    
    SelfCorrectionData scdat;
    
    scdat.atype = typeMap_[*atid];
    scdat.eFrame = Mat3x3d(eFrame);
    scdat.skippedCharge = *skippedCharge;
    scdat.pot = *pot;
    scdat.t = Vector3d(t);

    pair<AtomType*, AtomType*> key = make_pair(scdat.atype, scdat.atype);
    set<NonBondedInteraction*>::iterator it;

    for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it){
      if ((*it)->getFamily() == ELECTROSTATIC_FAMILY) {
        dynamic_cast<ElectrostaticInteraction*>(*it)->calcSelfCorrection(scdat);
      }
    }
        
    t[0] = scdat.t.x();
    t[1] = scdat.t.y();
    t[2] = scdat.t.z();

    return;    
  }


  RealType InteractionManager::getSuggestedCutoffRadius(int *atid) {
    if (!initialized_) initialize();
    
    AtomType* atype = typeMap_[*atid];

    pair<AtomType*, AtomType*> key = make_pair(atype, atype);
    set<NonBondedInteraction*>::iterator it;
    RealType cutoff = 0.0;
    
    for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it)
      cutoff = max(cutoff, (*it)->getSuggestedCutoffRadius(atype, atype));   
    return cutoff;    
  }

  RealType InteractionManager::getSuggestedCutoffRadius(AtomType* atype) {
    if (!initialized_) initialize();
    
    pair<AtomType*, AtomType*> key = make_pair(atype, atype);
    set<NonBondedInteraction*>::iterator it;
    RealType cutoff = 0.0;
    
    for (it = interactions_[key].begin(); it != interactions_[key].end(); ++it)
      cutoff = max(cutoff, (*it)->getSuggestedCutoffRadius(atype, atype));   
    return cutoff;    
  }

} //end namespace OpenMD

extern "C" {
  
#define fortranDoPrePair FC_FUNC(do_prepair, DO_PREPAIR)
#define fortranDoPreForce FC_FUNC(do_preforce, DO_PREFORCE)
#define fortranDoPair FC_FUNC(do_pair, DO_PAIR)
#define fortranDoSkipCorrection FC_FUNC(do_skip_correction, DO_SKIP_CORRECTION)
#define fortranDoSelfCorrection FC_FUNC(do_self_correction, DO_SELF_CORRECTION)
#define fortranGetCutoff FC_FUNC(get_cutoff, GET_CUTOFF)

  void fortranDoPrePair(int *atid1, int *atid2, RealType *rij, 
                        RealType *rho_i_at_j, RealType *rho_j_at_i) {
	    
    return OpenMD::InteractionManager::Instance()->doPrePair(atid1, atid2, rij,
                                                             rho_i_at_j,  
                                                             rho_j_at_i);
  }
  void fortranDoPreForce(int *atid, RealType *rho, RealType *frho,
                         RealType *dfrhodrho) {  
    
    return OpenMD::InteractionManager::Instance()->doPreForce(atid, rho, frho,
                                                              dfrhodrho);    
  }
  
  void fortranDoPair(int *atid1, int *atid2, RealType *d, RealType *r, 
                     RealType *r2, RealType *rcut, RealType *sw, 
                     RealType *vdwMult, RealType *electroMult, RealType *pot, 
                     RealType *vpair, RealType *f1, RealType *eFrame1, 
                     RealType *eFrame2, RealType *A1, RealType *A2, 
                     RealType *t1, RealType *t2, RealType *rho1, RealType *rho2,
                     RealType *dfrho1, RealType *dfrho2, RealType *fshift1, 
                     RealType *fshift2){
    
    return OpenMD::InteractionManager::Instance()->doPair(atid1, atid2, d, r, 
                                                          r2, rcut, sw, 
                                                          vdwMult, electroMult,
                                                          pot, vpair, f1, 
                                                          eFrame1, eFrame2, 
                                                          A1, A2, t1, t2, rho1,
                                                          rho2, dfrho1, dfrho2,
                                                          fshift1, fshift2);
  }
  
  void fortranDoSkipCorrection(int *atid1, int *atid2, RealType *d, 
                               RealType *r, RealType *skippedCharge1, 
                               RealType *skippedCharge2, RealType *sw, 
                               RealType *electroMult, RealType *pot, 
                               RealType *vpair, RealType *f1,
                               RealType *eFrame1, RealType *eFrame2,
                               RealType *t1, RealType *t2){
    
    return OpenMD::InteractionManager::Instance()->doSkipCorrection(atid1, 
                                                                    atid2, d, 
                                                                    r, 
                                                                    skippedCharge1,
                                                                    skippedCharge2,
                                                                    sw, electroMult, pot,
                                                                    vpair, f1, eFrame1, 
                                                                    eFrame2, t1, t2);
  }
  
  void fortranDoSelfCorrection(int *atid, RealType *eFrame, RealType *skippedCharge, 
                               RealType *pot, RealType *t) {
    
    return OpenMD::InteractionManager::Instance()->doSelfCorrection(atid, 
                                                                    eFrame, 
                                                                    skippedCharge,
                                                                    pot, t);
  }
  RealType fortranGetCutoff(int *atid) {    
    return OpenMD::InteractionManager::Instance()->getSuggestedCutoffRadius(atid);
  }
}
