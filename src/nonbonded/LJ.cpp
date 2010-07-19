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

#include <stdio.h>
#include <string.h>

#include <cmath>
#include "nonbonded/LJ.hpp"
#include "utils/simError.h"


namespace OpenMD {

  bool LJ::initialized_ = false;
  bool LJ::shiftedPot_ = false;
  bool LJ::shiftedFrc_ = false;
  ForceField* LJ::forceField_ = NULL;
  std::map<int, AtomType*> LJ::LJMap;
  std::map<std::pair<AtomType*, AtomType*>, LJInteractionData> LJ::MixingMap;
  
  LJ* LJ::_instance = NULL;

  LJ* LJ::Instance() {
    if (!_instance) {
      _instance = new LJ();
    }
    return _instance;
  }

  LJParam LJ::getLJParam(AtomType* atomType) {
    
    // Do sanity checking on the AtomType we were passed before
    // building any data structures:
    if (!atomType->isLennardJones()) {
      sprintf( painCave.errMsg,
               "LJ::getLJParam was passed an atomType (%s) that does not\n"
               "\tappear to be a Lennard-Jones atom.\n",
               atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();
    }
    
    GenericData* data = atomType->getPropertyByName("LennardJones");
    if (data == NULL) {
      sprintf( painCave.errMsg, "LJ::getLJParam could not find Lennard-Jones\n"
               "\tparameters for atomType %s.\n", atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError(); 
    }
    
    LJParamGenericData* ljData = dynamic_cast<LJParamGenericData*>(data);
    if (ljData == NULL) {
      sprintf( painCave.errMsg,
               "LJ::getLJParam could not convert GenericData to LJParam for\n"
               "\tatom type %s\n", atomType->getName().c_str());
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    
    return ljData->getData();
  }

  RealType LJ::getSigma(AtomType* atomType) {    
    LJParam ljParam = getLJParam(atomType);
    return ljParam.sigma;
  }

  RealType LJ::getSigma(int atid) { 
    std::map<int, AtomType*> :: const_iterator it;
    it = LJMap.find(atid);
    if (it == LJMap.end()) {
      sprintf( painCave.errMsg,
               "LJ::getSigma could not find atid %d in LJMap\n",
               (atid));
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    AtomType* atype = it->second;

    return getSigma(atype);    
  }

  RealType LJ::getSigma(AtomType* atomType1, AtomType* atomType2) {    
    RealType sigma1 = getSigma(atomType1);
    RealType sigma2 = getSigma(atomType2);
    
    ForceFieldOptions& fopts = forceField_->getForceFieldOptions();
    std::string DistanceMix = fopts.getDistanceMixingRule();
    toUpper(DistanceMix);

    if (DistanceMix == "GEOMETRIC") 
      return sqrt(sigma1 * sigma2);
    else 
      return 0.5 * (sigma1 + sigma2);
  }

  RealType LJ::getEpsilon(AtomType* atomType) {    
    LJParam ljParam = getLJParam(atomType);
    return ljParam.epsilon;
  }

  RealType LJ::getEpsilon(int atid) {    
    std::map<int, AtomType*> :: const_iterator it;
    it = LJMap.find(atid);
    if (it == LJMap.end()) {
      sprintf( painCave.errMsg,
               "LJ::getEpsilon could not find atid %d in LJMap\n",
               (atid));
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal = 1;
      simError();          
    }
    AtomType* atype = it->second;

    return getEpsilon(atype);    
  }

  RealType LJ::getEpsilon(AtomType* atomType1, AtomType* atomType2) {    
    RealType epsilon1 = getEpsilon(atomType1);
    RealType epsilon2 = getEpsilon(atomType2);
    return sqrt(epsilon1 * epsilon2);
  }

  void LJ::initialize() {    
    ForceField::AtomTypeContainer atomTypes = forceField_->getAtomTypes();
    ForceField::AtomTypeContainer::MapTypeIterator i;
    AtomType* at;

    for (at = atomTypes.beginType(i); at != NULL; 
         at = atomTypes.nextType(i)) {
      
      if (at->isLennardJones())
        addType(at);
    }

    ForceField::NonBondedInteractionTypeContainer nbiTypes = forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    NonBondedInteractionType* nbt;

    for (nbt = nbiTypes.beginType(j); nbt != NULL; 
         nbt = nbiTypes.nextType(j)) {
      
      if (nbt->isLennardJones()) {
        
        std::pair<AtomType*, AtomType*> atypes = nbt->getAtomTypes();
        
        GenericData* data = nbt->getPropertyByName("LennardJones");
        if (data == NULL) {
          sprintf( painCave.errMsg, "LJ::rebuildMixingMap could not find\n"
               "\tLennard-Jones parameters for %s - %s interaction.\n", 
                   atypes.first->getName().c_str(),
                   atypes.second->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError(); 
        }
    
        LJParamGenericData* ljData = dynamic_cast<LJParamGenericData*>(data);
        if (ljData == NULL) {
          sprintf( painCave.errMsg,
                   "LJ::rebuildMixingMap could not convert GenericData to\n"
                   "\tLJParam for %s - %s interaction.\n", 
                   atypes.first->getName().c_str(),
                   atypes.second->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
        
        LJParam ljParam = ljData->getData();

        RealType sigma = ljParam.sigma;
        RealType epsilon = ljParam.epsilon;

        addExplicitInteraction(atypes.first, atypes.second, sigma, epsilon);
      }
    }  
    initialized_ = true;
  }
      


  void LJ::addType(AtomType* atomType){
    RealType sigma1 = getSigma(atomType);
    RealType epsilon1 = getEpsilon(atomType);

    // add it to the map:
    AtomTypeProperties atp = atomType->getATP();    
    std::pair<std::map<int,AtomType*>::iterator,bool> ret;
    ret = LJMap.insert( std::pair<int, AtomType*>(atp.ident, atomType) );
    if (ret.second == false) {
      sprintf( painCave.errMsg,
               "LJ already had a previous entry with ident %d\n",
               atp.ident);
      painCave.severity = OPENMD_INFO;
      painCave.isFatal = 0;
      simError();         
    }
    
    // Now, iterate over all known types and add to the mixing map:
    
    std::map<int, AtomType*>::iterator it;
    for( it = LJMap.begin(); it != LJMap.end(); ++it) {
      
      AtomType* atype2 = (*it).second;

      LJInteractionData mixer;
      mixer.sigma = getSigma(atomType, atype2);
      mixer.epsilon = getEpsilon(atomType, atype2);
      mixer.sigmai = 1.0 / mixer.sigma;
      mixer.explicitlySet = false;

      std::pair<AtomType*, AtomType*> key1, key2;
      key1 = std::make_pair(atomType, atype2);
      key2 = std::make_pair(atype2, atomType);
      
      MixingMap[key1] = mixer;
      if (key2 != key1) {
        MixingMap[key2] = mixer;
      }
    }      
  }
  
  void LJ::addExplicitInteraction(AtomType* atype1, AtomType* atype2, RealType sigma, RealType epsilon){
    
    // in case these weren't already in the map
    addType(atype1);
    addType(atype2);

    LJInteractionData mixer;
    mixer.sigma = sigma;
    mixer.epsilon = epsilon;
    mixer.sigmai = 1.0 / mixer.sigma;
    mixer.explicitlySet = true;

    std::pair<AtomType*, AtomType*> key1, key2;
    key1 = std::make_pair(atype1, atype2);
    key2 = std::make_pair(atype2, atype1);
    
    MixingMap[key1] = mixer;
    if (key2 != key1) {
      MixingMap[key2] = mixer;
    }    
  }
 
  void LJ::calcForce(AtomType* at1, AtomType* at2, Vector3d d, 
                     RealType rij, RealType r2, RealType rcut, RealType sw, 
                     RealType vdwMult, RealType &vpair, RealType &pot, 
                     Vector3d &f1) {

    if (!initialized_) initialize();
    
    RealType ros;
    RealType rcos;
    RealType myPot = 0.0;
    RealType myPotC = 0.0;
    RealType myDeriv = 0.0;
    RealType myDerivC = 0.0;

    std::pair<AtomType*, AtomType*> key = std::make_pair(at1, at2);
    LJInteractionData mixer = MixingMap[key];

    RealType sigmai = mixer.sigmai;
    RealType epsilon = mixer.epsilon;
    
    ros = rij * sigmai;

    getLJfunc(ros, myPot, myDeriv);

    if (shiftedPot_) {
      rcos = rcut * sigmai;
      getLJfunc(rcos, myPotC, myDerivC);
      myDerivC = 0.0;
    } else if (LJ::shiftedFrc_) {
      rcos = rcut * sigmai;
      getLJfunc(rcos, myPotC, myDerivC);
      myPotC = myPotC + myDerivC * (rij - rcut) * sigmai;
    } else {
      myPotC = 0.0;
      myDerivC = 0.0;
    }

    RealType pot_temp = vdwMult * epsilon * (myPot - myPotC);
    vpair += pot_temp;

    RealType dudr = sw * vdwMult * epsilon * (myDeriv - myDerivC)*sigmai;
    
    pot += sw * pot_temp;
    f1 = d * dudr / rij;

    return;
  }

  void LJ::do_lj_pair(int *atid1, int *atid2, RealType *d, RealType *rij, 
                      RealType *r2, RealType *rcut, RealType *sw, 
                      RealType *vdwMult,
                      RealType *vpair, RealType *pot, RealType *f1) {

    if (!initialized_) initialize();
    
    AtomType* atype1 = LJMap[*atid1];
    AtomType* atype2 = LJMap[*atid2];
    
    Vector3d disp(d[0], d[1], d[2]);
    Vector3d frc(f1[0], f1[1], f1[2]);
    
    calcForce(atype1, atype2, disp, *rij, *r2, *rcut, *sw, *vdwMult, *vpair, 
              *pot, frc);
      
    f1[0] = frc.x();
    f1[1] = frc.y();
    f1[2] = frc.z();
    
    return;    
  }
  
  void LJ::getLJfunc(RealType r, RealType &pot, RealType &deriv) {

    RealType ri = 1.0 / r;
    RealType ri2 = ri * ri;
    RealType ri6 = ri2 * ri2 * ri2;
    RealType ri7 = ri6 * ri;
    RealType ri12 = ri6 * ri6;
    RealType ri13 = ri12 * ri;
    
    pot = 4.0 * (ri12 - ri6);
    deriv = 24.0 * (ri7 - 2.0 * ri13);

    return;
  }
  

  void LJ::setLJDefaultCutoff(RealType *thisRcut, int *sP, int *sF) {
    shiftedPot_ = (bool)(*sP);
    shiftedFrc_ = (bool)(*sF);
  }
}

extern "C" {
  
#define fortranGetSigma FC_FUNC(getsigma, GETSIGMA)
#define fortranGetEpsilon FC_FUNC(getepsilon, GETEPSILON)
#define fortranSetLJCutoff FC_FUNC(setljdefaultcutoff, SETLJDEFAULTCUTOFF)
#define fortranDoLJPair FC_FUNC(do_lj_pair, DO_LJ_PAIR)
  
  RealType fortranGetSigma(int* atid) {
    return OpenMD::LJ::Instance()->getSigma(*atid);
  }
  RealType fortranGetEpsilon(int* atid) {  
    return OpenMD::LJ::Instance()->getEpsilon(*atid);
  }
  void fortranSetLJCutoff(RealType *rcut, int *shiftedPot, int *shiftedFrc) {
    return OpenMD::LJ::Instance()->setLJDefaultCutoff(rcut, shiftedPot, 
                                                      shiftedFrc);
  }
  void fortranDoLJPair(int *atid1, int *atid2, RealType *d, RealType *rij, 
                       RealType *r2, RealType *rcut, RealType *sw, 
                       RealType *vdwMult, RealType* vpair, RealType* pot, 
                       RealType *f1){
    
    return OpenMD::LJ::Instance()->do_lj_pair(atid1, atid2, d, rij, r2, rcut,
                                              sw, vdwMult, vpair, pot, f1);
  }
}
