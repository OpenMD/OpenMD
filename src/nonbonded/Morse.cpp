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
#include "nonbonded/Morse.hpp"
#include "utils/simError.h"
#include "types/NonBondedInteractionType.hpp"

using namespace std;

namespace OpenMD {

  bool Morse::initialized_ = false;
  bool Morse::shiftedPot_ = false;
  bool Morse::shiftedFrc_ = false;
  ForceField* Morse::forceField_ = NULL;
  map<int, AtomType*> Morse::MorseMap;
  map<pair<AtomType*, AtomType*>, MorseInteractionData> Morse::MixingMap;
  map<string, MorseInteractionType> Morse::stringToEnumMap_;
  
  Morse* Morse::_instance = NULL;

  Morse* Morse::Instance() {
    if (!_instance) {
      _instance = new Morse();
    }
    return _instance;
  }

  void Morse::initialize() {    

    stringToEnumMap_["shiftedMorse"] = shiftedMorse;
    stringToEnumMap_["repulsiveMorse"] = repulsiveMorse;

    ForceField::NonBondedInteractionTypeContainer* nbiTypes = forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    NonBondedInteractionType* nbt;

    for (nbt = nbiTypes->beginType(j); nbt != NULL; 
         nbt = nbiTypes->nextType(j)) {
      
      if (nbt->isMorse()) {
        
        pair<AtomType*, AtomType*> atypes = nbt->getAtomTypes();
        
        GenericData* data = nbt->getPropertyByName("Morse");
        if (data == NULL) {
          sprintf( painCave.errMsg, "Morse::initialize could not find\n"
                   "\tMorse parameters for %s - %s interaction.\n", 
                   atypes.first->getName().c_str(),
                   atypes.second->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError(); 
        }
        
        MorseData* morseData = dynamic_cast<MorseData*>(data);
        if (morseData == NULL) {
          sprintf( painCave.errMsg,
                   "Morse::initialize could not convert GenericData to\n"
                   "\tMorseData for %s - %s interaction.\n", 
                   atypes.first->getName().c_str(),
                   atypes.second->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
        
        MorseParam morseParam = morseData->getData();

        RealType De = morseParam.De;
        RealType Re = morseParam.Re;
        RealType beta = morseParam.beta;
        string interactionType = morseParam.interactionType;

        toUpper(interactionType);
        map<string, MorseInteractionType>::iterator i;
        i = stringToEnumMap_.find(interactionType);
        if (i != stringToEnumMap_.end()) { 
          addExplicitInteraction(atypes.first, atypes.second, 
                                 De, Re, beta, i->second );
        } else {
          sprintf( painCave.errMsg,
                   "Morse::initialize found unknown Morse interaction type\n"
                   "\t(%s) for %s - %s interaction.\n", 
                   morseParam.interactionType.c_str(),
                   atypes.first->getName().c_str(),
                   atypes.second->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
      }
    }  
    initialized_ = true;
  }
      
  void Morse::addExplicitInteraction(AtomType* atype1, AtomType* atype2, 
                                     RealType De, RealType Re, RealType beta, 
                                     MorseInteractionType mit) {

    AtomTypeProperties atp = atype1->getATP();    
    MorseMap.insert( pair<int, AtomType*>(atp.ident, atype1) );

    atp = atype2->getATP();    
    MorseMap.insert( pair<int, AtomType*>(atp.ident, atype2) );

    MorseInteractionData mixer;
    mixer.De = De;
    mixer.Re = Re;
    mixer.beta = beta;
    mixer.interactionType = mit;

    pair<AtomType*, AtomType*> key1, key2;
    key1 = make_pair(atype1, atype2);
    key2 = make_pair(atype2, atype1);
    
    MixingMap[key1] = mixer;
    if (key2 != key1) {
      MixingMap[key2] = mixer;
    }    
  }
  
  void Morse::calcForce(AtomType* at1, AtomType* at2, Vector3d d, 
                        RealType r, RealType r2, RealType rcut, RealType sw, 
                        RealType vdwMult, RealType &vpair, RealType &pot, 
                        Vector3d &f1) {

    if (!initialized_) initialize();
    
    RealType myPot = 0.0;
    RealType myPotC = 0.0;
    RealType myDeriv = 0.0;
    RealType myDerivC = 0.0;

    pair<AtomType*, AtomType*> key = make_pair(at1, at2);
    MorseInteractionData mixer = MixingMap[key];

    RealType De = mixer.De;
    RealType Re = mixer.Re;
    RealType beta = mixer.beta;
    MorseInteractionType interactionType = mixer.interactionType;
   
    // V(r) = D_e exp(-a(r-re)(exp(-a(r-re))-2)
    
    RealType expt     = -beta*(r - Re);
    RealType expfnc   = exp(expt);
    RealType expfnc2  = expfnc*expfnc;

    RealType exptC = 0.0;
    RealType expfncC = 0.0;
    RealType expfnc2C = 0.0;

    if (Morse::shiftedPot_ || Morse::shiftedFrc_) {
      exptC     = -beta*(rcut - Re);
      expfncC   = exp(exptC);
      expfnc2C  = expfncC*expfncC;
    }


    switch(interactionType) {
    case shiftedMorse : {

      myPot  = De * (expfnc2  - 2.0 * expfnc);
      myDeriv   = 2.0 * De * beta * (expfnc - expfnc2);

      if (Morse::shiftedPot_) {
        myPotC = De * (expfnc2C - 2.0 * expfncC);
        myDerivC = 0.0;
      } else if (Morse::shiftedFrc_) {
        myPotC = De * (expfnc2C - 2.0 * expfncC);
        myDerivC  = 2.0 * De * beta * (expfnc2C - expfnc2C);
        myPotC += myDerivC * (r - rcut);
      } else {
        myPotC = 0.0;
        myDerivC = 0.0;
      }

      break;
    }
    case repulsiveMorse : {

      myPot  = De * expfnc2;
      myDeriv  = -2.0 * De * beta * expfnc2;

      if (Morse::shiftedPot_) {
        myPotC = De * expfnc2C;
        myDerivC = 0.0;
      } else if (Morse::shiftedFrc_) {
        myPotC = De * expfnc2C;
        myDerivC = -2.0 * De * beta * expfnc2C;
        myPotC += myDerivC * (r - rcut);
      } else {
        myPotC = 0.0;
        myDerivC = 0.0;
      }

      break;
    }
    }

    RealType pot_temp = vdwMult * (myPot - myPotC);
    vpair += pot_temp;

    RealType dudr = sw * vdwMult * (myDeriv - myDerivC);
    
    pot += sw * pot_temp;
    f1 = d * dudr / r;

    return;

  }

  void Morse::do_morse_pair(int *atid1, int *atid2, RealType *d, 
                            RealType *rij, RealType *r2, RealType *rcut, 
                            RealType *sw, RealType *vdwMult,
                            RealType *vpair, RealType *pot, RealType *f1) {

    if (!initialized_) initialize();
    
    AtomType* atype1 = MorseMap[*atid1];
    AtomType* atype2 = MorseMap[*atid2];
    
    Vector3d disp(d[0], d[1], d[2]);
    Vector3d frc(f1[0], f1[1], f1[2]);
    
    calcForce(atype1, atype2, disp, *rij, *r2, *rcut, *sw, *vdwMult, *vpair, 
              *pot, frc);
      
    f1[0] = frc.x();
    f1[1] = frc.y();
    f1[2] = frc.z();

    return;    
  }
    
  void Morse::setMorseDefaultCutoff(RealType *thisRcut, int *sP, int *sF) {
    shiftedPot_ = (bool)(*sP);
    shiftedFrc_ = (bool)(*sF);
  }
}

extern "C" {
  
#define fortranSetMorseCutoff FC_FUNC(setmorsedefaultcutoff, SETMORSEDEFAULTCUTOFF)
#define fortranDoMorsePair FC_FUNC(do_morse_pair, DO_MORSE_PAIR)
  
  void fortranSetMorseCutoff(RealType *rcut, int *shiftedPot, int *shiftedFrc) {
    return OpenMD::Morse::Instance()->setMorseDefaultCutoff(rcut, shiftedPot, 
                                                            shiftedFrc);
  }
  void fortranDoMorsePair(int *atid1, int *atid2, RealType *d, RealType *rij, 
                          RealType *r2, RealType *rcut, RealType *sw, 
                          RealType *vdwMult, RealType* vpair, RealType* pot, 
                          RealType *f1){
    
    return OpenMD::Morse::Instance()->do_morse_pair(atid1, atid2, d, rij, r2, 
                                                    rcut, sw, vdwMult, vpair, 
                                                    pot, f1);
  }
}
