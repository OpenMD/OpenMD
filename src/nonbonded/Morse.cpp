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
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */

#include <stdio.h>
#include <string.h>

#include <cmath>
#include "nonbonded/Morse.hpp"
#include "utils/simError.h"
#include "types/MorseInteractionType.hpp"

using namespace std;

namespace OpenMD {

  Morse::Morse() : name_("Morse"), initialized_(false), forceField_(NULL) {}
  
  void Morse::initialize() {    

    Mtypes.clear();
    Mtids.clear();
    MixingMap.clear();
    Mtids.resize( forceField_->getNAtomType(), -1);

    ForceField::NonBondedInteractionTypeContainer* nbiTypes = forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    ForceField::NonBondedInteractionTypeContainer::KeyType keys;
    NonBondedInteractionType* nbt;
    int mtid1, mtid2;

    for (nbt = nbiTypes->beginType(j); nbt != NULL; 
         nbt = nbiTypes->nextType(j)) {
      
      if (nbt->isMorse()) {
        keys = nbiTypes->getKeys(j);
        AtomType* at1 = forceField_->getAtomType(keys[0]);
        AtomType* at2 = forceField_->getAtomType(keys[1]);

        int atid1 = at1->getIdent();
        if (Mtids[atid1] == -1) {
          mtid1 = Mtypes.size();
          Mtypes.insert(atid1);
          Mtids[atid1] = mtid1;         
        }
        int atid2 = at2->getIdent();
        if (Mtids[atid2] == -1) {
          mtid2 = Mtypes.size();
          Mtypes.insert(atid2);
          Mtids[atid2] = mtid2;
        }
    
        MorseInteractionType* mit = dynamic_cast<MorseInteractionType*>(nbt);

        if (mit == NULL) {
          sprintf( painCave.errMsg,
                   "Morse::initialize could not convert NonBondedInteractionType\n"
                   "\tto MorseInteractionType for %s - %s interaction.\n", 
                   at1->getName().c_str(),
                   at2->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
        
        RealType De = mit->getD();
        RealType Re = mit->getR();
        RealType beta = mit->getBeta();
   
        MorseType variant = mit->getInteractionType();
        addExplicitInteraction(at1, at2, De, Re, beta, variant );
      }
    }  
    initialized_ = true;
  }
      
  void Morse::addExplicitInteraction(AtomType* atype1, AtomType* atype2, 
                                     RealType De, RealType Re, RealType beta, 
                                     MorseType mt) {

    MorseInteractionData mixer;
    mixer.De = De;
    mixer.Re = Re;
    mixer.beta = beta;
    mixer.variant = mt;

    int mtid1 = Mtids[atype1->getIdent()];
    int mtid2 = Mtids[atype2->getIdent()];
    int nM = Mtypes.size();

    MixingMap.resize(nM);
    MixingMap[mtid1].resize(nM);
    
    MixingMap[mtid1][mtid2] = mixer;
    if (mtid2 != mtid1) {
      MixingMap[mtid2].resize(nM);
      MixingMap[mtid2][mtid1] = mixer;
    }    
  }
  
  void Morse::calcForce(InteractionData &idat) {

    if (!initialized_) initialize();
    
    MorseInteractionData &mixer = MixingMap[Mtids[idat.atid1]][Mtids[idat.atid2]];

    RealType myPot = 0.0;
    RealType myPotC = 0.0;
    RealType myDeriv = 0.0;
    RealType myDerivC = 0.0;
    
    RealType De = mixer.De;
    RealType Re = mixer.Re;
    RealType beta = mixer.beta;
    MorseType variant = mixer.variant;
    
    // V(r) = D_e exp(-a(r-re)(exp(-a(r-re))-2)
    
    RealType expt     = -beta*( *(idat.rij) - Re);
    RealType expfnc   = exp(expt);
    RealType expfnc2  = expfnc*expfnc;
    
    RealType exptC = 0.0;
    RealType expfncC = 0.0;
    RealType expfnc2C = 0.0;
    
    if (idat.shiftedPot || idat.shiftedForce) {
      exptC     = -beta*( *(idat.rcut) - Re);
      expfncC   = exp(exptC);
      expfnc2C  = expfncC*expfncC;
    }
    
    
    switch(variant) {
    case mtShifted : {
      
      myPot  = De * (expfnc2  - 2.0 * expfnc);
      myDeriv   = 2.0 * De * beta * (expfnc - expfnc2);
      
      if (idat.shiftedPot) {
        myPotC = De * (expfnc2C - 2.0 * expfncC);
        myDerivC = 0.0;
      } else if (idat.shiftedForce) {
        myPotC = De * (expfnc2C - 2.0 * expfncC);
        myDerivC  = 2.0 * De * beta * (expfncC - expfnc2C);
        myPotC += myDerivC * ( *(idat.rij) - *(idat.rcut) );
      } else {
        myPotC = 0.0;
        myDerivC = 0.0;
      }
      
      break;
    }
    case mtRepulsive : {
        
      myPot  = De * expfnc2;
      myDeriv  = -2.0 * De * beta * expfnc2;
      
      if (idat.shiftedPot) {
        myPotC = De * expfnc2C;
        myDerivC = 0.0;
      } else if (idat.shiftedForce) {
        myPotC = De * expfnc2C;
        myDerivC = -2.0 * De * beta * expfnc2C;
        myPotC += myDerivC * ( *(idat.rij) - *(idat.rcut));
      } else {
        myPotC = 0.0;
        myDerivC = 0.0;
      }
      
      break;
    }
    case mtUnknown: {
      // don't know what to do so don't do anything
      break;
    }
    }
    
    RealType pot_temp = *(idat.vdwMult) * (myPot - myPotC);
    *(idat.vpair) += pot_temp;
    
    RealType dudr = *(idat.sw) * *(idat.vdwMult) * (myDeriv - myDerivC);
    
    (*(idat.pot))[VANDERWAALS_FAMILY] += *(idat.sw) * pot_temp;
    *(idat.f1) = *(idat.d) * dudr / *(idat.rij);
    
    return;    
  }

  RealType Morse::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();   
    
    int atid1 = atypes.first->getIdent();
    int atid2 = atypes.second->getIdent();
    int mtid1 = Mtids[atid1];
    int mtid2 = Mtids[atid2];
    
    if ( mtid1 == -1 || mtid2 == -1) return 0.0;
    else {      
      MorseInteractionData mixer = MixingMap[mtid1][mtid2];
      RealType Re = mixer.Re;
      RealType beta = mixer.beta;
      // This value of the r corresponds to an energy about 1.48% of 
      // the energy at the bottom of the Morse well.  For comparison, the
      // Lennard-Jones function is about 1.63% of it's minimum value at
      // a distance of 2.5 sigma.
      return (4.9 + beta * Re) / beta;
    }
  }
}

