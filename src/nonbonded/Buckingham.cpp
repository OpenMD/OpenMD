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

#include <cstdio>
#include <cstring>
#include <cmath>

#include "nonbonded/Buckingham.hpp"
#include "utils/simError.h"
#include "types/BuckinghamInteractionType.hpp"

using namespace std;

namespace OpenMD {

  Buckingham::Buckingham() : initialized_(false), forceField_(NULL), name_("Buckingham") {}
  
  void Buckingham::initialize() {    

    Btypes.clear();
    Btids.clear();
    MixingMap.clear();
    Btids.resize( forceField_->getNAtomType(), -1);

    ForceField::NonBondedInteractionTypeContainer* nbiTypes = forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    ForceField::NonBondedInteractionTypeContainer::KeyType keys;
    NonBondedInteractionType* nbt;
    int btid1, btid2;

    for (nbt = nbiTypes->beginType(j); nbt != NULL; 
         nbt = nbiTypes->nextType(j)) {
      
      if (nbt->isBuckingham()) {
        keys = nbiTypes->getKeys(j);
        AtomType* at1 = forceField_->getAtomType(keys[0]);
        if (at1 == NULL) {
          sprintf( painCave.errMsg,
                   "Buckingham::initialize could not find AtomType %s\n"
                   "\tto for for %s - %s interaction.\n",
                   keys[0].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }

        AtomType* at2 = forceField_->getAtomType(keys[1]);
        if (at2 == NULL) {
          sprintf( painCave.errMsg,
                   "Buckingham::initialize could not find AtomType %s\n"
                   "\tfor %s - %s nonbonded interaction.\n",
                   keys[1].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }

        int atid1 = at1->getIdent();
        if (Btids[atid1] == -1) {
          btid1 = Btypes.size();
          Btypes.insert(atid1);
          Btids[atid1] = btid1;
        }
        int atid2 = at2->getIdent();
        if (Btids[atid2] == -1) {
          btid2 = Btypes.size();
          Btypes.insert(atid2);
          Btids[atid2] = btid2;
        }
    
        BuckinghamInteractionType* bit = dynamic_cast<BuckinghamInteractionType*>(nbt);

        if (bit == NULL) {
          sprintf( painCave.errMsg,
                   "Buckingham::initialize could not convert NonBondedInteractionType\n"
                   "\tto BuckinghamInteractionType for %s - %s interaction.\n", 
                   at1->getName().c_str(),
                   at2->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
        
        RealType A = bit->getA();
        RealType B = bit->getB();
        RealType C = bit->getC();
   
        BuckinghamType variant = bit->getInteractionType();
        addExplicitInteraction(at1, at2, A, B, C, variant );
      }
    }  
    initialized_ = true;
  }
      
  void Buckingham::addExplicitInteraction(AtomType* atype1, AtomType* atype2, 
                                          RealType A, RealType B, RealType C, 
                                          BuckinghamType bt) {

    BuckinghamInteractionData mixer;
    mixer.A = A;
    mixer.B = B;
    mixer.C = C;
    mixer.variant = bt;

    int btid1 = Btids[atype1->getIdent()];
    int btid2 = Btids[atype2->getIdent()];
    int nB = Btypes.size();

    MixingMap.resize(nB);
    MixingMap[btid1].resize(nB);
    
    MixingMap[btid1][btid2] = mixer;
    if (btid2 != btid1) {
      MixingMap[btid2].resize(nB);
      MixingMap[btid2][btid1] = mixer;
    }    
  }

  void Buckingham::addExplicitInteraction(AtomType* atype1, AtomType* atype2, 
                                          RealType A, RealType B, RealType C,
                                          RealType sigma, RealType epsilon,
                                          BuckinghamType bt) {
    
    BuckinghamInteractionData mixer;
    mixer.A = A;
    mixer.B = B;
    mixer.C = C;
    mixer.sigma = sigma;
    mixer.epsilon = epsilon;
    mixer.variant = bt;

    int btid1 = Btids[atype1->getIdent()];
    int btid2 = Btids[atype2->getIdent()];
    int nB = Btypes.size();

    MixingMap.resize(nB);
    MixingMap[btid1].resize(nB);
    
    MixingMap[btid1][btid2] = mixer;
    if (btid2 != btid1) {
      MixingMap[btid2].resize(nB);
      MixingMap[btid2][btid1] = mixer;
    }    
  }

  void Buckingham::calcForce(InteractionData &idat) {

    if (!initialized_) initialize();
   
    BuckinghamInteractionData &mixer = MixingMap[Btids[idat.atid1]][Btids[idat.atid2]];

    RealType myPot = 0.0;
    RealType myPotC = 0.0;
    RealType myDeriv = 0.0;
    RealType myDerivC = 0.0;
    
    RealType A = mixer.A;
    RealType B = mixer.B;
    RealType C = mixer.C;
    RealType sigma = mixer.sigma;
    RealType epsilon = mixer.epsilon;
    BuckinghamType variant = mixer.variant;
    
    RealType expt     = -B * idat.rij;
    RealType expfnc   = exp(expt);
    RealType fnc6  = 1.0 / pow(idat.rij, 6);
    RealType fnc7  = fnc6 / idat.rij;
    
    RealType exptC = 0.0;
    RealType expfncC = 0.0;
    RealType fnc6C = 0.0;
    RealType fnc7C = 0.0;
    
    if (idat.shiftedPot || idat.shiftedForce) {
      exptC     = -B * idat.rcut;
      expfncC   =  exp(exptC);
      fnc6C     =  1.0 / pow( idat.rcut, 6);
      fnc7C     =  fnc6C / idat.rcut;
    }
    
    switch(variant) {
    case btTraditional: {
      
      // V(r) = A exp(-B*r) - C/r^6 
      myPot  = A*expfnc - C * fnc6;
      myDeriv = - A * B * expfnc + C * fnc7;
      
      if (idat.shiftedPot) {
        myPotC = A*expfncC - C * fnc6C ;
        myDerivC = 0.0;
      } else if (idat.shiftedForce) {
        myPotC = A * expfncC - C * fnc6C;
        myDerivC =  -A * B * expfncC + C * fnc7C;
        myPotC += myDerivC * ( idat.rij - idat.rcut );
      } else {
        myPotC = 0.0;
        myDerivC = 0.0;
      }
      break;
    }
    case btModified: {
      RealType s6 = pow(sigma, 6);
      RealType s7 = pow(sigma, 7);
      RealType fnc30 = pow(sigma / idat.rij, 30);
      RealType fnc31 = fnc30 * sigma / idat.rij;
      RealType fnc30C = 0.0;
      RealType fnc31C = 0.0;
      
      if (idat.shiftedPot || idat.shiftedForce) {
        fnc30C     =  pow( sigma / idat.rcut, 30);
        fnc31C     =  fnc30C * sigma / idat.rcut;
      }
      
      // V(r) = A exp(-B*r) - C/r^6 + 4 epsilon ((sigma/r)^30 - (sigma/r)^6)
      myPot  = A*expfnc - C * fnc6 + 4.0 * epsilon * ( fnc30 - s6*fnc6 );
      myDeriv = - A * B * expfnc + C * fnc7 + 4.0 * epsilon * (-30.0*fnc31 + 6.0*s7*fnc7) / sigma;
      
      if (idat.shiftedPot) {
        myPotC  = A*expfncC - C * fnc6C + 4.0 * epsilon * ( fnc30C - s6*fnc6C );
        myDerivC = 0.0;
      } else if (idat.shiftedForce) {
        myPotC  = A*expfncC - C * fnc6C + 4.0 * epsilon * ( fnc30C - s6*fnc6C );
        myDeriv = - A * B * expfncC + C * fnc7C + 4.0 * epsilon * (-30.0*fnc31C + 6.0*s7*fnc7C) / sigma;
        myPotC += myDerivC * ( idat.rij - idat.rcut );
      } else {
        myPotC = 0.0;
        myDerivC = 0.0;
      }

      break;
    }
    case btUnknown: {
      // don't know what to do so don't do anything
      break;
    }
    }
      
    RealType pot_temp = idat.vdwMult * (myPot - myPotC);
    idat.vpair += pot_temp;
    
    RealType dudr = idat.sw * idat.vdwMult * (myDeriv - myDerivC);
    
    idat.pot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;
    if (idat.isSelected)
      idat.selePot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;

    idat.f1 += idat.d * dudr / idat.rij;

    return;    
  }

  RealType Buckingham::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();   
    
    int atid1 = atypes.first->getIdent();
    int atid2 = atypes.second->getIdent();
    int btid1 = Btids[atid1];
    int btid2 = Btids[atid2];
    
    if ( btid1 == -1 || btid2 == -1) return 0.0;
    else {      
      // Uncomment if we ever want to query the simulated atoms types
      // for a suggested cutoff:
      //
      // BuckinghamInteractionData mixer = MixingMap[btid1][btid2];
      //
      // suggested cutoff for most implementations of the BKS potential are
      // around 1 nm (10 angstroms):
      return 10.0;
    }
  }
}

