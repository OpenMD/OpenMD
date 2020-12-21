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

#include "nonbonded/RepulsivePower.hpp"
#include "utils/simError.h"
#include "types/RepulsivePowerInteractionType.hpp"

using namespace std;

namespace OpenMD {

  RepulsivePower::RepulsivePower() : initialized_(false), forceField_(NULL), 
				     name_("RepulsivePower") {}
  
  void RepulsivePower::initialize() {    

    RPtypes.clear();
    RPtids.clear();
    MixingMap.clear();
    RPtids.resize( forceField_->getNAtomType(), -1);

    ForceField::NonBondedInteractionTypeContainer* nbiTypes = forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    ForceField::NonBondedInteractionTypeContainer::KeyType keys;
    NonBondedInteractionType* nbt;
    int rptid1, rptid2;

    for (nbt = nbiTypes->beginType(j); nbt != NULL; 
         nbt = nbiTypes->nextType(j)) {
      
      if (nbt->isRepulsivePower()) {
        keys = nbiTypes->getKeys(j);
        AtomType* at1 = forceField_->getAtomType(keys[0]);
        if (at1 == NULL) {
          sprintf( painCave.errMsg,
                   "RepulsivePower::initialize could not find AtomType %s\n"
                   "\tto for for %s - %s interaction.\n",
                   keys[0].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }

        AtomType* at2 = forceField_->getAtomType(keys[1]);
        if (at2 == NULL) {
          sprintf( painCave.errMsg,
                   "RepulsivePower::initialize could not find AtomType %s\n"
                   "\tfor %s - %s nonbonded interaction.\n",
                   keys[1].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }

        int atid1 = at1->getIdent();
        if (RPtids[atid1] == -1) {
          rptid1 = RPtypes.size();
          RPtypes.insert(atid1);
          RPtids[atid1] = rptid1;         
        }
        int atid2 = at2->getIdent();
        if (RPtids[atid2] == -1) {
          rptid2 = RPtypes.size();
          RPtypes.insert(atid2);
          RPtids[atid2] = rptid2;         
        }
          
        RepulsivePowerInteractionType* rpit = dynamic_cast<RepulsivePowerInteractionType*>(nbt);
        if (rpit == NULL) {
          sprintf( painCave.errMsg,
                   "RepulsivePower::initialize could not convert NonBondedInteractionType\n"
                   "\tto RepulsivePowerInteractionType for %s - %s interaction.\n", 
                   at1->getName().c_str(),
                   at2->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
        
        RealType sigma = rpit->getSigma();
        RealType epsilon = rpit->getEpsilon();
        int nRep = rpit->getNrep();

        addExplicitInteraction(at1, at2, sigma, epsilon, nRep);
      }
    }
    initialized_ = true;
  }
      
  void RepulsivePower::addExplicitInteraction(AtomType* atype1, 
                                              AtomType* atype2, 
                                              RealType sigma, 
                                              RealType epsilon, 
                                              int nRep) {

    RPInteractionData mixer;
    mixer.sigma = sigma;
    mixer.epsilon = epsilon;
    mixer.sigmai = 1.0 / mixer.sigma;
    mixer.nRep = nRep;

    int rptid1 = RPtids[atype1->getIdent()];
    int rptid2 = RPtids[atype2->getIdent()];
    int nRP = RPtypes.size();

    MixingMap.resize(nRP);
    MixingMap[rptid1].resize(nRP);
    
    MixingMap[rptid1][rptid2] = mixer;
    if (rptid2 != rptid1) {
      MixingMap[rptid2].resize(nRP);
      MixingMap[rptid2][rptid1] = mixer;
    }    
  }
  
  void RepulsivePower::calcForce(InteractionData &idat) {

    if (!initialized_) initialize();
    
    RPInteractionData &mixer = MixingMap[RPtids[idat.atid1]][RPtids[idat.atid2]];
    RealType sigmai = mixer.sigmai;
    RealType epsilon = mixer.epsilon;
    int nRep = mixer.nRep;
    
    RealType ros;
    RealType rcos;
    RealType myPot = 0.0;
    RealType myPotC = 0.0;
    RealType myDeriv = 0.0;
    RealType myDerivC = 0.0;
    
    ros = idat.rij * sigmai;     
    
    getNRepulsionFunc(ros, nRep, myPot, myDeriv);
    
    if (idat.shiftedPot) {
      rcos = idat.rcut * sigmai;
      getNRepulsionFunc(rcos, nRep, myPotC, myDerivC);
      myDerivC = 0.0;
    } else if (idat.shiftedForce) {
      rcos = idat.rcut * sigmai;
      getNRepulsionFunc(rcos, nRep, myPotC, myDerivC);
      myPotC = myPotC + myDerivC * (idat.rij - idat.rcut) * sigmai;
    } else {
      myPotC = 0.0;
      myDerivC = 0.0;        
    }
    
    RealType pot_temp = idat.vdwMult * epsilon * (myPot - myPotC);
    idat.vpair += pot_temp;
    
    RealType dudr = idat.sw * idat.vdwMult * epsilon * (myDeriv - 
							myDerivC)*sigmai;
    
    idat.pot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;
    if (idat.isSelected)
      idat.selePot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;

    idat.f1 += idat.d * dudr / idat.rij;
    
    return;
  }

  void RepulsivePower::getNRepulsionFunc(const RealType &r, int &n,
                                         RealType &pot, RealType &deriv) {

    RealType ri = 1.0 / r;
    RealType rin = pow(ri, n);
    RealType rin1 = rin * ri;

    pot = rin;
    deriv = -n * rin1;
    
    return;
  }
  

  RealType RepulsivePower::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {   
    if (!initialized_) initialize();   
    
    int atid1 = atypes.first->getIdent();
    int atid2 = atypes.second->getIdent();
    int rptid1 = RPtids[atid1];
    int rptid2 = RPtids[atid2];
    
    if (rptid1 == -1 || rptid2 == -1) return 0.0;
    else {      
      RPInteractionData mixer = MixingMap[rptid1][rptid2];
      return 2.5 * mixer.sigma;
    }
  }  
}

