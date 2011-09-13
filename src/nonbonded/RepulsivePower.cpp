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
#include "nonbonded/RepulsivePower.hpp"
#include "utils/simError.h"
#include "types/NonBondedInteractionType.hpp"

using namespace std;

namespace OpenMD {

  RepulsivePower::RepulsivePower() : name_("RepulsivePower"), 
                                     initialized_(false), forceField_(NULL) {}
  
  void RepulsivePower::initialize() {    

    ForceField::NonBondedInteractionTypeContainer* nbiTypes = forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    NonBondedInteractionType* nbt;

    for (nbt = nbiTypes->beginType(j); nbt != NULL; 
         nbt = nbiTypes->nextType(j)) {
      
      if (nbt->isRepulsivePower()) {
        
        pair<AtomType*, AtomType*> atypes = nbt->getAtomTypes();
        
        GenericData* data = nbt->getPropertyByName("RepulsivePower");
        if (data == NULL) {
          sprintf( painCave.errMsg, "RepulsivePower::initialize could not\n"
                   "\tfind RepulsivePower parameters for %s - %s interaction.\n", 
                   atypes.first->getName().c_str(),
                   atypes.second->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError(); 
        }
        
        RepulsivePowerData* rpData = dynamic_cast<RepulsivePowerData*>(data);
        if (rpData == NULL) {
          sprintf( painCave.errMsg,
                   "RepulsivePower::initialize could not convert GenericData\n"
                   "\tto RepulsivePowerData for %s - %s interaction.\n", 
                   atypes.first->getName().c_str(),
                   atypes.second->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }
        
        RepulsivePowerParam rpParam = rpData->getData();

        RealType sigma = rpParam.sigma;
        RealType epsilon = rpParam.epsilon;
        int nRep = rpParam.nRep;

        addExplicitInteraction(atypes.first, atypes.second, 
                               sigma, epsilon, nRep);
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

    pair<AtomType*, AtomType*> key1, key2;
    key1 = make_pair(atype1, atype2);
    key2 = make_pair(atype2, atype1);
    
    MixingMap[key1] = mixer;
    if (key2 != key1) {
      MixingMap[key2] = mixer;
    }    
  }
  
  void RepulsivePower::calcForce(InteractionData &idat) {

    if (!initialized_) initialize();
    
    map<pair<AtomType*, AtomType*>, RPInteractionData>::iterator it;
    it = MixingMap.find( idat.atypes );

    if (it != MixingMap.end()) {

      RPInteractionData mixer = (*it).second;
      RealType sigmai = mixer.sigmai;
      RealType epsilon = mixer.epsilon;
      int nRep = mixer.nRep;
      
      RealType ros;
      RealType rcos;
      RealType myPot = 0.0;
      RealType myPotC = 0.0;
      RealType myDeriv = 0.0;
      RealType myDerivC = 0.0;
     
      ros = *(idat.rij) * sigmai;     
      
      getNRepulsionFunc(ros, nRep, myPot, myDeriv);
      
      if (idat.shiftedPot) {
        rcos = *(idat.rcut) * sigmai;
        getNRepulsionFunc(rcos, nRep, myPotC, myDerivC);
        myDerivC = 0.0;
      } else if (idat.shiftedForce) {
        rcos = *(idat.rcut) * sigmai;
        getNRepulsionFunc(rcos, nRep, myPotC, myDerivC);
        myPotC = myPotC + myDerivC * (*(idat.rij) - *(idat.rcut)) * sigmai;
      } else {
        myPotC = 0.0;
        myDerivC = 0.0;        
      }

      RealType pot_temp = *(idat.vdwMult) * epsilon * (myPot - myPotC);
      *(idat.vpair) += pot_temp;
      
      RealType dudr = *(idat.sw) * *(idat.vdwMult) * epsilon * (myDeriv - 
                                                                myDerivC)*sigmai;      

      (*(idat.pot))[VANDERWAALS_FAMILY] += *(idat.sw) * pot_temp;
      *(idat.f1) = *(idat.d) * dudr / *(idat.rij);
    }
    return;
  }

  void RepulsivePower::getNRepulsionFunc(RealType r, int n, RealType &pot, RealType &deriv) {

    RealType ri = 1.0 / r;
    RealType rin = pow(ri, n);
    RealType rin1 = rin * ri;

    pot = rin;
    deriv = -n * rin1;
    
    return;
  }
  

  RealType RepulsivePower::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {   
    if (!initialized_) initialize();   
    map<pair<AtomType*, AtomType*>, RPInteractionData>::iterator it;
    it = MixingMap.find(atypes);
    if (it == MixingMap.end()) 
      return 0.0;
    else  {
      RPInteractionData mixer = (*it).second;
      return 2.5 * mixer.sigma;
    }
  }

}

