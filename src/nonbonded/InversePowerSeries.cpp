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

#include "nonbonded/InversePowerSeries.hpp"
#include "utils/simError.h"
#include "types/InversePowerSeriesInteractionType.hpp"

using namespace std;

namespace OpenMD {

  InversePowerSeries::InversePowerSeries() : initialized_(false),
                                             forceField_(NULL), 
                                             name_("InversePowerSeries") {}
  
  void InversePowerSeries::initialize() {    
    
    InversePowerSeriesTypes.clear();
    InversePowerSeriesTids.clear();
    MixingMap.clear();
    InversePowerSeriesTids.resize( forceField_->getNAtomType(), -1);

    ForceField::NonBondedInteractionTypeContainer* nbiTypes = forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    ForceField::NonBondedInteractionTypeContainer::KeyType keys;
    NonBondedInteractionType* nbt;
    int ipstid1, ipstid2;

    for (nbt = nbiTypes->beginType(j); nbt != NULL; 
         nbt = nbiTypes->nextType(j)) {
      
      if (nbt->isInversePowerSeries()) {
        keys = nbiTypes->getKeys(j);
        AtomType* at1 = forceField_->getAtomType(keys[0]);
        if (at1 == NULL) {
          sprintf( painCave.errMsg,
                   "InversePowerSeries::initialize could not find AtomType %s\n"
                   "\tto for for %s - %s interaction.\n",
                   keys[0].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }

        AtomType* at2 = forceField_->getAtomType(keys[1]);
        if (at2 == NULL) {
          sprintf( painCave.errMsg,
                   "InversePowerSeries::initialize could not find AtomType %s\n"
                   "\tfor %s - %s nonbonded interaction.\n",
                   keys[1].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }

        int atid1 = at1->getIdent();
        if (InversePowerSeriesTids[atid1] == -1) {
          ipstid1 = InversePowerSeriesTypes.size();
          InversePowerSeriesTypes.insert(atid1);
          InversePowerSeriesTids[atid1] = ipstid1;         
        }
        int atid2 = at2->getIdent();
        if (InversePowerSeriesTids[atid2] == -1) {
          ipstid2 = InversePowerSeriesTypes.size();
          InversePowerSeriesTypes.insert(atid2);
          InversePowerSeriesTids[atid2] = ipstid2;         
        }
          
        InversePowerSeriesInteractionType* ipsit = dynamic_cast<InversePowerSeriesInteractionType*>(nbt);
        if (ipsit == NULL) {
          sprintf( painCave.errMsg,
                   "InversePowerSeries::initialize could not convert NonBondedInteractionType\n"
                   "\tto InversePowerSeriesInteractionType for %s - %s interaction.\n", 
                   at1->getName().c_str(),
                   at2->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal = 1;
          simError();          
        }

        std::vector<int> powers = ipsit->getPowers();
        std::vector<RealType> coefficients = ipsit->getCoefficients();

        addExplicitInteraction(at1, at2, powers, coefficients);
      }
    }
    initialized_ = true;
  }
  
  void InversePowerSeries::addExplicitInteraction(AtomType* atype1, 
                                                  AtomType* atype2, 
                                                  std::vector<int> powers, 
                                                  std::vector<RealType> coefficients) {

    InversePowerSeriesInteractionData mixer;
    mixer.powers = powers;
    mixer.coefficients = coefficients;

    int ipstid1 = InversePowerSeriesTids[atype1->getIdent()];
    int ipstid2 = InversePowerSeriesTids[atype2->getIdent()];
    int nInversePowerSeries = InversePowerSeriesTypes.size();

    MixingMap.resize(nInversePowerSeries);
    MixingMap[ipstid1].resize(nInversePowerSeries);
    
    MixingMap[ipstid1][ipstid2] = mixer;
    if (ipstid2 != ipstid1) {
      MixingMap[ipstid2].resize(nInversePowerSeries);
      MixingMap[ipstid2][ipstid1] = mixer;
    }    
  }
  
  void InversePowerSeries::calcForce(InteractionData &idat) {

    if (!initialized_) initialize();
    
    InversePowerSeriesInteractionData &mixer = MixingMap[InversePowerSeriesTids[idat.atid1]][InversePowerSeriesTids[idat.atid2]];
    std::vector<int> powers = mixer.powers;
    std::vector<RealType> coefficients = mixer.coefficients;

    RealType myPot = 0.0;
    RealType myPotC = 0.0;
    RealType myDeriv = 0.0;
    RealType myDerivC = 0.0;

    RealType ri  = 1.0 / idat.rij;
    RealType ric = 1.0 / idat.rcut;

    RealType fn, fnc;
    
    for (unsigned int i = 0; i < powers.size(); i++) {
      fn  = coefficients[i] * pow(ri,  powers[i]);
      fnc = coefficients[i] * pow(ric, powers[i]);
      myPot  += fn;
      myPotC += fnc;
      myDeriv -= powers[i] * fn * ri;
      myDerivC -= powers[i] * fnc * ric;
    }
       
    if (idat.shiftedPot) {
      myDerivC = 0.0;
    } else if (idat.shiftedForce) {      
      myPotC = myPotC + myDerivC * ( idat.rij - idat.rcut );
    } else {
      myPotC = 0.0;
      myDerivC = 0.0;        
    }
    
    RealType pot_temp = idat.vdwMult * (myPot - myPotC);
    idat.vpair += pot_temp;
    
    RealType dudr = idat.sw * idat.vdwMult * (myDeriv -  myDerivC);
    
    idat.pot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;
    if (idat.isSelected)
      idat.selePot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;

    idat.f1 += idat.d * dudr / idat.rij;
    
    return;
  }
  
  RealType InversePowerSeries::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {   
    if (!initialized_) initialize();   
    
    int atid1 = atypes.first->getIdent();
    int atid2 = atypes.second->getIdent();
    int ipstid1 = InversePowerSeriesTids[atid1];
    int ipstid2 = InversePowerSeriesTids[atid2];
    
    if (ipstid1 == -1 || ipstid2 == -1) return 0.0;
    else {
      // This seems to work moderately well as a default.  There's no
      // inherent scale for 1/r^n interactions that we can standardize.
      // 12 angstroms seems to be a reasonably good guess for most
      // cases.
      return 12.0;
    }
  }  
}

