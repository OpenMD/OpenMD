/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "nonbonded/RepulsivePower.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>

#include "types/RepulsivePowerInteractionType.hpp"
#include "utils/simError.h"

using namespace std;

namespace OpenMD {

  RepulsivePower::RepulsivePower() :
      initialized_(false), forceField_(NULL), name_("RepulsivePower") {}

  void RepulsivePower::initialize() {
    RPtypes.clear();
    RPtids.clear();
    MixingMap.clear();
    RPtids.resize(forceField_->getNAtomType(), -1);

    ForceField::NonBondedInteractionTypeContainer* nbiTypes =
        forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    ForceField::NonBondedInteractionTypeContainer::KeyType keys;
    NonBondedInteractionType* nbt;
    int rptid1, rptid2;

    for (nbt = nbiTypes->beginType(j); nbt != NULL;
         nbt = nbiTypes->nextType(j)) {
      if (nbt->isRepulsivePower()) {
        keys          = nbiTypes->getKeys(j);
        AtomType* at1 = forceField_->getAtomType(keys[0]);
        if (at1 == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "RepulsivePower::initialize could not find AtomType %s\n"
                   "\tto for for %s - %s interaction.\n",
                   keys[0].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        AtomType* at2 = forceField_->getAtomType(keys[1]);
        if (at2 == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "RepulsivePower::initialize could not find AtomType %s\n"
                   "\tfor %s - %s nonbonded interaction.\n",
                   keys[1].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
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

        RepulsivePowerInteractionType* rpit =
            dynamic_cast<RepulsivePowerInteractionType*>(nbt);
        if (rpit == NULL) {
          snprintf(
              painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
              "RepulsivePower::initialize could not convert "
              "NonBondedInteractionType\n"
              "\tto RepulsivePowerInteractionType for %s - %s interaction.\n",
              at1->getName().c_str(), at2->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        RealType sigma   = rpit->getSigma();
        RealType epsilon = rpit->getEpsilon();
        int nRep         = rpit->getNrep();

        addExplicitInteraction(at1, at2, sigma, epsilon, nRep);
      }
    }
    initialized_ = true;
  }

  void RepulsivePower::addExplicitInteraction(AtomType* atype1,
                                              AtomType* atype2, RealType sigma,
                                              RealType epsilon, int nRep) {
    RPInteractionData mixer;
    mixer.sigma   = sigma;
    mixer.epsilon = epsilon;
    mixer.sigmai  = 1.0 / mixer.sigma;
    mixer.nRep    = nRep;

    int nRP   = RPtypes.size();
    int atid1 = atype1->getIdent();
    int atid2 = atype2->getIdent();
    int rptid1, rptid2;

    pair<set<int>::iterator, bool> ret;
    ret = RPtypes.insert(atid1);
    if (ret.second == false) {
      // already had this type in the MieMap, just get the mietid:
      rptid1 = RPtids[atid1];
    } else {
      // didn't already have it, so make a new one and assign it:
      rptid1        = nRP;
      RPtids[atid1] = nRP;
      nRP++;
    }

    ret = RPtypes.insert(atid2);
    if (ret.second == false) {
      // already had this type in the MieMap, just get the mietid:
      rptid2 = RPtids[atid2];
    } else {
      // didn't already have it, so make a new one and assign it:
      rptid2        = nRP;
      RPtids[atid2] = nRP;
      nRP++;
    }

    MixingMap.resize(nRP);
    MixingMap[rptid1].resize(nRP);

    MixingMap[rptid1][rptid2] = mixer;
    if (rptid2 != rptid1) {
      MixingMap[rptid2].resize(nRP);
      MixingMap[rptid2][rptid1] = mixer;
    }
  }

  void RepulsivePower::calcForce(InteractionData& idat) {
    if (!initialized_) initialize();

    RPInteractionData& mixer =
        MixingMap[RPtids[idat.atid1]][RPtids[idat.atid2]];
    RealType sigmai  = mixer.sigmai;
    RealType epsilon = mixer.epsilon;
    int nRep         = mixer.nRep;

    RealType ros;
    RealType rcos;
    RealType myPot    = 0.0;
    RealType myPotC   = 0.0;
    RealType myDeriv  = 0.0;
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
      myPotC   = 0.0;
      myDerivC = 0.0;
    }

    RealType pot_temp = idat.vdwMult * epsilon * (myPot - myPotC);
    idat.vpair += pot_temp;

    RealType dudr =
        idat.sw * idat.vdwMult * epsilon * (myDeriv - myDerivC) * sigmai;

    idat.pot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;
    if (idat.isSelected) idat.selePot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;

    idat.f1 += idat.d * dudr / idat.rij;

    return;
  }

  void RepulsivePower::getNRepulsionFunc(const RealType& r, int& n,
                                         RealType& pot, RealType& deriv) {
    RealType ri   = 1.0 / r;
    RealType rin  = pow(ri, n);
    RealType rin1 = rin * ri;

    pot   = rin;
    deriv = -n * rin1;

    return;
  }

  RealType RepulsivePower::getSuggestedCutoffRadius(
      pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();

    int atid1  = atypes.first->getIdent();
    int atid2  = atypes.second->getIdent();
    int rptid1 = RPtids[atid1];
    int rptid2 = RPtids[atid2];

    if (rptid1 == -1 || rptid2 == -1)
      return 0.0;
    else {
      RPInteractionData mixer = MixingMap[rptid1][rptid2];
      return 2.5 * mixer.sigma;
    }
  }
}  // namespace OpenMD
