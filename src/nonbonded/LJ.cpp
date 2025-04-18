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

#include "nonbonded/LJ.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>

#include "types/LennardJonesAdapter.hpp"
#include "types/LennardJonesInteractionType.hpp"
#include "utils/simError.h"

namespace OpenMD {

  LJ::LJ() : initialized_(false), forceField_(NULL), name_("LJ") {}

  RealType LJ::getSigma(AtomType* atomType1, AtomType* atomType2) {
    LennardJonesAdapter lja1 = LennardJonesAdapter(atomType1);
    LennardJonesAdapter lja2 = LennardJonesAdapter(atomType2);
    RealType sigma1          = lja1.getSigma();
    RealType sigma2          = lja2.getSigma();

    ForceFieldOptions& fopts = forceField_->getForceFieldOptions();
    string DistanceMix       = fopts.getDistanceMixingRule();
    toUpper(DistanceMix);

    if (DistanceMix == "GEOMETRIC")
      return sqrt(sigma1 * sigma2);
    else
      return 0.5 * (sigma1 + sigma2);
  }

  RealType LJ::getEpsilon(AtomType* atomType1, AtomType* atomType2) {
    LennardJonesAdapter lja1 = LennardJonesAdapter(atomType1);
    LennardJonesAdapter lja2 = LennardJonesAdapter(atomType2);

    RealType epsilon1 = lja1.getEpsilon();
    RealType epsilon2 = lja2.getEpsilon();
    return sqrt(epsilon1 * epsilon2);
  }

  void LJ::initialize() {
    LJtypes.clear();
    LJtids.clear();
    MixingMap.clear();
    nLJ_ = 0;

    LJtids.resize(forceField_->getNAtomType(), -1);

    AtomTypeSet::iterator at;
    for (at = simTypes_.begin(); at != simTypes_.end(); ++at) {
      if ((*at)->isLennardJones()) nLJ_++;
    }

    MixingMap.resize(nLJ_);

    for (at = simTypes_.begin(); at != simTypes_.end(); ++at) {
      if ((*at)->isLennardJones()) addType(*at);
    }

    ForceField::NonBondedInteractionTypeContainer* nbiTypes =
        forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    NonBondedInteractionType* nbt;
    ForceField::NonBondedInteractionTypeContainer::KeyType keys;

    for (nbt = nbiTypes->beginType(j); nbt != NULL;
         nbt = nbiTypes->nextType(j)) {
      if (nbt->isLennardJones()) {
        keys          = nbiTypes->getKeys(j);
        keys          = nbiTypes->getKeys(j);
        AtomType* at1 = forceField_->getAtomType(keys[0]);
        if (at1 == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "LennardJones::initialize could not find AtomType %s\n"
                   "\tto for for %s - %s interaction.\n",
                   keys[0].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        AtomType* at2 = forceField_->getAtomType(keys[1]);
        if (at2 == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "LennardJones::initialize could not find AtomType %s\n"
                   "\tfor %s - %s nonbonded interaction.\n",
                   keys[1].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        LennardJonesInteractionType* ljit =
            dynamic_cast<LennardJonesInteractionType*>(nbt);

        if (ljit == NULL) {
          snprintf(
              painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
              "LJ::initialize could not convert NonBondedInteractionType\n"
              "\tto LennardJonesInteractionType for %s - %s interaction.\n",
              at1->getName().c_str(), at2->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        RealType sigma   = ljit->getSigma();
        RealType epsilon = ljit->getEpsilon();
        addExplicitInteraction(at1, at2, sigma, epsilon);
      }
    }
    initialized_ = true;
  }

  void LJ::addType(AtomType* atomType) {
    // add it to the map:
    int atid  = atomType->getIdent();
    int ljtid = LJtypes.size();

    pair<set<int>::iterator, bool> ret;
    ret = LJtypes.insert(atid);
    if (ret.second == false) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "LJ already had a previous entry with ident %d\n", atid);
      painCave.severity = OPENMD_INFO;
      painCave.isFatal  = 0;
      simError();
    }

    // Check to make sure the 1/sigma won't cause problems later:
    RealType s = getSigma(atomType, atomType);
    if (fabs(s) < std::numeric_limits<RealType>::epsilon()) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "Lennard-Jones atom %s was defined with a sigma value (%f)\n"
               "\tthat was too close to zero.",
               atomType->getName().c_str(), s);
      painCave.severity = OPENMD_ERROR;
      painCave.isFatal  = 1;
      simError();
    }

    LJtids[atid] = ljtid;
    MixingMap[ljtid].resize(nLJ_);

    // Now, iterate over all known types and add to the mixing map:

    std::set<int>::iterator it;
    for (it = LJtypes.begin(); it != LJtypes.end(); ++it) {
      int ljtid2       = LJtids[(*it)];
      AtomType* atype2 = forceField_->getAtomType((*it));

      LJInteractionData mixer;
      mixer.sigma         = getSigma(atomType, atype2);
      mixer.epsilon       = getEpsilon(atomType, atype2);
      mixer.sigmai        = 1.0 / mixer.sigma;
      mixer.explicitlySet = false;
      MixingMap[ljtid2].resize(nLJ_);

      MixingMap[ljtid][ljtid2] = mixer;
      if (ljtid2 != ljtid) { MixingMap[ljtid2][ljtid] = mixer; }
    }
  }

  void LJ::addExplicitInteraction(AtomType* atype1, AtomType* atype2,
                                  RealType sigma, RealType epsilon) {
    LJInteractionData mixer;
    mixer.sigma         = sigma;
    mixer.epsilon       = epsilon;
    mixer.sigmai        = 1.0 / mixer.sigma;
    mixer.explicitlySet = true;

    int atid1 = atype1->getIdent();
    int atid2 = atype2->getIdent();

    int ljtid1, ljtid2;

    pair<set<int>::iterator, bool> ret;
    ret = LJtypes.insert(atid1);
    if (ret.second == false) {
      // already had this type in the LJMap, just get the ljtid:
      ljtid1 = LJtids[atid1];
    } else {
      // didn't already have it, so make a new one and assign it:
      ljtid1        = nLJ_;
      LJtids[atid1] = nLJ_;
      nLJ_++;
    }

    ret = LJtypes.insert(atid2);
    if (ret.second == false) {
      // already had this type in the LJMap, just get the ljtid:
      ljtid2 = LJtids[atid2];
    } else {
      // didn't already have it, so make a new one and assign it:
      ljtid2        = nLJ_;
      LJtids[atid2] = nLJ_;
      nLJ_++;
    }

    MixingMap.resize(nLJ_);
    MixingMap[ljtid1].resize(nLJ_);

    MixingMap[ljtid1][ljtid2] = mixer;
    if (ljtid2 != ljtid1) {
      MixingMap[ljtid2].resize(nLJ_);
      MixingMap[ljtid2][ljtid1] = mixer;
    }
  }

  void LJ::calcForce(InteractionData& idat) {
    if (!initialized_) initialize();

    LJInteractionData& mixer =
        MixingMap[LJtids[idat.atid1]][LJtids[idat.atid2]];

    RealType sigmai  = mixer.sigmai;
    RealType epsilon = mixer.epsilon;

    RealType ros;
    RealType rcos;
    RealType myPot    = 0.0;
    RealType myPotC   = 0.0;
    RealType myDeriv  = 0.0;
    RealType myDerivC = 0.0;

    ros = idat.rij * sigmai;

    getLJfunc(ros, myPot, myDeriv);

    if (idat.shiftedPot) {
      rcos = idat.rcut * sigmai;
      getLJfunc(rcos, myPotC, myDerivC);
      myDerivC = 0.0;
    } else if (idat.shiftedForce) {
      rcos = idat.rcut * sigmai;
      getLJfunc(rcos, myPotC, myDerivC);
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

  void LJ::getLJfunc(RealType r, RealType& pot, RealType& deriv) {
    RealType ri   = 1.0 / r;
    RealType ri2  = ri * ri;
    RealType ri6  = ri2 * ri2 * ri2;
    RealType ri7  = ri6 * ri;
    RealType ri12 = ri6 * ri6;
    RealType ri13 = ri12 * ri;

    pot   = 4.0 * (ri12 - ri6);
    deriv = 24.0 * (ri7 - 2.0 * ri13);

    return;
  }

  RealType LJ::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();

    int atid1  = atypes.first->getIdent();
    int atid2  = atypes.second->getIdent();
    int ljtid1 = LJtids[atid1];
    int ljtid2 = LJtids[atid2];

    if (ljtid1 == -1 || ljtid2 == -1)
      return 0.0;
    else {
      LJInteractionData mixer = MixingMap[ljtid1][ljtid2];
      return 2.5 * mixer.sigma;
    }
  }

}  // namespace OpenMD
