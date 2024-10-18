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

#include "nonbonded/Morse.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>

#include "types/MorseInteractionType.hpp"
#include "utils/simError.h"

using namespace std;

namespace OpenMD {

  Morse::Morse() : initialized_(false), forceField_(NULL), name_("Morse") {}

  void Morse::initialize() {
    Mtypes.clear();
    Mtids.clear();
    MixingMap.clear();
    nM_ = 0;

    Mtids.resize(forceField_->getNAtomType(), -1);

    ForceField::NonBondedInteractionTypeContainer* nbiTypes =
        forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    ForceField::NonBondedInteractionTypeContainer::KeyType keys;
    NonBondedInteractionType* nbt;

    for (nbt = nbiTypes->beginType(j); nbt != NULL;
         nbt = nbiTypes->nextType(j)) {
      if (nbt->isMorse()) {
        keys          = nbiTypes->getKeys(j);
        AtomType* at1 = forceField_->getAtomType(keys[0]);
        if (at1 == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Morse::initialize could not find AtomType %s\n"
                   "\tto for for %s - %s interaction.\n",
                   keys[0].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        AtomType* at2 = forceField_->getAtomType(keys[1]);
        if (at2 == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Morse::initialize could not find AtomType %s\n"
                   "\tfor %s - %s nonbonded interaction.\n",
                   keys[1].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        MorseInteractionType* mit = dynamic_cast<MorseInteractionType*>(nbt);

        if (mit == NULL) {
          snprintf(
              painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
              "Morse::initialize could not convert NonBondedInteractionType\n"
              "\tto MorseInteractionType for %s - %s interaction.\n",
              at1->getName().c_str(), at2->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        RealType De   = mit->getD();
        RealType Re   = mit->getR();
        RealType beta = mit->getBeta();

        MorseType variant = mit->getInteractionType();
        addExplicitInteraction(at1, at2, De, Re, beta, variant);
      }
    }
    initialized_ = true;
  }

  void Morse::addExplicitInteraction(AtomType* atype1, AtomType* atype2,
                                     RealType De, RealType Re, RealType beta,
                                     MorseType mt) {
    MorseInteractionData mixer;
    mixer.De      = De;
    mixer.Re      = Re;
    mixer.beta    = beta;
    mixer.variant = mt;

    int atid1 = atype1->getIdent();
    int atid2 = atype2->getIdent();

    int mtid1, mtid2;

    pair<set<int>::iterator, bool> ret;
    ret = Mtypes.insert(atid1);
    if (ret.second == false) {
      // already had this type in the Mtypes list, just get the mtid:
      mtid1 = Mtids[atid1];
    } else {
      // didn't already have it, so make a new one and assign it:
      mtid1        = nM_;
      Mtids[atid1] = nM_;
      nM_++;
    }
    ret = Mtypes.insert(atid2);
    if (ret.second == false) {
      // already had this type in the Mtypes list, just get the mtid:
      mtid2 = Mtids[atid2];
    } else {
      // didn't already have it, so make a new one and assign it:
      mtid2        = nM_;
      Mtids[atid2] = nM_;
      nM_++;
    }

    MixingMap.resize(nM_);
    MixingMap[mtid1].resize(nM_);
    MixingMap[mtid1][mtid2] = mixer;
    if (mtid2 != mtid1) {
      MixingMap[mtid2].resize(nM_);
      MixingMap[mtid2][mtid1] = mixer;
    }
  }

  void Morse::calcForce(InteractionData& idat) {
    if (!initialized_) initialize();

    MorseInteractionData& mixer =
        MixingMap[Mtids[idat.atid1]][Mtids[idat.atid2]];

    RealType myPot    = 0.0;
    RealType myPotC   = 0.0;
    RealType myDeriv  = 0.0;
    RealType myDerivC = 0.0;

    RealType De       = mixer.De;
    RealType Re       = mixer.Re;
    RealType beta     = mixer.beta;
    MorseType variant = mixer.variant;

    // V(r) = D_e exp(-a(r-re)(exp(-a(r-re))-2)

    RealType expt    = -beta * (idat.rij - Re);
    RealType expfnc  = exp(expt);
    RealType expfnc2 = expfnc * expfnc;

    RealType exptC    = 0.0;
    RealType expfncC  = 0.0;
    RealType expfnc2C = 0.0;

    if (idat.shiftedPot || idat.shiftedForce) {
      exptC    = -beta * (idat.rcut - Re);
      expfncC  = exp(exptC);
      expfnc2C = expfncC * expfncC;
    }

    switch (variant) {
    case mtShifted: {
      myPot   = De * (expfnc2 - 2.0 * expfnc);
      myDeriv = 2.0 * De * beta * (expfnc - expfnc2);

      if (idat.shiftedPot) {
        myPotC   = De * (expfnc2C - 2.0 * expfncC);
        myDerivC = 0.0;
      } else if (idat.shiftedForce) {
        myPotC   = De * (expfnc2C - 2.0 * expfncC);
        myDerivC = 2.0 * De * beta * (expfncC - expfnc2C);
        myPotC += myDerivC * (idat.rij - idat.rcut);
      } else {
        myPotC   = 0.0;
        myDerivC = 0.0;
      }

      break;
    }
    case mtRepulsive: {
      myPot   = De * expfnc2;
      myDeriv = -2.0 * De * beta * expfnc2;

      if (idat.shiftedPot) {
        myPotC   = De * expfnc2C;
        myDerivC = 0.0;
      } else if (idat.shiftedForce) {
        myPotC   = De * expfnc2C;
        myDerivC = -2.0 * De * beta * expfnc2C;
        myPotC += myDerivC * (idat.rij - idat.rcut);
      } else {
        myPotC   = 0.0;
        myDerivC = 0.0;
      }

      break;
    }
    case mtUnknown: {
      // don't know what to do so don't do anything
      break;
    }
    }

    RealType pot_temp = idat.vdwMult * (myPot - myPotC);
    idat.vpair += pot_temp;

    RealType dudr = idat.sw * idat.vdwMult * (myDeriv - myDerivC);

    idat.pot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;
    if (idat.isSelected) idat.selePot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;

    idat.f1 += idat.d * dudr / idat.rij;

    return;
  }

  RealType Morse::getSuggestedCutoffRadius(pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();

    int atid1 = atypes.first->getIdent();
    int atid2 = atypes.second->getIdent();
    int mtid1 = Mtids[atid1];
    int mtid2 = Mtids[atid2];

    if (mtid1 == -1 || mtid2 == -1)
      return 0.0;
    else {
      MorseInteractionData mixer = MixingMap[mtid1][mtid2];
      RealType Re                = mixer.Re;
      RealType beta              = mixer.beta;
      // This value of the r corresponds to an energy about 1.48% of
      // the energy at the bottom of the Morse well.  For comparison, the
      // Lennard-Jones function is about 1.63% of it's minimum value at
      // a distance of 2.5 sigma.
      return (4.9 + beta * Re) / beta;
    }
  }
}  // namespace OpenMD
