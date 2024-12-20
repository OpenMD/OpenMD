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

#include "nonbonded/Buckingham.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>

#include "types/BuckinghamInteractionType.hpp"
#include "utils/simError.h"

using namespace std;

namespace OpenMD {

  Buckingham::Buckingham() :
      initialized_(false), forceField_(NULL), name_("Buckingham") {}

  void Buckingham::initialize() {
    Btypes.clear();
    Btids.clear();
    MixingMap.clear();
    Btids.resize(forceField_->getNAtomType(), -1);

    ForceField::NonBondedInteractionTypeContainer* nbiTypes =
        forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    ForceField::NonBondedInteractionTypeContainer::KeyType keys;
    NonBondedInteractionType* nbt;
    int btid1, btid2;

    for (nbt = nbiTypes->beginType(j); nbt != NULL;
         nbt = nbiTypes->nextType(j)) {
      if (nbt->isBuckingham()) {
        keys          = nbiTypes->getKeys(j);
        AtomType* at1 = forceField_->getAtomType(keys[0]);
        if (at1 == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Buckingham::initialize could not find AtomType %s\n"
                   "\tto for for %s - %s interaction.\n",
                   keys[0].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        AtomType* at2 = forceField_->getAtomType(keys[1]);
        if (at2 == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Buckingham::initialize could not find AtomType %s\n"
                   "\tfor %s - %s nonbonded interaction.\n",
                   keys[1].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
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

        BuckinghamInteractionType* bit =
            dynamic_cast<BuckinghamInteractionType*>(nbt);

        if (bit == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "Buckingham::initialize could not convert "
                   "NonBondedInteractionType\n"
                   "\tto BuckinghamInteractionType for %s - %s interaction.\n",
                   at1->getName().c_str(), at2->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        RealType A = bit->getA();
        RealType B = bit->getB();
        RealType C = bit->getC();

        BuckinghamType variant = bit->getInteractionType();
        addExplicitInteraction(at1, at2, A, B, C, variant);
      }
    }
    initialized_ = true;
  }

  void Buckingham::addExplicitInteraction(AtomType* atype1, AtomType* atype2,
                                          RealType A, RealType B, RealType C,
                                          BuckinghamType bt) {
    BuckinghamInteractionData mixer;
    mixer.A       = A;
    mixer.B       = B;
    mixer.C       = C;
    mixer.variant = bt;

    int btid1 = Btids[atype1->getIdent()];
    int btid2 = Btids[atype2->getIdent()];
    int nB    = Btypes.size();

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
    mixer.A       = A;
    mixer.B       = B;
    mixer.C       = C;
    mixer.sigma   = sigma;
    mixer.epsilon = epsilon;
    mixer.variant = bt;

    int btid1 = Btids[atype1->getIdent()];
    int btid2 = Btids[atype2->getIdent()];
    int nB    = Btypes.size();

    MixingMap.resize(nB);
    MixingMap[btid1].resize(nB);

    MixingMap[btid1][btid2] = mixer;
    if (btid2 != btid1) {
      MixingMap[btid2].resize(nB);
      MixingMap[btid2][btid1] = mixer;
    }
  }

  void Buckingham::calcForce(InteractionData& idat) {
    if (!initialized_) initialize();

    BuckinghamInteractionData& mixer =
        MixingMap[Btids[idat.atid1]][Btids[idat.atid2]];

    RealType myPot    = 0.0;
    RealType myPotC   = 0.0;
    RealType myDeriv  = 0.0;
    RealType myDerivC = 0.0;

    RealType A             = mixer.A;
    RealType B             = mixer.B;
    RealType C             = mixer.C;
    RealType sigma         = mixer.sigma;
    RealType epsilon       = mixer.epsilon;
    BuckinghamType variant = mixer.variant;

    RealType expt   = -B * idat.rij;
    RealType expfnc = exp(expt);
    RealType fnc6   = 1.0 / pow(idat.rij, 6);
    RealType fnc7   = fnc6 / idat.rij;

    RealType exptC   = 0.0;
    RealType expfncC = 0.0;
    RealType fnc6C   = 0.0;
    RealType fnc7C   = 0.0;

    if (idat.shiftedPot || idat.shiftedForce) {
      exptC   = -B * idat.rcut;
      expfncC = exp(exptC);
      fnc6C   = 1.0 / pow(idat.rcut, 6);
      fnc7C   = fnc6C / idat.rcut;
    }

    switch (variant) {
    case btTraditional: {
      // V(r) = A exp(-B*r) - C/r^6
      myPot   = A * expfnc - C * fnc6;
      myDeriv = -A * B * expfnc + C * fnc7;

      if (idat.shiftedPot) {
        myPotC   = A * expfncC - C * fnc6C;
        myDerivC = 0.0;
      } else if (idat.shiftedForce) {
        myPotC   = A * expfncC - C * fnc6C;
        myDerivC = -A * B * expfncC + C * fnc7C;
        myPotC += myDerivC * (idat.rij - idat.rcut);
      } else {
        myPotC   = 0.0;
        myDerivC = 0.0;
      }
      break;
    }
    case btModified: {
      RealType s6     = pow(sigma, 6);
      RealType s7     = pow(sigma, 7);
      RealType fnc30  = pow(sigma / idat.rij, 30);
      RealType fnc31  = fnc30 * sigma / idat.rij;
      RealType fnc30C = 0.0;
      RealType fnc31C = 0.0;

      if (idat.shiftedPot || idat.shiftedForce) {
        fnc30C = pow(sigma / idat.rcut, 30);
        fnc31C = fnc30C * sigma / idat.rcut;
      }

      // V(r) = A exp(-B*r) - C/r^6 + 4 epsilon ((sigma/r)^30 - (sigma/r)^6)
      myPot   = A * expfnc - C * fnc6 + 4.0 * epsilon * (fnc30 - s6 * fnc6);
      myDeriv = -A * B * expfnc + C * fnc7 +
                4.0 * epsilon * (-30.0 * fnc31 + 6.0 * s7 * fnc7) / sigma;

      if (idat.shiftedPot) {
        myPotC =
            A * expfncC - C * fnc6C + 4.0 * epsilon * (fnc30C - s6 * fnc6C);
        myDerivC = 0.0;
      } else if (idat.shiftedForce) {
        myPotC =
            A * expfncC - C * fnc6C + 4.0 * epsilon * (fnc30C - s6 * fnc6C);
        myDeriv = -A * B * expfncC + C * fnc7C +
                  4.0 * epsilon * (-30.0 * fnc31C + 6.0 * s7 * fnc7C) / sigma;
        myPotC += myDerivC * (idat.rij - idat.rcut);
      } else {
        myPotC   = 0.0;
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
    if (idat.isSelected) idat.selePot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;

    idat.f1 += idat.d * dudr / idat.rij;

    return;
  }

  RealType Buckingham::getSuggestedCutoffRadius(
      pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();

    int atid1 = atypes.first->getIdent();
    int atid2 = atypes.second->getIdent();
    int btid1 = Btids[atid1];
    int btid2 = Btids[atid2];

    if (btid1 == -1 || btid2 == -1)
      return 0.0;
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
}  // namespace OpenMD
