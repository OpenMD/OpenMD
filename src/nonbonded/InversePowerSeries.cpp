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

#include "nonbonded/InversePowerSeries.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>

#include "types/InversePowerSeriesInteractionType.hpp"
#include "utils/simError.h"

using namespace std;

namespace OpenMD {

  InversePowerSeries::InversePowerSeries() :
      initialized_(false), forceField_(NULL), name_("InversePowerSeries") {}

  void InversePowerSeries::initialize() {
    InversePowerSeriesTypes.clear();
    InversePowerSeriesTids.clear();
    MixingMap.clear();
    InversePowerSeriesTids.resize(forceField_->getNAtomType(), -1);

    ForceField::NonBondedInteractionTypeContainer* nbiTypes =
        forceField_->getNonBondedInteractionTypes();
    ForceField::NonBondedInteractionTypeContainer::MapTypeIterator j;
    ForceField::NonBondedInteractionTypeContainer::KeyType keys;
    NonBondedInteractionType* nbt;
    int ipstid1, ipstid2;

    for (nbt = nbiTypes->beginType(j); nbt != NULL;
         nbt = nbiTypes->nextType(j)) {
      if (nbt->isInversePowerSeries()) {
        keys          = nbiTypes->getKeys(j);
        AtomType* at1 = forceField_->getAtomType(keys[0]);
        if (at1 == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "InversePowerSeries::initialize could not find AtomType %s\n"
                   "\tto for for %s - %s interaction.\n",
                   keys[0].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        AtomType* at2 = forceField_->getAtomType(keys[1]);
        if (at2 == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "InversePowerSeries::initialize could not find AtomType %s\n"
                   "\tfor %s - %s nonbonded interaction.\n",
                   keys[1].c_str(), keys[0].c_str(), keys[1].c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
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

        InversePowerSeriesInteractionType* ipsit =
            dynamic_cast<InversePowerSeriesInteractionType*>(nbt);
        if (ipsit == NULL) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "InversePowerSeries::initialize could not convert "
                   "NonBondedInteractionType\n"
                   "\tto InversePowerSeriesInteractionType for %s - %s "
                   "interaction.\n",
                   at1->getName().c_str(), at2->getName().c_str());
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
        }

        std::vector<int> powers            = ipsit->getPowers();
        std::vector<RealType> coefficients = ipsit->getCoefficients();

        addExplicitInteraction(at1, at2, powers, coefficients);
      }
    }
    initialized_ = true;
  }

  void InversePowerSeries::addExplicitInteraction(
      AtomType* atype1, AtomType* atype2, std::vector<int> powers,
      std::vector<RealType> coefficients) {
    InversePowerSeriesInteractionData mixer;
    mixer.powers       = powers;
    mixer.coefficients = coefficients;

    int ipstid1             = InversePowerSeriesTids[atype1->getIdent()];
    int ipstid2             = InversePowerSeriesTids[atype2->getIdent()];
    int nInversePowerSeries = InversePowerSeriesTypes.size();

    MixingMap.resize(nInversePowerSeries);
    MixingMap[ipstid1].resize(nInversePowerSeries);

    MixingMap[ipstid1][ipstid2] = mixer;
    if (ipstid2 != ipstid1) {
      MixingMap[ipstid2].resize(nInversePowerSeries);
      MixingMap[ipstid2][ipstid1] = mixer;
    }
  }

  void InversePowerSeries::calcForce(InteractionData& idat) {
    if (!initialized_) initialize();

    InversePowerSeriesInteractionData& mixer =
        MixingMap[InversePowerSeriesTids[idat.atid1]]
                 [InversePowerSeriesTids[idat.atid2]];
    std::vector<int> powers            = mixer.powers;
    std::vector<RealType> coefficients = mixer.coefficients;

    RealType myPot    = 0.0;
    RealType myPotC   = 0.0;
    RealType myDeriv  = 0.0;
    RealType myDerivC = 0.0;

    RealType ri  = 1.0 / idat.rij;
    RealType ric = 1.0 / idat.rcut;

    RealType fn, fnc;

    for (unsigned int i = 0; i < powers.size(); i++) {
      fn  = coefficients[i] * pow(ri, powers[i]);
      fnc = coefficients[i] * pow(ric, powers[i]);
      myPot += fn;
      myPotC += fnc;
      myDeriv -= powers[i] * fn * ri;
      myDerivC -= powers[i] * fnc * ric;
    }

    if (idat.shiftedPot) {
      myDerivC = 0.0;
    } else if (idat.shiftedForce) {
      myPotC = myPotC + myDerivC * (idat.rij - idat.rcut);
    } else {
      myPotC   = 0.0;
      myDerivC = 0.0;
    }

    RealType pot_temp = idat.vdwMult * (myPot - myPotC);
    idat.vpair += pot_temp;

    RealType dudr = idat.sw * idat.vdwMult * (myDeriv - myDerivC);

    idat.pot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;
    if (idat.isSelected) idat.selePot[VANDERWAALS_FAMILY] += idat.sw * pot_temp;

    idat.f1 += idat.d * dudr / idat.rij;

    return;
  }

  RealType InversePowerSeries::getSuggestedCutoffRadius(
      pair<AtomType*, AtomType*> atypes) {
    if (!initialized_) initialize();

    int atid1   = atypes.first->getIdent();
    int atid2   = atypes.second->getIdent();
    int ipstid1 = InversePowerSeriesTids[atid1];
    int ipstid2 = InversePowerSeriesTids[atid2];

    if (ipstid1 == -1 || ipstid2 == -1)
      return 0.0;
    else {
      // This seems to work moderately well as a default.  There's no
      // inherent scale for 1/r^n interactions that we can standardize.
      // 12 angstroms seems to be a reasonably good guess for most
      // cases.
      return 12.0;
    }
  }
}  // namespace OpenMD
