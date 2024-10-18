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

#include "nonbonded/InteractionManager.hpp"

#include "types/InversePowerSeriesInteractionType.hpp"
#include "types/LennardJonesInteractionType.hpp"
#include "types/MAWInteractionType.hpp"
#include "types/MieInteractionType.hpp"
#include "types/MorseInteractionType.hpp"
#include "types/RepulsivePowerInteractionType.hpp"

namespace OpenMD {

  InteractionManager::InteractionManager() {
    initialized_ = false;

    lj_                 = std::make_shared<LJ>();
    gb_                 = std::make_shared<GB>();
    sticky_             = std::make_shared<Sticky>();
    morse_              = std::make_shared<Morse>();
    repulsivePower_     = std::make_shared<RepulsivePower>();
    mie_                = std::make_shared<Mie>();
    eam_                = std::make_shared<EAM>();
    sc_                 = std::make_shared<SC>();
    electrostatic_      = std::make_shared<Electrostatic>();
    maw_                = std::make_shared<MAW>();
    inversePowerSeries_ = std::make_shared<InversePowerSeries>();
  }

  void InteractionManager::initialize() {
    if (initialized_) return;

    ForceField* forceField_ = info_->getForceField();

    lj_->setForceField(forceField_);
    gb_->setForceField(forceField_);
    sticky_->setForceField(forceField_);
    eam_->setForceField(forceField_);
    eam_->setElectrostatic(electrostatic_.get());
    sc_->setForceField(forceField_);
    morse_->setForceField(forceField_);
    electrostatic_->setSimInfo(info_);
    electrostatic_->setForceField(forceField_);
    maw_->setForceField(forceField_);
    repulsivePower_->setForceField(forceField_);
    mie_->setForceField(forceField_);
    inversePowerSeries_->setForceField(forceField_);

    ForceField::AtomTypeContainer* atomTypes = forceField_->getAtomTypes();
    int nTypes                               = atomTypes->size();
    sHash_.resize(nTypes);
    iHash_.resize(nTypes);
    interactions_.resize(nTypes);
    ForceField::AtomTypeContainer::MapTypeIterator i1, i2;
    AtomType* atype1;
    AtomType* atype2;
    int atid1, atid2;

    // We only need to worry about the types that are actually in the
    // simulation:

    AtomTypeSet atypes = info_->getSimulatedAtomTypes();

    lj_->setSimulatedAtomTypes(atypes);
    gb_->setSimulatedAtomTypes(atypes);
    sticky_->setSimulatedAtomTypes(atypes);
    eam_->setSimulatedAtomTypes(atypes);
    sc_->setSimulatedAtomTypes(atypes);
    morse_->setSimulatedAtomTypes(atypes);
    electrostatic_->setSimInfo(info_);
    electrostatic_->setSimulatedAtomTypes(atypes);
    maw_->setSimulatedAtomTypes(atypes);
    repulsivePower_->setSimulatedAtomTypes(atypes);
    mie_->setSimulatedAtomTypes(atypes);
    inversePowerSeries_->setSimulatedAtomTypes(atypes);

    AtomTypeSet::iterator at;
    set<NonBondedInteractionPtr>::iterator it;

    for (at = atypes.begin(); at != atypes.end(); ++at) {
      atype1 = *at;
      atid1  = atype1->getIdent();
      iHash_[atid1].resize(nTypes);
      interactions_[atid1].resize(nTypes);

      // add it to the map:
      pair<map<int, AtomType*>::iterator, bool> ret;
      ret = typeMap_.insert(pair<int, AtomType*>(atid1, atype1));
      if (ret.second == false) {
        snprintf(
            painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
            "InteractionManager already had a previous entry with ident %d\n",
            atype1->getIdent());
        painCave.severity = OPENMD_INFO;
        painCave.isFatal  = 0;
        simError();
      }

      if (atype1->isLennardJones()) { sHash_[atid1] |= LJ_INTERACTION; }
      if (atype1->isElectrostatic()) {
        sHash_[atid1] |= ELECTROSTATIC_INTERACTION;
      }
      if (atype1->isSticky()) { sHash_[atid1] |= STICKY_INTERACTION; }
      if (atype1->isStickyPower()) { sHash_[atid1] |= STICKY_INTERACTION; }
      if (atype1->isEAM()) { sHash_[atid1] |= EAM_INTERACTION; }
      if (atype1->isSC()) { sHash_[atid1] |= SC_INTERACTION; }
      if (atype1->isGayBerne()) { sHash_[atid1] |= GB_INTERACTION; }
    }
    // Now, iterate over all known types and add to the interaction map:

    map<int, AtomType*>::iterator it1, it2;
    for (it1 = typeMap_.begin(); it1 != typeMap_.end(); ++it1) {
      atype1 = (*it1).second;
      atid1  = atype1->getIdent();

      for (it2 = typeMap_.begin(); it2 != typeMap_.end(); ++it2) {
        atype2 = (*it2).second;
        atid2  = atype2->getIdent();

        iHash_[atid1][atid2] = 0;

        if (atype1->isLennardJones() && atype2->isLennardJones()) {
          interactions_[atid1][atid2].insert(lj_);
          iHash_[atid1][atid2] |= LJ_INTERACTION;
        }
        if (atype1->isElectrostatic() && atype2->isElectrostatic()) {
          interactions_[atid1][atid2].insert(electrostatic_);
          iHash_[atid1][atid2] |= ELECTROSTATIC_INTERACTION;
        }

        // A special case for calculating local fields:
        if (info_->getSimParams()->getOutputElectricField()) {
          if (atype1->isElectrostatic() || atype2->isElectrostatic()) {
            interactions_[atid1][atid2].insert(electrostatic_);
            iHash_[atid1][atid2] |= ELECTROSTATIC_INTERACTION;
          }
        }

        if (atype1->isSticky() && atype2->isSticky()) {
          interactions_[atid1][atid2].insert(sticky_);
          iHash_[atid1][atid2] |= STICKY_INTERACTION;
        }
        if (atype1->isStickyPower() && atype2->isStickyPower()) {
          interactions_[atid1][atid2].insert(sticky_);
          iHash_[atid1][atid2] |= STICKY_INTERACTION;
        }
        if (atype1->isEAM() && atype2->isEAM()) {
          interactions_[atid1][atid2].insert(eam_);
          iHash_[atid1][atid2] |= EAM_INTERACTION;
        }
        if (atype1->isSC() && atype2->isSC()) {
          interactions_[atid1][atid2].insert(sc_);
          iHash_[atid1][atid2] |= SC_INTERACTION;
        }
        if (atype1->isGayBerne() && atype2->isGayBerne()) {
          interactions_[atid1][atid2].insert(gb_);
          iHash_[atid1][atid2] |= GB_INTERACTION;
        }
        if ((atype1->isGayBerne() && atype2->isLennardJones()) ||
            (atype1->isLennardJones() && atype2->isGayBerne())) {
          interactions_[atid1][atid2].insert(gb_);
          iHash_[atid1][atid2] |= GB_INTERACTION;
        }

        // look for an explicitly-set non-bonded interaction type using the
        // two atom types.
        NonBondedInteractionType* nbiType =
            forceField_->getNonBondedInteractionType(atype1->getName(),
                                                     atype2->getName());

        if (nbiType != NULL) {
          bool vdwExplicit = false;
          bool metExplicit = false;
          // bool hbExplicit = false;

          if (nbiType->isLennardJones()) {
            // We found an explicit Lennard-Jones interaction.
            // override all other vdw entries for this pair of atom types:
            for (it = interactions_[atid1][atid2].begin();
                 it != interactions_[atid1][atid2].end();) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == VANDERWAALS_FAMILY) {
                iHash_[atid1][atid2] ^= (*it)->getHash();
                interactions_[atid1][atid2].erase(it++);
              } else {
                ++it;
              }
            }
            interactions_[atid1][atid2].insert(lj_);
            iHash_[atid1][atid2] |= LJ_INTERACTION;
            LennardJonesInteractionType* ljit =
                dynamic_cast<LennardJonesInteractionType*>(nbiType);
            lj_->addExplicitInteraction(atype1, atype2, ljit->getSigma(),
                                        ljit->getEpsilon());
            vdwExplicit = true;
          }

          if (nbiType->isMorse()) {
            if (vdwExplicit) {
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "InteractionManager::initialize found more than one "
                       "explicit \n"
                       "\tvan der Waals interaction for atom types %s - %s\n",
                       atype1->getName().c_str(), atype2->getName().c_str());
              painCave.severity = OPENMD_ERROR;
              painCave.isFatal  = 1;
              simError();
            }
            // We found an explicit Morse interaction.
            // override all other vdw entries for this pair of atom types:
            for (it = interactions_[atid1][atid2].begin();
                 it != interactions_[atid1][atid2].end();) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == VANDERWAALS_FAMILY) {
                iHash_[atid1][atid2] ^= (*it)->getHash();
                interactions_[atid1][atid2].erase(it++);
              } else {
                ++it;
              }
            }
            interactions_[atid1][atid2].insert(morse_);
            iHash_[atid1][atid2] |= MORSE_INTERACTION;
            MorseInteractionType* mit =
                dynamic_cast<MorseInteractionType*>(nbiType);
            morse_->addExplicitInteraction(atype1, atype2, mit->getD(),
                                           mit->getR(), mit->getBeta(),
                                           mit->getInteractionType());
            vdwExplicit = true;
          }

          if (nbiType->isRepulsivePower()) {
            if (vdwExplicit) {
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "InteractionManager::initialize found more than one "
                       "explicit \n"
                       "\tvan der Waals interaction for atom types %s - %s\n",
                       atype1->getName().c_str(), atype2->getName().c_str());
              painCave.severity = OPENMD_ERROR;
              painCave.isFatal  = 1;
              simError();
            }
            // We found an explicit RepulsivePower interaction.
            // override all other vdw entries for this pair of atom types:
            for (it = interactions_[atid1][atid2].begin();
                 it != interactions_[atid1][atid2].end();) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == VANDERWAALS_FAMILY) {
                iHash_[atid1][atid2] ^= (*it)->getHash();
                interactions_[atid1][atid2].erase(it++);
              } else {
                ++it;
              }
            }
            interactions_[atid1][atid2].insert(repulsivePower_);
            iHash_[atid1][atid2] |= REPULSIVEPOWER_INTERACTION;
            RepulsivePowerInteractionType* rpit =
                dynamic_cast<RepulsivePowerInteractionType*>(nbiType);

            repulsivePower_->addExplicitInteraction(
                atype1, atype2, rpit->getSigma(), rpit->getEpsilon(),
                rpit->getNrep());

            vdwExplicit = true;
          }

          if (nbiType->isMie()) {
            if (vdwExplicit) {
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "InteractionManager::initialize found more than one "
                       "explicit \n"
                       "\tvan der Waals interaction for atom types %s - %s\n",
                       atype1->getName().c_str(), atype2->getName().c_str());
              painCave.severity = OPENMD_ERROR;
              painCave.isFatal  = 1;
              simError();
            }
            // We found an explicit Mie interaction.
            // override all other vdw entries for this pair of atom types:
            for (it = interactions_[atid1][atid2].begin();
                 it != interactions_[atid1][atid2].end();) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == VANDERWAALS_FAMILY) {
                iHash_[atid1][atid2] ^= (*it)->getHash();
                interactions_[atid1][atid2].erase(it++);
              } else {
                ++it;
              }
            }
            interactions_[atid1][atid2].insert(mie_);
            iHash_[atid1][atid2] |= MIE_INTERACTION;
            MieInteractionType* mit =
                dynamic_cast<MieInteractionType*>(nbiType);

            mie_->addExplicitInteraction(atype1, atype2, mit->getSigma(),
                                         mit->getEpsilon(), mit->getNrep(),
                                         mit->getMatt());

            vdwExplicit = true;
          }

          if (nbiType->isEAMTable() || nbiType->isEAMZhou()) {
            // We found an explicit EAM interaction.
            // override all other metallic entries for this pair of atom types:
            for (it = interactions_[atid1][atid2].begin();
                 it != interactions_[atid1][atid2].end();) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == METALLIC_EMBEDDING_FAMILY) {
                iHash_[atid1][atid2] ^= (*it)->getHash();
                interactions_[atid1][atid2].erase(it++);
              } else {
                ++it;
              }
            }
            interactions_[atid1][atid2].insert(eam_);
            iHash_[atid1][atid2] |= EAM_INTERACTION;
            metExplicit = true;
          }

          if (nbiType->isSC()) {
            if (metExplicit) {
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "InteractionManager::initialize found more than one "
                       "explicit\n"
                       "\tmetallic interaction for atom types %s - %s\n",
                       atype1->getName().c_str(), atype2->getName().c_str());
              painCave.severity = OPENMD_ERROR;
              painCave.isFatal  = 1;
              simError();
            }
            // We found an explicit Sutton-Chen interaction.
            // override all other metallic entries for this pair of atom types:
            for (it = interactions_[atid1][atid2].begin();
                 it != interactions_[atid1][atid2].end();) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == METALLIC_EMBEDDING_FAMILY) {
                iHash_[atid1][atid2] ^= (*it)->getHash();
                interactions_[atid1][atid2].erase(it++);
              } else {
                ++it;
              }
            }
            interactions_[atid1][atid2].insert(sc_);
            iHash_[atid1][atid2] |= SC_INTERACTION;
            metExplicit = true;
          }

          if (nbiType->isMAW()) {
            if (vdwExplicit) {
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "InteractionManager::initialize found more than one "
                       "explicit\n"
                       "\tvan der Waals interaction for atom types %s - %s\n",
                       atype1->getName().c_str(), atype2->getName().c_str());
              painCave.severity = OPENMD_ERROR;
              painCave.isFatal  = 1;
              simError();
            }
            // We found an explicit MAW interaction.
            // override all other vdw entries for this pair of atom types:
            for (it = interactions_[atid1][atid2].begin();
                 it != interactions_[atid1][atid2].end();) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == VANDERWAALS_FAMILY) {
                iHash_[atid1][atid2] ^= (*it)->getHash();
                interactions_[atid1][atid2].erase(it++);
              } else {
                ++it;
              }
            }
            interactions_[atid1][atid2].insert(maw_);
            iHash_[atid1][atid2] |= MAW_INTERACTION;
            MAWInteractionType* mit =
                dynamic_cast<MAWInteractionType*>(nbiType);
            maw_->addExplicitInteraction(atype1, atype2, mit->getD(),
                                         mit->getBeta(), mit->getR(),
                                         mit->getCA1(), mit->getCB1());
            vdwExplicit = true;
          }

          if (nbiType->isInversePowerSeries()) {
            if (vdwExplicit) {
              snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                       "InteractionManager::initialize found more than one "
                       "explicit \n"
                       "\tvan der Waals interaction for atom types %s - %s\n",
                       atype1->getName().c_str(), atype2->getName().c_str());
              painCave.severity = OPENMD_ERROR;
              painCave.isFatal  = 1;
              simError();
            }
            // We found an explicit InversePowerSeries interaction.
            // override all other vdw entries for this pair of atom types:
            for (it = interactions_[atid1][atid2].begin();
                 it != interactions_[atid1][atid2].end();) {
              InteractionFamily ifam = (*it)->getFamily();
              if (ifam == VANDERWAALS_FAMILY) {
                iHash_[atid1][atid2] ^= (*it)->getHash();
                interactions_[atid1][atid2].erase(it++);
              } else {
                ++it;
              }
            }
            interactions_[atid1][atid2].insert(inversePowerSeries_);
            iHash_[atid1][atid2] |= INVERSEPOWERSERIES_INTERACTION;
            InversePowerSeriesInteractionType* ipsit =
                dynamic_cast<InversePowerSeriesInteractionType*>(nbiType);

            inversePowerSeries_->addExplicitInteraction(
                atype1, atype2, ipsit->getPowers(), ipsit->getCoefficients());
            vdwExplicit = true;
          }
        }
      }
    }

    // Make sure every pair of atom types in this simulation has a
    // non-bonded interaction.  If not, just inform the user.

    AtomTypeSet simTypes = info_->getSimulatedAtomTypes();
    AtomTypeSet::iterator bt;

    for (at = simTypes.begin(); at != simTypes.end(); ++at) {
      atype1 = (*at);
      atid1  = atype1->getIdent();
      for (bt = at; bt != simTypes.end(); ++bt) {
        atype2 = (*bt);
        atid2  = atype2->getIdent();

        if (interactions_[atid1][atid2].size() == 0) {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "InteractionManager could not find a matching non-bonded\n"
                   "\tinteraction for atom types %s - %s\n"
                   "\tProceeding without this interaction.\n",
                   atype1->getName().c_str(), atype2->getName().c_str());
          painCave.severity = OPENMD_INFO;
          painCave.isFatal  = 0;
          simError();
        }
      }
    }

    initialized_ = true;
  }

  void InteractionManager::setCutoffRadius(RealType rcut) {
    electrostatic_->setCutoffRadius(rcut);
    eam_->setCutoffRadius(rcut);
  }

  void InteractionManager::doPrePair(InteractionData& idat) {
    if (!initialized_) initialize();

    // excluded interaction, so just return
    if (idat.excluded) return;

    int& iHash = iHash_[idat.atid1][idat.atid2];

    if ((iHash & EAM_INTERACTION) != 0) eam_->calcDensity(idat);
    if ((iHash & SC_INTERACTION) != 0) sc_->calcDensity(idat);

    // set<NonBondedInteraction*>::iterator it;
    //
    // for (it = interactions_[ idat.atypes ].begin();
    //      it != interactions_[ idat.atypes ].end(); ++it){
    //   if ((*it)->getFamily() == METALLIC_EMBEDDING_FAMILY) {
    //     dynamic_cast<MetallicInteraction*>(*it)->calcDensity(idat);
    //   }
    // }

    return;
  }

  void InteractionManager::doPreForce(SelfData& sdat) {
    if (!initialized_) initialize();

    int& sHash = sHash_[sdat.atid];

    if ((sHash & EAM_INTERACTION) != 0) eam_->calcFunctional(sdat);
    if ((sHash & SC_INTERACTION) != 0) sc_->calcFunctional(sdat);

    // set<NonBondedInteraction*>::iterator it;
    //
    // for (it = interactions_[atid1][atid2].begin();
    //      it != interactions_[atid1][atid2].end(); ++it){
    //   if ((*it)->getFamily() == METALLIC_EMBEDDING_FAMILY) {
    //     dynamic_cast<MetallicInteraction*>(*it)->calcFunctional(sdat);
    //   }
    // }

    return;
  }

  void InteractionManager::doPair(InteractionData& idat) {
    if (!initialized_) initialize();

    int& iHash = iHash_[idat.atid1][idat.atid2];

    if ((iHash & ELECTROSTATIC_INTERACTION) != 0)
      electrostatic_->calcForce(idat);

    // electrostatics still has to worry about indirect
    // contributions from excluded pairs of atoms, but nothing else does:

    if (idat.excluded) return;

    if ((iHash & LJ_INTERACTION) != 0) lj_->calcForce(idat);
    if ((iHash & GB_INTERACTION) != 0) gb_->calcForce(idat);
    if ((iHash & STICKY_INTERACTION) != 0) sticky_->calcForce(idat);
    if ((iHash & MORSE_INTERACTION) != 0) morse_->calcForce(idat);
    if ((iHash & REPULSIVEPOWER_INTERACTION) != 0)
      repulsivePower_->calcForce(idat);
    if ((iHash & MIE_INTERACTION) != 0) mie_->calcForce(idat);
    if ((iHash & EAM_INTERACTION) != 0) eam_->calcForce(idat);
    if ((iHash & SC_INTERACTION) != 0) sc_->calcForce(idat);
    if ((iHash & MAW_INTERACTION) != 0) maw_->calcForce(idat);
    if ((iHash & INVERSEPOWERSERIES_INTERACTION) != 0)
      inversePowerSeries_->calcForce(idat);

    // set<NonBondedInteraction*>::iterator it;
    //
    // for (it = interactions_[ idat.atypes ].begin();
    //      it != interactions_[ idat.atypes ].end(); ++it) {
    //
    //   if (!idat.excluded || (*it)->getFamily() == ELECTROSTATIC_FAMILY) {
    //     (*it)->calcForce(idat);
    //   }
    // }

    return;
  }

  void InteractionManager::doSelfCorrection(SelfData& sdat) {
    if (!initialized_) initialize();

    int& sHash = sHash_[sdat.atid];

    if ((sHash & ELECTROSTATIC_INTERACTION) != 0) {
      electrostatic_->calcSelfCorrection(sdat);
    }

    // set<NonBondedInteraction*>::iterator it;
    //
    // for (it = interactions_[atid1][atid2].begin();
    //      it != interactions_[atid1][atid2].end(); ++it){
    //   if ((*it)->getFamily() == ELECTROSTATIC_FAMILY) {
    //     dynamic_cast<ElectrostaticInteraction*>(*it)->calcSelfCorrection(sdat);
    //   }
    // }

    return;
  }

  void InteractionManager::doSurfaceTerm(bool slabGeometry, int axis,
                                         RealType& pot) {
    if (!initialized_) initialize();
    electrostatic_->calcSurfaceTerm(slabGeometry, axis, pot);
  }

  void InteractionManager::doReciprocalSpaceSum(RealType& pot) {
    if (!initialized_) initialize();
    electrostatic_->ReciprocalSpaceSum(pot);
  }

  RealType InteractionManager::getSuggestedCutoffRadius(int* atid) {
    if (!initialized_) initialize();

    AtomType* atype = typeMap_[*atid];

    set<NonBondedInteractionPtr>::iterator it;
    RealType cutoff = 0.0;

    for (it = interactions_[*atid][*atid].begin();
         it != interactions_[*atid][*atid].end(); ++it) {
      cutoff =
          max(cutoff, (*it)->getSuggestedCutoffRadius(make_pair(atype, atype)));
    }
    return cutoff;
  }

  RealType InteractionManager::getSuggestedCutoffRadius(AtomType* atype) {
    if (!initialized_) initialize();

    int atid = atype->getIdent();

    set<NonBondedInteractionPtr>::iterator it;
    RealType cutoff = 0.0;

    for (it = interactions_[atid][atid].begin();
         it != interactions_[atid][atid].end(); ++it) {
      cutoff =
          max(cutoff, (*it)->getSuggestedCutoffRadius(make_pair(atype, atype)));
    }
    return cutoff;
  }
}  // namespace OpenMD
