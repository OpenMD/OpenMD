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

#include "types/AtomType.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "types/MultipoleAdapter.hpp"
#include "types/StickyAdapter.hpp"

using namespace std;

namespace OpenMD {
  AtomType::AtomType() {
    // initially, all atom types are their own base types.
    base_    = this;
    hasBase_ = false;

    // initialize to an error:
    ident_ = -1;

    // and massless:
    mass_                       = 0.0;
    myResponsibilities_["mass"] = false;
  }

  void AtomType::useBase(AtomType* base) {
    hasBase_ = true;
    base_    = base;
    base->addZig(this);
  }

  void AtomType::copyAllData(AtomType* orig) {
    // makes an exact replica of another atom type down to the
    // atom ID and any base attachments it may have.
    // use with caution!

    hasBase_ = orig->hasBase_;
    base_    = orig->base_;
    mass_    = orig->mass_;
    name_    = string(orig->name_);
    ident_   = orig->ident_;

    map<string, bool>::iterator i;
    ;
    map<string, RealType>::iterator j;

    for (i = orig->myResponsibilities_.begin();
         i != orig->myResponsibilities_.end(); ++i) {
      myResponsibilities_[(*i).first] = orig->myResponsibilities_[(*i).first];
    }

    for (j = orig->myValues_.begin(); j != orig->myValues_.end(); ++j) {
      myValues_[(*j).first] = orig->myValues_[(*j).first];
    }

    std::vector<std::shared_ptr<GenericData>> oprops = orig->getProperties();
    std::vector<std::shared_ptr<GenericData>>::iterator it;

    for (it = oprops.begin(); it != oprops.end(); ++it) {
      addProperty(*it);
    }
  }

  void AtomType::addProperty(std::shared_ptr<GenericData> genData) {
    myResponsibilities_[genData->getID()] = true;
    properties_.addProperty(genData);
  }

  void AtomType::removeProperty(const string& propName) {
    properties_.removeProperty(propName);
    myResponsibilities_[propName] = false;
  }

  std::vector<string> AtomType::getPropertyNames() {
    return properties_.getPropertyNames();
  }

  std::vector<std::shared_ptr<GenericData>> AtomType::getProperties() {
    return properties_.getProperties();
  }

  bool AtomType::hasProperty(const string& propName) {
    if (hasBase_ && !myResponsibilities_[propName]) {
      return base_->hasProperty(propName);
    } else
      return properties_.hasProperty(propName);
  }

  std::shared_ptr<GenericData> AtomType::getPropertyByName(
      const string& propName) {
    if (hasBase_ && !myResponsibilities_[propName]) {
      return base_->getPropertyByName(propName);
    } else
      return properties_.getPropertyByName(propName);
  }

  void AtomType::setMass(RealType m) {
    myResponsibilities_["mass"] = true;
    mass_                       = m;
  }

  RealType AtomType::getMass(void) {
    if (hasBase_ && !myResponsibilities_["mass"])
      return base_->getMass();
    else
      return mass_;
  }

  void AtomType::setIdent(int id) { ident_ = id; }

  int AtomType::getIdent() { return ident_; }

  void AtomType::setName(const string& name) { name_ = name; }

  string AtomType::getName() { return name_; }

  bool AtomType::isLennardJones() { return hasProperty("LJ"); }

  bool AtomType::isElectrostatic() { return isCharge() || isMultipole(); }

  bool AtomType::isEAM() { return hasProperty("EAM"); }

  bool AtomType::isCharge() { return isFixedCharge() || isFluctuatingCharge(); }

  bool AtomType::isDirectional() { return hasProperty("Directional"); }

  bool AtomType::isFluctuatingCharge() { return hasProperty("FlucQ"); }

  bool AtomType::isFixedCharge() { return hasProperty("Charge"); }

  bool AtomType::isDipole() {
    MultipoleAdapter ma = MultipoleAdapter(this);
    if (ma.isMultipole()) {
      return ma.isDipole();
    } else
      return false;
  }

  bool AtomType::isQuadrupole() {
    MultipoleAdapter ma = MultipoleAdapter(this);
    if (ma.isMultipole()) {
      return ma.isQuadrupole();
    } else
      return false;
  }

  bool AtomType::isMultipole() { return hasProperty("Multipole"); }

  bool AtomType::isGayBerne() { return hasProperty("GB"); }

  bool AtomType::isSticky() { return hasProperty("Sticky"); }

  bool AtomType::isStickyPower() {
    StickyAdapter sa = StickyAdapter(this);
    return sa.isStickyPower();
  }

  bool AtomType::isShape() { return hasProperty("Shape"); }

  bool AtomType::isSC() { return hasProperty("SC"); }

  bool AtomType::isMetal() { return isSC() || isEAM(); }

  std::vector<AtomType*> AtomType::allYourBase() {
    std::vector<AtomType*> myChain;

    if (hasBase_) {
      myChain = base_->allYourBase();
      myChain.insert(myChain.begin(), this);
    } else {
      myChain.push_back(this);
    }

    return myChain;
  }

}  // namespace OpenMD
