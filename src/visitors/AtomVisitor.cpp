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

#include "visitors/AtomVisitor.hpp"

#include <cstring>
#include <memory>

#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/GayBerneAdapter.hpp"
#include "types/MultipoleAdapter.hpp"

namespace OpenMD {

  BaseAtomVisitor::BaseAtomVisitor(SimInfo* info) : BaseVisitor() {
    storageLayout_ = info->getAtomStorageLayout();
  }

  void BaseAtomVisitor::visit(RigidBody*) {
    // vector<Atom*> myAtoms;
    // vector<Atom*>::iterator atomIter;

    // myAtoms = rb->getAtoms();

    // for(atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter)
    //  (*atomIter)->accept(this);
  }

  void BaseAtomVisitor::setVisited(Atom* atom) {
    std::shared_ptr<GenericData> data;
    data = atom->getPropertyByName("VISITED");

    // if visited property is not existed, add it as new property
    if (data == nullptr) {
      data = std::make_shared<GenericData>();
      data->setID("VISITED");
      atom->addProperty(data);
    }
  }

  bool BaseAtomVisitor::isVisited(Atom* atom) {
    std::shared_ptr<GenericData> data;
    data = atom->getPropertyByName("VISITED");
    return data == nullptr ? false : true;
  }

  //------------------------------------------------------------------------//

  void DefaultAtomVisitor::visit(Atom* atom) {
    std::shared_ptr<AtomData> atomData;
    std::shared_ptr<AtomInfo> atomInfo;
    AtomType* atype = atom->getAtomType();

    if (isVisited(atom)) return;

    atomInfo               = std::make_shared<AtomInfo>();
    atomInfo->atomTypeName = atom->getType();
    atomInfo->globalID     = atom->getGlobalIndex();
    atomInfo->pos          = atom->getPos();
    atomInfo->vel          = atom->getVel();
    atomInfo->frc          = atom->getFrc();
    atomInfo->vec          = V3Zero;
    atomInfo->hasVelocity  = true;
    atomInfo->hasForce     = true;
    atomInfo->hasGlobalID  = true;

    FixedChargeAdapter fca = FixedChargeAdapter(atype);
    if (fca.isFixedCharge()) {
      atomInfo->hasCharge = true;
      atomInfo->charge    = fca.getCharge();
    }

    FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
    if (fqa.isFluctuatingCharge()) {
      atomInfo->hasCharge = true;
      atomInfo->charge += atom->getFlucQPos();
    }

    if ((storageLayout_ & DataStorage::dslElectricField) &&
        (atype->isElectrostatic())) {
      atomInfo->hasElectricField = true;
      atomInfo->eField           = atom->getElectricField();
    }

    atomData = std::make_shared<AtomData>();
    atomData->setID("ATOMDATA");
    atomData->addAtomInfo(atomInfo);

    atom->addProperty(atomData);

    setVisited(atom);
  }

  void DefaultAtomVisitor::visit(DirectionalAtom* datom) {
    std::shared_ptr<AtomData> atomData;
    std::shared_ptr<AtomInfo> atomInfo;
    AtomType* atype = datom->getAtomType();

    if (isVisited(datom)) return;

    atomInfo               = std::make_shared<AtomInfo>();
    atomInfo->atomTypeName = datom->getType();
    atomInfo->globalID     = datom->getGlobalIndex();
    atomInfo->pos          = datom->getPos();
    atomInfo->vel          = datom->getVel();
    atomInfo->frc          = datom->getFrc();
    atomInfo->hasVelocity  = true;
    atomInfo->hasForce     = true;
    atomInfo->hasGlobalID  = true;

    FixedChargeAdapter fca = FixedChargeAdapter(atype);
    if (fca.isFixedCharge()) {
      atomInfo->hasCharge = true;
      atomInfo->charge    = fca.getCharge();
    }

    FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
    if (fqa.isFluctuatingCharge()) {
      atomInfo->hasCharge = true;
      atomInfo->charge += datom->getFlucQPos();
    }

    if ((storageLayout_ & DataStorage::dslElectricField) &&
        (atype->isElectrostatic())) {
      atomInfo->hasElectricField = true;
      atomInfo->eField           = datom->getElectricField();
    }

    GayBerneAdapter gba = GayBerneAdapter(atype);
    MultipoleAdapter ma = MultipoleAdapter(atype);

    if (gba.isGayBerne()) {
      atomInfo->hasVector = true;
      atomInfo->vec       = datom->getA().transpose() * V3Z;
    } else if (ma.isDipole()) {
      atomInfo->hasVector = true;
      atomInfo->vec       = datom->getDipole();
    } else if (ma.isQuadrupole()) {
      atomInfo->hasVector = true;
      atomInfo->vec       = datom->getA().transpose() * V3Z;
    }

    atomData = std::make_shared<AtomData>();
    atomData->setID("ATOMDATA");
    atomData->addAtomInfo(atomInfo);

    datom->addProperty(atomData);

    setVisited(datom);
  }

  const std::string DefaultAtomVisitor::toString() {
    char buffer[65535];
    std::string result;

    snprintf(
        buffer, 65535,
        "--------------------------------------------------------------\n");
    result += buffer;

    snprintf(buffer, 65535, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    snprintf(buffer, 65535,
             "Visitor Description: copy atom infomation into atom data\n");
    result += buffer;

    snprintf(
        buffer, 65535,
        "--------------------------------------------------------------\n");
    result += buffer;

    return result;
  }
}  // namespace OpenMD
