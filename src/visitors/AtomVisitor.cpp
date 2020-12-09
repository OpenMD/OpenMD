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
 
#include <cstring>
#include <memory>

#include "visitors/AtomVisitor.hpp"
#include "primitives/DirectionalAtom.hpp"
#include "primitives/RigidBody.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/MultipoleAdapter.hpp"
#include "types/GayBerneAdapter.hpp"

namespace OpenMD {

  BaseAtomVisitor::BaseAtomVisitor(SimInfo* info) : BaseVisitor() {
    storageLayout_ = info->getStorageLayout(); 
  }    
  
  void BaseAtomVisitor::visit(RigidBody *rb) {
    //vector<Atom*> myAtoms;
    //vector<Atom*>::iterator atomIter;

    //myAtoms = rb->getAtoms();

    //for(atomIter = myAtoms.begin(); atomIter != myAtoms.end(); ++atomIter)
    //  (*atomIter)->accept(this);
  }

  void BaseAtomVisitor::setVisited(Atom *atom) {
    std::shared_ptr<GenericData> data;
    data = atom->getPropertyByName("VISITED");

    //if visited property is not existed, add it as new property
    if (data == nullptr) {
      data = std::make_shared<GenericData>();
      data->setID("VISITED");
      atom->addProperty(data);
    }
  }

  bool BaseAtomVisitor::isVisited(Atom *atom) {
    std::shared_ptr<GenericData> data;
    data = atom->getPropertyByName("VISITED");
    return data == nullptr ? false : true;
  }

  //------------------------------------------------------------------------//
	
  void DefaultAtomVisitor::visit(Atom *atom) {
    std::shared_ptr<AtomData> atomData;
    AtomInfo* atomInfo;
    AtomType* atype = atom->getAtomType();
              
    if (isVisited(atom))
      return;
    
    atomInfo = new AtomInfo();
    atomInfo->atomTypeName = atom->getType();
    atomInfo->globalID = atom->getGlobalIndex();
    atomInfo->pos = atom->getPos();
    atomInfo->vel = atom->getVel();
    atomInfo->frc = atom->getFrc();
    atomInfo->vec = V3Zero;
    atomInfo->hasVelocity = true;
    atomInfo->hasForce = true;
    atomInfo->hasGlobalID = true;

        
    FixedChargeAdapter fca = FixedChargeAdapter(atype);
    if ( fca.isFixedCharge() ) {
      atomInfo->hasCharge = true;
      atomInfo->charge = fca.getCharge();
    }
          
    FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
    if ( fqa.isFluctuatingCharge() ) {
      atomInfo->hasCharge = true;
      atomInfo->charge += atom->getFlucQPos();
    }
    
    if ((storageLayout_ & DataStorage::dslElectricField) && 
        (atype->isElectrostatic())) {
      atomInfo->hasElectricField = true;
      atomInfo->eField = atom->getElectricField();
    }

    atomData = std::make_shared<AtomData>();
    atomData->setID("ATOMDATA");   
    atomData->addAtomInfo(atomInfo);
    
    atom->addProperty(atomData);
    
    setVisited(atom);
  }
  
  void DefaultAtomVisitor::visit(DirectionalAtom *datom) {
    std::shared_ptr<AtomData> atomData;
    AtomInfo *atomInfo;
    AtomType* atype = datom->getAtomType();

    if (isVisited(datom))
      return;
    
    atomInfo = new AtomInfo;
    atomInfo->atomTypeName = datom->getType();
    atomInfo->globalID = datom->getGlobalIndex();
    atomInfo->pos = datom->getPos();
    atomInfo->vel = datom->getVel();
    atomInfo->frc = datom->getFrc();
    atomInfo->hasVelocity = true;
    atomInfo->hasForce = true;
    atomInfo->hasGlobalID = true;


    FixedChargeAdapter fca = FixedChargeAdapter(atype);
    if ( fca.isFixedCharge() ) {
      atomInfo->hasCharge = true;
      atomInfo->charge = fca.getCharge();
    }
          
    FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
    if ( fqa.isFluctuatingCharge() ) {
      atomInfo->hasCharge = true;
      atomInfo->charge += datom->getFlucQPos();
    }

    if ((storageLayout_ & DataStorage::dslElectricField) && 
        (atype->isElectrostatic())) {
      atomInfo->hasElectricField = true;
      atomInfo->eField = datom->getElectricField();
    }

    GayBerneAdapter gba = GayBerneAdapter(atype);
    MultipoleAdapter ma = MultipoleAdapter(atype);
    
    if (gba.isGayBerne()) {
      atomInfo->hasVector = true;
      atomInfo->vec = datom->getA().transpose()*V3Z;
    } else if (ma.isDipole()) {
      atomInfo->hasVector = true;
      atomInfo->vec = datom->getDipole();
    } else if (ma.isQuadrupole()) {
      atomInfo->hasVector = true;
      atomInfo->vec = datom->getA().transpose()*V3Z;
    }

    atomData = std::make_shared<AtomData>();
    atomData->setID("ATOMDATA");   
    atomData->addAtomInfo(atomInfo);

    datom->addProperty(atomData);

    setVisited(datom);
  }

  const std::string DefaultAtomVisitor::toString() {
    char   buffer[65535];
    std::string result;

    sprintf(buffer,
            "--------------------------------------------------------------\n");
    result += buffer;

    sprintf(buffer, "Visitor name: %s\n", visitorName.c_str());
    result += buffer;

    sprintf(buffer,
            "Visitor Description: copy atom infomation into atom data\n");
    result += buffer;

    sprintf(buffer,
            "--------------------------------------------------------------\n");
    result += buffer;

    return result;
  }
} //namespace OpenMD
