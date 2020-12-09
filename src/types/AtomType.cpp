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
 
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>

#include "types/AtomType.hpp"
#include "types/StickyAdapter.hpp"
#include "types/MultipoleAdapter.hpp"

using namespace std;

namespace OpenMD {
  AtomType::AtomType(){

    // initially, all atom types are their own base types.
    base_ = this;
    hasBase_ = false;

    // initialize to an error:
    ident_ = -1;
    
    // and massless:
    mass_ = 0.0;    
    myResponsibilities_["mass"] = false;
  }
  
  void AtomType::useBase(AtomType* base) {
    hasBase_=true;
    base_ = base;
    base->addZig(this);
  }

  void AtomType::copyAllData(AtomType* orig) {

    // makes an exact replica of another atom type down to the
    // atom ID and any base attachments it may have.  
    // use with caution!
    
    hasBase_=orig->hasBase_;
    base_ = orig->base_;
    mass_ = orig->mass_;    
    name_ = string(orig->name_);
    ident_ = orig->ident_;

    map< string, bool>::iterator i;;
    map< string, RealType>::iterator j;

    for (i = orig->myResponsibilities_.begin(); 
         i != orig->myResponsibilities_.end(); ++i) {
      myResponsibilities_[(*i).first] = orig->myResponsibilities_[(*i).first];
    }
    
    for (j = orig->myValues_.begin(); j != orig->myValues_.end(); ++j) {
      myValues_[(*j).first] = orig->myValues_[(*j).first];
    }
    
    std::vector< std::shared_ptr<GenericData> > oprops = orig->getProperties();
    std::vector< std::shared_ptr<GenericData> >::iterator it;

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
  
  std::vector<std::shared_ptr<GenericData> > AtomType::getProperties() { 
    return properties_.getProperties(); 
  }

  bool AtomType::hasProperty(const string& propName) {
    if (hasBase_ && !myResponsibilities_[propName]){
      return base_->hasProperty(propName);}
    else
      return properties_.hasProperty(propName); 
  }  
  
  std::shared_ptr<GenericData> AtomType::getPropertyByName(const string& propName) {
    if (hasBase_ && !myResponsibilities_[propName]){
      return base_->getPropertyByName(propName);}
    else
      return properties_.getPropertyByName(propName); 
  }  

  void AtomType::setMass(RealType m) {
    myResponsibilities_["mass"] = true;
    mass_ = m;
  }
  
  RealType AtomType::getMass(void) {
    if (hasBase_ && !myResponsibilities_["mass"]) 
      return base_->getMass();
    else
      return mass_;
  }
  
  void AtomType::setIdent(int id) {
    ident_ = id;
  }
  
  int AtomType::getIdent() {
    return ident_;
  }
  
  void AtomType::setName(const string&name) {
    name_ = name;
  }
  
  string AtomType::getName() {
    return name_;
  }

  bool AtomType::isLennardJones() {
    return hasProperty("LJ");
  }

  bool AtomType::isElectrostatic() {
    return isCharge() || isMultipole();
  }
  
  bool AtomType::isEAM() {
    return hasProperty("EAM");
  }
  
  bool AtomType::isCharge() {
    return isFixedCharge() || isFluctuatingCharge();
  }

  bool AtomType::isDirectional() {
    return hasProperty("Directional");
  }

  bool AtomType::isFluctuatingCharge() {
    return hasProperty("FlucQ");
  }

  bool AtomType::isFixedCharge() {
    return hasProperty("Charge");
  }

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
      
  bool AtomType::isMultipole() {
    return hasProperty("Multipole");
  }

  bool AtomType::isGayBerne() {
    return hasProperty("GB");
  }

  bool AtomType::isSticky() {
    return hasProperty("Sticky");
  }

  bool AtomType::isStickyPower() {
    StickyAdapter sa = StickyAdapter(this);
    return sa.isStickyPower();
  }

  bool AtomType::isShape() {
    return hasProperty("Shape");
  }
  
  bool AtomType::isSC() {
    return hasProperty("SC");
  }

  bool AtomType::isMetal() {
    return isSC() || isEAM();
  }

  std::vector<AtomType* > AtomType::allYourBase() {   
    std::vector<AtomType* > myChain;   

    if(hasBase_){
      myChain = base_->allYourBase();
      myChain.insert(myChain.begin(), this);
    } else {
      myChain.push_back(this);
    }
    
    return myChain;
  }

}
