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
 
#include "applications/hydrodynamics/HydrodynamicsModelFactory.hpp"
#include "applications/hydrodynamics/HydrodynamicsModelCreator.hpp"
#include "applications/hydrodynamics/HydrodynamicsModel.hpp"
#include "brains/SimInfo.hpp"
#include "utils/MemoryUtils.hpp"

namespace OpenMD {

  //initialize instance of HydrodynamicsModelFactory
  HydrodynamicsModelFactory* HydrodynamicsModelFactory::instance_ = NULL;

  HydrodynamicsModelFactory::~HydrodynamicsModelFactory() {
    MemoryUtils::deletePointers(creatorMap_);
  }

  bool HydrodynamicsModelFactory::registerHydrodynamicsModel(HydrodynamicsModelCreator* creator) {
    return creatorMap_.insert(
			      CreatorMapType::value_type(creator->getIdent(), creator)).second;
  }

  bool HydrodynamicsModelFactory::unregisterHydrodynamicsModel(const std::string& id) {
    return creatorMap_.erase(id) == 1;
  }

  HydrodynamicsModel* HydrodynamicsModelFactory::createHydrodynamicsModel(const std::string& id, StuntDouble* sd, SimInfo* info) {
    CreatorMapType::iterator i = creatorMap_.find(id);
    if (i != creatorMap_.end()) {
      //invoke functor to create object
      return (i->second)->create(sd, info);
    } else {
      return NULL;
    }
  }

  std::vector<std::string> HydrodynamicsModelFactory::getIdents() {
    IdentVectorType idents;
    CreatorMapType::iterator i;

    for (i = creatorMap_.begin(); i != creatorMap_.end(); ++i) {
      idents.push_back(i->first);
    }
    
    return idents;
  }

  std::ostream& operator <<(std::ostream& o, HydrodynamicsModelFactory& factory) {
    HydrodynamicsModelFactory::IdentVectorType idents;
    HydrodynamicsModelFactory::IdentVectorIterator i;

    idents = factory.getIdents();

    o << "Avaliable type identifiers in this factory: " << std::endl;
    for (i = idents.begin(); i != idents.end(); ++i) {
      o << *i << std::endl;
    }

    return o;
  }

}

