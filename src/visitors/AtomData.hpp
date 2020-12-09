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
 
#ifndef VISITOR_ATOMDATA_HPP
#define VISITOR_ATOMDATA_HPP
#include <string>
#include <vector>

#include "utils/GenericData.hpp"
#include "math/Vector3.hpp"

namespace OpenMD {

  struct AtomInfo {
    AtomInfo() : pos(V3Zero), vec(V3Zero), vel(V3Zero), frc(V3Zero),
                 eField(V3Zero), charge(0.0),
                 hasCharge(false), hasVector(false), hasVelocity(false), 
                 hasForce(false), hasElectricField(false), hasGlobalID(false) {}
    
    std::string atomTypeName;
    int globalID;
    Vector3d pos;
    Vector3d vec;  
    Vector3d vel;  
    Vector3d frc;  
    Vector3d eField;
    RealType charge;
    bool hasCharge;
    bool hasVector;
    bool hasVelocity;
    bool hasForce;
    bool hasElectricField;
    bool hasGlobalID;
  };

  class AtomData : public GenericData{
  public:

    AtomData(const std::string& id = "ATOMDATA") : GenericData(id) {}
        
    void addAtomInfo(AtomInfo* info) {data.push_back(info);}

    void clearAllAtomInfo();

    AtomInfo* beginAtomInfo(std::vector<AtomInfo*>::iterator& i){
      i = data.begin();
      return i != data.end()? *i : NULL;
    }

    AtomInfo* nextAtomInfo(std::vector<AtomInfo*>::iterator& i){
      ++i;
      return i != data.end()? *i: NULL;
    }

    std::vector<AtomInfo*> getData() {return data;}

    int getSize() {return data.size();}

  protected:

    std::vector<AtomInfo*> data;
  };
}
#endif //VISITOR_ATOMDATA_HPP
