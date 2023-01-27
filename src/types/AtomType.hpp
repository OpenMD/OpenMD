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

#ifndef TYPES_ATOMTYPE_HPP

#define TYPES_ATOMTYPE_HPP

#include <config.h>

#include <set>
#include <string>

#include "utils/PropertyMap.hpp"

using namespace std;

namespace OpenMD {
  /**
   * @class AtomType
   * AtomType is what OpenMD looks to for unchanging data about an atom.
   * Things that belong to AtomType are universal properties (i.e. does
   * this atom have a Charge?  What is it's mass_?)  Dynamic properties of
   * an atom are not intended to be properties of an atom type
   */
  class AtomType {
  public:
    AtomType();

    virtual ~AtomType() {};

    virtual void useBase(AtomType* base);
    virtual void copyAllData(AtomType* orig);
    void setMass(RealType m);
    RealType getMass();
    void setIdent(int id);
    int getIdent();
    void setName(const string& name);
    string getName();
    std::vector<AtomType*> allYourBase();
    std::vector<AtomType*> allYourZIG() { return everyZIG; }
    void addZig(AtomType* at) { everyZIG.push_back(at); }

    // below functions are just forward functions
    /**
     * Adds property into property map
     * @param genData GenericData to be added into PropertyMap
     */
    void addProperty(std::shared_ptr<GenericData> genData);

    /**
     * Removes property from PropertyMap by name
     * @param propName the name of property to be removed
     */
    void removeProperty(const string& propName);

    /**
     * Returns all names of properties
     * @return all names of properties
     */
    std::vector<string> getPropertyNames();

    /**
     * Returns all of the properties in PropertyMap
     * @return all of the properties in PropertyMap
     */
    std::vector<std::shared_ptr<GenericData>> getProperties();

    /**
     * Checks if property is in this PropertyMap
     * @param propName name of property
     * @return boolean
     */
    bool hasProperty(const string& propName);

    /**
     * Returns property
     * @param propName name of property
     * @return a pointer point to property with propName. If no
     * property named propName exists, return NULL
     */
    std::shared_ptr<GenericData> getPropertyByName(const string& propName);

    bool isLennardJones();
    bool isElectrostatic();
    bool isEAM();
    bool isCharge();
    bool isFixedCharge();
    bool isDirectional();
    bool isDipole();
    bool isQuadrupole();
    bool isMultipole();
    bool isGayBerne();
    bool isSticky();
    bool isStickyPower();
    bool isShape();
    bool isSC();
    bool isMetal();
    bool isFluctuatingCharge();

  protected:
    int ident_;
    RealType mass_;
    string name_;
    bool hasBase_;  // All your base are belong to us
    AtomType* base_;
    vector<AtomType*> everyZIG;  // list of atom types which use us as a base
    map<string, bool> myResponsibilities_;
    map<string, RealType> myValues_;

  private:
    // prevent copy construction and copy assignment, since property
    // map contains pointers which can not be copied and managed
    // safely, except make the generic data at PropertyMap as copy on
    // write shared pointer
    AtomType(const AtomType&);
    AtomType& operator=(const AtomType& atomType);
    PropertyMap properties_;
  };

  class AtomTypeCompare {
  public:
    bool operator()(AtomType* lhs, AtomType* rhs) const {
      return lhs->getIdent() < rhs->getIdent();
    }
  };

  using AtomTypeSet = std::set<AtomType*, AtomTypeCompare>;
}  // namespace OpenMD

#endif
