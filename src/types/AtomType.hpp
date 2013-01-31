/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 * [1]  Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).             
 * [2]  Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).          
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 24107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
 
#ifndef TYPES_ATOMTYPE_HPP

#define TYPES_ATOMTYPE_HPP

#include <string>
#include "config.h"
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
    
    virtual ~AtomType() { } ;

    virtual void useBase(AtomType* base);
    virtual void copyAllData(AtomType* orig);   
    void setMass(RealType m);
    RealType getMass();        
    void setIdent(int id);    
    int getIdent();    
    void setName(const string&name);
    string getName();
    vector<AtomType*> allYourBase();
    vector<AtomType*> allYourZIG() {return everyZIG;}
    void addZig(AtomType* at) {everyZIG.push_back(at);}

    //below functions are just forward functions
    /**
     * Adds property into property map
     * @param genData GenericData to be added into PropertyMap
     */
    void addProperty(GenericData* genData);

    /**
     * Removes property from PropertyMap by name
     * @param propName the name of property to be removed
     */
    void removeProperty(const string& propName);

    /**
     * clear all of the properties
     */
    void clearProperties();

    /**
     * Returns all names of properties
     * @return all names of properties
     */
    vector<string> getPropertyNames();

    /**
     * Returns all of the properties in PropertyMap
     * @return all of the properties in PropertyMap
     */      
    vector<GenericData*> getProperties();

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
    GenericData* getPropertyByName(const string& propName);

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
    bool hasBase_; // All your base are belong to us
    AtomType* base_;
    vector< AtomType*> everyZIG;  // list of atom types which use us as a base
    map< string, bool> myResponsibilities_;
    map< string, RealType> myValues_;

  private:
    //prevent copy construction and copy assignment, since property
    //map contains pointers which can not be copied and managed
    //safely, except make the generic data at PropertyMap as copy on
    //write shared pointer
    AtomType(const AtomType&);
    AtomType& operator=(const AtomType& atomType);
    PropertyMap properties_;
  };  
}

#endif
