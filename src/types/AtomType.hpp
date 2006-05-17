/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#ifndef TYPES_ATOMTYPE_HPP

#define TYPES_ATOMTYPE_HPP

#include <string>
#include "config.h"
#include "utils/PropertyMap.hpp"

#define __C
#include "types/AtomTypeProperties.h"
#include "UseTheForce/DarkSide/atype_interface.h"

namespace oopse {
  /**
   * @class AtomType
   * AtomType is what OOPSE looks to for unchanging data about an atom.
   * Things that belong to AtomType are universal properties (i.e. does
   * this atom have a Charge?  What is it's mass_?)  Dynamic properties of
   * an atom are not intended to be properties of an atom type
   */
  class AtomType {
  public:

    AtomType();

    virtual ~AtomType() { } ;

    virtual void complete();

    /**
     * Finishes off the AtomType by communicating the logical portions of the
     * structure to the Fortran atype module
     */
    void makeFortranAtomType();
            
    void setMass(RealType m) {
      mass_ = m;
    }

    RealType getMass(void) {
      return mass_;
    }

    void setIdent(int id) {
      atp.ident = id;
    }

    int getIdent() {
      return atp.ident;
    }

    void setName(const std::string&name) {
      name_ = name;
    }

    std::string getName() {
      return name_;
    }

    void setLennardJones() {
      atp.is_LennardJones = 1;
    }

    bool isLennardJones() {
      return atp.is_LennardJones;
    }


    bool isElectrostatic() {
      return isCharge() || isMultipole();
    }

    void setEAM() {
      atp.is_EAM = 1;
    }

    bool isEAM() {
      return atp.is_EAM;
    }

    void setCharge() {
      atp.is_Charge = 1;
    }

    bool isCharge() {
      return atp.is_Charge;
    }

    bool isDirectional() {
      return atp.is_Directional;
    }

    bool isDipole() {
      return atp.is_Dipole;
    }

    bool isSplitDipole() {
      return atp.is_SplitDipole;
    }

    bool isQuadrupole() {
      return atp.is_Quadrupole;
    }
            
    bool isMultipole() {
      return isDipole() || isQuadrupole();
    }

    bool isGayBerne() {
      return atp.is_GayBerne;
    }

    bool isSticky() {
      return atp.is_Sticky;
    }
    
    bool isStickyPower() {
      return atp.is_StickyPower;
    }

    bool isShape() {
      return atp.is_Shape;
    }

    bool isSC() {
      return atp.is_SC;
    }
    
    void setSC() {
      atp.is_SC = 1;
    }
    
    
    bool isMEAM() {
      return atp.is_MEAM;
    }
    
    void setMEAM() {
      atp.is_MEAM = 1;
    }
    
    
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
    void removeProperty(const std::string& propName);

    /**
     * clear all of the properties
     */
    void clearProperties();

    /**
     * Returns all names of properties
     * @return all names of properties
     */
    std::vector<std::string> getPropertyNames();

    /**
     * Returns all of the properties in PropertyMap
     * @return all of the properties in PropertyMap
     */      
    std::vector<GenericData*> getProperties();

    /**
     * Returns property 
     * @param propName name of property
     * @return a pointer point to property with propName. If no property named propName
     * exists, return NULL
     */      
    GenericData* getPropertyByName(const std::string& propName);

  protected:

    AtomTypeProperties atp;
    RealType mass_;
    std::string name_;

  private:
    //prevent copy construction and copy assignment, since property map contains
    //pointers which can not be copied and managered safely, except make the generic data
    //at PropertyMap as copy on write shared pointer
    AtomType(const AtomType&);
    AtomType& operator=(const AtomType& atomType);

            
    PropertyMap properties_;
             
  };

  struct LJParam {
    RealType epsilon;
    RealType sigma;
    int soft_pot;
  };
  typedef SimpleTypeData<LJParam> LJParamGenericData;

  struct EAMParam {
    RealType latticeConstant;         
    int nrho;
    RealType drho;
    int nr;
    RealType dr;
    RealType rcut;
    std::vector<RealType> rvals; 
    std::vector<RealType> rhovals;
    std::vector<RealType> Frhovals;    
  };

  typedef SimpleTypeData<EAMParam> EAMParamGenericData;
  
  struct SCParam {
    RealType c;
    RealType m;
    RealType n;
    RealType alpha;
    RealType epsilon;
  };
  typedef SimpleTypeData<SCParam> SCParamGenericData;
  
  
  
  
}

#endif
