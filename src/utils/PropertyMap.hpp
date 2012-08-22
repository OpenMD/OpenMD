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
 
/**
 * @file PropertyMap.hpp
 * @author tlin
 * @date 09/21/2004
 * @time 9:20am
 * @version 1.0
 */
 
#ifndef UTIL_PROPERTYMAP_HPP
#define UTIL_PROPERTYMAP_HPP

#include <map>
#include <string>
#include <vector>

#include <utils/GenericData.hpp>

namespace OpenMD{

  /**
   * @class PropertyMap
   * PropertyMap class maintains a list of GenericData. Type of Property is actually GenericData.
   */
  class PropertyMap{
  public:

    /** trivial constructor */
    PropertyMap() {}

    /**
     * Virtual Destructor responsible for deleting all of the generc data in PropertyMap
     */
    virtual ~PropertyMap();

    /**
     * Adds property into property map
     *
     * @param genData GenericData to be added into PropertyMap
     *
     * @see #removeProperty
     * @see #clearProperties
     */
    void addProperty(GenericData* genData);

    /**
     * Removes property from PropertyMap by name
     *
     * @param propName the name of property to be removed
     *
     * @see #addProperty
     * @see #clearProperties
     */
    bool removeProperty(const std::string& propName);

    /**
     * clear all of the properties
     *
     * @see #addProperty
     * @see #removeProperty
     */
    void clearProperties();

    /**
     * Returns all names of properties
     *
     * @return all names of properties
     */
    std::vector<std::string> getPropertyNames();

    /**
     * Returns all of the properties in PropertyMap
     *
     * @return all of the properties in PropertyMap
     *
     * @see #getPropertyByName
     */      
    std::vector<GenericData*> getProperties();

    /**
     * Checks if property is in this PropertyMap
     *
     * @param propName name of property
     *
     * @return boolean
     *
     * @see #getProperties, #getPropertyByName
     */      
    bool hasProperty(const std::string& propName);

    /**
     * Returns property 
     *
     * @param propName name of property
     *
     * @return a pointer point to property with propName. If no property named propName
     * exists, return NULL
     *
     * @see #getProperties
     */      
    GenericData* getPropertyByName(const std::string& propName);
    //template<typename T = GenericData*> T getPropertyByName(std:std::string& propName);

  protected:
    std::map<std::string, GenericData*> propMap_;

  private:

    /** prevent copy constructing */
    PropertyMap(const PropertyMap&);

    /** prevent copy assignment */
    PropertyMap& operator=(const PropertyMap&);
  };

}// namespace OpenMD

#endif //UTIL_PROPERTYMAP_HPP
