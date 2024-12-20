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

/**
 * @file PropertyMap.hpp
 * @author tlin
 * @date 09/21/2004
 * @version 1.0
 */

#ifndef UTIL_PROPERTYMAP_HPP
#define UTIL_PROPERTYMAP_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <utils/GenericData.hpp>

namespace OpenMD {

  /**
   * @class PropertyMap
   * PropertyMap class maintains a list of GenericData. Type of Property is
   * actually GenericData.
   */
  class PropertyMap {
  public:
    /** trivial constructor */
    PropertyMap() {}

    /**
     * Virtual Destructor responsible for deleting all of the generc data in
     * PropertyMap
     */
    virtual ~PropertyMap() = default;

    /**
     * Adds property into property map
     *
     * @param genData GenericData to be added into PropertyMap
     *
     * @see #removeProperty
     * @see #clearProperties
     */
    void addProperty(std::shared_ptr<GenericData> genData);

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
    std::vector<std::shared_ptr<GenericData>> getProperties();

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
     * @return a pointer point to property with propName. If no property named
     * propName exists, return NULL
     *
     * @see #getProperties
     */
    std::shared_ptr<GenericData> getPropertyByName(const std::string& propName);
    // template<typename T = GenericData*> T getPropertyByName(std:std::string&
    // propName);

  protected:
    std::map<std::string, std::shared_ptr<GenericData>> propMap_;

  private:
    /** prevent copy constructing */
    PropertyMap(const PropertyMap&);

    /** prevent copy assignment */
    PropertyMap& operator=(const PropertyMap&);
  };
}  // namespace OpenMD

#endif  // UTIL_PROPERTYMAP_HPP
