/*
 * Copyright (C) 2000-2004  Object Oriented Parallel Simulation Engine (OOPSE) project
 * 
 * Contact: oopse@oopse.org
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
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

#include <GenericData.hpp>

namespace oopse{

  /**
   * @class PropertyMap
   * PropertyMap class maintains a list of GenericData. Type of Property is actually GenericData.
   */
  class PropertyMap{
    public:
      
      /** trivial constructor */
      PropertyMap();

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
      void removeProperty(std::string& propName);

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
       * Returns property 
       *
       * @param propName name of property
       *
       * @return a pointer point to property with propName. If no property named propName
       * exists, return NULL
       *
       * @see #getProperties
       */      
      GenericData* getPropertyByName(std:string& propName);
      template<typename T = GenericData*> T getPropertyByName(std:string& propName);
      
    protected:
      std::map<std::string, GenericData*> propMap_;

    private:

      /** prevent copy constructing */
      PropertyMap(const PropertyMap&);

      /** prevent copy assignment */
      PropertyMap& operator=(const PropertyMap&);
  };

}// namespace oopse
#endif //UTIL_PROPERTYMAP_HPP