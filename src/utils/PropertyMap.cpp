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
 
#include "PropertyMap.hpp"
#include <cassert>

PropertyMap::~PropertyMap(){
  clearProperties();
}


PropertyMap::addProperty(GenericData* genData){
  std::map<std::string, GenericData*>::iterator iter;

  assert( genData != NULL);

  iter = propMap_.find(genData->getID());

  if (iter == propMap_.end()){
    propMap_.insert(make_pair(genData->getID(), genData));
  }
  else
    logger.warn("");
  
}

void PropertyMap::removeProperty(std::string& propName){
  std::map<std::string, GenericData*>::iterator iter;

  iter = propMap_.find(propName);

  if (iter != propMap.end()){
    delete iter->second;
    propMap_.erase(iter);    
  }
  else
    logger.warn("Can not find property with name: " + propName);
}

void PropertyMap::clearProperties(){
  std::map<std::string, GenericData*>::iterator iter;

  for (iter = propMap_.begin(); iter != propMap_.end(); ++iter)
    delete iter->second;
  
  propMap_.clear();
}

std::vector<std::string> PropertyMap::getPropertyNames(){
  vector<string> propNames;
  std::map<std::string, GenericData*>::iterator iter;

  for (iter = propMap_.begin(); iter != propMap_.end(); ++iter)
    propNames.push_back(iter->first);

  return propNames;
}

std::vector<GenericData*> PropertyMap::getProperties(){
  vector<GenericData*> properties;
  std::map<std::string, GenericData*>::iterator iter;

  for (iter = propMap_.begin(); iter != propMap_.end(); ++iter)
    properties.push_back(iter->second);

  return properties;
}

GenericData* PropertyMap::getPropertyByName(std:string& propName){
  std::map<std::string, GenericData*>::iterator iter;

  iter = propMap_.find(propName);

  if (iter != propMap_.end())
    return iter->second;
  else
    return NULL;
}        
