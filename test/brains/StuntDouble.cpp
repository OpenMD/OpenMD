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
 * @file StuntDouble.cpp
 * @author    tlin
 * @date  10/22/2004
 * @version 1.0
 */ 

#include "primitives/StuntDouble.hpp"

namespace oopse {

StuntDouble::StuntDouble(ObjectType objType, DataStoragePointer storage) : 
    objType_(objType), storage_(storage),    
    linear_(false), linearAxis_(-1), globalIndex_(-1), localIndex_(-1), snapshotMan_(NULL){
}

StuntDouble::~StuntDouble() {

}

void StuntDouble::zeroForces() {
    setFrc(V3Zero);
}
void StuntDouble::addProperty(GenericData* genData) {
    properties_.addProperty(genData);  
}

void StuntDouble::removeProperty(std::string& propName) {
    properties_.removeProperty(propName);  
}

void StuntDouble::clearProperties() {
    properties_.clearProperties(); 
}

std::vector<std::string> StuntDouble::getPropertyNames() {
    return properties_.getPropertyNames();  
}
      
std::vector<GenericData*> StuntDouble::getProperties() { 
    return properties_.getProperties(); 
}

GenericData* StuntDouble::getPropertyByName(std::string& propName) {
    return properties_.getPropertyByName(propName); 
}


}
