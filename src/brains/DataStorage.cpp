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
  * @file DataStorage.cpp
  * @author tlin
  * @date 10/26/2004
  * @time 11:56am
  * @version 1.0
  */

#include "brains/DataStorage.hpp"


DataStorage::DataStorage() : size_(0), storageLayout_(0){

}

DataStorage::DataStorage(int size, int storageLayout) : size_(size){
    setStorageLayout(storageLayout);
    resize(size);
}

int DataStorage::getSize() {

    if (storageLayout_ & dslPosition && position.size() != size_) {
        //error
        std::cerr << "size does not match"<< std::endl;
    }

    if (storageLayout_ & dslVelocity && velocity.size() != size_) {
        //error
        std::cerr << "size does not match"<< std::endl;        
    }

    if (storageLayout_ & dslAmat && aMat.size() != size_) {
        //error
        std::cerr << "size does not match"<< std::endl;        
    }

    if (storageLayout_ & dslAngularMomentum && angularMomentum.size() != size_) {
        //error
        std::cerr << "size does not match"<< std::endl;        
    }

    if (storageLayout_ & dslUnitVector && unitVector.size() != size_) {
        //error
        std::cerr << "size does not match"<< std::endl;        
    }

    if (storageLayout_ & dslZAngle && zAngle.size() != size_) {
        //error
        std::cerr << "size does not match"<< std::endl;        
    }

    if (storageLayout_ & dslForce && force.size() != size_) {
        //error
        std::cerr << "size does not match"<< std::endl;        
    }

    if (storageLayout_ & dslTorque && torque.size() != size_) {
        //error
        std::cerr << "size does not match"<< std::endl;        
    }

    return size_;

}

void DataStorage::resize(int newSize) {

    if (storageLayout_ & dslPosition) {
        internalResize(position, newSize);
    }

    if (storageLayout_ & dslVelocity) {
        internalResize(velocity, newSize);
    }

    if (storageLayout_ & dslAmat) {
        internalResize(aMat, newSize);
    }

    if (storageLayout_ & dslAngularMomentum) {
        internalResize(angularMomentum, newSize);
    }

    if (storageLayout_ & dslUnitVector) {
        internalResize(unitVector, newSize);
    }
    
    if (storageLayout_ & dslZAngle) {
        internalResize(zAngle, newSize);
    }

    if (storageLayout_ & dslForce) {
        internalResize(force, newSize);
    }

    if (storageLayout_ & dslTorque) {
        internalResize(torque, newSize);
    }

    size_ = newSize;
}

void DataStorage::reserve(int size) {
    if (storageLayout_ & dslPosition) {
        position.reserve(size);
    } 

    if (storageLayout_ & dslVelocity) {
        velocity.reserve(size);
    } 

    if (storageLayout_ & dslAmat) {
        aMat.reserve(size);
    } 

    if (storageLayout_ & dslAngularMomentum) {
        angularMomentum.reserve(size);
    } 

    if (storageLayout_ & dslUnitVector) {
        unitVector.reserve(size);
    }
    
    if (storageLayout_ & dslZAngle) {
        zAngle.reserve(size);
    }

    if (storageLayout_ & dslForce) {
        force.reserve(size);
    } 

    if (storageLayout_ & dslTorque) {
        torque.reserve(size);
    }

}

void DataStorage::copy(int source, int num, int target) {
    if (num + target > size_ ) {
        //error
    }
    
    if (storageLayout_ & dslPosition) {
        interalCopy(position, source, num, target);
    } 

    if (storageLayout_ & dslVelocity) {
        interalCopy(velocity, source, num, target);
  } 

    if (storageLayout_ & dslAmat) {
        interalCopy(aMat, source, num, target);
   } 

    if (storageLayout_ & dslAngularMomentum) {
        interalCopy(angularMomentum, source, num, target);
    } 

    if (storageLayout_ & dslUnitVector) {
        interalCopy(unitVector, source, num, target);
    }
    
    if (storageLayout_ & dslZAngle) {
        interalCopy(zAngle, source, num, target);
    }

    if (storageLayout_ & dslForce) {
        interalCopy(force, source, num, target);
    } 

    if (storageLayout_ & dslTorque) {
        interalCopy(torque, source, num, target); 
    }
    

}

int DataStorage::getStorageLayout() {
    return storageLayout_;
}

void DataStorage::setStorageLayout(int layout) {
    storageLayout_ = layout;
    resize(size_);
}

double* DataStorage::getArrayPointer(int whichArray) {

    switch (whichArray) {
        case dslPosition:
            return internalGetArrayPointer(torque);
            break;
            
        case dslVelocity:
            return internalGetArrayPointer(velocity);
            break;
            
        case dslAmat:
            return internalGetArrayPointer(aMat);
            break;            
            
        case dslAngularMomentum:
            return internalGetArrayPointer(angularMomentum);
            break;
            
        case dslUnitVector:
            return internalGetArrayPointer(unitVector);
            break;
            
        case dslZAngle:
            return internalGetArrayPointer(zAngle);
            break;

        case dslForce:
            return internalGetArrayPointer(force);
            break;            

        case dslTorque:
            return internalGetArrayPointer(torque);
            break;  
            
        default:
            //error message
            return NULL;

    }
}    

double* DataStorage::internalGetArrayPointer(std::vector<Vector3d>& v) {
    if (v.size() == 0) {
        return NULL;
    } else {
        return v[0].getArrayPointer();
    }
}

double* DataStorage::internalGetArrayPointer(std::vector<RotMat3x3d>& v) {
    if (v.size() == 0) {
        return NULL;
    } else {
        return v[0].getArrayPointer();
    }

}

double* DataStorage::internalGetArrayPointer(std::vector<double>& v) {
    if (v.size() == 0) {
        return NULL;
    } else {
        return &(v[0]);
    }

}    

template<typename T>
void DataStorage::internalResize(std::vector<T>& v, int newSize){
    int oldSize = v.size();

    if (oldSize == newSize) {
        return;
    } else if (oldSize < newSize) {
        v.insert(v.end(), newSize-oldSize, T());
    } else {
        typename std::vector<T>::iterator i;
        i = v.begin();
        std::advance(i, newSize);        
        v.erase(i, v.end());
    }
}

template<typename T>
void DataStorage::interalCopy(std::vector<T>& v, int source,  int num, int target) {
    typename std::vector<T>::iterator first;
    typename std::vector<T>::iterator last;
    typename std::vector<T>::iterator result;

    first = v.begin();
    last = v.begin();
    result = v.begin();

    std::advance(first, source);
    //STL algorithm use half opened range
    std::advance(last, num + 1);
    std::advance(result, target );

    std::copy(first, last, result);
}
