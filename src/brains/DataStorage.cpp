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
 
 /**
  * @file DataStorage.cpp
  * @author tlin
  * @date 10/26/2004
  * @time 11:56am
  * @version 1.0
  */

#include "brains/DataStorage.hpp"

namespace oopse {

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

    if (storageLayout_ & dslElectroFrame && electroFrame.size() != size_) {
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

    if (storageLayout_ & dslElectroFrame) {
        internalResize(electroFrame, newSize);
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

    if (storageLayout_ & dslElectroFrame) {
        electroFrame.reserve(size);
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
        internalCopy(position, source, num, target);
    } 

    if (storageLayout_ & dslVelocity) {
        internalCopy(velocity, source, num, target);
  } 

    if (storageLayout_ & dslAmat) {
        internalCopy(aMat, source, num, target);
   } 

    if (storageLayout_ & dslAngularMomentum) {
        internalCopy(angularMomentum, source, num, target);
    } 

    if (storageLayout_ & dslElectroFrame) {
        internalCopy(electroFrame, source, num, target);
    }
    
    if (storageLayout_ & dslZAngle) {
        internalCopy(zAngle, source, num, target);
    }

    if (storageLayout_ & dslForce) {
        internalCopy(force, source, num, target);
    } 

    if (storageLayout_ & dslTorque) {
        internalCopy(torque, source, num, target); 
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
            return internalGetArrayPointer(position);
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
            
        case dslElectroFrame:
            return internalGetArrayPointer(electroFrame);
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
void DataStorage::internalCopy(std::vector<T>& v, int source,  int num, int target) {
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

int DataStorage::getBytesPerStuntDouble(int layout) {
    int bytes = 0;
    if (layout & dslPosition) {
        bytes += sizeof(Vector3d);
    } else if (layout & dslVelocity) {
        bytes += sizeof(Vector3d);
    } else if (layout & dslAmat) {
        bytes += sizeof(Mat3x3d);    
    } else if (layout & dslAngularMomentum) {
        bytes += sizeof(Vector3d);
    } else if (layout & dslElectroFrame) {
        bytes += sizeof(Mat3x3d);
    } else if (layout & dslZAngle) {
        bytes += sizeof(Vector3d);
    } else if (layout & dslForce) {
        bytes += sizeof(Vector3d);
    } else if (layout & dslTorque) {
        bytes += sizeof(Vector3d);
    }
    return bytes;
}

}
