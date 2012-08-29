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
 * @file DataStorage.cpp
 * @author tlin
 * @date 10/26/2004
 * @time 11:56am
 * @version 1.0
 */

#include "brains/DataStorage.hpp"
using namespace std;
namespace OpenMD {

  DataStorage::DataStorage() : size_(0), storageLayout_(0){

  }

  DataStorage::DataStorage(int size, int storageLayout) : size_(size){
    setStorageLayout(storageLayout);
    resize(size);
  }

  int DataStorage::getSize() {

    if (storageLayout_ & dslPosition && position.size() != size_) {
      //error
      cerr << "size does not match"<< endl;
    }

    if (storageLayout_ & dslVelocity && velocity.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslForce && force.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslAmat && aMat.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslAngularMomentum && angularMomentum.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslTorque && torque.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslParticlePot && particlePot.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslDensity && density.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslFunctional && functional.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslFunctionalDerivative && functionalDerivative.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslDipole && dipole.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslQuadrupole && quadrupole.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslElectricField && electricField.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslSkippedCharge && skippedCharge.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslFlucQPosition && flucQPos.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslFlucQVelocity && flucQVel.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
    }

    if (storageLayout_ & dslFlucQForce && flucQFrc.size() != size_) {
      //error
      cerr << "size does not match"<< endl;        
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

    if (storageLayout_ & dslForce) {
      internalResize(force, newSize);
    }

    if (storageLayout_ & dslAmat) {
      internalResize(aMat, newSize);
    }

    if (storageLayout_ & dslAngularMomentum) {
      internalResize(angularMomentum, newSize);
    }
    
    if (storageLayout_ & dslTorque) {
      internalResize(torque, newSize);
    }

    if (storageLayout_ & dslParticlePot) {
      internalResize(particlePot, newSize);
    }

    if (storageLayout_ & dslDensity) {
      internalResize(density, newSize);
    }

    if (storageLayout_ & dslFunctional) {
      internalResize(functional, newSize);
    }

    if (storageLayout_ & dslFunctionalDerivative) {
      internalResize(functionalDerivative, newSize);
    }

    if (storageLayout_ & dslDipole) {
      internalResize(dipole, newSize);
    }

    if (storageLayout_ & dslQuadrupole) {
      internalResize(quadrupole, newSize);
    }

    if (storageLayout_ & dslElectricField) {
      internalResize(electricField, newSize);
    }

    if (storageLayout_ & dslSkippedCharge) {
      internalResize(skippedCharge, newSize);
    }

    if (storageLayout_ & dslFlucQPosition) {
      internalResize(flucQPos, newSize);
    }

    if (storageLayout_ & dslFlucQVelocity) {
      internalResize(flucQVel, newSize);
    }

    if (storageLayout_ & dslFlucQForce) {
      internalResize(flucQFrc, newSize);
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

    if (storageLayout_ & dslForce) {
      force.reserve(size);
    } 

    if (storageLayout_ & dslAmat) {
      aMat.reserve(size);
    } 

    if (storageLayout_ & dslAngularMomentum) {
      angularMomentum.reserve(size);
    } 

    if (storageLayout_ & dslTorque) {
      torque.reserve(size);
    }
    
    if (storageLayout_ & dslParticlePot) {
      particlePot.reserve(size);
    }

    if (storageLayout_ & dslDensity) {
      density.reserve(size);
    }

    if (storageLayout_ & dslFunctional) {
      functional.reserve(size);
    }

    if (storageLayout_ & dslFunctionalDerivative) {
      functionalDerivative.reserve(size);
    }

    if (storageLayout_ & dslDipole) {
      dipole.reserve(size);
    }

    if (storageLayout_ & dslQuadrupole) {
      quadrupole.reserve(size);
    }

    if (storageLayout_ & dslElectricField) {
      electricField.reserve(size);
    }

    if (storageLayout_ & dslSkippedCharge) {
      skippedCharge.reserve(size);
    }

    if (storageLayout_ & dslFlucQPosition) {
      flucQPos.reserve(size);
    }

    if (storageLayout_ & dslFlucQVelocity) {
      flucQVel.reserve(size);
    }

    if (storageLayout_ & dslFlucQForce) {
      flucQFrc.reserve(size);
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

    if (storageLayout_ & dslForce) {
      internalCopy(force, source, num, target);
    } 

    if (storageLayout_ & dslAmat) {
      internalCopy(aMat, source, num, target);
    } 

    if (storageLayout_ & dslAngularMomentum) {
      internalCopy(angularMomentum, source, num, target);
    } 

    if (storageLayout_ & dslTorque) {
      internalCopy(torque, source, num, target); 
    }

    if (storageLayout_ & dslParticlePot) {
      internalCopy(particlePot, source, num, target); 
    }

    if (storageLayout_ & dslDensity) {
      internalCopy(density, source, num, target);
    }

    if (storageLayout_ & dslFunctional) {
      internalCopy(functional, source, num, target);
    }

    if (storageLayout_ & dslFunctionalDerivative) {
      internalCopy(functionalDerivative, source, num, target);
    }

    if (storageLayout_ & dslDipole) {
      internalCopy(dipole, source, num, target);
    }

    if (storageLayout_ & dslQuadrupole) {
      internalCopy(quadrupole, source, num, target);
    }

    if (storageLayout_ & dslElectricField) {
      internalCopy(electricField, source, num, target);
    }

    if (storageLayout_ & dslSkippedCharge) {
      internalCopy(skippedCharge, source, num, target);
    }

    if (storageLayout_ & dslFlucQPosition) {
      internalCopy(flucQPos, source, num, target);
    }

    if (storageLayout_ & dslFlucQVelocity) {
      internalCopy(flucQVel, source, num, target);
    }
    if (storageLayout_ & dslFlucQForce) {
      internalCopy(flucQFrc, source, num, target);
    }
  }

  int DataStorage::getStorageLayout() {
    return storageLayout_;
  }

  void DataStorage::setStorageLayout(int layout) {
    storageLayout_ = layout;
    resize(size_);
  }

  RealType* DataStorage::getArrayPointer(int whichArray) {

    switch (whichArray) {
    case dslPosition:
      return internalGetArrayPointer(position);
      break;
            
    case dslVelocity:
      return internalGetArrayPointer(velocity);
      break;
            
    case dslForce:
      return internalGetArrayPointer(force);
      break;            

    case dslAmat:
      return internalGetArrayPointer(aMat);
      break;            
            
    case dslAngularMomentum:
      return internalGetArrayPointer(angularMomentum);
      break;
            
    case dslTorque:
      return internalGetArrayPointer(torque);
      break;

    case dslParticlePot:
      return internalGetArrayPointer(particlePot);
      break;

    case dslDensity:
      return internalGetArrayPointer(density);
      break;

    case dslFunctional:
      return internalGetArrayPointer(functional);
      break;

    case dslFunctionalDerivative:
      return internalGetArrayPointer(functionalDerivative);
      break;

    case dslDipole:
      return internalGetArrayPointer(dipole);
      break;

    case dslQuadrupole:
      return internalGetArrayPointer(quadrupole);
      break;

    case dslElectricField:
      return internalGetArrayPointer(electricField);
      break;

    case dslSkippedCharge:
      return internalGetArrayPointer(skippedCharge);
      break;

    case dslFlucQPosition:
      return internalGetArrayPointer(flucQPos);
      break;

    case dslFlucQVelocity:
      return internalGetArrayPointer(flucQVel);
      break;
           
    case dslFlucQForce:
      return internalGetArrayPointer(flucQFrc);
      break;
           
    default:
      //error message
      return NULL;

    }
  }    

  RealType* DataStorage::internalGetArrayPointer(std::vector<Vector3d>& v) {
    if (v.size() == 0) {
      return NULL;
    } else {
      return v[0].getArrayPointer();
    }
  }

  RealType* DataStorage::internalGetArrayPointer(std::vector<Mat3x3d>& v) {
    if (v.size() == 0) {
      return NULL;
    } else {
      return v[0].getArrayPointer();
    }

  }

  RealType* DataStorage::internalGetArrayPointer(std::vector<RealType>& v) {
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
    }
    if (layout & dslVelocity) {
      bytes += sizeof(Vector3d);
    }
    if (layout & dslForce) {
      bytes += sizeof(Vector3d);
    }
    if (layout & dslAmat) {
      bytes += sizeof(RotMat3x3d);    
    }
    if (layout & dslAngularMomentum) {
      bytes += sizeof(Vector3d);
    }
    if (layout & dslTorque) {
      bytes += sizeof(Vector3d);
    }
    if (layout & dslParticlePot) {
      bytes += sizeof(RealType);
    }
    if (layout & dslDensity) {
      bytes += sizeof(RealType);
    }
    if (layout & dslFunctional) {
      bytes += sizeof(RealType);
    }
    if (layout & dslFunctionalDerivative) {
      bytes += sizeof(RealType);
    }
    if (layout & dslDipole) {
      bytes += sizeof(Vector3d);
    }
    if (layout & dslQuadrupole) {
      bytes += sizeof(Mat3x3d);
    }
    if (layout & dslElectricField) {
      bytes += sizeof(Vector3d);
    }
    if (layout & dslSkippedCharge) {
      bytes += sizeof(RealType);
    }
    if (layout & dslFlucQPosition) {
      bytes += sizeof(RealType);
    }
    if (layout & dslFlucQVelocity) {
      bytes += sizeof(RealType);
    }
    if (layout & dslFlucQForce) {
      bytes += sizeof(RealType);
    }

    return bytes;
  }

}
