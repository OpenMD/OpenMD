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
 * @file DataStorage.cpp
 * @author tlin
 * @date 10/26/2004
 * @version 1.0
 */

#include "brains/DataStorage.hpp"
using namespace std;
namespace OpenMD {

  DataStorage::DataStorage() : size_(0), storageLayout_(0) {}

  DataStorage::DataStorage(std::size_t size, int storageLayout) : size_(size) {
    setStorageLayout(storageLayout);
    resize(size);
  }

  std::size_t DataStorage::getSize() {
    if (storageLayout_ & dslPosition && position.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslVelocity && velocity.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslForce && force.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslAmat && aMat.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslAngularMomentum &&
        angularMomentum.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslTorque && torque.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslParticlePot && particlePot.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslDensity && density.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslFunctional && functional.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslFunctionalDerivative &&
        functionalDerivative.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslDipole && dipole.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslQuadrupole && quadrupole.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslElectricField && electricField.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslSkippedCharge && skippedCharge.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslFlucQPosition && flucQPos.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslFlucQVelocity && flucQVel.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslFlucQForce && flucQFrc.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    if (storageLayout_ & dslSitePotential && sitePotential.size() != size_) {
      // error
      cerr << "size does not match" << endl;
    }

    return size_;
  }

  void DataStorage::resize(std::size_t newSize) {
    if (storageLayout_ & dslPosition) { internalResize(position, newSize); }

    if (storageLayout_ & dslVelocity) { internalResize(velocity, newSize); }

    if (storageLayout_ & dslForce) { internalResize(force, newSize); }

    if (storageLayout_ & dslAmat) { internalResize(aMat, newSize); }

    if (storageLayout_ & dslAngularMomentum) {
      internalResize(angularMomentum, newSize);
    }

    if (storageLayout_ & dslTorque) { internalResize(torque, newSize); }

    if (storageLayout_ & dslParticlePot) {
      internalResize(particlePot, newSize);
    }

    if (storageLayout_ & dslDensity) { internalResize(density, newSize); }

    if (storageLayout_ & dslFunctional) { internalResize(functional, newSize); }

    if (storageLayout_ & dslFunctionalDerivative) {
      internalResize(functionalDerivative, newSize);
    }

    if (storageLayout_ & dslDipole) { internalResize(dipole, newSize); }

    if (storageLayout_ & dslQuadrupole) { internalResize(quadrupole, newSize); }

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

    if (storageLayout_ & dslFlucQForce) { internalResize(flucQFrc, newSize); }

    if (storageLayout_ & dslSitePotential) {
      internalResize(sitePotential, newSize);
    }

    size_ = newSize;
  }

  void DataStorage::reserve(std::size_t size) {
    if (storageLayout_ & dslPosition) { position.reserve(size); }

    if (storageLayout_ & dslVelocity) { velocity.reserve(size); }

    if (storageLayout_ & dslForce) { force.reserve(size); }

    if (storageLayout_ & dslAmat) { aMat.reserve(size); }

    if (storageLayout_ & dslAngularMomentum) { angularMomentum.reserve(size); }

    if (storageLayout_ & dslTorque) { torque.reserve(size); }

    if (storageLayout_ & dslParticlePot) { particlePot.reserve(size); }

    if (storageLayout_ & dslDensity) { density.reserve(size); }

    if (storageLayout_ & dslFunctional) { functional.reserve(size); }

    if (storageLayout_ & dslFunctionalDerivative) {
      functionalDerivative.reserve(size);
    }

    if (storageLayout_ & dslDipole) { dipole.reserve(size); }

    if (storageLayout_ & dslQuadrupole) { quadrupole.reserve(size); }

    if (storageLayout_ & dslElectricField) { electricField.reserve(size); }

    if (storageLayout_ & dslSkippedCharge) { skippedCharge.reserve(size); }

    if (storageLayout_ & dslFlucQPosition) { flucQPos.reserve(size); }

    if (storageLayout_ & dslFlucQVelocity) { flucQVel.reserve(size); }

    if (storageLayout_ & dslFlucQForce) { flucQFrc.reserve(size); }

    if (storageLayout_ & dslSitePotential) { sitePotential.reserve(size); }
  }

  void DataStorage::copy(int source, std::size_t num, std::size_t target) {
    if (num + target > size_) {
      // error
    }

    if (storageLayout_ & dslPosition) {
      internalCopy(position, source, num, target);
    }

    if (storageLayout_ & dslVelocity) {
      internalCopy(velocity, source, num, target);
    }

    if (storageLayout_ & dslForce) { internalCopy(force, source, num, target); }

    if (storageLayout_ & dslAmat) { internalCopy(aMat, source, num, target); }

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

    if (storageLayout_ & dslSitePotential) {
      internalCopy(sitePotential, source, num, target);
    }
  }

  int DataStorage::getStorageLayout() { return storageLayout_; }

  void DataStorage::setStorageLayout(int layout) {
    storageLayout_ = layout;
    resize(size_);
  }

  RealType* DataStorage::getArrayPointer(int whichArray) {
    switch (whichArray) {
    case dslPosition:
      return internalGetArrayPointer(position);

    case dslVelocity:
      return internalGetArrayPointer(velocity);

    case dslForce:
      return internalGetArrayPointer(force);

    case dslAmat:
      return internalGetArrayPointer(aMat);

    case dslAngularMomentum:
      return internalGetArrayPointer(angularMomentum);

    case dslTorque:
      return internalGetArrayPointer(torque);

    case dslParticlePot:
      return internalGetArrayPointer(particlePot);

    case dslDensity:
      return internalGetArrayPointer(density);

    case dslFunctional:
      return internalGetArrayPointer(functional);

    case dslFunctionalDerivative:
      return internalGetArrayPointer(functionalDerivative);

    case dslDipole:
      return internalGetArrayPointer(dipole);

    case dslQuadrupole:
      return internalGetArrayPointer(quadrupole);

    case dslElectricField:
      return internalGetArrayPointer(electricField);

    case dslSkippedCharge:
      return internalGetArrayPointer(skippedCharge);

    case dslFlucQPosition:
      return internalGetArrayPointer(flucQPos);

    case dslFlucQVelocity:
      return internalGetArrayPointer(flucQVel);

    case dslFlucQForce:
      return internalGetArrayPointer(flucQFrc);

    case dslSitePotential:
      return internalGetArrayPointer(sitePotential);

    default:
      // error message
      return NULL;
    }
  }

  RealType* DataStorage::internalGetArrayPointer(std::vector<Vector3d>& v) {
    if (v.empty()) {
      return NULL;
    } else {
      return v[0].getArrayPointer();
    }
  }

  RealType* DataStorage::internalGetArrayPointer(std::vector<Mat3x3d>& v) {
    if (v.empty()) {
      return NULL;
    } else {
      return v[0].getArrayPointer();
    }
  }

  RealType* DataStorage::internalGetArrayPointer(std::vector<RealType>& v) {
    if (v.empty()) {
      return NULL;
    } else {
      return &(v[0]);
    }
  }

  template<typename T>
  void DataStorage::internalResize(std::vector<T>& v, std::size_t newSize) {
    std::size_t oldSize = v.size();

    if (oldSize == newSize) {
      return;
    } else if (oldSize < newSize) {
      v.insert(v.end(), newSize - oldSize, T());
    } else {
      typename std::vector<T>::iterator i;
      i = v.begin();
      std::advance(i, newSize);
      v.erase(i, v.end());
    }
  }

  template<typename T>
  void DataStorage::internalCopy(std::vector<T>& v, int source, std::size_t num,
                                 std::size_t target) {
    typename std::vector<T>::iterator first;
    typename std::vector<T>::iterator last;
    typename std::vector<T>::iterator result;

    first  = v.begin();
    last   = v.begin();
    result = v.begin();

    std::advance(first, source);
    // STL algorithm use half opened range
    std::advance(last, num + 1);
    std::advance(result, target);

    std::copy(first, last, result);
  }

  std::size_t DataStorage::getBytesPerStuntDouble(int layout) {
    std::size_t bytes = 0;
    if (layout & dslPosition) { bytes += sizeof(Vector3d); }
    if (layout & dslVelocity) { bytes += sizeof(Vector3d); }
    if (layout & dslForce) { bytes += sizeof(Vector3d); }
    if (layout & dslAmat) { bytes += sizeof(RotMat3x3d); }
    if (layout & dslAngularMomentum) { bytes += sizeof(Vector3d); }
    if (layout & dslTorque) { bytes += sizeof(Vector3d); }
    if (layout & dslParticlePot) { bytes += sizeof(RealType); }
    if (layout & dslDensity) { bytes += sizeof(RealType); }
    if (layout & dslFunctional) { bytes += sizeof(RealType); }
    if (layout & dslFunctionalDerivative) { bytes += sizeof(RealType); }
    if (layout & dslDipole) { bytes += sizeof(Vector3d); }
    if (layout & dslQuadrupole) { bytes += sizeof(Mat3x3d); }
    if (layout & dslElectricField) { bytes += sizeof(Vector3d); }
    if (layout & dslSkippedCharge) { bytes += sizeof(RealType); }
    if (layout & dslFlucQPosition) { bytes += sizeof(RealType); }
    if (layout & dslFlucQVelocity) { bytes += sizeof(RealType); }
    if (layout & dslFlucQForce) { bytes += sizeof(RealType); }
    if (layout & dslSitePotential) { bytes += sizeof(RealType); }

    return bytes;
  }

}  // namespace OpenMD
