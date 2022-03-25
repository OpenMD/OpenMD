/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

/**
 * @file Vector.hpp
 * @author Teng Lin
 * @date 09/14/2004
 * @version 1.0
 */

#ifndef BRAINS_DATASTORAGE_HPP
#define BRAINS_DATASTORAGE_HPP

#include <vector>

#include <math/SquareMatrix3.hpp>
#include <math/Vector3.hpp>

using namespace std;
namespace OpenMD {
  /**
   * @class DataStorage
   * @warning do not try to insert element into (or ease element from) private
   * member data of DataStorage directly.
   * @todo DataStorage may need refactoring. Every vector can inherit from the
   * same base class which will make it easy to maintain
   */
  class DataStorage {
  public:
    enum {
      dslPosition             = 1,
      dslVelocity             = 2,
      dslForce                = 4,
      dslAmat                 = 8,
      dslAngularMomentum      = 16,
      dslTorque               = 32,
      dslParticlePot          = 64,
      dslDensity              = 128,
      dslFunctional           = 256,
      dslFunctionalDerivative = 512,
      dslDipole               = 1024,
      dslQuadrupole           = 2048,
      dslElectricField        = 4096,
      dslSkippedCharge        = 8192,
      dslFlucQPosition        = 16384,
      dslFlucQVelocity        = 32768,
      dslFlucQForce           = 65536,
      dslSitePotential        = 131072
    };

    DataStorage();
    DataStorage(std::size_t size, int storageLayout = 0);
    /** return the size of this DataStorage. */
    std::size_t getSize();
    /**
     * Changes the size of this DataStorage.
     * @param newSize new size of this DataStorage
     */
    void resize(std::size_t newSize);
    /**
     * Reallocates memory manually.
     *
     * The main reason for using reserve() is efficiency if you know
     * the capacity to which your vector must eventually grow,
     * then it is usually more efficient to allocate that memory all
     * at once.
     */
    void reserve(std::size_t size);
    /**
     * Copies data inside DataStorage class.
     *
     * Copy function actually calls copy for every vector in
     * DataStorage class.  One Precondition of copy is that
     * target is not within the range [source, soruce + num]
     *
     * @param source
     * @param num number of element to be moved
     * @param target
     */
    void copy(int source, std::size_t num, std::size_t target);
    /** Returns the storage layout  */
    int getStorageLayout();
    /** Sets the storage layout  */
    void setStorageLayout(int layout);
    /** Returns the pointer of internal array */
    RealType* getArrayPointer(int whichArray);

    vector<Vector3d> position;        /** position array */
    vector<Vector3d> velocity;        /** velocity array */
    vector<Vector3d> force;           /** force array */
    vector<RotMat3x3d> aMat;          /** rotation matrix array */
    vector<Vector3d> angularMomentum; /** angular momentum array (body-fixed) */
    vector<Vector3d> torque;          /** torque array */
    vector<RealType> particlePot;     /** particle potential arrray */
    vector<RealType> density;         /** electron density */
    vector<RealType> functional;      /** density functional */
    vector<RealType> functionalDerivative; /** derivative of functional */
    vector<Vector3d> dipole;               /** space-frame dipole vector */
    vector<Mat3x3d> quadrupole;            /** space-frame quadrupole tensor */
    vector<Vector3d> electricField;        /** local electric field */
    vector<RealType>
        skippedCharge; /** charge skipped during normal pairwise calculation */
    vector<RealType> flucQPos;      /** fluctuating charges */
    vector<RealType> flucQVel;      /** fluctuating charge velocities */
    vector<RealType> flucQFrc;      /** fluctuating charge forces */
    vector<RealType> sitePotential; /** electrostatic site potentials */

    static std::size_t getBytesPerStuntDouble(int layout);

  private:
    RealType* internalGetArrayPointer(vector<Vector3d>& v);
    RealType* internalGetArrayPointer(vector<Mat3x3d>& v);
    RealType* internalGetArrayPointer(vector<RealType>& v);

    template<typename T>
    void internalResize(std::vector<T>& v, std::size_t newSize);

    template<typename T>
    void internalCopy(std::vector<T>& v, int source, std::size_t num,
                      std::size_t target);

    std::size_t size_;
    int storageLayout_;
  };
}  // namespace OpenMD

#endif  // BRAINS_DATASTORAGE_HPP
