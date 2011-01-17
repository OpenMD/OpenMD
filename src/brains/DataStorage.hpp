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
 * [4]  Vardeman & Gezelter, in progress (2009).                        
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
#include <math/Vector3.hpp>
#include <math/SquareMatrix3.hpp>

namespace OpenMD {
  /**
   * @class DataStorage
   * @warning do not try to insert element into (or ease element from) private member data 
   * of DataStorage directly.
   * @todo DataStorage may need refactoring. Every std::vector can inherit from the same base class
   * which will make it easy to maintain
   */
  class DataStorage {
  public:

    enum{
      dslPosition = 1,
      dslVelocity = 2,
      dslAmat = 4, 
      dslAngularMomentum = 8,
      dslElectroFrame = 16,
      dslZAngle = 32,
      dslForce = 64, 
      dslTorque = 128,
      dslParticlePot = 256,
      dslDensity = 512,
      dslFunctional = 1024,
      dslFunctionalDerivative = 2048,
      dslElectricField = 4096
    };


    DataStorage();
    DataStorage(int size, int storageLayout = 8191);
    /** return the size of this DataStorage. */
    int getSize();
    /**
     * Changes the size of this DataStorage.
     * @param size new size of this DataStorage
     */
    void resize(int newSize);
    /**
     * Reallocates memory manually. 
     *
     * The main reason for using reserve() is efficiency if you know
     * the capacity to which your std::vector must eventually grow,
     * then it is usually more efficient to allocate that memory all
     * at once.
     */
    void reserve(int size);
    /**
     * Copies data inside DataStorage class.
     *
     * Copy function actually call std::copy for every std::vector in
     * DataStorage class.  One Precondition of std::copy is that
     * target is not within the range [soruce, soruce + num]
     *
     * @param souce 
     * @param num number of element to be moved
     * @param target
     */
    void copy(int source, int num, int target);
    /** Returns the storage layout  */
    int getStorageLayout();
    /** Sets the storage layout  */
    void setStorageLayout(int layout);
    /** Returns the pointer of internal array */
    RealType *getArrayPointer(int whichArray);

    std::vector<Vector3d> position;        /** position array */
    std::vector<Vector3d> velocity;        /** velocity array */
    std::vector<RotMat3x3d> aMat;          /** rotation matrix array */
    std::vector<Vector3d> angularMomentum; /** angular momentum array
                                               (body-fixed) */
    std::vector<Mat3x3d> electroFrame;     /** the lab frame unit vector array*/
    std::vector<RealType> zAngle;          /** z-angle array */        
    std::vector<Vector3d> force;           /** force array */
    std::vector<Vector3d> torque;          /** torque array */
    std::vector<RealType> particlePot;     /** particle potential arrray */
    std::vector<RealType> density;         /** electron density */
    std::vector<RealType> functional;      /** density functional */
    std::vector<RealType> functionalDerivative; /** derivative of functional */
    std::vector<Vector3d> electricField;   /** local electric field */

    static int getBytesPerStuntDouble(int layout);

  private:

    RealType* internalGetArrayPointer(std::vector<Vector3d>& v);            
    RealType* internalGetArrayPointer(std::vector<RotMat3x3d>& v);
    RealType* internalGetArrayPointer(std::vector<RealType>& v);

            
    template<typename T>
    void internalResize(std::vector<T>& v, int newSize);

    template<typename T>
    void internalCopy(std::vector<T>& v, int source,  int num, int target);
            
    int size_;
    int storageLayout_;

  };

}
#endif //BRAINS_DATASTORAGE_HPP
