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

using namespace oopse;



    //forward declaration
    class Snapshot;
    class StuntDouble;
    class DataStorageTestCase;
    /**
     * @class DataStorage
     * @warning do not try to insert element into (or ease element from) private member data 
     * of DataStorage directly.
     * @todo DataStorage may need refactorying. Every std::vector can inherit from the same base class
     * which will make it easy to maintain
     */
    class DataStorage {
        public:

            enum{
                dslPosition = 1,
                dslVelocity = 2,
                dslAmat = 4, 
                dslAngularMomentum = 8,
                dslUnitVector = 16,
                dslZAngle = 32,
                dslForce = 64, 
                dslTorque = 128
            };


            DataStorage();
            DataStorage(int size, int storageLayout = 0x11111111);
            /** return the size of this DataStorage. */
            int getSize();
            /**
             * Changes the size of this DataStorage.
             * @param size new size of this DataStorage
             */
            void resize(int newSize);
            /**
             * Reallocates memory manually. The main reason for using reserve() is efficiency
             * if you know the capacity to which your std::vector must eventually grow, then it is usually more
             * efficient to allocate that memory all at once.
             */
            void reserve(int size);
            /**
             * Copies data inside DataStorage class.
             * Copy function actually call std::copy for every std::vector in DataStorage class. 
             * One Precondition of std::copy is that target is not within the range [soruce, soruce + num]
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
            double *getArrayPointer(int whichArray);
            friend class StuntDouble;
            friend class DataStorageTestCase;
        private:


            double* internalGetArrayPointer(std::vector<Vector3d>& v);
            
            double* internalGetArrayPointer(std::vector<RotMat3x3d>& v);
            double* internalGetArrayPointer(std::vector<double>& v);
            
            template<typename T>
            void internalResize(std::vector<T>& v, int newSize);

            template<typename T>
            void interalCopy(std::vector<T>& v, int source,  int num, int target);
            
            int size_;
            int storageLayout_;
            std::vector<Vector3d> position;               /** position array */
            std::vector<Vector3d> velocity;               /** velocity array */
            std::vector<RotMat3x3d> aMat;            /** rotation matrix array */
            std::vector<Vector3d> angularMomentum;/** velocity array */
            std::vector<Vector3d> unitVector;                /** the lab frame unit std::vector array*/
            std::vector<double> zAngle;              /** z -angle array */        
            std::vector<Vector3d> force;               /** force array */
            std::vector<Vector3d> torque;               /** torque array */
    };


#endif //BRAINS_DATASTORAGE_HPP
