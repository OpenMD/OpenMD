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
  * @file Snapshot.hpp
  * @author tlin
  * @date 10/20/2004
  * @time 23:56am
  * @version 1.0
  */
  
#ifndef BRAINS_SNAPSHOT_HPP
#define BRAINS_SNAPSHOT_HPP

#include <vector>

#include "math/Vector3.hpp"
#include "math/SquareMatrix3.hpp"

using namespace std;

namespace oopse{


    //forward declaration
    class Snapshot;
    class SnapshotManager;
    class StuntDouble;
    
    /**
     * @class DataStorage
     * @warning do not try to insert element into (or ease element from) private member data 
     * of DataStorage directly.
     */
    class DataStorage {
        public:

            enum{
                slPosition = 1,
                slVelocity = 2,
                slAmat = 4, 
                slAngularMomentum = 8,
                slUnitVector = 16,
                slZAngle = 32,
                slForce = 64, 
                slTorque = 128
            };
            
            DataStorage(int size);

            /** return the size of this DataStorage */
            int size();
            void resize(int size);
            void reserve(int size);

            void move(int source, int num, int target);
            int getStorageLayout();
            void setStorageLayout(int layout);

            double *getArrayPointer(int );

            friend class StuntDouble;

        private:
            int size_;
            int storageLayout_;
            vector<Vector3d> position;               /** position array */
            vector<Vector3d> velocity;               /** velocity array */
            vector<RotMat3x3d> aMat;            /** rotation matrix array */
            vector<Vector3d> angularMomentum;/** velocity array */
            vector<Vector3d> unitVector;                /** the lab frame unit vector array*/
            vector<double> zAngle;              /** z -angle array */        
            vector<Vector3d> force;               /** force array */
            vector<Vector3d> torque;               /** torque array */
    };

    /**
     * @class Snapshot Snapshot.hpp "brains/Snapshot.hpp"
     * @brief Snapshot class is a repository class for storing dynamic data during 
     *  Simulation
     * Every snapshot class will contain one DataStorage  for atoms and one DataStorage
     *  for rigid bodies.
     * @see SimData
     */
    class Snapshot {
        public:
            
            Snapshot(int nAtoms, int nRigidbodies) {
                atomData.resize(nAtoms);
                rigidbodyData.resize(nRigidbodies);
            }

            Snapshot(const Snapshot& s);

            Snapshot& operator =(const Snapshot& s);
            
            /** Returns the id of this Snapshot */
            int getID() {
                return id_;
            }

            /** Sets the id of this Snapshot */
            void setID(int id) {
                id_ = id;
            }

            //template<typename T>
            //static typename T::ElemPointerType getArrayPointer(vector<T>& v) {
            //    return v[0]->getArrayPointer();
            //}

            int getSize() {
                return atomData.size() + rigidbodyData.size();
            }

            /** Returns the number of atoms */
            int getNumberOfAtoms() {
                return atomData.size();
            }

            /** Returns the number of rigid bodies */
            int getNumberOfRigidBodies() {
                return rigidbodyData.size();
            }

            /** Returns the H-Matrix */
            Mat3x3d getHmat() {
                return hmat_;
            }

            /** Sets the H-Matrix */
            void setHamt(const Mat3x3d& m) {
                hmat_ = m;
                invHmat_ = hmat_.inverse();
            }

            /** Returns the inverse H-Matrix */
            Mat3x3d getInvHmat() {
                return invHmat_
            }

            double getTimeStamp() {
                return timeStamp_;
            }

            void setTimeStamp(double timeStamp) {
                timeStamp_ =timeStamp;
            }


            DataStorage atomData;
            DataStorage rigidbodyData;
            
            friend class StuntDouble;
            
        private:
            double timeStamp_;
            Mat3x3d hmat_;
            Mat3x3d invHmat_;
            int id_; /**< identification number of the snapshot */
    };

    typedef DataStorage (Snapshot::*DataStoragePointer); 
}
#endif //BRAINS_SNAPSHOT_HPP
