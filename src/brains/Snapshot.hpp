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

    enum DataStorageLayout {
        dslPosition = 1,
        dslVelocity = 2,
        dslAMat = 4, 
        dslAngularMomentum = 8,
        dslUnitVector = 16,
        dslZAngle = 32,
        dslForce = 64, 
        dslTorque = 128
    };
    //forward declaration
    class Snapshot;
    class SnapshotManager;
    class StuntDouble;
    
    /**
     * @class DataStorage
     * @brief
     * @warning do not try to insert element into (or ease element from) private member data 
     * of DataStorage directly.
     */
    class DataStorage {
        public:
            DataStorage() {};
            DataStorage(size_t size);

            int size();
            void resize(size_t size);
            void reserve(size_t size);
            
            friend class Snapshot;
            friend class StuntDouble;
        private:
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

            static double* getArrayPointer(vector<Vector3d>& v) {
                assert(v.size() > 0);
                return v[0].getArrayPointer();
            }
            
            static double* getArrayPointer(vector<RotMat3x3d>& v) {
                assert(v.size() > 0);
                return v[0].getArrayPointer();
            }
            
            static double* getArrayPointer(vector<double>& v) {
                assert(v.size() > 0);
                return &(v[0]);
            }

            int getSize() {
                return atomData.size() + rigidbodyData.size();
            }
            
            int getSizeOfAtoms() {
                return atomData.size();
            }
            
            int getSizeOfRigidBodies() {
                return rigidbodyData.size();
            }
            
            DataStorage atomData;
            DataStorage rigidbodyData;

            friend class StuntDouble;
        private:

            int id_; /**< identification number of the snapshot */
    };

    typedef DataStorage (Snapshot::*DataStoragePointer); 
}
#endif //BRAINS_SNAPSHOT_HPP
