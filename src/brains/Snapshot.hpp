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

#include "brains/DataStorage.hpp"

using namespace std;

namespace oopse{

    class StuntDouble;
    
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
                return atomData.getSize() + rigidbodyData.getSize();
            }

            /** Returns the number of atoms */
            int getNumberOfAtoms() {
                return atomData.getSize();
            }

            /** Returns the number of rigid bodies */
            int getNumberOfRigidBodies() {
                return rigidbodyData.getSize();
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
