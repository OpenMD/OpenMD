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
 #ifndef _STUNTDOUBLE_HPP_
 #define _STUNTDOUBLE_HPP_

#include <vector>

#include <core/BaseVisitor.hpp>
#include <math/Quaternion.hpp>
#include <math/Mat3x3d.hpp>
#include <math/Vector3d.hpp>
#include <util/PropertyMap.hpp>

namespace oopse{

  /**
   * The base class for atom types. Typically, atom types are used to descibe the behavior of an atom
   * of a element in differnet enviroment. In OOPSE, atom types are also used to represent the coarse-
   * grainded atom
   * @author    tlin
   * @date  09/08/2004
   * @version 1.0
   */

  /*
   * Design Decision:
   * the member data of stuntdouble will be stored outside of the class
   */
   class StuntDouble{
    public:    

        /**
         * Returns the global index of this stuntdouble.
         * @return  the global index of this stuntdouble 
         */
        int getGlobalIndex() {
            return globalIndex_;
        }

        /**
         * Sets the global index of this stuntdouble.
         * @param new global index to be set
         */
        void setGlobalIndex(int index) {
            return globalIndex_;
        }
        
        /** 
         * Returns the local index of this stuntdouble 
         * @return the local index of this stuntdouble
         */
        int getLocalIndex() {
            return localIndex_;
        }

        /**
         * Sets the local index of this stuntdouble
         * @param index new index to be set
         */        
        void setLocalIndex(int index) {
            localIndex_ = index;
        }

       /**
         * Returns  the previous position of this stuntdouble
         * @return the position of this stuntdouble
         */    
        void Vector3d getPrevPos() {
            return (snapshotMan_->getPrevSnapshot())->storage->position[localId_];
        }
       
        /**
         * Returns  the current position of this stuntdouble
         * @return the position of this stuntdouble
         */    
        void Vector3d getPos() {
            return (snapshotMan_->getCurrentSnapshot())->storage->position[localId_];
        }

       /**
         * Returns  the position of this stuntdouble in specified snapshot 
         *
         * @return the position of this stuntdouble
         * @param snapshotNo
         */    
        Vector3d getPos(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage->position[localId_];
        }

       /**
         * Sets  the previous position of this stuntdouble
         *
         * @param pos  new position 
         * @see #getPos
         */         
       void setPrevPos(const Vector3d& pos) {
            (snapshotMan_->getPrevSnapshot())->storage->position[localId_] = pos;
       }
       
       /**
         * Sets  the current position of this stuntdouble
         * @param pos  new position 
         */         
        void setPos(const Vector3d& pos) {
            (snapshotMan_->getCurrentSnapshot())->storage->position[localId_] = pos;
        }

       /**
         * Sets  the position of this stuntdouble in specified snapshot
         *
         * @param pos position to be set 
         * @param snapshotNo 
         * @see #getPos
         */         
        void setPos(const Vector3d& pos, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage->position[localId_] = pos;
        }
       
       /**
         * Returns  the previous velocity of this stuntdouble
         * @return the velocity of this stuntdouble
         */    
        Vector3d getPrevVel() {
            return (snapshotMan_->getPrevSnapshot())->storage->velocity[localId_];
        }
       
        /**
         * Returns  the current velocity of this stuntdouble
         * @return the velocity of this stuntdouble
         */    
        Vector3d getVel() {
            return (snapshotMan_->getCurrentSnapshot())->storage->velocity[localId_];
        }

       /**
         * Returns  the velocity of this stuntdouble in specified snapshot 
         *
         * @return the velocity of this stuntdouble
         * @param snapshotNo
         */    
         Vector3d getVel(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage->velocity[localId_];
        }

       /**
         * Sets  the previous velocity of this stuntdouble
         *
         * @param vel  new velocity 
         * @see #getVel
         */         
       void setPrevVel(const Vector3d& vel) {
            (snapshotMan_->getPrevSnapshot())->storage->velocity[localId_] = vel;
       }
       
       /**
         * Sets  the current velocity of this stuntdouble
         * @param vel  new velocity 
         */         
        void setVel(const Vector3d& vel) {
            (snapshotMan_->getCurrentSnapshot())->storage->velocity[localId_] = vel;
        }

       /**
         * Sets  the velocity of this stuntdouble in specified snapshot
         *
         * @param vel velocity to be set 
         * @param snapshotNo 
         * @see #getVel
         */         
        void setVel(const Vector3d& vel, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage->velocity[localId_] = vel;
        }

       /**
         * Returns  the previous rotation matrix of this stuntdouble
         * @return the rotation matrix of this stuntdouble
         */    
        RotMat3x3d getPrevA() {
            return (snapshotMan_->getPrevSnapshot())->storage->aMat[localId_];
        }
       
        /**
         * Returns  the current rotation matrix of this stuntdouble
         * @return the rotation matrix of this stuntdouble
         */    
        RotMat3x3d getA() {
            return (snapshotMan_->getCurrentSnapshot())->storage->aMat[localId_];
        }

       /**
         * Returns  the rotation matrix of this stuntdouble in specified snapshot 
         *
         * @return the rotation matrix of this stuntdouble
         * @param snapshotNo
         */    
         RotMat3x3d getA(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage->aMat[localId_];
        }

       /**
         * Sets  the previous rotation matrix of this stuntdouble
         *
         * @param a  new rotation matrix 
         * @see #getA
         */         
       void setPrevA(const RotMat3x3d& a) {
            (snapshotMan_->getPrevSnapshot())->storage->aMat[localId_] = a;
       }
       
       /**
         * Sets  the current rotation matrix of this stuntdouble
         * @param a  new rotation matrix 
         */         
        void setA(const RotMat3x3d& a) {
            (snapshotMan_->getCurrentSnapshot())->storage->aMat[localId_] = a;
        }

       /**
         * Sets  the rotation matrix of this stuntdouble in specified snapshot
         *
         * @param a rotation matrix to be set 
         * @param snapshotNo 
         * @see #getA
         */         
        void setA(const RotMat3x3d& a, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage->a[localId_] = a;
        }       

       /**
         * Returns  the previous angular momentum of this stuntdouble
         * @return the angular momentum of this stuntdouble
         */    
        Vector3d getPrevJ() {
            return (snapshotMan_->getPrevSnapshot())->storage->angularMomentum[localId_];
        }
       
        /**
         * Returns  the current angular momentum of this stuntdouble
         * @return the angular momentum of this stuntdouble
         */    
        Vector3d getJ() {
            return (snapshotMan_->getCurrentSnapshot())->storage->angularMomentum[localId_];
        }

       /**
         * Returns  the angular momentum of this stuntdouble in specified snapshot 
         *
         * @return the angular momentum of this stuntdouble
         * @param snapshotNo
         */    
         Vector3d getJ(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage->angularMomentum[localId_];
        }

       /**
         * Sets  the previous angular momentum of this stuntdouble
         *
         * @param angMom  new angular momentum 
         * @see #getJ
         */         
       void setPrevJ(const Vector3d& angMom) {
            (snapshotMan_->getPrevSnapshot())->storage->angularMomentum[localId_] = angMom;
       }
       
       /**
         * Sets  the current angular momentum of this stuntdouble
         * @param angMom  new angular momentum 
         */         
        void setJ(const Vector3d& angMom) {
            (snapshotMan_->getCurrentSnapshot())->storage->angularMomentum[localId_] = angMom;
        }

       /**
         * Sets the angular momentum of this stuntdouble in specified snapshot
         *
         * @param angMom angular momentum to be set 
         * @param snapshotNo 
         * @see #getJ
         */         
        void setJ(const Vector3d& angMom, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage->angularMomentum[localId_] = angMom;
        }
         
       /**
         * Returns  the previous quaternion of this stuntdouble
         * @return the quaternion of this stuntdouble
         */    
        Quat4d getPrevQ() {
            return (snapshotMan_->getPrevSnapshot())->storage->aMat[localId_].toQuaternion();
        }
       
        /**
         * Returns  the current quaternion of this stuntdouble
         * @return the quaternion of this stuntdouble
         */    
        Quat4d getQ() {
            return (snapshotMan_->getCurrentSnapshot())->storage->aMat[localId_].toQuaternion();
        }

       /**
         * Returns  the quaternion of this stuntdouble in specified snapshot 
         *
         * @return the quaternion of this stuntdouble
         * @param snapshotNo
         */    
         Quat4d getQ(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage->aMat[localId_].toQuaternion();
        }

       /**
         * Sets  the previous quaternion of this stuntdouble
         *
         * @param q  new quaternion 
         * @see #getQ
         */         
       void setPrevQ(const Quat4d& q) {
            (snapshotMan_->getPrevSnapshot())->storage->aMat[localId_] = q;
       }
       
       /**
         * Sets  the current quaternion of this stuntdouble
         * @param q  new quaternion 
         */         
        void setQ(const Quat4d& q) {
            (snapshotMan_->getCurrentSnapshot())->storage->aMat[localId_] = q;
        }

       /**
         * Sets  the quaternion of this stuntdouble in specified snapshot
         *
         * @param q quaternion to be set 
         * @param snapshotNo 
         * @see #getQ
         */         
        void setQ(const Quat4d& q, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage->aMat[localId_] = q;
        }
       
       /**
         * Returns the force of this stuntdouble
         *
         * @return the quaternion of this stuntdouble
         *
         * @see #setFrc
         * @see #addFrc
         */
        virtual Vector3d getFrc() = 0;

       /**
         * Sets the force of this stuntdouble
         *
         * @param frc new force
         *
         * @see #getFrc
         * @see #addFrc
         */     
        virtual void setFrc(Vector3d& frc) = 0;

       /**
         * Adds the force into this stuntdouble
         *
         * @param frc  force to be added
         *
         * @see #getFrc
         * @see #setFrc
         */          
        virtual void addFrc(Vector3d& frc) = 0;

       /**
         * Returns the torque of this stuntdouble
         *
         * @return the torque of this stuntdouble
         *
         * @see #setTrq
         * @see #addTrq
         */
        virtual Vector3d getTrq() = 0;

       /**
         * Sets the torque of this stuntdouble
         *
         * @param trq new torque
         *
         * @see #getTrq
         * @see #addTrq
         */      
        virtual void setTrq(Vector3d& trq) = 0;

       /**
         * Adds the torque into this stuntdouble
         *
         * @param trq  torque to be added
         *
         * @see #getTrq
         * @see #setTrq
         */
        virtual void addTrq(Vector3d& trq) = 0;

       /**
         * Returns the inertia tensor of this stuntdouble
         *
         * @return the inertia tensor of this stuntdouble
         *
         * @see #setI
         */ 
        virtual Mat3x3d getI() = 0;

       /**
         * Sets the inertia tensor of this stuntdouble
         *
         * @param trq new inertia tensor
         *
         * @see #getI
         */      
        virtual void setI(Mat3x3d& I) = 0;

       /**
         * Returns the gradient of this stuntdouble
         *
         * @return the inertia tensor of this stuntdouble
         *
         * @see #setI
         */ 
        virtual std::vector<double> getGrad() = 0;

       /**
         * Returns the euler angles of this stuntdouble
         * <p>
         *   <ol>
         *   <il>
         *   <il>
         *   <il>
         *   </ol>
         * </p>
         * @return the euler angles of this stuntdouble
         *
         * @see #setEuler
         */ 
        virtual Vector3d getEuler() = 0;

       /**
         * Sets the euler angles of this stuntdouble
         *
         * @param e new euler angles
         * <p>
         *   <ol>
         *   <il>e[0] = phi </il>
         *   <il>e[1] = theta <il>
         *   <il>e[2] = psi </il>
         *   </ol>
         * </p>
         *
         * @see #setEuler
         */ 
        virtual void setEuler(Vector3d& e) = 0;
        
        virtual bool isLinear() {return false;}

       /**
         * Returns the linear axis of the rigidbody, atom and directional atom will always return -1
         *
         * @return the linear axis of the rigidbody
         * 
         * @see #isLinear
         */ 
        virtual int linearAxis() {return -1;}

       /**
        * Returns the zangle
        *
        * @return the zangle
        *
        * @see #setZangle
        * @see #addZangle
        */
        virtual double   getZangle() = 0;

       /**
        *
        */
        virtual void   setZangle(double zAngle) = 0;

       /**
        *
        */
        virtual void   addZangle(double zAngle) = 0;

       /**
         * 
         * <p>
         * The purpose of the Visitor Pattern is to encapsulate an operation that you want to perform on
         * the elements of a data structure. In this way, you can change the operation being performed 
         * on a structure without the need of changing the classes of the elements that you are operating
         * on. Using a Visitor pattern allows you to decouple the classes for the data structure and the 
         * algorithms used upon them
         * </p>
         * @param v visitor
         */      
        virtual void accept(BaseVisitor* v) = 0;

        StorageMethod getStorageMethd(){return stMethod;}
        void setStorageMethod(StorageMethod method) {stMethod = method;}

        //below functions are just forward functions
        /**
         * Adds property into property map
         *
         * @param genData GenericData to be added into PropertyMap
         *
         * @see #removeProperty
         * @see #clearProperties
         */
        void addProperty(GenericData* genData) {  properties.addProperty(genData);  }

        /**
         * Removes property from PropertyMap by name
         *
         * @param propName the name of property to be removed
         *
         * @see #addProperty
         * @see #clearProperties
         */
        void removeProperty(std::string& propName) {  properties.removeProperty();  }

        /**
         * clear all of the properties
         *
         * @see #addProperty
         * @see #removeProperty
         */
        void clearProperties() {  properties.clearProperties(); }

        /**
         * Returns all names of properties
         *
         * @return all names of properties
         */
        std::vector<std::string> getPropertyNames() {return properties.getPropertyNames();  }

        /**
         * Returns all of the properties in PropertyMap
         *
         * @return all of the properties in PropertyMap
         *
         * @see #getPropertyByName
         */      
        std::vector<GenericData*> getProperties() { return properties.getProperties(); }

        /**
         * Returns property 
         *
         * @param propName name of property
         *
         * @return a pointer point to property with propName. If no property named propName
         * exists, return NULL
         *
         * @see #getProperties
         */      
        GenericData* getPropertyByName(std:string& propName) {  return properties.getPropertyByName(propName); }

      private:
        
        int globalIndex_;
        int localIndex_;
        PropertyMap properties;
        DataStoragePointer storage_;
        SnapshotManager* snapshotMan_;
 };

}//end namespace oopse
#endif //ifndef _STUNTDOUBLE_HPP_
