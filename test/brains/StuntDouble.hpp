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
 * @file StuntDouble.hpp
 * @author    tlin
 * @date  10/22/2004
 * @version 1.0
 */ 
   
 #ifndef _STUNTDOUBLE_HPP_
 #define _STUNTDOUBLE_HPP_

#include <vector>

#include <core/BaseVisitor.hpp>
#include <math/Quaternion.hpp>
#include <math/SquareMatrix3.hpp>
#include <math/Vector3.hpp>
#include <utils/PropertyMap.hpp>

namespace oopse{


   
  /**
   * @class StuntDouble StuntDouble.hpp "Primitives/StuntDouble.hpp"
   * @brief 
   * StuntDouble is a very strange idea.  A StuntDouble stands in for
   * some object that can be manipulated by the Integrators or
   * Minimizers.  Some of the manipulable objects are Atoms, some are
   * DirectionalAtoms, and some are RigidBodies.  StuntDouble
   * provides an interface for the Integrators and Minimizers to use,
   * and does some preliminary sanity checking so that the program
   * doesn't try to do something stupid like torque an Atom
   * @note the dynamoc data of stuntdouble will be stored outside of the class
   */
   class StuntDouble{
   public:    

        enum ObjectType{
            otAtom,
            otDAtom,
            otRigidBody
        };

        virtual ~StuntDouble();
        
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
         * Tests if this stuntdouble is an atom 
         * @return true is this stuntdouble is an atom(or a directional atom), return false otherwise
         */
        bool isAtom(){
            return objType_ == otAtom || objType_ == otDAtom;
        }

        /** 
         * Tests if this stuntdouble is an directional atom 
         * @return true if this stuntdouble is an directional atom, return false otherwise
         */
        bool isDirectionalAtom(){
            return objType_ == otDAtom;
        }

        /**
         * Tests if this stuntdouble is a rigid body. 
         * @return true if this stuntdouble is a rigid body, otherwise return false
         */
        bool isRigidBody(){
            return objType_ == otRigidBody;
        }

        /**
         * Tests if this stuntdouble is a directional one. 
         * @return true is this stuntdouble is a directional atom or a rigid body, return false otherwise
         */
        bool isDirectional(){
            return isDirectionalAtom() || isRigidBody();
        }

       /**
         * Returns the previous position of this stuntdouble
         * @return the position of this stuntdouble
         */    
        void Vector3d getPrevPos() {
            return (snapshotMan_->getPrevSnapshot())->storage_->position[localIndex_];
        }
       
        /**
         * Returns the current position of this stuntdouble
         * @return the position of this stuntdouble
         */    
        void Vector3d getPos() {
            return (snapshotMan_->getCurrentSnapshot())->storage_->position[localIndex_];
        }

       /**
         * Returns the position of this stuntdouble in specified snapshot 
         * @return the position of this stuntdouble
         * @param snapshotNo
         */    
        Vector3d getPos(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage_->position[localIndex_];
        }

       /**
         * Sets  the previous position of this stuntdouble
         * @param pos  new position 
         * @see #getPos
         */         
       void setPrevPos(const Vector3d& pos) {
            (snapshotMan_->getPrevSnapshot())->storage_->position[localIndex_] = pos;
       }
       
       /**
         * Sets  the current position of this stuntdouble
         * @param pos  new position 
         */         
        void setPos(const Vector3d& pos) {
            (snapshotMan_->getCurrentSnapshot())->storage_->position[localIndex_] = pos;
        }

       /**
         * Sets  the position of this stuntdouble in specified snapshot
         * @param pos position to be set 
         * @param snapshotNo 
         * @see #getPos
         */         
        void setPos(const Vector3d& pos, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->position[localIndex_] = pos;
        }
       
       /**
         * Returns the previous velocity of this stuntdouble
         * @return the velocity of this stuntdouble
         */    
        Vector3d getPrevVel() {
            return (snapshotMan_->getPrevSnapshot())->storage_->velocity[localIndex_];
        }
       
        /**
         * Returns the current velocity of this stuntdouble
         * @return the velocity of this stuntdouble
         */    
        Vector3d getVel() {
            return (snapshotMan_->getCurrentSnapshot())->storage_->velocity[localIndex_];
        }

       /**
         * Returns the velocity of this stuntdouble in specified snapshot 
         * @return the velocity of this stuntdouble
         * @param snapshotNo
         */    
         Vector3d getVel(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage_->velocity[localIndex_];
        }

       /**
         * Sets  the previous velocity of this stuntdouble
         * @param vel  new velocity 
         * @see #getVel
         */         
       void setPrevVel(const Vector3d& vel) {
            (snapshotMan_->getPrevSnapshot())->storage_->velocity[localIndex_] = vel;
       }
       
       /**
         * Sets  the current velocity of this stuntdouble
         * @param vel  new velocity 
         */         
        void setVel(const Vector3d& vel) {
            (snapshotMan_->getCurrentSnapshot())->storage_->velocity[localIndex_] = vel;
        }

       /**
         * Sets  the velocity of this stuntdouble in specified snapshot
         * @param vel velocity to be set 
         * @param snapshotNo 
         * @see #getVel
         */         
        void setVel(const Vector3d& vel, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->velocity[localIndex_] = vel;
        }

       /**
         * Returns the previous rotation matrix of this stuntdouble
         * @return the rotation matrix of this stuntdouble
         */    
        RotMat3x3d getPrevA() {
            return (snapshotMan_->getPrevSnapshot())->storage_->aMat[localIndex_];
        }
       
        /**
         * Returns the current rotation matrix of this stuntdouble
         * @return the rotation matrix of this stuntdouble
         */    
        RotMat3x3d getA() {
            return (snapshotMan_->getCurrentSnapshot())->storage_->aMat[localIndex_];
        }

       /**
         * Returns the rotation matrix of this stuntdouble in specified snapshot 
         *
         * @return the rotation matrix of this stuntdouble
         * @param snapshotNo
         */    
         RotMat3x3d getA(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage_->aMat[localIndex_];
        }

       /**
         * Sets  the previous rotation matrix of this stuntdouble
         * @param a  new rotation matrix 
         * @see #getA
         */         
       virtual void setPrevA(const RotMat3x3d& a) {
            (snapshotMan_->getPrevSnapshot())->storage_->aMat[localIndex_] = a;
       }
       
       /**
         * Sets  the current rotation matrix of this stuntdouble
         * @param a  new rotation matrix 
         */         
        virtual void setA(const RotMat3x3d& a) {
            (snapshotMan_->getCurrentSnapshot())->storage_->aMat[localIndex_] = a;
        }

       /**
         * Sets  the rotation matrix of this stuntdouble in specified snapshot
         * @param a rotation matrix to be set 
         * @param snapshotNo 
         * @see #getA
         */         
        virtual void setA(const RotMat3x3d& a, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->aMat[localIndex_] = a;
        }       

       /**
         * Returns the previous angular momentum of this stuntdouble (body-fixed).
         * @return the angular momentum of this stuntdouble
         */    
        Vector3d getPrevJ() {
            return (snapshotMan_->getPrevSnapshot())->storage_->angularMomentum[localIndex_];
        }
       
        /**
         * Returns the current angular momentum of this stuntdouble (body -fixed).
         * @return the angular momentum of this stuntdouble
         */    
        Vector3d getJ() {
            return (snapshotMan_->getCurrentSnapshot())->storage_->angularMomentum[localIndex_];
        }

       /**
         * Returns the angular momentum of this stuntdouble in specified snapshot (body-fixed).
         * @return the angular momentum of this stuntdouble
         * @param snapshotNo
         */    
         Vector3d getJ(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage_->angularMomentum[localIndex_];
        }

       /**
         * Sets  the previous angular momentum of this stuntdouble (body-fixed).
         * @param angMom  new angular momentum 
         * @see #getJ
         */         
       void setPrevJ(const Vector3d& angMom) {
            (snapshotMan_->getPrevSnapshot())->storage_->angularMomentum[localIndex_] = angMom;
       }
       
       /**
         * Sets  the current angular momentum of this stuntdouble (body-fixed).
         * @param angMom  new angular momentum 
         */         
        void setJ(const Vector3d& angMom) {
            (snapshotMan_->getCurrentSnapshot())->storage_->angularMomentum[localIndex_] = angMom;
        }

       /**
         * Sets the angular momentum of this stuntdouble in specified snapshot(body-fixed).
         * @param angMom angular momentum to be set 
         * @param snapshotNo 
         * @see #getJ
         */         
        void setJ(const Vector3d& angMom, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->angularMomentum[localIndex_] = angMom;
        }
         
       /**
         * Returns the previous quaternion of this stuntdouble
         * @return the quaternion of this stuntdouble
         */    
        Quat4d getPrevQ() {
            return (snapshotMan_->getPrevSnapshot())->storage_->aMat[localIndex_].toQuaternion();
        }
       
        /**
         * Returns the current quaternion of this stuntdouble
         * @return the quaternion of this stuntdouble
         */    
        Quat4d getQ() {
            return (snapshotMan_->getCurrentSnapshot())->storage_->aMat[localIndex_].toQuaternion();
        }

       /**
         * Returns the quaternion of this stuntdouble in specified snapshot 
         * @return the quaternion of this stuntdouble
         * @param snapshotNo
         */    
         Quat4d getQ(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage_->aMat[localIndex_].toQuaternion();
        }

       /**
         * Sets  the previous quaternion of this stuntdouble
         * @param q  new quaternion 
         * @note actual storage data is rotation matrix
         */         
       void setPrevQ(const Quat4d& q) {
            (snapshotMan_->getPrevSnapshot())->storage_->aMat[localIndex_] = q;
       }
       
       /**
         * Sets  the current quaternion of this stuntdouble
         * @param q  new quaternion 
         * @note actual storage data is rotation matrix
         */         
        void setQ(const Quat4d& q) {
            (snapshotMan_->getCurrentSnapshot())->storage_->aMat[localIndex_] = q;
        }

       /**
         * Sets  the quaternion of this stuntdouble in specified snapshot
         *
         * @param q quaternion to be set 
         * @param snapshotNo 
         * @note actual storage data is rotation matrix
         */         
        void setQ(const Quat4d& q, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->aMat[localIndex_] = q;
        }

       /**
         * Returns the previous euler angles of this stuntdouble
         * @return the euler angles of this stuntdouble
         */    
        Vector3d getPrevEuler() {
            return (snapshotMan_->getPrevSnapshot())->storage_->aMat[localIndex_].toEulerAngles();
        }
       
        /**
         * Returns the current euler angles of this stuntdouble
         * @return the euler angles of this stuntdouble
         */    
        Vector3d getEuler() {
            return (snapshotMan_->getCurrentSnapshot())->storage_->aMat[localIndex_].toEulerAngles();
        }

       /**
         * Returns the euler angles of this stuntdouble in specified snapshot.
         * @return the euler angles of this stuntdouble
         * @param snapshotNo
         */    
         Vector3d getEuler(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage_->aMat[localIndex_].toEulerAngles();
        }

       /**
         * Sets  the previous euler angles of this stuntdouble.
         * @param euler  new euler angles 
         * @see #getEuler
         * @note actual storage data is rotation matrix         
         */         
       void setPrevEuler(const Vector3d& euler) {
            (snapshotMan_->getPrevSnapshot())->storage_->aMat[localIndex_] = euler;
       }
       
       /**
         * Sets  the current euler angles of this stuntdouble
         * @param euler  new euler angles 
         */         
        void setEuler(const Vector3d& euler) {
            (snapshotMan_->getCurrentSnapshot())->storage_->aMat[localIndex_] = euler;
        }

       /**
         * Sets  the euler angles  of this stuntdouble in specified snapshot
         *
         * @param euler euler angles to be set 
         * @param snapshotNo 
         * @note actual storage data is rotation matrix                  
         */         
        void setEuler(const Vector3d& euler, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->aMat[localIndex_] = euler;
        }
       
       /**
         * Returns the previous unit vectors of this stuntdouble
         * @return the unit vectors of this stuntdouble
         */    
        Vector3d getPrevUnitVector() {
            return (snapshotMan_->getPrevSnapshot())->storage_->unitVector[localIndex_];
        }
       
        /**
         * Returns the current unit vectors of this stuntdouble
         * @return the unit vectors of this stuntdouble
         */    
        Vector3d getUnitVector() {
            return (snapshotMan_->getCurrentSnapshot())->storage_->unitVector[localIndex_];
        }

       /**
         * Returns the unit vectors of this stuntdouble in specified snapshot 
         *
         * @return the unit vectors of this stuntdouble
         * @param snapshotNo
         */    
         Vector3d getUnitVector(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage_->unitVector[localIndex_];
        }

       /**
         * Returns the previous force of this stuntdouble
         * @return the force of this stuntdouble
         */    
        Vector3d getPrevFrc() {
            return (snapshotMan_->getPrevSnapshot())->storage_->force[localIndex_];
        }
       
        /**
         * Returns the current force of this stuntdouble
         * @return the force of this stuntdouble
         */    
        Vector3d getFrc() {
            return (snapshotMan_->getCurrentSnapshot())->storage_->force[localIndex_];
        }

       /**
         * Returns the force of this stuntdouble in specified snapshot 
         *
         * @return the force of this stuntdouble
         * @param snapshotNo
         */    
         Vector3d getFrc(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage_->force[localIndex_];
        }

       /**
         * Sets  the previous force of this stuntdouble
         *
         * @param frc  new force 
         * @see #getFrc
         */         
       void setPrevFrc(const Vector3d& frc) {
            (snapshotMan_->getPrevSnapshot())->storage_->force[localIndex_] = frc;
       }
       
       /**
         * Sets  the current force of this stuntdouble
         * @param frc  new force 
         */         
        void setFrc(const Vector3d& frc) {
            (snapshotMan_->getCurrentSnapshot())->storage_->force[localIndex_] = frc;
        }

       /**
         * Sets  the force of this stuntdouble in specified snapshot
         *
         * @param frc force to be set 
         * @param snapshotNo 
         * @see #getFrc
         */         
        void setFrc(const Vector3d& frc, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->force[localIndex_] = frc;
        }

       /**
         * Adds force into the previous force of this stuntdouble
         *
         * @param frc  new force 
         * @see #getFrc
         */         
       void setPrevFrc(const Vector3d& frc) {
            (snapshotMan_->getPrevSnapshot())->storage_->force[localIndex_] += frc;
       }
       
       /**
         * Adds force into the current force of this stuntdouble
         * @param frc  new force 
         */         
        void setFrc(const Vector3d& frc) {
            (snapshotMan_->getCurrentSnapshot())->storage_->force[localIndex_] += frc;
        }

       /**
         * Adds force into the force of this stuntdouble in specified snapshot
         *
         * @param frc force to be set 
         * @param snapshotNo 
         * @see #getFrc
         */         
        void setFrc(const Vector3d& frc, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->force[localIndex_] += frc;
        }

       /**
         * Returns the previous torque of this stuntdouble
         * @return the torque of this stuntdouble
         */    
        Vector3d getPrevTrq() {
            return (snapshotMan_->getPrevSnapshot())->storage_->torque[localIndex_];
        }
       
        /**
         * Returns the current torque of this stuntdouble
         * @return the torque of this stuntdouble
         */    
        Vector3d getTrq() {
            return (snapshotMan_->getCurrentSnapshot())->storage_->torque[localIndex_];
        }

       /**
         * Returns the torque of this stuntdouble in specified snapshot 
         *
         * @return the torque of this stuntdouble
         * @param snapshotNo
         */    
         Vector3d getTrq(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage_->torque[localIndex_];
        }

       /**
         * Sets  the previous torque of this stuntdouble
         *
         * @param trq  new torque 
         * @see #getTrq
         */         
       void setPrevTrq(const Vector3d& trq) {
            (snapshotMan_->getPrevSnapshot())->storage_->torque[localIndex_] = trq;
       }
       
       /**
         * Sets  the current torque of this stuntdouble
         * @param trq  new torque 
         */         
        void setTrq(const Vector3d& trq) {
            (snapshotMan_->getCurrentSnapshot())->storage_->torque[localIndex_] = trq;
        }

       /**
         * Sets  the torque of this stuntdouble in specified snapshot
         *
         * @param trq torque to be set 
         * @param snapshotNo 
         * @see #getTrq
         */         
        void setTrq(const Vector3d& trq, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->torque[localIndex_] = trq;
        }

       /**
         * Adds torque into the previous torque of this stuntdouble
         *
         * @param trq  new torque 
         * @see #getTrq
         */         
       void addPrevTrq(const Vector3d& trq) {
            (snapshotMan_->getPrevSnapshot())->storage_->torque[localIndex_] += trq;
       }
       
       /**
         * Adds torque into the current torque of this stuntdouble
         * @param trq  new torque 
         */         
        void addTrq(const Vector3d& trq) {
            (snapshotMan_->getCurrentSnapshot())->storage_->torque[localIndex_] += trq;
        }

       /**
         * Adds torque into the torque of this stuntdouble in specified snapshot
         *
         * @param trq torque to be add 
         * @param snapshotNo 
         * @see #getTrq
         */         
        void addTrq(const Vector3d& trq, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->torque[localIndex_] += trq;
        }       


       /**
         * Returns the previous z-angle of this stuntdouble
         * @return the z-angle of this stuntdouble
         */    
        Vector3d getPrevZangle() {
            return (snapshotMan_->getPrevSnapshot())->storage_->zAngle[localIndex_];
        }
       
        /**
         * Returns the current z-angle of this stuntdouble
         * @return the z-angle of this stuntdouble
         */    
        Vector3d getZangle() {
            return (snapshotMan_->getCurrentSnapshot())->storage_->zAngle[localIndex_];
        }

       /**
         * Returns the z-angle of this stuntdouble in specified snapshot
         * @return the z-angle of this stuntdouble
         * @param snapshotNo
         */    
         Vector3d getZangle(int snapshotNo) {
            return (snapshotMan_->getSnapshot(snapShotNo))->storage_->zAngle[localIndex_];
        }

       /**
         * Sets  the previous z-angle of this stuntdouble
         * @param angle  new z-angle 
         * @see #getZangle
         */         
       void setPrevZangle(const Vector3d& angle) {
            (snapshotMan_->getPrevSnapshot())->storage_->zAngle[localIndex_] = angle;
       }
       
       /**
         * Sets  the current z-angle of this stuntdouble
         * @param angle  new z-angle 
         */         
        void setZangle(const Vector3d& angle) {
            (snapshotMan_->getCurrentSnapshot())->storage_->zAngle[localIndex_] = angle;
        }

       /**
         * Sets  the z-angle of this stuntdouble in specified snapshot
         * @param angle z-angle to be set 
         * @param snapshotNo 
         * @see #getZangle
         */         
        void setZangle(const Vector3d& angle, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->zAngle[localIndex_] = angle;
        }

       /**
         * Adds z-angle into the previous z-angle of this stuntdouble
         * @param angle  new z-angle 
         * @see #getZangle
         */         
       void addPrevZangle(const Vector3d& angle) {
            (snapshotMan_->getPrevSnapshot())->storage_->zAngle[localIndex_] += angle;
       }
       
       /**
         * Adds z-angle into the current z-angle of this stuntdouble
         * @param angle  new z-angle 
         */         
        void addZangle(const Vector3d& angle) {
            (snapshotMan_->getCurrentSnapshot())->storage_->zAngle[localIndex_] += angle;
        }

       /**
         * Adds z-angle into the z-angle of this stuntdouble in specified snapshot
         * @param angle z-angle to be add 
         * @param snapshotNo 
         * @see #getZangle
         */         
        void addZangle(const Vector3d& angle, int snapshotNo) {
            (snapshotMan_->getSnapshot(snapshotNo))->storage_->zAngle[localIndex_] += angle;
        }       

       /**
         * Returns the inertia tensor of this stuntdouble
         * @return the inertia tensor of this stuntdouble
         */ 
        virtual Mat3x3d getI() = 0;

       /**
         * Returns the gradient of this stuntdouble
         * @return the gradient of this stuntdouble
         */ 
        virtual std::vector<double> getGrad() = 0;

       /**
         * Tests the  if this stuntdouble is a  linear rigidbody
         *
         * @return true if this stuntdouble is a  linear rigidbody, otherwise return false
         * @note atom and directional atom will always return false
         * 
         * @see #linearAxis
         */         
        bool isLinear() {
            return linear_;
        }

       /**
         * Returns the linear axis of the rigidbody, atom and directional atom will always return -1
         *
         * @return the linear axis of the rigidbody
         * 
         * @see #isLinear
         */ 
        int linearAxis() {
            return linearAxis_;
        }

        /** Returns the mass of this stuntdouble */
        double getMass() {
            return mass_;
        }

        /**
         * Sets the mass of this stuntdoulbe
         * @param mass the mass to be set
         */         
        void setMass(double mass) {
            mass_ = mass;
        }

        /** Returns the name of this stuntdouble */
        std::string getType();
        
        /** Sets the name of this stuntdouble*/
        void setType(const std::string& name);

        /**
         * Converts a lab fixed vector to a body fixed vector.
         * @return body fixed vector
         * @param v lab fixed vector
         */
        Vector3d lab2Body(const Vector3d& v);

        /**
         * Converts a body fixed vector to a lab fixed vector.
         * @return corresponding lab fixed vector
         * @param v body fixed vector
         */
        Vector3d body2Lab(const Vector3d& v);
        /**
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

        //below functions are just forward functions
        /**
         * Adds property into property map
         * @param genData GenericData to be added into PropertyMap
         */
        void addProperty(GenericData* genData);

        /**
         * Removes property from PropertyMap by name
         * @param propName the name of property to be removed
         */
        void removeProperty(std::string& propName);

        /**
         * clear all of the properties
         */
        void clearProperties();

        /**
         * Returns all names of properties
         * @return all names of properties
         */
        std::vector<std::string> getPropertyNames();

        /**
         * Returns all of the properties in PropertyMap
         * @return all of the properties in PropertyMap
         */      
        std::vector<GenericData*> getProperties();

        /**
         * Returns property 
         * @param propName name of property
         * @return a pointer point to property with propName. If no property named propName
         * exists, return NULL
         */      
        GenericData* getPropertyByName(std:string& propName);

    protected:
        
        StuntDouble();
        
        StuntDouble(const StuntDouble& sd);
        StuntDouble& operator=(const StuntDouble& sd);

        ObjectType objType_;

        bool linear_;
        int linearAxis_;        

        DataStoragePointer storage_;
        SnapshotManager* snapshotMan_;
        
    private:
        
        int globalIndex_;
        int localIndex_;

        std::string name_;

        double mass_;
        
        PropertyMap properties_;
    };

}//end namespace oopse
#endif //ifndef _STUNTDOUBLE_HPP_
