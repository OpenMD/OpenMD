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

#ifndef CONSTRAINTS_CONTRAINTELEM_HPP
#define CONSTRAINTS_CONTRAINTELEM_HPP

#include "primitives/StuntDouble.hpp"
#include "utils/GenericData.hpp"
#include "utils/simError.h"
namespace OpenMD {
  
  /**
   * @class ConstraintElem ConstraintElem.hpp "constraints/ConstraintElem.hpp"
   * An adapter class of StuntDouble which is used at constraint algorithm
   */
  
  class ConstraintElem{
  public:
    ConstraintElem(StuntDouble* sd) : sd_(sd) {
      GenericData* movedData = sd_->getPropertyByName("Moved");
      if (movedData !=NULL) { //if generic data with keyword "moved" exists, assign it to moved_
	moved_ = dynamic_cast<BoolGenericData*>(movedData);
	if (moved_ == NULL) {
	  sprintf(painCave.errMsg,
		  "Generic Data with keyword Moved exists, however, it can not be casted to a BoolGenericData.\n");
	  painCave.isFatal = 1;;
	  simError();
	}
      }else { //if generic data with keyword "moved" does not exist, create it
	moved_ = new BoolGenericData("Moved");
	sd_->addProperty(moved_);
      }
      
      GenericData* movingData = sd_->getPropertyByName("Moving");
      if (movingData !=NULL) {
	moving_ = dynamic_cast<BoolGenericData*>(movingData);
	if (moving_ == NULL) {
	  sprintf(painCave.errMsg,
		  "Generic Data with keyword Moving exists, however, it can not be casted to a BoolGenericData.\n");
	  painCave.isFatal = 1;;
	  simError();
	}
      }else {
	moving_ = new BoolGenericData("Moving");
	sd_->addProperty(moving_);                
      }
      
    }
    
    bool getMoved() { return moved_->getData(); }
    void setMoved(bool moved) { moved_->setData(moved);}
    
    bool getMoving() { return moving_->getData(); }
    void setMoving(bool moving) { moving_->setData(moving); }
    
    StuntDouble* getStuntDouble() { return sd_; }
    
    /**
     * Returns the global index of this stuntRealType.
     * @return  the global index of this stuntdouble 
     */
    int getGlobalIndex() {
      return sd_->getGlobalIndex();
    }
    
    /**
     * Sets the global index of this stuntRealType.
     * @param new global index to be set
     */
    void setGlobalIndex(int index) {
      sd_->setGlobalIndex(index);
    }
    
    /** 
     * Returns the local index of this stuntdouble 
     * @return the local index of this stuntdouble
     */
    int getLocalIndex() {
      return sd_->getLocalIndex();
    }
    
    /**
     * Sets the local index of this stuntdouble
     * @param index new index to be set
     */        
    void setLocalIndex(int index) {
      sd_->setLocalIndex(index);
    }
    
    /**
     * Sets the Snapshot Manager of this stuntdouble
     */
    void setSnapshotManager(SnapshotManager* sman) {
      sd_->setSnapshotManager(sman);
    }
    
    /**
     * Tests if this stuntdouble is an atom 
     * @return true is this stuntdouble is an atom(or a directional atom), return false otherwise
     */
    bool isAtom(){
      return sd_->isAtom();
    }
    
    /** 
     * Tests if this stuntdouble is an directional atom 
     * @return true if this stuntdouble is an directional atom, return false otherwise
     */
    bool isDirectionalAtom(){
      return sd_->isDirectional();
    }
    
    /**
     * Tests if this stuntdouble is a rigid body. 
     * @return true if this stuntdouble is a rigid body, otherwise return false
     */
    bool isRigidBody(){
      return sd_->isRigidBody();
    }
    
    /**
     * Tests if this stuntdouble is a directional one. 
     * @return true is this stuntdouble is a directional atom or a rigid body, return false otherwise
     */
    bool isDirectional(){
      return sd_->isDirectional();
    }
    
    /**
     * Returns the previous position of this stuntdouble
     * @return the position of this stuntdouble
     */    
    Vector3d getPrevPos() {
      return sd_->getPrevPos();
    }
    
    /**
     * Returns the current position of this stuntdouble
     * @return the position of this stuntdouble
     */    
    Vector3d getPos() {
      return sd_->getPos();
    }
    
    /**
     * Returns the position of this stuntdouble in specified snapshot 
     * @return the position of this stuntdouble
     * @param snapshotNo
     */    
    Vector3d getPos(int snapshotNo) {
      return sd_->getPos(snapshotNo);
    }
    
    /**
     * Sets  the previous position of this stuntdouble
     * @param pos  new position 
     * @see #getPos
     */         
    void setPrevPos(const Vector3d& pos) {
      sd_->setPrevPos(pos);
    }
    
    /**
     * Sets  the current position of this stuntdouble
     * @param pos  new position 
     */         
    void setPos(const Vector3d& pos) {
      sd_->setPos(pos);
    }
    
    /**
     * Sets  the position of this stuntdouble in specified snapshot
     * @param pos position to be set 
     * @param snapshotNo 
     * @see #getPos
     */         
    void setPos(const Vector3d& pos, int snapshotNo) {
      sd_->setPos(pos, snapshotNo);
    }
    
    /**
     * Returns the previous velocity of this stuntdouble
     * @return the velocity of this stuntdouble
     */    
    Vector3d getPrevVel() {
      return sd_->getPrevVel();
    }
    
    /**
     * Returns the current velocity of this stuntdouble
     * @return the velocity of this stuntdouble
     */    
    Vector3d getVel() {
      return sd_->getVel();
    }
    
    /**
     * Returns the velocity of this stuntdouble in specified snapshot 
     * @return the velocity of this stuntdouble
     * @param snapshotNo
     */    
    
    Vector3d getVel(int snapshotNo) {
      return sd_->getVel(snapshotNo);
    }
    
    /**
     * Sets  the previous velocity of this stuntdouble
     * @param vel  new velocity 
     * @see #getVel
     */         
    void setPrevVel(const Vector3d& vel) {
      sd_->setPrevVel(vel);
    }
    
    /**
     * Sets  the current velocity of this stuntdouble
     * @param vel  new velocity 
     */         
    void setVel(const Vector3d& vel) {
      sd_->setVel(vel);
    }
    
    /**
     * Sets  the velocity of this stuntdouble in specified snapshot
     * @param vel velocity to be set 
     * @param snapshotNo 
     * @see #getVel
     */         
    void setVel(const Vector3d& vel, int snapshotNo) {
      sd_->setVel(vel, snapshotNo);
    }
    
    /**
     * Returns the previous rotation matrix of this stuntdouble
     * @return the rotation matrix of this stuntdouble
     */    
    RotMat3x3d getPrevA() {
      return sd_->getPrevA();
    }
    
    /**
     * Returns the current rotation matrix of this stuntdouble
     * @return the rotation matrix of this stuntdouble
     */    
    RotMat3x3d getA() {
      return sd_->getA();
    }
    
    /**
     * Returns the rotation matrix of this stuntdouble in specified snapshot 
     *
     * @return the rotation matrix of this stuntdouble
     * @param snapshotNo
     */    
    RotMat3x3d getA(int snapshotNo) {
      return sd_->getA(snapshotNo);
    }
    
    /**
     * Sets  the previous rotation matrix of this stuntdouble
     * @param a  new rotation matrix 
     * @see #getA
     */         
    void setPrevA(const RotMat3x3d& a) {
      sd_->setPrevA(a);
    }
    
    /**
     * Sets  the current rotation matrix of this stuntdouble
     * @param a  new rotation matrix 
     */         
    void setA(const RotMat3x3d& a) {
      sd_->setA(a);
    }
    
    /**
     * Sets  the rotation matrix of this stuntdouble in specified snapshot
     * @param a rotation matrix to be set 
     * @param snapshotNo 
     * @see #getA
     */         
    void setA(const RotMat3x3d& a, int snapshotNo) {
      sd_->setA(a, snapshotNo);
    }       
    
    /**
     * Returns the previous angular momentum of this stuntdouble (body-fixed).
     * @return the angular momentum of this stuntdouble
     */    
    Vector3d getPrevJ() {
      return sd_->getPrevJ();
    }
    
    /**
     * Returns the current angular momentum of this stuntdouble (body -fixed).
     * @return the angular momentum of this stuntdouble
     */    
    Vector3d getJ() {
      return sd_->getJ();
    }
    
    /**
     * Returns the angular momentum of this stuntdouble in specified snapshot (body-fixed).
     * @return the angular momentum of this stuntdouble
     * @param snapshotNo
     */    
    Vector3d getJ(int snapshotNo) {
      return sd_->getJ(snapshotNo);
    }
    
    /**
     * Sets  the previous angular momentum of this stuntdouble (body-fixed).
     * @param angMom  new angular momentum 
     * @see #getJ
     */         
    void setPrevJ(const Vector3d& angMom) {
      sd_->setPrevJ(angMom);
    }
    
    /**
     * Sets  the current angular momentum of this stuntdouble (body-fixed).
     * @param angMom  new angular momentum 
     */         
    void setJ(const Vector3d& angMom) {
      sd_->setJ(angMom);
    }
    
    /**
     * Sets the angular momentum of this stuntdouble in specified snapshot(body-fixed).
     * @param angMom angular momentum to be set 
     * @param snapshotNo 
     * @see #getJ
     */         
    void setJ(const Vector3d& angMom, int snapshotNo) {
      sd_->setJ(angMom, snapshotNo);
    }
    
    /**
     * Returns the previous quaternion of this stuntdouble
     * @return the quaternion of this stuntdouble
     */    
    Quat4d getPrevQ() {
      return sd_->getPrevQ();
    }
    
    /**
     * Returns the current quaternion of this stuntdouble
     * @return the quaternion of this stuntdouble
     */    
    Quat4d getQ() {
      return sd_->getQ();
    }
    
    /**
     * Returns the quaternion of this stuntdouble in specified snapshot 
     * @return the quaternion of this stuntdouble
     * @param snapshotNo
     */    
    Quat4d getQ(int snapshotNo) {
      return sd_->getQ(snapshotNo);
    }
    
    /**
     * Sets  the previous quaternion of this stuntdouble
     * @param q  new quaternion 
     * @note actual storage data is rotation matrix
     */         
    void setPrevQ(const Quat4d& q) {
      sd_->setPrevQ(q);
    }
    
    /**
     * Sets  the current quaternion of this stuntdouble
     * @param q  new quaternion 
     * @note actual storage data is rotation matrix
     */         
    void setQ(const Quat4d& q) {
      sd_->setQ(q);
    }
    
    /**
     * Sets  the quaternion of this stuntdouble in specified snapshot
     *
     * @param q quaternion to be set 
     * @param snapshotNo 
     * @note actual storage data is rotation matrix
     */         
    void setQ(const Quat4d& q, int snapshotNo) {
      sd_->setQ(q, snapshotNo);
    }
    
    /**
     * Returns the previous euler angles of this stuntdouble
     * @return the euler angles of this stuntdouble
     */    
    Vector3d getPrevEuler() {
      return sd_->getPrevEuler();
    }
    
    /**
     * Returns the current euler angles of this stuntdouble
     * @return the euler angles of this stuntdouble
     */    
    Vector3d getEuler() {
      return sd_->getEuler();
    }
    
    /**
     * Returns the euler angles of this stuntdouble in specified snapshot.
     * @return the euler angles of this stuntdouble
     * @param snapshotNo
     */    
    Vector3d getEuler(int snapshotNo) {
      return sd_->getEuler(snapshotNo);
    }
    
    /**
     * Sets  the previous euler angles of this stuntRealType.
     * @param euler  new euler angles 
     * @see #getEuler
     * @note actual storage data is rotation matrix         
     */         
    void setPrevEuler(const Vector3d& euler) {
      sd_->setPrevEuler(euler);
    }
    
    /**
     * Sets  the current euler angles of this stuntdouble
     * @param euler  new euler angles 
     */         
    void setEuler(const Vector3d& euler) {
      sd_->setEuler(euler);
    }
    
    /**
     * Sets  the euler angles  of this stuntdouble in specified snapshot
     *
     * @param euler euler angles to be set 
     * @param snapshotNo 
     * @note actual storage data is rotation matrix                  
     */         
    void setEuler(const Vector3d& euler, int snapshotNo) {
      sd_->setEuler(euler, snapshotNo);
    }
    
    /**
     * Returns the previous unit vectors of this stuntdouble
     * @return the unit vectors of this stuntdouble
     */    
    RotMat3x3d getPrevElectroFrame() {
      return sd_->getPrevElectroFrame();
    }
    
    /**
     * Returns the current unit vectors of this stuntdouble
     * @return the unit vectors of this stuntdouble
     */    
    RotMat3x3d getElectroFrame() {
      return sd_->getElectroFrame();
    }
    
    /**
     * Returns the unit vectors of this stuntdouble in specified snapshot 
     *
     * @return the unit vectors of this stuntdouble
     * @param snapshotNo
     */    
    RotMat3x3d getElectroFrame(int snapshotNo) {
      return sd_->getElectroFrame(snapshotNo);
    }
    
    /**
     * Returns the previous force of this stuntdouble
     * @return the force of this stuntdouble
     */    
    Vector3d getPrevFrc() {
      return sd_->getPrevFrc();
    }
    
    /**
     * Returns the current force of this stuntdouble
     * @return the force of this stuntdouble
     */    
    Vector3d getFrc() {
      return sd_->getFrc();
    }
    
    /**
     * Returns the force of this stuntdouble in specified snapshot 
     *
     * @return the force of this stuntdouble
     * @param snapshotNo
     */    
    Vector3d getFrc(int snapshotNo) {
      return sd_->getFrc(snapshotNo);
    }
    
    /**
     * Sets  the previous force of this stuntdouble
     *
     * @param frc  new force 
     * @see #getFrc
     */         
    void setPrevFrc(const Vector3d& frc) {
      sd_->setPrevFrc(frc);
    }
    
    /**
     * Sets  the current force of this stuntdouble
     * @param frc  new force 
     */         
    void setFrc(const Vector3d& frc) {
      sd_->setFrc(frc);
    }
    
    /**
     * Sets  the force of this stuntdouble in specified snapshot
     *
     * @param frc force to be set 
     * @param snapshotNo 
     * @see #getFrc
     */         
    void setFrc(const Vector3d& frc, int snapshotNo) {
      sd_->setFrc(frc, snapshotNo);
    }
    
    /**
     * Adds force into the previous force of this stuntdouble
     *
     * @param frc  new force 
     * @see #getFrc
     */         
    void addPrevFrc(const Vector3d& frc) {
      sd_->addPrevFrc(frc);
    }
    
    /**
     * Adds force into the current force of this stuntdouble
     * @param frc  new force 
     */         
    void addFrc(const Vector3d& frc) {
      sd_->addFrc(frc);
    }
    
    /**
     * Adds force into the force of this stuntdouble in specified snapshot
     *
     * @param frc force to be set 
     * @param snapshotNo 
     * @see #getFrc
     */         
    void addFrc(const Vector3d& frc, int snapshotNo) {
      sd_->addFrc(frc, snapshotNo);
    }
    
    /**
     * Returns the previous torque of this stuntdouble
     * @return the torque of this stuntdouble
     */    
    Vector3d getPrevTrq() {
      return sd_->getPrevTrq();
    }
    
    /**
     * Returns the current torque of this stuntdouble
     * @return the torque of this stuntdouble
     */    
    Vector3d getTrq() {
      return sd_->getTrq();
    }
    
    /**
     * Returns the torque of this stuntdouble in specified snapshot 
     *
     * @return the torque of this stuntdouble
     * @param snapshotNo
     */    
    Vector3d getTrq(int snapshotNo) {
      return sd_->getTrq(snapshotNo);
    }
    
    /**
     * Sets  the previous torque of this stuntdouble
     *
     * @param trq  new torque 
     * @see #getTrq
     */         
    void setPrevTrq(const Vector3d& trq) {
      sd_->setPrevTrq(trq);
    }
    
    /**
     * Sets  the current torque of this stuntdouble
     * @param trq  new torque 
     */         
    void setTrq(const Vector3d& trq) {
      sd_->setTrq(trq);
    }
    
    /**
     * Sets  the torque of this stuntdouble in specified snapshot
     *
     * @param trq torque to be set 
     * @param snapshotNo 
     * @see #getTrq
     */         
    void setTrq(const Vector3d& trq, int snapshotNo) {
      sd_->setTrq(trq, snapshotNo);
    }
    
    /**
     * Adds torque into the previous torque of this stuntdouble
     *
     * @param trq  new torque 
     * @see #getTrq
     */         
    void addPrevTrq(const Vector3d& trq) {
      sd_->addPrevTrq(trq);
    }
    
    /**
     * Adds torque into the current torque of this stuntdouble
     * @param trq  new torque 
     */         
    void addTrq(const Vector3d& trq) {
      sd_->addTrq(trq);
    }
    
    /**
     * Adds torque into the torque of this stuntdouble in specified snapshot
     *
     * @param trq torque to be add 
     * @param snapshotNo 
     * @see #getTrq
     */         
    void addTrq(const Vector3d& trq, int snapshotNo) {
      sd_->addTrq(trq, snapshotNo);
    }       
    
    /** Set the force of this stuntdouble to zero */
    void zeroForcesAndTorques() {
      sd_->zeroForcesAndTorques();
    }
    /**
     * Returns the inertia tensor of this stuntdouble
     * @return the inertia tensor of this stuntdouble
     */ 
    Mat3x3d getI() {
      return sd_->getI();
    }
    
    /**
     * Returns the gradient of this stuntdouble
     * @return the gradient of this stuntdouble
     */ 
    std::vector<RealType> getGrad() {
      return sd_->getGrad();
    }
    
    /**
     * Tests the  if this stuntdouble is a  linear rigidbody
     *
     * @return true if this stuntdouble is a  linear rigidbody, otherwise return false
     * @note atom and directional atom will always return false
     * 
     * @see #linearAxis
     */         
    bool isLinear() {
      return sd_->isLinear();
    }
    
    /**
     * Returns the linear axis of the rigidbody, atom and directional atom will always return -1
     *
     * @return the linear axis of the rigidbody
     * 
     * @see #isLinear
     */ 
    int linearAxis() {
      return sd_->linearAxis();
    }
    
    /** Returns the mass of this stuntdouble */
    RealType getMass() {
      return sd_->getMass();
    }
    
    /**
     * Sets the mass of this stuntdoulbe
     * @param mass the mass to be set
     */         
    void setMass(RealType mass) {
      sd_->setMass(mass);
    }
    
    /** Returns the name of this stuntdouble */
    std::string getType() {
      return sd_->getType();
    }
    
    /** Sets the name of this stuntRealType*/
    void setType(const std::string& name) {
      sd_->setType(name);
    }
    
    /**
     * Converts a lab fixed vector to a body fixed vector.
     * @return body fixed vector
     * @param v lab fixed vector
     */
    Vector3d lab2Body(const Vector3d& v) {
      return sd_->lab2Body(v);
    }
    
    /**
     * Converts a body fixed vector to a lab fixed vector.
     * @return corresponding lab fixed vector
     * @param v body fixed vector
     */
    Vector3d body2Lab(const Vector3d& v){
      return sd_->body2Lab(v);
    }
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
    void accept(BaseVisitor* v) {
      sd_->accept(v);
    }
    
    //below functions are just forward functions
    /**
     * Adds property into property map
     * @param genData GenericData to be added into PropertyMap
     */
    void addProperty(GenericData* genData){
      sd_->addProperty(genData);
    }
    
    /**
     * Removes property from PropertyMap by name
     * @param propName the name of property to be removed
     */
    void removeProperty(const std::string& propName) {
      sd_->removeProperty(propName);
    }
    
    /**
     * clear all of the properties
     */
    void clearProperties() {
      sd_->clearProperties();
    }
    
    /**
     * Returns all names of properties
     * @return all names of properties
     */
    std::vector<std::string> getPropertyNames() {
      return sd_->getPropertyNames();
    }
    
    /**
     * Returns all of the properties in PropertyMap
     * @return all of the properties in PropertyMap
     */      
    std::vector<GenericData*> getProperties() {
      return sd_->getProperties();
    }
    
    /**
     * Returns property 
     * @param propName name of property
     * @return a pointer point to property with propName. If no property named propName
     * exists, return NULL
     */      
    GenericData* getPropertyByName(const std::string& propName) {
      return sd_->getPropertyByName(propName);
    }
    
  private:
    StuntDouble* sd_;
    BoolGenericData* moved_;
    BoolGenericData* moving_;
  };
  
}

#endif
