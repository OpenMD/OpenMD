/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
#include "brains/Stats.hpp"
#include "UseTheForce/DarkSide/simulation_interface.h"

namespace oopse{

  /**
   * @class Snapshot Snapshot.hpp "brains/Snapshot.hpp"
   * @brief Snapshot class is a repository class for storing dynamic data during 
   *  Simulation
   * Every snapshot class will contain one DataStorage  for atoms and one DataStorage
   *  for rigid bodies.
   */
  class Snapshot {
  public:
            
    Snapshot(int nAtoms, int nRigidbodies) : atomData(nAtoms), rigidbodyData(nRigidbodies),
					     currentTime_(0), orthoRhombic_(0), chi_(0.0), integralOfChiDt_(0.0), eta_(0.0), id_(-1) {

    }

    Snapshot(int nAtoms, int nRigidbodies, int storageLayout) 
      : atomData(nAtoms, storageLayout), rigidbodyData(nRigidbodies, storageLayout),
	currentTime_(0), orthoRhombic_(0), chi_(0.0), integralOfChiDt_(0.0), eta_(0.0), id_(-1) {

      }
            
    /** Returns the id of this Snapshot */
    int getID() {
      return id_;
    }

    /** Sets the id of this Snapshot */
    void setID(int id) {
      id_ = id;
    }

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
    void setHmat(const Mat3x3d& m);
            
    RealType getVolume() {
      return hmat_.determinant();
    }

    /** Returns the inverse H-Matrix */
    Mat3x3d getInvHmat() {
      return invHmat_;
    }

    /** Wrapping the vector according to periodic boundary condition*/
    void wrapVector(Vector3d& v);

            
    RealType getTime() {
      return currentTime_;
    }

    void increaseTime(RealType dt) {
      setTime(getTime() + dt);
    }

    void setTime(RealType time) {
      currentTime_ =time;
      //time at statData is redundant
      statData[Stats::TIME] = currentTime_;
    }

    RealType getChi() {
      return chi_;
    }

    void setChi(RealType chi) {
      chi_ = chi;
    }

    RealType getIntegralOfChiDt() {
      return integralOfChiDt_;
    }

    void setIntegralOfChiDt(RealType integralOfChiDt) {
      integralOfChiDt_ = integralOfChiDt;
    }
            
    Mat3x3d getEta() {
      return eta_;
    }

    void setEta(const Mat3x3d& eta) {
      eta_ = eta;
    }
            
    DataStorage atomData;
    DataStorage rigidbodyData;
    Stats statData;
            
  private:
    RealType currentTime_;

    Mat3x3d hmat_;
    Mat3x3d invHmat_;
    int orthoRhombic_;

    RealType chi_;
    RealType integralOfChiDt_;
    Mat3x3d eta_;
            
    int id_; /**< identification number of the snapshot */
  };

  typedef DataStorage (Snapshot::*DataStoragePointer); 
}
#endif //BRAINS_SNAPSHOT_HPP
