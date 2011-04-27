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

namespace OpenMD{

  /**
   * @class Snapshot Snapshot.hpp "brains/Snapshot.hpp"
   * @brief Snapshot class is a repository class for storing dynamic data during 
   *  Simulation
   * Every snapshot class will contain one DataStorage for atoms and one DataStorage
   *  for rigid bodies.
   */
  class Snapshot {
  public:
            
    Snapshot(int nAtoms, int nRigidbodies, 
             int nCutoffGroups) : atomData(nAtoms), 
                                  rigidbodyData(nRigidbodies),
                                  cgData(nCutoffGroups, DataStorage::dslPosition),
                                  currentTime_(0), 
                                  orthoTolerance_(1e-6), 
                                  orthoRhombic_(0), 
                                  chi_(0.0), 
                                  integralOfChiDt_(0.0), 
                                  eta_(0.0), id_(-1), hasCOM_(false), 
                                  hasVolume_(false), volume_(0.0) {

    }

    Snapshot(int nAtoms, int nRigidbodies, int nCutoffGroups, 
             int storageLayout) : atomData(nAtoms, storageLayout), 
                                  rigidbodyData(nRigidbodies, storageLayout),
                                  cgData(nCutoffGroups, DataStorage::dslPosition),
                                  currentTime_(0), orthoTolerance_(1e-6), 
                                  orthoRhombic_(0), chi_(0.0), 
                                  integralOfChiDt_(0.0), eta_(0.0), id_(-1), 
                                  hasCOM_(false), hasVolume_(false),
                                  volume_(0.0)  {
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

    /** Returns the number of rigid bodies */
    int getNumberOfCutoffGroups() {
      return cgData.getSize();
    }

    /** Returns the H-Matrix */
    Mat3x3d getHmat() {
      return hmat_;
    }

    /** Sets the H-Matrix */
    void setHmat(const Mat3x3d& m);
            
    RealType getVolume() {
      if (hasVolume_){
        return volume_;
      }else{
        return hmat_.determinant();
      }
    }

    void setVolume(RealType volume){
      hasVolume_=true;
      volume_ = volume;
    }

    /** Returns the inverse H-Matrix */
    Mat3x3d getInvHmat() {
      return invHmat_;
    }

    /** Wrapping the vector according to periodic boundary condition*/
    void wrapVector(Vector3d& v);
    Vector3d getCOM();
    Vector3d getCOMvel();
    Vector3d getCOMw();
            
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
            

    void setOrthoTolerance(RealType orthoTolerance) {
      orthoTolerance_ = orthoTolerance;
    }

    Mat3x3d getEta() {
      return eta_;
    }

    void setEta(const Mat3x3d& eta) {
      eta_ = eta;
    }

    bool hasCOM() {
      return hasCOM_;
    }

    void setCOMprops(const Vector3d& COM, const Vector3d& COMvel, const Vector3d& COMw) {
      COM_ = COM;
      COMvel_ = COMvel;
      COMw_ = COMw;
      hasCOM_ = true;
    }

    Vector3d getAtomPosByIindex(int iIndex) {
#ifdef IS_MPI
      return atomIData.position[iIndex];
#else
      return atomData.position[iIndex];
#endif
    }
    Vector3d getAtomPosByJindex(int jIndex) {
#ifdef IS_MPI
      return atomJData.position[jIndex];
#else
      return atomData.position[jIndex];
#endif
    }

    Vector3d getCutoffGroupPosByIindex(int iIndex) {
#ifdef IS_MPI
      return cgIData.position[iIndex];
#else
      return cgData.position[iIndex];
#endif
    }
    Vector3d getCutoffGroupPosByJindex(int jIndex) {
#ifdef IS_MPI
      return cgJData.position[jIndex];
#else
      return cgData.position[jIndex];
#endif
    }

    DataStorage atomData;
    DataStorage rigidbodyData;
    DataStorage cgData;
    Stats statData;

#ifdef IS_MPI
    DataStorage atomIData;
    DataStorage atomJData;
    DataStorage cgIData;
    DataStorage cgJData;
#endif
   
            
  private:
    RealType currentTime_;

    Mat3x3d hmat_;
    Mat3x3d invHmat_;
    RealType orthoTolerance_;
    int orthoRhombic_;
    RealType volume_;

    RealType chi_;
    RealType integralOfChiDt_;
    Mat3x3d eta_;
    Vector3d COM_;
    Vector3d COMvel_;
    Vector3d COMw_;
    int id_; /**< identification number of the snapshot */
    bool hasCOM_;
    bool hasVolume_;
            
  };

  typedef DataStorage (Snapshot::*DataStoragePointer); 
}
#endif //BRAINS_SNAPSHOT_HPP
