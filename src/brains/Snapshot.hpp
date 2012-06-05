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

  struct FrameData {
    int id;                   /**< identification number of the snapshot */
    RealType currentTime;     /**< current time */
    Mat3x3d hmat;             /**< axes of the periodic box in matrix form */
    Mat3x3d invHmat;          /**< the inverse of the Hmat matrix */
    bool orthoRhombic;        /**< is this an orthorhombic periodic box? */
    RealType volume;          /**< total volume of this frame */
    RealType pressure;        /**< pressure of this frame */
    RealType totalEnergy;     /**< total energy of this frame */
    RealType kineticEnergy;   /**< kinetic energy of this frame */
    RealType potentialEnergy; /**< potential energy of this frame */
    RealType temperature;     /**< temperature of this frame */
    RealType chi;             /**< thermostat velocity */
    RealType integralOfChiDt; /**< the actual thermostat */
    RealType electronicTemperature; /**< temperature of the electronic degrees of freedom */
    RealType chiQ;            /**< fluctuating charge thermostat velocity */
    RealType integralOfChiQDt; /**< the actual fluctuating charge thermostat */
    Mat3x3d eta;              /**< barostat matrix */
    Vector3d COM;             /**< location of center of mass */
    Vector3d COMvel;          /**< system center of mass velocity */
    Vector3d COMw;            /**< system center of mass angular velocity */
    Mat3x3d stressTensor;     /**< stress tensor */
    Mat3x3d pressureTensor;   /**< pressure tensor */
    Vector3d systemDipole;    /**< total system dipole moment */
    Vector3d conductiveHeatFlux; /**< heat flux vector (conductive only) */
  };


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
                                  orthoTolerance_(1e-6), hasCOM_(false), hasVolume_(false){
      
      frameData.id = -1;                   
      frameData.currentTime = 0;     
      frameData.hmat = Mat3x3d(0.0);             
      frameData.invHmat = Mat3x3d(0.0);          
      frameData.orthoRhombic = false;        
      frameData.volume = 0.0;          
      frameData.pressure = 0.0;        
      frameData.totalEnergy = 0.0;     
      frameData.kineticEnergy = 0.0;   
      frameData.potentialEnergy = 0.0; 
      frameData.temperature = 0.0;     
      frameData.chi = 0.0;             
      frameData.integralOfChiDt = 0.0; 
      frameData.electronicTemperature = 0.0;
      frameData.chiQ = 0.0;            
      frameData.integralOfChiQDt = 0.0; 
      frameData.eta = Mat3x3d(0.0);              
      frameData.COM = V3Zero;             
      frameData.COMvel = V3Zero;          
      frameData.COMw = V3Zero;            
      frameData.stressTensor = Mat3x3d(0.0);              
      frameData.pressureTensor = Mat3x3d(0.0);   
      frameData.systemDipole = Vector3d(0.0);            
      frameData.conductiveHeatFlux = Vector3d(0.0, 0.0, 0.0);
    }

    Snapshot(int nAtoms, int nRigidbodies, int nCutoffGroups, 
             int storageLayout) : atomData(nAtoms, storageLayout), 
                                  rigidbodyData(nRigidbodies, storageLayout),
                                  cgData(nCutoffGroups, DataStorage::dslPosition),
                                  orthoTolerance_(1e-6),
                                  hasCOM_(false),
                                  hasVolume_(false) {
      frameData.id = -1;                   
      frameData.currentTime = 0;     
      frameData.hmat = Mat3x3d(0.0);             
      frameData.invHmat = Mat3x3d(0.0);          
      frameData.orthoRhombic = false;        
      frameData.volume = 0.0;          
      frameData.pressure = 0.0;        
      frameData.totalEnergy = 0.0;     
      frameData.kineticEnergy = 0.0;   
      frameData.potentialEnergy = 0.0; 
      frameData.temperature = 0.0;     
      frameData.chi = 0.0;             
      frameData.integralOfChiDt = 0.0; 
      frameData.electronicTemperature = 0.0;
      frameData.chiQ = 0.0;            
      frameData.integralOfChiQDt = 0.0; 
      frameData.eta = Mat3x3d(0.0);              
      frameData.COM = V3Zero;             
      frameData.COMvel = V3Zero;          
      frameData.COMw = V3Zero;            
      frameData.stressTensor = Mat3x3d(0.0);              
      frameData.pressureTensor = Mat3x3d(0.0);   
      frameData.systemDipole = V3Zero;            
      frameData.conductiveHeatFlux = Vector3d(0.0, 0.0, 0.0);            
    }
    
    /** Returns the id of this Snapshot */
    int getID() {
      return frameData.id;
    }

    /** Sets the id of this Snapshot */
    void setID(int id) {
      frameData.id = id;
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
      return frameData.hmat;
    }

    /** Sets the H-Matrix */
    void setHmat(const Mat3x3d& m);
            
    RealType getVolume() {
      if (hasVolume_){
        return frameData.volume;
      }else{
        return frameData.hmat.determinant();
      }
    }

    void setVolume(RealType volume){
      hasVolume_=true;
      frameData.volume = volume;
    }

    /** Returns the inverse H-Matrix */
    Mat3x3d getInvHmat() {
      return frameData.invHmat;
    }

    /** Wrapping the vector according to periodic boundary condition*/
    void wrapVector(Vector3d& v);
    /** Scaling a vector to multiples of the periodic box */
    Vector3d scaleVector(Vector3d &v);


    Vector3d getCOM();
    Vector3d getCOMvel();
    Vector3d getCOMw();
            
    RealType getTime() {
      return frameData.currentTime;
    }

    void increaseTime(RealType dt) {
      setTime(getTime() + dt);
    }

    void setTime(RealType time) {
      frameData.currentTime =time;
      //time at statData is redundant
      statData[Stats::TIME] = frameData.currentTime;
    }

    RealType getChi() {
      return frameData.chi;
    }

    void setChi(RealType chi) {
      frameData.chi = chi;
    }

    RealType getIntegralOfChiDt() {
      return frameData.integralOfChiDt;
    }

    void setIntegralOfChiDt(RealType integralOfChiDt) {
      frameData.integralOfChiDt = integralOfChiDt;
    }
            
    RealType getChiElectronic() {
      return frameData.chiQ;
    }

    void setChiElectronic(RealType chiQ) {
      frameData.chiQ = chiQ;
    }

    RealType getIntegralOfChiElectronicDt() {
      return frameData.integralOfChiQDt;
    }

    void setIntegralOfChiElectronicDt(RealType integralOfChiQDt) {
      frameData.integralOfChiQDt = integralOfChiQDt;
    }
            

    void setOrthoTolerance(RealType orthoTolerance) {
      orthoTolerance_ = orthoTolerance;
    }

    Mat3x3d getEta() {
      return frameData.eta;
    }

    void setEta(const Mat3x3d& eta) {
      frameData.eta = eta;
    }

    Mat3x3d getStressTensor() {
      return frameData.stressTensor;
    }
        
    void setStressTensor(const Mat3x3d& stressTensor) {
      frameData.stressTensor = stressTensor;
    }

    Vector3d getConductiveHeatFlux() {
      return frameData.conductiveHeatFlux;
    }
        
    void setConductiveHeatFlux(const Vector3d& heatFlux) {
      frameData.conductiveHeatFlux = heatFlux;
    }

    bool hasCOM() {
      return hasCOM_;
    }

    void setCOMprops(const Vector3d& COM, const Vector3d& COMvel, const Vector3d& COMw) {
      frameData.COM = COM;
      frameData.COMvel = COMvel;
      frameData.COMw = COMw;
      hasCOM_ = true;
    }

    DataStorage atomData;
    DataStorage rigidbodyData;
    DataStorage cgData;
    FrameData frameData;
    Stats statData;

  private:
    RealType orthoTolerance_;
    bool hasCOM_;
    bool hasVolume_;    
  };

  typedef DataStorage (Snapshot::*DataStoragePointer); 
}
#endif //BRAINS_SNAPSHOT_HPP
