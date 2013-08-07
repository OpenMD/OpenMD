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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
  
/**
 * @file Snapshot.cpp
 * @author tlin
 * @date 11/11/2004
 * @version 1.0
 */

#include "brains/Snapshot.hpp"
#include "utils/NumericConstant.hpp"
#include "utils/simError.h"
#include "utils/Utility.hpp"
#include <cstdio>

namespace OpenMD {

  Snapshot::Snapshot(int nAtoms, int nRigidbodies, int nCutoffGroups) : 
    atomData(nAtoms), rigidbodyData(nRigidbodies),
    cgData(nCutoffGroups, DataStorage::dslPosition), 
    orthoTolerance_(1e-6) {
    
    frameData.id = -1;                   
    frameData.currentTime = 0;     
    frameData.hmat = Mat3x3d(0.0);             
    frameData.invHmat = Mat3x3d(0.0);          
    frameData.orthoRhombic = false;        
    frameData.bondPotential = 0.0;      
    frameData.bendPotential = 0.0;      
    frameData.torsionPotential = 0.0;   
    frameData.inversionPotential = 0.0; 
    frameData.lrPotentials = potVec(0.0);
    frameData.reciprocalPotential = 0.0;
    frameData.excludedPotentials = potVec(0.0); 
    frameData.restraintPotential = 0.0; 
    frameData.rawPotential = 0.0;   
    frameData.xyArea = 0.0;
    frameData.volume = 0.0;          
    frameData.thermostat = make_pair(0.0, 0.0);
    frameData.electronicThermostat = make_pair(0.0, 0.0);
    frameData.barostat = Mat3x3d(0.0);              
    frameData.stressTensor = Mat3x3d(0.0);              
    frameData.conductiveHeatFlux = Vector3d(0.0, 0.0, 0.0);

    clearDerivedProperties();
  }
  
  Snapshot::Snapshot(int nAtoms, int nRigidbodies, int nCutoffGroups, 
                     int storageLayout) : 
    atomData(nAtoms, storageLayout), 
    rigidbodyData(nRigidbodies, storageLayout),
    cgData(nCutoffGroups, DataStorage::dslPosition),
    orthoTolerance_(1e-6) {
    
    frameData.id = -1;                   
    frameData.currentTime = 0;     
    frameData.hmat = Mat3x3d(0.0);             
    frameData.invHmat = Mat3x3d(0.0);      
    frameData.bBox = Mat3x3d(0.0);             
    frameData.invBbox = Mat3x3d(0.0);
    frameData.orthoRhombic = false;        
    frameData.bondPotential = 0.0;      
    frameData.bendPotential = 0.0;      
    frameData.torsionPotential = 0.0;   
    frameData.inversionPotential = 0.0; 
    frameData.lrPotentials = potVec(0.0);
    frameData.reciprocalPotential = 0.0;
    frameData.excludedPotentials = potVec(0.0); 
    frameData.restraintPotential = 0.0; 
    frameData.rawPotential = 0.0;       
    frameData.xyArea = 0.0;
    frameData.volume = 0.0;          
    frameData.thermostat = make_pair(0.0, 0.0);
    frameData.electronicThermostat = make_pair(0.0, 0.0);
    frameData.barostat = Mat3x3d(0.0);              
    frameData.stressTensor = Mat3x3d(0.0);              
    frameData.conductiveHeatFlux = Vector3d(0.0, 0.0, 0.0);

    clearDerivedProperties();
  }

  void Snapshot::clearDerivedProperties() {
    frameData.totalEnergy = 0.0;     
    frameData.translationalKinetic = 0.0;   
    frameData.rotationalKinetic = 0.0;   
    frameData.kineticEnergy = 0.0;   
    frameData.potentialEnergy = 0.0; 
    frameData.shortRangePotential = 0.0;
    frameData.longRangePotential = 0.0; 
    frameData.pressure = 0.0;        
    frameData.temperature = 0.0;
    frameData.pressureTensor = Mat3x3d(0.0);   
    frameData.systemDipole = Vector3d(0.0);            
    frameData.convectiveHeatFlux = Vector3d(0.0, 0.0, 0.0);
    frameData.electronicTemperature = 0.0;
    frameData.COM = V3Zero;             
    frameData.COMvel = V3Zero;          
    frameData.COMw = V3Zero;  

    hasTotalEnergy = false;         
    hasTranslationalKineticEnergy = false;       
    hasRotationalKineticEnergy = false;       
    hasKineticEnergy = false;       
    hasShortRangePotential = false;
    hasLongRangePotential = false;
    hasPotentialEnergy = false;   
    hasXYarea = false;
    hasVolume = false;         
    hasPressure = false;       
    hasTemperature = false;    
    hasElectronicTemperature = false;
    hasCOM = false;
    hasCOMvel = false;
    hasCOMw = false;
    hasPressureTensor = false;    
    hasSystemDipole = false;      
    hasConvectiveHeatFlux = false;  
    hasInertiaTensor = false;
    hasGyrationalVolume = false;   
    hasHullVolume = false;
    hasConservedQuantity = false; 
    hasBoundingBox = false;
  }

  /** Returns the id of this Snapshot */
  int Snapshot::getID() {
    return frameData.id;
  }
  
  /** Sets the id of this Snapshot */
  void Snapshot::setID(int id) {
    frameData.id = id;
  }
  
  int Snapshot::getSize() {
    return atomData.getSize() + rigidbodyData.getSize();
  }
  
  /** Returns the number of atoms */
  int Snapshot::getNumberOfAtoms() {
    return atomData.getSize();
  }
  
  /** Returns the number of rigid bodies */
  int Snapshot::getNumberOfRigidBodies() {
    return rigidbodyData.getSize();
  }
  
  /** Returns the number of rigid bodies */
  int Snapshot::getNumberOfCutoffGroups() {
    return cgData.getSize();
  }
  
  /** Returns the H-Matrix */
  Mat3x3d Snapshot::getHmat() {
    return frameData.hmat;
  }

  /** Sets the H-Matrix */  
  void Snapshot::setHmat(const Mat3x3d& m) {
    hasVolume = false;
    frameData.hmat = m;
    frameData.invHmat = frameData.hmat.inverse();
    
    //determine whether the box is orthoTolerance or not
    bool oldOrthoRhombic = frameData.orthoRhombic;
    
    RealType smallDiag = fabs(frameData.hmat(0, 0));
    if(smallDiag > fabs(frameData.hmat(1, 1))) smallDiag = fabs(frameData.hmat(1, 1));
    if(smallDiag > fabs(frameData.hmat(2, 2))) smallDiag = fabs(frameData.hmat(2, 2));    
    RealType tol = smallDiag * orthoTolerance_;

    frameData.orthoRhombic = true;

    for (int i = 0; i < 3; i++ ) {
      for (int j = 0 ; j < 3; j++) {
	if (i != j) {
	  if (frameData.orthoRhombic) {
	    if ( fabs(frameData.hmat(i, j)) >= tol)
	      frameData.orthoRhombic = false;
	  }        
	}
      }
    }
    
    if( oldOrthoRhombic != frameData.orthoRhombic){
      
      // It is finally time to suppress these warnings once and for
      // all.  They were annoying and not very informative.

      // if( frameData.orthoRhombic ) {
      //   sprintf( painCave.errMsg,
      //   	 "OpenMD is switching from the default Non-Orthorhombic\n"
      //   	 "\tto the faster Orthorhombic periodic boundary computations.\n"
      //   	 "\tThis is usually a good thing, but if you want the\n"
      //   	 "\tNon-Orthorhombic computations, make the orthoBoxTolerance\n"
      //   	 "\tvariable ( currently set to %G ) smaller.\n",
      //   	 orthoTolerance_);
      //   painCave.severity = OPENMD_INFO;
      //   simError();
      // }
      // else {
      //   sprintf( painCave.errMsg,
      //   	 "OpenMD is switching from the faster Orthorhombic to the more\n"
      //   	 "\tflexible Non-Orthorhombic periodic boundary computations.\n"
      //   	 "\tThis is usually because the box has deformed under\n"
      //   	 "\tNPTf integration. If you want to live on the edge with\n"
      //   	 "\tthe Orthorhombic computations, make the orthoBoxTolerance\n"
      //   	 "\tvariable ( currently set to %G ) larger.\n",
      //   	 orthoTolerance_);
      //   painCave.severity = OPENMD_WARNING;
      //   simError();
      // }
    }    
  }
  
  /** Returns the inverse H-Matrix */
  Mat3x3d Snapshot::getInvHmat() {
    return frameData.invHmat;
  }

  /** Returns the Bounding Box */
  Mat3x3d Snapshot::getBoundingBox() {
    return frameData.bBox;
  }

  /** Sets the Bounding Box */  
  void Snapshot::setBoundingBox(const Mat3x3d& m) {
    frameData.bBox = m;
    frameData.invBbox = frameData.bBox.inverse();
    hasBoundingBox = true;
  }

  /** Returns the inverse Bounding Box */
  Mat3x3d Snapshot::getInvBoundingBox() {
    return frameData.invBbox;
  }

  RealType Snapshot::getXYarea() {
    if (!hasXYarea) {
      Vector3d x = frameData.hmat.getColumn(0);
      Vector3d y = frameData.hmat.getColumn(1);
      frameData.xyArea = cross(x,y).length();
      hasXYarea = true;
    }
    return frameData.xyArea;
  }

  RealType Snapshot::getVolume() {
    if (!hasVolume) {
      frameData.volume = frameData.hmat.determinant();
      hasVolume = true;
    }
    return frameData.volume;
  }

  void Snapshot::setVolume(RealType vol) {
    hasVolume = true;
    frameData.volume = vol;
  }


  /** Wrap a vector according to periodic boundary conditions */
  void Snapshot::wrapVector(Vector3d& pos) {
    
    if( !frameData.orthoRhombic ) {
      Vector3d scaled = frameData.invHmat * pos;
      for (int i = 0; i < 3; i++) {
        scaled[i] -= roundMe( scaled[i] );        
      }
      // calc the wrapped real coordinates from the wrapped scaled coordinates
      pos = frameData.hmat * scaled;
    } else {
      RealType scaled;
      for (int i=0; i<3; i++) {      
        scaled = pos[i] * frameData.invHmat(i,i);
        scaled -= roundMe( scaled );
        pos[i] = scaled * frameData.hmat(i,i);
      }
    }
  }

  /** Scaling a vector to multiples of the periodic box */
  inline Vector3d Snapshot::scaleVector(Vector3d& pos) {   
    
    Vector3d scaled;

    if( !frameData.orthoRhombic )
      scaled = frameData.invHmat * pos;
    else {
      // calc the scaled coordinates.
      for (int i=0; i<3; i++) 
        scaled[i] = pos[i] * frameData.invHmat(i, i);
    }

    return scaled;
  }

  void Snapshot::setCOM(const Vector3d& com) {
    frameData.COM = com;
    hasCOM = true;
  }
  
  void Snapshot::setCOMvel(const Vector3d& comVel) {
    frameData.COMvel = comVel;
    hasCOMvel = true;
  }
  
  void Snapshot::setCOMw(const Vector3d& comw) {
    frameData.COMw = comw;
    hasCOMw = true;
  }
  
  Vector3d Snapshot::getCOM() {
    return frameData.COM;
  }
  
  Vector3d Snapshot::getCOMvel() {
    return frameData.COMvel;
  }
  
  Vector3d Snapshot::getCOMw() {
    return frameData.COMw;
  }
  
  RealType Snapshot::getTime() {
    return frameData.currentTime;
  }
  
  void Snapshot::increaseTime(RealType dt) {
    setTime(getTime() + dt);
  }
  
  void Snapshot::setTime(RealType time) {
    frameData.currentTime = time;
  }
 
  void Snapshot::setBondPotential(RealType bp) {
    frameData.bondPotential = bp;
  }
  
  void Snapshot::setBendPotential(RealType bp) {
    frameData.bendPotential = bp;
  }
  
  void Snapshot::setTorsionPotential(RealType tp) {
    frameData.torsionPotential = tp;
  }
  
  void Snapshot::setInversionPotential(RealType ip) {
    frameData.inversionPotential = ip;
  }


  RealType Snapshot::getBondPotential() {
    return frameData.bondPotential;
  }
  RealType Snapshot::getBendPotential() {
    return frameData.bendPotential;
  }
  RealType Snapshot::getTorsionPotential() {
    return frameData.torsionPotential;
  }
  RealType Snapshot::getInversionPotential() {
    return frameData.inversionPotential;
  }

  RealType Snapshot::getShortRangePotential() {
    if (!hasShortRangePotential) {
      frameData.shortRangePotential = frameData.bondPotential;
      frameData.shortRangePotential += frameData.bendPotential;
      frameData.shortRangePotential += frameData.torsionPotential;
      frameData.shortRangePotential += frameData.inversionPotential;
      hasShortRangePotential = true;
    }
    return frameData.shortRangePotential;
  }

  void Snapshot::setReciprocalPotential(RealType rp){
    frameData.reciprocalPotential = rp;
  }

  RealType Snapshot::getReciprocalPotential() {
    return frameData.reciprocalPotential;
  }

  void Snapshot::setLongRangePotential(potVec lrPot) {
    frameData.lrPotentials = lrPot;
  }
    
  RealType Snapshot::getLongRangePotential() {
    if (!hasLongRangePotential) {
      for (int i = 0; i < N_INTERACTION_FAMILIES; i++) {
        frameData.longRangePotential += frameData.lrPotentials[i];
      }
      frameData.longRangePotential += frameData.reciprocalPotential;
      hasLongRangePotential = true;
    }   
    return frameData.longRangePotential;
  }

  potVec Snapshot::getLongRangePotentials() {
    return frameData.lrPotentials;
  }

  RealType Snapshot::getPotentialEnergy() {
    if (!hasPotentialEnergy) {
      frameData.potentialEnergy = this->getLongRangePotential();
      frameData.potentialEnergy += this->getShortRangePotential();
      hasPotentialEnergy = true;
    }
    return frameData.potentialEnergy;
  }
    
  void Snapshot::setExcludedPotentials(potVec exPot) {
    frameData.excludedPotentials = exPot;
  }

  potVec Snapshot::getExcludedPotentials() {
    return frameData.excludedPotentials;
  }
      
  void Snapshot::setRestraintPotential(RealType rp) {
    frameData.restraintPotential = rp;
  }
  
  RealType Snapshot::getRestraintPotential() {
    return frameData.restraintPotential;
  }
  
  void Snapshot::setRawPotential(RealType rp) {
    frameData.rawPotential = rp;
  }
  
  RealType Snapshot::getRawPotential() {
    return frameData.rawPotential;
  }

  RealType Snapshot::getTranslationalKineticEnergy() {
    return frameData.translationalKinetic;
  }

  RealType Snapshot::getRotationalKineticEnergy() {
    return frameData.rotationalKinetic;
  }

  RealType Snapshot::getKineticEnergy() {
    return frameData.kineticEnergy;
  }

  void Snapshot::setTranslationalKineticEnergy(RealType tke) {
    hasTranslationalKineticEnergy = true;
    frameData.translationalKinetic = tke;
  }

  void Snapshot::setRotationalKineticEnergy(RealType rke) {
    hasRotationalKineticEnergy = true;
    frameData.rotationalKinetic = rke;
  }

  void Snapshot::setKineticEnergy(RealType ke) {
    hasKineticEnergy = true;
    frameData.kineticEnergy = ke;
  }

  RealType Snapshot::getTotalEnergy() {
    return frameData.totalEnergy;
  }

  void Snapshot::setTotalEnergy(RealType te) {
    hasTotalEnergy = true;
    frameData.totalEnergy = te;
  }

  RealType Snapshot::getConservedQuantity() {
    return frameData.conservedQuantity;
  }

  void Snapshot::setConservedQuantity(RealType cq) {
    hasConservedQuantity = true;
    frameData.conservedQuantity = cq;
  }

  RealType Snapshot::getTemperature() {
    return frameData.temperature;
  }

  void Snapshot::setTemperature(RealType temp) {
    hasTemperature = true;
    frameData.temperature = temp;
  }

  RealType Snapshot::getElectronicTemperature() {
    return frameData.electronicTemperature;
  }

  void Snapshot::setElectronicTemperature(RealType eTemp) {
    hasElectronicTemperature = true;
    frameData.electronicTemperature = eTemp;
  }

  RealType Snapshot::getPressure() {
    return frameData.pressure;
  }

  void Snapshot::setPressure(RealType pressure) {
    hasPressure = true;
    frameData.pressure = pressure;
  }

  Mat3x3d Snapshot::getPressureTensor() {
    return frameData.pressureTensor;
  }


  void Snapshot::setPressureTensor(const Mat3x3d& pressureTensor) {
    hasPressureTensor = true;
    frameData.pressureTensor = pressureTensor;
  }

  void Snapshot::setStressTensor(const Mat3x3d& stressTensor) {
    frameData.stressTensor = stressTensor;
  }

  Mat3x3d  Snapshot::getStressTensor() {
    return frameData.stressTensor;
  }

  void Snapshot::setConductiveHeatFlux(const Vector3d& chf) {
    frameData.conductiveHeatFlux = chf;
  }

  Vector3d Snapshot::getConductiveHeatFlux() {
    return frameData.conductiveHeatFlux;
  }
  
  Vector3d Snapshot::getConvectiveHeatFlux() {
    return frameData.convectiveHeatFlux;
  }

  void Snapshot::setConvectiveHeatFlux(const Vector3d& chf) {    
    hasConvectiveHeatFlux = true;
    frameData.convectiveHeatFlux = chf;
  }

  Vector3d Snapshot::getHeatFlux() {
    // BE CAREFUL WITH UNITS
    return getConductiveHeatFlux() + getConvectiveHeatFlux();
  }

  Vector3d Snapshot::getSystemDipole() {
    return frameData.systemDipole;
  }

  void Snapshot::setSystemDipole(const Vector3d& bd) {    
    hasSystemDipole = true;
    frameData.systemDipole = bd;
  }

  void Snapshot::setThermostat(const pair<RealType, RealType>& thermostat) {
    frameData.thermostat = thermostat;
  }

  pair<RealType, RealType> Snapshot::getThermostat() {
    return frameData.thermostat;
  }

  void Snapshot::setElectronicThermostat(const pair<RealType, RealType>& eTherm) {
    frameData.electronicThermostat = eTherm;
  }

  pair<RealType, RealType> Snapshot::getElectronicThermostat() {
    return frameData.electronicThermostat;
  }
 
  void Snapshot::setBarostat(const Mat3x3d& barostat) {
    frameData.barostat = barostat;
  }

  Mat3x3d Snapshot::getBarostat() {
    return frameData.barostat;
  }

  void Snapshot::setInertiaTensor(const Mat3x3d& inertiaTensor) {
    frameData.inertiaTensor = inertiaTensor;
    hasInertiaTensor = true;
  }

  Mat3x3d Snapshot::getInertiaTensor() {
    return frameData.inertiaTensor;
  }

  void Snapshot::setGyrationalVolume(const RealType gyrationalVolume) {
    frameData.gyrationalVolume = gyrationalVolume;
    hasGyrationalVolume = true;
  }

  RealType Snapshot::getGyrationalVolume() {
    return frameData.gyrationalVolume;
  }

  void Snapshot::setHullVolume(const RealType hullVolume) {
    frameData.hullVolume = hullVolume;
    hasHullVolume = true;
  }

  RealType Snapshot::getHullVolume() {
    return frameData.hullVolume;
  }

  void Snapshot::setOrthoTolerance(RealType ot) {
    orthoTolerance_ = ot;
  }
}
