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
  
#ifndef BRAINS_SNAPSHOT_HPP
#define BRAINS_SNAPSHOT_HPP

#include <vector>

#include "brains/DataStorage.hpp"
#include "nonbonded/NonBondedInteraction.hpp"
#include "brains/Stats.hpp"

using namespace std;
namespace OpenMD{

  /**
   * FrameData is a structure for holding system-wide dynamic data
   * about the simulation.
   */
  
  struct FrameData {
    int id;                       /**< identification number of the snapshot */
    RealType currentTime;         /**< current time */
    Mat3x3d  hmat;                /**< axes of the periodic box in matrix form */
    Mat3x3d  invHmat;             /**< the inverse of the Hmat matrix */
    bool     orthoRhombic;        /**< is this an orthorhombic periodic box? */
    RealType totalEnergy;         /**< total energy of this frame */
    RealType translationalKinetic; /**< translational kinetic energy of this frame */
    RealType rotationalKinetic;   /**< rotational kinetic energy of this frame */
    RealType kineticEnergy;       /**< kinetic energy of this frame */
    RealType potentialEnergy;     /**< potential energy of this frame */
    RealType shortRangePotential; /**< short-range contributions to the potential*/
    RealType longRangePotential;  /**< long-range contributions to the potential */
    RealType bondPotential;       /**< bonded contribution to the potential */
    RealType bendPotential;       /**< angle-bending contribution to the potential */
    RealType torsionPotential;    /**< dihedral (torsion angle) contribution to the potential */
    RealType inversionPotential;  /**< inversion (planarity) contribution to the potential */
    potVec   lrPotentials;        /**< breakdown of long-range potentials by family */
    potVec   excludedPotentials;  /**< breakdown of excluded potentials by family */
    RealType restraintPotential;  /**< potential energy of restraints */
    RealType rawPotential;        /**< unrestrained potential energy (when restraints are applied) */
    RealType volume;              /**< total volume of this frame */
    RealType pressure;            /**< pressure of this frame */
    RealType temperature;         /**< temperature of this frame */
    pair<RealType, RealType> thermostat;    /**< thermostat variables */
    RealType electronicTemperature; /**< temperature of the electronic degrees of freedom */
    pair<RealType, RealType> electronicThermostat; /**< thermostat variables for electronic degrees of freedom */
    Mat3x3d  barostat;            /**< barostat matrix */
    Vector3d COM;                 /**< location of system center of mass */
    Vector3d COMvel;              /**< system center of mass velocity */
    Vector3d COMw;                /**< system center of mass angular velocity */
    Mat3x3d  inertiaTensor;       /**< inertia tensor for entire system */
    RealType gyrationalVolume;    /**< gyrational volume for entire system */
    RealType hullVolume;          /**< hull volume for entire system */
    Mat3x3d  stressTensor;        /**< stress tensor */
    Mat3x3d  pressureTensor;      /**< pressure tensor */
    Vector3d systemDipole;        /**< total system dipole moment */
    Vector3d conductiveHeatFlux;  /**< heat flux vector (conductive only) */
    Vector3d convectiveHeatFlux;  /**< heat flux vector (convective only) */
    RealType conservedQuantity;   /**< anything conserved by the integrator */
  };


  /**
   * @class Snapshot 
   * @brief The Snapshot class is a repository storing dynamic data during a
   * Simulation.  Every Snapshot contains FrameData (for global information)
   * as well as DataStorage (one for Atoms, one for RigidBodies, and one for 
   * CutoffGroups).
   */
  class Snapshot {

  public:            
    Snapshot(int nAtoms, int nRigidbodies, int nCutoffGroups);
    Snapshot(int nAtoms, int nRigidbodies, int nCutoffGroups, int storageLayout);    
    /** Returns the id of this Snapshot */
    int      getID();
    /** Sets the id of this Snapshot */
    void     setID(int id);

    /** sets the state of the computed properties to false */
    void     clearDerivedProperties();

    int      getSize();
    /** Returns the number of atoms */
    int      getNumberOfAtoms();
    /** Returns the number of rigid bodies */
    int      getNumberOfRigidBodies();
    /** Returns the number of rigid bodies */
    int      getNumberOfCutoffGroups();

    /** Returns the H-Matrix */
    Mat3x3d  getHmat();
    /** Sets the H-Matrix */
    void     setHmat(const Mat3x3d& m);
    /** Returns the inverse H-Matrix */
    Mat3x3d  getInvHmat();
            
    RealType getVolume();
    void     setVolume(const RealType vol);

    /** Wrapping the vector according to periodic boundary condition*/
    void     wrapVector(Vector3d& v);

    /** Scaling a vector to multiples of the periodic box */
    Vector3d scaleVector(Vector3d &v);

    void     setCOM(const Vector3d &com);
    void     setCOMvel(const Vector3d &comVel);
    void     setCOMw(const Vector3d &comw);

    Vector3d getCOM();
    Vector3d getCOMvel();
    Vector3d getCOMw();
            
    RealType getTime();
    void     increaseTime(const RealType dt);
    void     setTime(const RealType time);

    void     setBondPotential(const RealType bp);
    void     setBendPotential(const RealType bp);
    void     setTorsionPotential(const RealType tp);
    void     setInversionPotential(const RealType ip);
    RealType getBondPotential();
    RealType getBendPotential();
    RealType getTorsionPotential();
    RealType getInversionPotential();

    RealType getShortRangePotential();

    void     setLongRangePotential(const potVec lrPot);
    RealType getLongRangePotential();
    potVec   getLongRangePotentials();

    void     setExcludedPotentials(const potVec exPot);
    potVec   getExcludedPotentials();
   
    void     setRestraintPotential(const RealType rp);
    RealType getRestraintPotential();

    void     setRawPotential(const RealType rp);
    RealType getRawPotential();

    RealType getPotentialEnergy();
    RealType getKineticEnergy();
    RealType getTranslationalKineticEnergy();
    RealType getRotationalKineticEnergy();
    void     setKineticEnergy(const RealType ke);
    void     setTranslationalKineticEnergy(const RealType tke);
    void     setRotationalKineticEnergy(const RealType rke);
    RealType getTotalEnergy();
    void     setTotalEnergy(const RealType te);
    RealType getConservedQuantity();
    void     setConservedQuantity(const RealType cq);
    RealType getTemperature();
    void     setTemperature(const RealType temp);
    RealType getElectronicTemperature();
    void     setElectronicTemperature(const RealType eTemp);
    RealType getPressure();
    void     setPressure(const RealType pressure);

    Mat3x3d  getPressureTensor();
    void     setPressureTensor(const Mat3x3d& pressureTensor);

    Mat3x3d  getStressTensor();
    void     setStressTensor(const Mat3x3d& stressTensor);

    Vector3d getConductiveHeatFlux();
    void     setConductiveHeatFlux(const Vector3d& chf);

    Vector3d getConvectiveHeatFlux();
    void     setConvectiveHeatFlux(const Vector3d& chf);

    Vector3d getHeatFlux();
    
    Vector3d getSystemDipole();
    void     setSystemDipole(const Vector3d& bd);

    pair<RealType, RealType> getThermostat();
    void setThermostat(const pair<RealType, RealType>& thermostat);

    pair<RealType, RealType> getElectronicThermostat();
    void setElectronicThermostat(const pair<RealType, RealType>& eThermostat);
            
    Mat3x3d  getBarostat();
    void     setBarostat(const Mat3x3d& barostat);

    Mat3x3d  getInertiaTensor();
    void     setInertiaTensor(const Mat3x3d& inertiaTensor);

    RealType getGyrationalVolume();
    void     setGyrationalVolume(const RealType gv);

    RealType getHullVolume();
    void     setHullVolume(const RealType hv);
    
    void     setOrthoTolerance(RealType orthoTolerance);

    DataStorage atomData;
    DataStorage rigidbodyData;
    DataStorage cgData;
    FrameData   frameData;

    bool hasTotalEnergy;         
    bool hasTranslationalKineticEnergy;    
    bool hasRotationalKineticEnergy;    
    bool hasKineticEnergy;    
    bool hasShortRangePotential;
    bool hasLongRangePotential;
    bool hasPotentialEnergy;     
    bool hasVolume;         
    bool hasPressure;       
    bool hasTemperature;    
    bool hasElectronicTemperature;
    bool hasCOM;             
    bool hasCOMvel;
    bool hasCOMw;
    bool hasPressureTensor;    
    bool hasSystemDipole;    
    bool hasConvectiveHeatFlux;
    bool hasInertiaTensor;
    bool hasGyrationalVolume;
    bool hasHullVolume;
    bool hasConservedQuantity;

  private:
    RealType orthoTolerance_;
    
  };

  typedef DataStorage (Snapshot::*DataStoragePointer); 
}
#endif //BRAINS_SNAPSHOT_HPP
