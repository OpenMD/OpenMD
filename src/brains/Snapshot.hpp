/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#ifndef BRAINS_SNAPSHOT_HPP
#define BRAINS_SNAPSHOT_HPP

#include <vector>

#include "brains/DataStorage.hpp"
#include "brains/Stats.hpp"
#include "nonbonded/NonBondedInteraction.hpp"

using namespace std;
namespace OpenMD {

  /**
   * Storage specific to the SPF-RNEMD method that allows for a new simulation
   *  to pick up where an old one left off.
   *
   * \note used primarily in a \c std::shared_ptr
   */
  struct SPFData {
    Vector3d pos {V3Zero}; /**< location to place a selected molecule */
    RealType lambda {0.0}; /**< how much of the molecule has been transferred */
    int globalID {-1};     /**< which molecule have we selected */

    /**
     * Reset member variables to their defaults. Prefer this to resetting a
     *  \c std::shared_ptr<SPFData> or allocating/deallocating more memory.
     */
    void clear() {
      pos      = V3Zero;
      lambda   = 0.0;
      globalID = -1;
    }
  };

  /**
   * FrameData is a structure for holding system-wide dynamic data
   * about the simulation.
   */
  struct FrameData {
    int id;               /**< identification number of the snapshot */
    RealType currentTime; /**< current time */
    Mat3x3d hmat;         /**< axes of the periodic box in matrix form */
    Mat3x3d invHmat;      /**< the inverse of the Hmat matrix */
    Mat3x3d bBox;         /**< axes of a bounding box in matrix form */
    Mat3x3d invBbox;      /**< the inverse of the bounding box */
    bool usePBC;          /**< are we using a periodic box? */
    bool orthoRhombic;    /**< is this an orthorhombic periodic box? */
    RealType totalEnergy; /**< total energy of this frame */
    RealType
        translationalKinetic; /**< translational kinetic energy of this frame */
    RealType rotationalKinetic; /**< rotational kinetic energy of this frame */
    RealType electronicKinetic; /**< electronic kinetic energy of this frame */
    RealType kineticEnergy;     /**< kinetic energy of this frame */
    RealType potentialEnergy;   /**< potential energy of this frame */
    RealType
        shortRangePotential; /**< short-range contributions to the potential*/
    RealType
        longRangePotential; /**< long-range contributions to the potential */
    RealType reciprocalPotential; /**< reciprocal-space contributions to the
                                     potential */
    RealType
        surfacePotential;   /**< surface-term contributions to the potential */
    RealType bondPotential; /**< bonded contribution to the potential */
    RealType bendPotential; /**< angle-bending contribution to the potential */
    RealType torsionPotential; /**< dihedral (torsion angle) contribution to the
                                  potential */
    RealType inversionPotential; /**< inversion (planarity) contribution to the
                                    potential */
    potVec lrPotentials;    /**< breakdown of long-range potentials by family */
    RealType selfPotential; /**< potential energy of self interactions */
    potVec selfPotentials;  /**< breakdown of self interactions by family */
    RealType
        excludedPotential; /**< potential energy excluded from atomic forces */
    potVec
        excludedPotentials; /**< breakdown of excluded potentials by family */
    RealType restraintPotential; /**< potential energy of restraints */
    RealType rawPotential; /**< unrestrained potential energy (when restraints
                              are applied) */
    potVec selectionPotentials; /**< potential of selected stuntDoubles */
    RealType xyArea;            /**< XY area of this frame */
    RealType xzArea;            /**< XZ area of this frame */
    RealType yzArea;            /**< YZ area of this frame */
    RealType volume;            /**< total volume of this frame */
    RealType pressure;          /**< pressure of this frame */
    RealType temperature;       /**< temperature of this frame */
    pair<RealType, RealType> thermostat; /**< thermostat variables */
    RealType electronicTemperature; /**< temperature of the electronic degrees
                                       of freedom */
    RealType netCharge;             /**< total net charge in the system */
    RealType chargeMomentum;        /**< total charge momentum in the system */
    pair<RealType, RealType>
        electronicThermostat;  /**< thermostat variables for electronic degrees
                                  of freedom */
    Mat3x3d barostat;          /**< barostat matrix */
    Vector3d COM;              /**< location of system center of mass */
    Vector3d COMvel;           /**< system center of mass velocity */
    Vector3d COMw;             /**< system center of mass angular velocity */
    Mat3x3d inertiaTensor;     /**< inertia tensor for entire system */
    RealType gyrationalVolume; /**< gyrational volume for entire system */
    RealType hullVolume;       /**< hull volume for entire system */
    Mat3x3d virialTensor;      /**< virial tensor */
    Mat3x3d pressureTensor;    /**< pressure tensor */
    Vector3d systemDipole;     /**< total system dipole moment */
    Mat3x3d systemQuadrupole;  /**< total system quadrupole moment */
    Vector3d conductiveHeatFlux; /**< heat flux vector (conductive only) */
    Vector3d convectiveHeatFlux; /**< heat flux vector (convective only) */
    RealType conservedQuantity;  /**< anything conserved by the integrator */
    std::shared_ptr<SPFData> spfData {
        nullptr}; /**< parameters for restarting an SPF simulation */
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
    Snapshot(int nAtoms, int nRigidbodies, int nCutoffGroups, bool usePBC);
    Snapshot(int nAtoms, int nRigidbodies, int nCutoffGroups,
             int atomStorageLayout, int rigidBodyStorageLayout,
             int cutoffGroupStorageLayout, bool usePBC);
    /** Returns the id of this Snapshot */
    int getID();
    /** Sets the id of this Snapshot */
    void setID(int id);

    /** sets the state of the computed properties to false */
    void clearDerivedProperties();

    int getSize();
    /** Returns the number of atoms */
    int getNumberOfAtoms();
    /** Returns the number of rigid bodies */
    int getNumberOfRigidBodies();
    /** Returns the number of rigid bodies */
    int getNumberOfCutoffGroups();
    /** Returns the number of bytes in a FrameData structure */
    static int getFrameDataSize();

    /** Returns the H-Matrix */
    Mat3x3d getHmat();
    /** Sets the H-Matrix */
    void setHmat(const Mat3x3d& m);
    /** Returns the inverse H-Matrix */
    Mat3x3d getInvHmat();

    /** Returns the Bounding Box */
    Mat3x3d getBoundingBox();
    /** Sets the Bounding Box */
    void setBoundingBox(const Mat3x3d& m);
    /** Returns the inverse Bounding Box*/
    Mat3x3d getInvBoundingBox();

    RealType getVolume();
    RealType getXYarea();
    RealType getXZarea();
    RealType getYZarea();
    void setVolume(const RealType vol);

    /** Wrapping the vector according to periodic boundary condition*/
    void wrapVector(Vector3d& v);

    /** Scaling a vector to multiples of the periodic box */
    Vector3d scaleVector(Vector3d& v);

    void setCOM(const Vector3d& com);
    void setCOMvel(const Vector3d& comVel);
    void setCOMw(const Vector3d& comw);

    Vector3d getCOM();
    Vector3d getCOMvel();
    Vector3d getCOMw();

    RealType getTime();
    void increaseTime(const RealType dt);
    void setTime(const RealType time);

    void setBondPotential(const RealType bp);
    void setBendPotential(const RealType bp);
    void setTorsionPotential(const RealType tp);
    void setInversionPotential(const RealType ip);
    RealType getBondPotential();
    RealType getBendPotential();
    RealType getTorsionPotential();
    RealType getInversionPotential();

    RealType getShortRangePotential();

    void setLongRangePotentials(const potVec lrPot);
    RealType getLongRangePotential();
    potVec getLongRangePotentials();

    void setReciprocalPotential(const RealType rp);
    RealType getReciprocalPotential();

    void setSurfacePotential(const RealType sp);
    RealType getSurfacePotential();

    void setSelfPotentials(const potVec sp);
    RealType getSelfPotential();
    potVec getSelfPotentials();

    void setExcludedPotentials(const potVec exPot);
    potVec getExcludedPotentials();
    RealType getExcludedPotential();

    void setRestraintPotential(const RealType rp);
    RealType getRestraintPotential();

    void setRawPotential(const RealType rp);
    RealType getRawPotential();

    void setSelectionPotentials(const potVec selPot);
    potVec getSelectionPotentials();

    RealType getPotentialEnergy();
    void setPotentialEnergy(const RealType pe);
    RealType getKineticEnergy();
    RealType getTranslationalKineticEnergy();
    RealType getRotationalKineticEnergy();
    RealType getElectronicKineticEnergy();
    void setKineticEnergy(const RealType ke);
    void setTranslationalKineticEnergy(const RealType tke);
    void setRotationalKineticEnergy(const RealType rke);
    void setElectronicKineticEnergy(const RealType eke);
    RealType getTotalEnergy();
    void setTotalEnergy(const RealType te);
    RealType getConservedQuantity();
    void setConservedQuantity(const RealType cq);
    RealType getTemperature();
    void setTemperature(const RealType temp);
    RealType getElectronicTemperature();
    void setElectronicTemperature(const RealType eTemp);
    RealType getNetCharge();
    void setNetCharge(const RealType nChg);
    RealType getChargeMomentum();
    void setChargeMomentum(const RealType cMom);

    RealType getPressure();
    void setPressure(const RealType pressure);

    Mat3x3d getPressureTensor();
    void setPressureTensor(const Mat3x3d& pressureTensor);

    Mat3x3d getVirialTensor();
    void setVirialTensor(const Mat3x3d& virialTensor);

    Vector3d getConductiveHeatFlux();
    void setConductiveHeatFlux(const Vector3d& chf);

    Vector3d getConvectiveHeatFlux();
    void setConvectiveHeatFlux(const Vector3d& chf);

    Vector3d getHeatFlux();

    Vector3d getSystemDipole();
    void setSystemDipole(const Vector3d& bd);

    Mat3x3d getSystemQuadrupole();
    void setSystemQuadrupole(const Mat3x3d& bq);

    pair<RealType, RealType> getThermostat();
    void setThermostat(const pair<RealType, RealType>& thermostat);

    pair<RealType, RealType> getElectronicThermostat();
    void setElectronicThermostat(const pair<RealType, RealType>& eThermostat);

    Mat3x3d getBarostat();
    void setBarostat(const Mat3x3d& barostat);

    std::shared_ptr<SPFData> getSPFData();
    void setSPFData(std::shared_ptr<SPFData> data);

    Mat3x3d getInertiaTensor();
    void setInertiaTensor(const Mat3x3d& inertiaTensor);

    RealType getGyrationalVolume();
    void setGyrationalVolume(const RealType gv);

    RealType getHullVolume();
    void setHullVolume(const RealType hv);

    void setOrthoTolerance(RealType orthoTolerance);

    DataStorage atomData;
    DataStorage rigidbodyData;
    DataStorage cgData;
    FrameData frameData;

    bool hasTotalEnergy;
    bool hasTranslationalKineticEnergy;
    bool hasRotationalKineticEnergy;
    bool hasElectronicKineticEnergy;
    bool hasKineticEnergy;
    bool hasShortRangePotential;
    bool hasLongRangePotential;
    bool hasExcludedPotential;
    bool hasSelfPotential;
    bool hasPotentialEnergy;
    bool hasXYarea;
    bool hasXZarea;
    bool hasYZarea;
    bool hasVolume;
    bool hasPressure;
    bool hasTemperature;
    bool hasElectronicTemperature;
    bool hasNetCharge;
    bool hasChargeMomentum;
    bool hasCOM;
    bool hasCOMvel;
    bool hasCOMw;
    bool hasPressureTensor;
    bool hasSystemDipole;
    bool hasSystemQuadrupole;
    bool hasConvectiveHeatFlux;
    bool hasInertiaTensor;
    bool hasGyrationalVolume;
    bool hasHullVolume;
    bool hasBoundingBox;

  private:
    RealType orthoTolerance_;
  };

  typedef DataStorage(Snapshot::*DataStoragePointer);

  // translation from typedef into using:
  // using DataStoragePointer = DataStorage Snapshot::*;
}  // namespace OpenMD

#endif  // BRAINS_SNAPSHOT_HPP
