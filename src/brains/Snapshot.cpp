/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

/**
 * @file Snapshot.cpp
 * @author tlin
 * @date 11/11/2004
 * @version 1.0
 */

#include "brains/Snapshot.hpp"

#include <cstdio>

#include "utils/Utility.hpp"
#include "utils/simError.h"

namespace OpenMD {

  Snapshot::Snapshot(int nAtoms, int nRigidbodies, int nCutoffGroups,
                     bool usePBC) :
      atomData(nAtoms),
      rigidbodyData(nRigidbodies),
      cgData(nCutoffGroups, DataStorage::dslPosition), orthoTolerance_(1e-6) {
    frameData.id                   = -1;
    frameData.currentTime          = 0;
    frameData.hmat                 = Mat3x3d(0.0);
    frameData.invHmat              = Mat3x3d(0.0);
    frameData.orthoRhombic         = false;
    frameData.usePBC               = usePBC;
    frameData.bondPotential        = 0.0;
    frameData.bendPotential        = 0.0;
    frameData.torsionPotential     = 0.0;
    frameData.inversionPotential   = 0.0;
    frameData.lrPotentials         = potVec(0.0);
    frameData.surfacePotential     = 0.0;
    frameData.reciprocalPotential  = 0.0;
    frameData.selfPotentials       = potVec(0.0);
    frameData.excludedPotentials   = potVec(0.0);
    frameData.restraintPotential   = 0.0;
    frameData.rawPotential         = 0.0;
    frameData.xyArea               = 0.0;
    frameData.xzArea               = 0.0;
    frameData.yzArea               = 0.0;
    frameData.volume               = 0.0;
    frameData.thermostat           = make_pair(0.0, 0.0);
    frameData.electronicThermostat = make_pair(0.0, 0.0);
    frameData.barostat             = Mat3x3d(0.0);
    frameData.virialTensor         = Mat3x3d(0.0);
    frameData.conductiveHeatFlux   = Vector3d(0.0, 0.0, 0.0);
    clearDerivedProperties();
  }

  Snapshot::Snapshot(int nAtoms, int nRigidbodies, int nCutoffGroups,
                     int storageLayout, bool usePBC) :
      atomData(nAtoms, storageLayout),
      rigidbodyData(nRigidbodies, storageLayout),
      cgData(nCutoffGroups, DataStorage::dslPosition), orthoTolerance_(1e-6) {
    frameData.id                   = -1;
    frameData.currentTime          = 0;
    frameData.hmat                 = Mat3x3d(0.0);
    frameData.invHmat              = Mat3x3d(0.0);
    frameData.bBox                 = Mat3x3d(0.0);
    frameData.invBbox              = Mat3x3d(0.0);
    frameData.orthoRhombic         = false;
    frameData.usePBC               = usePBC;
    frameData.bondPotential        = 0.0;
    frameData.bendPotential        = 0.0;
    frameData.torsionPotential     = 0.0;
    frameData.inversionPotential   = 0.0;
    frameData.lrPotentials         = potVec(0.0);
    frameData.surfacePotential     = 0.0;
    frameData.reciprocalPotential  = 0.0;
    frameData.selfPotentials       = potVec(0.0);
    frameData.excludedPotentials   = potVec(0.0);
    frameData.restraintPotential   = 0.0;
    frameData.rawPotential         = 0.0;
    frameData.xyArea               = 0.0;
    frameData.xzArea               = 0.0;
    frameData.yzArea               = 0.0;
    frameData.volume               = 0.0;
    frameData.thermostat           = make_pair(0.0, 0.0);
    frameData.electronicThermostat = make_pair(0.0, 0.0);
    frameData.barostat             = Mat3x3d(0.0);
    frameData.virialTensor         = Mat3x3d(0.0);
    frameData.conductiveHeatFlux   = Vector3d(0.0, 0.0, 0.0);

    clearDerivedProperties();
  }

  void Snapshot::clearDerivedProperties() {
    frameData.totalEnergy           = 0.0;
    frameData.translationalKinetic  = 0.0;
    frameData.rotationalKinetic     = 0.0;
    frameData.electronicKinetic     = 0.0;
    frameData.kineticEnergy         = 0.0;
    frameData.potentialEnergy       = 0.0;
    frameData.shortRangePotential   = 0.0;
    frameData.longRangePotential    = 0.0;
    frameData.excludedPotential     = 0.0;
    frameData.selfPotential         = 0.0;
    frameData.pressure              = 0.0;
    frameData.temperature           = 0.0;
    frameData.pressureTensor        = Mat3x3d(0.0);
    frameData.systemDipole          = Vector3d(0.0);
    frameData.systemQuadrupole      = Mat3x3d(0.0);
    frameData.convectiveHeatFlux    = Vector3d(0.0, 0.0, 0.0);
    frameData.electronicTemperature = 0.0;
    frameData.netCharge             = 0.0;
    frameData.chargeMomentum        = 0.0;
    frameData.COM                   = V3Zero;
    frameData.COMvel                = V3Zero;
    frameData.COMw                  = V3Zero;

    hasTotalEnergy                = false;
    hasTranslationalKineticEnergy = false;
    hasRotationalKineticEnergy    = false;
    hasElectronicKineticEnergy    = false;
    hasKineticEnergy              = false;
    hasShortRangePotential        = false;
    hasLongRangePotential         = false;
    hasExcludedPotential          = false;
    hasSelfPotential              = false;
    hasPotentialEnergy            = false;
    hasXYarea                     = false;
    hasXZarea                     = false;
    hasYZarea                     = false;
    hasVolume                     = false;
    hasPressure                   = false;
    hasTemperature                = false;
    hasElectronicTemperature      = false;
    hasNetCharge                  = false;
    hasChargeMomentum             = false;
    hasCOM                        = false;
    hasCOMvel                     = false;
    hasCOMw                       = false;
    hasPressureTensor             = false;
    hasSystemDipole               = false;
    hasSystemQuadrupole           = false;
    hasConvectiveHeatFlux         = false;
    hasInertiaTensor              = false;
    hasGyrationalVolume           = false;
    hasHullVolume                 = false;
    hasBoundingBox                = false;
  }

  /** Returns the id of this Snapshot */
  int Snapshot::getID() { return frameData.id; }

  /** Sets the id of this Snapshot */
  void Snapshot::setID(int id) { frameData.id = id; }

  int Snapshot::getSize() {
    return atomData.getSize() + rigidbodyData.getSize();
  }

  /** Returns the number of atoms */
  int Snapshot::getNumberOfAtoms() { return atomData.getSize(); }

  /** Returns the number of rigid bodies */
  int Snapshot::getNumberOfRigidBodies() { return rigidbodyData.getSize(); }

  /** Returns the number of rigid bodies */
  int Snapshot::getNumberOfCutoffGroups() { return cgData.getSize(); }

  /** Returns the number of bytes in a FrameData structure */
  int Snapshot::getFrameDataSize() { return sizeof(FrameData); }

  /** Returns the H-Matrix */
  Mat3x3d Snapshot::getHmat() { return frameData.hmat; }

  /** Sets the H-Matrix */
  void Snapshot::setHmat(const Mat3x3d& m) {
    hasVolume         = false;
    frameData.hmat    = m;
    frameData.invHmat = frameData.hmat.inverse();

    // determine whether the box is orthoTolerance or not
    bool oldOrthoRhombic = frameData.orthoRhombic;

    RealType smallDiag = fabs(frameData.hmat(0, 0));
    if (smallDiag > fabs(frameData.hmat(1, 1)))
      smallDiag = fabs(frameData.hmat(1, 1));
    if (smallDiag > fabs(frameData.hmat(2, 2)))
      smallDiag = fabs(frameData.hmat(2, 2));
    RealType tol = smallDiag * orthoTolerance_;

    frameData.orthoRhombic = true;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (i != j) {
          if (frameData.orthoRhombic) {
            if (fabs(frameData.hmat(i, j)) >= tol)
              frameData.orthoRhombic = false;
          }
        }
      }
    }

    if (oldOrthoRhombic != frameData.orthoRhombic) {
      // It is finally time to suppress these warnings once and for
      // all.  They were annoying and not very informative.

      // if( frameData.orthoRhombic ) {
      //   sprintf( painCave.errMsg,
      //   	 "OpenMD is switching from the default Non-Orthorhombic\n"
      //   	 "\tto the faster Orthorhombic periodic boundary
      //   computations.\n"
      //   	 "\tThis is usually a good thing, but if you want the\n"
      //   	 "\tNon-Orthorhombic computations, make the orthoBoxTolerance\n"
      //   	 "\tvariable ( currently set to %G ) smaller.\n",
      //   	 orthoTolerance_);
      //   painCave.severity = OPENMD_INFO;
      //   simError();
      // }
      // else {
      //   sprintf( painCave.errMsg,
      //   	 "OpenMD is switching from the faster Orthorhombic to the
      //   more\n"
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
  Mat3x3d Snapshot::getInvHmat() { return frameData.invHmat; }

  /** Returns the Bounding Box */
  Mat3x3d Snapshot::getBoundingBox() { return frameData.bBox; }

  /** Sets the Bounding Box */
  void Snapshot::setBoundingBox(const Mat3x3d& m) {
    frameData.bBox    = m;
    frameData.invBbox = frameData.bBox.inverse();
    hasBoundingBox    = true;
  }

  /** Returns the inverse Bounding Box */
  Mat3x3d Snapshot::getInvBoundingBox() { return frameData.invBbox; }

  RealType Snapshot::getXYarea() {
    if (!hasXYarea) {
      Vector3d x       = frameData.hmat.getColumn(0);
      Vector3d y       = frameData.hmat.getColumn(1);
      frameData.xyArea = cross(x, y).length();
      hasXYarea        = true;
    }
    return frameData.xyArea;
  }

  RealType Snapshot::getXZarea() {
    if (!hasXZarea) {
      Vector3d x       = frameData.hmat.getColumn(0);
      Vector3d z       = frameData.hmat.getColumn(2);
      frameData.xzArea = cross(x, z).length();
      hasXZarea        = true;
    }
    return frameData.xzArea;
  }

  RealType Snapshot::getYZarea() {
    if (!hasYZarea) {
      Vector3d y       = frameData.hmat.getColumn(1);
      Vector3d z       = frameData.hmat.getColumn(2);
      frameData.yzArea = cross(y, z).length();
      hasYZarea        = true;
    }
    return frameData.yzArea;
  }

  RealType Snapshot::getVolume() {
    if (!hasVolume) {
      frameData.volume = frameData.hmat.determinant();
      hasVolume        = true;
    }
    return frameData.volume;
  }

  void Snapshot::setVolume(RealType vol) {
    hasVolume        = true;
    frameData.volume = vol;
  }

  /** Wrap a vector according to periodic boundary conditions */
  void Snapshot::wrapVector(Vector3d& pos) {
    if (!frameData.usePBC) return;

    if (!frameData.orthoRhombic) {
      Vector3d scaled = frameData.invHmat * pos;
      for (int i = 0; i < 3; i++) {
        scaled[i] -= roundMe(scaled[i]);
      }
      // calc the wrapped real coordinates from the wrapped scaled coordinates
      pos = frameData.hmat * scaled;
    } else {
      RealType scaled;
      for (int i = 0; i < 3; i++) {
        scaled = pos[i] * frameData.invHmat(i, i);
        scaled -= roundMe(scaled);
        pos[i] = scaled * frameData.hmat(i, i);
      }
    }
  }

  /** Scaling a vector to multiples of the periodic box */
  inline Vector3d Snapshot::scaleVector(Vector3d& pos) {
    Vector3d scaled;

    if (!frameData.orthoRhombic)
      scaled = frameData.invHmat * pos;
    else {
      // calc the scaled coordinates.
      for (int i = 0; i < 3; i++)
        scaled[i] = pos[i] * frameData.invHmat(i, i);
    }

    return scaled;
  }

  void Snapshot::setCOM(const Vector3d& com) {
    frameData.COM = com;
    hasCOM        = true;
  }

  void Snapshot::setCOMvel(const Vector3d& comVel) {
    frameData.COMvel = comVel;
    hasCOMvel        = true;
  }

  void Snapshot::setCOMw(const Vector3d& comw) {
    frameData.COMw = comw;
    hasCOMw        = true;
  }

  Vector3d Snapshot::getCOM() { return frameData.COM; }

  Vector3d Snapshot::getCOMvel() { return frameData.COMvel; }

  Vector3d Snapshot::getCOMw() { return frameData.COMw; }

  RealType Snapshot::getTime() { return frameData.currentTime; }

  void Snapshot::increaseTime(RealType dt) { setTime(getTime() + dt); }

  void Snapshot::setTime(RealType time) { frameData.currentTime = time; }

  void Snapshot::setBondPotential(RealType bp) {
    frameData.bondPotential = bp;
    hasShortRangePotential  = false;
    hasPotentialEnergy      = false;
    hasTotalEnergy          = false;
  }

  void Snapshot::setBendPotential(RealType bp) {
    frameData.bendPotential = bp;
    hasShortRangePotential  = false;
    hasPotentialEnergy      = false;
    hasTotalEnergy          = false;
  }

  void Snapshot::setTorsionPotential(RealType tp) {
    frameData.torsionPotential = tp;
    hasShortRangePotential     = false;
    hasPotentialEnergy         = false;
    hasTotalEnergy             = false;
  }

  void Snapshot::setInversionPotential(RealType ip) {
    frameData.inversionPotential = ip;
    hasShortRangePotential       = false;
    hasPotentialEnergy           = false;
    hasTotalEnergy               = false;
  }

  RealType Snapshot::getBondPotential() { return frameData.bondPotential; }

  RealType Snapshot::getBendPotential() { return frameData.bendPotential; }

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
      hasPotentialEnergy     = false;
      hasTotalEnergy         = false;
    }
    return frameData.shortRangePotential;
  }

  void Snapshot::setSurfacePotential(RealType sp) {
    frameData.surfacePotential = sp;
    hasLongRangePotential      = false;
    hasPotentialEnergy         = false;
    hasTotalEnergy             = false;
  }

  RealType Snapshot::getSurfacePotential() {
    return frameData.surfacePotential;
  }

  void Snapshot::setReciprocalPotential(RealType rp) {
    frameData.reciprocalPotential = rp;
    hasLongRangePotential         = false;
    hasPotentialEnergy            = false;
  }

  RealType Snapshot::getReciprocalPotential() {
    return frameData.reciprocalPotential;
  }

  void Snapshot::setSelfPotentials(potVec sp) {
    frameData.selfPotentials = sp;
    hasSelfPotential         = false;
    hasPotentialEnergy       = false;
    hasTotalEnergy           = false;
  }

  potVec Snapshot::getSelfPotentials() { return frameData.selfPotentials; }

  RealType Snapshot::getSelfPotential() {
    if (!hasSelfPotential) {
      for (int i = 0; i < N_INTERACTION_FAMILIES; i++) {
        frameData.selfPotential += frameData.selfPotentials[i];
      }
      hasSelfPotential   = true;
      hasPotentialEnergy = false;
      hasTotalEnergy     = false;
    }
    return frameData.selfPotential;
  }

  void Snapshot::setLongRangePotentials(potVec lrPot) {
    frameData.lrPotentials = lrPot;
    hasLongRangePotential  = false;
    hasPotentialEnergy     = false;
    hasTotalEnergy         = false;
  }

  RealType Snapshot::getLongRangePotential() {
    if (!hasLongRangePotential) {
      for (int i = 0; i < N_INTERACTION_FAMILIES; i++) {
        frameData.longRangePotential += frameData.lrPotentials[i];
      }
      frameData.longRangePotential += frameData.reciprocalPotential;
      frameData.longRangePotential += frameData.surfacePotential;
      hasLongRangePotential = true;
      hasPotentialEnergy    = false;
      hasTotalEnergy        = false;
    }
    return frameData.longRangePotential;
  }

  potVec Snapshot::getLongRangePotentials() { return frameData.lrPotentials; }

  RealType Snapshot::getPotentialEnergy() {
    if (!hasPotentialEnergy) {
      frameData.potentialEnergy = this->getLongRangePotential();
      frameData.potentialEnergy += this->getShortRangePotential();
      frameData.potentialEnergy += this->getSelfPotential();
      frameData.potentialEnergy += this->getExcludedPotential();
      hasPotentialEnergy = true;
      hasTotalEnergy     = false;
    }
    return frameData.potentialEnergy;
  }

  void Snapshot::setPotentialEnergy(const RealType pe) {
    frameData.potentialEnergy = pe;
    hasPotentialEnergy        = true;
    hasTotalEnergy            = false;
  }

  void Snapshot::setExcludedPotentials(potVec exPot) {
    frameData.excludedPotentials = exPot;
    hasExcludedPotential         = false;
    hasPotentialEnergy           = false;
    hasTotalEnergy               = false;
  }

  potVec Snapshot::getExcludedPotentials() {
    return frameData.excludedPotentials;
  }

  RealType Snapshot::getExcludedPotential() {
    if (!hasExcludedPotential) {
      for (int i = 0; i < N_INTERACTION_FAMILIES; i++) {
        frameData.excludedPotential += frameData.excludedPotentials[i];
      }
      hasExcludedPotential = true;
      hasPotentialEnergy   = false;
      hasTotalEnergy       = false;
    }
    return frameData.excludedPotential;
  }

  void Snapshot::setRestraintPotential(RealType rp) {
    frameData.restraintPotential = rp;
  }

  RealType Snapshot::getRestraintPotential() {
    return frameData.restraintPotential;
  }

  void Snapshot::setRawPotential(RealType rp) { frameData.rawPotential = rp; }

  RealType Snapshot::getRawPotential() { return frameData.rawPotential; }

  void Snapshot::setSelectionPotentials(potVec selPot) {
    frameData.selectionPotentials = selPot;
  }

  potVec Snapshot::getSelectionPotentials() {
    return frameData.selectionPotentials;
  }

  RealType Snapshot::getTranslationalKineticEnergy() {
    return frameData.translationalKinetic;
  }

  RealType Snapshot::getRotationalKineticEnergy() {
    return frameData.rotationalKinetic;
  }

  RealType Snapshot::getElectronicKineticEnergy() {
    return frameData.electronicKinetic;
  }

  RealType Snapshot::getKineticEnergy() { return frameData.kineticEnergy; }

  void Snapshot::setTranslationalKineticEnergy(RealType tke) {
    frameData.translationalKinetic = tke;
    hasTranslationalKineticEnergy  = true;
    hasKineticEnergy               = false;
    hasTotalEnergy                 = false;
  }

  void Snapshot::setRotationalKineticEnergy(RealType rke) {
    frameData.rotationalKinetic = rke;
    hasRotationalKineticEnergy  = true;
    hasKineticEnergy            = false;
    hasTotalEnergy              = false;
  }

  void Snapshot::setElectronicKineticEnergy(RealType eke) {
    frameData.electronicKinetic = eke;
    hasElectronicKineticEnergy  = true;
    hasKineticEnergy            = false;
    hasTotalEnergy              = false;
  }

  void Snapshot::setKineticEnergy(RealType ke) {
    frameData.kineticEnergy = ke;
    hasKineticEnergy        = true;
    hasTotalEnergy          = false;
  }

  RealType Snapshot::getTotalEnergy() { return frameData.totalEnergy; }

  void Snapshot::setTotalEnergy(RealType te) {
    frameData.totalEnergy = te;
    hasTotalEnergy        = true;
  }

  RealType Snapshot::getConservedQuantity() {
    return frameData.conservedQuantity;
  }

  void Snapshot::setConservedQuantity(RealType cq) {
    frameData.conservedQuantity = cq;
  }

  RealType Snapshot::getTemperature() { return frameData.temperature; }

  void Snapshot::setTemperature(RealType temp) {
    hasTemperature        = true;
    frameData.temperature = temp;
  }

  RealType Snapshot::getElectronicTemperature() {
    return frameData.electronicTemperature;
  }

  void Snapshot::setElectronicTemperature(RealType eTemp) {
    hasElectronicTemperature        = true;
    frameData.electronicTemperature = eTemp;
  }

  RealType Snapshot::getNetCharge() { return frameData.netCharge; }

  void Snapshot::setNetCharge(RealType nChg) {
    hasNetCharge        = true;
    frameData.netCharge = nChg;
  }

  RealType Snapshot::getChargeMomentum() { return frameData.chargeMomentum; }

  void Snapshot::setChargeMomentum(RealType cMom) {
    hasChargeMomentum        = true;
    frameData.chargeMomentum = cMom;
  }

  RealType Snapshot::getPressure() { return frameData.pressure; }

  void Snapshot::setPressure(RealType pressure) {
    hasPressure        = true;
    frameData.pressure = pressure;
  }

  Mat3x3d Snapshot::getPressureTensor() { return frameData.pressureTensor; }

  void Snapshot::setPressureTensor(const Mat3x3d& pressureTensor) {
    hasPressureTensor        = true;
    frameData.pressureTensor = pressureTensor;
  }

  void Snapshot::setVirialTensor(const Mat3x3d& virialTensor) {
    frameData.virialTensor = virialTensor;
  }

  Mat3x3d Snapshot::getVirialTensor() { return frameData.virialTensor; }

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
    hasConvectiveHeatFlux        = true;
    frameData.convectiveHeatFlux = chf;
  }

  Vector3d Snapshot::getHeatFlux() {
    // BE CAREFUL WITH UNITS
    return getConductiveHeatFlux() + getConvectiveHeatFlux();
  }

  Vector3d Snapshot::getSystemDipole() { return frameData.systemDipole; }

  void Snapshot::setSystemDipole(const Vector3d& bd) {
    hasSystemDipole        = true;
    frameData.systemDipole = bd;
  }

  Mat3x3d Snapshot::getSystemQuadrupole() { return frameData.systemQuadrupole; }

  void Snapshot::setSystemQuadrupole(const Mat3x3d& bq) {
    hasSystemQuadrupole        = true;
    frameData.systemQuadrupole = bq;
  }

  void Snapshot::setThermostat(const pair<RealType, RealType>& thermostat) {
    frameData.thermostat = thermostat;
  }

  pair<RealType, RealType> Snapshot::getThermostat() {
    return frameData.thermostat;
  }

  void Snapshot::setElectronicThermostat(
      const pair<RealType, RealType>& eTherm) {
    frameData.electronicThermostat = eTherm;
  }

  pair<RealType, RealType> Snapshot::getElectronicThermostat() {
    return frameData.electronicThermostat;
  }

  void Snapshot::setBarostat(const Mat3x3d& barostat) {
    frameData.barostat = barostat;
  }

  Mat3x3d Snapshot::getBarostat() { return frameData.barostat; }

  void Snapshot::setInertiaTensor(const Mat3x3d& inertiaTensor) {
    frameData.inertiaTensor = inertiaTensor;
    hasInertiaTensor        = true;
  }

  Mat3x3d Snapshot::getInertiaTensor() { return frameData.inertiaTensor; }

  void Snapshot::setGyrationalVolume(const RealType gyrationalVolume) {
    frameData.gyrationalVolume = gyrationalVolume;
    hasGyrationalVolume        = true;
  }

  RealType Snapshot::getGyrationalVolume() {
    return frameData.gyrationalVolume;
  }

  void Snapshot::setHullVolume(const RealType hullVolume) {
    frameData.hullVolume = hullVolume;
    hasHullVolume        = true;
  }

  RealType Snapshot::getHullVolume() { return frameData.hullVolume; }

  void Snapshot::setOrthoTolerance(RealType ot) { orthoTolerance_ = ot; }
}  // namespace OpenMD
