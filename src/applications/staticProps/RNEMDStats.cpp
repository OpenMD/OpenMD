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

#include "applications/staticProps/RNEMDStats.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "applications/staticProps/SpatialStatistics.hpp"
#include "brains/DataStorage.hpp"
#include "brains/SimInfo.hpp"
#include "io/Globals.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector3.hpp"
#include "primitives/Atom.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/RigidBody.hpp"
#include "primitives/StuntDouble.hpp"
#include "rnemd/RNEMDParameters.hpp"
#include "types/AtomType.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/Accumulator.hpp"
#include "utils/AccumulatorView.hpp"
#include "utils/BaseAccumulator.hpp"
#include "utils/Constants.hpp"
#include "utils/StringUtils.hpp"

using namespace OpenMD::Utils;

namespace OpenMD {

  RNEMDZ::RNEMDZ(SimInfo* info, const std::string& filename,
                 const std::string& sele, int nzbins, int axis) :
      SlabStatistics(info, filename, sele, nzbins, axis) {
    setOutputName(getPrefix(filename) + ".rnemdZ");

    evaluator_.loadScriptString(sele);
    seleMan_.setSelectionSet(evaluator_.evaluate());
    AtomTypeSet osTypes = seleMan_.getSelectedAtomTypes();
    std::copy(osTypes.begin(), osTypes.end(), std::back_inserter(outputTypes_));

    // Pre-load the OutputData
    data_.resize(RNEMDZ::ENDINDEX);

    OutputData z;
    z.units        = "Angstroms";
    z.title        = axisLabel_;
    z.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      z.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());
    data_[Z] = std::move(z);

    OutputData temperature;
    temperature.units        = "K";
    temperature.title        = "Temperature";
    temperature.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      temperature.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());
    data_[TEMPERATURE] = std::move(temperature);

    OutputData velocity;
    velocity.units        = "angstroms/fs";
    velocity.title        = "Velocity";
    velocity.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      velocity.accumulator.push_back(
          std::make_unique<AccumulatorView<Vector3dAccumulator>>());
    data_[VELOCITY] = std::move(velocity);

    OutputData density;
    density.units        = "g cm^-3";
    density.title        = "Density";
    density.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      density.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());
    data_[DENSITY] = std::move(density);

    OutputData activity;
    activity.units        = "unitless";
    activity.title        = "Activity";
    activity.dataHandling = DataHandling::Average;
    unsigned int nTypes   = outputTypes_.size();
    // Only do activities if we have atoms in the selection
    if (nTypes > 0) {
      for (unsigned int i = 0; i < nBins_; i++)
        activity.accumulator.push_back(
            std::make_unique<AccumulatorView<StdVectorAccumulator>>());
      data_[ACTIVITY] = std::move(activity);
    }

    OutputData eField;
    eField.units        = "kcal/mol/angstroms/e";
    eField.title        = "Electric Field";
    eField.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      eField.accumulator.push_back(
          std::make_unique<AccumulatorView<Vector3dAccumulator>>());

    OutputData ePot;
    ePot.units        = "kcal/mol/e";
    ePot.title        = "Electrostatic Potential";
    ePot.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      ePot.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());

    OutputData charge;
    charge.units        = "e";
    charge.title        = "Charge";
    charge.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      charge.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());

    OutputData chargeVelocity;
    chargeVelocity.units        = "e/fs";
    chargeVelocity.title        = "Charge_Velocity";
    chargeVelocity.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      chargeVelocity.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());

    outputMask_.set(Z);
    outputMask_.set(TEMPERATURE);
    outputMask_.set(VELOCITY);
    outputMask_.set(DENSITY);
    outputMask_.set(ACTIVITY);

    int atomStorageLayout        = info_->getAtomStorageLayout();
    int rigidBodyStorageLayout   = info->getRigidBodyStorageLayout();
    int cutoffGroupStorageLayout = info->getCutoffGroupStorageLayout();

    if (atomStorageLayout & DataStorage::dslElectricField) {
      outputMask_.set(ELECTRICFIELD);
      outputMask_.set(ELECTROSTATICPOTENTIAL);

      data_[ELECTRICFIELD]          = std::move(eField);
      data_[ELECTROSTATICPOTENTIAL] = std::move(ePot);
    }

    if (info_->usesElectrostaticAtoms() ||
        atomStorageLayout & DataStorage::dslFlucQPosition) {
      outputMask_.set(CHARGE);

      data_[CHARGE] = std::move(charge);
    }

    if (atomStorageLayout & DataStorage::dslFlucQVelocity) {
      outputMask_.set(CHARGEVELOCITY);

      data_[CHARGEVELOCITY] = std::move(chargeVelocity);
    }
  }

  void RNEMDZ::processFrame(int istep) {
    SlabStatistics::processFrame(istep);

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    int binNo;
    int typeIndex(-1);
    RealType mass;
    Vector3d vel;
    RealType KE;
    RealType q = 0.0;
    RealType w = 0.0;
    Vector3d eField;

    std::vector<RealType> binMass(nBins_, 0.0);
    std::vector<Vector3d> binP(nBins_, V3Zero);
    std::vector<RealType> binCharge(nBins_, 0.0);
    std::vector<RealType> binChargeVelocity(nBins_, 0.0);
    std::vector<RealType> binKE(nBins_, 0.0);
    std::vector<Vector3d> binEField(nBins_, V3Zero);
    std::vector<int> binDOF(nBins_, 0);
    std::vector<int> binCount(nBins_, 0);
    std::vector<std::vector<int>> binTypeCounts;
    std::vector<int> binEFieldCount(nBins_, 0);

    if (outputMask_[ACTIVITY]) {
      binTypeCounts.resize(nBins_);
      for (unsigned int i = 0; i < nBins_; i++) {
        binTypeCounts[i].resize(outputTypes_.size(), 0);
      }
    }

    SimInfo::MoleculeIterator miter;
    std::vector<StuntDouble*>::iterator iiter;
    std::vector<AtomType*>::iterator at;
    Molecule* mol;
    StuntDouble* sd;
    AtomType* atype;

    for (mol = info_->beginMolecule(miter); mol != NULL;
         mol = info_->nextMolecule(miter)) {
      for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
           sd = mol->nextIntegrableObject(iiter)) {
        if (seleMan_.isSelected(sd)) {
          Vector3d pos = sd->getPos();
          binNo        = getBin(pos);

          mass = sd->getMass();
          vel  = sd->getVel();
          KE   = 0.5 * mass * vel.lengthSquare();

          if (outputMask_[ACTIVITY]) {
            typeIndex = -1;
            if (sd->isAtom()) {
              atype = static_cast<Atom*>(sd)->getAtomType();
              at = std::find(outputTypes_.begin(), outputTypes_.end(), atype);
              if (at != outputTypes_.end()) {
                typeIndex = std::distance(outputTypes_.begin(), at);
              }
            }
          }

          if (binNo >= 0 && binNo < int(nBins_)) {
            binCount[binNo]++;
            binMass[binNo] += mass;
            binP[binNo] += mass * vel;
            binKE[binNo] += KE;
            binDOF[binNo] += 3;

            if (outputMask_[ACTIVITY] && typeIndex != -1)
              binTypeCounts[binNo][typeIndex]++;

            if (outputMask_[CHARGE] || outputMask_[CHARGEVELOCITY]) {
              if (sd->isAtom()) {
                AtomType* atomType     = static_cast<Atom*>(sd)->getAtomType();
                FixedChargeAdapter fca = FixedChargeAdapter(atomType);
                if (fca.isFixedCharge()) { q = fca.getCharge(); }
                FluctuatingChargeAdapter fqa =
                    FluctuatingChargeAdapter(atomType);
                if (fqa.isFluctuatingCharge()) {
                  q += sd->getFlucQPos();
                  w += sd->getFlucQVel();
                }

                if (outputMask_[CHARGE]) binCharge[binNo] += q;
                if (outputMask_[CHARGEVELOCITY]) binChargeVelocity[binNo] += w;
              } else if (sd->isRigidBody()) {
                RigidBody* rb = static_cast<RigidBody*>(sd);
                std::vector<Atom*>::iterator ai;
                Atom* atom;
                for (atom = rb->beginAtom(ai); atom != NULL;
                     atom = rb->nextAtom(ai)) {
                  binNo                  = getBin(atom->getPos());
                  AtomType* atomType     = atom->getAtomType();
                  FixedChargeAdapter fca = FixedChargeAdapter(atomType);
                  if (fca.isFixedCharge()) { q = fca.getCharge(); }

                  FluctuatingChargeAdapter fqa =
                      FluctuatingChargeAdapter(atomType);
                  if (fqa.isFluctuatingCharge()) {
                    q += sd->getFlucQPos();
                    w += sd->getFlucQVel();
                  }

                  if (outputMask_[CHARGE]) binCharge[binNo] += q;
                  if (outputMask_[CHARGEVELOCITY])
                    binChargeVelocity[binNo] += w;
                }
              }
            }

            if (sd->isDirectional()) {
              Vector3d angMom = sd->getJ();
              Mat3x3d Ia      = sd->getI();
              if (sd->isLinear()) {
                int i = sd->linearAxis();
                int j = (i + 1) % 3;
                int k = (i + 2) % 3;
                binKE[binNo] += 0.5 * (angMom[j] * angMom[j] / Ia(j, j) +
                                       angMom[k] * angMom[k] / Ia(k, k));
                binDOF[binNo] += 2;
              } else {
                binKE[binNo] += 0.5 * (angMom[0] * angMom[0] / Ia(0, 0) +
                                       angMom[1] * angMom[1] / Ia(1, 1) +
                                       angMom[2] * angMom[2] / Ia(2, 2));
                binDOF[binNo] += 3;
              }
            }
          }
        }

        // Calculate the electric field (kcal/mol/e/Angstrom) for all atoms in
        // the box
        if (outputMask_[ELECTRICFIELD]) {
          if (sd->isRigidBody()) {
            RigidBody* rb = static_cast<RigidBody*>(sd);
            std::vector<Atom*>::iterator ai;
            Atom* atom;
            for (atom = rb->beginAtom(ai); atom != NULL;
                 atom = rb->nextAtom(ai)) {
              binNo  = getBin(atom->getPos());
              eField = atom->getElectricField();
              binEFieldCount[binNo]++;
              binEField[binNo] += eField;
            }
          } else {
            eField = sd->getElectricField();
            binNo  = getBin(sd->getPos());

            binEFieldCount[binNo]++;
            binEField[binNo] += eField;
          }
        }
      }
      if (seleMan_.isSelected(mol)) {
        Vector3d pos    = mol->getCom();
        binNo           = getBin(pos);
        int constraints = mol->getNConstraintPairs();
        binDOF[binNo] -= constraints;
      }
    }

    for (unsigned int i = 0; i < nBins_; i++) {
      RealType temp(0.0), ePot(0.0);
      Vector3d vel(0.0), eField(0.0);
      RealType z, den(0.0), binVolume(0.0), dz(0.0);
      std::vector<RealType> nden(outputTypes_.size(), 0.0);

      z = (((RealType)i + 0.5) / (RealType)nBins_) * hmat_(axis_, axis_);
      data_[Z].accumulator[i]->add(z);

      binVolume = volume_ / nBins_;
      dz        = hmat_(axis_, axis_) / (RealType)nBins_;

      // The calculations of the following properties are done regardless
      //   of whether or not the selected species are present in the bin
      if (outputMask_[ELECTRICFIELD] && binEFieldCount[i] > 0) {
        eField = binEField[i] / RealType(binEFieldCount[i]);
        data_[ELECTRICFIELD].accumulator[i]->add(eField);
      }

      if (outputMask_[ELECTROSTATICPOTENTIAL] && binEFieldCount[i] > 0) {
        ePot += eField[axis_] * dz;
        data_[ELECTROSTATICPOTENTIAL].accumulator[i]->add(ePot);
      }

      // For the following properties, zero should be added if the selected
      //   species is not present in the bin
      if (outputMask_[DENSITY]) {
        den = binMass[i] * Constants::densityConvert / binVolume;
        data_[DENSITY].accumulator[i]->add(den);
      }

      if (outputMask_[ACTIVITY]) {
        for (unsigned int j = 0; j < outputTypes_.size(); j++) {
          nden[j] = (binTypeCounts[i][j] / binVolume) *
                    Constants::concentrationConvert;
        }
        data_[ACTIVITY].accumulator[i]->add(nden);
      }

      if (binCount[i] > 0) {
        // The calculations of the following properties are undefined if
        //   the selected species is not found in the bin
        if (outputMask_[VELOCITY]) {
          vel = binP[i] / binMass[i];
          data_[VELOCITY].accumulator[i]->add(vel);
        }

        if (outputMask_[TEMPERATURE]) {
          temp = 2.0 * binKE[i] /
                 (binDOF[i] * Constants::kb * Constants::energyConvert);
          data_[TEMPERATURE].accumulator[i]->add(temp);
        }

        if (outputMask_[CHARGE])
          data_[CHARGE].accumulator[i]->add(binCharge[i]);

        if (outputMask_[CHARGEVELOCITY])
          data_[CHARGEVELOCITY].accumulator[i]->add(binChargeVelocity[i]);
      }
    }
  }

  RNEMDR::RNEMDR(SimInfo* info, const std::string& filename,
                 const std::string& sele, const std::string& comsele,
                 int nrbins) :
      ShellStatistics(info, filename, sele, comsele, nrbins) {
    setOutputName(getPrefix(filename) + ".rnemdR");

    // Pre-load the OutputData
    data_.resize(RNEMDR::ENDINDEX);

    OutputData r;
    r.units        = "Angstroms";
    r.title        = "R";
    r.dataHandling = DataHandling::Average;
    for (int i = 0; i < nBins_; i++)
      r.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());
    data_[R] = std::move(r);

    OutputData temperature;
    temperature.units        = "K";
    temperature.title        = "Temperature";
    temperature.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      temperature.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());
    data_[TEMPERATURE] = std::move(temperature);

    OutputData angularVelocity;
    angularVelocity.units        = "angstroms/fs";
    angularVelocity.title        = "Velocity";
    angularVelocity.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      angularVelocity.accumulator.push_back(
          std::make_unique<AccumulatorView<Vector3dAccumulator>>());
    data_[ANGULARVELOCITY] = std::move(angularVelocity);

    OutputData density;
    density.units        = "g cm^-3";
    density.title        = "Density";
    density.dataHandling = DataHandling::Average;
    for (unsigned int i = 0; i < nBins_; i++)
      density.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());
    data_[DENSITY] = std::move(density);

    outputMask_.set(R);
    outputMask_.set(TEMPERATURE);
    outputMask_.set(ANGULARVELOCITY);
    outputMask_.set(DENSITY);
  }

  void RNEMDR::processFrame(int istep) {
    ShellStatistics::processFrame(istep);

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    int binNo;
    RealType mass;
    Vector3d vel;
    Vector3d rPos;
    RealType KE;
    Vector3d L;
    Mat3x3d I;
    RealType r2;

    std::vector<int> binCount(nBins_, 0);
    std::vector<RealType> binMass(nBins_, 0.0);
    std::vector<Vector3d> binP(nBins_, V3Zero);
    std::vector<RealType> binOmega(nBins_, 0.0);
    std::vector<Vector3d> binL(nBins_, V3Zero);
    std::vector<Mat3x3d> binI(nBins_);
    std::vector<RealType> binKE(nBins_, 0.0);
    std::vector<int> binDOF(nBins_, 0);

    SimInfo::MoleculeIterator miter;
    std::vector<StuntDouble*>::iterator iiter;
    std::vector<AtomType*>::iterator at;
    Molecule* mol;
    StuntDouble* sd;

    // loop over the selected atoms:
    for (mol = info_->beginMolecule(miter); mol != NULL;
         mol = info_->nextMolecule(miter)) {
      for (sd = mol->beginIntegrableObject(iiter); sd != NULL;
           sd = mol->nextIntegrableObject(iiter)) {
        if (seleMan_.isSelected(sd)) {
          // figure out where that object is:
          binNo = getBin(sd->getPos());

          if (binNo >= 0 && binNo < int(nBins_)) {
            mass = sd->getMass();
            vel  = sd->getVel();
            rPos = sd->getPos() - coordinateOrigin_;
            KE   = 0.5 * mass * vel.lengthSquare();
            L    = mass * cross(rPos, vel);
            I    = outProduct(rPos, rPos) * mass;
            r2   = rPos.lengthSquare();
            I(0, 0) += mass * r2;
            I(1, 1) += mass * r2;
            I(2, 2) += mass * r2;

            binCount[binNo]++;
            binMass[binNo] += mass;
            binP[binNo] += mass * vel;
            binKE[binNo] += KE;
            binI[binNo] += I;
            binL[binNo] += L;
            binDOF[binNo] += 3;

            if (sd->isDirectional()) {
              Vector3d angMom = sd->getJ();
              Mat3x3d Ia      = sd->getI();
              if (sd->isLinear()) {
                int i = sd->linearAxis();
                int j = (i + 1) % 3;
                int k = (i + 2) % 3;
                binKE[binNo] += 0.5 * (angMom[j] * angMom[j] / Ia(j, j) +
                                       angMom[k] * angMom[k] / Ia(k, k));
                binDOF[binNo] += 2;
              } else {
                binKE[binNo] += 0.5 * (angMom[0] * angMom[0] / Ia(0, 0) +
                                       angMom[1] * angMom[1] / Ia(1, 1) +
                                       angMom[2] * angMom[2] / Ia(2, 2));
                binDOF[binNo] += 3;
              }
            }
          }
        }
      }
      if (seleMan_.isSelected(mol)) {
        Vector3d pos    = mol->getCom();
        binNo           = getBin(pos);
        int constraints = mol->getNConstraintPairs();
        binDOF[binNo] -= constraints;
      }
    }

    for (unsigned int i = 0; i < nBins_; i++) {
      RealType r, rinner, router, den(0.0), binVolume(0.0), temp(0.0);
      Vector3d omega(0.0);

      r      = (((RealType)i + 0.5) * binWidth_);
      rinner = (RealType)i * binWidth_;
      router = (RealType)(i + 1) * binWidth_;
      binVolume =
          (4.0 * Constants::PI * (pow(router, 3) - pow(rinner, 3))) / 3.0;

      data_[R].accumulator[i]->add(r);

      // For the following properties, zero should be added if the selected
      //   species is not present in the bin
      den = binMass[i] * Constants::densityConvert / binVolume;
      data_[DENSITY].accumulator[i]->add(den);

      if (binDOF[i] > 0) {
        // The calculations of the following properties are undefined if
        //   the selected species is not found in the bin
        omega = binI[i].inverse() * binL[i];
        data_[ANGULARVELOCITY].accumulator[i]->add(omega);

        temp = 2.0 * binKE[i] /
               (binDOF[i] * Constants::kb * Constants::energyConvert);
        data_[TEMPERATURE].accumulator[i]->add(temp);
      }
    }
  }

  RNEMDRTheta::RNEMDRTheta(SimInfo* info, const std::string& filename,
                           const std::string& sele, const std::string& comsele,
                           int nrbins, int nangleBins) :
      ShellStatistics(info, filename, sele, comsele, nrbins),
      nAngleBins_(nangleBins) {
    Globals* simParams                  = info->getSimParams();
    RNEMD::RNEMDParameters* rnemdParams = simParams->getRNEMDParameters();
    bool hasAngularMomentumFluxVector =
        rnemdParams->haveAngularMomentumFluxVector();

    if (hasAngularMomentumFluxVector) {
      std::vector<RealType> amf = rnemdParams->getAngularMomentumFluxVector();

      if (amf.size() != 3) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "RNEMDRTheta: Incorrect number of parameters specified for "
                 "angularMomentumFluxVector.\n"
                 "\tthere should be 3 parameters, but %lu were specified.\n",
                 amf.size());
        painCave.isFatal = 1;
        simError();
      }
      fluxVector_.x() = amf[0];
      fluxVector_.y() = amf[1];
      fluxVector_.z() = amf[2];
    } else {
      std::string fluxStr = rnemdParams->getFluxType();

      if (fluxStr.find("Lx") != std::string::npos) {
        fluxVector_ = V3X;
      } else if (fluxStr.find("Ly") != std::string::npos) {
        fluxVector_ = V3Y;
      } else {
        fluxVector_ = V3Z;
      }
    }

    fluxVector_.normalize();

    setOutputName(getPrefix(filename) + ".rnemdRTheta");

    // Pre-load the OutputData
    r_.units = "Angstroms";
    r_.title = "R";
    for (int i = 0; i < nBins_; i++)
      r_.accumulator.push_back(
          std::make_unique<AccumulatorView<RealAccumulator>>());

    angularVelocity_.units = "1/fs";
    angularVelocity_.title = "Projected Angular Velocity";
    for (unsigned int i = 0; i < nBins_; i++) {
      angularVelocity_.accumulator.push_back(
          std::make_unique<AccumulatorView<StdVectorAccumulator>>());
    }
  }

  std::pair<int, int> RNEMDRTheta::getBins(Vector3d pos) {
    std::pair<int, int> result;

    Vector3d rPos     = pos - coordinateOrigin_;
    RealType cosAngle = dot(rPos, fluxVector_) / rPos.length();

    result.first  = int(rPos.length() / binWidth_);
    result.second = int((nAngleBins_)*0.5 * (cosAngle + 1.0));
    return result;
  }

  void RNEMDRTheta::processFrame(int istep) {
    ShellStatistics::processFrame(istep);

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    StuntDouble* sd;
    int i;

    std::vector<std::vector<int>> binCount(nBins_);
    std::vector<std::vector<Mat3x3d>> binI(nBins_);
    std::vector<std::vector<Vector3d>> binL(nBins_);

    for (std::size_t i {}; i < nBins_; ++i) {
      binCount[i].resize(nAngleBins_);
      binI[i].resize(nAngleBins_);
      binL[i].resize(nAngleBins_);
    }

    // loop over the selected atoms:
    for (sd = seleMan_.beginSelected(i); sd != NULL;
         sd = seleMan_.nextSelected(i)) {
      // figure out where that object is:
      std::pair<int, int> bins = getBins(sd->getPos());

      if (bins.first >= 0 && bins.first < int(nBins_)) {
        if (bins.second >= 0 && bins.second < nAngleBins_) {
          Vector3d rPos = sd->getPos() - coordinateOrigin_;
          Vector3d vel  = sd->getVel();
          RealType m    = sd->getMass();
          Vector3d L    = m * cross(rPos, vel);
          Mat3x3d I(0.0);
          I           = outProduct(rPos, rPos) * m;
          RealType r2 = rPos.lengthSquare();
          I(0, 0) += m * r2;
          I(1, 1) += m * r2;
          I(2, 2) += m * r2;

          binI[bins.first][bins.second] += I;
          binL[bins.first][bins.second] += L;
          binCount[bins.first][bins.second]++;
        }
      }
    }

    for (unsigned int i = 0; i < nBins_; i++) {
      RealType r = (((RealType)i + 0.5) * binWidth_);
      r_.accumulator[i]->add(r);

      std::vector<RealType> projections(nAngleBins_);

      for (int j = 0; j < nAngleBins_; j++) {
        Vector3d omega(0.0);

        if (binCount[i][j] > 0) { omega = binI[i][j].inverse() * binL[i][j]; }

        // RealType omegaProj = dot(omega, fluxVector_);
        projections[j] = dot(omega, fluxVector_);
      }

      angularVelocity_.accumulator[i]->add(projections);
    }
  }

  void RNEMDRTheta::writeOutput() {
    std::ofstream outStream(outputFilename_.c_str());

    if (outStream.is_open()) {
      // write title
      outStream << "# SPATIAL STATISTICS\n";
      outStream << "#nBins = " << nBins_ << "\t binWidth = " << binWidth_
                << " maxR = " << nBins_ * binWidth_ << "\n";
      outStream << "#fluxVector = " << fluxVector_ << "\tBins = " << nAngleBins_
                << "\n";
      outStream << "#\t" << angularVelocity_.title << "("
                << angularVelocity_.units << ")\t\t";

      outStream << std::endl;

      outStream.precision(8);

      for (unsigned int i = 0; i < nBins_; i++) {
        std::size_t n {r_.accumulator[i]->getCount()};

        if (n != 0) {
          std::string message =
              "StaticAnalyser detected a numerical error writing: " +
              angularVelocity_.title + " for bin " + std::to_string(i);

          angularVelocity_.accumulator[i]->writeData(outStream, message);
        }

        outStream << std::endl;
      }
    }
  }
}  // namespace OpenMD
