/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "rnemd/RNEMD.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/ForceManager.hpp"
#include "brains/Thermo.hpp"
#include "io/Globals.hpp"
#include "math/ConvexHull.hpp"
#include "math/Polynomial.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector.hpp"
#include "math/Vector3.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "rnemd/RNEMDParameters.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/Accumulator.hpp"
#include "utils/Constants.hpp"
#include "utils/Tuple.hpp"

namespace OpenMD {
  namespace RNEMD {

    RNEMD::RNEMD(SimInfo* info, ForceManager* forceMan) :
        info_(info), forceMan_(forceMan), evaluator_(info_), seleMan_(info_),
        evaluatorA_(info_), evaluatorB_(info_), seleManA_(info_),
        seleManB_(info_), commonA_(info_), commonB_(info_),
        outputEvaluator_(info_), outputSeleMan_(info_) {
      trialCount_     = 0;
      failTrialCount_ = 0;
      failRootCount_  = 0;

      Globals* simParams           = info->getSimParams();
      RNEMDParameters* rnemdParams = simParams->getRNEMDParameters();

      usePeriodicBoundaryConditions_ =
          simParams->getUsePeriodicBoundaryConditions();

      doRNEMD_ = rnemdParams->getUseRNEMD();
      if (!doRNEMD_) return;

      // Determine Flux Type
      std::map<std::string, RNEMDFluxType> stringToFluxType;

      stringToFluxType["KE"]         = rnemdKE;
      stringToFluxType["Px"]         = rnemdPx;
      stringToFluxType["Py"]         = rnemdPy;
      stringToFluxType["Pz"]         = rnemdPz;
      stringToFluxType["Pvector"]    = rnemdPvector;
      stringToFluxType["Lx"]         = rnemdLx;
      stringToFluxType["Ly"]         = rnemdLy;
      stringToFluxType["Lz"]         = rnemdLz;
      stringToFluxType["Lvector"]    = rnemdLvector;
      stringToFluxType["Particle"]   = rnemdParticle;
      stringToFluxType["KE+Px"]      = rnemdKePx;
      stringToFluxType["KE+Py"]      = rnemdKePy;
      stringToFluxType["KE+Pvector"] = rnemdKePvector;
      stringToFluxType["KE+Lx"]      = rnemdKeLx;
      stringToFluxType["KE+Ly"]      = rnemdKeLy;
      stringToFluxType["KE+Lz"]      = rnemdKeLz;
      stringToFluxType["KE+Lvector"] = rnemdKeLvector;

      if (rnemdParams->haveFluxType()) {
        rnemdFluxTypeLabel_ = rnemdParams->getFluxType();
        rnemdFluxType_ = stringToFluxType.find(rnemdFluxTypeLabel_)->second;
      } else {
        std::string allowedFluxTypes;
        int currentLineLength = 0;

        for (std::map<std::string, RNEMDFluxType>::iterator fluxStrIter =
                 stringToFluxType.begin();
             fluxStrIter != stringToFluxType.end(); ++fluxStrIter) {
          allowedFluxTypes += fluxStrIter->first + ", ";
          currentLineLength += fluxStrIter->first.length() + 2;

          if (currentLineLength >= 50) {
            allowedFluxTypes += "\n\t\t";
            currentLineLength = 0;
          }
        }

        allowedFluxTypes.erase(allowedFluxTypes.length() - 2, 2);

        sprintf(
            painCave.errMsg,
            "RNEMD: No fluxType was set in the omd file. This parameter\n"
            "\tmust be set to use RNEMD, and can take any of these values:\n"
            "\t\t%s.\n",
            allowedFluxTypes.c_str());
        painCave.isFatal  = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      }

      // Determine Privileged Axis
      const std::string privAxis = rnemdParams->getPrivilegedAxis();

      if (privAxis == "x") {
        rnemdAxisLabel_      = "x";
        rnemdPrivilegedAxis_ = rnemdX;
      } else if (privAxis == "y") {
        rnemdAxisLabel_      = "y";
        rnemdPrivilegedAxis_ = rnemdY;
      } else {
        rnemdAxisLabel_      = "z";
        rnemdPrivilegedAxis_ = rnemdZ;
      }

      runTime_    = simParams->getRunTime();
      statusTime_ = simParams->getStatusTime();

      rnemdObjectSelection_ = rnemdParams->getObjectSelection();

      bool hasSlabWidth     = rnemdParams->haveSlabWidth();
      bool hasSlabACenter   = rnemdParams->haveSlabACenter();
      bool hasSlabBCenter   = rnemdParams->haveSlabBCenter();
      bool hasSphereARadius = rnemdParams->haveSphereARadius();
      bool hasSphereBRadius = rnemdParams->haveSphereBRadius();

      hasSelectionA_ = rnemdParams->haveSelectionA();
      hasSelectionB_ = rnemdParams->haveSelectionB();

      hasDividingArea_ = rnemdParams->haveDividingArea();
      dividingArea_    = rnemdParams->getDividingArea();

      bool hasCoordinateOrigin = rnemdParams->haveCoordinateOrigin();
      bool hasOutputFileName   = rnemdParams->haveOutputFileName();
      bool hasOutputFields     = rnemdParams->haveOutputFields();
      bool hasOutputSelection  = rnemdParams->haveOutputSelection();

      if (hasOutputSelection) {
        outputSelection_ = rnemdParams->getOutputSelection();
      } else {
        outputSelection_ = rnemdObjectSelection_;
      }

      if (hasCoordinateOrigin) {
        std::vector<RealType> co = rnemdParams->getCoordinateOrigin();
        if (co.size() != 3) {
          sprintf(painCave.errMsg,
                  "RNEMD: Incorrect number of parameters specified for "
                  "coordinateOrigin.\n"
                  "\tthere should be 3 parameters, but %lu were specified.\n",
                  co.size());
          painCave.isFatal = 1;
          simError();
        }
        coordinateOrigin_.x() = co[0];
        coordinateOrigin_.y() = co[1];
        coordinateOrigin_.z() = co[2];
      } else {
        coordinateOrigin_ = V3Zero;
      }

      outputEvaluator_.loadScriptString(outputSelection_);
      outputSeleMan_.setSelectionSet(outputEvaluator_.evaluate());
      std::set<AtomType*> osTypes = outputSeleMan_.getSelectedAtomTypes();
      std::copy(osTypes.begin(), osTypes.end(),
                std::back_inserter(outputTypes_));

      nBins_    = rnemdParams->getOutputBins();
      binWidth_ = rnemdParams->getOutputBinWidth();

      data_.resize(RNEMD::ENDINDEX);

      OutputData z;
      z.units    = "Angstroms";
      z.title    = rnemdAxisLabel_;
      z.dataType = "RealType";
      z.accumulator.reserve(nBins_);
      for (unsigned int i = 0; i < nBins_; i++)
        z.accumulator.push_back(new Accumulator());
      data_[Z]        = z;
      outputMap_["Z"] = Z;

      OutputData r;
      r.units    = "Angstroms";
      r.title    = "R";
      r.dataType = "RealType";
      r.accumulator.reserve(nBins_);
      for (unsigned int i = 0; i < nBins_; i++)
        r.accumulator.push_back(new Accumulator());
      data_[R]        = r;
      outputMap_["R"] = R;

      OutputData temperature;
      temperature.units    = "K";
      temperature.title    = "Temperature";
      temperature.dataType = "RealType";
      temperature.accumulator.reserve(nBins_);
      for (unsigned int i = 0; i < nBins_; i++)
        temperature.accumulator.push_back(new Accumulator());
      data_[TEMPERATURE]        = temperature;
      outputMap_["TEMPERATURE"] = TEMPERATURE;

      OutputData velocity;
      velocity.units    = "angstroms/fs";
      velocity.title    = "Velocity";
      velocity.dataType = "Vector3d";
      velocity.accumulator.reserve(nBins_);
      for (unsigned int i = 0; i < nBins_; i++)
        velocity.accumulator.push_back(new VectorAccumulator());
      data_[VELOCITY]        = velocity;
      outputMap_["VELOCITY"] = VELOCITY;

      OutputData angularVelocity;
      angularVelocity.units    = "angstroms^2/fs";
      angularVelocity.title    = "AngularVelocity";
      angularVelocity.dataType = "Vector3d";
      angularVelocity.accumulator.reserve(nBins_);
      for (unsigned int i = 0; i < nBins_; i++)
        angularVelocity.accumulator.push_back(new VectorAccumulator());
      data_[ANGULARVELOCITY]        = angularVelocity;
      outputMap_["ANGULARVELOCITY"] = ANGULARVELOCITY;

      OutputData density;
      density.units    = "g cm^-3";
      density.title    = "Density";
      density.dataType = "RealType";
      density.accumulator.reserve(nBins_);
      for (unsigned int i = 0; i < nBins_; i++)
        density.accumulator.push_back(new Accumulator());
      data_[DENSITY]        = density;
      outputMap_["DENSITY"] = DENSITY;

      OutputData activity;
      activity.units      = "unitless";
      activity.title      = "Activity";
      activity.dataType   = "Array2d";
      unsigned int nTypes = outputTypes_.size();
      activity.accumulatorArray2d.resize(nBins_);
      for (unsigned int i = 0; i < nBins_; i++) {
        activity.accumulatorArray2d[i].resize(nTypes);
        for (unsigned int j = 0; j < nTypes; j++) {
          activity.accumulatorArray2d[i][j] = new Accumulator();
        }
      }
      data_[ACTIVITY]        = (activity);
      outputMap_["ACTIVITY"] = ACTIVITY;

      OutputData eField;
      eField.units    = "kcal/mol/angstroms/e";
      eField.title    = "Electric Field";
      eField.dataType = "Vector3d";
      eField.accumulator.reserve(nBins_);
      for (unsigned int i = 0; i < nBins_; i++)
        eField.accumulator.push_back(new VectorAccumulator());
      data_[ELECTRICFIELD]        = eField;
      outputMap_["ELECTRICFIELD"] = ELECTRICFIELD;

      OutputData ePot;
      ePot.units    = "kcal/mol/e";
      ePot.title    = "Electrostatic Potential";
      ePot.dataType = "RealType";
      ePot.accumulator.reserve(nBins_);
      for (unsigned int i = 0; i < nBins_; i++)
        ePot.accumulator.push_back(new Accumulator());
      data_[ELECTROSTATICPOTENTIAL]        = ePot;
      outputMap_["ELECTROSTATICPOTENTIAL"] = ELECTROSTATICPOTENTIAL;

      if (hasOutputFields) {
        parseOutputFileFormat(rnemdParams->getOutputFields());
      } else {
        if (usePeriodicBoundaryConditions_)
          outputMask_.set(Z);
        else
          outputMask_.set(R);
        switch (rnemdFluxType_) {
        case rnemdKE:
        case rnemdRotKE:
        case rnemdFullKE:
          outputMask_.set(TEMPERATURE);
          break;
        case rnemdPx:
        case rnemdPy:
          outputMask_.set(VELOCITY);
          break;
        case rnemdPz:
        case rnemdPvector:
          outputMask_.set(VELOCITY);
          outputMask_.set(DENSITY);
          break;
        case rnemdLx:
        case rnemdLy:
        case rnemdLz:
        case rnemdLvector:
          outputMask_.set(ANGULARVELOCITY);
          break;
        case rnemdKeLx:
        case rnemdKeLy:
        case rnemdKeLz:
        case rnemdKeLvector:
          outputMask_.set(TEMPERATURE);
          outputMask_.set(ANGULARVELOCITY);
          break;
        case rnemdKePx:
        case rnemdKePy:
          outputMask_.set(TEMPERATURE);
          outputMask_.set(VELOCITY);
          break;
        case rnemdKePvector:
          outputMask_.set(TEMPERATURE);
          outputMask_.set(VELOCITY);
          outputMask_.set(DENSITY);
          break;
        default:
          break;
        }
      }

      if (hasOutputFileName) {
        rnemdFileName_ = rnemdParams->getOutputFileName();
      } else {
        rnemdFileName_ = getPrefix(info->getFinalConfigFileName()) + ".rnemd";
      }

      // Exchange time should not be less than time step and should be a
      // multiple of dt
      exchangeTime_  = rnemdParams->getExchangeTime();
      RealType dt    = simParams->getDt();
      RealType newET = std::ceil(exchangeTime_ / dt) * dt;

      if (std::fabs(newET - exchangeTime_) > 1e-6) {
        sprintf(painCave.errMsg,
                "RNEMD: The exchangeTime was reset to %lf,\n"
                "\t\twhich is a multiple of dt, %lf.\n",
                newET, dt);
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_WARNING;
        simError();
        exchangeTime_ = newET;
      }

      currentSnap_ = info->getSnapshotManager()->getCurrentSnapshot();
      hmat_        = currentSnap_->getHmat();

      // Set up the slab selection logic
      std::ostringstream selectionAstream;
      std::ostringstream selectionBstream;

      if (hasSelectionA_) {
        selectionA_ = rnemdParams->getSelectionA();
      } else {
        if (usePeriodicBoundaryConditions_) {
          slabWidth_ =
              hasSlabWidth ?
                  rnemdParams->getSlabWidth() :
                  hmat_(rnemdPrivilegedAxis_, rnemdPrivilegedAxis_) / 10.0;

          slabACenter_ = hasSlabACenter ? rnemdParams->getSlabACenter() : 0.0;

          selectionA_ = this->setSelection(slabACenter_);
        } else {
          if (hasSphereARadius)
            sphereARadius_ = rnemdParams->getSphereARadius();
          else {
            // use an initial guess to the size of the inner slab to be 1/10 the
            // radius of an approximately spherical hull:
            Thermo thermo(info);
            RealType hVol = thermo.getHullVolume();
            sphereARadius_ =
                0.1 * pow((3.0 * hVol / (4.0 * Constants::PI)), 1.0 / 3.0);
          }
          selectionAstream << "select r < " << sphereARadius_;
          selectionA_ = selectionAstream.str();
        }
      }

      if (hasSelectionB_) {
        selectionB_ = rnemdParams->getSelectionB();
      } else {
        if (usePeriodicBoundaryConditions_) {
          slabWidth_ =
              hasSlabWidth ?
                  rnemdParams->getSlabWidth() :
                  hmat_(rnemdPrivilegedAxis_, rnemdPrivilegedAxis_) / 10.0;

          slabBCenter_ =
              hasSlabBCenter ?
                  rnemdParams->getSlabBCenter() :
                  hmat_(rnemdPrivilegedAxis_, rnemdPrivilegedAxis_) / 2.0;

          selectionB_ = this->setSelection(slabBCenter_);
        } else {
          if (hasSphereBRadius) {
            sphereBRadius_ = rnemdParams->getSphereBRadius();
            selectionBstream << "select r > " << sphereBRadius_;
            selectionB_ = selectionBstream.str();
          } else {
            selectionB_    = "select hull";
            hasSelectionB_ = true;
          }
        }
      }

      // Static object evaluators
      evaluator_.loadScriptString(rnemdObjectSelection_);
      if (!evaluator_.isDynamic())
        seleMan_.setSelectionSet(evaluator_.evaluate());

      evaluatorA_.loadScriptString(selectionA_);
      if (!evaluatorA_.isDynamic())
        seleManA_.setSelectionSet(evaluatorA_.evaluate());

      evaluatorB_.loadScriptString(selectionB_);
      if (!evaluatorB_.isDynamic())
        seleManB_.setSelectionSet(evaluatorB_.evaluate());

      // Do some sanity checking
      int selectionCount = seleMan_.getSelectionCount();
      int nIntegrable    = info->getNGlobalIntegrableObjects();

      if (selectionCount > nIntegrable) {
        sprintf(painCave.errMsg,
                "RNEMD: The current objectSelection,\n"
                "\t\t%s\n"
                "\thas resulted in %d selected objects.  However,\n"
                "\tthe total number of integrable objects in the system\n"
                "\tis only %d.  This is almost certainly not what you want\n"
                "\tto do.  A likely cause of this is forgetting the _RB_0\n"
                "\tselector in the selection script!\n",
                rnemdObjectSelection_.c_str(), selectionCount, nIntegrable);
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_WARNING;
        simError();
      }
    }

    RNEMD::~RNEMD() {
      if (!doRNEMD_) return;
#ifdef IS_MPI
      if (worldRank == 0) {
#endif
        writeOutputFile();

        rnemdFile_.close();
#ifdef IS_MPI
      }
#endif

      for (auto& data : data_) {
        if (!data.accumulatorArray2d.empty())
          Utils::deletePointers(data.accumulatorArray2d);
        else
          Utils::deletePointers(data.accumulator);
      }
    }

    void RNEMD::getStarted() {
      if (!doRNEMD_) return;
      collectData();
      writeOutputFile();
    }

    void RNEMD::doRNEMD() {
      if (!doRNEMD_) return;
      trialCount_++;
      hmat_ = currentSnap_->getHmat();

      // dynamic object evaluators:
      evaluator_.loadScriptString(rnemdObjectSelection_);
      if (evaluator_.isDynamic())
        seleMan_.setSelectionSet(evaluator_.evaluate());

      evaluatorA_.loadScriptString(selectionA_);
      if (evaluatorA_.isDynamic())
        seleManA_.setSelectionSet(evaluatorA_.evaluate());

      evaluatorB_.loadScriptString(selectionB_);
      if (evaluatorB_.isDynamic())
        seleManB_.setSelectionSet(evaluatorB_.evaluate());

      commonA_ = seleManA_ & seleMan_;
      commonB_ = seleManB_ & seleMan_;

      // Target exchange quantities (in each exchange) = flux * dividingArea *
      // dt flux = target flux dividingArea = smallest dividing surface between
      // the two regions dt = exchange time interval

      RealType area = getDefaultDividingArea();

      kineticTarget_  = kineticFlux_ * exchangeTime_ * area;
      momentumTarget_ = momentumFluxVector_ * exchangeTime_ * area;
      angularMomentumTarget_ =
          angularMomentumFluxVector_ * exchangeTime_ * area;
      particleTarget_ = particleFlux_ * exchangeTime_ * area;

      if (std::fabs(particleTarget_) > 1.0) {
        sprintf(painCave.errMsg,
                "RNEMD: The current particleFlux,\n"
                "\t\t%f\n"
                "\thas resulted in a target particle exchange of %f.\n"
                "\tThis is equivalent to moving more than one particle\n"
                "\tduring each exchange.  Please reduce your particleFlux.\n",
                particleFlux_, particleTarget_);
        painCave.isFatal  = 1;
        painCave.severity = OPENMD_ERROR;
        simError();
      }

      this->doRNEMDImpl(commonA_, commonB_);
    }

    void RNEMD::collectData() {
      if (!doRNEMD_) return;
      currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
      hmat_        = currentSnap_->getHmat();

      // collectData() can be called more frequently than the doRNEMD(), so use
      // the computed area from the last exchange time:
      RealType area = getDefaultDividingArea();
      areaAccumulator_.add(area);

      Vector3d u = angularMomentumFluxVector_;
      u.normalize();

      if (outputEvaluator_.isDynamic()) {
        outputSeleMan_.setSelectionSet(outputEvaluator_.evaluate());
      }

      int binNo;
      int typeIndex(-1);
      RealType mass;
      Vector3d vel;
      Vector3d rPos;
      RealType KE;
      Vector3d L;
      Mat3x3d I;
      RealType r2;
      Vector3d eField;

      vector<RealType> binMass(nBins_, 0.0);
      vector<Vector3d> binP(nBins_, V3Zero);
      vector<RealType> binOmega(nBins_, 0.0);
      vector<Vector3d> binL(nBins_, V3Zero);
      vector<Mat3x3d> binI(nBins_);
      vector<RealType> binKE(nBins_, 0.0);
      vector<Vector3d> binEField(nBins_, V3Zero);
      vector<int> binDOF(nBins_, 0);
      vector<int> binCount(nBins_, 0);
      vector<int> binEFieldCount(nBins_, 0);
      vector<vector<int>> binTypeCounts;

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
          if (outputSeleMan_.isSelected(sd)) {
            Vector3d pos = sd->getPos();
            binNo        = getBin(pos);

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
              binI[binNo] += I;
              binL[binNo] += L;
              binDOF[binNo] += 3;

              if (outputMask_[ACTIVITY] && typeIndex != -1)
                binTypeCounts[binNo][typeIndex]++;

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
      }

#ifdef IS_MPI
      for (unsigned int i = 0; i < nBins_; i++) {
        MPI_Allreduce(MPI_IN_PLACE, &binCount[i], 1, MPI_INT, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &binMass[i], 1, MPI_REALTYPE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, binP[i].getArrayPointer(), 3, MPI_REALTYPE,
                      MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, binL[i].getArrayPointer(), 3, MPI_REALTYPE,
                      MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, binI[i].getArrayPointer(), 9, MPI_REALTYPE,
                      MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &binKE[i], 1, MPI_REALTYPE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &binDOF[i], 1, MPI_INT, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &binEFieldCount[i], 1, MPI_INT, MPI_SUM,
                      MPI_COMM_WORLD);

        if (outputMask_[ELECTRICFIELD]) {
          MPI_Allreduce(MPI_IN_PLACE, binEField[i].getArrayPointer(), 3,
                        MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
        }
        if (outputMask_[ACTIVITY]) {
          MPI_Allreduce(MPI_IN_PLACE, &binTypeCounts[i][0], outputTypes_.size(),
                        MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        }
      }
#endif

      Vector3d omega;
      RealType z, r, temp, binVolume, den(0.0), dz(0.0);
      std::vector<RealType> nden(outputTypes_.size(), 0.0);
      RealType boxVolume = currentSnap_->getVolume();
      RealType ePot(0.0);

      for (unsigned int i = 0; i < nBins_; i++) {
        if (usePeriodicBoundaryConditions_) {
          z = (((RealType)i + 0.5) / (RealType)nBins_) *
              hmat_(rnemdPrivilegedAxis_, rnemdPrivilegedAxis_);
          dynamic_cast<Accumulator*>(data_[Z].accumulator[i])->add(z);

          binVolume = boxVolume / nBins_;
          dz        = hmat_(rnemdPrivilegedAxis_, rnemdPrivilegedAxis_) /
               (RealType)nBins_;
        } else {
          r = (((RealType)i + 0.5) * binWidth_);
          dynamic_cast<Accumulator*>(data_[R].accumulator[i])->add(r);

          RealType rinner = (RealType)i * binWidth_;
          RealType router = (RealType)(i + 1) * binWidth_;
          binVolume =
              (4.0 * Constants::PI * (pow(router, 3) - pow(rinner, 3))) / 3.0;
        }

        // The calculations of the following properties are done regardless
        //   of whether or not the selected species are present in the bin
        if (outputMask_[ELECTRICFIELD] && binEFieldCount[i] > 0) {
          eField = binEField[i] / RealType(binEFieldCount[i]);
          dynamic_cast<VectorAccumulator*>(data_[ELECTRICFIELD].accumulator[i])
              ->add(eField);
        }

        if (outputMask_[ELECTROSTATICPOTENTIAL]) {
          if (usePeriodicBoundaryConditions_ && binEFieldCount[i] > 0) {
            ePot += eField[rnemdPrivilegedAxis_] * dz;
            dynamic_cast<Accumulator*>(
                data_[ELECTROSTATICPOTENTIAL].accumulator[i])
                ->add(ePot);
          }
        }

        // For the following properties, zero should be added if the selected
        //   species is not present in the bin
        if (outputMask_[DENSITY]) {
          den = binMass[i] * Constants::densityConvert / binVolume;
          dynamic_cast<Accumulator*>(data_[DENSITY].accumulator[i])->add(den);
        }

        if (outputMask_[ACTIVITY]) {
          for (unsigned int j = 0; j < outputTypes_.size(); j++) {
            nden[j] = (binTypeCounts[i][j] / binVolume) *
                      Constants::concentrationConvert;
            dynamic_cast<Accumulator*>(data_[ACTIVITY].accumulatorArray2d[i][j])
                ->add(nden[j]);
          }
        }

        if (binCount[i] > 0) {
          // The calculations of the following properties are meaningless if
          //   the selected species is not found in the bin
          if (outputMask_[VELOCITY]) {
            vel = binP[i] / binMass[i];
            dynamic_cast<VectorAccumulator*>(data_[VELOCITY].accumulator[i])
                ->add(vel);
          }

          if (outputMask_[ANGULARVELOCITY]) {
            omega = binI[i].inverse() * binL[i];
            dynamic_cast<VectorAccumulator*>(
                data_[ANGULARVELOCITY].accumulator[i])
                ->add(omega);
          }

          if (outputMask_[TEMPERATURE]) {
            temp = 2.0 * binKE[i] /
                   (binDOF[i] * Constants::kb * Constants::energyConvert);
            dynamic_cast<Accumulator*>(data_[TEMPERATURE].accumulator[i])
                ->add(temp);
          }
        }
      }

      hasData_ = true;
    }

    void RNEMD::writeOutputFile() {
      if (!doRNEMD_) return;
      if (!hasData_) return;

#ifdef IS_MPI
      // If we're the primary node, should we print out the results
      int worldRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

      if (worldRank == 0) {
#endif
        rnemdFile_.open(rnemdFileName_.c_str(),
                        std::ios::out | std::ios::trunc);

        if (!rnemdFile_) {
          sprintf(painCave.errMsg, "Could not open \"%s\" for RNEMD output.\n",
                  rnemdFileName_.c_str());
          painCave.isFatal = 1;
          simError();
        }

        RealType time    = currentSnap_->getTime();
        RealType avgArea = areaAccumulator_.getAverage();

        RealType Jz(0.0);
        Vector3d JzP(V3Zero);
        Vector3d JzL(V3Zero);
        RealType Jpart(0.0);

        if (time >= info_->getSimParams()->getDt()) {
          Jz  = kineticExchange_ / (time * avgArea) / Constants::energyConvert;
          JzP = momentumExchange_ / (time * avgArea);
          JzL = angularMomentumExchange_ / (time * avgArea);
          Jpart = particleExchange_ / (time * avgArea);
        }

        rnemdFile_
            << "#######################################################\n";
        rnemdFile_ << "# RNEMD {\n";
        rnemdFile_ << "#    exchangeMethod  = \"" << rnemdMethodLabel_
                   << "\";\n";
        rnemdFile_ << "#    fluxType  = \"" << rnemdFluxTypeLabel_ << "\";\n";

        if (usePeriodicBoundaryConditions_)
          rnemdFile_ << "#    privilegedAxis = " << rnemdAxisLabel_ << ";\n";

        rnemdFile_ << "#    exchangeTime = " << exchangeTime_ << ";\n";
        rnemdFile_ << "#    objectSelection = \"" << rnemdObjectSelection_
                   << "\";\n";
        rnemdFile_ << "#    selectionA = \"" << selectionA_ << "\";\n";
        rnemdFile_ << "#    selectionB = \"" << selectionB_ << "\";\n";
        rnemdFile_ << "#    outputSelection = \"" << outputSelection_
                   << "\";\n";
        rnemdFile_ << "# }\n";
        rnemdFile_
            << "#######################################################\n";
        rnemdFile_ << "# RNEMD report:\n";
        rnemdFile_ << "#      running time = " << time << " fs\n";
        rnemdFile_ << "# Target flux:\n";
        rnemdFile_ << "#           kinetic = "
                   << kineticFlux_ / Constants::energyConvert
                   << " (kcal/mol/A^2/fs)\n";
        rnemdFile_ << "#          momentum = " << momentumFluxVector_
                   << " (amu/A/fs^2)\n";
        rnemdFile_ << "#  angular momentum = " << angularMomentumFluxVector_
                   << " (amu/A^2/fs^2)\n";
        rnemdFile_ << "#          particle = " << particleFlux_
                   << " (particles/A^2/fs)\n";

        rnemdFile_ << "# Target one-time exchanges:\n";
        rnemdFile_ << "#          kinetic = "
                   << kineticTarget_ / Constants::energyConvert
                   << " (kcal/mol)\n";
        rnemdFile_ << "#          momentum = " << momentumTarget_
                   << " (amu*A/fs)\n";
        rnemdFile_ << "#  angular momentum = " << angularMomentumTarget_
                   << " (amu*A^2/fs)\n";
        rnemdFile_ << "#          particle = " << particleTarget_
                   << " (particles)\n";
        rnemdFile_ << "# Actual exchange totals:\n";
        rnemdFile_ << "#          kinetic = "
                   << kineticExchange_ / Constants::energyConvert
                   << " (kcal/mol)\n";
        rnemdFile_ << "#          momentum = " << momentumExchange_
                   << " (amu*A/fs)\n";
        rnemdFile_ << "#  angular momentum = " << angularMomentumExchange_
                   << " (amu*A^2/fs)\n";
        rnemdFile_ << "#         particles = " << particleExchange_
                   << " (particles)\n";

        rnemdFile_ << "# Actual flux:\n";
        rnemdFile_ << "#          kinetic = " << Jz << " (kcal/mol/A^2/fs)\n";
        rnemdFile_ << "#          momentum = " << JzP << " (amu/A/fs^2)\n";
        rnemdFile_ << "#  angular momentum = " << JzL << " (amu/A^2/fs^2)\n";
        rnemdFile_ << "#          particle = " << Jpart
                   << " (particles/A^2/fs)\n";
        rnemdFile_ << "# Exchange statistics:\n";
        rnemdFile_ << "#               attempted = " << trialCount_ << "\n";
        rnemdFile_ << "#                  failed = " << failTrialCount_ << "\n";
        if (rnemdMethodLabel_ == "NIVS") {
          rnemdFile_ << "#  NIVS root-check errors = " << failRootCount_
                     << "\n";
        }
        rnemdFile_
            << "#######################################################\n";

        // write title
        rnemdFile_ << "#";
        for (unsigned int i = 0; i < outputMask_.size(); ++i) {
          if (outputMask_[i]) {
            rnemdFile_ << "\t" << data_[i].title << "(" << data_[i].units
                       << ")";
            // add some extra tabs for column alignment
            if (data_[i].dataType == "Vector3d") rnemdFile_ << "\t\t";
            if (data_[i].dataType == "Array2d") {
              rnemdFile_ << "(";
              for (unsigned int j = 0;
                   j < data_[i].accumulatorArray2d[0].size(); j++) {
                rnemdFile_ << outputTypes_[j]->getName() << "\t";
              }
              rnemdFile_ << ")\t";
            }
          }
        }
        rnemdFile_ << std::endl;

        rnemdFile_.precision(8);

        for (unsigned int j = 0; j < nBins_; j++) {
          for (unsigned int i = 0; i < outputMask_.size(); ++i) {
            if (outputMask_[i]) {
              if (data_[i].dataType == "RealType")
                writeReal(i, j);
              else if (data_[i].dataType == "Vector3d")
                writeVector(i, j);
              else if (data_[i].dataType == "Array2d")
                writeArray(i, j);
              else {
                sprintf(painCave.errMsg,
                        "RNEMD found an unknown data type for: %s ",
                        data_[i].title.c_str());
                painCave.isFatal = 1;
                simError();
              }
            }
          }
          rnemdFile_ << std::endl;
        }

        rnemdFile_
            << "#######################################################\n";
        rnemdFile_
            << "# 95% confidence intervals in those quantities follow:\n";
        rnemdFile_
            << "#######################################################\n";

        for (unsigned int j = 0; j < nBins_; j++) {
          rnemdFile_ << "#";
          for (unsigned int i = 0; i < outputMask_.size(); ++i) {
            if (outputMask_[i]) {
              if (data_[i].dataType == "RealType")
                writeRealErrorBars(i, j);
              else if (data_[i].dataType == "Vector3d")
                writeVectorErrorBars(i, j);
              else if (data_[i].dataType == "Array2d")
                writeArrayErrorBars(i, j);
              else {
                sprintf(painCave.errMsg,
                        "RNEMD found an unknown data type for: %s ",
                        data_[i].title.c_str());
                painCave.isFatal = 1;
                simError();
              }
            }
          }
          rnemdFile_ << std::endl;
        }

        rnemdFile_.flush();
        rnemdFile_.close();

#ifdef IS_MPI
      }
#endif
    }

    void RNEMD::setKineticFlux(RealType kineticFlux) {
      // convert the kcal / mol / Angstroms^2 / fs values in the md file
      // into  amu / fs^3:
      kineticFlux_ = kineticFlux * Constants::energyConvert;
    }

    void RNEMD::setParticleFlux(RealType particleFlux) {
      particleFlux_ = particleFlux;
    }

    void RNEMD::setMomentumFluxVector(
        const std::vector<RealType>& momentumFluxVector) {
      if (momentumFluxVector.size() != 3) {
        sprintf(painCave.errMsg,
                "RNEMD: Incorrect number of parameters specified for "
                "momentumFluxVector.\n"
                "\tthere should be 3 parameters, but %lu were specified.\n",
                momentumFluxVector.size());
        painCave.isFatal = 1;
        simError();
      }

      momentumFluxVector_.x() = momentumFluxVector[0];
      momentumFluxVector_.y() = momentumFluxVector[1];
      momentumFluxVector_.z() = momentumFluxVector[2];
    }

    void RNEMD::setAngularMomentumFluxVector(
        const std::vector<RealType>& angularMomentumFluxVector) {
      if (angularMomentumFluxVector.size() != 3) {
        sprintf(painCave.errMsg,
                "RNEMD: Incorrect number of parameters specified for "
                "angularMomentumFluxVector.\n"
                "\tthere should be 3 parameters, but %lu were specified.\n",
                angularMomentumFluxVector.size());
        painCave.isFatal = 1;
        simError();
      }

      angularMomentumFluxVector_.x() = angularMomentumFluxVector[0];
      angularMomentumFluxVector_.y() = angularMomentumFluxVector[1];
      angularMomentumFluxVector_.z() = angularMomentumFluxVector[2];
    }

    void RNEMD::parseOutputFileFormat(const std::string& format) {
      if (!doRNEMD_) return;
      StringTokenizer tokenizer(format, " ,;|\t\n\r");

      while (tokenizer.hasMoreTokens()) {
        std::string token(tokenizer.nextToken());
        toUpper(token);
        OutputMapType::iterator i = outputMap_.find(token);
        if (i != outputMap_.end()) {
          outputMask_.set(i->second);
        } else {
          sprintf(painCave.errMsg,
                  "RNEMD::parseOutputFileFormat: %s is not a recognized\n"
                  "\toutputFileFormat keyword.\n",
                  token.c_str());
          painCave.isFatal  = 0;
          painCave.severity = OPENMD_ERROR;
          simError();
        }
      }
    }

    std::string RNEMD::setSelection(RealType& slabCenter) {
      bool printSlabCenterWarning {false};

      Vector3d tempSlabCenter {V3Zero};
      tempSlabCenter[rnemdPrivilegedAxis_] = slabCenter;

      RealType hmat_2 = hmat_(rnemdPrivilegedAxis_, rnemdPrivilegedAxis_) / 2.0;

      if (slabCenter > hmat_2) {
        currentSnap_->wrapVector(tempSlabCenter);
        printSlabCenterWarning = true;
      } else if (slabCenter < -hmat_2) {
        currentSnap_->wrapVector(tempSlabCenter);
        printSlabCenterWarning = true;
      }

      if (printSlabCenterWarning) {
        sprintf(painCave.errMsg,
                "The given slab center was set to %0.2f. In the wrapped "
                "coordinates\n"
                "\t[-Hmat/2, +Hmat/2], this has been remapped to %0.2f.\n",
                slabCenter, tempSlabCenter[rnemdPrivilegedAxis_]);
        painCave.isFatal  = 0;
        painCave.severity = OPENMD_WARNING;
        simError();

        slabCenter = tempSlabCenter[rnemdPrivilegedAxis_];
      }

      Vector3d leftSlab {V3Zero};
      const RealType& leftSlabBoundary = leftSlab[rnemdPrivilegedAxis_];
      leftSlab[rnemdPrivilegedAxis_]   = slabCenter - 0.5 * slabWidth_;
      currentSnap_->wrapVector(leftSlab);

      Vector3d rightSlab {V3Zero};
      const RealType& rightSlabBoundary = rightSlab[rnemdPrivilegedAxis_];
      rightSlab[rnemdPrivilegedAxis_]   = slabCenter + 0.5 * slabWidth_;
      currentSnap_->wrapVector(rightSlab);

      std::ostringstream selectionStream;

      selectionStream << "select wrapped" << rnemdAxisLabel_
                      << " >= " << leftSlabBoundary;

      if (leftSlabBoundary > rightSlabBoundary)
        selectionStream << " || wrapped" << rnemdAxisLabel_ << " < "
                        << rightSlabBoundary;
      else
        selectionStream << " && wrapped" << rnemdAxisLabel_ << " < "
                        << rightSlabBoundary;

      return selectionStream.str();
    }

    RealType RNEMD::getDefaultDividingArea() {
      if (hasDividingArea_) return dividingArea_;

      Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

      if (hasSelectionA_) {
        if (evaluatorA_.hasSurfaceArea()) {
          areaA_   = evaluatorA_.getSurfaceArea();
          volumeA_ = evaluatorA_.getVolume();
        } else {
          int isd;
          StuntDouble* sd;
          std::vector<StuntDouble*> aSites;
          seleManA_.setSelectionSet(evaluatorA_.evaluate());
          for (sd = seleManA_.beginSelected(isd); sd != NULL;
               sd = seleManA_.nextSelected(isd)) {
            aSites.push_back(sd);
          }
#if defined(HAVE_QHULL)
          ConvexHull* surfaceMeshA = new ConvexHull();
          surfaceMeshA->computeHull(aSites);
          areaA_   = surfaceMeshA->getArea();
          volumeA_ = surfaceMeshA->getVolume();
          delete surfaceMeshA;
#else
          sprintf(
              painCave.errMsg,
              "RNEMD::getDividingArea : Hull calculation is not possible\n"
              "\twithout libqhull. Please rebuild OpenMD with qhull enabled.");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
#endif
        }
      } else {
        if (usePeriodicBoundaryConditions_) {
          // in periodic boundaries, the surface area is twice the
          // area of the current box, normal to the privileged axis:
          switch (rnemdPrivilegedAxis_) {
          case rnemdX:
            areaA_ = 2.0 * snap->getYZarea();
            break;
          case rnemdY:
            areaA_ = 2.0 * snap->getXZarea();
            break;
          case rnemdZ:
          default:
            areaA_ = 2.0 * snap->getXYarea();
          }

          volumeA_ = areaA_ * slabWidth_;
        } else {
          // in non-periodic simulations, without explicitly setting
          // selections, the sphere radius sets the surface area of the
          // dividing surface:
          areaA_   = 4.0 * Constants::PI * std::pow(sphereARadius_, 2);
          volumeA_ = 4.0 * Constants::PI * std::pow(sphereARadius_, 3) / 3.0;
        }
      }

      if (hasSelectionB_) {
        if (evaluatorB_.hasSurfaceArea()) {
          areaB_   = evaluatorB_.getSurfaceArea();
          volumeB_ = evaluatorB_.getVolume();
        } else {
          int isd;
          StuntDouble* sd;
          std::vector<StuntDouble*> bSites;
          seleManB_.setSelectionSet(evaluatorB_.evaluate());
          for (sd = seleManB_.beginSelected(isd); sd != NULL;
               sd = seleManB_.nextSelected(isd)) {
            bSites.push_back(sd);
          }

#if defined(HAVE_QHULL)
          ConvexHull* surfaceMeshB = new ConvexHull();
          surfaceMeshB->computeHull(bSites);
          areaB_   = surfaceMeshB->getArea();
          volumeB_ = surfaceMeshB->getVolume();
          delete surfaceMeshB;
#else
          sprintf(
              painCave.errMsg,
              "RNEMD::getDividingArea : Hull calculation is not possible\n"
              "\twithout libqhull. Please rebuild OpenMD with qhull enabled.");
          painCave.severity = OPENMD_ERROR;
          painCave.isFatal  = 1;
          simError();
#endif
        }
      } else {
        if (usePeriodicBoundaryConditions_) {
          // in periodic boundaries, the surface area is twice the
          // area of the current box, normal to the privileged axis:
          switch (rnemdPrivilegedAxis_) {
          case rnemdX:
            areaB_ = 2.0 * snap->getYZarea();
            break;
          case rnemdY:
            areaB_ = 2.0 * snap->getXZarea();
            break;
          case rnemdZ:
          default:
            areaB_ = 2.0 * snap->getXYarea();
          }

          volumeB_ = areaB_ * slabWidth_;
        } else {
          // in non-periodic simulations, without explicitly setting
          // selections, but if a sphereBradius has been set, just use that:
          areaB_ = 4.0 * Constants::PI * pow(sphereBRadius_, 2);
          Thermo thermo(info_);
          RealType hVol = thermo.getHullVolume();
          volumeB_ = hVol - 4.0 * Constants::PI * pow(sphereBRadius_, 3) / 3.0;
        }
      }

      dividingArea_    = min(areaA_, areaB_);
      hasDividingArea_ = true;
      return dividingArea_;
    }

    int RNEMD::getBin(Vector3d pos) {
      if (usePeriodicBoundaryConditions_) {
        currentSnap_->wrapVector(pos);

        return int(nBins_ *
                   (pos[rnemdPrivilegedAxis_] /
                        hmat_(rnemdPrivilegedAxis_, rnemdPrivilegedAxis_) +
                    0.5)) %
               nBins_;
      } else {
        Vector3d rPos = pos - coordinateOrigin_;
        return int(rPos.length() / binWidth_);
      }
    }

    void RNEMD::writeReal(int index, unsigned int bin) {
      if (!doRNEMD_) return;
      assert(index >= 0 && index < ENDINDEX);
      assert(bin < nBins_);

      RealType s;
      std::size_t count = data_[index].accumulator[bin]->count();

      if (count == 0) {
        rnemdFile_ << "\t";
      } else {
        dynamic_cast<Accumulator*>(data_[index].accumulator[bin])
            ->getAverage(s);

        if (std::isinf(s) || std::isnan(s)) {
          sprintf(painCave.errMsg,
                  "RNEMD detected a numerical error writing: %s for bin %u",
                  data_[index].title.c_str(), bin);
          painCave.isFatal = 1;
          simError();
        } else {
          rnemdFile_ << "\t" << s;
        }
      }
    }

    void RNEMD::writeVector(int index, unsigned int bin) {
      if (!doRNEMD_) return;
      assert(index >= 0 && index < ENDINDEX);
      assert(bin < nBins_);

      Vector3d s;
      std::size_t count = data_[index].accumulator[bin]->count();

      if (count == 0) {
        rnemdFile_ << "\t\t\t";
      } else {
        dynamic_cast<VectorAccumulator*>(data_[index].accumulator[bin])
            ->getAverage(s);

        if (std::isinf(s[0]) || std::isnan(s[0]) || std::isinf(s[1]) ||
            std::isnan(s[1]) || std::isinf(s[2]) || std::isnan(s[2])) {
          sprintf(painCave.errMsg,
                  "RNEMD detected a numerical error writing: %s for bin %u",
                  data_[index].title.c_str(), bin);
          painCave.isFatal = 1;
          simError();
        } else {
          rnemdFile_ << "\t" << s[0] << "\t" << s[1] << "\t" << s[2];
        }
      }
    }

    void RNEMD::writeArray(int index, unsigned int bin) {
      if (!doRNEMD_) return;
      assert(index >= 0 && index < ENDINDEX);
      assert(bin < nBins_);

      RealType s;
      std::size_t columns = data_[index].accumulatorArray2d[0].size();

      for (std::size_t i = 0; i < columns; i++) {
        std::size_t count = data_[index].accumulatorArray2d[bin][i]->count();

        if (count == 0) {
          rnemdFile_ << "\t";
        } else {
          dynamic_cast<Accumulator*>(data_[index].accumulatorArray2d[bin][i])
              ->getAverage(s);

          if (std::isinf(s) || std::isnan(s)) {
            sprintf(painCave.errMsg,
                    "RNEMD detected a numerical error writing: %s for bin %u, "
                    "column %u",
                    data_[index].title.c_str(), bin,
                    static_cast<unsigned int>(i));
            painCave.isFatal = 1;
            simError();
          } else {
            rnemdFile_ << "\t" << s;
          }
        }
      }
    }

    void RNEMD::writeRealErrorBars(int index, unsigned int bin) {
      if (!doRNEMD_) return;
      assert(index >= 0 && index < ENDINDEX);
      assert(bin < nBins_);

      RealType s;
      std::size_t count = data_[index].accumulator[bin]->count();

      if (count == 0) {
        rnemdFile_ << "\t";
      } else {
        dynamic_cast<Accumulator*>(data_[index].accumulator[bin])
            ->get95percentConfidenceInterval(s);

        if (std::isinf(s) || std::isnan(s)) {
          sprintf(painCave.errMsg,
                  "RNEMD detected a numerical error writing: %s std. dev. for "
                  "bin %u",
                  data_[index].title.c_str(), bin);
          painCave.isFatal = 1;
          simError();
        } else {
          rnemdFile_ << "\t" << s;
        }
      }
    }

    void RNEMD::writeVectorErrorBars(int index, unsigned int bin) {
      if (!doRNEMD_) return;
      assert(index >= 0 && index < ENDINDEX);
      assert(bin < nBins_);

      Vector3d s;
      std::size_t count = data_[index].accumulator[bin]->count();

      if (count == 0) {
        rnemdFile_ << "\t\t\t";
      } else {
        dynamic_cast<VectorAccumulator*>(data_[index].accumulator[bin])
            ->get95percentConfidenceInterval(s);

        if (std::isinf(s[0]) || std::isnan(s[0]) || std::isinf(s[1]) ||
            std::isnan(s[1]) || std::isinf(s[2]) || std::isnan(s[2])) {
          sprintf(painCave.errMsg,
                  "RNEMD detected a numerical error writing: %s std. dev. for "
                  "bin %u",
                  data_[index].title.c_str(), bin);
          painCave.isFatal = 1;
          simError();
        } else {
          rnemdFile_ << "\t" << s[0] << "\t" << s[1] << "\t" << s[2];
        }
      }
    }

    void RNEMD::writeArrayErrorBars(int index, unsigned int bin) {
      if (!doRNEMD_) return;
      assert(index >= 0 && index < ENDINDEX);
      assert(bin < nBins_);

      RealType s;
      std::size_t columns = data_[index].accumulatorArray2d[0].size();

      for (std::size_t i = 0; i < columns; i++) {
        std::size_t count = data_[index].accumulatorArray2d[bin][i]->count();

        if (count == 0)
          rnemdFile_ << "\t";
        else {
          dynamic_cast<Accumulator*>(data_[index].accumulatorArray2d[bin][i])
              ->get95percentConfidenceInterval(s);

          if (std::isinf(s) || std::isnan(s)) {
            sprintf(
                painCave.errMsg,
                "RNEMD detected a numerical error writing: %s std. dev. for "
                "bin %u, "
                "column %u",
                data_[index].title.c_str(), bin, static_cast<unsigned int>(i));
            painCave.isFatal = 1;
            simError();
          } else {
            rnemdFile_ << "\t" << s;
          }
        }
      }
    }
  }  // namespace RNEMD
}  // namespace OpenMD
