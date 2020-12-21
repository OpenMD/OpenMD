/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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

#include <algorithm>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef IS_MPI
#include <mpi.h>
#endif

#include "brains/Thermo.hpp"
#include "brains/ForceManager.hpp"
#include "io/Globals.hpp"
#include "math/ConvexHull.hpp"
#include "math/Polynomial.hpp"
#include "math/SquareMatrix3.hpp"
#include "math/Vector.hpp"
#include "math/Vector3.hpp"
#include "primitives/Molecule.hpp"
#include "primitives/StuntDouble.hpp"
#include "rnemd/RNEMD.hpp"
#include "rnemd/RNEMDParameters.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/Accumulator.hpp"
#include "utils/Constants.hpp"
#include "utils/Tuple.hpp"

#define HONKING_LARGE_VALUE 1.0e10

using namespace std;

namespace OpenMD {

  RNEMD::RNEMD(SimInfo* info) : info_(info),
				evaluator_(info_), seleMan_(info_),
				evaluatorA_(info_), seleManA_(info_),
				evaluatorB_(info_), seleManB_(info_),
				commonA_(info_), commonB_(info_),
				outputEvaluator_(info_), outputSeleMan_(info_),
				usePeriodicBoundaryConditions_(info_->getSimParams()->getUsePeriodicBoundaryConditions()),
				hasDividingArea_(false), hasData_(false) {

    trialCount_ = 0;
    failTrialCount_ = 0;
    failRootCount_ = 0;

    Globals* simParams = info->getSimParams();
    RNEMDParameters* rnemdParams = simParams->getRNEMDParameters();

    doRNEMD_ = rnemdParams->getUseRNEMD();
    if (!doRNEMD_) return;

    stringToMethod_["Swap"]  = rnemdSwap;
    stringToMethod_["NIVS"]  = rnemdNIVS;
    stringToMethod_["VSS"]   = rnemdVSS;

    const std::string methStr = rnemdParams->getMethod();
    rnemdMethod_ = stringToMethod_.find(methStr)->second;

    stringToFluxType_["KE"]  = rnemdKE;
    stringToFluxType_["Px"]  = rnemdPx;
    stringToFluxType_["Py"]  = rnemdPy;
    stringToFluxType_["Pz"]  = rnemdPz;
    stringToFluxType_["Pvector"]  = rnemdPvector;
    stringToFluxType_["Lx"]  = rnemdLx;
    stringToFluxType_["Ly"]  = rnemdLy;
    stringToFluxType_["Lz"]  = rnemdLz;
    stringToFluxType_["Lvector"]  = rnemdLvector;
    stringToFluxType_["Current"]  = rnemdCurrent;
    stringToFluxType_["Single"]   = rnemdSingle;
    stringToFluxType_["KE+Px"]  = rnemdKePx;
    stringToFluxType_["KE+Py"]  = rnemdKePy;
    stringToFluxType_["KE+Pvector"]  = rnemdKePvector;
    stringToFluxType_["KE+Current"]  = rnemdKeCurrent;
    stringToFluxType_["KE+Lx"]  = rnemdKeLx;
    stringToFluxType_["KE+Ly"]  = rnemdKeLy;
    stringToFluxType_["KE+Lz"]  = rnemdKeLz;
    stringToFluxType_["KE+Lvector"]  = rnemdKeLvector;

    bool hasFluxType = rnemdParams->haveFluxType();
    std::string fluxStr;

    if (hasFluxType) {
      fluxStr = rnemdParams->getFluxType();
      rnemdFluxType_ = stringToFluxType_.find(fluxStr)->second;

    } else {
      std::string allowedFluxTypes;
      int currentLineLength = 0;

      for (std::map<std::string, RNEMDFluxType>::iterator fluxStrIter = stringToFluxType_.begin();
	   fluxStrIter != stringToFluxType_.end(); ++fluxStrIter) {
	allowedFluxTypes += fluxStrIter->first + ", ";
	currentLineLength += fluxStrIter->first.length() + 2;

	if (currentLineLength >= 50) {
	  allowedFluxTypes += "\n\t\t";
	  currentLineLength = 0;
	}
      }

      allowedFluxTypes.erase(allowedFluxTypes.length() - 2, 2);

      sprintf(painCave.errMsg,
	      "RNEMD: No fluxType was set in the omd file. This parameter\n"
	      "\tmust be set to use RNEMD, and can take any of these values:\n"
	      "\t\t%s.\n",
	      allowedFluxTypes.c_str());
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    stringToPrivilegedAxis_["x"] = rnemdX;
    stringToPrivilegedAxis_["y"] = rnemdY;
    stringToPrivilegedAxis_["z"] = rnemdZ;

    const std::string privAxis = rnemdParams->getPrivilegedAxis();
    rnemdPrivilegedAxis_ = stringToPrivilegedAxis_.find(privAxis)->second;

    runTime_ = simParams->getRunTime();
    statusTime_ = simParams->getStatusTime();

    rnemdObjectSelection_ = rnemdParams->getObjectSelection();

    bool hasKineticFlux = rnemdParams->haveKineticFlux();
    bool hasMomentumFlux = rnemdParams->haveMomentumFlux();
    bool hasCurrentDensity = rnemdParams->haveCurrentDensity();
    bool hasMomentumFluxVector = rnemdParams->haveMomentumFluxVector();
    bool hasAngularMomentumFlux = rnemdParams->haveAngularMomentumFlux();
    bool hasAngularMomentumFluxVector = rnemdParams->haveAngularMomentumFluxVector();
    hasSelectionA_ = rnemdParams->haveSelectionA();
    hasSelectionB_ = rnemdParams->haveSelectionB();
    bool hasSlabWidth = rnemdParams->haveSlabWidth();
    bool hasSlabACenter = rnemdParams->haveSlabACenter();
    bool hasSlabBCenter = rnemdParams->haveSlabBCenter();
    bool hasSphereARadius = rnemdParams->haveSphereARadius();
    hasSphereBRadius_ = rnemdParams->haveSphereBRadius();

    hasDividingArea_ = rnemdParams->haveDividingArea();
    dividingArea_ = rnemdParams->getDividingArea();

    bool hasCoordinateOrigin = rnemdParams->haveCoordinateOrigin();
    bool hasOutputFileName = rnemdParams->haveOutputFileName();
    bool hasOutputFields = rnemdParams->haveOutputFields();

    hasOutputSelection_ = rnemdParams->haveOutputSelection();
    if (hasOutputSelection_) {
      outputSelection_ = rnemdParams->getOutputSelection();
    } else {
      outputSelection_ = rnemdObjectSelection_;
    }

    switch(rnemdPrivilegedAxis_) {
    case 0:
      rnemdAxisLabel_ = "x";
      break;
    case 1:
      rnemdAxisLabel_ = "y";
      break;
    case 2:
    default:
      rnemdAxisLabel_ = "z";
      break;
    }

    bool methodFluxMismatch = false;
    bool hasCorrectFlux = false;
    switch(rnemdMethod_) {
    case rnemdSwap:
      switch (rnemdFluxType_) {
      case rnemdKE:
        hasCorrectFlux = hasKineticFlux;
        break;
      case rnemdPx:
      case rnemdPy:
      case rnemdPz:
        hasCorrectFlux = hasMomentumFlux;
        break;
      default :
        methodFluxMismatch = true;
        break;
      }
      break;
    case rnemdNIVS:
      switch (rnemdFluxType_) {
      case rnemdKE:
      case rnemdRotKE:
      case rnemdFullKE:
        hasCorrectFlux = hasKineticFlux;
        break;
      case rnemdPx:
      case rnemdPy:
      case rnemdPz:
        hasCorrectFlux = hasMomentumFlux;
        break;
      case rnemdKePx:
      case rnemdKePy:
        hasCorrectFlux = hasMomentumFlux && hasKineticFlux;
        break;
      default:
        methodFluxMismatch = true;
        break;
      }
      break;
    case rnemdVSS:
      switch (rnemdFluxType_) {
      case rnemdKE:
      case rnemdRotKE:
      case rnemdFullKE:
        hasCorrectFlux = hasKineticFlux;
        break;
      case rnemdPx:
      case rnemdPy:
      case rnemdPz:
        hasCorrectFlux = hasMomentumFlux;
        break;
      case rnemdLx:
      case rnemdLy:
      case rnemdLz:
        hasCorrectFlux = hasAngularMomentumFlux;
        break;
      case rnemdPvector:
        hasCorrectFlux = hasMomentumFluxVector;
        break;
      case rnemdLvector:
        hasCorrectFlux = hasAngularMomentumFluxVector;
        break;
      case rnemdCurrent:
      case rnemdSingle:
        hasCorrectFlux = hasCurrentDensity;
        break;
      case rnemdKePx:
      case rnemdKePy:
        hasCorrectFlux = hasMomentumFlux && hasKineticFlux;
        break;
      case rnemdKeLx:
      case rnemdKeLy:
      case rnemdKeLz:
        hasCorrectFlux = hasAngularMomentumFlux && hasKineticFlux;
        break;
      case rnemdKeCurrent:
        hasCorrectFlux = hasCurrentDensity && hasKineticFlux;
        break;
      case rnemdKePvector:
        hasCorrectFlux = hasMomentumFluxVector && hasKineticFlux;
        break;
      case rnemdKeLvector:
        hasCorrectFlux = hasAngularMomentumFluxVector && hasKineticFlux;
        break;
      default:
        methodFluxMismatch = true;
        break;
      }
    default:
      break;
    }

    if (methodFluxMismatch) {
      sprintf(painCave.errMsg,
              "RNEMD: The current method,\n"
              "\t\t%s\n"
              "\tcannot be used with the current flux type, %s\n",
              methStr.c_str(), fluxStr.c_str());
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }
    if (!hasCorrectFlux) {
      sprintf(painCave.errMsg,
              "RNEMD: The current method, %s, and flux type, %s,\n"
              "\tdid not have the correct flux value specified. Options\n"
              "\tinclude: kineticFlux, momentumFlux, angularMomentumFlux,\n"
              "\tmomentumFluxVector, angularMomentumFluxVector, and\n"
              "\tcurrentDensity\n",
              methStr.c_str(), fluxStr.c_str());
      painCave.isFatal = 1;
      painCave.severity = OPENMD_ERROR;
      simError();
    }

    if (hasKineticFlux) {
      // convert the kcal / mol / Angstroms^2 / fs values in the md file
      // into  amu / fs^3:
      kineticFlux_ = rnemdParams->getKineticFlux()
        * Constants::energyConvert;
    } else {
      kineticFlux_ = 0.0;
    }
    if (hasMomentumFluxVector) {
      std::vector<RealType> mf = rnemdParams->getMomentumFluxVector();
      if (mf.size() != 3) {
        sprintf(painCave.errMsg,
                "RNEMD: Incorrect number of parameters specified for momentumFluxVector.\n"
                "\tthere should be 3 parameters, but %lu were specified.\n",
                mf.size());
        painCave.isFatal = 1;
        simError();
      }
      momentumFluxVector_.x() = mf[0];
      momentumFluxVector_.y() = mf[1];
      momentumFluxVector_.z() = mf[2];
    } else {
      momentumFluxVector_ = V3Zero;
      if (hasMomentumFlux) {
        RealType momentumFlux = rnemdParams->getMomentumFlux();
        switch (rnemdFluxType_) {
        case rnemdPx:
          momentumFluxVector_.x() = momentumFlux;
          break;
        case rnemdPy:
          momentumFluxVector_.y() = momentumFlux;
          break;
        case rnemdPz:
          momentumFluxVector_.z() = momentumFlux;
          break;
        case rnemdKePx:
          momentumFluxVector_.x() = momentumFlux;
          break;
        case rnemdKePy:
          momentumFluxVector_.y() = momentumFlux;
          break;
        default:
          break;
        }
      }
    }
    if (hasAngularMomentumFluxVector) {
      std::vector<RealType> amf = rnemdParams->getAngularMomentumFluxVector();
      if (amf.size() != 3) {
        sprintf(painCave.errMsg,
                "RNEMD: Incorrect number of parameters specified for angularMomentumFluxVector.\n"
                "\tthere should be 3 parameters, but %lu were specified.\n",
                amf.size());
        painCave.isFatal = 1;
        simError();
      }
      angularMomentumFluxVector_.x() = amf[0];
      angularMomentumFluxVector_.y() = amf[1];
      angularMomentumFluxVector_.z() = amf[2];
    } else {
      angularMomentumFluxVector_ = V3Zero;
      if (hasAngularMomentumFlux) {
        RealType angularMomentumFlux = rnemdParams->getAngularMomentumFlux();
        switch (rnemdFluxType_) {
        case rnemdLx:
          angularMomentumFluxVector_.x() = angularMomentumFlux;
          break;
        case rnemdLy:
          angularMomentumFluxVector_.y() = angularMomentumFlux;
          break;
        case rnemdLz:
          angularMomentumFluxVector_.z() = angularMomentumFlux;
          break;
        case rnemdKeLx:
          angularMomentumFluxVector_.x() = angularMomentumFlux;
          break;
        case rnemdKeLy:
          angularMomentumFluxVector_.y() = angularMomentumFlux;
          break;
        case rnemdKeLz:
          angularMomentumFluxVector_.z() = angularMomentumFlux;
          break;
        default:
          break;
        }
      }
    }
    if (hasCurrentDensity) {
      // convert the Amp m^-2 values in the omd file into electrons fs^-1 A^-2:
      currentDensity_ = rnemdParams->getCurrentDensity()
        * Constants::currentDensityConvert;
    } else {
      currentDensity_ = 0.0;
    }

    if (hasCoordinateOrigin) {
      std::vector<RealType> co = rnemdParams->getCoordinateOrigin();
      if (co.size() != 3) {
        sprintf(painCave.errMsg,
                "RNEMD: Incorrect number of parameters specified for coordinateOrigin.\n"
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
    std::copy(osTypes.begin(), osTypes.end(), std::back_inserter(outputTypes_));

    nBins_ = rnemdParams->getOutputBins();
    binWidth_ = rnemdParams->getOutputBinWidth();

    data_.resize(RNEMD::ENDINDEX);

    OutputData z;
    z.units =  "Angstroms";
    z.title =  rnemdAxisLabel_;
    z.dataType = "RealType";
    z.accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      z.accumulator.push_back( new Accumulator() );
    data_[Z] = z;
    outputMap_["Z"] =  Z;

    OutputData r;
    r.units =  "Angstroms";
    r.title =  "R";
    r.dataType = "RealType";
    r.accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      r.accumulator.push_back( new Accumulator() );
    data_[R] = r;
    outputMap_["R"] =  R;

    OutputData temperature;
    temperature.units =  "K";
    temperature.title =  "Temperature";
    temperature.dataType = "RealType";
    temperature.accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      temperature.accumulator.push_back( new Accumulator() );
    data_[TEMPERATURE] = temperature;
    outputMap_["TEMPERATURE"] =  TEMPERATURE;

    OutputData velocity;
    velocity.units = "angstroms/fs";
    velocity.title =  "Velocity";
    velocity.dataType = "Vector3d";
    velocity.accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      velocity.accumulator.push_back( new VectorAccumulator() );
    data_[VELOCITY] = velocity;
    outputMap_["VELOCITY"] = VELOCITY;

    OutputData angularVelocity;
    angularVelocity.units = "angstroms^2/fs";
    angularVelocity.title =  "AngularVelocity";
    angularVelocity.dataType = "Vector3d";
    angularVelocity.accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      angularVelocity.accumulator.push_back( new VectorAccumulator() );
    data_[ANGULARVELOCITY] = angularVelocity;
    outputMap_["ANGULARVELOCITY"] = ANGULARVELOCITY;

    OutputData density;
    density.units =  "g cm^-3";
    density.title =  "Density";
    density.dataType = "RealType";
    density.accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      density.accumulator.push_back( new Accumulator() );
    data_[DENSITY] = density;
    outputMap_["DENSITY"] =  DENSITY;

    OutputData activity;
    activity.units = "unitless";
    activity.title =  "Activity";
    activity.dataType = "Array2d";
    unsigned int nTypes = outputTypes_.size();
    activity.accumulatorArray2d.resize(nBins_);
    for (unsigned int i = 0; i < nBins_; i++) {
      activity.accumulatorArray2d[i].resize(nTypes);
      for (unsigned int j = 0 ; j < nTypes; j++) {
        activity.accumulatorArray2d[i][j] = new Accumulator();
      }
    }
    data_[ACTIVITY] = (activity);
    outputMap_["ACTIVITY"] =  ACTIVITY;

    OutputData eField;
    eField.units =  "kcal/mol/angstroms/e";
    eField.title =  "Electric Field";
    eField.dataType = "Vector3d";
    eField.accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      eField.accumulator.push_back( new VectorAccumulator() );
    data_[ELECTRICFIELD] = eField;
    outputMap_["ELECTRICFIELD"] =  ELECTRICFIELD;

    OutputData ePot;
    ePot.units =  "kcal/mol/e";
    ePot.title =  "Electrostatic Potential";
    ePot.dataType = "RealType";
    ePot.accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      ePot.accumulator.push_back( new Accumulator() );
    data_[ELECTROSTATICPOTENTIAL] = ePot;
    outputMap_["ELECTROSTATICPOTENTIAL"] =  ELECTROSTATICPOTENTIAL;

    OutputData currentDensity;
    currentDensity.units = "e/angstroms^2/fs";
    currentDensity.title =  "Current Density";
    currentDensity.dataType = "Array2d";
    nTypes = outputTypes_.size();
    currentDensity.accumulatorArray2d.resize(nBins_);
    for (unsigned int i = 0; i < nBins_; i++) {
      currentDensity.accumulatorArray2d[i].resize(nTypes);
      for (unsigned int j = 0 ; j < nTypes; j++) {
        currentDensity.accumulatorArray2d[i][j] = new Accumulator();
      }
    }
    data_[CURRENTDENSITY] = (currentDensity);
    outputMap_["CURRENTDENSITY"] =  CURRENTDENSITY;

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
      case rnemdCurrent:
      case rnemdSingle:
      case rnemdKeCurrent:
        outputMask_.set(TEMPERATURE);
        outputMask_.set(VELOCITY);
        outputMask_.set(DENSITY);
        outputMask_.set(ACTIVITY);
        outputMask_.set(ELECTRICFIELD);
        outputMask_.set(ELECTROSTATICPOTENTIAL);
	outputMask_.set(CURRENTDENSITY);
      default:
        break;
      }
    }

    if (hasOutputFileName) {
      rnemdFileName_ = rnemdParams->getOutputFileName();
    } else {
      rnemdFileName_ = getPrefix(info->getFinalConfigFileName()) + ".rnemd";
    }

    // Exchange time should not be less than time step and should be a multiple of dt
    exchangeTime_ = rnemdParams->getExchangeTime();
    RealType dt = info->getSimParams()->getDt();
    RealType newET = std::ceil(exchangeTime_ / dt) * dt;

    if ( std::fabs(newET - exchangeTime_) > 1e-6 ) {
      sprintf(painCave.errMsg,
	      "RNEMD: The exchangeTime was reset to %lf,\n"
	      "\t\twhich is a multiple of dt, %lf.\n",
	      newET, dt);
      painCave.isFatal = 0;
      painCave.severity = OPENMD_WARNING;
      simError();
      exchangeTime_ = newET;
    }

    currentSnap_ = info->getSnapshotManager()->getCurrentSnapshot();
    hmat_ = currentSnap_->getHmat();

    // total exchange sums are zeroed out at the beginning:
    kineticExchange_ = 0.0;
    momentumExchange_ = V3Zero;
    particleFlux_h_ = V3Zero;
    particleFlux_c_ = V3Zero;
    angularMomentumExchange_ = V3Zero;

    std::ostringstream selectionAstream;
    std::ostringstream selectionBstream;

    if (rnemdFluxType_ == rnemdSingle) {
      if (hasSelectionA_) {
        selectionA_ = rnemdParams->getSelectionA();
      } else {
        if (usePeriodicBoundaryConditions_) {
          if (hasSlabWidth)
            slabWidth_ = rnemdParams->getSlabWidth();
          else
            slabWidth_ = hmat_(rnemdPrivilegedAxis_,rnemdPrivilegedAxis_) / 2.0;

          if (hasSlabACenter) {
            slabACenter_ = rnemdParams->getSlabACenter();
          } else {
            slabACenter_ = 0.0;
          }
          selectionAstream << "select wrappedz > "
                           << slabACenter_ - 0.5*slabWidth_
                           <<  " && wrappedz < "
                           << -slabACenter_ + 0.5*slabWidth_;
          selectionA_ = selectionAstream.str();
        }
      }

      if (hasSelectionB_) {
        sprintf(painCave.errMsg,
                "RNEMD: The rnemdSingle flux type only allows for a selectionA to be specified.\n");
        painCave.isFatal = 0;
        painCave.severity = OPENMD_WARNING;
        simError();
      } else selectionB_ = std::string("select none");
    } else {
      if (hasSelectionA_) {
        selectionA_ = rnemdParams->getSelectionA();
      } else {
        if (usePeriodicBoundaryConditions_) {
          if (hasSlabWidth)
            slabWidth_ = rnemdParams->getSlabWidth();
          else
            slabWidth_ = hmat_(rnemdPrivilegedAxis_,rnemdPrivilegedAxis_) / 10.0;

          if (hasSlabACenter) slabACenter_ = rnemdParams->getSlabACenter();
          else slabACenter_ = 0.0;

          selectionA_ = this->setSelection(slabACenter_);
        } else {
          if (hasSphereARadius)
            sphereARadius_ = rnemdParams->getSphereARadius();
          else {
            // use an initial guess to the size of the inner slab to be 1/10 the
            // radius of an approximately spherical hull:
            Thermo thermo(info);
            RealType hVol = thermo.getHullVolume();
            sphereARadius_ = 0.1 * pow((3.0 * hVol / (4.0 * Constants::PI)), 1.0/3.0);
          }
          selectionAstream << "select r < " << sphereARadius_;
          selectionA_ = selectionAstream.str();
        }
      }

      if (hasSelectionB_) {
        selectionB_ = rnemdParams->getSelectionB();
      } else {
        if (usePeriodicBoundaryConditions_) {
          if (hasSlabWidth)
            slabWidth_ = rnemdParams->getSlabWidth();
          else
            slabWidth_ = hmat_(rnemdPrivilegedAxis_,rnemdPrivilegedAxis_) / 10.0;

          if (hasSlabBCenter) slabBCenter_ = rnemdParams->getSlabBCenter();
          else slabBCenter_ = hmat_(rnemdPrivilegedAxis_,rnemdPrivilegedAxis_) / 2.0;

          selectionB_ = this->setSelection(slabBCenter_);

        } else {
          if (hasSphereBRadius_) {
            sphereBRadius_ = rnemdParams->getSphereBRadius();
            selectionBstream << "select r > " << sphereBRadius_;
            selectionB_ = selectionBstream.str();
          } else {
            selectionB_ = "select hull";
            BisHull_ = true;
            hasSelectionB_ = true;
          }
        }
      }
    }

    // static object evaluators:
    evaluator_.loadScriptString(rnemdObjectSelection_);
    if (!evaluator_.isDynamic())
      seleMan_.setSelectionSet(evaluator_.evaluate());

    evaluatorA_.loadScriptString(selectionA_);
    if (!evaluatorA_.isDynamic())
      seleManA_.setSelectionSet(evaluatorA_.evaluate());

    evaluatorB_.loadScriptString(selectionB_);
    if (!evaluatorB_.isDynamic())
      seleManB_.setSelectionSet(evaluatorB_.evaluate());

    // do some sanity checking
    int selectionCount = seleMan_.getSelectionCount();
    int nIntegrable = info->getNGlobalIntegrableObjects();

    if (selectionCount > nIntegrable) {
      sprintf(painCave.errMsg,
              "RNEMD: The current objectSelection,\n"
              "\t\t%s\n"
              "\thas resulted in %d selected objects.  However,\n"
              "\tthe total number of integrable objects in the system\n"
              "\tis only %d.  This is almost certainly not what you want\n"
              "\tto do.  A likely cause of this is forgetting the _RB_0\n"
              "\tselector in the selection script!\n",
              rnemdObjectSelection_.c_str(),
              selectionCount, nIntegrable);
      painCave.isFatal = 0;
      painCave.severity = OPENMD_WARNING;
      simError();
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
              "slabAcenter was set to %0.2f. In the wrapped coordinates\n"
              "\t[-Hmat/2, +Hmat/2], this has been remapped to %0.2f.\n",
              slabCenter, tempSlabCenter[rnemdPrivilegedAxis_]);
      painCave.isFatal = 0;
      painCave.severity = OPENMD_WARNING;
      simError();

      slabCenter = tempSlabCenter[rnemdPrivilegedAxis_];
    }

    RealType leftSlabBoundary = slabCenter - 0.5*slabWidth_;
    RealType rightSlabBoundary = slabCenter + 0.5*slabWidth_;

    std::ostringstream selectionStream;

    selectionStream << "select wrappedz > " << leftSlabBoundary;

    if ( (leftSlabBoundary < -hmat_2) || (rightSlabBoundary > hmat_2) )
      selectionStream <<  " || wrappedz < " << -leftSlabBoundary;
    else
      selectionStream <<  " && wrappedz < " << rightSlabBoundary;

    return selectionStream.str();
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
      if ( !data.accumulatorArray2d.empty() )
        MemoryUtils::deletePointers(data.accumulatorArray2d);
      else
        MemoryUtils::deletePointers(data.accumulator);
    }
  }

  void RNEMD::doSwap(SelectionManager& smanA, SelectionManager& smanB) {

    if (!doRNEMD_) return;
    int selei;
    int selej;

    StuntDouble* sd;

    RealType min_val(0.0);
    int min_found = 0;
    StuntDouble* min_sd = NULL;

    RealType max_val(0.0);
    int max_found = 0;
    StuntDouble* max_sd = NULL;

    for (sd = seleManA_.beginSelected(selei); sd != NULL;
         sd = seleManA_.nextSelected(selei)) {

      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();
      RealType value(0.0);

      switch(rnemdFluxType_) {
      case rnemdKE :

        value = mass * vel.lengthSquare();

        if (sd->isDirectional()) {
          Vector3d angMom = sd->getJ();
          Mat3x3d I = sd->getI();

          if (sd->isLinear()) {
            int i = sd->linearAxis();
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;
            value += angMom[j] * angMom[j] / I(j, j) +
              angMom[k] * angMom[k] / I(k, k);
          } else {
            value += angMom[0]*angMom[0]/I(0, 0)
              + angMom[1]*angMom[1]/I(1, 1)
              + angMom[2]*angMom[2]/I(2, 2);
          }
        } // angular momenta exchange enabled
        value *= 0.5;
        break;
      case rnemdPx :
        value = mass * vel[0];
        break;
      case rnemdPy :
        value = mass * vel[1];
        break;
      case rnemdPz :
        value = mass * vel[2];
        break;
      default :
        break;
      }
      if (!max_found) {
        max_val = value;
        max_sd = sd;
        max_found = 1;
      } else {
        if (max_val < value) {
          max_val = value;
          max_sd = sd;
        }
      }
    }

    for (sd = seleManB_.beginSelected(selej); sd != NULL;
         sd = seleManB_.nextSelected(selej)) {

      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();
      RealType value(0.0);

      switch(rnemdFluxType_) {
      case rnemdKE :

        value = mass * vel.lengthSquare();

        if (sd->isDirectional()) {
          Vector3d angMom = sd->getJ();
          Mat3x3d I = sd->getI();

          if (sd->isLinear()) {
            int i = sd->linearAxis();
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;
            value += angMom[j] * angMom[j] / I(j, j) +
              angMom[k] * angMom[k] / I(k, k);
          } else {
            value += angMom[0]*angMom[0]/I(0, 0)
              + angMom[1]*angMom[1]/I(1, 1)
              + angMom[2]*angMom[2]/I(2, 2);
          }
        } // angular momenta exchange enabled
        value *= 0.5;
        break;
      case rnemdPx :
        value = mass * vel[0];
        break;
      case rnemdPy :
        value = mass * vel[1];
        break;
      case rnemdPz :
        value = mass * vel[2];
        break;
      default :
        break;
      }

      if (!min_found) {
        min_val = value;
        min_sd = sd;
        min_found = 1;
      } else {
        if (min_val > value) {
          min_val = value;
          min_sd = sd;
        }
      }
    }

#ifdef IS_MPI
    int worldRank;
    MPI_Comm_rank( MPI_COMM_WORLD, &worldRank);

    int my_min_found = min_found;
    int my_max_found = max_found;

    // Even if we didn't find a minimum, did someone else?
    MPI_Allreduce(&my_min_found, &min_found, 1, MPI_INT, MPI_LOR,
                  MPI_COMM_WORLD);
    // Even if we didn't find a maximum, did someone else?
    MPI_Allreduce(&my_max_found, &max_found, 1, MPI_INT, MPI_LOR,
                  MPI_COMM_WORLD);
#endif

    if (max_found && min_found) {

#ifdef IS_MPI
      struct {
	RealType val;
	int rank;
      } max_vals, min_vals;

      if (my_min_found) {
        min_vals.val = min_val;
      } else {
        min_vals.val = HONKING_LARGE_VALUE;
      }
      min_vals.rank = worldRank;

      // Who had the minimum?
      MPI_Allreduce(&min_vals, &min_vals,
                    1, MPI_REALTYPE_INT, MPI_MINLOC, MPI_COMM_WORLD);
      min_val = min_vals.val;

      if (my_max_found) {
        max_vals.val = max_val;
      } else {
        max_vals.val = -HONKING_LARGE_VALUE;
      }
      max_vals.rank = worldRank;

      // Who had the maximum?
      MPI_Allreduce(&max_vals, &max_vals,
                    1, MPI_REALTYPE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      max_val = max_vals.val;
#endif

      if (min_val < max_val) {

#ifdef IS_MPI
        if (max_vals.rank == worldRank && min_vals.rank == worldRank) {
          // I have both maximum and minimum, so proceed like a single
          // processor version:
#endif

          Vector3d min_vel = min_sd->getVel();
          Vector3d max_vel = max_sd->getVel();
          RealType temp_vel;

          switch(rnemdFluxType_) {
          case rnemdKE :
            min_sd->setVel(max_vel);
            max_sd->setVel(min_vel);
	    if (min_sd->isDirectional() && max_sd->isDirectional()) {
              Vector3d min_angMom = min_sd->getJ();
              Vector3d max_angMom = max_sd->getJ();
              min_sd->setJ(max_angMom);
              max_sd->setJ(min_angMom);
	    }// angular momenta exchange enabled
	    // assumes same rigid body identity
            break;
          case rnemdPx :
            temp_vel = min_vel.x();
            min_vel.x() = max_vel.x();
            max_vel.x() = temp_vel;
            min_sd->setVel(min_vel);
            max_sd->setVel(max_vel);
            break;
          case rnemdPy :
            temp_vel = min_vel.y();
            min_vel.y() = max_vel.y();
            max_vel.y() = temp_vel;
            min_sd->setVel(min_vel);
            max_sd->setVel(max_vel);
            break;
          case rnemdPz :
            temp_vel = min_vel.z();
            min_vel.z() = max_vel.z();
            max_vel.z() = temp_vel;
            min_sd->setVel(min_vel);
            max_sd->setVel(max_vel);
            break;
          default :
            break;
          }

#ifdef IS_MPI
          // the rest of the cases only apply in parallel simulations:
        } else if (max_vals.rank == worldRank) {
          // I had the max, but not the minimum

          Vector3d min_vel;
          Vector3d max_vel = max_sd->getVel();
          MPI_Status status;

          // point-to-point swap of the velocity vector
          MPI_Sendrecv(max_vel.getArrayPointer(), 3, MPI_REALTYPE,
                       min_vals.rank, 0,
                       min_vel.getArrayPointer(), 3, MPI_REALTYPE,
                       min_vals.rank, 0, MPI_COMM_WORLD, &status);

          switch(rnemdFluxType_) {
          case rnemdKE :
            max_sd->setVel(min_vel);
            // angular momenta exchange enabled
            if (max_sd->isDirectional()) {
              Vector3d min_angMom;
              Vector3d max_angMom = max_sd->getJ();

              // point-to-point swap of the angular momentum vector
              MPI_Sendrecv(max_angMom.getArrayPointer(), 3,
                           MPI_REALTYPE, min_vals.rank, 1,
                           min_angMom.getArrayPointer(), 3,
                           MPI_REALTYPE, min_vals.rank, 1,
                           MPI_COMM_WORLD, &status);

              max_sd->setJ(min_angMom);
	    }
            break;
          case rnemdPx :
            max_vel.x() = min_vel.x();
            max_sd->setVel(max_vel);
            break;
          case rnemdPy :
            max_vel.y() = min_vel.y();
            max_sd->setVel(max_vel);
            break;
          case rnemdPz :
            max_vel.z() = min_vel.z();
            max_sd->setVel(max_vel);
            break;
          default :
            break;
          }
        } else if (min_vals.rank == worldRank) {
          // I had the minimum but not the maximum:

          Vector3d max_vel;
          Vector3d min_vel = min_sd->getVel();
          MPI_Status status;

          // point-to-point swap of the velocity vector
          MPI_Sendrecv(min_vel.getArrayPointer(), 3, MPI_REALTYPE,
                       max_vals.rank, 0,
                       max_vel.getArrayPointer(), 3, MPI_REALTYPE,
                       max_vals.rank, 0, MPI_COMM_WORLD, &status);

          switch(rnemdFluxType_) {
          case rnemdKE :
            min_sd->setVel(max_vel);
            // angular momenta exchange enabled
            if (min_sd->isDirectional()) {
              Vector3d min_angMom = min_sd->getJ();
              Vector3d max_angMom;

              // point-to-point swap of the angular momentum vector
              MPI_Sendrecv(min_angMom.getArrayPointer(), 3,
                           MPI_REALTYPE, max_vals.rank, 1,
                           max_angMom.getArrayPointer(), 3,
                           MPI_REALTYPE, max_vals.rank, 1,
                           MPI_COMM_WORLD, &status);

              min_sd->setJ(max_angMom);
            }
            break;
          case rnemdPx :
            min_vel.x() = max_vel.x();
            min_sd->setVel(min_vel);
            break;
          case rnemdPy :
            min_vel.y() = max_vel.y();
            min_sd->setVel(min_vel);
            break;
          case rnemdPz :
            min_vel.z() = max_vel.z();
            min_sd->setVel(min_vel);
            break;
          default :
            break;
          }
        }
#endif

        switch(rnemdFluxType_) {
        case rnemdKE:
          kineticExchange_ += max_val - min_val;
          break;
        case rnemdPx:
          momentumExchange_.x() += max_val - min_val;
          break;
        case rnemdPy:
          momentumExchange_.y() += max_val - min_val;
          break;
        case rnemdPz:
          momentumExchange_.z() += max_val - min_val;
          break;
        default:
          break;
        }
      } else {
        sprintf(painCave.errMsg,
                "RNEMD::doSwap exchange NOT performed because min_val > max_val\n");
        painCave.isFatal = 0;
        painCave.severity = OPENMD_INFO;
        simError();
        failTrialCount_++;
      }
    } else {
      sprintf(painCave.errMsg,
              "RNEMD::doSwap exchange NOT performed because selected object\n"
	      "\twas not present in at least one of the two slabs.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      failTrialCount_++;
    }
  }

  void RNEMD::doNIVS(SelectionManager& smanA, SelectionManager& smanB) {

    if (!doRNEMD_) return;
    int selei;
    int selej;

    StuntDouble* sd;

    vector<StuntDouble*> hotBin, coldBin;

    RealType Phx = 0.0;
    RealType Phy = 0.0;
    RealType Phz = 0.0;
    RealType Khx = 0.0;
    RealType Khy = 0.0;
    RealType Khz = 0.0;
    RealType Khw = 0.0;
    RealType Pcx = 0.0;
    RealType Pcy = 0.0;
    RealType Pcz = 0.0;
    RealType Kcx = 0.0;
    RealType Kcy = 0.0;
    RealType Kcz = 0.0;
    RealType Kcw = 0.0;

    for (sd = smanA.beginSelected(selei); sd != NULL;
         sd = smanA.nextSelected(selei)) {

      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);


      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();

      hotBin.push_back(sd);
      Phx += mass * vel.x();
      Phy += mass * vel.y();
      Phz += mass * vel.z();
      Khx += mass * vel.x() * vel.x();
      Khy += mass * vel.y() * vel.y();
      Khz += mass * vel.z() * vel.z();
      if (sd->isDirectional()) {
        Vector3d angMom = sd->getJ();
        Mat3x3d I = sd->getI();
        if (sd->isLinear()) {
          int i = sd->linearAxis();
          int j = (i + 1) % 3;
          int k = (i + 2) % 3;
          Khw += angMom[j] * angMom[j] / I(j, j) +
            angMom[k] * angMom[k] / I(k, k);
        } else {
          Khw += angMom[0]*angMom[0]/I(0, 0)
            + angMom[1]*angMom[1]/I(1, 1)
            + angMom[2]*angMom[2]/I(2, 2);
        }
      }
    }
    for (sd = smanB.beginSelected(selej); sd != NULL;
         sd = smanB.nextSelected(selej)) {
      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();

      coldBin.push_back(sd);
      Pcx += mass * vel.x();
      Pcy += mass * vel.y();
      Pcz += mass * vel.z();
      Kcx += mass * vel.x() * vel.x();
      Kcy += mass * vel.y() * vel.y();
      Kcz += mass * vel.z() * vel.z();
      if (sd->isDirectional()) {
        Vector3d angMom = sd->getJ();
        Mat3x3d I = sd->getI();
        if (sd->isLinear()) {
          int i = sd->linearAxis();
          int j = (i + 1) % 3;
          int k = (i + 2) % 3;
          Kcw += angMom[j] * angMom[j] / I(j, j) +
            angMom[k] * angMom[k] / I(k, k);
        } else {
          Kcw += angMom[0]*angMom[0]/I(0, 0)
            + angMom[1]*angMom[1]/I(1, 1)
            + angMom[2]*angMom[2]/I(2, 2);
        }
      }
    }

    Khx *= 0.5;
    Khy *= 0.5;
    Khz *= 0.5;
    Khw *= 0.5;
    Kcx *= 0.5;
    Kcy *= 0.5;
    Kcz *= 0.5;
    Kcw *= 0.5;

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &Phx, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Phy, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Phz, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Pcx, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Pcy, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Pcz, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &Khx, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Khy, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Khz, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Khw, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &Kcx, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kcy, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kcz, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kcw, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif

    // solve coldBin coeff's first
    RealType px = Pcx / Phx;
    RealType py = Pcy / Phy;
    RealType pz = Pcz / Phz;
    RealType c(0.0), x(0.0), y(0.0), z(0.0);
    bool successfulScale = false;
    if ((rnemdFluxType_ == rnemdFullKE) ||
	(rnemdFluxType_ == rnemdRotKE)) {
      // may need sanity check Khw & Kcw > 0

      if (rnemdFluxType_ == rnemdFullKE) {
	c = 1.0 - kineticTarget_ / (Kcx + Kcy + Kcz + Kcw);
      } else {
	c = 1.0 - kineticTarget_ / Kcw;
      }

      if ((c > 0.81) && (c < 1.21)) {	// restrict scaling coefficients
	c = sqrt(c);

	RealType w = 0.0;
	if (rnemdFluxType_ ==  rnemdFullKE) {
	  x = 1.0 + px * (1.0 - c);
	  y = 1.0 + py * (1.0 - c);
	  z = 1.0 + pz * (1.0 - c);
	  /* more complicated way
	     w = 1.0 + (Kcw - Kcw * c * c - (c * c * (Kcx + Kcy + Kcz
	     + Khx * px * px + Khy * py * py + Khz * pz * pz)
	     - 2.0 * c * (Khx * px * (1.0 + px) + Khy * py * (1.0 + py)
	     + Khz * pz * (1.0 + pz)) + Khx * px * (2.0 + px)
	     + Khy * py * (2.0 + py) + Khz * pz * (2.0 + pz)
	     - Kcx - Kcy - Kcz)) / Khw; the following is simpler
	  */
	  if ((std::fabs(x - 1.0) < 0.1) && (std::fabs(y - 1.0) < 0.1) &&
	      (std::fabs(z - 1.0) < 0.1)) {
	    w = 1.0 + (kineticTarget_
                       + Khx * (1.0 - x * x) + Khy * (1.0 - y * y)
		       + Khz * (1.0 - z * z)) / Khw;
	  }// no need to calculate w if x, y or z is out of range
	} else {
	  w = 1.0 + kineticTarget_ / Khw;
	}
	if ((w > 0.81) && (w < 1.21)) {// restrict scaling coefficients
	  // if w is in the right range, so should be x, y, z.
	  std::vector<StuntDouble*>::iterator sdi;
	  Vector3d vel;
	  for (sdi = coldBin.begin(); sdi != coldBin.end(); ++sdi) {
	    if (rnemdFluxType_ == rnemdFullKE) {
	      vel = (*sdi)->getVel() * c;
	      (*sdi)->setVel(vel);
	    }
	    if ((*sdi)->isDirectional()) {
	      Vector3d angMom = (*sdi)->getJ() * c;
	      (*sdi)->setJ(angMom);
	    }
	  }
	  w = sqrt(w);
	  for (sdi = hotBin.begin(); sdi != hotBin.end(); ++sdi) {
	    if (rnemdFluxType_ == rnemdFullKE) {
	      vel = (*sdi)->getVel();
	      vel.x() *= x;
	      vel.y() *= y;
	      vel.z() *= z;
	      (*sdi)->setVel(vel);
	    }
	    if ((*sdi)->isDirectional()) {
	      Vector3d angMom = (*sdi)->getJ() * w;
	      (*sdi)->setJ(angMom);
	    }
	  }
	  successfulScale = true;
	  kineticExchange_ += kineticTarget_;
	}
      }
    } else {
      RealType a000(0.0), a110(0.0), c0(0.0);
      RealType a001(0.0), a111(0.0), b01(0.0), b11(0.0), c1(0.0);
      switch(rnemdFluxType_) {
      case rnemdKE :
	/* used hotBin coeff's & only scale x & y dimensions
	   RealType px = Phx / Pcx;
	   RealType py = Phy / Pcy;
	   a110 = Khy;
	   c0 = - Khx - Khy - kineticTarget_;
	   a000 = Khx;
	   a111 = Kcy * py * py;
	   b11 = -2.0 * Kcy * py * (1.0 + py);
	   c1 = Kcy * py * (2.0 + py) + Kcx * px * ( 2.0 + px) + kineticTarget_;
	   b01 = -2.0 * Kcx * px * (1.0 + px);
	   a001 = Kcx * px * px;
	*/
	// scale all three dimensions, let c_x = c_y
	a000 = Kcx + Kcy;
	a110 = Kcz;
	c0 = kineticTarget_ - Kcx - Kcy - Kcz;
	a001 = Khx * px * px + Khy * py * py;
	a111 = Khz * pz * pz;
	b01 = -2.0 * (Khx * px * (1.0 + px) + Khy * py * (1.0 + py));
	b11 = -2.0 * Khz * pz * (1.0 + pz);
	c1 = Khx * px * (2.0 + px) + Khy * py * (2.0 + py)
	  + Khz * pz * (2.0 + pz) - kineticTarget_;
	break;
      case rnemdPx :
	c = 1 - momentumTarget_.x() / Pcx;
	a000 = Kcy;
	a110 = Kcz;
	c0 = Kcx * c * c - Kcx - Kcy - Kcz;
	a001 = py * py * Khy;
	a111 = pz * pz * Khz;
	b01 = -2.0 * Khy * py * (1.0 + py);
	b11 = -2.0 * Khz * pz * (1.0 + pz);
	c1 = Khy * py * (2.0 + py) + Khz * pz * (2.0 + pz)
	  + Khx * (fastpow(c * px - px - 1.0, 2) - 1.0);
	break;
      case rnemdPy :
	c = 1 - momentumTarget_.y() / Pcy;
	a000 = Kcx;
	a110 = Kcz;
	c0 = Kcy * c * c - Kcx - Kcy - Kcz;
	a001 = px * px * Khx;
	a111 = pz * pz * Khz;
	b01 = -2.0 * Khx * px * (1.0 + px);
	b11 = -2.0 * Khz * pz * (1.0 + pz);
	c1 = Khx * px * (2.0 + px) + Khz * pz * (2.0 + pz)
	  + Khy * (fastpow(c * py - py - 1.0, 2) - 1.0);
	break;
      case rnemdPz :// we don't really do this, do we?
	c = 1 - momentumTarget_.z() / Pcz;
	a000 = Kcx;
	a110 = Kcy;
	c0 = Kcz * c * c - Kcx - Kcy - Kcz;
	a001 = px * px * Khx;
	a111 = py * py * Khy;
	b01 = -2.0 * Khx * px * (1.0 + px);
	b11 = -2.0 * Khy * py * (1.0 + py);
	c1 = Khx * px * (2.0 + px) + Khy * py * (2.0 + py)
	  + Khz * (fastpow(c * pz - pz - 1.0, 2) - 1.0);
	break;
      default :
	break;
      }

      RealType v1 = a000 * a111 - a001 * a110;
      RealType v2 = a000 * b01;
      RealType v3 = a000 * b11;
      RealType v4 = a000 * c1 - a001 * c0;
      RealType v8 = a110 * b01;
      RealType v10 = - b01 * c0;

      RealType u0 = v2 * v10 - v4 * v4;
      RealType u1 = -2.0 * v3 * v4;
      RealType u2 = -v2 * v8 - v3 * v3 - 2.0 * v1 * v4;
      RealType u3 = -2.0 * v1 * v3;
      RealType u4 = - v1 * v1;
      // rescale coefficients
      RealType maxAbs = fabs(u0);
      if (maxAbs < fabs(u1)) maxAbs = fabs(u1);
      if (maxAbs < fabs(u2)) maxAbs = fabs(u2);
      if (maxAbs < fabs(u3)) maxAbs = fabs(u3);
      if (maxAbs < fabs(u4)) maxAbs = fabs(u4);
      u0 /= maxAbs;
      u1 /= maxAbs;
      u2 /= maxAbs;
      u3 /= maxAbs;
      u4 /= maxAbs;
      // max_element(start, end) is also available.
      Polynomial<RealType> poly; // same as DoublePolynomial poly;
      poly.setCoefficient(4, u4);
      poly.setCoefficient(3, u3);
      poly.setCoefficient(2, u2);
      poly.setCoefficient(1, u1);
      poly.setCoefficient(0, u0);
      vector<RealType> realRoots = poly.FindRealRoots();

      vector<RealType>::iterator ri;
      RealType r1, r2, alpha0;
      vector<pair<RealType,RealType> > rps;
      for (ri = realRoots.begin(); ri !=realRoots.end(); ++ri) {
	r2 = *ri;
	// Check to see if FindRealRoots() gave the right answer:
	if ( fabs(u0 + r2 * (u1 + r2 * (u2 + r2 * (u3 + r2 * u4)))) > 1e-6 ) {
	  sprintf(painCave.errMsg,
		  "RNEMD Warning: polynomial solve seems to have an error!");
	  painCave.isFatal = 0;
	  simError();
	  failRootCount_++;
	}
	// Might not be useful w/o rescaling coefficients
	alpha0 = -c0 - a110 * r2 * r2;
	if (alpha0 >= 0.0) {
	  r1 = sqrt(alpha0 / a000);
	  if (fabs(c1 + r1 * (b01 + r1 * a001) + r2 * (b11 + r2 * a111))
	      < 1e-6)
	    { rps.push_back(make_pair(r1, r2)); }
	  if (r1 > 1e-6) { // r1 non-negative
	    r1 = -r1;
	    if (fabs(c1 + r1 * (b01 + r1 * a001) + r2 * (b11 + r2 * a111))
		< 1e-6)
	      { rps.push_back(make_pair(r1, r2)); }
	  }
	}
      }
      // Consider combining together the part for solving for the pair
      // w/ the searching for the best solution part so that we don't
      // need the pairs vector:
      if (!rps.empty()) {
	RealType smallestDiff = HONKING_LARGE_VALUE;
	RealType diff(0.0);
	std::pair<RealType,RealType> bestPair = std::make_pair(1.0, 1.0);
	std::vector<std::pair<RealType,RealType> >::iterator rpi;
	for (rpi = rps.begin(); rpi != rps.end(); ++rpi) {
	  r1 = (*rpi).first;
	  r2 = (*rpi).second;
	  switch(rnemdFluxType_) {
	  case rnemdKE :
	    diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcx, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcy, 2);
	    break;
	  case rnemdPx :
	    diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcy, 2);
	    break;
	  case rnemdPy :
	    diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcz/Kcx, 2);
	    break;
	  case rnemdPz :
	    diff = fastpow(1.0 - r1, 2) + fastpow(1.0 - r2, 2)
	      + fastpow(r1 * r1 / r2 / r2 - Kcy/Kcx, 2);
	  default :
	    break;
	  }
	  if (diff < smallestDiff) {
	    smallestDiff = diff;
	    bestPair = *rpi;
	  }
	}
#ifdef IS_MPI
	if (worldRank == 0) {
#endif
	  // sprintf(painCave.errMsg,
	  //         "RNEMD: roots r1= %lf\tr2 = %lf\n",
	  //         bestPair.first, bestPair.second);
	  // painCave.isFatal = 0;
	  // painCave.severity = OPENMD_INFO;
	  // simError();
#ifdef IS_MPI
	}
#endif

	switch(rnemdFluxType_) {
	case rnemdKE :
	  x = bestPair.first;
	  y = bestPair.first;
	  z = bestPair.second;
	  break;
	case rnemdPx :
	  x = c;
	  y = bestPair.first;
	  z = bestPair.second;
	  break;
	case rnemdPy :
	  x = bestPair.first;
	  y = c;
	  z = bestPair.second;
	  break;
	case rnemdPz :
	  x = bestPair.first;
	  y = bestPair.second;
	  z = c;
	  break;
	default :
	  break;
	}
	vector<StuntDouble*>::iterator sdi;
	Vector3d vel;
	for (sdi = coldBin.begin(); sdi != coldBin.end(); ++sdi) {
	  vel = (*sdi)->getVel();
	  vel.x() *= x;
	  vel.y() *= y;
	  vel.z() *= z;
	  (*sdi)->setVel(vel);
	}
	// convert to hotBin coefficient
	x = 1.0 + px * (1.0 - x);
	y = 1.0 + py * (1.0 - y);
	z = 1.0 + pz * (1.0 - z);
	for (sdi = hotBin.begin(); sdi != hotBin.end(); ++sdi) {
	  vel = (*sdi)->getVel();
	  vel.x() *= x;
	  vel.y() *= y;
	  vel.z() *= z;
	  (*sdi)->setVel(vel);
	}
	successfulScale = true;
        switch(rnemdFluxType_) {
        case rnemdKE :
          kineticExchange_ += kineticTarget_;
	  break;
	case rnemdPx :
	case rnemdPy :
	case rnemdPz :
          momentumExchange_ += momentumTarget_;
	  break;
	default :
	  break;
	}
      }
    }
    if (successfulScale != true) {
      sprintf(painCave.errMsg,
              "RNEMD::doNIVS exchange NOT performed - roots that solve\n"
              "\tthe constraint equations may not exist or there may be\n"
              "\tno selected objects in one or both slabs.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      failTrialCount_++;
    }
  }

  void RNEMD::doVSS(SelectionManager& smanA, SelectionManager& smanB) {

    if (!doRNEMD_) return;
    int selei;
    int selej;

    StuntDouble* sd;

    vector<StuntDouble*> hotBin, coldBin;

    Vector3d Ph(V3Zero);
    Vector3d Lh(V3Zero);
    RealType Mh = 0.0;
    Mat3x3d Ih(0.0);
    RealType Kh = 0.0;
    std::size_t count_h = 0;
    Vector3d Pc(V3Zero);
    Vector3d Lc(V3Zero);
    RealType Mc = 0.0;
    Mat3x3d Ic(0.0);
    RealType Kc = 0.0;
    std::size_t count_c = 0;

    // Constraints can be on only the linear or angular momentum, but
    // not both.  Usually, the user will specify which they want, but
    // in case they don't, the use of periodic boundaries should make
    // the choice for us.
    bool doLinearPart = false;
    bool doAngularPart = false;

    switch (rnemdFluxType_) {
    case rnemdPx:
    case rnemdPy:
    case rnemdPz:
    case rnemdPvector:
    case rnemdKePx:
    case rnemdKePy:
    case rnemdKePvector:
      doLinearPart = true;
      break;
    case rnemdLx:
    case rnemdLy:
    case rnemdLz:
    case rnemdLvector:
    case rnemdKeLx:
    case rnemdKeLy:
    case rnemdKeLz:
    case rnemdKeLvector:
      doAngularPart = true;
      break;
    case rnemdKE:
    case rnemdRotKE:
    case rnemdFullKE:
    default:
      if (usePeriodicBoundaryConditions_)
        doLinearPart = true;
      else
        doAngularPart = true;
      break;
    }

    for (sd = smanA.beginSelected(selei); sd != NULL;
         sd = smanA.nextSelected(selei)) {

      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();
      Vector3d rPos = sd->getPos() - coordinateOrigin_;
      RealType r2;

      hotBin.push_back(sd);
      Ph += mass * vel;
      Mh += mass;
      Kh += mass * vel.lengthSquare();
      Lh += mass * cross(rPos, vel);
      Ih -= outProduct(rPos, rPos) * mass;
      r2 = rPos.lengthSquare();
      Ih(0, 0) += mass * r2;
      Ih(1, 1) += mass * r2;
      Ih(2, 2) += mass * r2;
      ++count_h;

      if (rnemdFluxType_ == rnemdFullKE) {
        if (sd->isDirectional()) {
          Vector3d angMom = sd->getJ();
          Mat3x3d I = sd->getI();
          if (sd->isLinear()) {
            int i = sd->linearAxis();
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;
            Kh += angMom[j] * angMom[j] / I(j, j) +
              angMom[k] * angMom[k] / I(k, k);
          } else {
            Kh += angMom[0] * angMom[0] / I(0, 0) +
              angMom[1] * angMom[1] / I(1, 1) +
              angMom[2] * angMom[2] / I(2, 2);
          }
        }
      }
    }
    for (sd = smanB.beginSelected(selej); sd != NULL;
         sd = smanB.nextSelected(selej)) {

      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();
      Vector3d rPos = sd->getPos() - coordinateOrigin_;
      RealType r2;

      coldBin.push_back(sd);
      Pc += mass * vel;
      Mc += mass;
      Kc += mass * vel.lengthSquare();
      Lc += mass * cross(rPos, vel);
      Ic -= outProduct(rPos, rPos) * mass;
      r2 = rPos.lengthSquare();
      Ic(0, 0) += mass * r2;
      Ic(1, 1) += mass * r2;
      Ic(2, 2) += mass * r2;
      ++count_c;

      if (rnemdFluxType_ == rnemdFullKE) {
        if (sd->isDirectional()) {
          Vector3d angMom = sd->getJ();
          Mat3x3d I = sd->getI();
          if (sd->isLinear()) {
            int i = sd->linearAxis();
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;
            Kc += angMom[j] * angMom[j] / I(j, j) +
              angMom[k] * angMom[k] / I(k, k);
          } else {
            Kc += angMom[0] * angMom[0] / I(0, 0) +
              angMom[1] * angMom[1] / I(1, 1) +
              angMom[2] * angMom[2] / I(2, 2);
          }
        }
      }
    }

    Kh *= 0.5;
    Kc *= 0.5;

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &Ph[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Pc[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Lh[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Lc[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Mh, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kh, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Mc, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kc, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, Ih.getArrayPointer(), 9,
                  MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, Ic.getArrayPointer(), 9,
                  MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &count_h, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &count_c, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    
#endif


    Vector3d ac, acrec, bc, bcrec;
    Vector3d ah, ahrec, bh, bhrec;

    bool successfulExchange = false;
    if ((Mh > 0.0) && (Mc > 0.0)) {	// both slabs are not empty
      Vector3d vc = Pc / Mc;
      ac = -momentumTarget_ / Mc + vc;
      acrec = -momentumTarget_ / Mc;

      // We now need the inverse of the inertia tensor to calculate the
      // angular velocity of the cold slab;
      Mat3x3d Ici = Ic.inverse();
      Vector3d omegac = Ici * Lc;
      bc  = -(Ici * angularMomentumTarget_) + omegac;
      bcrec = bc - omegac;

      RealType cNumerator = Kc - kineticTarget_;
      if (doLinearPart)
        cNumerator -= 0.5 * Mc * ac.lengthSquare();

      if (doAngularPart)
        cNumerator -= 0.5 * ( dot(bc, Ic * bc));

      if (cNumerator > 0.0) {

        RealType cDenominator = Kc;

        if (doLinearPart)
          cDenominator -= 0.5 * Mc * vc.lengthSquare();

        if (doAngularPart)
          cDenominator -= 0.5*(dot(omegac, Ic * omegac));

	if (cDenominator > 0.0) {
	  RealType c = sqrt(cNumerator / cDenominator);
	  if ((c > 0.9) && (c < 1.1)) {// restrict scaling coefficients

	    Vector3d vh = Ph / Mh;
            ah = momentumTarget_ / Mh + vh;
            ahrec = momentumTarget_ / Mh;

            // We now need the inverse of the inertia tensor to
            // calculate the angular velocity of the hot slab;
            Mat3x3d Ihi = Ih.inverse();
            Vector3d omegah = Ihi * Lh;
            bh  = (Ihi * angularMomentumTarget_) + omegah;
            bhrec = bh - omegah;

            RealType hNumerator = Kh + kineticTarget_;
            if (doLinearPart)
              hNumerator -= 0.5 * Mh * ah.lengthSquare();

            if (doAngularPart)
              hNumerator -= 0.5 * ( dot(bh, Ih * bh));

            if (hNumerator > 0.0) {

              RealType hDenominator = Kh;
              if (doLinearPart)
                hDenominator -= 0.5 * Mh * vh.lengthSquare();
              if (doAngularPart)
                hDenominator -= 0.5*(dot(omegah, Ih * omegah));

	      if (hDenominator > 0.0) {
		RealType h = sqrt(hNumerator / hDenominator);
		if ((h > 0.9) && (h < 1.1)) {

		  vector<StuntDouble*>::iterator sdi;
		  Vector3d vel;
                  Vector3d rPos;

		  for (sdi = coldBin.begin(); sdi != coldBin.end(); ++sdi) {
		    // vel = (*sdi)->getVel();
                    rPos = (*sdi)->getPos() - coordinateOrigin_;
                    if (doLinearPart)
                      vel = ((*sdi)->getVel() - vc) * c + ac;
                    if (doAngularPart)
                      vel = ((*sdi)->getVel() - cross(omegac, rPos)) * c + cross(bc, rPos);

		    (*sdi)->setVel(vel);
		    if (rnemdFluxType_ == rnemdFullKE) {
		      if ((*sdi)->isDirectional()) {
			Vector3d angMom = (*sdi)->getJ() * c;
			(*sdi)->setJ(angMom);
		      }
		    }
		  }
		  for (sdi = hotBin.begin(); sdi != hotBin.end(); ++sdi) {
		    // vel = (*sdi)->getVel();
                    rPos = (*sdi)->getPos() - coordinateOrigin_;
                    if (doLinearPart)
                      vel = ((*sdi)->getVel() - vh) * h + ah;
                    if (doAngularPart)
                      vel = ((*sdi)->getVel() - cross(omegah, rPos)) * h + cross(bh, rPos);

		    (*sdi)->setVel(vel);
		    if (rnemdFluxType_ == rnemdFullKE) {
		      if ((*sdi)->isDirectional()) {
			Vector3d angMom = (*sdi)->getJ() * h;
			(*sdi)->setJ(angMom);
		      }
		    }
		  }
		  successfulExchange = true;
		  kineticExchange_ += kineticTarget_;
                  momentumExchange_ += momentumTarget_;
                  particleFlux_h_ += smanA.getSelectionCount() / volumeA_ * momentumTarget_ / Mh;
                  particleFlux_c_ += smanB.getSelectionCount() / volumeB_ * momentumTarget_ / Mc;
                  angularMomentumExchange_ += angularMomentumTarget_;
		}
	      }
	    }
	  }
	}
      }
    }
    if (successfulExchange != true) {
      sprintf(painCave.errMsg,
              "RNEMD::doVSS exchange NOT performed - roots that solve\n"
              "\tthe constraint equations may not exist or there may be\n"
              "\tno selected objects in one or both slabs.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();
      failTrialCount_++;
    }
  }

  void RNEMD::doVSSCurrent(SelectionManager& smanA, SelectionManager& smanB) {

    if (!doRNEMD_) return;
    int selei;
    int selej;

    StuntDouble* sd;
    AtomType* atype;

    vector<StuntDouble*> hotBin, coldBin;

    Vector3d Ph(V3Zero);
    RealType Mh = 0.0;
    RealType Kh = 0.0;
    RealType Khz = 0.0;
    RealType MQhp = 0.0;
    RealType MQhn = 0.0;
    RealType MQvhp = 0.0;
    RealType MQvhn = 0.0;
    RealType MQ2hp = 0.0;
    RealType MQ2hn = 0.0;
    RealType Q2hp = 0.0;
    RealType Q2hn = 0.0;

    Vector3d Pc(V3Zero);
    RealType Mc = 0.0;
    RealType Kc = 0.0;
    RealType Kcz = 0.0;
    RealType MQcp = 0.0;
    RealType MQcn = 0.0;
    RealType MQvcp = 0.0;
    RealType MQvcn = 0.0;
    RealType MQ2cp = 0.0;
    RealType MQ2cn = 0.0;
    RealType Q2cp = 0.0;
    RealType Q2cn = 0.0;
    RealType Volc = 0.0;

    RealType Jc_total = 0.0;
    RealType Jc_cation = 0.0;
    RealType Jc_anion = 0.0;

    if (!usePeriodicBoundaryConditions_) {
      sprintf(painCave.errMsg,
              "RNEMD::doVSSCurrent can't be used without periodic boundary\n"
              "\tconditions enabled.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();
    }

    for (sd = smanA.beginSelected(selei); sd != NULL;
         sd = smanA.nextSelected(selei)) {

      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();

      hotBin.push_back(sd);
      Ph += mass * vel;
      Mh += mass;
      Kh += mass * vel.lengthSquare();
      Khz += mass * vel.z() * vel.z();

      RealType charge = 0.0;

      if (sd->isAtom()) {
        atype = static_cast<Atom*>(sd)->getAtomType();
        FixedChargeAdapter fca = FixedChargeAdapter(atype);
        if ( fca.isFixedCharge() ) charge = fca.getCharge();
        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
        if ( fqa.isFluctuatingCharge() ) charge += sd->getFlucQPos();
      }

      if (charge > 0) {
        MQhp += mass * charge;
        MQvhp += mass * vel.z() * charge;
        MQ2hp += mass * charge * charge;
        Q2hp += charge * charge;
      } else {
        MQhn += mass * charge;
        MQvhn += mass * vel.z() * charge;
        MQ2hn += mass * charge * charge;
        Q2hn += charge * charge;
      }
    }


    for (sd = smanB.beginSelected(selej); sd != NULL;
         sd = smanB.nextSelected(selej)) {

      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();

      coldBin.push_back(sd);
      Pc += mass * vel;
      Mc += mass;
      Kc += mass * vel.lengthSquare();
      Kcz += mass * vel.z() * vel.z();

      RealType charge = 0.0;

      if (sd->isAtom()) {
        atype = static_cast<Atom*>(sd)->getAtomType();
        FixedChargeAdapter fca = FixedChargeAdapter(atype);
        if ( fca.isFixedCharge() ) charge = fca.getCharge();
        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
        if ( fqa.isFluctuatingCharge() ) charge += sd->getFlucQPos();
      }

      if (charge > 0) {
        MQcp += mass * charge;
        MQvcp += mass * vel.z() * charge;
        MQ2cp += mass * charge * charge;
        Q2cp += charge * charge;
      } else {
        MQcn += mass * charge;
        MQvcn += mass * vel.z() * charge;
        MQ2cn += mass * charge * charge;
        Q2cn += charge * charge;
      }
    }

    Volc = volumeB_;

    Kh *= 0.5;
    Kc *= 0.5;
    Khz *= 0.5;
    Kcz *= 0.5;

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &Ph[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Pc[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &Mh, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kh, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Khz, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Mc, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kc, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Kcz, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &MQcp, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &MQcn, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &MQvcp, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &MQvcn, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &MQ2cp, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &MQ2cn, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Q2cp, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Q2cn, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &MQhp, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &MQhn, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &MQvhp, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &MQvhn, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &MQ2hp, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &MQ2hn, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Q2hp, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Q2hn, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif

    RealType alphac = 0.0;
    RealType betac = 0.0;
    RealType alphah = 0.0;
    RealType betah = 0.0;

    bool successfulExchange = false;
    if ( ( (MQcp > 0.0) && (MQcn < 0.0) && (MQhp > 0.0) && (MQhn < 0.0) ) || ( (MQcp > 0.0) && (MQhp > 0.0) ) || ( (MQcn < 0.0) && (MQhn < 0.0) ) ) {
      Vector3d vc = Pc / Mc;

      // units of velocity (Angstrom/fs):
      if (MQcn == 0.0) {
        alphac = currentDensity_ * Volc / Q2cp;
        alphah = -alphac * MQcp / MQhp;
        betac = 0.0;
        betah = 0.0;
      } else if (MQcp == 0.0) {
        alphac = 0.0;
        alphah = 0.0;
        betac = currentDensity_ * Volc / Q2cn;
        betah = -betac * MQcn / MQhn;
      } else {
        alphac = (currentDensity_ * Volc) / (Q2cp - Q2cn * (MQcp/MQcn));
        alphah = -alphac * MQcp / MQhp;
        betac = -alphac * MQcp / MQcn;
        betah = -betac * MQcn / MQhn;
      }

      Jc_total = (Q2cp * alphac + Q2cn * betac) / Volc;
      Jc_cation = Q2cp * alphac / Volc;
      Jc_anion = Q2cn * betac / Volc;

      RealType cNumerator = Kc - kineticTarget_;
      cNumerator -= 0.5 * Mc * vc.x()*vc.x();
      cNumerator -= 0.5 * Mc * vc.y()*vc.y();
      cNumerator -= Kcz;
      cNumerator -= MQvcp * alphac;
      cNumerator -= MQvcn * betac;
      cNumerator -= 0.5 * MQ2cp * alphac * alphac;
      cNumerator -= 0.5 * MQ2cn * betac * betac;

      RealType cDenominator = Kc - Kcz;
      cDenominator -= 0.5 * Mc * (vc.x()*vc.x() + vc.y()*vc.y());

      if (cNumerator/cDenominator > 0.0) {
        RealType c = sqrt(cNumerator / cDenominator);

        if ((c > 0.9) && (c < 1.1)) { // restrict scaling coefficients
          Vector3d vh = Ph / Mh;

          RealType hNumerator = Kh + kineticTarget_;
          hNumerator -= 0.5 * Mh * vh.x()*vh.x();
          hNumerator -= 0.5 * Mh * vh.y()*vh.y();
          hNumerator -= Khz;
          hNumerator -= MQvhp * alphah;
          hNumerator -= MQvhn * betah;
          hNumerator -= 0.5 * MQ2hp * alphah * alphah;
          hNumerator -= 0.5 * MQ2hn * betah * betah;

          RealType hDenominator = Kh - Khz;
          hDenominator -= 0.5 * Mh * (vh.x()*vh.x() + vh.y()*vh.y());

          if (hNumerator/hDenominator > 0.0) {
            RealType h = sqrt(hNumerator / hDenominator);

            if ((h > 0.9) && (h < 1.1)) {

              vector<StuntDouble*>::iterator sdi;
              Vector3d vel;
              Vector3d rPos;

              for (sdi = coldBin.begin(); sdi != coldBin.end(); ++sdi) {
                vel = (*sdi)->getVel();
                vel.x() = (vel.x() - vc.x()) * c + vc.x();
                vel.y() = (vel.y() - vc.y()) * c + vc.y();
                RealType q = 0.0;
                if ((*sdi)->isAtom()) {
                  atype = static_cast<Atom*>(*sdi)->getAtomType();
                  FixedChargeAdapter fca = FixedChargeAdapter(atype);
                  if ( fca.isFixedCharge() ) q = fca.getCharge();
                  FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
                  if ( fqa.isFluctuatingCharge() ) q += (*sdi)->getFlucQPos();
                  if (q > 0)
                    vel.z() += q * alphac;
                  else
                    vel.z() += q * betac;
                }
                (*sdi)->setVel(vel);
              }
              for (sdi = hotBin.begin(); sdi != hotBin.end(); ++sdi) {
                vel = (*sdi)->getVel();
                vel.x() = (vel.x() - vh.x()) * h + vh.x();
                vel.y() = (vel.y() - vh.y()) * h + vh.y();
                RealType q = 0.0;
                if ((*sdi)->isAtom()) {
                  atype = static_cast<Atom*>(*sdi)->getAtomType();
                  FixedChargeAdapter fca = FixedChargeAdapter(atype);
                  if ( fca.isFixedCharge() ) q = fca.getCharge();
                  FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
                  if ( fqa.isFluctuatingCharge() ) q += (*sdi)->getFlucQPos();
                  if (q > 0)
                    vel.z() += q * alphah;
                  else
                    vel.z() += q * betah;
                }
                (*sdi)->setVel(vel);
              }
              successfulExchange = true;
              kineticExchange_ += kineticTarget_;
            }
          }
        }
      }
    }
    if (successfulExchange != true) {
      Jc_total = 0.0;
      Jc_cation = 0.0;
      Jc_anion = 0.0;
      failTrialCount_++;
    }

    Jc_totalAccumulator_.add(Jc_total);
    Jc_cationAccumulator_.add(Jc_cation);
    Jc_anionAccumulator_.add(Jc_anion);
  }

  void RNEMD::doVSSSingle(SelectionManager& smanA) {

    if (!doRNEMD_) return;
    int selei;

    StuntDouble* sd;
    AtomType* atype;

    vector<StuntDouble*> hotBin;

    Vector3d P(V3Zero);
    RealType M = 0.0;
    RealType K = 0.0;
    RealType Mp = 0.0;
    RealType Mn = 0.0;
    RealType Qp = 0.0;
    RealType Qn = 0.0;
    RealType Pp = 0.0;
    RealType Pn = 0.0;
    RealType Vol = 0.0;

    RealType Jc_total = 0.0;
    RealType Jc_cation = 0.0;
    RealType Jc_anion = 0.0;

    if (!usePeriodicBoundaryConditions_) {
      sprintf(painCave.errMsg,
              "RNEMD::doVSSSingle can't be used without periodic boundary\n"
              "\tconditions enabled.\n");
      painCave.isFatal = 0;
      painCave.severity = OPENMD_INFO;
      simError();
    }

    for (sd = smanA.beginSelected(selei); sd != NULL;
         sd = smanA.nextSelected(selei)) {

      Vector3d pos = sd->getPos();

      // wrap the stuntdouble's position back into the box:

      if (usePeriodicBoundaryConditions_)
        currentSnap_->wrapVector(pos);

      RealType mass = sd->getMass();
      Vector3d vel = sd->getVel();

      hotBin.push_back(sd);
      P += mass * vel;
      M += mass;
      K += mass * vel.lengthSquare();

      RealType charge = 0.0;

      if (sd->isAtom()) {
        atype = static_cast<Atom*>(sd)->getAtomType();
        FixedChargeAdapter fca = FixedChargeAdapter(atype);
        if ( fca.isFixedCharge() ) charge = fca.getCharge();
        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
        if ( fqa.isFluctuatingCharge() ) charge += sd->getFlucQPos();

        if (charge > 0.0) {
          Mp += mass;
          Qp += charge;
          Pp += mass * vel.z();
        } else if (charge < 0.0) {
          Mn += mass;
          Qn += charge;
          Pn += mass * vel.z();
        }
      }
    }

    Vol = volumeA_;

    K *= 0.5;

#ifdef IS_MPI
    MPI_Allreduce(MPI_IN_PLACE, &P[0], 3, MPI_REALTYPE, MPI_SUM,
                  MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &M, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &K, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &Mp, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Mn, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Qp, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Qn, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Pp, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Pn, 1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif

    RealType alpha = 0.0;
    RealType beta = 0.0;

    bool successfulExchange = false;
    if ( (Mp > 0.0) && (Mn > 0.0) ) {
      Vector3d va = P / M;

      // units of velocity (Angstrom/fs):
      alpha = (currentDensity_ * Vol) / (Qp - Qn * (Mp/Mn));
      beta = -alpha * Mp / Mn;

      Jc_total = (Qp * alpha + Qn * beta) / Vol;
      Jc_cation = Qp * alpha / Vol;
      Jc_anion = Qn * beta / Vol;

      // Quadratic in a: A a^2 + B a + C = 0
      RealType A = K;
      A -= 0.5 * M * va.lengthSquare();

      RealType B = Pp * alpha;
      B += Pn * beta;

      RealType C = 0.5 * M * va.lengthSquare();
      C += 0.5 * Mp * alpha * alpha;
      C += 0.5 * Mn * beta * beta;
      C -= K;

      RealType insideSqrt = B * B - 4 * A * C;
      RealType a = 1.0;

      if (fabs(A) > std::numeric_limits<RealType>::epsilon()) {
        if (insideSqrt >= 0.0) {
          RealType root1 = (-B + sqrt(insideSqrt)) / (2 * A);
          RealType root2 = (-B - sqrt(insideSqrt)) / (2 * A);

          if ( (1.0 - root1) >= 0 ) {
            if ( ((1.0 - root2) >= 0) && ((1.0 - root1) > (1.0 - root2)) ) {
              a = root2;
            } else {
              a = root1;
            }
          } else if ( (1.0 - root2) >= 0 ) {
            a = root2;
          }
        } else
          a = 0.0;
      } else if (fabs(B) > std::numeric_limits<RealType>::epsilon()) {
        a = - C / B;
      } else a = 0.0;

      if ( (a > 0.9) && (a <= 1.0) ) { // restrict scaling coefficients
        vector<StuntDouble*>::iterator sdi;
        Vector3d vel;
        Vector3d rPos;

        for (sdi = hotBin.begin(); sdi != hotBin.end(); ++sdi) {
          vel = ((*sdi)->getVel() - va) * a + va;
          RealType q = 0.0;
          if ((*sdi)->isAtom()) {
            atype = static_cast<Atom*>(*sdi)->getAtomType();
            FixedChargeAdapter fca = FixedChargeAdapter(atype);
            if ( fca.isFixedCharge() )
              q = fca.getCharge();
            FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atype);
            if ( fqa.isFluctuatingCharge() )
              q += (*sdi)->getFlucQPos();
            if (q > 0.0)
              vel.z() += alpha;
            else if (q < 0.0)
              vel.z() += beta;
          }
          (*sdi)->setVel(vel);
        }
        successfulExchange = true;
      }
    }
    if (successfulExchange != true) {
      // sprintf(painCave.errMsg,
      //         "RNEMD::doVSSCurrent exchange NOT performed - roots that solve\n"
      //         "\tthe constraint equations may not exist or there may be\n"
      //         "\tno selected objects in one or both slabs.\n");
      // painCave.isFatal = 0;
      // painCave.severity = OPENMD_INFO;
      // simError();
      Jc_total = 0.0;
      Jc_cation = 0.0;
      Jc_anion = 0.0;
      failTrialCount_++;
    }

    Jc_totalAccumulator_.add(Jc_total);
    Jc_cationAccumulator_.add(Jc_cation);
    Jc_anionAccumulator_.add(Jc_anion);
  }

  int RNEMD::getBin(Vector3d pos) {

    if (usePeriodicBoundaryConditions_) {
      currentSnap_->wrapVector(pos);
      return int(nBins_ * (pos[rnemdPrivilegedAxis_] /
			   hmat_(rnemdPrivilegedAxis_, rnemdPrivilegedAxis_)
			   + 0.5)) % nBins_;
    } else {
      Vector3d rPos = pos - coordinateOrigin_;
      return int(rPos.length() / binWidth_);
    }
  }

  RealType RNEMD::getDividingArea() {

    if (hasDividingArea_) return dividingArea_;

    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();

    if (hasSelectionA_) {

      if (evaluatorA_.hasSurfaceArea()) {
        areaA_ = evaluatorA_.getSurfaceArea();
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
        areaA_ = surfaceMeshA->getArea();
        volumeA_ = surfaceMeshA->getVolume();
        delete surfaceMeshA;
#else
        sprintf( painCave.errMsg,
                 "RNEMD::getDividingArea : Hull calculation is not possible\n"
                 "\twithout libqhull. Please rebuild OpenMD with qhull enabled.");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
#endif
      }

    } else {
      if (usePeriodicBoundaryConditions_) {
        // in periodic boundaries, the surface area is twice the
        // area of the current box, normal to the privileged axis:
        switch(rnemdPrivilegedAxis_) {
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
        areaA_ = 4.0 * Constants::PI * std::pow(sphereARadius_, 2);
        volumeA_ = 4.0 * Constants::PI * std::pow(sphereARadius_, 3) / 3.0;
      }
    }

    if (hasSelectionB_) {
      if (evaluatorB_.hasSurfaceArea()) {
        areaB_ = evaluatorB_.getSurfaceArea();
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
        areaB_ = surfaceMeshB->getArea();
        volumeB_ = surfaceMeshB->getVolume();
        delete surfaceMeshB;
#else
        sprintf( painCave.errMsg,
                 "RNEMD::getDividingArea : Hull calculation is not possible\n"
                 "\twithout libqhull. Please rebuild OpenMD with qhull enabled.");
        painCave.severity = OPENMD_ERROR;
        painCave.isFatal = 1;
        simError();
#endif
      }

    } else {
      if (usePeriodicBoundaryConditions_) {
        // in periodic boundaries, the surface area is twice the
        // area of the current box, normal to the privileged axis:
        switch(rnemdPrivilegedAxis_) {
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

    dividingArea_ = min(areaA_, areaB_);
    hasDividingArea_ = true;
    return dividingArea_;
  }

  void RNEMD::doRNEMD() {

    if (!doRNEMD_) return;
    trialCount_++;

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

    // Target exchange quantities (in each exchange) = dividingArea * dt * flux
    // dt = exchange time interval
    // flux = target flux
    // dividingArea = smallest dividing surface between the two regions

    RealType area;

    if (info_->getSimParams()->getRNEMDParameters()->haveDividingArea()) {
      hasDividingArea_ = true;
      area = info_->getSimParams()->getRNEMDParameters()->getDividingArea();
    } else {
      hasDividingArea_ = false;
      area = getDividingArea();
    }

    kineticTarget_ = kineticFlux_ * exchangeTime_ * area;
    momentumTarget_ = momentumFluxVector_ * exchangeTime_ * area;
    angularMomentumTarget_ = angularMomentumFluxVector_ * exchangeTime_ * area;

    switch(rnemdMethod_) {
    case rnemdSwap:
      doSwap(commonA_, commonB_);
      break;
    case rnemdNIVS:
      doNIVS(commonA_, commonB_);
      break;
    case rnemdVSS:
      switch (rnemdFluxType_) {
      case rnemdCurrent:
      case rnemdKeCurrent:
        doVSSCurrent(commonA_, commonB_);
        break;
      case rnemdSingle:
        doVSSSingle(commonA_);
        break;
      default:
        doVSS(commonA_, commonB_);
        break;
      }
      break;
    case rnemdUnkownMethod:
    default :
      break;
    }
  }

  void RNEMD::collectData() {

    if (!doRNEMD_) return;
    currentSnap_ = info_->getSnapshotManager()->getCurrentSnapshot();
    hmat_ = currentSnap_->getHmat();

    // collectData can be called more frequently than the doRNEMD, so use the
    // computed area from the last exchange time:
    RealType area = getDividingArea();
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
    RealType charge = 0.0;

    vector<RealType> binMass(nBins_, 0.0);
    vector<Vector3d> binP(nBins_, V3Zero);
    vector<RealType> binOmega(nBins_, 0.0);
    vector<Vector3d> binL(nBins_, V3Zero);
    vector<Mat3x3d>  binI(nBins_);
    vector<RealType> binKE(nBins_, 0.0);
    vector<Vector3d> binEField(nBins_, V3Zero);
    vector<int> binDOF(nBins_, 0);
    vector<int> binCount(nBins_, 0);
    vector<int> binEFieldCount(nBins_, 0);
    vector<vector<int> > binTypeCounts;
    vector<vector<int> > binChargedTypeCounts;
    vector<vector<RealType> > binJc;

    if (outputMask_[ACTIVITY]) {
      binTypeCounts.resize(nBins_);
      for (unsigned int i = 0; i < nBins_; i++) {
        binTypeCounts[i].resize(outputTypes_.size(), 0);
      }
    }

    if (outputMask_[CURRENTDENSITY]) {
      binJc.resize(nBins_);
      binChargedTypeCounts.resize(nBins_);
      for (unsigned int i = 0; i < nBins_; i++) {
	binJc[i].resize(outputTypes_.size(), 0);
        binChargedTypeCounts[i].resize(outputTypes_.size(), 0);
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
	  binNo = getBin(pos);

	  mass = sd->getMass();
	  vel = sd->getVel();
	  rPos = sd->getPos() - coordinateOrigin_;
	  KE = 0.5 * mass * vel.lengthSquare();
	  L = mass * cross(rPos, vel);
	  I = outProduct(rPos, rPos) * mass;
	  r2 = rPos.lengthSquare();
	  I(0, 0) += mass * r2;
	  I(1, 1) += mass * r2;
	  I(2, 2) += mass * r2;

	  if (outputMask_[ACTIVITY] || outputMask_[CURRENTDENSITY]) {
	    typeIndex = -1;
	    if (sd->isAtom()) {
	      atype = static_cast<Atom*>(sd)->getAtomType();
	      at = std::find(outputTypes_.begin(), outputTypes_.end(), atype);
	      if (at != outputTypes_.end()) {
		typeIndex = std::distance(outputTypes_.begin(), at);
	      }
	    }
	  }

	  if (binNo >= 0 && binNo < int(nBins_))  {
	    binCount[binNo]++;
	    binMass[binNo] += mass;
	    binP[binNo] += mass*vel;
	    binKE[binNo] += KE;
	    binI[binNo] += I;
	    binL[binNo] += L;
	    binDOF[binNo] += 3;

	    if (outputMask_[ACTIVITY] && typeIndex != -1)
	      binTypeCounts[binNo][typeIndex]++;

	    if (outputMask_[CURRENTDENSITY] && typeIndex != -1) {
	      if (sd->isAtom()) {
		AtomType* atomType = static_cast<Atom*>(sd)->getAtomType();
		FixedChargeAdapter fca = FixedChargeAdapter(atomType);
		if (fca.isFixedCharge())
		  charge = fca.getCharge();

		FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
		if (fqa.isFluctuatingCharge())
		  charge += sd->getFlucQPos();
	      }

	      binJc[binNo][typeIndex] += charge * vel[rnemdPrivilegedAxis_];
	      binChargedTypeCounts[binNo][typeIndex]++;
	    }

	    if (sd->isDirectional()) {
	      Vector3d angMom = sd->getJ();
	      Mat3x3d Ia = sd->getI();
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

	// Calculate the electric field (kcal/mol/e/Angstrom) for all atoms in the box
	if (outputMask_[ELECTRICFIELD]) {
	  if (sd->isRigidBody()) {
	    RigidBody* rb = static_cast<RigidBody*>(sd);
	    std::vector<Atom*>::iterator ai;
	    Atom* atom;
	    for (atom = rb->beginAtom(ai); atom != NULL;
		 atom = rb->nextAtom(ai)) {

	      binNo = getBin(atom->getPos());
	      eField = atom->getElectricField();

	      binEFieldCount[binNo]++;
	      binEField[binNo] += eField;
	    }
	  } else {
	    eField = sd->getElectricField();
	    binNo = getBin(sd->getPos());

	    binEFieldCount[binNo]++;
	    binEField[binNo] += eField;
	  }
	}
      }
    }

#ifdef IS_MPI

    for (unsigned int i = 0; i < nBins_; i++) {

      MPI_Allreduce(MPI_IN_PLACE, &binCount[i],
                    1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &binMass[i],
                    1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, binP[i].getArrayPointer(),
                    3, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, binL[i].getArrayPointer(),
                    3, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, binI[i].getArrayPointer(),
                    9, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &binKE[i],
                    1, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &binDOF[i],
                    1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &binEFieldCount[i],
                    1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      if (outputMask_[ELECTRICFIELD]) {
        MPI_Allreduce(MPI_IN_PLACE, binEField[i].getArrayPointer(),
                      3, MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
      }
      if (outputMask_[ACTIVITY]) {
        MPI_Allreduce(MPI_IN_PLACE, &binTypeCounts[i][0],
                      outputTypes_.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      }
      if (outputMask_[CURRENTDENSITY]) {
        MPI_Allreduce(MPI_IN_PLACE, &binJc[i][0],
                      outputTypes_.size(), MPI_REALTYPE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &binChargedTypeCounts[i][0],
                      outputTypes_.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      }
    }

#endif

    Vector3d omega;
    RealType z, r, temp, binVolume, den(0.0), dz(0.0);
    std::vector<RealType> nden(outputTypes_.size(), 0.0);
    std::vector<RealType> Jc(outputTypes_.size(), 0.0);
    RealType boxVolume = currentSnap_->getVolume();
    RealType ePot(0.0);

    for (unsigned int i = 0; i < nBins_; i++) {

      if (usePeriodicBoundaryConditions_) {
        z = (((RealType)i + 0.5) / (RealType)nBins_) * hmat_(rnemdPrivilegedAxis_,rnemdPrivilegedAxis_);
        dynamic_cast<Accumulator *>(data_[Z].accumulator[i])->add(z);

        binVolume = boxVolume / nBins_;
        dz = hmat_(rnemdPrivilegedAxis_, rnemdPrivilegedAxis_) / (RealType)nBins_;
      } else {
        r = (((RealType)i + 0.5) * binWidth_);
        dynamic_cast<Accumulator *>(data_[R].accumulator[i])->add(r);

        RealType rinner = (RealType)i * binWidth_;
        RealType router = (RealType)(i+1) * binWidth_;
        binVolume = (4.0 * Constants::PI * (pow(router,3) - pow(rinner,3))) / 3.0;
      }

      // The calculations of the following properties are done regardless
      //   of whether or not the selected species are present in the bin
      if (outputMask_[ELECTRICFIELD] && binEFieldCount[i] > 0) {
	eField = binEField[i] / RealType(binEFieldCount[i]);
	dynamic_cast<VectorAccumulator *>(data_[ELECTRICFIELD].accumulator[i])->add(eField);
      }

      if (outputMask_[ELECTROSTATICPOTENTIAL]) {
	if (usePeriodicBoundaryConditions_ && binEFieldCount[i] > 0) {
	  ePot += eField[rnemdPrivilegedAxis_] * dz;
	  dynamic_cast<Accumulator *>(data_[ELECTROSTATICPOTENTIAL].accumulator[i])->add(ePot);
	}
      }

      // For the following properties, zero should be added if the selected
      //   species is not present in the bin
      if (outputMask_[DENSITY]) {
	den = binMass[i] * Constants::densityConvert / binVolume;
        dynamic_cast<Accumulator *>(data_[DENSITY].accumulator[i])->add(den);
      }

      if (outputMask_[ACTIVITY]) {
        for (unsigned int j = 0; j < outputTypes_.size(); j++) {
          nden[j] = (binTypeCounts[i][j] / binVolume)
            * Constants::concentrationConvert;
          dynamic_cast<Accumulator *>(data_[ACTIVITY].accumulatorArray2d[i][j])->add(nden[j]);
        }
      }

      if (binCount[i] > 0) {
	// The calculations of the following properties are meaningless if
	//   the selected species is not found in the bin
        if (outputMask_[VELOCITY]) {
          vel = binP[i] / binMass[i];
          dynamic_cast<VectorAccumulator *>(data_[VELOCITY].accumulator[i])->add(vel);
        }

        if (outputMask_[ANGULARVELOCITY]) {
          omega = binI[i].inverse() * binL[i];
          dynamic_cast<VectorAccumulator *>(data_[ANGULARVELOCITY].accumulator[i])->add(omega);
        }

        if (outputMask_[TEMPERATURE]) {
          temp = 2.0 * binKE[i] / (binDOF[i] * Constants::kb *
                                   Constants::energyConvert);
          dynamic_cast<Accumulator *>(data_[TEMPERATURE].accumulator[i])->add(temp);
        }

	if (outputMask_[CURRENTDENSITY]) {
	  for (unsigned int j = 0; j < outputTypes_.size(); j++) {
            if (binChargedTypeCounts[i][j] > 0) {
	      Jc[j] = binJc[i][j] / binVolume;
	      dynamic_cast<Accumulator *>(data_[CURRENTDENSITY].accumulatorArray2d[i][j])->add(Jc[j]);
            }
	  }
	}
      }
    }
    hasData_ = true;
  }

  void RNEMD::getStarted() {

    if (!doRNEMD_) return;
    if (info_->getSimParams()->getRNEMDParameters()->haveDividingArea())
      hasDividingArea_ = true;
    else
      hasDividingArea_ = false;
    collectData();
    writeOutputFile();
  }

  void RNEMD::parseOutputFileFormat(const std::string& format) {

    if (!doRNEMD_) return;
    StringTokenizer tokenizer(format, " ,;|\t\n\r");

    while(tokenizer.hasMoreTokens()) {
      std::string token(tokenizer.nextToken());
      toUpper(token);
      OutputMapType::iterator i = outputMap_.find(token);
      if (i != outputMap_.end()) {
        outputMask_.set(i->second);
      } else {
        sprintf( painCave.errMsg,
                 "RNEMD::parseOutputFileFormat: %s is not a recognized\n"
                 "\toutputFileFormat keyword.\n", token.c_str() );
        painCave.isFatal = 0;
        painCave.severity = OPENMD_ERROR;
        simError();
      }
    }
  }

  void RNEMD::writeOutputFile() {

    if (!doRNEMD_) return;
    if (!hasData_) return;

#ifdef IS_MPI
    // If we're the primary node, should we print out the results
    int worldRank;
    MPI_Comm_rank( MPI_COMM_WORLD, &worldRank);

    if (worldRank == 0) {
#endif
      rnemdFile_.open(rnemdFileName_.c_str(), std::ios::out | std::ios::trunc );

      if( !rnemdFile_ ) {
        sprintf( painCave.errMsg,
                 "Could not open \"%s\" for RNEMD output.\n",
                 rnemdFileName_.c_str());
        painCave.isFatal = 1;
        simError();
      }

      RealType time = currentSnap_->getTime();
      RealType avgArea = areaAccumulator_.getAverage();

      RealType avgJc_total = Jc_totalAccumulator_.getAverage();
      RealType avgJc_cation = Jc_cationAccumulator_.getAverage();
      RealType avgJc_anion = Jc_anionAccumulator_.getAverage();

      RealType Jz(0.0);
      Vector3d JzP(V3Zero);
      Vector3d JzL(V3Zero);

      if (time >= info_->getSimParams()->getDt()) {
        Jz = kineticExchange_ / (time * avgArea)
          / Constants::energyConvert;
        JzP = momentumExchange_ / (time * avgArea);
        JzL = angularMomentumExchange_ / (time * avgArea);
      }

      rnemdFile_ << "#######################################################\n";
      rnemdFile_ << "# RNEMD {\n";

      map<string, RNEMDMethod>::iterator mi;
      for(mi = stringToMethod_.begin(); mi != stringToMethod_.end(); ++mi) {
        if ( (*mi).second == rnemdMethod_)
          rnemdFile_ << "#    exchangeMethod  = \"" << (*mi).first << "\";\n";
      }
      map<string, RNEMDFluxType>::iterator fi;
      for(fi = stringToFluxType_.begin(); fi != stringToFluxType_.end(); ++fi) {
        if ( (*fi).second == rnemdFluxType_)
          rnemdFile_ << "#    fluxType  = \"" << (*fi).first << "\";\n";
      }
      if (usePeriodicBoundaryConditions_)
        rnemdFile_ << "#    privilegedAxis = " << rnemdAxisLabel_ << ";\n";
      rnemdFile_ << "#    exchangeTime = " << exchangeTime_ << ";\n";

      rnemdFile_ << "#    objectSelection = \""
                 << rnemdObjectSelection_ << "\";\n";
      rnemdFile_ << "#    selectionA = \"" << selectionA_ << "\";\n";
      rnemdFile_ << "#    selectionB = \"" << selectionB_ << "\";\n";
      rnemdFile_ << "#    outputSelection = \"" << outputSelection_ << "\";\n";
      rnemdFile_ << "# }\n";
      rnemdFile_ << "#######################################################\n";
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
      rnemdFile_ << "#   current density = " << currentDensity_
                 << " (electrons/A^2/fs)\n";
      rnemdFile_ << "# Target one-time exchanges:\n";
      rnemdFile_ << "#          kinetic = "
                 << kineticTarget_ / Constants::energyConvert
                 << " (kcal/mol)\n";
      rnemdFile_ << "#          momentum = " << momentumTarget_
                 << " (amu*A/fs)\n";
      rnemdFile_ << "#  angular momentum = " << angularMomentumTarget_
                 << " (amu*A^2/fs)\n";
      rnemdFile_ << "# Actual exchange totals:\n";
      rnemdFile_ << "#          kinetic = "
                 << kineticExchange_ / Constants::energyConvert
                 << " (kcal/mol)\n";
      rnemdFile_ << "#          momentum = " << momentumExchange_
                 << " (amu*A/fs)\n";
      if (rnemdFluxType_ == rnemdPvector && rnemdAxisLabel_ == "z") {
        rnemdFile_ << "#        part (hot) = " << particleFlux_h_
                 << " (particles/A^2/fs)\n";
        rnemdFile_ << "#        part (cold) = " << particleFlux_c_
                 << " (particles/A^2/fs)\n";
      }
      rnemdFile_ << "#  angular momentum = " << angularMomentumExchange_
                 << " (amu*A^2/fs)\n";
      rnemdFile_ << "# Actual flux:\n";
      rnemdFile_ << "#          kinetic = " << Jz
                 << " (kcal/mol/A^2/fs)\n";
      rnemdFile_ << "#          momentum = " << JzP
                 << " (amu/A/fs^2)\n";
      rnemdFile_ << "#  angular momentum = " << JzL
                 << " (amu/A^2/fs^2)\n";
      if ( (rnemdFluxType_ == rnemdCurrent)   ||
           (rnemdFluxType_ == rnemdKeCurrent) ||
           (rnemdFluxType_ == rnemdSingle) ) {
        rnemdFile_ << "#   Total current density = " << avgJc_total
                   << " (electrons/A^2/fs)\n";
        rnemdFile_ << "#   cation current density = " << avgJc_cation
                   << " (electrons/A^2/fs)\n";
        rnemdFile_ << "#   anion current density = " << avgJc_anion
                   << " (electrons/A^2/fs)\n";
      }
      rnemdFile_ << "# Exchange statistics:\n";
      rnemdFile_ << "#               attempted = " << trialCount_ << "\n";
      rnemdFile_ << "#                  failed = " << failTrialCount_ << "\n";
      if (rnemdMethod_ == rnemdNIVS) {
        rnemdFile_ << "#  NIVS root-check errors = "
                   << failRootCount_ << "\n";
      }
      rnemdFile_ << "#######################################################\n";

      // write title
      rnemdFile_ << "#";
      for (unsigned int i = 0; i < outputMask_.size(); ++i) {
        if (outputMask_[i]) {
          rnemdFile_ << "\t" << data_[i].title <<
            "(" << data_[i].units << ")";
          // add some extra tabs for column alignment
          if (data_[i].dataType == "Vector3d") rnemdFile_ << "\t\t";
          if (data_[i].dataType == "Array2d") {
            rnemdFile_ << "(";
            for (unsigned int j = 0;
                 j <  data_[i].accumulatorArray2d[0].size(); j++) {
              rnemdFile_<< outputTypes_[j]->getName() << "\t";
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
              writeReal(i,j);
            else if (data_[i].dataType == "Vector3d")
              writeVector(i,j);
            else if (data_[i].dataType == "Array2d")
              writeArray(i, j);
            else {
              sprintf( painCave.errMsg,
                       "RNEMD found an unknown data type for: %s ",
                       data_[i].title.c_str());
              painCave.isFatal = 1;
              simError();
            }
          }
        }
        rnemdFile_ << std::endl;

      }

      rnemdFile_ << "#######################################################\n";
      rnemdFile_ << "# 95% confidence intervals in those quantities follow:\n";
      rnemdFile_ << "#######################################################\n";


      for (unsigned int j = 0; j < nBins_; j++) {
        rnemdFile_ << "#";
        for (unsigned int i = 0; i < outputMask_.size(); ++i) {
          if (outputMask_[i]) {
            if (data_[i].dataType == "RealType")
              writeRealErrorBars(i,j);
            else if (data_[i].dataType == "Vector3d")
              writeVectorErrorBars(i,j);
            else if (data_[i].dataType == "Array2d")
              writeArrayErrorBars(i,j);
            else {
              sprintf( painCave.errMsg,
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

  void RNEMD::writeReal(int index, unsigned int bin) {

    if (!doRNEMD_) return;
    assert(index >= 0 && index < ENDINDEX);
    assert(bin < nBins_);

    RealType s;
    std::size_t count = data_[index].accumulator[bin]->count();

    if (count == 0) {
      rnemdFile_ << "\t";
    } else {
      dynamic_cast<Accumulator *>(data_[index].accumulator[bin])->getAverage(s);

      if ( std::isinf(s) || std::isnan(s) ) {
      	sprintf( painCave.errMsg,
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
      dynamic_cast<VectorAccumulator*>(data_[index].accumulator[bin])->getAverage(s);

      if ( std::isinf(s[0]) || std::isnan(s[0]) ||
	   std::isinf(s[1]) || std::isnan(s[1]) ||
	   std::isinf(s[2]) || std::isnan(s[2]) ) {
      	sprintf( painCave.errMsg,
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
        dynamic_cast<Accumulator*>(data_[index].accumulatorArray2d[bin][i])->getAverage(s);

        if ( std::isinf(s) || std::isnan(s) ) {
          sprintf( painCave.errMsg,
                   "RNEMD detected a numerical error writing: %s for bin %u, column %u",
                   data_[index].title.c_str(), bin, static_cast<unsigned int>(i));
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
      dynamic_cast<Accumulator *>(data_[index].accumulator[bin])->get95percentConfidenceInterval(s);

      if ( std::isinf(s) || std::isnan(s) ) {
      	sprintf( painCave.errMsg,
               	 "RNEMD detected a numerical error writing: %s std. dev. for bin %u",
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
      dynamic_cast<VectorAccumulator*>(data_[index].accumulator[bin])->get95percentConfidenceInterval(s);

      if ( std::isinf(s[0]) || std::isnan(s[0]) ||
	   std::isinf(s[1]) || std::isnan(s[1]) ||
	   std::isinf(s[2]) || std::isnan(s[2]) ) {
      	sprintf( painCave.errMsg,
               	 "RNEMD detected a numerical error writing: %s std. dev. for bin %u",
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
        dynamic_cast<Accumulator *>(data_[index].accumulatorArray2d[bin][i])->get95percentConfidenceInterval(s);

        if ( std::isinf(s) || std::isnan(s) ) {
          sprintf( painCave.errMsg,
                   "RNEMD detected a numerical error writing: %s std. dev. for bin %u, column %u",
                   data_[index].title.c_str(), bin, static_cast<unsigned int>(i));
          painCave.isFatal = 1;
          simError();
        } else {
          rnemdFile_ << "\t" << s;
        }
      }
    }
  }
}
