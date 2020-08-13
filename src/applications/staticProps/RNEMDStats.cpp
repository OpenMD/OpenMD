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
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "applications/staticProps/RNEMDStats.hpp"
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
#include "utils/Constants.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  RNEMDZ::RNEMDZ(SimInfo* info, const std::string& filename,
		 const std::string& sele, int nzbins, int axis)
    : SlabStatistics(info, filename, sele, nzbins, axis) {

    setOutputName(getPrefix(filename) + ".rnemdZ");

    evaluator_.loadScriptString(sele);
    seleMan_.setSelectionSet(evaluator_.evaluate());
    std::set<AtomType*> osTypes = seleMan_.getSelectedAtomTypes();
    std::copy(osTypes.begin(), osTypes.end(), std::back_inserter(outputTypes_));

    data_.resize(RNEMDZ::ENDINDEX);

    temperature = new OutputData;
    temperature->units =  "K";
    temperature->title =  "Temperature";
    temperature->dataType = odtReal;
    temperature->dataHandling = odhAverage;
    temperature->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      temperature->accumulator.push_back( new Accumulator() );
    addOutputDataAt(temperature, TEMPERATURE);

    velocity = new OutputData;
    velocity->units = "angstroms/fs";
    velocity->title =  "Velocity";
    velocity->dataType = odtVector3;
    velocity->dataHandling = odhAverage;
    velocity->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      velocity->accumulator.push_back( new VectorAccumulator() );
    addOutputDataAt(velocity, VELOCITY);

    density = new OutputData;
    density->units =  "g cm^-3";
    density->title =  "Density";
    density->dataType = odtReal;
    density->dataHandling = odhAverage;
    density->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      density->accumulator.push_back( new Accumulator() );
    addOutputDataAt(density, DENSITY);

    activity = new OutputData;;
    activity->units = "unitless";
    activity->title =  "Activity";
    activity->dataType = odtArray2d;
    activity->dataHandling = odhAverage;
    unsigned int nTypes = outputTypes_.size();
    // don't do activities if we don't have any atoms in the selection
    if (nTypes > 0) {
      activity->accumulatorArray2d.resize(nBins_);
      for (unsigned int i = 0; i < nBins_; i++) {
	activity->accumulatorArray2d[i].resize(nTypes);
	for (unsigned int j = 0 ; j < nTypes; j++) {
	  activity->accumulatorArray2d[i][j] = new Accumulator();
	}
      }
      for (unsigned int j = 0 ; j < nTypes; j++)
	activity->columnNames.push_back( outputTypes_[j]->getName() );
      addOutputDataAt(activity, ACTIVITY);
    }

    eField = new OutputData;
    eField->units =  "kcal/mol/angstroms/e";
    eField->title =  "Electric Field";
    eField->dataType = odtVector3;
    eField->dataHandling = odhAverage;
    eField->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      eField->accumulator.push_back( new VectorAccumulator() );

    ePot = new OutputData;
    ePot->units =  "kcal/mol/e";
    ePot->title =  "Electrostatic Potential";
    ePot->dataType = odtReal;
    ePot->dataHandling = odhAverage;
    ePot->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      ePot->accumulator.push_back( new Accumulator() );

    charge = new OutputData;
    charge->units =  "e";
    charge->title =  "Charge";
    charge->dataType = odtReal;
    charge->dataHandling = odhAverage;
    charge->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      charge->accumulator.push_back( new Accumulator() );

    chargeVelocity = new OutputData;
    chargeVelocity->units =  "e/fs";
    chargeVelocity->title =  "Charge_Velocity";
    chargeVelocity->dataType = odtReal;
    chargeVelocity->dataHandling = odhAverage;
    chargeVelocity->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      chargeVelocity->accumulator.push_back( new Accumulator() );

    outputMask_.set(TEMPERATURE);
    outputMask_.set(VELOCITY);
    outputMask_.set(DENSITY);
    outputMask_.set(ACTIVITY);

    int storageLayout = info_->getStorageLayout();

    if (storageLayout & DataStorage::dslElectricField) {
      outputMask_.set(ELECTRICFIELD);
      outputMask_.set(ELECTROSTATICPOTENTIAL);
      addOutputDataAt(eField, ELECTRICFIELD);
      addOutputDataAt(ePot, ELECTROSTATICPOTENTIAL);
    }

    if (info_->usesElectrostaticAtoms() ||
	storageLayout & DataStorage::dslFlucQPosition) {
      outputMask_.set(CHARGE);
      addOutputDataAt(charge, CHARGE);
    }

    if (storageLayout & DataStorage::dslFlucQVelocity) {
      outputMask_.set(CHARGEVELOCITY);
      addOutputDataAt(chargeVelocity, CHARGEVELOCITY);
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
    std::vector< std::vector<int> > binTypeCounts;
    std::vector<int> binEFieldCount(nBins_, 0);

    if (outputMask_[ACTIVITY]) {
      binTypeCounts.resize(nBins_);
      for (unsigned int i = 0; i < nBins_; i++) {
	binTypeCounts[i].resize( outputTypes_.size(), 0);
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
	  binNo = getBin(pos);

	  mass = sd->getMass();
	  vel = sd->getVel();
	  KE = 0.5 * mass * vel.lengthSquare();

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

	  if ( binNo >= 0 && binNo < int(nBins_) ) {
	    binCount[binNo]++;
	    binMass[binNo] += mass;
	    binP[binNo] += mass * vel;
	    binKE[binNo] += KE;
	    binDOF[binNo] += 3;

	    if ( outputMask_[ACTIVITY] && typeIndex != -1 )
	      binTypeCounts[binNo][typeIndex]++;

	    if ( outputMask_[CHARGE] || outputMask_[CHARGEVELOCITY] ) {
	      if (sd->isAtom()) {
		AtomType* atomType = static_cast<Atom*>(sd)->getAtomType();
		FixedChargeAdapter fca = FixedChargeAdapter(atomType);
		if ( fca.isFixedCharge() ) {
		  q = fca.getCharge();
		}
		FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
		if ( fqa.isFluctuatingCharge() ) {
		  q += sd->getFlucQPos();
		  w += sd->getFlucQVel();
		}

		if (outputMask_[CHARGE])
		  binCharge[binNo] += q;
		if (outputMask_[CHARGEVELOCITY])
		  binChargeVelocity[binNo] += w;
	      } else if (sd->isRigidBody()) {
		RigidBody* rb = static_cast<RigidBody*>(sd);
		std::vector<Atom*>::iterator ai;
		Atom* atom;
		for (atom = rb->beginAtom(ai); atom != NULL;
		     atom = rb->nextAtom(ai)) {

		  binNo = getBin(atom->getPos());
		  AtomType* atomType = atom->getAtomType();
		  FixedChargeAdapter fca = FixedChargeAdapter(atomType);
		  if ( fca.isFixedCharge() ) {
		    q = fca.getCharge();
		  }

		  FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
		  if ( fqa.isFluctuatingCharge() ) {
		    q += sd->getFlucQPos();
		    w += sd->getFlucQVel();
		  }

		  if (outputMask_[CHARGE])
		    binCharge[binNo] += q;
		  if (outputMask_[CHARGEVELOCITY])
		    binChargeVelocity[binNo] += w;
		}
	      }
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

    for (unsigned int i = 0; i < nBins_; i++) {

      RealType temp(0.0), ePot(0.0);
      Vector3d vel(0.0), eField(0.0);
      RealType z, den(0.0), binVolume(0.0), dz(0.0);
      std::vector<RealType> nden(outputTypes_.size(), 0.0);

      z = (((RealType)i + 0.5) / (RealType)nBins_) * hmat_(axis_, axis_);
      dynamic_cast<Accumulator *>(data_[Z]->accumulator[i])->add(z);

      binVolume = volume_ / nBins_;
      dz = hmat_(axis_, axis_) / (RealType)nBins_;

      // The calculations of the following properties are done regardless
      //   of whether or not the selected species are present in the bin
      if (outputMask_[ELECTRICFIELD] && binEFieldCount[i] > 0) {
	eField = binEField[i] / RealType(binEFieldCount[i]);
	dynamic_cast<VectorAccumulator *>(data_[ELECTRICFIELD]->accumulator[i])->add(eField);
      }

      if (outputMask_[ELECTROSTATICPOTENTIAL] && binEFieldCount[i] > 0) {
	ePot += eField[axis_] * dz;
	dynamic_cast<Accumulator *>(data_[ELECTROSTATICPOTENTIAL]->accumulator[i])->add(ePot);
      }

      // For the following properties, zero should be added if the selected
      //   species is not present in the bin
      if (outputMask_[DENSITY]) {
	den = binMass[i] * Constants::densityConvert / binVolume;
	dynamic_cast<Accumulator *>(data_[DENSITY]->accumulator[i])->add(den);
      }

      if (outputMask_[ACTIVITY]) {
	for (unsigned int j = 0; j < outputTypes_.size(); j++) {
	  nden[j] = (binTypeCounts[i][j] / binVolume)
	    * Constants::concentrationConvert;
	  dynamic_cast<Accumulator *>(data_[ACTIVITY]->accumulatorArray2d[i][j])->add(nden[j]);
	}
      }

      if (binCount[i] > 0) {
	// The calculations of the following properties are undefined if
	//   the selected species is not found in the bin
	if (outputMask_[VELOCITY]) {
	  vel = binP[i] / binMass[i];
	  dynamic_cast<VectorAccumulator *>(data_[VELOCITY]->accumulator[i])->add(vel);
	}

	if (outputMask_[TEMPERATURE]) {
	  temp = 2.0 * binKE[i] / (binDOF[i] * Constants::kb
				   * Constants::energyConvert);
	  dynamic_cast<Accumulator *>(data_[TEMPERATURE]->accumulator[i])->add(temp);
	}

	if (outputMask_[CHARGE])
	  dynamic_cast<Accumulator *>(data_[CHARGE]->accumulator[i])->add(binCharge[i]);
	if (outputMask_[CHARGEVELOCITY])
	  dynamic_cast<Accumulator *>(data_[CHARGEVELOCITY]->accumulator[i])->add(binChargeVelocity[i]);
      }
    }
  }


  RNEMDR::RNEMDR(SimInfo* info, const std::string& filename,
		 const std::string& sele, int nrbins)
    : ShellStatistics(info, filename, sele, nrbins) {

    setOutputName(getPrefix(filename) + ".rnemdR");

    data_.resize(RNEMDR::ENDINDEX);

    temperature = new OutputData;
    temperature->units =  "K";
    temperature->title =  "Temperature";
    temperature->dataType = odtReal;
    temperature->dataHandling = odhAverage;
    temperature->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      temperature->accumulator.push_back( new Accumulator() );
    addOutputDataAt(temperature, TEMPERATURE);

    angularVelocity = new OutputData;
    angularVelocity->units = "angstroms/fs";
    angularVelocity->title =  "Velocity";
    angularVelocity->dataType = odtVector3;
    angularVelocity->dataHandling = odhAverage;
    angularVelocity->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      angularVelocity->accumulator.push_back( new VectorAccumulator() );
    addOutputDataAt(angularVelocity, ANGULARVELOCITY);

    density = new OutputData;
    density->units =  "g cm^-3";
    density->title =  "Density";
    density->dataType = odtReal;
    density->dataHandling = odhAverage;
    density->accumulator.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++)
      density->accumulator.push_back( new Accumulator() );
    addOutputDataAt(density, DENSITY);
  }

  void RNEMDR::processFrame(int istep) {

    ShellStatistics::processFrame(istep);

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    StuntDouble* sd;
    int i;

    RealType mass;
    Vector3d vel;
    Vector3d rPos;
    RealType KE;
    Vector3d L;
    Mat3x3d I;
    RealType r2;

    std::vector<int> binDOF(nBins_, 0);
    std::vector<int> binCount(nBins_, 0);
    std::vector<RealType> binMass(nBins_, 0.0);
    std::vector<Vector3d> binP(nBins_, V3Zero);
    std::vector<RealType> binOmega(nBins_, 0.0);
    std::vector<Vector3d> binL(nBins_, V3Zero);
    std::vector<Mat3x3d>  binI(nBins_);
    std::vector<RealType> binKE(nBins_, 0.0);

    // loop over the selected atoms:
    for (sd = seleMan_.beginSelected(i); sd != NULL;
	 sd = seleMan_.nextSelected(i)) {

      // figure out where that object is:
      int binNo = getBin( sd->getPos() );

      if (binNo >= 0 && binNo < int(nBins_))  {

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

	binCount[binNo]++;
	binMass[binNo] += mass;
	binP[binNo] += mass*vel;
	binKE[binNo] += KE;
	binI[binNo] += I;
	binL[binNo] += L;
	binDOF[binNo] += 3;

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

    for (unsigned int i = 0; i < nBins_; i++) {

      RealType r, rinner, router, den(0.0), binVolume(0.0), temp(0.0);
      Vector3d omega(0.0);

      r = (((RealType)i + 0.5) * binWidth_);
      rinner = (RealType)i * binWidth_;
      router = (RealType)(i+1) * binWidth_;
      binVolume = (4.0 * Constants::PI * (pow(router,3) - pow(rinner,3))) / 3.0;

      dynamic_cast<Accumulator *>(data_[R]->accumulator[i])->add(r);

      // For the following properties, zero should be added if the selected
      //   species is not present in the bin
      den = binMass[i] * Constants::densityConvert / binVolume;
      dynamic_cast<Accumulator *>(data_[DENSITY]->accumulator[i])->add(den);

      if (binDOF[i] > 0) {
	// The calculations of the following properties are undefined if
	//   the selected species is not found in the bin
	omega = binI[i].inverse() * binL[i];
	dynamic_cast<VectorAccumulator *>(data_[ANGULARVELOCITY]->accumulator[i])->add(omega);

	temp = 2.0 * binKE[i] / (binDOF[i] * Constants::kb
				 * Constants::energyConvert);
	dynamic_cast<Accumulator *>(data_[TEMPERATURE]->accumulator[i])->add(temp);
      }
    }
  }


  RNEMDRTheta::RNEMDRTheta(SimInfo* info, const std::string& filename,
			   const std::string& sele, int nrbins, int nangleBins)
    : ShellStatistics(info, filename, sele, nrbins), nAngleBins_(nangleBins) {

    Globals* simParams = info->getSimParams();
    RNEMDParameters* rnemdParams = simParams->getRNEMDParameters();
    bool hasAngularMomentumFluxVector = rnemdParams->haveAngularMomentumFluxVector();

    if (hasAngularMomentumFluxVector) {
      std::vector<RealType> amf = rnemdParams->getAngularMomentumFluxVector();

      if (amf.size() != 3) {
	sprintf(painCave.errMsg,
		"RNEMDRTheta: Incorrect number of parameters specified for angularMomentumFluxVector.\n"
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

    data_.resize(RNEMDRTheta::ENDINDEX);

    angularVelocity = new OutputData;
    angularVelocity->units = "angstroms/fs";
    angularVelocity->title =  "Velocity";
    angularVelocity->dataType = odtVector3;
    angularVelocity->dataHandling = odhAverage;
    angularVelocity->accumulatorArray2d.reserve(nBins_);
    for (unsigned int i = 0; i < nBins_; i++) {
      angularVelocity->accumulatorArray2d[i].reserve(nAngleBins_);
      for (int j = 0 ; j < nAngleBins_; j++) {
	angularVelocity->accumulatorArray2d[i][j] = new Accumulator();
      }
    }
    addOutputDataAt(angularVelocity, ANGULARVELOCITY);
  }

  std::pair<int,int> RNEMDRTheta::getBins(Vector3d pos) {
    std::pair<int,int> result;

    Vector3d rPos = pos - coordinateOrigin_;
    RealType cosAngle= dot(rPos, fluxVector_) / rPos.length();

    result.first = int(rPos.length() / binWidth_);
    result.second = int( (nAngleBins_ - 1) * 0.5 * (cosAngle + 1.0) );
    return result;
  }

  void RNEMDRTheta::processFrame(int istep) {

    ShellStatistics::processFrame(istep);

    if (evaluator_.isDynamic()) {
      seleMan_.setSelectionSet(evaluator_.evaluate());
    }

    StuntDouble* sd;
    int i;

    std::vector< std::vector<int> > binCount;
    std::vector< std::vector<Mat3x3d> >  binI;
    std::vector< std::vector<Vector3d> > binL;


    // loop over the selected atoms:
    for (sd = seleMan_.beginSelected(i); sd != NULL;
	 sd = seleMan_.nextSelected(i)) {

      // figure out where that object is:
      std::pair<int,int> bins = getBins( sd->getPos() );

      if ( bins.first >= 0 && bins.first < int(nBins_) )  {
	if ( bins.second >= 0 && bins.second < nAngleBins_ ) {

	  Vector3d rPos = sd->getPos() - coordinateOrigin_;
	  Vector3d vel = sd->getVel();
	  RealType m = sd->getMass();
	  Vector3d L = m * cross(rPos, vel);
	  Mat3x3d I(0.0);
	  I = outProduct(rPos, rPos) * m;
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

      for (int j = 0; j < nAngleBins_; j++) {

	Vector3d omega(0.0);
	if (binCount[i][j] > 0) {
	  omega = binI[i][j].inverse() * binL[i][j];
	}

	RealType omegaProj = dot(omega, fluxVector_);

	dynamic_cast<Accumulator *>(data_[R]->accumulator[i])->add(r);
	dynamic_cast<Accumulator *>(data_[ANGULARVELOCITY]->accumulatorArray2d[i][j])->add(omegaProj);
      }
    }
  }

  void RNEMDRTheta::writeOutput() {

    std::vector<OutputData*>::iterator i;
    OutputData* outputData;

    std::ofstream outStream(outputFilename_.c_str());
    if (outStream.is_open()) {

      //write title
      outStream << "# SPATIAL STATISTICS\n";
      outStream << "#";

      for(outputData = beginOutputData(i); outputData;
	  outputData = nextOutputData(i)) {

	outStream << "\t" << outputData->title <<
	  "(" << outputData->units << ")";
	// add some extra tabs for column alignment
	if (outputData->dataType == odtVector3) outStream << "\t\t";
      }

      outStream << std::endl;

      outStream.precision(8);

      for (unsigned int j = 0; j < nBins_; j++) {

	for(outputData = beginOutputData(i); outputData;
	    outputData = nextOutputData(i)) {

	  int n = outputData->accumulator[j]->count();
	  if (n != 0) {
	    writeData( outStream, outputData, j );
	  }
	}
	outStream << std::endl;
      }
    }
  }
}
