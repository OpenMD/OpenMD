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

#include "brains/Stats.hpp"
#include "brains/Thermo.hpp"
#include "utils/MemoryUtils.hpp"

#include <sstream>
#include <iomanip>

namespace OpenMD {

  Stats::Stats(SimInfo* info) : info_(info), isInit_(false) {

    if (!isInit_) {
      init();
      isInit_ = true;
    }
  }

  Stats::~Stats() {
    for (auto& data : data_) {
      if ( !data.accumulatorArray2d.empty() )
        MemoryUtils::deletePointers(data.accumulatorArray2d);
      else
        delete data.accumulator;
    }

    data_.clear();
    statsMap_.clear();
  }

  void Stats::init() {

    data_.resize(Stats::ENDINDEX);

    StatsData time;
    time.units =  "fs";
    time.title =  "Time";
    time.dataType = "RealType";
    time.accumulator = new Accumulator();
    data_[TIME] = time;
    statsMap_["TIME"] = TIME;

    StatsData total_energy;
    total_energy.units =  "kcal/mol";
    total_energy.title =  "Total Energy";
    total_energy.dataType = "RealType";
    total_energy.accumulator = new Accumulator();
    data_[TOTAL_ENERGY] = total_energy;
    statsMap_["TOTAL_ENERGY"] =  TOTAL_ENERGY;

    StatsData potential_energy;
    potential_energy.units =  "kcal/mol";
    potential_energy.title =  "Potential Energy";
    potential_energy.dataType = "RealType";
    potential_energy.accumulator = new Accumulator();
    data_[POTENTIAL_ENERGY] = potential_energy;
    statsMap_["POTENTIAL_ENERGY"] =  POTENTIAL_ENERGY;

    StatsData kinetic_energy;
    kinetic_energy.units =  "kcal/mol";
    kinetic_energy.title =  "Kinetic Energy";
    kinetic_energy.dataType = "RealType";
    kinetic_energy.accumulator = new Accumulator();
    data_[KINETIC_ENERGY] = kinetic_energy;
    statsMap_["KINETIC_ENERGY"] =  KINETIC_ENERGY;

    StatsData temperature;
    temperature.units =  "K";
    temperature.title =  "Temperature";
    temperature.dataType = "RealType";
    temperature.accumulator = new Accumulator();
    data_[TEMPERATURE] = temperature;
    statsMap_["TEMPERATURE"] =  TEMPERATURE;

    StatsData pressure;
    pressure.units =  "atm";
    pressure.title =  "Pressure";
    pressure.dataType = "RealType";
    pressure.accumulator = new Accumulator();
    data_[PRESSURE] = pressure;
    statsMap_["PRESSURE"] =  PRESSURE;


    StatsData volume;
    volume.units =  "A^3";
    volume.title =  "Volume";
    volume.dataType = "RealType";
    volume.accumulator = new Accumulator();
    data_[VOLUME] = volume;
    statsMap_["VOLUME"] =  VOLUME;

    StatsData hullvolume;
    hullvolume.units =  "A^3";
    hullvolume.title =  "Hull Volume";
    hullvolume.dataType = "RealType";
    hullvolume.accumulator = new Accumulator();
    data_[HULLVOLUME] = hullvolume;
    statsMap_["HULLVOLUME"] =  HULLVOLUME;

    StatsData gyrvolume;
    gyrvolume.units =  "A^3";
    gyrvolume.title =  "Gyrational Volume";
    gyrvolume.dataType = "RealType";
    gyrvolume.accumulator = new Accumulator();
    data_[GYRVOLUME] = gyrvolume;
    statsMap_["GYRVOLUME"] =  GYRVOLUME;

    StatsData conserved_quantity;
    conserved_quantity.units =  "kcal/mol";
    conserved_quantity.title =  "Conserved Quantity";
    conserved_quantity.dataType = "RealType";
    conserved_quantity.accumulator = new Accumulator();
    data_[CONSERVED_QUANTITY] = conserved_quantity;
    statsMap_["CONSERVED_QUANTITY"] =  CONSERVED_QUANTITY;

    StatsData translational_kinetic;
    translational_kinetic.units =  "kcal/mol";
    translational_kinetic.title =  "Translational Kinetic";
    translational_kinetic.dataType = "RealType";
    translational_kinetic.accumulator = new Accumulator();
    data_[TRANSLATIONAL_KINETIC] = translational_kinetic;
    statsMap_["TRANSLATIONAL_KINETIC"] =  TRANSLATIONAL_KINETIC;

    StatsData rotational_kinetic;
    rotational_kinetic.units =  "kcal/mol";
    rotational_kinetic.title =  "Rotational Kinetic";
    rotational_kinetic.dataType = "RealType";
    rotational_kinetic.accumulator = new Accumulator();
    data_[ROTATIONAL_KINETIC] = rotational_kinetic;
    statsMap_["ROTATIONAL_KINETIC"] =  ROTATIONAL_KINETIC;

    StatsData electronic_kinetic;
    electronic_kinetic.units =  "kcal/mol";
    electronic_kinetic.title =  "Electronic Kinetic";
    electronic_kinetic.dataType = "RealType";
    electronic_kinetic.accumulator = new Accumulator();
    data_[ELECTRONIC_KINETIC] = electronic_kinetic;
    statsMap_["ELECTRONIC_KINETIC"] =  ELECTRONIC_KINETIC;

    StatsData long_range_potential;
    long_range_potential.units =  "kcal/mol";
    long_range_potential.title =  "Long Range Potential";
    long_range_potential.dataType = "RealType";
    long_range_potential.accumulator = new Accumulator();
    data_[LONG_RANGE_POTENTIAL] = long_range_potential;
    statsMap_["LONG_RANGE_POTENTIAL"] =  LONG_RANGE_POTENTIAL;

    StatsData vanderwaals_potential;
    vanderwaals_potential.units =  "kcal/mol";
    vanderwaals_potential.title =  "van der waals Potential";
    vanderwaals_potential.dataType = "RealType";
    vanderwaals_potential.accumulator = new Accumulator();
    data_[VANDERWAALS_POTENTIAL] = vanderwaals_potential;
    statsMap_["VANDERWAALS_POTENTIAL"] =  VANDERWAALS_POTENTIAL;

    StatsData electrostatic_potential;
    electrostatic_potential.units =  "kcal/mol";
    electrostatic_potential.title =  "Electrostatic Potential";
    electrostatic_potential.dataType = "RealType";
    electrostatic_potential.accumulator = new Accumulator();
    data_[ELECTROSTATIC_POTENTIAL] = electrostatic_potential;
    statsMap_["ELECTROSTATIC_POTENTIAL"] =  ELECTROSTATIC_POTENTIAL;

    StatsData metallic_potential;
    metallic_potential.units =  "kcal/mol";
    metallic_potential.title =  "Metallic Potential";
    metallic_potential.dataType = "RealType";
    metallic_potential.accumulator = new Accumulator();
    data_[METALLIC_POTENTIAL] = metallic_potential;
    statsMap_["METALLIC_POTENTIAL"] =  METALLIC_POTENTIAL;

    StatsData metallic_embedding;
    metallic_embedding.units =  "kcal/mol";
    metallic_embedding.title =  "Metallic Embedding";
    metallic_embedding.dataType = "RealType";
    metallic_embedding.accumulator = new Accumulator();
    data_[METALLIC_EMBEDDING] = metallic_embedding;
    statsMap_["METALLIC_EMBEDDING"] =  METALLIC_EMBEDDING;

    StatsData metallic_pair;
    metallic_pair.units =  "kcal/mol";
    metallic_pair.title =  "Metallic Pair";
    metallic_pair.dataType = "RealType";
    metallic_pair.accumulator = new Accumulator();
    data_[METALLIC_PAIR] = metallic_pair;
    statsMap_["METALLIC_PAIR"] =  METALLIC_PAIR;

    StatsData hydrogenbonding_potential;
    hydrogenbonding_potential.units =  "kcal/mol";
    hydrogenbonding_potential.title =  "Hydrogen Bonding Pot.";
    hydrogenbonding_potential.dataType = "RealType";
    hydrogenbonding_potential.accumulator = new Accumulator();
    data_[HYDROGENBONDING_POTENTIAL] = hydrogenbonding_potential;
    statsMap_["HYDROGENBONDING_POTENTIAL"] =  HYDROGENBONDING_POTENTIAL;

    StatsData reciprocal_potential;
    reciprocal_potential.units =  "kcal/mol";
    reciprocal_potential.title =  "Reciprocal Space Pot.";
    reciprocal_potential.dataType = "RealType";
    reciprocal_potential.accumulator = new Accumulator();
    data_[RECIPROCAL_POTENTIAL] = reciprocal_potential;
    statsMap_["RECIPROCAL_POTENTIAL"] =  RECIPROCAL_POTENTIAL;

    StatsData surface_potential;
    surface_potential.units =  "kcal/mol";
    surface_potential.title =  "Surface Potential";
    surface_potential.dataType = "RealType";
    surface_potential.accumulator = new Accumulator();
    data_[SURFACE_POTENTIAL] = surface_potential;
    statsMap_["SURFACE_POTENTIAL"] =  SURFACE_POTENTIAL;

    StatsData short_range_potential;
    short_range_potential.units =  "kcal/mol";
    short_range_potential.title =  "Short Range Potential";
    short_range_potential.dataType = "RealType";
    short_range_potential.accumulator = new Accumulator();
    data_[SHORT_RANGE_POTENTIAL] = short_range_potential;
    statsMap_["SHORT_RANGE_POTENTIAL"] =  SHORT_RANGE_POTENTIAL;

    StatsData bond_potential;
    bond_potential.units =  "kcal/mol";
    bond_potential.title =  "Bond Potential";
    bond_potential.dataType = "RealType";
    bond_potential.accumulator = new Accumulator();
    data_[BOND_POTENTIAL] = bond_potential;
    statsMap_["BOND_POTENTIAL"] =  BOND_POTENTIAL;

    StatsData bend_potential;
    bend_potential.units =  "kcal/mol";
    bend_potential.title =  "Bend Potential";
    bend_potential.dataType = "RealType";
    bend_potential.accumulator = new Accumulator();
    data_[BEND_POTENTIAL] = bend_potential;
    statsMap_["BEND_POTENTIAL"] =  BEND_POTENTIAL;

    StatsData dihedral_potential;
    dihedral_potential.units =  "kcal/mol";
    dihedral_potential.title =  "Dihedral Potential";
    dihedral_potential.dataType = "RealType";
    dihedral_potential.accumulator = new Accumulator();
    data_[DIHEDRAL_POTENTIAL] = dihedral_potential;
    statsMap_["DIHEDRAL_POTENTIAL"] =  DIHEDRAL_POTENTIAL;

    StatsData inversion_potential;
    inversion_potential.units =  "kcal/mol";
    inversion_potential.title =  "Inversion Potential";
    inversion_potential.dataType = "RealType";
    inversion_potential.accumulator = new Accumulator();
    data_[INVERSION_POTENTIAL] = inversion_potential;
    statsMap_["INVERSION_POTENTIAL"] =  INVERSION_POTENTIAL;

    StatsData vraw;
    vraw.units =  "kcal/mol";
    vraw.title =  "Raw Potential";
    vraw.dataType = "RealType";
    vraw.accumulator = new Accumulator();
    data_[RAW_POTENTIAL] = vraw;
    statsMap_["RAW_POTENTIAL"] =  RAW_POTENTIAL;

    StatsData vrestraint;
    vrestraint.units =  "kcal/mol";
    vrestraint.title =  "Restraint Potential";
    vrestraint.dataType = "RealType";
    vrestraint.accumulator = new Accumulator();
    data_[RESTRAINT_POTENTIAL] = vrestraint;
    statsMap_["RESTRAINT_POTENTIAL"] =  RESTRAINT_POTENTIAL;

    StatsData vexcluded;
    vexcluded.units =  "kcal/mol";
    vexcluded.title =  "Excluded Potential";
    vexcluded.dataType = "RealType";
    vexcluded.accumulator = new Accumulator();
    data_[EXCLUDED_POTENTIAL] = vexcluded;
    statsMap_["EXCLUDED_POTENTIAL"] =  EXCLUDED_POTENTIAL;

    StatsData pressure_tensor;
    pressure_tensor.units =  "amu/fs^2/A";
    pressure_tensor.title =  "Pressure Tensor";
    pressure_tensor.dataType = "Mat3x3d";
    pressure_tensor.accumulator = new MatrixAccumulator();
    data_[PRESSURE_TENSOR] = pressure_tensor;
    statsMap_["PRESSURE_TENSOR"] =  PRESSURE_TENSOR;

    //virial tensor added
    StatsData virial_tensor;
    virial_tensor.units =  "kcal/mol";
    virial_tensor.title =  "Virial Tensor";
    virial_tensor.dataType = "Mat3x3d";
    virial_tensor.accumulator = new MatrixAccumulator();
    data_[VIRIAL_TENSOR] = virial_tensor;
    statsMap_["VIRIAL_TENSOR"] =  VIRIAL_TENSOR;

    StatsData system_dipole;
    system_dipole.units =  "C m";
    system_dipole.title =  "System Dipole";
    system_dipole.dataType = "Vector3d";
    system_dipole.accumulator = new VectorAccumulator();
    data_[SYSTEM_DIPOLE] = system_dipole;
    statsMap_["SYSTEM_DIPOLE"] =  SYSTEM_DIPOLE;

    StatsData system_quadrupole;
    system_quadrupole.units =  "C m^2";
    system_quadrupole.title =  "System Quadrupole";
    system_quadrupole.dataType = "Mat3x3d";
    system_quadrupole.accumulator = new MatrixAccumulator();
    data_[SYSTEM_QUADRUPOLE] = system_quadrupole;
    statsMap_["SYSTEM_QUADRUPOLE"] =  SYSTEM_QUADRUPOLE;

    StatsData tagged_pair_distance;
    tagged_pair_distance.units =  "A";
    tagged_pair_distance.title =  "Tagged Pair Distance";
    tagged_pair_distance.dataType = "RealType";
    tagged_pair_distance.accumulator = new Accumulator();
    data_[TAGGED_PAIR_DISTANCE] = tagged_pair_distance;
    statsMap_["TAGGED_PAIR_DISTANCE"] =  TAGGED_PAIR_DISTANCE;

    StatsData shadowh;
    shadowh.units =  "kcal/mol";
    shadowh.title =  "Shadow Hamiltonian";
    shadowh.dataType = "RealType";
    shadowh.accumulator = new Accumulator();
    data_[SHADOWH] = shadowh;
    statsMap_["SHADOWH"] =  SHADOWH;

    StatsData helfandmoment;
    helfandmoment.units =  "A*kcal/mol";
    helfandmoment.title =  "Thermal Helfand Moment";
    helfandmoment.dataType = "Vector3d";
    helfandmoment.accumulator = new VectorAccumulator();
    data_[HELFANDMOMENT] = helfandmoment;
    statsMap_["HELFANDMOMENT"] = HELFANDMOMENT;

    StatsData heatflux;
    heatflux.units = "amu/fs^3";
    heatflux.title =  "Heat Flux";
    heatflux.dataType = "Vector3d";
    heatflux.accumulator = new VectorAccumulator();
    data_[HEATFLUX] = heatflux;
    statsMap_["HEATFLUX"] = HEATFLUX;

    StatsData electronic_temperature;
    electronic_temperature.units = "K";
    electronic_temperature.title =  "Electronic Temperature";
    electronic_temperature.dataType = "RealType";
    electronic_temperature.accumulator = new Accumulator();
    data_[ELECTRONIC_TEMPERATURE] = electronic_temperature;
    statsMap_["ELECTRONIC_TEMPERATURE"] = ELECTRONIC_TEMPERATURE;

    StatsData com;
    com.units =  "A";
    com.title =  "Center of Mass";
    com.dataType = "Vector3d";
    com.accumulator = new VectorAccumulator();
    data_[COM] = com;
    statsMap_["COM"] =  COM;

    StatsData comVel;
    comVel.units =  "A/fs";
    comVel.title =  "COM Velocity";
    comVel.dataType = "Vector3d";
    comVel.accumulator = new VectorAccumulator();
    data_[COM_VELOCITY] = comVel;
    statsMap_["COM_VELOCITY"] =  COM_VELOCITY;

    StatsData angMom;
    angMom.units =  "amu A^2/fs";
    angMom.title =  "Angular Momentum";
    angMom.dataType = "Vector3d";
    angMom.accumulator = new VectorAccumulator();
    data_[ANGULAR_MOMENTUM] = angMom;
    statsMap_["ANGULAR_MOMENTUM"] =  ANGULAR_MOMENTUM;

    StatsData potSelection;
    potSelection.units =  "kcal/mol";
    potSelection.title =  "Selection Potentials";
    potSelection.dataType = "potVec";
    potSelection.accumulator = new PotVecAccumulator();
    data_[POTENTIAL_SELECTION] = potSelection;
    statsMap_["POTENTIAL_SELECTION"] =  POTENTIAL_SELECTION;

    StatsData netCharge;
    netCharge.units = "e";
    netCharge.title =  "Net Charge";
    netCharge.dataType = "RealType";
    netCharge.accumulator = new Accumulator();
    data_[NET_CHARGE] = netCharge;
    statsMap_["NET_CHARGE"] = NET_CHARGE;

    StatsData chargeMomentum;
    chargeMomentum.units = "kcal fs / e / mol";
    chargeMomentum.title =  "Charge Momentum";
    chargeMomentum.dataType = "RealType";
    chargeMomentum.accumulator = new Accumulator();
    data_[CHARGE_MOMENTUM] = chargeMomentum;
    statsMap_["CHARGE_MOMENTUM"] = CHARGE_MOMENTUM;

    StatsData currentDensity;
    currentDensity.units = "amps / m^2";
    currentDensity.title = "Current Density: total, then by AtomType";
    currentDensity.dataType = "Array2d";
    unsigned int nCols = 3 * (info_->getSimulatedAtomTypes().size()) + 3;
    currentDensity.accumulatorArray2d.resize(nCols);
    for (unsigned int j = 0 ; j < nCols; j++) {
      currentDensity.accumulatorArray2d[j] = new Accumulator();
    }
    data_[CURRENT_DENSITY] = currentDensity;
    statsMap_["CURRENT_DENSITY"] = CURRENT_DENSITY;

    // Now, set some defaults in the mask:

    Globals* simParams = info_->getSimParams();
    std::string statFileFormatString = simParams->getStatFileFormat();
    parseStatFileFormat(statFileFormatString);

    // if we're doing a thermodynamic integration, we'll want the raw
    // potential as well as the full potential:

    if (simParams->getUseThermodynamicIntegration())
      statsMask_.set(RAW_POTENTIAL);

    // if we've got restraints turned on, we'll also want a report of the
    // total harmonic restraints
    if (simParams->getUseRestraints()){
      statsMask_.set(RESTRAINT_POTENTIAL);
    }

    if (simParams->havePrintPressureTensor() &&
	simParams->getPrintPressureTensor()){
      statsMask_.set(PRESSURE_TENSOR);
    }
    
    if (simParams->havePrintVirialTensor() &&
        simParams->getPrintVirialTensor()){
      statsMask_.set(VIRIAL_TENSOR);
    }

    // Why do we have both of these?
    if (simParams->getAccumulateBoxDipole()) {
      statsMask_.set(SYSTEM_DIPOLE);
    }
    if (info_->getCalcBoxDipole()){
      statsMask_.set(SYSTEM_DIPOLE);
    }

    // Why do we have both of these?
    if (simParams->getAccumulateBoxQuadrupole()) {
      statsMask_.set(SYSTEM_QUADRUPOLE);
    }
    if (info_->getCalcBoxQuadrupole()){
      statsMask_.set(SYSTEM_QUADRUPOLE);
    }

    if (simParams->havePrintHeatFlux()) {
      if (simParams->getPrintHeatFlux()){
        statsMask_.set(HEATFLUX);
      }
    }

    if (simParams->haveTaggedAtomPair() &&
        simParams->havePrintTaggedPairDistance()) {
      if (simParams->getPrintTaggedPairDistance()) {
        statsMask_.set(TAGGED_PAIR_DISTANCE);
      }
    }

    if (simParams->havePotentialSelection()) {
      statsMask_.set(POTENTIAL_SELECTION);
    }
  }

  int Stats::getPrecision() {
    Globals* simParams = info_->getSimParams();
    int statFilePrecision = simParams->getStatFilePrecision();
    return statFilePrecision;
  }

  void Stats::parseStatFileFormat(const std::string& format) {
    StringTokenizer tokenizer(format, " ,;|\t\n\r");

    while(tokenizer.hasMoreTokens()) {
      std::string token(tokenizer.nextToken());
      toUpper(token);
      StatsMapType::iterator i = statsMap_.find(token);
      if (i != statsMap_.end()) {
        statsMask_.set(i->second);
      } else {
        sprintf( painCave.errMsg,
                 "Stats::parseStatFileFormat: %s is not a recognized\n"
                 "\tstatFileFormat keyword.\n", token.c_str() );
        painCave.isFatal = 0;
        painCave.severity = OPENMD_ERROR;
        simError();
      }
    }
  }

  std::string Stats::getTitle(int index) {
    assert(index >=0 && index < ENDINDEX);
    return data_[index].title;
  }

  std::string Stats::getUnits(int index) {
    assert(index >=0 && index < ENDINDEX);
    return data_[index].units;
  }

  std::string Stats::getDataType(int index) {
    assert(index >=0 && index < ENDINDEX);
    return data_[index].dataType;
  }

  void Stats::collectStats(){
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    Thermo thermo(info_);

    for (unsigned int i = 0; i < statsMask_.size(); ++i) {
      if (statsMask_[i]) {
        switch (i) {
        case TIME:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getTime());
          break;
        case KINETIC_ENERGY:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getKinetic());
          break;
        case POTENTIAL_ENERGY:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getPotential());
          break;
        case TOTAL_ENERGY:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getTotalEnergy());
          break;
        case TEMPERATURE:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getTemperature());
          break;
        case PRESSURE:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getPressure());
          break;
        case VOLUME:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getVolume());
          break;
        case CONSERVED_QUANTITY:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getConservedQuantity());
          break;
        case PRESSURE_TENSOR:
          dynamic_cast<MatrixAccumulator *>(data_[i].accumulator)->add(thermo.getPressureTensor());
          break;
          //virial Tensor
        case VIRIAL_TENSOR:
          dynamic_cast<MatrixAccumulator *>(data_[i].accumulator)->add(snap->getVirialTensor());
          break;

        case SYSTEM_DIPOLE:
          dynamic_cast<VectorAccumulator *>(data_[i].accumulator)->add(thermo.getSystemDipole());
          break;
        case SYSTEM_QUADRUPOLE:
          dynamic_cast<MatrixAccumulator *>(data_[i].accumulator)->add(thermo.getSystemQuadrupole());
          break;
        case HEATFLUX:
          dynamic_cast<VectorAccumulator *>(data_[i].accumulator)->add(thermo.getHeatFlux());
          break;
        case HULLVOLUME:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getHullVolume());
          break;
        case GYRVOLUME:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getGyrationalVolume());
          break;
        case TRANSLATIONAL_KINETIC:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getTranslationalKinetic());
          break;
        case ROTATIONAL_KINETIC:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getRotationalKinetic());
          break;
        case ELECTRONIC_KINETIC:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getElectronicKinetic());
          break;
        case LONG_RANGE_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getLongRangePotential());
          break;
        case VANDERWAALS_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getLongRangePotentials()[VANDERWAALS_FAMILY]);
          break;
        case ELECTROSTATIC_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getLongRangePotentials()[ELECTROSTATIC_FAMILY] +
                                                                 snap->getSelfPotentials()[ELECTROSTATIC_FAMILY]);
          break;
        case METALLIC_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getSelfPotentials()[METALLIC_EMBEDDING_FAMILY] +
                                                                 snap->getLongRangePotentials()[METALLIC_PAIR_FAMILY]);
          break;
        case METALLIC_EMBEDDING:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getSelfPotentials()[METALLIC_EMBEDDING_FAMILY]);
          break;
        case METALLIC_PAIR:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getLongRangePotentials()[METALLIC_PAIR_FAMILY]);
          break;
        case HYDROGENBONDING_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getLongRangePotentials()[HYDROGENBONDING_FAMILY]);
          break;
        case RECIPROCAL_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getReciprocalPotential());
          break;
        case SURFACE_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getSurfacePotential());
          break;
        case SHORT_RANGE_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getShortRangePotential());
          break;
        case BOND_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getBondPotential());
          break;
        case BEND_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getBendPotential());
          break;
        case DIHEDRAL_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getTorsionPotential());
          break;
        case INVERSION_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getInversionPotential());
          break;
        case RAW_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getRawPotential());
          break;
        case RESTRAINT_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getRestraintPotential());
          break;
        case EXCLUDED_POTENTIAL:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(snap->getExcludedPotential());
          break;
        case TAGGED_PAIR_DISTANCE:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getTaggedAtomPairDistance());
          break;
        case ELECTRONIC_TEMPERATURE:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getElectronicTemperature());
          break;
        case COM:
          dynamic_cast<VectorAccumulator *>(data_[i].accumulator)->add(thermo.getCom());
          break;
        case COM_VELOCITY:
          dynamic_cast<VectorAccumulator *>(data_[i].accumulator)->add(thermo.getComVel());
          break;
        case ANGULAR_MOMENTUM:
          dynamic_cast<VectorAccumulator *>(data_[i].accumulator)->add(thermo.getAngularMomentum());
          break;
        case POTENTIAL_SELECTION:
          dynamic_cast<PotVecAccumulator *>(data_[i].accumulator)->add(thermo.getSelectionPotentials());
          break;
        case NET_CHARGE:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getNetCharge());
          break;
        case CHARGE_MOMENTUM:
          dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getChargeMomentum());
          break;
        case CURRENT_DENSITY:
          unsigned int k = 0;
          std::vector<Vector3d> Jc = thermo.getCurrentDensity();
          for (unsigned int j = 0; j < Jc.size(); ++j) {
            dynamic_cast<Accumulator *>(data_[i].accumulatorArray2d[k++])->add(Jc[j][0]);
            dynamic_cast<Accumulator *>(data_[i].accumulatorArray2d[k++])->add(Jc[j][1]);
            dynamic_cast<Accumulator *>(data_[i].accumulatorArray2d[k++])->add(Jc[j][2]);
          }
          break;

          /*
            case SHADOWH:
            dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getShadowHamiltionian());
            break;
            case HELFANDMOMENT:
            dynamic_cast<Accumulator *>(data_[i].accumulator)->add(thermo.getHelfandMoment());
            break;
          */
        }
      }
    }
  }

  int Stats::getIntData(int index) {
    assert(index >=0 && index < ENDINDEX);
    RealType value;
    dynamic_cast<Accumulator *>(data_[index].accumulator)->getLastValue(value);
    return (int) value;
  }
  RealType Stats::getRealData(int index) {
    assert(index >=0 && index < ENDINDEX);
    RealType value(0.0);
    dynamic_cast<Accumulator *>(data_[index].accumulator)->getLastValue(value);
    return value;
  }
  Vector3d Stats::getVectorData(int index) {
    assert(index >=0 && index < ENDINDEX);
    Vector3d value;
    dynamic_cast<VectorAccumulator*>(data_[index].accumulator)->getLastValue(value);
    return value;
  }
  potVec Stats::getPotVecData(int index) {
    assert(index >=0 && index < ENDINDEX);
    potVec value;
    dynamic_cast<PotVecAccumulator*>(data_[index].accumulator)->getLastValue(value);
    return value;
  }
  Mat3x3d Stats::getMatrixData(int index) {
    assert(index >=0 && index < ENDINDEX);
    Mat3x3d value;
    dynamic_cast<MatrixAccumulator*>(data_[index].accumulator)->getLastValue(value);
    return value;
  }
  std::vector<RealType> Stats::getArrayData(int index) {
    assert(index >=0 && index < ENDINDEX);
    std::vector<RealType> value;
    RealType v;
    for (unsigned int i = 0; i < data_[index].accumulatorArray2d.size(); ++i) {
      dynamic_cast<Accumulator*>(data_[index].accumulatorArray2d[i])->getLastValue(v);
      value.push_back(v);
    }
    return value;
  }

  int Stats::getIntAverage(int index) {
    assert(index >=0 && index < ENDINDEX);
    RealType value;
    dynamic_cast<Accumulator *>(data_[index].accumulator)->getAverage(value);
    return (int) value;
  }
  RealType Stats::getRealAverage(int index) {
    assert(index >=0 && index < ENDINDEX);
    RealType value(0.0);
    dynamic_cast<Accumulator *>(data_[index].accumulator)->getAverage(value);
    return value;
  }
  Vector3d Stats::getVectorAverage(int index) {
    assert(index >=0 && index < ENDINDEX);
    Vector3d value;
    dynamic_cast<VectorAccumulator*>(data_[index].accumulator)->getAverage(value);
    return value;
  }
  potVec Stats::getPotVecAverage(int index) {
    assert(index >=0 && index < ENDINDEX);
    potVec value;
    dynamic_cast<PotVecAccumulator*>(data_[index].accumulator)->getAverage(value);
    return value;
  }
  Mat3x3d Stats::getMatrixAverage(int index) {
    assert(index >=0 && index < ENDINDEX);
    Mat3x3d value;
    dynamic_cast<MatrixAccumulator*>(data_[index].accumulator)->getAverage(value);
    return value;
  }
  std::vector<RealType> Stats::getArrayAverage(int index) {
    assert(index >=0 && index < ENDINDEX);
    std::vector<RealType> value;
    RealType v;
    for (unsigned int i = 0; i < data_[index].accumulatorArray2d.size(); ++i) {
      dynamic_cast<Accumulator*>(data_[index].accumulatorArray2d[i])->getAverage(v);
      value.push_back(v);
    }
    return value;
  }


  int Stats::getIntError(int index) {
    assert(index >=0 && index < ENDINDEX);
    RealType value;
    dynamic_cast<Accumulator *>(data_[index].accumulator)->get95percentConfidenceInterval(value);
    return (int) value;
  }
  RealType Stats::getRealError(int index) {
    assert(index >=0 && index < ENDINDEX);
    RealType value(0.0);
    dynamic_cast<Accumulator *>(data_[index].accumulator)->get95percentConfidenceInterval(value);
    return value;
  }
  Vector3d Stats::getVectorError(int index) {
    assert(index >=0 && index < ENDINDEX);
    Vector3d value;
    dynamic_cast<VectorAccumulator*>(data_[index].accumulator)->get95percentConfidenceInterval(value);
    return value;
  }
  potVec Stats::getPotVecError(int index) {
    assert(index >=0 && index < ENDINDEX);
    potVec value;
    dynamic_cast<PotVecAccumulator*>(data_[index].accumulator)->get95percentConfidenceInterval(value);
    return value;
  }
  Mat3x3d Stats::getMatrixError(int index) {
    assert(index >=0 && index < ENDINDEX);
    Mat3x3d value;
    dynamic_cast<MatrixAccumulator*>(data_[index].accumulator)->get95percentConfidenceInterval(value);
    return value;
  }
  std::vector<RealType> Stats::getArrayError(int index) {
    assert(index >=0 && index < ENDINDEX);
    std::vector<RealType> value;
    RealType v;
    for (unsigned int i = 0; i < data_[index].accumulatorArray2d.size(); ++i) {
      dynamic_cast<Accumulator*>(data_[index].accumulatorArray2d[i])->get95percentConfidenceInterval(v);
      value.push_back(v);
    }
    return value;
  }


  Stats::StatsBitSet Stats::getStatsMask() {
    return statsMask_;
  }
  Stats::StatsMapType Stats::getStatsMap() {
    return statsMap_;
  }
  void Stats::setStatsMask(Stats::StatsBitSet mask) {
    statsMask_ = mask;
  }

  std::string Stats::getStatsReport() {
    std::stringstream report;
#if defined (_MSC_VER)
    std::string pm = " +/- ";
    std::string luc = "[";
    std::string lex = "|";
    std::string llc = "[";
    std::string ruc = "]";
    std::string rex = "|";
    std::string rlc = "]";
#else
    std::string pm = "  \u00B1  ";
    std::string luc = "\u23A1";
    std::string lex = "\u23A2";
    std::string llc = "\u23A3";
    std::string ruc = "\u23A4";
    std::string rex = "\u23A5";
    std::string rlc = "\u23A6";
#endif

    int nSamp = dynamic_cast<Accumulator *>(data_[TIME].accumulator)->count();

    std::string head(79, '#');
    report << head << std::endl;
    report << "# Status Report:" << std::string(62, ' ') << "#" << std::endl;
    report << "# " << right << setw(24) << "Total Time:";
    report << setw(12) << getRealData(TIME);
    report << " " << setw(17) << left << getUnits(TIME)
           << "                      #" << std::endl;
    report << "# " << right << setw(24) << "Number of Samples:";
    report << setw(12) << nSamp;
    report << "                                        #" << std::endl;

    for (unsigned int i = 0; i < statsMask_.size(); ++i) {
      if (statsMask_[i] && i != TIME) {

        if (getDataType(i) == "RealType") {
          report << "# " << right << setw(23) << getTitle(i) << ":";
          report << right << setw(12) << getRealAverage(i);
          report << pm << left << setw(12) << getRealError(i);
          report << " " << left << setw(17) << getUnits(i);
          report << "     #" << std::endl;

        }
        else if (getDataType(i) == "Vector3d") {
          Vector3d s = getVectorAverage(i);
          Vector3d e = getVectorError(i);

          report << "#                       ";
          report << luc << right << setw(12) << s(0) << ruc << "     ";
          report << luc << right <<  setw(12) << e(0);
          report << ruc << "                    #" << std::endl;

          report << "# " << right << setw(23) << getTitle(i) << ":";

          report << lex << right << setw(12) << s(1) << rex << pm;
          report << lex << right << setw(12) << e(1) << rex << " ";
          report << left << setw(17) << getUnits(i) << "  #";
          report << std::endl;

          report << "#                       ";
          report << llc << right << setw(12) << s(2) << rlc << "     ";
          report << llc << right <<  setw(12) << e(2);
          report << rlc << "                     #" << std::endl;

        }
        else if (getDataType(i) == "potVec") {
          potVec s = getPotVecAverage(i);
          potVec e = getPotVecError(i);

          report << "# " << right << setw(23) << getTitle(i);
          report << ":                                                    #";
          report << std::endl;

          for (unsigned int j = 1; j < N_INTERACTION_FAMILIES; j++) {
            switch (j) {
            case VANDERWAALS_FAMILY:
              report << "# " << right << setw(24) << "van der Waals:";
              break;
            case ELECTROSTATIC_FAMILY:
              report << "# " << right << setw(24) << "Electrostatic:";
              break;
            case METALLIC_EMBEDDING:
              report << "# " << right << setw(24) << "Metallic Embedding:";
              break;
            case METALLIC_PAIR:
              report << "# " << right << setw(24) << "Metallic Pair:";
              break;
            case HYDROGENBONDING_FAMILY:
              report << "# " << right << setw(24) << "Hydrogen Bonding:";
              break;
            case BONDED_FAMILY:
              report << "# " << right << setw(24) << "Bonded (1-2,1-3,1-4):";
              break;
            default:
              report << "# " << right << setw(24) << "Unknown:";
              break;
            }
            report << right << setw(12) << s[j];
            report << pm << left << setw(12) << e[j];
            report << " " << left << setw(17) << getUnits(i) << "     #";
            report << std::endl;
          }

        }
        else if (getDataType(i) == "Mat3x3d") {
          Mat3x3d s = getMatrixAverage(i);
          Mat3x3d e = getMatrixError(i);

          report << "#                         ";
          report << luc << right << setw(12) << s(0,0) << " ";
          report << right << setw(12) << s(0,1) << " ";
          report << right << setw(12) << s(0,2) << ruc << "            #";
          report << std::endl;

          report << "# " << right << setw(23) << getTitle(i) << ":";

          report << lex << right << setw(12) << s(1,0) << " ";
          report << right << setw(12) << s(1,1) << " ";
          report << right << setw(12) << s(1,2) << rex << " ";
          report << left <<  setw(11) << getUnits(i) << "#" << std::endl;

          report << "#                         ";
          report << llc << right << setw(12) << s(2,0) << " ";
          report << right << setw(12) << s(2,1) << " ";
          report << right << setw(12) << s(2,2) << rlc << "            #";
          report << std::endl;

          report << "#                                 ";
          report << luc << right << setw(12) << e(0,0) << " ";
          report << right << setw(12) << e(0,1) << " ";
          report << right << setw(12) << e(0,2) << ruc << "    #";
          report << std::endl;

          report << "#                            " << pm ;

          report << lex << right << setw(12) << e(1,0) << " ";
          report << right << setw(12) << e(1,1) << " ";
          report << right << setw(12) << e(1,2) << rex << "    #";
          report << std::endl;

          report << "#                                 ";
          report << llc << right << setw(12) << e(2,0) << " ";
          report << right << setw(12) << e(2,1) << " ";
          report << right << setw(12) << e(2,2) << rlc << "    #";
          report << std::endl;

        }
      }
    }
    report << head << std::endl;

    return report.str();
  }

}
