/*
 * Copyright (c) 2005, 2009 The University of Notre Dame. All Rights Reserved.
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
 * @file Stats.cpp
 * @author tlin
 * @date 11/04/2004
 * @time 14:26am
 * @version 1.0
 */

#include "brains/Stats.hpp"
#include "brains/Thermo.hpp"

namespace OpenMD {

  Stats::Stats(SimInfo* info) : isInit_(false), info_(info) {   

    if (!isInit_) {
      init();
      isInit_ = true;
    }
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

    StatsData hydrogenbonding_potential;
    hydrogenbonding_potential.units =  "kcal/mol";
    hydrogenbonding_potential.title =  "Metallic Potential";    
    hydrogenbonding_potential.dataType = "RealType";
    hydrogenbonding_potential.accumulator = new Accumulator();
    data_[HYDROGENBONDING_POTENTIAL] = hydrogenbonding_potential;
    statsMap_["HYDROGENBONDING_POTENTIAL"] =  HYDROGENBONDING_POTENTIAL;

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

    StatsData pressure_tensor;
    pressure_tensor.units =  "amu*fs^-2*Ang^-1";
    pressure_tensor.title =  "Ptensor";
    pressure_tensor.dataType = "Mat3x3d";
    pressure_tensor.accumulator = new MatrixAccumulator();
    data_[PRESSURE_TENSOR] = pressure_tensor;
    statsMap_["PRESSURE_TENSOR"] =  PRESSURE_TENSOR;

    StatsData system_dipole;
    system_dipole.units =  "C*m";
    system_dipole.title =  "System Dipole";
    system_dipole.dataType = "Vector3d";
    system_dipole.accumulator = new VectorAccumulator();
    data_[SYSTEM_DIPOLE] = system_dipole;
    statsMap_["SYSTEM_DIPOLE"] =  SYSTEM_DIPOLE;

    StatsData tagged_pair_distance;
    tagged_pair_distance.units =  "Ang";
    tagged_pair_distance.title =  "Tagged_Pair_Distance";
    tagged_pair_distance.dataType = "RealType";
    tagged_pair_distance.accumulator = new Accumulator();
    data_[TAGGED_PAIR_DISTANCE] = tagged_pair_distance;
    statsMap_["TAGGED_PAIR_DISTANCE"] =  TAGGED_PAIR_DISTANCE;

    StatsData rnemd_exchange_total;
    rnemd_exchange_total.units =  "Variable";
    rnemd_exchange_total.title =  "RNEMD_exchange_total";
    rnemd_exchange_total.dataType = "RealType";
    rnemd_exchange_total.accumulator = new Accumulator();
    data_[RNEMD_EXCHANGE_TOTAL] = rnemd_exchange_total;
    statsMap_["RNEMD_EXCHANGE_TOTAL"] =  RNEMD_EXCHANGE_TOTAL;

    StatsData shadowh;
    shadowh.units =  "kcal/mol";
    shadowh.title =  "Shadow Hamiltonian";
    shadowh.dataType = "RealType";
    shadowh.accumulator = new Accumulator();
    data_[SHADOWH] = shadowh;
    statsMap_["SHADOWH"] =  SHADOWH;

    StatsData helfandmoment;
    helfandmoment.units =  "Ang*kcal/mol";
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

    // Why do we have both of these?
    if (simParams->getAccumulateBoxDipole()) {
      statsMask_.set(SYSTEM_DIPOLE);
    }
    if (info_->getCalcBoxDipole()){
      statsMask_.set(SYSTEM_DIPOLE);
    }

    if (simParams->havePrintHeatFlux()) {
      if (simParams->getPrintHeatFlux()){
        statsMask_.set(HEATFLUX);
      }
    }    
    
    
    if (simParams->haveTaggedAtomPair() && simParams->havePrintTaggedPairDistance()) {
      if (simParams->getPrintTaggedPairDistance()) {
        statsMask_.set(TAGGED_PAIR_DISTANCE);
      }
    }
    
    if (simParams->getRNEMDParameters()->getUseRNEMD()) {
      statsMask_.set(RNEMD_EXCHANGE_TOTAL);
    }



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
    Globals* simParams = info_->getSimParams();
    Snapshot* snap = info_->getSnapshotManager()->getCurrentSnapshot();
    Thermo thermo(info_);
   
    for (int i = 0; i < statsMask_.size(); ++i) {
      if (statsMask_[i]) {
        switch (i) {
        case TIME:
          data_[i].accumulator->add(snap->getTime());
          break;
        case KINETIC_ENERGY:
          data_[i].accumulator->add(thermo.getKinetic());
          break;
        case POTENTIAL_ENERGY:
          data_[i].accumulator->add(thermo.getPotential());
          break;
        case TOTAL_ENERGY:
          data_[i].accumulator->add(thermo.getTotalEnergy());
          break;
        case TEMPERATURE:
          data_[i].accumulator->add(thermo.getTemperature());
          break;
        case PRESSURE:
          data_[i].accumulator->add(thermo.getPressure());
          break;
        case VOLUME:
          data_[i].accumulator->add(thermo.getVolume());
          break;
        case CONSERVED_QUANTITY:
          data_[i].accumulator->add(snap->getConservedQuantity());
          break;
        case PRESSURE_TENSOR:
          dynamic_cast<MatrixAccumulator *>(data_[i].accumulator)->add(thermo.getPressureTensor());
          break;
        case SYSTEM_DIPOLE:
          dynamic_cast<VectorAccumulator *>(data_[i].accumulator)->add(thermo.getSystemDipole());
          break;
        case HEATFLUX:
          dynamic_cast<VectorAccumulator *>(data_[i].accumulator)->add(thermo.getHeatFlux());
          break;
        case HULLVOLUME:
          data_[i].accumulator->add(thermo.getHullVolume());
          break;
        case GYRVOLUME:
          data_[i].accumulator->add(thermo.getGyrationalVolume());
          break;
        case TRANSLATIONAL_KINETIC:
          data_[i].accumulator->add(thermo.getTranslationalKinetic());
          break;
        case ROTATIONAL_KINETIC:
          data_[i].accumulator->add(thermo.getRotationalKinetic());
          break;
        case LONG_RANGE_POTENTIAL:
          data_[i].accumulator->add(snap->getLongRangePotential());
          break;
        case VANDERWAALS_POTENTIAL:
          data_[i].accumulator->add(snap->getLongRangePotentials()[VANDERWAALS_FAMILY]);
          break;
        case ELECTROSTATIC_POTENTIAL:
          data_[i].accumulator->add(snap->getLongRangePotentials()[ELECTROSTATIC_FAMILY]);
          break;
        case METALLIC_POTENTIAL:
          data_[i].accumulator->add(snap->getLongRangePotentials()[METALLIC_FAMILY]);
          break;
        case HYDROGENBONDING_POTENTIAL:
          data_[i].accumulator->add(snap->getLongRangePotentials()[HYDROGENBONDING_FAMILY]);
          break;
        case SHORT_RANGE_POTENTIAL:
          data_[i].accumulator->add(snap->getShortRangePotential());
          break;
        case BOND_POTENTIAL:
          data_[i].accumulator->add(snap->getBondPotential());
          break;
        case BEND_POTENTIAL:
          data_[i].accumulator->add(snap->getBendPotential());
          break;
        case DIHEDRAL_POTENTIAL:
          data_[i].accumulator->add(snap->getTorsionPotential());
          break;
        case INVERSION_POTENTIAL:
          data_[i].accumulator->add(snap->getInversionPotential());
          break;
        case RAW_POTENTIAL:
          data_[i].accumulator->add(snap->getRawPotential());
          break;
        case RESTRAINT_POTENTIAL:
          data_[i].accumulator->add(snap->getRestraintPotential());
          break;
        case TAGGED_PAIR_DISTANCE:
          data_[i].accumulator->add(thermo.getTaggedAtomPairDistance());
          break;
          /*
        case RNEMD_EXCHANGE_TOTAL:
          data_[i].accumulator->add(thermo.get_RNEMD_exchange_total());
          break;
        case SHADOWH:
          data_[i].accumulator->add(thermo.getShadowHamiltionian());
          break;
        case HELFANDMOMENT:
          data_[i].accumulator->add(thermo.getHelfandMoment());
          break;
          */
        case ELECTRONIC_TEMPERATURE:
          data_[i].accumulator->add(thermo.getElectronicTemperature());
          break; 
        }
      }
    }
  }

  int Stats::getIntData(int index) { 
    assert(index >=0 && index < ENDINDEX);
    RealType value;
    data_[index].accumulator->getLastValue(value);
    return (int) value;
  }
  RealType Stats::getRealData(int index) {
    assert(index >=0 && index < ENDINDEX);
    RealType value(0.0);
    data_[index].accumulator->getLastValue(value);
    return value;
  }
  Vector3d Stats::getVectorData(int index) {
    assert(index >=0 && index < ENDINDEX);
    Vector3d value;
    dynamic_cast<VectorAccumulator*>(data_[index].accumulator)->getLastValue(value);
    return value;
  }
  Mat3x3d Stats::getMatrixData(int index) {
    assert(index >=0 && index < ENDINDEX);
    Mat3x3d value;
    dynamic_cast<MatrixAccumulator*>(data_[index].accumulator)->getLastValue(value);
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

}
