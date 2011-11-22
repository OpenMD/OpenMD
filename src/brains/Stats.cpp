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

namespace OpenMD {

  bool Stats::isInit_ = false;
  std::string Stats::title_[Stats::ENDINDEX - Stats::BEGININDEX];
  std::string Stats::units_[Stats::ENDINDEX - Stats::BEGININDEX];
  Stats::StatsMapType Stats::statsMap;
  Stats::Stats() {

    if (!isInit_) {
      init();
      isInit_ = true;
    }

  }

  void Stats::init() {

    Stats::title_[TIME] = "Time";
    Stats::title_[TOTAL_ENERGY] = "Total Energy";
    Stats::title_[POTENTIAL_ENERGY] = "Potential Energy";
    Stats::title_[KINETIC_ENERGY] = "Kinetic Energy";
    Stats::title_[TEMPERATURE] = "Temperature";
    Stats::title_[PRESSURE] = "Pressure";
    Stats::title_[VOLUME] = "Volume";
    Stats::title_[HULLVOLUME] = "Hull Volume";
    Stats::title_[GYRVOLUME] = "Gyrational Volume";
    Stats::title_[CONSERVED_QUANTITY] = "Conserved Quantity";             
    Stats::title_[TRANSLATIONAL_KINETIC] = "Translational Kinetic";
    Stats::title_[ROTATIONAL_KINETIC] = "Rotational Kinetic";
    Stats::title_[LONG_RANGE_POTENTIAL] = "Long Range Potential";
    Stats::title_[SHORT_RANGE_POTENTIAL] = "Short Range Potential";
    Stats::title_[VANDERWAALS_POTENTIAL] = "van der waals Potential";
    Stats::title_[ELECTROSTATIC_POTENTIAL] = "Electrostatic Potential";    
    Stats::title_[BOND_POTENTIAL] = "Bond Potential";
    Stats::title_[BEND_POTENTIAL] = "Bend Potential";
    Stats::title_[DIHEDRAL_POTENTIAL] = "Dihedral Potential";
    Stats::title_[INVERSION_POTENTIAL] = "Inversion Potential";
    Stats::title_[VRAW] = "Raw Potential";
    Stats::title_[VHARM] = "Harmonic Potential";
    Stats::title_[SHADOWH] = "Shadow Hamiltonian";
    Stats::title_[PRESSURE_TENSOR_XX] = "P_xx";
    Stats::title_[PRESSURE_TENSOR_XY] = "P_xy";
    Stats::title_[PRESSURE_TENSOR_XZ] = "P_xz";
    Stats::title_[PRESSURE_TENSOR_YX] = "P_yx";
    Stats::title_[PRESSURE_TENSOR_YY] = "P_yy";
    Stats::title_[PRESSURE_TENSOR_YZ] = "P_yz";
    Stats::title_[PRESSURE_TENSOR_ZX] = "P_zx";
    Stats::title_[PRESSURE_TENSOR_ZY] = "P_zy";
    Stats::title_[PRESSURE_TENSOR_ZZ] = "P_zz";
    Stats::title_[BOX_DIPOLE_X] = "box dipole x";
    Stats::title_[BOX_DIPOLE_Y] = "box dipole y";
    Stats::title_[BOX_DIPOLE_Z] = "box dipole z";
    Stats::title_[TAGGED_PAIR_DISTANCE] = "Tagged_Pair_Distance";
    Stats::title_[RNEMD_EXCHANGE_TOTAL] = "RNEMD_exchange_total";
    
    Stats::units_[TIME] = "fs";
    Stats::units_[TOTAL_ENERGY] = "kcal/mol";
    Stats::units_[POTENTIAL_ENERGY] = "kcal/mol";
    Stats::units_[KINETIC_ENERGY] = "kcal/mol";
    Stats::units_[TEMPERATURE] = "K";
    Stats::units_[PRESSURE] = "atm";
    Stats::units_[VOLUME] = "A^3";
    Stats::units_[HULLVOLUME] = "A^3";
    Stats::units_[GYRVOLUME] = "A^3";
    Stats::units_[CONSERVED_QUANTITY] = "kcal/mol";             
    Stats::units_[TRANSLATIONAL_KINETIC] = "kcal/mol";
    Stats::units_[ROTATIONAL_KINETIC] = "kcal/mol";
    Stats::units_[LONG_RANGE_POTENTIAL] = "kcal/mol";
    Stats::units_[SHORT_RANGE_POTENTIAL] = "kcal/mol";
    Stats::units_[VANDERWAALS_POTENTIAL] = "kcal/mol";
    Stats::units_[ELECTROSTATIC_POTENTIAL] = "kcal/mol";
    Stats::units_[BOND_POTENTIAL] = "kcal/mol";
    Stats::units_[BEND_POTENTIAL] = "kcal/mol";
    Stats::units_[DIHEDRAL_POTENTIAL] = "kcal/mol";
    Stats::units_[INVERSION_POTENTIAL] = "kcal/mol";
    Stats::units_[VRAW] = "kcal/mol";
    Stats::units_[VHARM] = "kcal/mol";
    Stats::units_[SHADOWH] = "kcal/mol";
    Stats::units_[PRESSURE_TENSOR_XX] = "amu*fs^-2*Ang^-1";
    Stats::units_[PRESSURE_TENSOR_XY] = "amu*fs^-2*Ang^-1";
    Stats::units_[PRESSURE_TENSOR_XZ] = "amu*fs^-2*Ang^-1";
    Stats::units_[PRESSURE_TENSOR_YX] = "amu*fs^-2*Ang^-1";
    Stats::units_[PRESSURE_TENSOR_YY] = "amu*fs^-2*Ang^-1";
    Stats::units_[PRESSURE_TENSOR_YZ] = "amu*fs^-2*Ang^-1";
    Stats::units_[PRESSURE_TENSOR_ZX] = "amu*fs^-2*Ang^-1";
    Stats::units_[PRESSURE_TENSOR_ZY] = "amu*fs^-2*Ang^-1";
    Stats::units_[PRESSURE_TENSOR_ZZ] = "amu*fs^-2*Ang^-1";
    Stats::units_[BOX_DIPOLE_X] = "C*m";
    Stats::units_[BOX_DIPOLE_Y] = "C*m";
    Stats::units_[BOX_DIPOLE_Z] = "C*m";
    Stats::units_[TAGGED_PAIR_DISTANCE] = "Ang";
    Stats::units_[RNEMD_EXCHANGE_TOTAL] = "Variable";

    Stats::statsMap.insert(StatsMapType::value_type("TIME", TIME));
    Stats::statsMap.insert(StatsMapType::value_type("TOTAL_ENERGY", TOTAL_ENERGY));
    Stats::statsMap.insert(StatsMapType::value_type("POTENTIAL_ENERGY", POTENTIAL_ENERGY));
    Stats::statsMap.insert(StatsMapType::value_type("KINETIC_ENERGY", KINETIC_ENERGY));
    Stats::statsMap.insert(StatsMapType::value_type("TEMPERATURE", TEMPERATURE));
    Stats::statsMap.insert(StatsMapType::value_type("PRESSURE", PRESSURE));
    Stats::statsMap.insert(StatsMapType::value_type("VOLUME", VOLUME));
    Stats::statsMap.insert(StatsMapType::value_type("HULLVOLUME", HULLVOLUME));
    Stats::statsMap.insert(StatsMapType::value_type("GYRVOLUME", GYRVOLUME));
    Stats::statsMap.insert(StatsMapType::value_type("CONSERVED_QUANTITY", CONSERVED_QUANTITY));
    Stats::statsMap.insert(StatsMapType::value_type("TRANSLATIONAL_KINETIC", TRANSLATIONAL_KINETIC));
    Stats::statsMap.insert(StatsMapType::value_type("ROTATIONAL_KINETIC", ROTATIONAL_KINETIC));
    Stats::statsMap.insert(StatsMapType::value_type("LONG_RANGE_POTENTIAL", LONG_RANGE_POTENTIAL));
    Stats::statsMap.insert(StatsMapType::value_type("SHORT_RANGE_POTENTIAL", SHORT_RANGE_POTENTIAL));
    Stats::statsMap.insert(StatsMapType::value_type("VANDERWAALS_POTENTIAL", VANDERWAALS_POTENTIAL));
    Stats::statsMap.insert(StatsMapType::value_type("ELECTROSTATIC_POTENTIAL", ELECTROSTATIC_POTENTIAL));
    Stats::statsMap.insert(StatsMapType::value_type("BOND_POTENTIAL", BOND_POTENTIAL));
    Stats::statsMap.insert(StatsMapType::value_type("BEND_POTENTIAL", BEND_POTENTIAL));
    Stats::statsMap.insert(StatsMapType::value_type("DIHEDRAL_POTENTIAL", DIHEDRAL_POTENTIAL));
    Stats::statsMap.insert(StatsMapType::value_type("INVERSION_POTENTIAL", INVERSION_POTENTIAL));
    Stats::statsMap.insert(StatsMapType::value_type("VRAW", VRAW));    
    Stats::statsMap.insert(StatsMapType::value_type("VHARM", VHARM));    
    Stats::statsMap.insert(StatsMapType::value_type("PRESSURE_TENSOR_XX", PRESSURE_TENSOR_XX));    
    Stats::statsMap.insert(StatsMapType::value_type("PRESSURE_TENSOR_XY", PRESSURE_TENSOR_XY));    
    Stats::statsMap.insert(StatsMapType::value_type("PRESSURE_TENSOR_XZ", PRESSURE_TENSOR_XZ));    
    Stats::statsMap.insert(StatsMapType::value_type("PRESSURE_TENSOR_YX", PRESSURE_TENSOR_YX));    
    Stats::statsMap.insert(StatsMapType::value_type("PRESSURE_TENSOR_YY", PRESSURE_TENSOR_YY));    
    Stats::statsMap.insert(StatsMapType::value_type("PRESSURE_TENSOR_YZ", PRESSURE_TENSOR_YZ));    
    Stats::statsMap.insert(StatsMapType::value_type("PRESSURE_TENSOR_ZX", PRESSURE_TENSOR_ZX));    
    Stats::statsMap.insert(StatsMapType::value_type("PRESSURE_TENSOR_ZY", PRESSURE_TENSOR_ZY));    
    Stats::statsMap.insert(StatsMapType::value_type("PRESSURE_TENSOR_ZZ", PRESSURE_TENSOR_ZZ));    
    Stats::statsMap.insert(StatsMapType::value_type("BOX_DIPOLE_X", BOX_DIPOLE_X));    
    Stats::statsMap.insert(StatsMapType::value_type("BOX_DIPOLE_Y", BOX_DIPOLE_Y));    
    Stats::statsMap.insert(StatsMapType::value_type("BOX_DIPOLE_Z", BOX_DIPOLE_Z));    
    Stats::statsMap.insert(StatsMapType::value_type("TAGGED_PAIR_DISTANCE", TAGGED_PAIR_DISTANCE));    
    Stats::statsMap.insert(StatsMapType::value_type("RNEMD_EXCHANGE_TOTAL", RNEMD_EXCHANGE_TOTAL));    
    Stats::statsMap.insert(StatsMapType::value_type("SHADOWH", SHADOWH));    
  }

}
