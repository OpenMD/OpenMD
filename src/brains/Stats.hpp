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
 * @file Stats.hpp
 * @author tlin
 * @date 11/04/2004
 * @time 23:56am
 * @version 1.0
 */

#ifndef BRAINS_STATS_HPP
#define BRAINS_STATS_HPP

#include <string>
#include <map>

#include "math/SquareMatrix3.hpp"
namespace OpenMD {

  /**
   * @class Stats Stats.hpp "brains/Stats.hpp"
   */
  class Stats{
  public:
    enum StatsIndex {
      BEGININDEX = 0,  //internal use
      TIME = BEGININDEX,            
      TOTAL_ENERGY,
      POTENTIAL_ENERGY,
      KINETIC_ENERGY,
      TEMPERATURE,
      PRESSURE,
      VOLUME,
      HULLVOLUME,
      GYRVOLUME,
      CONSERVED_QUANTITY,             
      TRANSLATIONAL_KINETIC,
      ROTATIONAL_KINETIC,
      LONG_RANGE_POTENTIAL,   
      SHORT_RANGE_POTENTIAL,
      VANDERWAALS_POTENTIAL,
      ELECTROSTATIC_POTENTIAL,
      BOND_POTENTIAL,
      BEND_POTENTIAL,
      DIHEDRAL_POTENTIAL,
      INVERSION_POTENTIAL,
      VRAW,
      VHARM,
      PRESSURE_TENSOR_XX,
      PRESSURE_TENSOR_XY,
      PRESSURE_TENSOR_XZ,
      PRESSURE_TENSOR_YX,
      PRESSURE_TENSOR_YY,
      PRESSURE_TENSOR_YZ,
      PRESSURE_TENSOR_ZX,
      PRESSURE_TENSOR_ZY,
      PRESSURE_TENSOR_ZZ,
      BOX_DIPOLE_X,
      BOX_DIPOLE_Y,
      BOX_DIPOLE_Z,
      TAGGED_PAIR_DISTANCE,
      RNEMD_EXCHANGE_TOTAL,
      SHADOWH,
      ENDINDEX  //internal use
    };

    Stats();
    const RealType& operator [](int index) const {
      assert(index >=0 && index < ENDINDEX);
      return data_[index];
    }

    RealType& operator [](int index){
      assert(index >=0 && index < ENDINDEX);            
      return data_[index];
    }
        
    static std::string getTitle(int index) {
      assert(index >=0 && index < ENDINDEX);
      return title_[index];
    }

    static std::string getUnits(int index) {
      assert(index >=0 && index < ENDINDEX);
      return units_[index];
    }

    Mat3x3d getTau() {
      return tau_;
    }
        
    void setTau(const Mat3x3d& tau) {
      tau_ = tau;
    }

    typedef std::map<std::string, Stats::StatsIndex> StatsMapType;
    static  StatsMapType statsMap;
  
  private:
    static void init();
    static bool isInit_;
    RealType data_[ENDINDEX - BEGININDEX];
    static std::string title_[ENDINDEX - BEGININDEX];
    static std::string units_[ENDINDEX - BEGININDEX];
    Mat3x3d tau_;
  };



} //end namespace OpenMD
#endif //BRAINS_STATS_HPP
