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

#ifndef BRAINS_STATS_HPP
#define BRAINS_STATS_HPP

#include <string>
#include <map>
#include <bitset>

#include "math/SquareMatrix3.hpp"
#include "utils/Accumulator.hpp"
#include "brains/SimInfo.hpp"

using namespace std;
namespace OpenMD {

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
      ELECTRONIC_KINETIC,
      LONG_RANGE_POTENTIAL,
      VANDERWAALS_POTENTIAL,
      ELECTROSTATIC_POTENTIAL,
      METALLIC_POTENTIAL,
      METALLIC_EMBEDDING,
      METALLIC_PAIR,
      HYDROGENBONDING_POTENTIAL,
      RECIPROCAL_POTENTIAL,
      SURFACE_POTENTIAL,
      SHORT_RANGE_POTENTIAL,
      BOND_POTENTIAL,
      BEND_POTENTIAL,
      DIHEDRAL_POTENTIAL,
      INVERSION_POTENTIAL,
      RAW_POTENTIAL,
      RESTRAINT_POTENTIAL,
      EXCLUDED_POTENTIAL,
      PRESSURE_TENSOR,
      VIRIAL_TENSOR,
      SYSTEM_DIPOLE,
      SYSTEM_QUADRUPOLE,
      TAGGED_PAIR_DISTANCE,
      SHADOWH,
      HELFANDMOMENT,
      HEATFLUX,
      ELECTRONIC_TEMPERATURE,
      COM,
      COM_VELOCITY,
      ANGULAR_MOMENTUM,
      POTENTIAL_SELECTION,
      NET_CHARGE,
      CHARGE_MOMENTUM,
      CURRENT_DENSITY,
      ENDINDEX  //internal use
    };

    struct StatsData {
      std::string title;
      std::string units;
      std::string dataType;
      BaseAccumulator* accumulator;
      std::vector<BaseAccumulator*> accumulatorArray2d;
    };

    typedef bitset<ENDINDEX-BEGININDEX> StatsBitSet;
    typedef map<std::string, StatsIndex> StatsMapType;

    Stats(SimInfo* info);
    virtual ~Stats();
    void parseStatFileFormat(const std::string& format);
    int getPrecision();
    void collectStats();
    std::string getStatsReport();

    StatsBitSet  getStatsMask();
    StatsMapType getStatsMap();
    void         setStatsMask(StatsBitSet mask);

    string    getTitle(int index);
    string    getUnits(int index);
    string    getDataType(int index);

    int       getIntData(int index);
    RealType  getRealData(int index);
    Vector3d  getVectorData(int index);
    potVec    getPotVecData(int index);
    Mat3x3d   getMatrixData(int index);
    std::vector<RealType> getArrayData(int index);

    int       getIntAverage(int index);
    RealType  getRealAverage(int index);
    Vector3d  getVectorAverage(int index);
    potVec    getPotVecAverage(int index);
    Mat3x3d   getMatrixAverage(int index);
    std::vector<RealType> getArrayAverage(int index);

    int       getIntError(int index);
    RealType  getRealError(int index);
    Vector3d  getVectorError(int index);
    potVec    getPotVecError(int index);
    Mat3x3d   getMatrixError(int index);
    std::vector<RealType> getArrayError(int index);

  private:
    SimInfo* info_ {nullptr};
    void init();
    bool isInit_;
    std::vector<StatsData> data_;
    StatsBitSet statsMask_;
    StatsMapType statsMap_;
  };
}
#endif
