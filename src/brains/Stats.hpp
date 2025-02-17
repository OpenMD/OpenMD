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
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#ifndef BRAINS_STATS_HPP
#define BRAINS_STATS_HPP

#include <bitset>
#include <map>
#include <string>

#include "brains/SimInfo.hpp"
#include "math/SquareMatrix3.hpp"
#include "utils/OldAccumulator.hpp"

using namespace std;
namespace OpenMD {

  class Stats {
  public:
    enum StatsIndex {
      BEGININDEX = 0,  // internal use
      TIME       = BEGININDEX,
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
      ENDINDEX  // internal use
    };

    struct StatsData {
      std::string title;
      std::string units;
      std::string dataType;
      BaseAccumulator* accumulator;
      std::vector<BaseAccumulator*> accumulatorArray2d;
    };

    using StatsBitSet  = std::bitset<ENDINDEX - BEGININDEX>;
    using StatsMapType = std::map<std::string, StatsIndex>;

    Stats(SimInfo* info);
    virtual ~Stats();
    void parseStatFileFormat(const std::string& format);
    int getPrecision();
    void collectStats();
    std::string getStatsReport();

    StatsBitSet getStatsMask();
    StatsMapType getStatsMap();
    void setStatsMask(StatsBitSet mask);

    string getTitle(int index);
    string getUnits(int index);
    string getDataType(int index);

    int getIntData(int index);
    RealType getRealData(int index);
    Vector3d getVectorData(int index);
    potVec getPotVecData(int index);
    Mat3x3d getMatrixData(int index);
    std::vector<RealType> getArrayData(int index);

    int getIntAverage(int index);
    RealType getRealAverage(int index);
    Vector3d getVectorAverage(int index);
    potVec getPotVecAverage(int index);
    Mat3x3d getMatrixAverage(int index);
    std::vector<RealType> getArrayAverage(int index);

    int getIntError(int index);
    RealType getRealError(int index);
    Vector3d getVectorError(int index);
    potVec getPotVecError(int index);
    Mat3x3d getMatrixError(int index);
    std::vector<RealType> getArrayError(int index);

  private:
    SimInfo* info_ {nullptr};
    void init();
    bool isInit_;
    std::vector<StatsData> data_;
    StatsBitSet statsMask_;
    StatsMapType statsMap_;
  };
}  // namespace OpenMD

#endif
