/*
 * Copyright (c) 2004-2022, The University of Notre Dame. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
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

#ifndef APPLICATIONS_STATICPROPS_RNEMDSTATS_HPP
#define APPLICATIONS_STATICPROPS_RNEMDSTATS_HPP

#include <bitset>
#include <string>
#include <utility>
#include <vector>

#include "applications/staticProps/SpatialStatistics.hpp"
#include "brains/SimInfo.hpp"
#include "math/Vector3.hpp"
#include "primitives/StuntDouble.hpp"
#include "types/AtomType.hpp"

namespace OpenMD {

  class RNEMDZ : public SlabStatistics {
  public:
    RNEMDZ(SimInfo* info, const std::string& filename, const std::string& sele,
           int nzbins, int axis = 2);
    void processFrame(int frame);
    void processStuntDouble(StuntDouble*, int) {}

  protected:
    enum OutputFields {
      BEGININDEX = 0,
      Z          = BEGININDEX,
      TEMPERATURE,
      VELOCITY,
      DENSITY,
      ACTIVITY,
      ELECTRICFIELD,
      ELECTROSTATICPOTENTIAL,
      CHARGE,
      CHARGEVELOCITY,
      ENDINDEX
    };

    using OutputBitSet = std::bitset<ENDINDEX - BEGININDEX>;

    std::vector<AtomType*> outputTypes_;

    OutputData* temperature;
    OutputData* velocity;
    OutputData* density;
    OutputData* activity;
    OutputData* eField;
    OutputData* ePot;
    OutputData* charge;
    OutputData* chargeVelocity;

    OutputBitSet outputMask_;
  };

  class RNEMDR : public ShellStatistics {
  public:
    RNEMDR(SimInfo* info, const std::string& filename, const std::string& sele,
           int nrbins);
    void processFrame(int frame);
    void processStuntDouble(StuntDouble*, int) {}

  protected:
    enum OutputFields {
      BEGININDEX = 0,
      R          = BEGININDEX,
      TEMPERATURE,
      ANGULARVELOCITY,
      DENSITY,
      ENDINDEX
    };

    using OutputBitSet = std::bitset<ENDINDEX - BEGININDEX>;

    OutputData* temperature;
    OutputData* angularVelocity;
    OutputData* density;

    OutputBitSet outputMask_;
  };

  class RNEMDRTheta : public ShellStatistics {
  public:
    RNEMDRTheta(SimInfo* info, const std::string& filename,
                const std::string& sele, int nrbins, int nanglebins);
    void processFrame(int frame);
    void processStuntDouble(StuntDouble*, int) {}
    std::pair<int, int> getBins(Vector3d pos);
    void writeOutput();

  protected:
    enum OutputFields {
      BEGININDEX = 0,
      R          = BEGININDEX,
      ANGULARVELOCITY,
      ENDINDEX
    };

    using OutputBitSet = std::bitset<ENDINDEX - BEGININDEX>;

    int nAngleBins_;
    Vector3d fluxVector_;

    OutputData* angularVelocity;

    OutputBitSet outputMask_;
  };
}  // namespace OpenMD

#endif
