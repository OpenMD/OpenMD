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

#ifndef IO_FORCEFIELDOPTIONS_HPP
#define IO_FORCEFIELDOPTIONS_HPP

#include "io/ParamConstraint.hpp"
#include "types/DataHolder.hpp"
#include "utils/ParameterManager.hpp"
#include "utils/StringUtils.hpp"
#include "utils/simError.h"

namespace OpenMD {

  class ForceFieldOptions : public DataHolder {
    DeclareParameter(Name, std::string);
    DeclareParameter(vdWtype, std::string);
    DeclareParameter(DistanceMixingRule, std::string);
    DeclareParameter(DistanceType, std::string);
    DeclareParameter(EnergyMixingRule, std::string);
    DeclareParameter(EnergyUnitScaling, RealType);
    DeclareParameter(MetallicEnergyUnitScaling, RealType);
    DeclareParameter(FluctuatingChargeEnergyUnitScaling, RealType);
    DeclareParameter(DistanceUnitScaling, RealType);
    DeclareParameter(AngleUnitScaling, RealType);
    DeclareParameter(ChargeUnitScaling, RealType);
    DeclareParameter(OxidationStateScaling, RealType);
    DeclareParameter(TorsionAngleConvention, std::string);
    DeclareParameter(vdw12scale, RealType);
    DeclareParameter(vdw13scale, RealType);
    DeclareParameter(vdw14scale, RealType);
    DeclareParameter(BondForceConstantScaling, RealType);
    DeclareParameter(BendForceConstantScaling, RealType);
    DeclareParameter(electrostatic12scale, RealType);
    DeclareParameter(electrostatic13scale, RealType);
    DeclareParameter(electrostatic14scale, RealType);
    DeclareParameter(GayBerneMu, RealType);
    DeclareParameter(GayBerneNu, RealType);
    DeclareParameter(EAMMixingMethod, std::string);
    DeclareParameter(DelayedParameterCalculation, bool);

  public:
    ForceFieldOptions();
    ForceFieldOptions(const ForceFieldOptions&);
    ForceFieldOptions& operator=(const ForceFieldOptions&);

    void validateOptions() {
      CheckParameter(vdWtype, isEqualIgnoreCase(std::string("Lennard-Jones")));
      CheckParameter(DistanceMixingRule,
                     isEqualIgnoreCase(std::string("arithmetic")) ||
                         isEqualIgnoreCase(std::string("geometric")) ||
                         isEqualIgnoreCase(std::string("cubic")));
      CheckParameter(DistanceType, isEqualIgnoreCase(std::string("sigma")) ||
                                       isEqualIgnoreCase(std::string("Rmin")));
      CheckParameter(EnergyMixingRule,
                     isEqualIgnoreCase(std::string("arithmetic")) ||
                         isEqualIgnoreCase(std::string("geometric")) ||
                         isEqualIgnoreCase(std::string("hhg")));
      CheckParameter(TorsionAngleConvention,
                     isEqualIgnoreCase(std::string("180_is_trans")) ||
                         isEqualIgnoreCase(std::string("0_is_trans")));
      CheckParameter(EAMMixingMethod,
                     isEqualIgnoreCase(std::string("Johnson")) ||
                         isEqualIgnoreCase(std::string("Daw")) ||
                         isEqualIgnoreCase(std::string("DREAM1")) ||
                         isEqualIgnoreCase(std::string("DREAM2")));
    }

    bool setData(const std::string& keyword, const std::string& value) {
      bool result(false);
      ParamMap::iterator i = parameters_.find(keyword);
      if (i != parameters_.end()) {
        if (isInteger(value)) {
          int ival = lexi_cast<int>(value);
          result   = i->second->setData(ival);
        } else if (isType<RealType>(value)) {
          RealType dval = lexi_cast<RealType>(value);
          result        = i->second->setData(dval);
        } else {
          result = i->second->setData(value);
        }

        if (!result) {
          sprintf(painCave.errMsg,
                  "Unrecognized data type for keyword: %s = %s\n",
                  keyword.c_str(), value.c_str());
          painCave.isFatal = 1;
          simError();
        }
      } else {
        sprintf(painCave.errMsg, "%s is an unrecognized keyword\n",
                keyword.c_str());
        painCave.isFatal = 0;
        simError();
      }

      return result;
    }

  private:
    using ParamMap = std::map<std::string, ParameterBase*>;
    ParamMap parameters_;
  };
}  // namespace OpenMD

#endif
