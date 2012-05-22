/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
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
 
#ifndef IO_FORCEFIELDOPTIONS_HPP
#define IO_FORCEFIELDOPTIONS_HPP
#include "utils/simError.h"
#include "types/DataHolder.hpp"
#include "utils/ParameterManager.hpp"
#include "utils/StringUtils.hpp"
#include "io/ParamConstraint.hpp"

namespace OpenMD {
  
  class ForceFieldOptions : public DataHolder {
    DeclareParameter(Name, std::string);
    DeclareParameter(vdWtype, std::string);
    DeclareParameter(DistanceMixingRule, std::string);
    DeclareParameter(DistanceType, std::string);
    DeclareParameter(EnergyMixingRule, std::string);
    DeclareParameter(CutoffPolicy, std::string);
    DeclareParameter(EnergyUnitScaling, RealType);
    DeclareParameter(MetallicEnergyUnitScaling, RealType);
    DeclareParameter(DistanceUnitScaling, RealType);
    DeclareParameter(AngleUnitScaling, RealType);
    DeclareParameter(TorsionAngleConvention, std::string);
    DeclareParameter(vdw12scale, RealType);
    DeclareParameter(vdw13scale, RealType);
    DeclareParameter(vdw14scale, RealType);
    DeclareParameter(electrostatic12scale, RealType);
    DeclareParameter(electrostatic13scale, RealType);
    DeclareParameter(electrostatic14scale, RealType);
    DeclareParameter(GayBerneMu, RealType);
    DeclareParameter(GayBerneNu, RealType);
    DeclareParameter(EAMMixingMethod, std::string);
    
  public:
    ForceFieldOptions();
    ForceFieldOptions(const ForceFieldOptions&);
    ForceFieldOptions& operator = (const ForceFieldOptions&);
    
    void validateOptions() {
      CheckParameter(vdWtype, isEqualIgnoreCase(std::string("Lennard-Jones")));
      CheckParameter(DistanceMixingRule, isEqualIgnoreCase(std::string("arithmetic")) || isEqualIgnoreCase(std::string("geometric")) || isEqualIgnoreCase(std::string("cubic")));
      CheckParameter(DistanceType, isEqualIgnoreCase(std::string("sigma")) || isEqualIgnoreCase(std::string("Rmin")));
      CheckParameter(EnergyMixingRule, isEqualIgnoreCase(std::string("arithmetic")) || isEqualIgnoreCase(std::string("geometric")) || isEqualIgnoreCase(std::string("hhg")));
      CheckParameter(TorsionAngleConvention, isEqualIgnoreCase(std::string("180_is_trans")) || isEqualIgnoreCase(std::string("0_is_trans")));
      CheckParameter(CutoffPolicy, isEqualIgnoreCase(std::string("MIX")) || isEqualIgnoreCase(std::string("MAX")) || isEqualIgnoreCase(std::string("TRADITIONAL")));
      CheckParameter(EAMMixingMethod, isEqualIgnoreCase(std::string("JOHNSON")) || isEqualIgnoreCase(std::string("DAW")));
   }
    
    bool setData(const std::string& keyword, const std::string& value) {
      bool result;
      ParamMap::iterator i =parameters_.find(keyword);
      if (i != parameters_.end()) {
        if(isInteger(value)){
          int ival = lexi_cast<int>(value);
          result = i->second->setData(ival);
        }      
        else if (isType<RealType>(value)){
          RealType dval = lexi_cast<RealType>(value);
          result = i->second->setData(dval);
        } else{
          result = i->second->setData(value);
        }

        if (!result) {
          sprintf(painCave.errMsg,  
                  "Unrecognized data type for keyword: %s = %s\n",
                  keyword.c_str(), value.c_str() );
          painCave.isFatal = 1;
          simError();                  
        }
      } else {
        sprintf(painCave.errMsg,  "%s is an unrecognized keyword\n", keyword.c_str() );
        painCave.isFatal = 0;
        simError();        
      }
      
      return result;
    }

  private:
    typedef std::map<std::string, ParameterBase*> ParamMap;
    ParamMap parameters_;                  
  };
  
}
#endif
