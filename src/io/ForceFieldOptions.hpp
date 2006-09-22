/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
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
 */
 
#ifndef IO_FORCEFIELDOPTIONS_HPP
#define IO_FORCEFIELDOPTIONS_HPP
#include "utils/simError.h"
#include "utils/ParameterManager.hpp"
#include "utils/StringUtils.hpp"
#include "io/ParamConstraint.hpp"
#define __C
#include "UseTheForce/fForceOptions.h"

namespace oopse {
  
  class ForceFieldOptions {
    DeclareParameter(Name, std::string);
    DeclareParameter(vdWtype, std::string);
    DeclareParameter(DistanceMixingRule, std::string);
    DeclareParameter(DistanceType, std::string);
    DeclareParameter(EnergyMixingRule, std::string);
    DeclareParameter(CutoffPolicy, std::string);
    DeclareParameter(EnergyUnitScaling, RealType);
    DeclareParameter(DistanceUnitScaling, RealType);
    DeclareParameter(AngleUnitScaling, RealType);
    DeclareParameter(TorsionAngleConvention, std::string);
    DeclareParameter(vdw14scale, RealType);
    DeclareParameter(electrostatic14scale, RealType);
    DeclareParameter(GayBerneMu, RealType);
    DeclareParameter(GayBerneNu, RealType);
    
  public:
    ForceFieldOptions();
    ForceFieldOptions(const ForceFieldOptions&);
    ForceFieldOptions& operator = (const ForceFieldOptions&);
    
    void validateOptions() {
      CheckParameter(vdWtype, isEqualIgnoreCase(std::string("Lennard-Jones")));
      CheckParameter(DistanceMixingRule, isEqualIgnoreCase(std::string("arithmetic")) || isEqualIgnoreCase(std::string("geometric")) || isEqualIgnoreCase(std::string("cubic")));
      CheckParameter(DistanceType, isEqualIgnoreCase(std::string("sigma")) || isEqualIgnoreCase(std::string("Rmin")));
      CheckParameter(EnergyMixingRule, isEqualIgnoreCase(std::string("arithmetic")) || isEqualIgnoreCase(std::string("geometric")) || isEqualIgnoreCase(std::string("hhg")));
      CheckParameter(TorsionAngleConvention, isEqualIgnoreCase(std::string("180 is trans")) || isEqualIgnoreCase(std::string("0 is trans")));
      CheckParameter(CutoffPolicy, isEqualIgnoreCase(std::string("MIX")) || isEqualIgnoreCase(std::string("MAX")) || isEqualIgnoreCase(std::string("TRADITIONAL")));
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
      } else {
        sprintf(painCave.errMsg,  "%s is an unrecognized keyword\n", keyword.c_str() );
        painCave.isFatal = 0;
        simError();        
      }
      
      return result;
    }

    void makeFortranOptions(ForceOptions & fortranForceOptions);
  private:
    typedef std::map<std::string, ParameterBase*> ParamMap;
    ParamMap parameters_;                  
  };
  
}
#endif
