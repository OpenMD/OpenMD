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
 * [3]  Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).          
 * [4]  Kuang & Gezelter,  J. Chem. Phys. 133, 164101 (2010).
 * [5]  Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 */
#define __OPENMD_C
#include "io/ForceFieldOptions.hpp"

namespace OpenMD {

    ForceFieldOptions::ForceFieldOptions() {
      DefineOptionalParameter(Name, "Name");
      DefineOptionalParameterWithDefaultValue(vdWtype, "vdWtype", "Lennard-Jones");
      DefineOptionalParameterWithDefaultValue(DistanceMixingRule, "DistanceMixingRule", "arithmetic");
      DefineOptionalParameterWithDefaultValue(DistanceType, "DistanceType", "sigma");
      DefineOptionalParameterWithDefaultValue(EnergyMixingRule, "EnergyMixingRule", "geometric");
      DefineOptionalParameterWithDefaultValue(EnergyUnitScaling, "EnergyUnitScaling", 1.0);
      DefineOptionalParameterWithDefaultValue(MetallicEnergyUnitScaling, "MetallicEnergyUnitScaling", 1.0);
      DefineOptionalParameterWithDefaultValue(DistanceUnitScaling, "DistanceUnitScaling", 1.0);
      DefineOptionalParameterWithDefaultValue(AngleUnitScaling, "AngleUnitScaling", 1.0);
      DefineOptionalParameterWithDefaultValue(TorsionAngleConvention, "TorsionAngleConvention", "180_is_trans");
      DefineOptionalParameter(vdw12scale, "vdW-12-scale");
      DefineOptionalParameter(vdw13scale, "vdW-13-scale");
      DefineOptionalParameter(vdw14scale, "vdW-14-scale");
      DefineOptionalParameter(electrostatic12scale, "electrostatic-12-scale");
      DefineOptionalParameter(electrostatic13scale, "electrostatic-13-scale");
      DefineOptionalParameter(electrostatic14scale, "electrostatic-14-scale");

      // DefineOptionalParameterWithDefaultValue(vdw12scale, "vdW-12-scale", 0.0);
      // DefineOptionalParameterWithDefaultValue(vdw13scale, "vdW-13-scale", 0.0);
      // DefineOptionalParameterWithDefaultValue(vdw14scale, "vdW-14-scale", 0.0);
      // DefineOptionalParameterWithDefaultValue(electrostatic12scale, "electrostatic-12-scale", 0.0);
      // DefineOptionalParameterWithDefaultValue(electrostatic13scale, "electrostatic-13-scale", 0.0);
      // DefineOptionalParameterWithDefaultValue(electrostatic14scale, "electrostatic-14-scale", 0.0);
      DefineOptionalParameterWithDefaultValue(GayBerneMu, "GayBerneMu", 2.0);
      DefineOptionalParameterWithDefaultValue(GayBerneNu, "GayBerneNu", 1.0);
      DefineOptionalParameterWithDefaultValue(EAMMixingMethod, "EAMMixingMethod", "Johnson");

      deprecatedKeywords_.insert("cutoffPolicy");
    }
}
