#include "io/ForceFieldOptions.hpp"
namespace oopse {
    ForceFieldOptions::ForceFieldOptions() {
      std::cout << "\n";
      DefineOptionalParameter(Name, "Name");
      DefineOptionalParameterWithDefaultValue(vdWtype, "vdWtype", "Lennard-Jones");
      DefineOptionalParameterWithDefaultValue(DistanceMixingRule, "DistanceMixingRule", "arithmetic");
      DefineOptionalParameterWithDefaultValue(DistanceType, "DistanceType", "sigma");
      DefineOptionalParameterWithDefaultValue(EnergyMixingRule, "EnergyMixingRule", "geometric");
      DefineOptionalParameterWithDefaultValue(EnergyUnitScaling, "EnergyUnitScaling", 1.0);
      DefineOptionalParameterWithDefaultValue(DistanceUnitScaling, "DistanceUnitScaling", 1.0);
      DefineOptionalParameterWithDefaultValue(AngleUnitScaling, "AngleUnitScaling", 1.0);
      DefineOptionalParameterWithDefaultValue(TorsionAngleConvention, "TorsionAngleConvention", "180 is trans");
      DefineOptionalParameterWithDefaultValue(vdw14scale, "vdW-14-scale", 0.0);
      DefineOptionalParameterWithDefaultValue(electrostatic14scale, "electrostatic-14-scale", 0.0);
      DefineOptionalParameterWithDefaultValue(dielectric, "dielectric", 1.0);
    }
}
