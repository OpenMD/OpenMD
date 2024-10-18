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

#include "flucq/FluctuatingChargeParameters.hpp"

#include <cstdlib>
#include <cstring>
#include <iostream>

namespace OpenMD {
  FluctuatingChargeParameters::FluctuatingChargeParameters() {
    DefineOptionalParameterWithDefaultValue(Propagator, "propagator", "Damped");
    DefineOptionalParameterWithDefaultValue(DoInitialOptimization,
                                            "doInitialOptimization", true);
    DefineOptionalParameterWithDefaultValue(ChargeOptimizationMethod, "method",
                                            "CG");
    DefineOptionalParameterWithDefaultValue(Friction, "friction", 1600.0);
    DefineOptionalParameterWithDefaultValue(Tolerance, "tolerance", 1.0e-5);
    DefineOptionalParameterWithDefaultValue(MaxIterations, "maxIterations",
                                            100);
    DefineOptionalParameterWithDefaultValue(TargetTemp, "targetTemp", 1.0e-6);
    DefineOptionalParameterWithDefaultValue(TauThermostat, "tauThermostat",
                                            10.0);
    DefineOptionalParameterWithDefaultValue(DragCoefficient, "dragCoefficient",
                                            0.01);
    DefineOptionalParameterWithDefaultValue(ConstrainRegions,
                                            "constrainRegions", false);
    DefineOptionalParameterWithDefaultValue(InitialStepSize, "initialStepSize",
                                            0.1);
  }

  void FluctuatingChargeParameters::validate() {
    CheckParameter(Propagator,
                   isEqualIgnoreCase("Damped") || isEqualIgnoreCase("NVT") ||
                       isEqualIgnoreCase("Langevin") ||
                       isEqualIgnoreCase("Minimizer") ||
                       isEqualIgnoreCase("NVE") || isEqualIgnoreCase("Exact"));
    CheckParameter(Friction, isNonNegative());
    CheckParameter(Tolerance, isPositive());
    CheckParameter(ChargeOptimizationMethod, isEqualIgnoreCase("SD") ||
                                                 isEqualIgnoreCase("CG") ||
                                                 isEqualIgnoreCase("BFGS"));
    CheckParameter(MaxIterations, isPositive());
    CheckParameter(TargetTemp, isNonNegative());
    CheckParameter(TauThermostat, isPositive());
    CheckParameter(DragCoefficient, isPositive());
    RealType zero = 0.0;
    RealType one  = 1.0;
    CheckParameter(InitialStepSize,
                   isGreaterThan(zero) && isLessThanOrEqualTo(one));
  }

}  // namespace OpenMD
