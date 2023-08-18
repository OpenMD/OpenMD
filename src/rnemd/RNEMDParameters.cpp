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

#include "rnemd/RNEMDParameters.hpp"

#include <cstdlib>
#include <cstring>
#include <iostream>

namespace OpenMD::RNEMD {

  RNEMDParameters::RNEMDParameters() {
    DefineOptionalParameterWithDefaultValue(UseRNEMD, "useRNEMD", false);
    DefineOptionalParameterWithDefaultValue(ObjectSelection, "objectSelection",
                                            "select all");

    DefineOptionalParameterWithDefaultValue(Method, "method", "VSS");
    DefineOptionalParameter(FluxType, "fluxType");

    DefineOptionalParameterWithDefaultValue(ExchangeTime, "exchangeTime",
                                            100.0);
    DefineOptionalParameter(KineticFlux, "kineticFlux");
    DefineOptionalParameter(ParticleFlux, "particleFlux");
    DefineOptionalParameter(MomentumFlux, "momentumFlux");
    DefineOptionalParameter(MomentumFluxVector, "momentumFluxVector");
    DefineOptionalParameter(AngularMomentumFlux, "angularMomentumFlux");
    DefineOptionalParameter(AngularMomentumFluxVector,
                            "angularMomentumFluxVector");
    DefineOptionalParameter(SlabWidth, "slabWidth");
    DefineOptionalParameter(SlabACenter, "slabAcenter");
    DefineOptionalParameter(SlabBCenter, "slabBcenter");
    DefineOptionalParameter(SphereARadius, "sphereAradius");
    DefineOptionalParameter(SphereBRadius, "sphereBradius");
    DefineOptionalParameter(SelectionA, "selectionA");
    DefineOptionalParameter(SelectionB, "selectionB");
    DefineOptionalParameter(CoordinateOrigin, "coordinateOrigin");
    DefineOptionalParameter(OutputFileName, "outputFileName");
    DefineOptionalParameterWithDefaultValue(OutputBins, "outputBins", 20);
    DefineOptionalParameterWithDefaultValue(OutputBinWidth, "outputBinWidth",
                                            2.0);
    DefineOptionalParameter(OutputSelection, "outputSelection");
    DefineOptionalParameter(OutputFields, "outputFields");
    DefineOptionalParameter(DividingArea, "dividingArea");
    DefineOptionalParameterWithDefaultValue(PrivilegedAxis, "privilegedAxis",
                                            "z");
    DefineOptionalParameterWithDefaultValue(SPFScalingPower, "spfScalingPower",
                                            3);
    DefineOptionalParameterWithDefaultValue(SPFUniformKineticScaling,
                                            "spfUniformKineticScaling",
                                            false);
  }

  void RNEMDParameters::validate() {
    CheckParameter(ExchangeTime, isPositive());
    CheckParameter(OutputBins, isPositive());
    CheckParameter(OutputBinWidth, isPositive());
    CheckParameter(Method,
                   isEqualIgnoreCase("Swap") || isEqualIgnoreCase("NIVS") ||
                       isEqualIgnoreCase("VSS") || isEqualIgnoreCase("SPF"));
    CheckParameter(
        FluxType,
        isEqualIgnoreCase("KE") || isEqualIgnoreCase("Px") ||
            isEqualIgnoreCase("Py") || isEqualIgnoreCase("Pz") ||
            isEqualIgnoreCase("Lx") || isEqualIgnoreCase("Ly") ||
            isEqualIgnoreCase("Lz") || isEqualIgnoreCase("Pvector") ||
            isEqualIgnoreCase("Lvector") || isEqualIgnoreCase("KE+Px") ||
            isEqualIgnoreCase("KE+Py") || isEqualIgnoreCase("KE+Lx") ||
            isEqualIgnoreCase("KE+Ly") || isEqualIgnoreCase("KE+Lz") ||
            isEqualIgnoreCase("KE+Pvector") ||
            isEqualIgnoreCase("KE+Lvector") || isEqualIgnoreCase("Particle")||
            isEqualIgnoreCase("Particle+KE"));
    CheckParameter(PrivilegedAxis, isEqualIgnoreCase("x") ||
                                       isEqualIgnoreCase("y") ||
                                       isEqualIgnoreCase("z"));
  }

  bool RNEMDParameters::requiresElectricField() {
    static bool wasParsed {false};

    if (!wasParsed) {
      StringTokenizer tokenizer(getOutputFields(), " ,;|\t\n\r");

      while (tokenizer.hasMoreTokens()) {
        std::string token(tokenizer.nextToken());
        toUpper(token);

        if (token == "ELECTRICFIELD" || token == "ELECTROSTATICPOTENTIAL") {
          calculateElectricField_ = true;
          break;
        }
      }

      wasParsed = true;
    }

    return calculateElectricField_;
  }
}  // namespace OpenMD::RNEMD
