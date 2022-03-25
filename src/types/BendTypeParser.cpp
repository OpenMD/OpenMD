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
#include "types/BendTypeParser.hpp"

#include <string>

#include "io/BendTypesSectionParser.hpp"
#include "types/CosineBendType.hpp"
#include "types/CosineSeriesBendType.hpp"
#include "types/CubicBendType.hpp"
#include "types/HarmonicBendType.hpp"
#include "types/HarmonicSineBendType.hpp"
#include "types/PolynomialBendType.hpp"
#include "types/QuarticBendType.hpp"
#include "types/SDKBendType.hpp"
#include "types/UreyBradleyBendType.hpp"
#include "utils/Constants.hpp"
#include "utils/OpenMDException.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  BendTypeParser::BendTypeParser() {
    stringToEnumMap_["Harmonic"]     = btHarmonic;
    stringToEnumMap_["GhostBend"]    = btGhostBend;
    stringToEnumMap_["UreyBradley"]  = btUreyBradley;
    stringToEnumMap_["Cubic"]        = btCubic;
    stringToEnumMap_["Quartic"]      = btQuartic;
    stringToEnumMap_["Polynomial"]   = btPolynomial;
    stringToEnumMap_["Cosine"]       = btCosine;
    stringToEnumMap_["SDK"]          = btSDK;
    stringToEnumMap_["CosineSeries"] = btCosineSeries;
    stringToEnumMap_["HarmonicSine"] = btHarmonicSine;
  }

  BendType* BendTypeParser::parseTypeAndPars(const std::string& type,
                                             std::vector<RealType> pars) {
    std::string line(type);

    std::vector<RealType>::iterator it;
    for (it = pars.begin(); it != pars.end(); ++it) {
      line.append("\t");
      line.append(std::to_string(*it));
    }
    // Assume all overrides know about our functional forms:
    return parseLine(line, 1.0);
  }

  BendType* BendTypeParser::parseLine(const std::string& line,
                                      RealType kScale) {
    StringTokenizer tokenizer(line);
    BendType* bendType = NULL;
    int nTokens        = tokenizer.countTokens();
    RealType theta0, ktheta;
    RealType deg2rad = Constants::PI / 180.0;

    if (nTokens < 2) {
      throw OpenMDException("BendTypeParser: Not enough tokens");
    }

    BendTypeEnum bt = getBendTypeEnum(tokenizer.nextToken());
    nTokens -= 1;

    // switch is a nightmare to maintain
    switch (bt) {
    case btHarmonic:

      if (nTokens < 2) {
        throw OpenMDException("BendTypeParser: Not enough tokens");
      } else {
        theta0   = tokenizer.nextTokenAsDouble() * deg2rad;
        ktheta   = tokenizer.nextTokenAsDouble() * kScale;
        bendType = new HarmonicBendType(theta0, ktheta);
      }
      break;
    case btGhostBend:
      if (nTokens < 2) {
        throw OpenMDException("BendTypeParser: Not enough tokens");
      } else {
        theta0   = tokenizer.nextTokenAsDouble() * deg2rad;
        ktheta   = tokenizer.nextTokenAsDouble() * kScale;
        bendType = new HarmonicBendType(theta0, ktheta);
      }
      break;

    case btUreyBradley:
      if (nTokens < 4) {
        throw OpenMDException("BendTypeParser: Not enough tokens");
      } else {
        theta0       = tokenizer.nextTokenAsDouble() * deg2rad;
        ktheta       = tokenizer.nextTokenAsDouble();
        RealType s0  = tokenizer.nextTokenAsDouble();
        RealType kub = tokenizer.nextTokenAsDouble();
        bendType     = new UreyBradleyBendType(theta0, ktheta, s0, kub);
      }
      break;

    case btCubic:
      if (nTokens < 5) {
        throw OpenMDException("BendTypeParser: Not enough tokens");
      } else {
        theta0      = tokenizer.nextTokenAsDouble() * deg2rad;
        RealType k3 = tokenizer.nextTokenAsDouble();
        RealType k2 = tokenizer.nextTokenAsDouble();
        RealType k1 = tokenizer.nextTokenAsDouble();
        RealType k0 = tokenizer.nextTokenAsDouble();

        bendType = new CubicBendType(theta0, k3, k2, k1, k0);
      }
      break;

    case btQuartic:
      if (nTokens < 6) {
        throw OpenMDException("BendTypeParser: Not enough tokens");
      } else {
        theta0      = tokenizer.nextTokenAsDouble() * deg2rad;
        RealType k4 = tokenizer.nextTokenAsDouble();
        RealType k3 = tokenizer.nextTokenAsDouble();
        RealType k2 = tokenizer.nextTokenAsDouble();
        RealType k1 = tokenizer.nextTokenAsDouble();
        RealType k0 = tokenizer.nextTokenAsDouble();

        bendType = new QuarticBendType(theta0, k4, k3, k2, k1, k0);
      }
      break;

    case btPolynomial: {
      theta0 = tokenizer.nextTokenAsDouble() * deg2rad;
      nTokens -= 1;
      if (nTokens < 2 || nTokens % 2 != 0) {
        throw OpenMDException("BendTypeParser: Not enough tokens");
      } else {
        int nPairs = nTokens / 2;
        int power;
        RealType coefficient;
        PolynomialBendType* pbt = new PolynomialBendType(theta0);

        for (int i = 0; i < nPairs; ++i) {
          power       = tokenizer.nextTokenAsInt();
          coefficient = tokenizer.nextTokenAsDouble();
          pbt->setCoefficient(power, coefficient);
        }
      }
    } break;

    case btCosine:

      if (nTokens < 2) {
        throw OpenMDException("BendTypeParser: Not enough tokens");
      } else {
        theta0   = tokenizer.nextTokenAsDouble() * deg2rad;
        ktheta   = tokenizer.nextTokenAsDouble();
        bendType = new CosineBendType(theta0, ktheta);
      }
      break;

    case btSDK:
      if (nTokens < 6) {
        throw OpenMDException("BendTypeParser: Not enough tokens");
      } else {
        theta0           = tokenizer.nextTokenAsDouble() * deg2rad;
        ktheta           = tokenizer.nextTokenAsDouble();
        RealType sigma   = tokenizer.nextTokenAsDouble();
        RealType epsilon = tokenizer.nextTokenAsDouble();
        int nRep         = tokenizer.nextTokenAsInt();
        int mAtt         = tokenizer.nextTokenAsInt();
        bendType = new SDKBendType(theta0, ktheta, sigma, epsilon, nRep, mAtt);
      }
      break;

    case btCosineSeries:

      if (nTokens < 2) {
        throw OpenMDException("BendTypeParser: Not enough tokens");
      } else {
        theta0   = tokenizer.nextTokenAsDouble() * deg2rad;
        ktheta   = tokenizer.nextTokenAsDouble() * kScale;
        bendType = new CosineSeriesBendType(theta0, ktheta);
      }
      break;

    case btHarmonicSine:

      if (nTokens < 1) {
        throw OpenMDException("BendTypeParser: Not enough tokens");
      } else {
        ktheta   = tokenizer.nextTokenAsDouble() * kScale;
        bendType = new HarmonicSineBendType(ktheta);
      }
      break;

    case btUnknown:
    default:
      throw OpenMDException("BendTypeParser: Unknown Bend Type");
    }

    return bendType;
  }

  BendTypeParser::BendTypeEnum BendTypeParser::getBendTypeEnum(
      const std::string& str) {
    std::map<std::string, BendTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);

    return i == stringToEnumMap_.end() ? btUnknown : i->second;
  }

}  // namespace OpenMD
