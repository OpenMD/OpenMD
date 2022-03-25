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

#include "types/InversionTypeParser.hpp"

#include <string>

#include "types/AmberImproperTorsionType.hpp"
#include "types/HarmonicInversionType.hpp"
#include "types/ImproperCosineInversionType.hpp"
#include "types/PolynomialInversionType.hpp"
#include "utils/Constants.hpp"
// #include "types/CentralAtomHeightInversionType.hpp"
// #include "types/DreidingInversionType.hpp"
#include "utils/OpenMDException.hpp"
#include "utils/StringTokenizer.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  InversionTypeParser::InversionTypeParser() {
    stringToEnumMap_["AmberImproper"]     = itAmberImproper;
    stringToEnumMap_["ImproperCosine"]    = itImproperCosine;
    stringToEnumMap_["Harmonic"]          = itHarmonic;
    stringToEnumMap_["CentralAtomHeight"] = itCentralAtomHeight;
    stringToEnumMap_["Dreiding"]          = itDreiding;
  }

  InversionType* InversionTypeParser::parseTypeAndPars(
      const std::string& type, std::vector<RealType> pars) {
    std::string line(type);

    std::vector<RealType>::iterator it;
    for (it = pars.begin(); it != pars.end(); ++it) {
      line.append("\t");
      line.append(std::to_string(*it));
    }

    return parseLine(line);
  }

  InversionType* InversionTypeParser::parseLine(const std::string& line) {
    StringTokenizer tokenizer(line);
    InversionType* inversionType = NULL;

    int nTokens = tokenizer.countTokens();

    if (nTokens < 1) {
      throw OpenMDException("InversionTypeParser: Not enough tokens");
    }

    InversionTypeEnum it = getInversionTypeEnum(tokenizer.nextToken());

    nTokens -= 1;

    switch (it) {
    case itImproperCosine:

      if (nTokens < 3 || nTokens % 3 != 0) {
        throw OpenMDException("InversionTypeParser: Not enough tokens");
      } else {
        int nSets = nTokens / 3;

        std::vector<ImproperCosineInversionParameter> parameters;
        for (int i = 0; i < nSets; ++i) {
          ImproperCosineInversionParameter currParam;
          currParam.kchi  = tokenizer.nextTokenAsDouble();
          currParam.n     = tokenizer.nextTokenAsInt();
          currParam.delta = tokenizer.nextTokenAsDouble() / 180.0 *
                            Constants::PI;  // convert to rad
          parameters.push_back(currParam);
        }
        inversionType = new ImproperCosineInversionType(parameters);
      }

      break;

    case itAmberImproper:

      if (nTokens < 1) {
        throw OpenMDException("InversionTypeParser: Not enough tokens");
      } else {
        RealType v2   = tokenizer.nextTokenAsDouble();
        inversionType = new AmberImproperTorsionType(v2);
      }

      break;

    case itHarmonic:
      if (nTokens < 2) {
        throw OpenMDException("InversionTypeParser: Not enough tokens");
      } else {
        // Most inversion don't have specific angle information since
        // they are cosine polynomials.  This one is different,
        // however.  To match our other force field files
        // (particularly for bends):
        //
        // d0 should be read in kcal / mol / degrees^2
        // phi0 should be read in degrees

        RealType degreesPerRadian = 180.0 / Constants::PI;

        // convert to kcal / mol / radians^2
        RealType d0 = tokenizer.nextTokenAsDouble() * pow(degreesPerRadian, 2);

        // convert to radians
        RealType phi0 = tokenizer.nextTokenAsDouble() / degreesPerRadian;
        inversionType = new HarmonicInversionType(d0, phi0);
      }
      break;

      /*
    case itCentralAtomHeight :
      if (nTokens < 1) {
        throw OpenMDException("InversionTypeParser: Not enough tokens");
      } else {

        RealType k = tokenizer.nextTokenAsDouble();

        inversionType = new CentralAtomHeightInversionType(k);
      }
      break;
    case itDreiding :
      if (nTokens < 3) {
        throw OpenMDException("InversionTypeParser: Not enough tokens");
      } else {

        RealType k = tokenizer.nextTokenAsDouble();

        inversionType = new CentralAtomHeightInversionType(k);
      }
      break;
      */
    case itUnknown:
    default:
      throw OpenMDException("InversionTypeParser: Unknown Inversion Type");
    }
    return inversionType;
  }

  InversionTypeParser::InversionTypeEnum
      InversionTypeParser::getInversionTypeEnum(const std::string& str) {
    std::map<std::string, InversionTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);

    return i == stringToEnumMap_.end() ? itUnknown : i->second;
  }

}  // namespace OpenMD
