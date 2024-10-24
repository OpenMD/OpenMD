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

#include "io/NonBondedInteractionsSectionParser.hpp"

#include "brains/ForceField.hpp"
#include "types/AtomType.hpp"
#include "types/BuckinghamInteractionType.hpp"
#include "types/EAMInteractionType.hpp"
#include "types/InversePowerSeriesInteractionType.hpp"
#include "types/LennardJonesInteractionType.hpp"
#include "types/MAWInteractionType.hpp"
#include "types/MieInteractionType.hpp"
#include "types/MorseInteractionType.hpp"
#include "types/RepulsivePowerInteractionType.hpp"
#include "utils/simError.h"

namespace OpenMD {

  NonBondedInteractionsSectionParser::NonBondedInteractionsSectionParser(
      ForceFieldOptions& options) :
      options_(options) {
    setSectionName("NonBondedInteractions");

    stringToEnumMap_["MAW"]                = MAW;
    stringToEnumMap_["ShiftedMorse"]       = ShiftedMorse;
    stringToEnumMap_["LennardJones"]       = LennardJones;
    stringToEnumMap_["RepulsiveMorse"]     = RepulsiveMorse;
    stringToEnumMap_["RepulsivePower"]     = RepulsivePower;
    stringToEnumMap_["Mie"]                = Mie;
    stringToEnumMap_["Buckingham"]         = Buckingham;
    stringToEnumMap_["EAMTable"]           = EAMTable;
    stringToEnumMap_["EAMZhou"]            = EAMZhou;
    stringToEnumMap_["EAMOxides"]          = EAMOxides;
    stringToEnumMap_["InversePowerSeries"] = InversePowerSeries;
  }

  void NonBondedInteractionsSectionParser::parseLine(ForceField& ff,
                                                     const std::string& line,
                                                     int lineNo) {
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();
    if (nTokens < 3) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "NonBondedInteractionsSectionParser Error: Not enough tokens at "
               "line %d\n",
               lineNo);
      painCave.isFatal = 1;
      simError();
    }

    meus_ = options_.getMetallicEnergyUnitScaling();
    eus_  = options_.getEnergyUnitScaling();
    dus_  = options_.getDistanceUnitScaling();

    std::string at1   = tokenizer.nextToken();
    std::string at2   = tokenizer.nextToken();
    std::string itype = tokenizer.nextToken();

    NonBondedInteractionTypeEnum nbit = getNonBondedInteractionTypeEnum(itype);
    nTokens -= 3;
    NonBondedInteractionType* interactionType = NULL;

    // switch is a nightmare to maintain
    switch (nbit) {
    case MAW:
      if (nTokens != 5) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "NonBondedInteractionsSectionParser Error: Token number "
                 "mismatch at line "
                 "%d. 8 tokens expected.  \n",
                 lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType r_e    = dus_ * tokenizer.nextTokenAsDouble();
        RealType D_e    = eus_ * tokenizer.nextTokenAsDouble();
        RealType beta   = tokenizer.nextTokenAsDouble() / dus_;
        RealType ca1    = tokenizer.nextTokenAsDouble();
        RealType cb1    = tokenizer.nextTokenAsDouble();
        interactionType = new MAWInteractionType(D_e, beta, r_e, ca1, cb1);
      }
      break;

    case ShiftedMorse:
      if (nTokens != 3) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "NonBondedInteractionsSectionParser Error: Token number "
                 "mismatch at line "
                 "%d. 6 tokens expected.  \n",
                 lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType r0     = dus_ * tokenizer.nextTokenAsDouble();
        RealType D0     = eus_ * tokenizer.nextTokenAsDouble();
        RealType beta0  = tokenizer.nextTokenAsDouble() / dus_;
        interactionType = new MorseInteractionType(D0, beta0, r0, mtShifted);
      }
      break;

    case RepulsiveMorse:
      if (nTokens != 3) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "NonBondedInteractionsSectionParser Error: Token number "
                 "mismatch at line "
                 "%d. 6 tokens expected.  \n",
                 lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType r0     = dus_ * tokenizer.nextTokenAsDouble();
        RealType D0     = eus_ * tokenizer.nextTokenAsDouble();
        RealType beta0  = tokenizer.nextTokenAsDouble() / dus_;
        interactionType = new MorseInteractionType(D0, beta0, r0, mtRepulsive);
      }
      break;

    case LennardJones:
      if (nTokens != 2) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "NonBondedInteractionsSectionParser Error: Token number "
                 "mismatch at line "
                 "%d. 5 tokens expected.  \n",
                 lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType sigma   = dus_ * tokenizer.nextTokenAsDouble();
        RealType epsilon = eus_ * tokenizer.nextTokenAsDouble();
        interactionType  = new LennardJonesInteractionType(sigma, epsilon);
      }
      break;

    case RepulsivePower:
      if (nTokens < 3) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "NonBondedInteractionsSectionParser Error: Token number "
                 "mismatch at line "
                 "%d. 6 tokens expected.  \n",
                 lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType sigma   = dus_ * tokenizer.nextTokenAsDouble();
        RealType epsilon = eus_ * tokenizer.nextTokenAsDouble();
        int nRep         = tokenizer.nextTokenAsInt();
        interactionType =
            new RepulsivePowerInteractionType(sigma, epsilon, nRep);
      }
      break;

    case Mie:
      if (nTokens != 4) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "NonBondedInteractionsSectionParser Error: Token number "
                 "mismatch at line "
                 "%d. 7 tokens expected.  \n",
                 lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType sigma   = dus_ * tokenizer.nextTokenAsDouble();
        RealType epsilon = eus_ * tokenizer.nextTokenAsDouble();
        int nRep         = tokenizer.nextTokenAsInt();
        int mAtt         = tokenizer.nextTokenAsInt();
        interactionType  = new MieInteractionType(sigma, epsilon, nRep, mAtt);
      }
      break;

    case Buckingham:
      if (nTokens < 4) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "NonBondedInteractionsSectionParser Error: Not enough tokens "
                 "at line %d\n",
                 lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        std::string btype = tokenizer.nextToken();
        toUpper(btype);

        RealType A = eus_ * tokenizer.nextTokenAsDouble();
        RealType B = tokenizer.nextTokenAsDouble() / dus_;
        RealType C =
            tokenizer.nextTokenAsDouble();  // should also have a scaling
        RealType sigma   = 0.0;
        RealType epsilon = 0.0;

        if (btype.compare("MODIFIED")) {
          sigma           = dus_ * tokenizer.nextTokenAsDouble();
          epsilon         = eus_ * tokenizer.nextTokenAsDouble();
          interactionType = new BuckinghamInteractionType(A, B, C, sigma,
                                                          epsilon, btModified);

        } else if (btype.compare("TRADITIONAL")) {
          interactionType =
              new BuckinghamInteractionType(A, B, C, btTraditional);
        } else {
          snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                   "NonBondedInteractionsSectionParser Error: Unknown "
                   "Buckingham Type at "
                   "line %d\n",
                   lineNo);
          painCave.isFatal = 1;
          simError();
        }
      }
      break;

    case EAMZhou:
      if (nTokens != 7) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "NonBondedInteractionsSectionParser Error: Token number "
                 "mismatch at line "
                 "%d. 10 tokens expected.  \n",
                 lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType re    = dus_ * tokenizer.nextTokenAsDouble();
        RealType alpha = tokenizer.nextTokenAsDouble();
        RealType beta  = tokenizer.nextTokenAsDouble();
        // Because EAM is a metallic potential, we'll use the metallic
        // unit scaling for these two parameters
        RealType A      = meus_ * tokenizer.nextTokenAsDouble();
        RealType B      = meus_ * tokenizer.nextTokenAsDouble();
        RealType kappa  = tokenizer.nextTokenAsDouble();
        RealType lambda = tokenizer.nextTokenAsDouble();

        interactionType =
            new EAMInteractionType(re, alpha, beta, A, B, kappa, lambda);
      }
      break;

    case EAMOxides:
      if (nTokens != 5) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "NonBondedInteractionsSectionParser Error: Token number "
                 "mismatch at line "
                 "%d. 8 tokens expected.  \n",
                 lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType re    = dus_ * tokenizer.nextTokenAsDouble();
        RealType alpha = tokenizer.nextTokenAsDouble();
        // Because EAM is a metallic potential, we'll use the metallic
        // unit scaling for these two parameters
        RealType A  = meus_ * tokenizer.nextTokenAsDouble();
        RealType Ci = tokenizer.nextTokenAsDouble();
        RealType Cj = tokenizer.nextTokenAsDouble();

        interactionType = new EAMInteractionType(re, alpha, A, Ci, Cj);
      }
      break;

    case InversePowerSeries:
      if (nTokens < 2 || nTokens % 2 != 0) {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "NonBondedInteractionsSectionParser Error: Not enough tokens "
                 "at line %d\n",
                 lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        std::vector<std::pair<int, RealType>> series;
        int nPairs = nTokens / 2;
        int power;
        RealType coefficient;

        for (int i = 0; i < nPairs; ++i) {
          power       = tokenizer.nextTokenAsInt();
          coefficient = tokenizer.nextTokenAsDouble() * eus_ * pow(dus_, power);
          series.push_back(std::make_pair(power, coefficient));
        }
        interactionType = new InversePowerSeriesInteractionType(series);
      }
      break;

    case Unknown:
    default:
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "NonBondedInteractionsSectionParser Error: Unknown Interaction "
               "Type at "
               "line %d\n",
               lineNo);
      painCave.isFatal = 1;
      simError();

      break;
    }

    if (interactionType != NULL) {
      ff.addNonBondedInteractionType(at1, at2, interactionType);
    }
  }

  NonBondedInteractionsSectionParser::NonBondedInteractionTypeEnum
      NonBondedInteractionsSectionParser::getNonBondedInteractionTypeEnum(
          const std::string& str) {
    std::map<std::string, NonBondedInteractionTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);

    return i == stringToEnumMap_.end() ? Unknown : i->second;
  }

}  // namespace OpenMD
