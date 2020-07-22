/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
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
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "io/NonBondedInteractionsSectionParser.hpp"
#include "types/AtomType.hpp"
#include "types/MorseInteractionType.hpp"
#include "types/MAWInteractionType.hpp"
#include "types/LennardJonesInteractionType.hpp"
#include "types/RepulsivePowerInteractionType.hpp"
#include "types/MieInteractionType.hpp"
#include "types/BuckinghamInteractionType.hpp"
#include "types/EAMInteractionType.hpp"
#include "types/InversePowerSeriesInteractionType.hpp"
#include "brains/ForceField.hpp"
#include "utils/simError.h"
namespace OpenMD {

  NonBondedInteractionsSectionParser::NonBondedInteractionsSectionParser(ForceFieldOptions& options) : options_(options){
    setSectionName("NonBondedInteractions");

    stringToEnumMap_["MAW"] =  MAW;
    stringToEnumMap_["ShiftedMorse"] =  ShiftedMorse;
    stringToEnumMap_["LennardJones"] = LennardJones;
    stringToEnumMap_["RepulsiveMorse"] = RepulsiveMorse;
    stringToEnumMap_["RepulsivePower"] = RepulsivePower;
    stringToEnumMap_["Mie"] = Mie;
    stringToEnumMap_["Buckingham"] = Buckingham;
    stringToEnumMap_["EAMTable"] = EAMTable;
    stringToEnumMap_["EAMZhou"] = EAMZhou;
    stringToEnumMap_["InversePowerSeries"] = InversePowerSeries;

  }

  void NonBondedInteractionsSectionParser::parseLine(ForceField& ff,const std::string& line, int lineNo){
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();

    if (nTokens < 3) {
      sprintf(painCave.errMsg,
              "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();
    }

    meus_ = options_.getMetallicEnergyUnitScaling();
    eus_  = options_.getEnergyUnitScaling();
    dus_  = options_.getDistanceUnitScaling();

    std::string at1 = tokenizer.nextToken();
    std::string at2 = tokenizer.nextToken();
    std::string itype = tokenizer.nextToken();

    NonBondedInteractionTypeEnum nbit = getNonBondedInteractionTypeEnum(itype);
    nTokens -= 3;
    NonBondedInteractionType* interactionType = NULL;

    //switch is a nightmare to maintain
    switch(nbit) {
    case MAW :
      if (nTokens < 5) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType r_e = dus_ * tokenizer.nextTokenAsDouble();
        RealType D_e = eus_ * tokenizer.nextTokenAsDouble();
        RealType beta = tokenizer.nextTokenAsDouble() / dus_ ;
        RealType ca1 = tokenizer.nextTokenAsDouble();
        RealType cb1 = tokenizer.nextTokenAsDouble();
        interactionType = new MAWInteractionType(D_e, beta, r_e, ca1, cb1);
      }
      break;

    case ShiftedMorse :
      if (nTokens < 3) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType r0 = dus_ * tokenizer.nextTokenAsDouble();
        RealType D0 = eus_ * tokenizer.nextTokenAsDouble();
        RealType beta0 = tokenizer.nextTokenAsDouble() / dus_;
        interactionType = new MorseInteractionType(D0, beta0, r0, mtShifted);
      }
      break;

    case RepulsiveMorse :
      if (nTokens < 3) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType r0 = dus_ * tokenizer.nextTokenAsDouble();
        RealType D0 = eus_ * tokenizer.nextTokenAsDouble();
        RealType beta0 = tokenizer.nextTokenAsDouble() / dus_;
        interactionType = new MorseInteractionType(D0, beta0, r0, mtRepulsive);
      }
      break;

    case LennardJones :
      if (nTokens < 2) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType sigma = dus_ * tokenizer.nextTokenAsDouble();
        RealType epsilon = eus_ * tokenizer.nextTokenAsDouble();
        interactionType = new LennardJonesInteractionType(sigma, epsilon);
      }
      break;

    case RepulsivePower :
      if (nTokens < 3) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType sigma = dus_ * tokenizer.nextTokenAsDouble();
        RealType epsilon = eus_ * tokenizer.nextTokenAsDouble();
        int nRep = tokenizer.nextTokenAsInt();
        interactionType = new RepulsivePowerInteractionType(sigma, epsilon,
                                                            nRep);
      }
      break;

    case Mie :
      if (nTokens < 4) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType sigma = dus_ * tokenizer.nextTokenAsDouble();
        RealType epsilon = eus_ * tokenizer.nextTokenAsDouble();
        int nRep = tokenizer.nextTokenAsInt();
        int mAtt = tokenizer.nextTokenAsInt();
        interactionType = new MieInteractionType(sigma, epsilon, nRep, mAtt);
      }
      break;

    case Buckingham :
      if (nTokens < 4) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        std::string btype = tokenizer.nextToken();
        toUpper(btype);

        RealType A = eus_ * tokenizer.nextTokenAsDouble();
        RealType B = tokenizer.nextTokenAsDouble() / dus_;
        RealType C = tokenizer.nextTokenAsDouble(); // should also have a scaling
        RealType sigma = 0.0;
        RealType epsilon = 0.0;

        if (btype.compare("MODIFIED")) {
          sigma = dus_ * tokenizer.nextTokenAsDouble();
          epsilon = eus_ * tokenizer.nextTokenAsDouble();
          interactionType = new BuckinghamInteractionType(A, B, C, sigma, epsilon, btModified);

        } else if(btype.compare("TRADITIONAL")) {
          interactionType = new BuckinghamInteractionType(A, B, C, btTraditional);
        } else {

          sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Unknown Buckingham Type at line %d\n",
                  lineNo);
          painCave.isFatal = 1;
          simError();
        }

      }
      break;

    case EAMZhou :
      if (nTokens < 7) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {

        RealType re = dus_ * tokenizer.nextTokenAsDouble();
        RealType alpha = tokenizer.nextTokenAsDouble();
        RealType beta = tokenizer.nextTokenAsDouble();
        // Because EAM is a metallic potential, we'll use the metallic
        // unit scaling for these two parameters
        RealType A = meus_ * tokenizer.nextTokenAsDouble();
        RealType B = meus_ * tokenizer.nextTokenAsDouble();
        RealType kappa = tokenizer.nextTokenAsDouble();
        RealType lambda = tokenizer.nextTokenAsDouble();

        interactionType = new EAMInteractionType(re, alpha, beta, A, B,
                                                 kappa, lambda);
      }
      break;

    case InversePowerSeries :
      if (nTokens < 2 || nTokens % 2 != 0) {
        sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {

        std::vector<std::pair<int,RealType> > series;
        int nPairs = nTokens / 2;
        int power;
        RealType coefficient;
        
        for (int i = 0; i < nPairs; ++i) {
          power = tokenizer.nextTokenAsInt();
          coefficient = tokenizer.nextTokenAsDouble() * eus_ * pow(dus_, power);
          series.push_back(std::make_pair( power, coefficient));          
        }
        interactionType = new InversePowerSeriesInteractionType(series);
        
      }
      break;

    case Unknown :
    default:
      sprintf(painCave.errMsg, "NonBondedInteractionsSectionParser Error: Unknown Interaction Type at line %d\n",
	      lineNo);
      painCave.isFatal = 1;
      simError();

      break;

    }

    if (interactionType != NULL) {
      ff.addNonBondedInteractionType(at1, at2, interactionType);
    }

  }

  NonBondedInteractionsSectionParser::NonBondedInteractionTypeEnum NonBondedInteractionsSectionParser::getNonBondedInteractionTypeEnum(const std::string& str) {
    std::map<std::string, NonBondedInteractionTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);

    return i == stringToEnumMap_.end() ? Unknown : i->second;
  }

} //end namespace OpenMD
