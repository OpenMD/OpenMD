/*
 * Copyright (c) 2012 The University of Notre Dame. All Rights Reserved.
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

#include "io/FluctuatingChargeAtomTypesSectionParser.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "io/ForceFieldOptions.hpp"
#include "types/AtomType.hpp"
#include "brains/ForceField.hpp"
#include "utils/Tuple.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {

  FluctuatingChargeAtomTypesSectionParser::FluctuatingChargeAtomTypesSectionParser(ForceFieldOptions& options) : options_(options) {
    setSectionName("FluctuatingChargeAtomTypes");

    stringToEnumMap_["Hardness"] =  fqtHardness;
    stringToEnumMap_["MultipleMinima"] =  fqtMultipleMinima;
    stringToEnumMap_["EAM"] =  fqtEAM;
  }

  void FluctuatingChargeAtomTypesSectionParser::parseLine(ForceField& ff,
                                                          const string& line,
                                                          int lineNo){
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();


    if (nTokens < 3)  {
      sprintf(painCave.errMsg,
              "FluctuatingChargeAtomTypesSectionParser Error: "
              "Not enough tokens at line %d\n",
              lineNo);
      painCave.isFatal = 1;
      simError();
    }

    eus_  = options_.getEnergyUnitScaling();

    string atomTypeName = tokenizer.nextToken();
    AtomType* atomType = ff.getAtomType(atomTypeName);
    if (atomType != NULL) {
        FixedChargeAdapter fca = FixedChargeAdapter(atomType);

	// All fluctuating charges are charges, and if we haven't
	// already set values for the charge, then start with zero.
	if (! fca.isFixedCharge()) {
	  RealType charge = 0.0;
	  fca.makeFixedCharge(charge);
	}

      } else {
      sprintf(painCave.errMsg,
              "FluctuatingChargeAtomTypesSectionParser Error: Atom Type [%s] "
              "has not been created yet\n", atomTypeName.c_str());
      painCave.isFatal = 1;
      simError();
    }


    RealType chargeMass = tokenizer.nextTokenAsDouble();
    FluctuatingTypeEnum fqt = getFluctuatingTypeEnum(tokenizer.nextToken());

    nTokens -= 3;

    switch(fqt) {

    case fqtHardness:
      // For Rick, Stuart, Berne style fluctuating charges, there's a
      // self charge potential defined by electronegativity and
      // hardness.  On molecular structures, the slater-type overlap
      // integral is used to compute the hardness.

      // atomTypeName, chargeMass, Hardness, electronegativity,
      //   hardness (Jii), slaterN, slaterZeta

      if (nTokens < 4) {
        sprintf(painCave.errMsg,
                "FluctuatingChargeAtomTypesSectionParser Error: "
                "Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType chi = tokenizer.nextTokenAsDouble();
        RealType Jii = tokenizer.nextTokenAsDouble();
        int slaterN = tokenizer.nextTokenAsInt();
        RealType slaterZeta = tokenizer.nextTokenAsDouble();

        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
        fqa.makeFluctuatingCharge(chargeMass, chi, Jii, slaterN, slaterZeta);
      }
      break;

    case fqtMultipleMinima:
      if (nTokens < 4 || nTokens % 3 != 1) {
        sprintf(painCave.errMsg,
                "FluctuatingChargeAtomTypesSectionParser Error: "
                "Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();

      } else {
        std::vector<tuple3<RealType, RealType, RealType> > diabaticStates;
        RealType curvature = tokenizer.nextTokenAsDouble();
        RealType coupling = tokenizer.nextTokenAsDouble();
        nTokens -= 2;
        int nStates = nTokens / 3;
        RealType charge;
        RealType ionizationEnergy;
        for (int i = 0; i < nStates; ++i) {
          charge = tokenizer.nextTokenAsDouble();
          ionizationEnergy = eus_ * tokenizer.nextTokenAsDouble();
          diabaticStates.push_back( make_tuple3( charge, ionizationEnergy,
                                                 curvature));
        }

        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
        fqa.makeFluctuatingCharge(chargeMass, coupling, diabaticStates);
      }
      break;

      case fqtEAM:
        if (nTokens < 5 || nTokens % 3 != 2) {
          sprintf(painCave.errMsg,
                  "FluctuatingChargeAtomTypesSectionParser Error: "
                  "Not enough tokens at line %d\n",
                  lineNo);
          painCave.isFatal = 1;
          simError();
        } else {
          RealType nValence = tokenizer.nextTokenAsDouble();
          nTokens -= 1;

          std::vector<tuple3<RealType, RealType, RealType> > diabaticStates;
          RealType coupling = tokenizer.nextTokenAsDouble();
          nTokens -= 1;
          int nStates = nTokens / 3;
          RealType charge;
          RealType ionizationEnergy;
          RealType curvature;
          for (int i = 0; i < nStates; ++i) {
            charge = tokenizer.nextTokenAsDouble();
            ionizationEnergy = eus_ * tokenizer.nextTokenAsDouble();
            curvature = tokenizer.nextTokenAsDouble();
            diabaticStates.push_back( make_tuple3( charge, ionizationEnergy,
                                                   curvature ));
          }

          FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
          fqa.makeFluctuatingCharge(chargeMass, nValence, coupling,
                                    diabaticStates);
        }
        break;


    case fqtUnknown:
    default:
      sprintf(painCave.errMsg,
              "FluctuatingChargeAtomTypesSectionParser Error: "
              "Unknown Fluctuating Charge Type at line %d\n",
              lineNo);
      painCave.isFatal = 1;
      simError();
      break;

    }
  }

  FluctuatingChargeAtomTypesSectionParser::FluctuatingTypeEnum FluctuatingChargeAtomTypesSectionParser::getFluctuatingTypeEnum(const std::string& str) {
    std::map<std::string, FluctuatingTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);

    return i == stringToEnumMap_.end() ? fqtUnknown : i->second;
  }

}
