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

#include "io/FluctuatingChargeAtomTypesSectionParser.hpp"

#include "brains/ForceField.hpp"
#include "io/ForceFieldOptions.hpp"
#include "types/AtomType.hpp"
#include "types/FixedChargeAdapter.hpp"
#include "types/FluctuatingChargeAdapter.hpp"
#include "utils/simError.h"

using namespace std;
namespace OpenMD {

  FluctuatingChargeAtomTypesSectionParser::
      FluctuatingChargeAtomTypesSectionParser(ForceFieldOptions& options) :
      options_(options) {
    setSectionName("FluctuatingChargeAtomTypes");

    stringToEnumMap_["Hardness"] = fqtHardness;
    stringToEnumMap_["EAM"]      = fqtEAM;
    stringToEnumMap_["EAMPoly"]  = fqtEAMPolynomial;
    stringToEnumMap_["DREAM2"]   = fqtDREAM2;
  }

  void FluctuatingChargeAtomTypesSectionParser::parseLine(ForceField& ff,
                                                          const string& line,
                                                          int lineNo) {
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();

    if (nTokens < 3) {
      sprintf(painCave.errMsg,
              "FluctuatingChargeAtomTypesSectionParser Error: "
              "Not enough tokens at line %d\n",
              lineNo);
      painCave.isFatal = 1;
      simError();
    }

    // conversion to kcal / mol
    eus_ = options_.getFluctuatingChargeEnergyUnitScaling();
    // conversion to electrons
    cus_ = options_.getChargeUnitScaling();
    // conversion to angstroms
    dus_ = options_.getDistanceUnitScaling();
    // oxidation state scaling
    oss_ = options_.getOxidationStateScaling();

    // electronegativity conversion
    RealType chius = eus_ / cus_;
    // curvature (hardness) conversion
    RealType curvus = eus_ / (cus_ * cus_);

    string atomTypeName = tokenizer.nextToken();
    AtomType* atomType  = ff.getAtomType(atomTypeName);
    if (atomType != NULL) {
      FixedChargeAdapter fca = FixedChargeAdapter(atomType);

      // All fluctuating charges are charges, and if we haven't
      // already set values for the charge, then start with zero.
      if (!fca.isFixedCharge()) {
        RealType charge = 0.0;
        fca.makeFixedCharge(charge);
      }

    } else {
      sprintf(painCave.errMsg,
              "FluctuatingChargeAtomTypesSectionParser Error: Atom Type [%s] "
              "has not been created yet\n",
              atomTypeName.c_str());
      painCave.isFatal = 1;
      simError();
    }

    RealType chargeMass     = tokenizer.nextTokenAsDouble();
    FluctuatingTypeEnum fqt = getFluctuatingTypeEnum(tokenizer.nextToken());

    nTokens -= 3;

    switch (fqt) {
    case fqtHardness:
      // For Rick, Stuart, Berne style fluctuating charges, there's a
      // self charge potential defined by electronegativity and
      // hardness. On molecular structures, the slater-type overlap
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
        RealType chi        = chius * tokenizer.nextTokenAsDouble();
        RealType Jii        = curvus * tokenizer.nextTokenAsDouble();
        int slaterN         = tokenizer.nextTokenAsInt();
        RealType slaterZeta = tokenizer.nextTokenAsDouble() / dus_;

        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
        fqa.makeFluctuatingCharge(chargeMass, chi, Jii, slaterN, slaterZeta);
      }
      break;

    case fqtEAMPolynomial:

      // For DR-EAM v1, oxidation state scaling is built in to
      // fluctuating charge parameters and nValence.

      if (nTokens < 3 || nTokens % 2 != 1) {
        sprintf(painCave.errMsg,
                "FluctuatingChargeAtomTypesSectionParser Error: "
                "Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        RealType nValence = tokenizer.nextTokenAsDouble();
        nTokens -= 1;

        DoublePolynomial vself;

        int nPairs = nTokens / 2;
        int power;
        RealType coefficient;

        for (int i = 0; i < nPairs; ++i) {
          power       = tokenizer.nextTokenAsInt();
          coefficient = tokenizer.nextTokenAsDouble() * eus_ / pow(cus_, power);
          vself.setCoefficient(power, coefficient);
        }

        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);
        fqa.makeFluctuatingCharge(chargeMass, nValence, vself);
      }
      break;

    case fqtDREAM2:

      // For DR-EAM v2, oxidation state scaling (oss_) is a force
      // field parameter, fluctuating charge parameters, nValence, and
      // nMobile in the force field file should assume unit charge
      // oxidation state fits to ionization states and electron
      // affinities.

      if (nTokens < 3 || nTokens % 2 != 0) {
        sprintf(painCave.errMsg,
                "FluctuatingChargeAtomTypesSectionParser Error: "
                "Not enough tokens at line %d\n",
                lineNo);
        painCave.isFatal = 1;
        simError();
      } else {
        // Neutral atom valence count
        RealType nValence = tokenizer.nextTokenAsDouble() * oss_;
        nTokens -= 1;
        // Mobile electron count
        RealType nMobile = tokenizer.nextTokenAsDouble() * oss_;
        nTokens -= 1;

        DoublePolynomial vself;

        int nPairs = nTokens / 2;
        int power;
        RealType coefficient;

        for (int i = 0; i < nPairs; ++i) {
          power       = tokenizer.nextTokenAsInt();
          coefficient = oss_ * oss_ * tokenizer.nextTokenAsDouble() * eus_ /
                        pow(oss_ * cus_, power);

          vself.setCoefficient(power, coefficient);
        }

        FluctuatingChargeAdapter fqa = FluctuatingChargeAdapter(atomType);

        fqa.makeFluctuatingCharge(chargeMass, nValence, nMobile, vself);
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

  FluctuatingChargeAtomTypesSectionParser::FluctuatingTypeEnum
      FluctuatingChargeAtomTypesSectionParser::getFluctuatingTypeEnum(
          const std::string& str) {
    std::map<std::string, FluctuatingTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);

    return i == stringToEnumMap_.end() ? fqtUnknown : i->second;
  }

}  // namespace OpenMD
