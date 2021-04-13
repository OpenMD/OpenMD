/*
 * Copyright (c) 2004-2021 The University of Notre Dame. All Rights Reserved.
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

#include "types/BondTypeParser.hpp"

#include "types/CubicBondType.hpp"
#include "types/FixedBondType.hpp"
#include "types/HarmonicBondType.hpp"
#include "types/MorseBondType.hpp"
#include "types/PolynomialBondType.hpp"
#include "types/QuarticBondType.hpp"
#include "types/ShiftedMieBondType.hpp"
#include "utils/OpenMDException.hpp"
#include "utils/StringTokenizer.hpp"
#include "utils/StringUtils.hpp"

namespace OpenMD {

  BondTypeParser::BondTypeParser() {
    stringToEnumMap_["Fixed"]      = btFixed;
    stringToEnumMap_["Harmonic"]   = btHarmonic;
    stringToEnumMap_["Cubic"]      = btCubic;
    stringToEnumMap_["Quartic"]    = btQuartic;
    stringToEnumMap_["Polynomial"] = btPolynomial;
    stringToEnumMap_["Morse"]      = btMorse;
    stringToEnumMap_["ShiftedMie"] = btShiftedMie;
  }

  BondType* BondTypeParser::parseTypeAndPars(const std::string& type,
                                             std::vector<RealType> pars) {
    std::string line(type);

    std::vector<RealType>::iterator it;
    for (it = pars.begin(); it != pars.end(); ++it) {
      line.append("\t");
      line.append(OpenMD::to_string(*it));
    }
    // assume all overrides know about our functional forms:
    return parseLine(line, 1.0);
  }

  BondType* BondTypeParser::parseLine(const std::string& line,
                                      RealType kScale) {
    StringTokenizer tokenizer(line);
    BondType* bondType = NULL;
    int nTokens        = tokenizer.countTokens();

    if (nTokens < 1) {
      throw OpenMDException("BondTypeParser: Not enough tokens");
    }

    BondTypeEnum bt = getBondTypeEnum(tokenizer.nextToken());
    nTokens -= 1;

    switch (bt) {
    case btFixed:
      if (nTokens < 1) {
        throw OpenMDException("BondTypeParser: Not enough tokens");
      } else {
        RealType b0 = tokenizer.nextTokenAsDouble();
        bondType    = new FixedBondType(b0);
      }
      break;

    case btHarmonic:
      if (nTokens < 2) {
        throw OpenMDException("BondTypeParser: Not enough tokens");
      } else {
        RealType b0 = tokenizer.nextTokenAsDouble();
        RealType kb = tokenizer.nextTokenAsDouble();
        kb *= kScale;
        bondType = new HarmonicBondType(b0, kb);
      }
      break;

    case btCubic:
      if (nTokens < 5) {
        throw OpenMDException("BondTypeParser: Not enough tokens");
      } else {
        RealType b0 = tokenizer.nextTokenAsDouble();
        RealType k3 = tokenizer.nextTokenAsDouble();
        RealType k2 = tokenizer.nextTokenAsDouble();
        RealType k1 = tokenizer.nextTokenAsDouble();
        RealType k0 = tokenizer.nextTokenAsDouble();

        bondType = new CubicBondType(b0, k3, k2, k1, k0);
      }
      break;

    case btQuartic:
      if (nTokens < 6) {
        throw OpenMDException("BondTypeParser: Not enough tokens");
      } else {
        RealType b0 = tokenizer.nextTokenAsDouble();
        RealType k4 = tokenizer.nextTokenAsDouble();
        RealType k3 = tokenizer.nextTokenAsDouble();
        RealType k2 = tokenizer.nextTokenAsDouble();
        RealType k1 = tokenizer.nextTokenAsDouble();
        RealType k0 = tokenizer.nextTokenAsDouble();

        bondType = new QuarticBondType(b0, k4, k3, k2, k1, k0);
      }
      break;

    case btPolynomial:

      if (nTokens < 3 || nTokens % 2 != 1) {
        throw OpenMDException("BondTypeParser: Not enough tokens");
      } else {
        RealType b0 = tokenizer.nextTokenAsDouble();
        nTokens -= 1;

        int nPairs = nTokens / 2;
        int power;
        RealType coefficient;
        PolynomialBondType* pbt = new PolynomialBondType(b0);

        for (int i = 0; i < nPairs; ++i) {
          power       = tokenizer.nextTokenAsInt();
          coefficient = tokenizer.nextTokenAsDouble();
          pbt->setCoefficient(power, coefficient);
        }
      }

      break;

    case btMorse:
      if (nTokens < 3) {
        throw OpenMDException("BondTypeParser: Not enough tokens");
      } else {
        RealType b0   = tokenizer.nextTokenAsDouble();
        RealType D    = tokenizer.nextTokenAsDouble();
        RealType beta = tokenizer.nextTokenAsDouble();
        bondType      = new MorseBondType(b0, D, beta);
      }
      break;

    case btShiftedMie:
      if (nTokens < 4) {
        throw OpenMDException("BondTypeParser: Not enough tokens");
      } else {
        RealType sigma   = tokenizer.nextTokenAsDouble();
        RealType epsilon = tokenizer.nextTokenAsDouble();
        int nRep         = tokenizer.nextTokenAsInt();
        int mAtt         = tokenizer.nextTokenAsInt();

        bondType = new ShiftedMieBondType(sigma, epsilon, nRep, mAtt);
      }
      break;

    case btUnknown:
    default:
      throw OpenMDException("BondTypeParser: Unknown Bond Type");
    }

    return bondType;
  }

  BondTypeParser::BondTypeEnum BondTypeParser::getBondTypeEnum(
      const std::string& str) {
    std::map<std::string, BondTypeEnum>::iterator i;
    i = stringToEnumMap_.find(str);

    return i == stringToEnumMap_.end() ? btUnknown : i->second;
  }

}  // namespace OpenMD
