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

#include "io/UFFAtomTypesSectionParser.hpp"

#include "brains/ForceField.hpp"
#include "types/UFFAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  UFFAtomTypesSectionParser::UFFAtomTypesSectionParser(
      ForceFieldOptions& options) :
      options_(options) {
    setSectionName("UFFAtomTypes");
  }

  void UFFAtomTypesSectionParser::parseLine(ForceField& ff,
                                            const std::string& line,
                                            int lineNo) {
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();

    // in UFFAtomTypesSectionParser, a line at least contains 12 tokens:
    // Atom r1 theta0 x1 D1 zeta Z1 Vi Uj Xi Hard Radius

    if (nTokens < 12) {
      snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
               "UFFAtomTypesSectionParser Error: "
               "Not enough tokens at line %d\n",
               lineNo);
      painCave.isFatal = 1;
      simError();
    } else {
      std::string atomTypeName = tokenizer.nextToken();
      AtomType* atomType       = ff.getAtomType(atomTypeName);

      if (atomType != NULL) {
        UFFAdapter uffa = UFFAdapter(atomType);

        RealType r1     = tokenizer.nextTokenAsDouble();
        RealType theta0 = tokenizer.nextTokenAsDouble();
        RealType x1     = tokenizer.nextTokenAsDouble();
        RealType D1     = tokenizer.nextTokenAsDouble();
        RealType zeta   = tokenizer.nextTokenAsDouble();
        RealType Z1     = tokenizer.nextTokenAsDouble();
        RealType Vi     = tokenizer.nextTokenAsDouble();
        RealType Uj     = tokenizer.nextTokenAsDouble();
        RealType Xi     = tokenizer.nextTokenAsDouble();
        RealType Hard   = tokenizer.nextTokenAsDouble();
        RealType Radius = tokenizer.nextTokenAsDouble();

        r1 *= options_.getDistanceUnitScaling();
        theta0 *= options_.getAngleUnitScaling();
        x1 *= options_.getDistanceUnitScaling();
        D1 *= options_.getEnergyUnitScaling();
        Z1 *= options_.getChargeUnitScaling();
        Vi *= options_.getEnergyUnitScaling();
        Uj *= options_.getEnergyUnitScaling();

        uffa.makeUFF(r1, theta0, x1, D1, zeta, Z1, Vi, Uj, Xi, Hard, Radius);

      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "UFFAtomTypesSectionParser Error: Atom Type [%s]"
                 " has not been created yet\n",
                 atomTypeName.c_str());
        painCave.isFatal = 1;
        simError();
      }
    }
  }

  void UFFAtomTypesSectionParser::validateSection(ForceField&) {
    // ForceField::AtomTypeContainer* atomTypes = ff.getAtomTypes();
    // ForceField::AtomTypeContainer::MapTypeIterator i;
    // AtomType* at;

    // std::vector<AtomType*> uffTypes;

    // for (at = atomTypes->beginType(i); at != NULL; at =
    // atomTypes->nextType(i))
    // {

    //   UFFAdapter uffa = UFFAdapter(at);
    //   if (uffa.isUFF())
    //     uffTypes.push_back(at);
    // }

    // for (std::size_t i = 0; i < uffTypes.size(); ++i) {
    //   RealType ri = uffTypes[i]->getR1();
    //   RealType chiI = uffTypes[i]->getXi();
    //   RealType Ra = uffTypes[i]->getX1();
    //   RealType ka = uffTypes[i]->getD1();

    //   for (std::size_t j = 0; j < uffTypes.size(); ++j) {
    //     RealType rj = uffTypes[j]->getR1();
    //     RealType chiJ = uffTypes[j]->getXi();
    //     RealType Rb = uffTypes[j]->getX1();
    //     RealType kb = uffTypes[j]->getD1();

    //     // Precompute the equilibrium bond distance
    //     // From equation 3
    //     rbo = -0.1332*(ri+rj)*log(bondorder);
    //     // From equation 4
    //     ren = ri*rj*(pow((sqrt(chiI) - sqrt(chiJ)),2.0)) / (chiI*ri +
    //     chiJ*rj);
    //     // From equation 2
    //     // NOTE: See http://towhee.sourceforge.net/forcefields/uff.html
    //     // There is a typo in the published paper
    //     rij = ri + rj + rbo - ren;

    //     kab = sqrt(ka * kb);

    //     // ka now represents the xij in equation 20 -- the expected vdw
    //     distance kaSquared = (Ra * Rb); ka = sqrt(kaSquared);
  }
}  // namespace OpenMD
