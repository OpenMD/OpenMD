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

#include "io/GayBerneAtomTypesSectionParser.hpp"

#include "brains/ForceField.hpp"
#include "types/AtomType.hpp"
#include "types/GayBerneAdapter.hpp"
#include "utils/simError.h"

namespace OpenMD {

  GayBerneAtomTypesSectionParser::GayBerneAtomTypesSectionParser(
      ForceFieldOptions& options) :
      options_(options) {
    setSectionName("GayBerneAtomTypes");
  }

  void GayBerneAtomTypesSectionParser::parseLine(ForceField& ff,
                                                 const std::string& line,
                                                 int lineNo) {
    StringTokenizer tokenizer(line);
    int nTokens = tokenizer.countTokens();
    // in GayBerneAtomTypesSection, a line at least contains 7 tokens
    // atomTypeName   d    l    eps_X  eps_S   eps_E   dw
    if (nTokens < 7) {
      snprintf(
          painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
          "GayBerneAtomTypesSectionParser: Not enough tokens at line %d\n"
          "\tPlease note that GB atoms now require separate specification of "
          "epsilon\n"
          "\tvalues for cross (X), Side-by-Side (S), and End-to-End (E) for "
          "each\n"
          "\tellipsoid type.\n",
          lineNo);
      painCave.isFatal = 1;
      simError();
    } else {
      RealType eus_            = options_.getEnergyUnitScaling();
      RealType dus_            = options_.getDistanceUnitScaling();
      std::string atomTypeName = tokenizer.nextToken();
      AtomType* atomType       = ff.getAtomType(atomTypeName);
      if (atomType != NULL) {
        RealType GB_d     = dus_ * tokenizer.nextTokenAsDouble();
        RealType GB_l     = dus_ * tokenizer.nextTokenAsDouble();
        RealType GB_eps_X = eus_ * tokenizer.nextTokenAsDouble();
        RealType GB_eps_S = eus_ * tokenizer.nextTokenAsDouble();
        RealType GB_eps_E = eus_ * tokenizer.nextTokenAsDouble();
        RealType GB_dw    = tokenizer.nextTokenAsDouble();

        GayBerneAdapter gba = GayBerneAdapter(atomType);
        gba.makeGayBerne(GB_d, GB_l, GB_eps_X, GB_eps_S, GB_eps_E, GB_dw);

      } else {
        snprintf(painCave.errMsg, MAX_SIM_ERROR_MSG_LENGTH,
                 "GayBerneAtomTypesSectionParser: Can not find matching "
                 "AtomType %s\n"
                 "\tfor this GayBerne atom type\n",
                 atomTypeName.c_str());
        painCave.isFatal = 1;
        simError();
      }
    }
  }
}  // namespace OpenMD
