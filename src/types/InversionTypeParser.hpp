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

#ifndef TYPES_INVERSIONTYPEPARSER_HPP
#define TYPES_INVERSIONTYPEPARSER_HPP

#include <map>
#include <string>
#include <vector>

#include "types/InversionType.hpp"

namespace OpenMD {

  /**
   * @class InversionTypeParser InversionTypeParser.hpp
   * "types/InversionTypeParser.hpp"
   */
  class InversionTypeParser {
  public:
    InversionTypeParser();
    InversionType* parseLine(const std::string& line);
    InversionType* parseTypeAndPars(const std::string& type,
                                    std::vector<RealType> pars);

  private:
    // Inversion types vary by force field:  In this description,
    // I is the central atom, while J, K, and L are atoms directly
    // bonded to I (but not to each other):

    // Amber uses a special bond (IL) as the hinge between the planes
    // IJL and IKL (the central atom I & peripheral atom L are common
    // in both planes).  It then applies a cosine series much like
    // other torsional forms.

    // Gromacs & Charmm use an improper torsion that is harmonic
    // in the angle between the IJK and JKL planes (Central atom is I, but
    // this atom appears in only one of the plane definitions.

    // MM2 and Tripos use a planarity definition of the central atom (I)
    // distance from the plane made by the other three atoms (JKL).

    // The Dreiding force field uses a complicated umbrella inversion
    // form.

    enum InversionTypeEnum {
      itImproperCosine,
      itHarmonic,
      itCentralAtomHeight,
      itAmberImproper,
      itDreiding,
      itUnknown
    };

    InversionTypeEnum getInversionTypeEnum(const std::string& str);
    std::map<std::string, InversionTypeEnum> stringToEnumMap_;
  };
}  // namespace OpenMD

#endif
