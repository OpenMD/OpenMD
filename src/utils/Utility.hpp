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

#ifndef UTILS_UTILITY_HPP
#define UTILS_UTILITY_HPP

#include <config.h>

#include <cmath>
#include <string>
#include <vector>

namespace OpenMD {
  inline RealType roundMe(const RealType& x) {
    return (x >= 0) ? std::floor(x + 0.5) : std::ceil(x - 0.5);
  }

  /**
   * @brief iteratively replace the sequence with wild cards
   * @return true if more combination sequence is avaliable, otherwise
   * return true
   * @param cont iterator container, if expect the whole series of
   * combinations, pass an empty iterator container. The user should
   * not modify this iterator container
   * @param sequence the whole sequence used to generate combination
   * @param result a possible combination sequence which is set on return
   * @param wildCard the wild card string. Its value is "X" by default
   * @note since next_combination never returns an empty sequence,
   * replaceWildCard will not generate one special combination, which
   * is n identical wild cards (n is equal to the size of the passing
   * sequence)
   *
   * @code
   * vector<string> sv;
   * vector<vector<string>::iterator> sic;
   * vector<string> resultString;
   * sv.push_back("H");
   * sv.push_back("C");
   * sv.push_back("N");

   * while (replaceWithWildCard(sic, sv, resultString)) {
   *     for(vector<string>::iterator i = resultString.begin();
   *          i != resultString.end(); ++i) {
   *         cout << *i << "\t";
   *     }
   *     cout << endl;
   * }
   * //output
   * //H X X
   * //X C X
   * //X X N
   * //H C X
   * //H X N
   * //X C N
   * //H C N
   * @endcode
   */
  bool replaceWithWildCard(
      std::vector<std::vector<std::string>::iterator>& cont,
      std::vector<std::string>& sequence, std::vector<std::string>& result,
      const std::string& wildCard = "X");
}  // namespace OpenMD

#endif
