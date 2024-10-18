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

#ifndef OPENMD_UTILS_TRIM_HPP
#define OPENMD_UTILS_TRIM_HPP

#include <algorithm>
#include <locale>
#include <string>

namespace OpenMD::Utils {

  /**
   * Remove all leading spaces in-place.
   * @param str An input sequence
   * @param loc The locale to use, defaults to std::locale
   */
  inline void trimLeft(std::string& str,
                       const std::locale& loc = std::locale()) {
    str.erase(str.begin(),
              std::find_if(str.begin(), str.end(),
                           [&loc](auto ch) { return !std::isspace(ch, loc); }));
  }

  /**
   * Remove all trailing spaces in-place.
   * @param str An input sequence
   * @param loc The locale to use, defaults to std::locale
   */
  inline void trimRight(std::string& str,
                        const std::locale& loc = std::locale()) {
    str.erase(std::find_if(str.rbegin(), str.rend(),
                           [&loc](auto ch) { return !std::isspace(ch, loc); })
                  .base(),
              str.end());
  }

  /**
   * Remove all leading and trailing spaces in-place
   * @param str An input sequence
   * @param loc The locale to use, defaults to std::locale
   */
  inline void trim(std::string& str, const std::locale& loc = std::locale()) {
    trimLeft(str, loc);
    trimRight(str, loc);
  }

  /**
   * Remove all leading spaces from the input.
   * @return A trimmed copy of the input
   * @param str An input sequence
   * @param loc The locale to use, defaults to std::locale
   */
  inline std::string trimLeftCopy(std::string str,
                                  const std::locale& loc = std::locale()) {
    trimLeft(str, loc);
    return str;
  }

  /**
   * Remove all trailing spaces from the input.
   * @return A trimmed copy of the input
   * @param str An input sequence
   * @param loc The locale to use, defaults to std::locale
   */
  inline std::string trimRightCopy(std::string str,
                                   const std::locale& loc = std::locale()) {
    trimRight(str, loc);
    return str;
  }

  /**
   * Remove all leading and trailing spaces from the input.
   * @return A trimmed copy of the input
   * @param str An input sequence
   * @param loc The locale to use, defaults to std::locale
   */
  inline std::string trimCopy(std::string str,
                              const std::locale& loc = std::locale()) {
    trim(str, loc);
    return str;
  }
}  // namespace OpenMD::Utils

#endif  // OPENMD_UTILS_TRIM_HPP
