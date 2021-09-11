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
