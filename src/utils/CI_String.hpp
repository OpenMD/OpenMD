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

#ifndef OPENMD_UTILS_CI_STRING_HPP
#define OPENMD_UTILS_CI_STRING_HPP

#include <cctype>
#include <cstddef>
#include <iostream>
#include <string>

namespace OpenMD::Utils {

  struct ci_char_traits : public std::char_traits<char> {
    static int compare(const char* s1_, const char* s2_,
                       std::size_t count_) noexcept {
      while (count_-- != 0) {
        if (std::toupper(*s1_) < std::toupper(*s2_)) return -1;
        if (std::toupper(*s1_) > std::toupper(*s2_)) return 1;

        ++s1_;
        ++s2_;
      }

      return 0;
    }

    static const char* find(const char* p_, std::size_t count_,
                            const char& ch_) noexcept {
      const auto CH {std::toupper(ch_)};

      while (count_-- != 0) {
        if (std::toupper(*p_) == CH) return p_;

        ++p_;
      }

      return nullptr;
    }

    static bool eq(char a, char b) noexcept {
      return std::toupper(a) == std::toupper(b);
    }

    static bool lt(char a, char b) noexcept {
      return std::toupper(a) < std::toupper(b);
    }
  };

  template<class OutputTraits, class InputTraits>
  inline std::basic_string<char, OutputTraits> traits_cast(
      const std::basic_string<char, InputTraits>& input) noexcept {
    std::basic_string<char, OutputTraits> ouput {input.data(), input.size()};
    return ouput;
  }

  inline std::ostream& operator<<(
      std::ostream& os, const std::basic_string<char, ci_char_traits>& str) {
    os << str.c_str();
    return os;
  }

  // Type alias for case-insensitive strings
  using ci_string = std::basic_string<char, ci_char_traits>;
}  // namespace OpenMD::Utils

#endif
