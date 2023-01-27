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

/**
 * @file StringUtils.hpp
 * @author Dan Gezelter
 * @date 10/18/2004
 * @version 1.0
 */

#ifndef UTILS_STRINGUTILS_HPP
#define UTILS_STRINGUTILS_HPP

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace OpenMD {

  /**
   * Converts a string to UPPER CASE
   * @param S
   */
  std::string UpperCase(const std::string& S);

  /**
   * Finds the location of the string "begin <startText>" in an input stream.
   * @param theStream
   * @param startText
   *
   * @return the line number of the block within the theStream
   */
  int findBegin(std::istream& theStream, const char* startText);

  /**
   * discovers whether or not the line contains the "end" token
   *
   * @param line  The line to test
   *
   * @return int  (==1 if the line has "end", ==0 if not).
   */
  int isEndLine(char* line);

  bool CaseInsensitiveEquals(char ch1, char ch2);

  size_t CaseInsensitiveFind(const std::string& str1, const std::string& str2);

  /**
   * Convert a variable to a string
   * @tparam T data type
   * @param v data to be converted
   * @return a string
   */
  template<typename T>
  std::string toString(const T& v) {
    std::ostringstream oss;
    if (!(oss << v)) { std::cerr << "toString Error" << std::endl; }
    return oss.str();
  }

  template<typename T>
  T lexi_cast(const std::string& str) {
    T result;
    std::istringstream iss(str);
    if (!(iss >> result)) { std::cerr << "lexi_cast Error" << std::endl; }
    return result;
  }

  template<typename T>
  bool isType(const std::string& str) {
    T result;
    std::istringstream iss(str);
    bool ret = true;
    if (!(iss >> result)) { ret = false; }
    return ret;
  }

  bool isInteger(const std::string& str);

  std::string OpenMD_itoa(int value, unsigned int base = 10);

  /**@todo need implementation */
  std::string getPrefix(const std::string& str);

  template<class ContainerType>
  std::string containerToString(const ContainerType& cont) {
    std::ostringstream oss;
    oss << "(";
    typename ContainerType::const_iterator i = cont.begin();
    if (i != cont.end()) {
      oss << *i;
      ++i;
    }
    for (; i != cont.end(); ++i) {
      oss << ", ";
      oss << *i;
    }
    oss << ")";
    return oss.str();
  }
}  // namespace OpenMD

#endif
