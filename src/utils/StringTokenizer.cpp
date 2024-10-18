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

#include "utils/StringTokenizer.hpp"

#include <iostream>
#include <iterator>
#include <sstream>

namespace OpenMD {

  StringTokenizer::StringTokenizer(const std::string& str,
                                   const std::string& delim) :
      tokenString_(str),
      delim_(delim), returnTokens_(false), currentPos_(tokenString_.begin()),
      end_(tokenString_.end()) {}

  StringTokenizer::StringTokenizer(std::string::const_iterator& first,
                                   std::string::const_iterator& last,
                                   const std::string& delim) :
      tokenString_(first, last),
      delim_(delim), returnTokens_(false), currentPos_(tokenString_.begin()),
      end_(tokenString_.end()) {}

  StringTokenizer::StringTokenizer(const std::string& str,
                                   const std::string& delim,
                                   bool returnTokens) :
      tokenString_(str),
      delim_(delim), returnTokens_(returnTokens),
      currentPos_(tokenString_.begin()), end_(tokenString_.end()) {}

  bool StringTokenizer::isDelimiter(const char c) {
    return delim_.find(c) == std::string::npos ? false : true;
  }

  int StringTokenizer::countTokens() {
    std::string::const_iterator tmpIter = currentPos_;
    int numToken                        = 0;

    while (true) {
      // skip delimiter first
      while (tmpIter != end_ && isDelimiter(*tmpIter)) {
        ++tmpIter;

        if (returnTokens_) {
          // if delimiter is consider as token
          ++numToken;
        }
      }

      if (tmpIter == end_) { break; }

      // encount a token here
      while (tmpIter != end_ && !isDelimiter(*tmpIter)) {
        ++tmpIter;
      }

      ++numToken;
    }

    return numToken;
  }

  bool StringTokenizer::hasMoreTokens() {
    if (currentPos_ == end_) {
      return false;
    } else if (returnTokens_) {
      return true;
    } else {
      std::string::const_iterator i = currentPos_;

      // walk through the remaining string to check whether it contains
      // non-delimeter or not
      while (i != end_ && isDelimiter(*i)) {
        ++i;
      }

      return i != end_ ? true : false;
    }
  }

  std::string StringTokenizer::nextToken() {
    std::string result;

    if (currentPos_ != end_) {
      std::insert_iterator<std::string> insertIter(result, result.begin());

      while (currentPos_ != end_ && isDelimiter(*currentPos_)) {
        if (returnTokens_) {
          *insertIter++ = *currentPos_++;
          return result;
        }

        ++currentPos_;
      }

      while (currentPos_ != end_ && !isDelimiter(*currentPos_)) {
        *insertIter++ = *currentPos_++;
      }
    }

    return result;
  }

  void StringTokenizer::skipToken() {
    if (currentPos_ != end_) {
      while (currentPos_ != end_ && isDelimiter(*currentPos_)) {
        if (returnTokens_) {
          currentPos_++;
          return;
        }

        ++currentPos_;
      }

      while (currentPos_ != end_ && !isDelimiter(*currentPos_)) {
        currentPos_++;
      }
    }
  }

  bool StringTokenizer::nextTokenAsBool() {
    std::string token = nextToken();
    std::istringstream iss(token);
    bool result;

    if (iss >> result) {
      return result;
    } else {
      std::cerr << "unable to convert " << token << " to a bool" << std::endl;
      return false;
    }
  }

  // Since libstdc++(GCC 3.2) has an i/ostream::operator>>/<<(streambuf*) bug
  // (Bug 9318) Instead of using iostream facility, we use C library
  int StringTokenizer::nextTokenAsInt() {
    std::string token = nextToken();

    return atoi(token.c_str());
  }

  float StringTokenizer::nextTokenAsFloat() {
    std::string token = nextToken();
    convertFortranNumber(token);
    return (float)(atof(token.c_str()));
  }

  RealType StringTokenizer::nextTokenAsDouble() {
    std::string token = nextToken();
    convertFortranNumber(token);
    return atof(token.c_str());
  }

  std::string StringTokenizer::peekNextToken() {
    std::string result;
    std::string::const_iterator tmpIter = currentPos_;

    if (tmpIter != end_) {
      std::insert_iterator<std::string> insertIter(result, result.begin());

      while (tmpIter != end_ && isDelimiter(*tmpIter)) {
        if (returnTokens_) {
          *insertIter++ = *tmpIter++;
          return result;
        }

        ++tmpIter;
      }

      while (tmpIter != end_ && !isDelimiter(*tmpIter)) {
        *insertIter++ = *tmpIter++;
      }
    }

    return result;
  }

  std::vector<std::string> StringTokenizer::getAllTokens() {
    std::vector<std::string> tokens;
    while (hasMoreTokens()) {
      tokens.push_back(nextToken());
    }
    return tokens;
  }

  void StringTokenizer::convertFortranNumber(std::string& fortranNumber) {
    std::string::iterator i;
    for (i = fortranNumber.begin(); i != fortranNumber.end(); ++i) {
      if (*i == 'd' || *i == 'D') { *i = 'E'; }
    }
  }

  std::string StringTokenizer::getRemainingString() const {
    std::string result;
    std::string::const_iterator tmpIter = currentPos_;
    if (tmpIter != end_) {
      std::insert_iterator<std::string> insertIter(result, result.begin());

      while (tmpIter != end_) {
        *insertIter++ = *tmpIter++;
      }
    }

    return result;
  }
}  // namespace OpenMD
