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

#include "utils/StringUtils.hpp"

#include <config.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <string>

#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strdup _strdup
#define strtoull _strtoui64
#endif

namespace OpenMD {
  std::string stripComments(const std::string& S) {
    unsigned int n = S.length();
    std::string res;

    // Flags to indicate that single line and multpile line comments
    // have started or not.
    bool s_cmt = false;
    bool m_cmt = false;

    // Traverse the line
    for (unsigned int i = 0; i < n; i++) {
      // If single line comment flag is on, then check for end of it
      if (s_cmt == true && S[i] == '\n') s_cmt = false;

      // If multiple line comment is on, then check for end of it
      else if (m_cmt == true && S[i] == '*' && S[i + 1] == '/')
        m_cmt = false, i++;

      // If this character is in a comment, ignore it
      else if (s_cmt || m_cmt)
        continue;

      // Check for beginning of comments and set the approproate flags
      else if ((S[i] == '/' && S[i + 1] == '/') || (S[i] == '#'))
        s_cmt = true, i++;
      else if (S[i] == '/' && S[i + 1] == '*')
        m_cmt = true, i++;

      // If current character is a non-comment character, append it to res
      else
        res += S[i];
    }
    return res;
  }

  std::string UpperCase(const std::string& S) {
    std::string uc = S;
    unsigned int n = uc.size();
    for (unsigned int j = 0; j < n; j++) {
      char sj = uc[j];
      if (sj >= 'a' && sj <= 'z') uc[j] = (char)(sj - ('a' - 'A'));
    }
    return uc;
  }

  int findBegin(std::istream& theStream, const char* startText) {
    const int MAXLEN = 1024;
    char readLine[MAXLEN];
    int foundText = 0;
    int lineNum;
    char* the_token;

    // rewind the stream
    theStream.seekg(0, std::ios::beg);
    lineNum = 0;

    if (!theStream.eof()) {
      theStream.getline(readLine, MAXLEN);
      lineNum++;
    } else {
      printf("Error fast forwarding stream: stream is empty.\n");
      return -1;
    }

    while (!foundText) {
      if (theStream.eof()) {
        printf("Error fast forwarding stream at line %d: "
               "stream ended unexpectedly.\n",
               lineNum);
        return -1;
      }

      the_token = strtok(readLine, " ,;\t");
      if (the_token != NULL)
        if (!strcasecmp("begin", the_token)) {
          the_token = strtok(NULL, " ,;\t");
          if (the_token != NULL) {
            foundText = !strcasecmp(startText, the_token);
          }
        }

      if (!foundText) {
        if (!theStream.eof()) {
          theStream.getline(readLine, MAXLEN);
          lineNum++;
        } else {
          printf("Error fast forwarding stream at line %d: "
                 "stream ended unexpectedly.\n",
                 lineNum);
          return -1;
        }
      }
    }
    return lineNum;
  }

  int isEndLine(char* line) {
    char* working_line;
    char* foo;

    working_line = strdup(line);

    foo = strtok(working_line, " ,;\t");

    if (foo != NULL) {
      if (!strcasecmp(foo, "end")) {
        free(working_line);
        return 1;
      }
    }

    free(working_line);
    return 0;
  }

  std::string OpenMD_itoa(int value, unsigned int base) {
    const char digitMap[] = "0123456789abcdef";
    std::string buf;

    if (base == 0 || base > 16) { return buf; }

    if (value == 0) {
      buf = "0";
      return buf;
    }

    // Take care negative int:

    std::string sign;
    int _value = value;
    if (value < 0) {
      _value = -value;
      sign   = "-";
    }

    // Translating number to string with base:
    for (int i = 30; _value && i; --i) {
      buf = digitMap[_value % base] + buf;
      _value /= base;
    }
    return sign.append(buf);
  }

  std::string getPrefix(const std::string& str) {
    return str.substr(0, str.rfind('.'));
  }

  bool isInteger(const std::string& str) {
    bool result = false;

    std::string::const_iterator i = str.begin();
    if (i != str.end() && (*i == '+' || *i == '-' || std::isdigit(*i))) {
      ++i;
      while (i != str.end() && std::isdigit(*i))
        ++i;
      if (i == str.end()) result = true;
    }

    return result;
  }

  bool CaseInsensitiveEquals(const char ch1, const char ch2) {
    return std::toupper((unsigned char)ch1) == std::toupper((unsigned char)ch2);
  }

  size_t CaseInsensitiveFind(const std::string& str1, const std::string& str2) {
    std::string::const_iterator pos =
        std::search(str1.begin(), str1.end(), str2.begin(), str2.end(),
                    CaseInsensitiveEquals);
    if (pos == str1.end())
      return std::string::npos;
    else
      return pos - str1.begin();
  }
}  // namespace OpenMD
