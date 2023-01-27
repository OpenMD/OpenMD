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

#include "io/SectionParser.hpp"

#include "utils/Trim.hpp"

namespace OpenMD {

  void SectionParser::parse(std::istream& input, ForceField& ff, int lineNo) {
    const int bufferSize = 65535;
    char buffer[bufferSize];
    std::string line, filteredLine;
    while (input.getline(buffer, bufferSize)) {
      ++lineNo;
      line = Utils::trimLeftCopy(buffer);
      // a line begins with "//" is comment
      // let's also call lines starting with # and ! as comments
      if (isEndSection(line)) {
        break;
      } else if (line.empty() ||
                 (line.size() >= 2 && line[0] == '/' && line[1] == '/') ||
                 (line.size() >= 1 && line[0] == '#') ||
                 (line.size() >= 1 && line[0] == '!')) {
        continue;
      } else {
        filteredLine = stripComments(line);
        parseLine(ff, filteredLine, lineNo);
      }
    }
  }

  std::string SectionParser::stripComments(const std::string& line) {
    unsigned int n = line.length();
    std::string res;

    // Flags to indicate that single line and multpile line comments
    // have started or not.
    bool s_cmt = false;
    bool m_cmt = false;

    // Traverse the line
    for (unsigned int i = 0; i < n; i++) {
      // If single line comment flag is on, then check for end of it
      if (s_cmt == true && line[i] == '\n') s_cmt = false;

      // If multiple line comment is on, then check for end of it
      else if (m_cmt == true && line[i] == '*' && line[i + 1] == '/')
        m_cmt = false, i++;

      // If this character is in a comment, ignore it
      else if (s_cmt || m_cmt)
        continue;

      // Check for beginning of comments and set the approproate flags
      else if (line[i] == '/' && line[i + 1] == '/')
        s_cmt = true, i++;
      else if (line[i] == '/' && line[i + 1] == '*')
        m_cmt = true, i++;

      // If current character is a non-comment character, append it to res
      else
        res += line[i];
    }
    return res;
  }

  bool SectionParser::isEndSection(const std::string& line) {
    StringTokenizer tokenizer(line);

    if (tokenizer.countTokens() >= 2) {
      std::string keyword = tokenizer.nextToken();

      if (keyword != "end") { return false; }

      std::string section = tokenizer.nextToken();
      if (section == sectionName_) {
        return true;
      } else {
        return false;
      }

    } else {
      return false;
    }
  }

}  // namespace OpenMD
